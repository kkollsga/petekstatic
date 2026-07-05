//! [`StaticModelTemplate`] — the Monte-Carlo regeneration seam (SPEC §7a, graph
//! `decision_staticmodel_regen_seam`, RATIFIED 2026-07-03 with peteksim's
//! amendments). Built once, it holds the **reusable** state — fixed lattice
//! geometry, control topology, layering scheme, and the warm-start surface — so
//! [`StaticModelTemplate::realize`] is the *cheap* per-realization call: a
//! warm-started solve from the previous converged surface instead of a cold
//! 20k-sweep relaxation (measured ~14× at 50×50).
//!
//! ## Kernel constraint (`decision_gridder_kernel_unification`)
//! The warm-start chain lives **entirely in petekTools kernel space**: the
//! template bootstraps its base surface with [`KernelSurface::flat`] (a constant
//! field is a fixed point of both kernels) and every subsequent surface comes
//! from `solve_surface_seeded` (the petekTools kernel). A cold `solve_surface`
//! output cannot enter the chain — the [`KernelSurface`] newtype makes that a
//! compile error. Parity gates at this seam are per-node, not aggregate.
//!
//! ## Ratified amendments honoured here
//! 1. `StaticModelTemplate: Send` + `StaticModel: Send` (Sync not required) —
//!    compile-checked in `lib.rs`; consumers shard one-template-per-worker.
//! 2. [`RealizationDraw`]: `#[non_exhaustive]` + `::new`, `pub seed_index: u64`,
//!    `Clone + Debug`, typed structural `Option`, **fvf excluded** (see `draw.rs`).
//! 3. Perf: fixed geometry stays in the template; grid-per-realization is its own
//!    budget (the 422µs/100k figure is the sampler's only). The additive
//!    [`StaticModelTemplate::realize_into`] buffer-recycling variant is **now
//!    implemented** (`task_petekstatic_realize_into`): it recycles the passed
//!    model's ZCORN + cube allocations (~100 MB/draw at 1M cells) in place, so the
//!    MC drivers (`run_structured_mc`/`_parallel`) allocate once per shard, not per
//!    draw — [`StaticModelTemplate::realize`] is `new-empty-then-realize_into`.
//! 4. `&mut self` on `realize` (the warm-start chain is serial per template) and
//!    `Result` per draw (a warm-started solve can fail to converge).

use crate::builder::{
    condition_scatter, repair_interface, resolve_stack_surfaces, spacing, stack_boundary_ring,
    stack_framework_from_names, substitute_tie_datums, validate_stack, warm_surface,
    wireframe_structure, BuildOpts, HorizonSource, HorizonStack, StackFrame, WellTie,
};
use crate::draw::{PerturbationField, RealizationDraw};
use crate::model::{Georef, StaticModel};
use crate::pipeline::{McMode, PropertyPipeline, PropertyReport};
use crate::population::{override_zone_priors, populate, PetroSample};
use crate::provenance::{
    BuildWarning, InterfaceRepair, PopulationMode, Provenance, StackProvenance, WellTieRecord,
    ZoneProvenance,
};
use crate::spec::BuildSpec;
use crate::trend::TrendSurface;
use crate::zones::ZoneTable;
use petekstatic_error::StaticError;
use petektools::geostat::sgs_unconditional;
use petektools::Lattice;
use srs_grid::{Grid, Property};
use srs_gridder::ExtrapolationPolicy;
use srs_gridder::{
    layer_grid_stack_into, solve_surface_seeded, Conformity, Control, KernelSurface, LayerScratch,
    SolveOpts, Surface, ZoneLayerSpec,
};
use srs_volumetrics::{validate_fraction, validate_positive, ConstantPriors, NTG, PORO, SW};
use srs_wireframe::{Contact, ContactKind, Hardness, Wireframe};
use std::borrow::Cow;

/// Bounded SGS search-neighbourhood node cap for the structural perturbation fields
/// (`decision_structural_uncertainty_isochore`) — the same default the property
/// pipelines use (`pipeline::DEFAULT_MAX_NEIGHBOURS`).
const STRUCTURAL_MAX_NEIGHBOURS: usize = 16;

/// Golden-ratio seed mixer (shared convention with `pipeline::SEED_GOLDEN`), so a
/// per-horizon perturbation seed is independent yet reproducible.
const SEED_GOLDEN: u64 = 0x9E37_79B9_7F4A_7C15;

/// A per-horizon structural-perturbation seed: `draw.seed_index` salted by the
/// horizon index (0 = top depth field, `z + 1` = zone `z`'s isochore field). Mutually
/// independent across horizons, bit-reproducible per draw seed.
fn horizon_seed(base: u64, salt: u64) -> u64 {
    base ^ salt.wrapping_add(1).wrapping_mul(SEED_GOLDEN)
}

/// Generate a node-major (`jp * nx + ip`) correlated structural perturbation field on
/// the areal node lattice for one horizon of a draw
/// (`decision_structural_uncertainty_isochore`): an unconditional Gaussian random
/// field (petekTools [`sgs_unconditional`]) with marginal `N(0, sd_m²)` and the
/// variogram's spatial continuity. `sd_m <= 0` (or non-finite) → an all-zero field
/// (a no-op). The search neighbourhood is derived from the variogram range (radius
/// `= max(range·1.5, spacing·4)`), matching the property-SGS default.
#[allow(clippy::too_many_arguments)]
fn perturbation_field(
    field: &PerturbationField,
    origin_x: f64,
    origin_y: f64,
    dx: f64,
    dy: f64,
    nx: usize,
    ny: usize,
    seed: u64,
) -> Result<Vec<f64>, StaticError> {
    let n = nx * ny;
    if !(field.sd_m.is_finite() && field.sd_m > 0.0) {
        return Ok(vec![0.0; n]);
    }
    let lattice = Lattice::regular(origin_x, origin_y, dx, dy, nx, ny);
    let spacing = dx.max(dy);
    let radius = (field.variogram.range * 1.5).max(spacing * 4.0);
    let arr = sgs_unconditional(
        &lattice,
        0.0,
        field.sd_m * field.sd_m,
        &field.variogram,
        STRUCTURAL_MAX_NEIGHBOURS,
        radius,
        seed,
    )
    .map_err(|e| StaticError::Grid(format!("structural perturbation field failed: {e}")))?;
    // `sgs_unconditional` returns an `(ncol, nrow)` array indexed `[ip, jp]`; flatten
    // to the node-major `jp * nx + ip` layout `Surface::offset_by_field` expects.
    let mut out = vec![0.0; n];
    for jp in 0..ny {
        for ip in 0..nx {
            out[jp * nx + ip] = arr[[ip, jp]];
        }
    }
    Ok(out)
}

/// Pin a perturbation field to **zero at every well-tie node** whose record names
/// `horizon` (`decision_structural_uncertainty_isochore`: perturbation → 0 at tie
/// wells). Radius-0 locality — a tie fixes exactly its node (matching the default
/// [`crate::TieMethod::Replace`]), leaving the correlated field to relax away over
/// the variogram range around it. Ordering is preserved (a zeroed thickness node
/// stays at its non-negative base isochore).
fn pin_ties(field: &mut [f64], well_ties: &[WellTieRecord], horizon: &str, nx: usize, ny: usize) {
    for rec in well_ties {
        if rec.ip < nx && rec.jp < ny && rec.residuals.iter().any(|r| r.horizon == horizon) {
            field[rec.jp * nx + rec.ip] = 0.0;
        }
    }
}

/// The **fixed** resolved state of a stack-aware template
/// ([`StaticModelTemplate::from_horizon_stack`]): the horizon surfaces resolved once
/// (their controls do not vary per draw, so nor do the surfaces), the per-zone
/// layering specs + names, the static per-zone contacts, the once-built framework
/// (Arc-shared horizons; only its contacts vary per draw), and the interface-repair
/// record. A stack realization varies only the areal footprint (spacing), the
/// per-zone contacts, and the per-zone property levels — never the geometry topology
/// — so the surfaces here are shared across every draw (bit-determinism by
/// construction).
#[derive(Debug, Clone)]
struct StackState {
    surfaces: Vec<Surface>,
    horizon_names: Vec<String>,
    zone_names: Vec<String>,
    /// Per-zone display colours (viewer hint), parallel to `zone_names`; carried
    /// from the source [`crate::StackZone::color`] onto each realization's zone table.
    zone_colors: Vec<Option<String>>,
    zone_specs: Vec<ZoneLayerSpec>,
    /// Static per-zone contacts from the stack (used when a draw does not override a
    /// zone), parallel to `zone_names`.
    zone_contacts: Vec<Vec<Contact>>,
    interface_repairs: Vec<InterfaceRepair>,
    /// Areal node spacing (m) from the nominal `opts.area_m2` footprint — the
    /// metre scale for [`crate::TieSettings`] radius locality (ties are
    /// draw-invariant; a per-draw area does not re-evaluate them).
    tie_dxdy: (f64, f64),
    /// Framework built once with the union of static contacts; per draw its
    /// `contacts` are replaced (Arc horizons are shared, O(1) clone).
    framework: Wireframe,
    /// The source [`HorizonStack`] (tied data substituted in by
    /// [`StaticModelTemplate::with_well_ties`]) — retained so ties and a policy
    /// change re-run the SAME stack resolution the deterministic builder uses
    /// (cumulative isochores, merge collapse, extrapolation policy).
    source: HorizonStack,
}

/// A per-property geostatistical pipeline attached to the template, with its
/// Monte-Carlo mode (`decision_mc_composition`) and — for [`McMode::LevelShift`] —
/// the once-propagated pattern + report cached at the first realization.
#[derive(Debug, Clone)]
struct PropertyModel {
    pipeline: PropertyPipeline,
    mode: McMode,
    /// `LevelShift`: the once-propagated pattern cube + its report, computed at the
    /// first `realize` and reused (shifted) thereafter. `None` until then / for
    /// `Resimulate`.
    cached: Option<(Vec<f64>, PropertyReport)>,
}

/// A **zone-scoped** geostatistical pipeline attached to the template (the MC analog
/// of [`crate::StaticModelBuilder::with_zone_property`]): the pipeline is run
/// **restricted to its named zone's `k`-range** each realization, over the per-zone
/// constant priors, with its Monte-Carlo mode and — for [`McMode::LevelShift`] — the
/// once-propagated cube cached at the first realization. Without this the zoned MC
/// realizes zone-piped zones from the ZONE PRIORS only, ignoring the staged
/// upscale+SGS cubes (`question_zoned_mc_zone_pipe_parity`).
#[derive(Debug, Clone)]
struct ZonePropertyModel {
    /// The stack zone this pipeline populates (matched against the stack's zone names).
    zone: String,
    pipeline: PropertyPipeline,
    mode: McMode,
    /// `LevelShift`: the once-propagated full cube pattern + its report, computed at
    /// the first `realize` and reused (its zone slice shifted) thereafter. Only the
    /// zone's `k`-range cells of the pattern are ever read. `None` until then /
    /// for `Resimulate`.
    cached: Option<(Vec<f64>, PropertyReport)>,
}

/// The nominal per-node gross-thickness field extracted from a supplied `Base`
/// horizon at template build (`decision_template_gross_scaling`): `node_dz`
/// carries the real relief *shape* (row-major `jp * nx + ip`); the per-draw
/// `gross_height_m` scales its *level* by `gross_height_m / mean`.
#[derive(Debug, Clone)]
struct GrossField {
    node_dz: Vec<f64>,
    mean: f64,
}

/// The reusable MC-regeneration template. One template per worker (it is `Send`;
/// the warm-start chain makes `realize` serial per template).
#[derive(Debug, Clone)]
pub struct StaticModelTemplate {
    ni: usize,
    nj: usize,
    nk: usize,
    conformity: Conformity,
    solve_opts: SolveOpts,
    top_controls: Vec<Control>,
    /// Real base relief as a nominal gross field (`None` = no Base horizon →
    /// constant per-draw offset, the backward-compatible fallback).
    gross_field: Option<GrossField>,
    logs: Option<Vec<PetroSample>>,
    /// Optional areal trend multiplier field applied per-realization to NTG
    /// (and optionally φ) after population — lateral shape only.
    trend: Option<TrendSurface>,
    framework: Wireframe,
    /// Entry-point default provenance label (`"wireframe"` / `"horizon-stack"`),
    /// used when the spec does not name one.
    default_inputs_ref: &'static str,
    /// The declarative build configuration — the SAME [`BuildSpec`] the
    /// deterministic builder consumes (`task_petekstatic_spec_mirror`); every
    /// `with_*` setter is thin sugar mutating it, and
    /// [`StaticModelTemplate::with_spec`] installs a whole one.
    spec: BuildSpec,
    /// The warm-start chain state — always petekTools kernel space.
    seed: KernelSurface,
    /// Per-property geostatistical pipelines run each realization after the base
    /// population (P5), each with its Monte-Carlo mode (`decision_mc_composition`).
    properties: Vec<PropertyModel>,
    /// Zone-scoped geostatistical pipelines (`with_zone_property`), each run
    /// restricted to its named zone's `k`-range after the per-zone constant priors
    /// on the stack path (the MC analog of the builder's per-zone population). Empty
    /// on the classic 2-surface path.
    zone_properties: Vec<ZonePropertyModel>,
    /// Per-horizon well-tie records (`with_well_ties`) — the tie residuals recorded
    /// when the stack surfaces were tied to the wells at template construction. The
    /// geometry is already tied (draw-invariant), so this only rides onto each
    /// realization's provenance at zero per-draw cost. Empty when untied.
    well_ties: Vec<WellTieRecord>,
    /// Reusable per-pillar layering scratch for [`StaticModelTemplate::realize_into`]
    /// — grown once and refilled each draw so the layering intermediates are not
    /// re-allocated per realization. Not logical state: a template clone starts with
    /// empty scratch (see [`LayerScratch`]), so one-template-per-worker sharding does
    /// not copy it.
    scratch: LayerScratch,
    /// The resolved multi-zone stack state for a stack-aware template
    /// ([`StaticModelTemplate::from_horizon_stack`]). `Some` routes `realize_into`
    /// down the stack path (per-zone contacts + per-zone property levels over a fixed
    /// multi-horizon framework); `None` is the classic 2-surface path.
    stack: Option<StackState>,
}

impl StaticModelTemplate {
    /// Build the template from a constraining wireframe: extract the fixed
    /// structure (lattice, control topology + nominal depths, layering) and
    /// converge the warm-start base surface **in kernel space** (flat bootstrap
    /// at the mean control depth, then one seeded solve).
    ///
    /// `opts` carries the fixed geometry (nk, conformity, nominal gridder
    /// settings) and default scalars; the per-draw scalars in a
    /// [`RealizationDraw`] override area/height/contact/priors per realization.
    ///
    /// # Errors
    /// [`StaticError`] if the wireframe has no `Top` horizon / contact, the
    /// lattice is degenerate, or the bootstrap solve fails.
    pub fn new(wf: &Wireframe, opts: BuildOpts) -> Result<Self, StaticError> {
        let s = wireframe_structure(wf, &opts)?;
        let (ni, nj) = (s.ni, s.nj);
        let top_controls = s.top_controls;
        // Structural solve via the SAME kernel-space path the deterministic
        // builder uses (`warm_surface`), so a mean-gross realization reproduces
        // the build (R2).
        let seed = warm_surface(ni + 1, nj + 1, &top_controls)?;

        // A supplied Base horizon becomes the nominal gross field g(x,y): solve
        // it once in kernel space and take the per-node thickness against the
        // nominal top (`decision_template_gross_scaling`). Per draw,
        // `gross_height_m` rescales g's level around its mean; without a Base
        // the draw keeps today's constant offset.
        let gross_field = match &s.base_controls {
            None => None,
            Some(bc) => {
                let base = warm_surface(ni + 1, nj + 1, bc)?;
                let mut node_dz = Vec::with_capacity((ni + 1) * (nj + 1));
                for jp in 0..=nj {
                    for ip in 0..=ni {
                        node_dz.push(base.z(ip, jp) - seed.z(ip, jp));
                    }
                }
                let mean = node_dz.iter().sum::<f64>() / node_dz.len() as f64;
                if !(mean.is_finite() && mean > 0.0) {
                    return Err(StaticError::Grid(format!(
                        "Base horizon gross field must have a positive mean, got {mean}"
                    )));
                }
                Some(GrossField { node_dz, mean })
            }
        };

        Ok(Self {
            ni,
            nj,
            nk: opts.nk,
            conformity: opts.conformity,
            solve_opts: opts.solve_opts,
            top_controls,
            gross_field,
            logs: None,
            trend: None,
            framework: wf.clone(),
            default_inputs_ref: "wireframe",
            spec: BuildSpec::default(),
            seed,
            properties: Vec::new(),
            zone_properties: Vec::new(),
            well_ties: Vec::new(),
            scratch: LayerScratch::new(),
            stack: None,
        })
    }

    /// Build a **stack-aware** MC template from an ordered [`HorizonStack`] (P8,
    /// `task_petekstatic_multizone_2`): the multi-zone analog of
    /// [`StaticModelTemplate::new`]. It resolves the horizon surfaces **once** (their
    /// controls do not vary per draw), repairs consecutive-horizon ordering, and
    /// caches the per-zone layering + the framework. [`StaticModelTemplate::realize_into`]
    /// then varies, per draw, only the areal footprint (spacing), the **per-zone
    /// contacts** (each [`crate::ZoneDraw`] owns an optional GOC/OWC; a contactless
    /// zone contributes GRV, zero HC), and the **per-zone property levels** — never
    /// the geometry topology, so realizations are bit-deterministic by construction
    /// and `realize_into` recycles buffers exactly as on the 2-surface path.
    ///
    /// `opts` supplies the fallback priors + gridder settings; `opts.nk` /
    /// `conformity` / `gross_height_m` are unused (per-zone layering governs). The
    /// `with_min_thickness_m` / `with_clamp_base_to_top` / `with_collapse_below_m` /
    /// `with_georef` / `with_inputs_ref` builders apply as usual.
    ///
    /// # Errors
    /// [`StaticError`] as [`crate::StaticModelBuilder::from_horizon_stack`] (bad stack
    /// shape, a degenerate/mismatched surface, a tops-only horizon with no mapped
    /// horizon above, or a solve failure).
    pub fn from_horizon_stack(stack: HorizonStack, opts: BuildOpts) -> Result<Self, StaticError> {
        let (ni, nj) = validate_stack(&stack)?;
        let (nx, ny) = (ni + 1, nj + 1);

        // Resolve every horizon surface by building down via non-negative isochores
        // (the merge-safe construction; see [`resolve_stack_surfaces`]) — byte-for-byte
        // the deterministic builder's resolution, so `realize` reproduces it.
        let policy = ExtrapolationPolicy::default();
        let mut surfaces = resolve_stack_surfaces(&stack, nx, ny, policy)?;
        let is_derived: Vec<bool> = stack
            .horizons
            .iter()
            .map(|h| matches!(h.source, HorizonSource::TopsOnly(_)))
            .collect();

        // Per-interface order-repair (fixed for the template's lifetime), honouring
        // mapped-over-derived precedence. A no-op on a consistent isochore-built
        // stack; retained as the per-draw crossing guard.
        let mut interface_repairs = Vec::new();
        for i in 0..surfaces.len() - 1 {
            let (upper, lower) = repair_interface(
                &surfaces[i],
                &surfaces[i + 1],
                None, // template default: order-repair to zero thickness (no floor)
                true, // clamp a crossing rather than erroring at build (per-draw stable)
                is_derived[i],
                is_derived[i + 1],
                i,
                &mut interface_repairs,
            )?;
            surfaces[i] = upper;
            surfaces[i + 1] = lower;
        }

        let horizon_names: Vec<String> = stack.horizons.iter().map(|h| h.name.clone()).collect();
        let zone_contacts: Vec<Vec<Contact>> = stack
            .zone_layers
            .iter()
            .map(|z| z.contacts.clone())
            .collect();
        let zone_specs: Vec<ZoneLayerSpec> = stack
            .zone_layers
            .iter()
            .map(|z| ZoneLayerSpec {
                conformity: z.conformity,
                requested_nk: z.nk,
            })
            .collect();
        let all_contacts: Vec<Contact> = zone_contacts.iter().flatten().copied().collect();
        // Placeholder unit-square ring; the real world outline is resolved per
        // realization once `with_georef` / `with_boundary` are known (finding 2).
        let framework = stack_framework_from_names(
            &horizon_names,
            &surfaces,
            all_contacts,
            nx,
            ny,
            stack_boundary_ring(None, None, nx, ny),
        );

        let zone_names: Vec<String> = stack.zone_layers.iter().map(|z| z.name.clone()).collect();
        let zone_colors: Vec<Option<String>> =
            stack.zone_layers.iter().map(|z| z.color.clone()).collect();
        let stack_state = StackState {
            surfaces,
            horizon_names,
            zone_names,
            zone_colors,
            zone_specs,
            zone_contacts,
            interface_repairs,
            tie_dxdy: spacing(opts.area_m2, ni, nj),
            framework: framework.clone(),
            source: stack,
        };

        Ok(Self {
            ni,
            nj,
            nk: opts.nk,
            conformity: opts.conformity,
            solve_opts: opts.solve_opts,
            top_controls: Vec::new(),
            gross_field: None,
            logs: None,
            trend: None,
            framework,
            default_inputs_ref: "horizon-stack",
            spec: BuildSpec::default(),
            seed: KernelSurface::flat(nx, ny, 0.0),
            properties: Vec::new(),
            zone_properties: Vec::new(),
            well_ties: Vec::new(),
            scratch: LayerScratch::new(),
            stack: Some(stack_state),
        })
    }

    /// Build a **stack-aware MC template** from a stack that may carry raw scatter
    /// horizons ([`HorizonSource::Scatter`]). Each scatter horizon is conditioned
    /// onto the `frame` lattice (the engine owns the gridding; genuine voids `NaN`),
    /// then resolved through the **same** merge-safe path as
    /// [`Self::from_horizon_stack`] — byte-for-byte the deterministic
    /// [`crate::StaticModelBuilder::from_scatter_stack`], so `realize` reproduces
    /// it. The `frame` georef is registered for the realized world frames.
    ///
    /// # Errors
    /// [`StaticError`] as [`Self::from_horizon_stack`] /
    /// [`crate::StaticModelBuilder::from_scatter_stack`].
    pub fn from_scatter_stack(
        mut stack: HorizonStack,
        opts: BuildOpts,
        frame: StackFrame,
    ) -> Result<Self, StaticError> {
        condition_scatter(&mut stack, &frame)?;
        let g = frame.georef;
        Ok(Self::from_horizon_stack(stack, opts)?.with_georef(
            g.origin_x,
            g.origin_y,
            g.spacing_x,
            g.spacing_y,
        ))
    }

    /// Opt into clamping a crossed base (base above top) to zero gross at the
    /// offending columns per realization, instead of erroring. Default is to error
    /// ([`StaticError::CrossedSurfaces`]); this zeroes only those columns (R1).
    #[must_use]
    pub fn with_clamp_base_to_top(mut self, clamp: bool) -> Self {
        self.spec.clamp_base_to_top = clamp;
        self
    }

    /// Opt into **post-gridding order-repair** per realization: where the scaled
    /// base sits less than `min_thickness_m` below the top (a thin or crossed
    /// column, e.g. gridder overshoot at a thin margin), pull the base **down** to
    /// exactly `top + min_thickness_m`, preserving the top (the better-constrained
    /// seismic pick). Off by default (the crossing guard stays the default); when
    /// enabled, each realization's [`Provenance::warnings`] records a
    /// [`BuildWarning::ThinColumnsRepaired`] with the repaired-node count and the
    /// worst violation. Takes precedence over [`Self::with_clamp_base_to_top`] (R-c).
    #[must_use]
    pub fn with_min_thickness_m(mut self, min_thickness_m: f64) -> Self {
        self.spec.min_thickness_m = Some(min_thickness_m);
        self
    }

    /// Opt into the **cell-collapse threshold** applied each realization — the MC
    /// analog of [`crate::StaticModelBuilder::with_collapse_below_m`]: after
    /// layering, any cell thinner than `collapse_below_m` collapses (its sliver
    /// merged into a thicker zone-interior neighbour; volume-conserving). Off by
    /// default; a realization records [`BuildWarning::CellsCollapsed`] when it bites.
    #[must_use]
    pub fn with_collapse_below_m(mut self, collapse_below_m: f64) -> Self {
        self.spec.collapse_below_m = Some(collapse_below_m);
        self
    }

    /// Attach positioned petro samples (TVD, φ, Sw) so every realization
    /// populates from upscaled logs (draw priors fill uncovered cells).
    #[must_use]
    pub fn with_logs(mut self, samples: Vec<PetroSample>) -> Self {
        self.logs = if samples.is_empty() {
            None
        } else {
            // Sort by TVD once (reused across every realization) so population
            // binary-searches each cell's depth range instead of scanning (V2).
            Some(crate::population::sort_by_tvd(samples))
        };
        self
    }

    /// Attach an areal trend multiplier field applied to every realization's NTG
    /// (and φ if flagged) — the template analog of
    /// [`crate::StaticModelBuilder::with_areal_trend`]. Lateral shape only; the
    /// per-draw prior gives the level. See [`TrendSurface`].
    ///
    /// **Superseded (deprecation-tracked)** — the same interim status as the builder
    /// method it mirrors: the fuller path is [`crate::PropertyPipeline`] +
    /// [`crate::Gaussian::with_trend`] (collocated cokriging). Retained for the
    /// simple post-population multiplier case (`task_petekstatic_organize`).
    #[must_use]
    pub fn with_areal_trend(mut self, trend: TrendSurface) -> Self {
        self.trend = Some(trend);
        self
    }

    /// Attach a per-property geostatistical pipeline run every realization, in the
    /// default [`McMode::LevelShift`] mode: the field is propagated **once** (first
    /// realization) and reused with only the draw's per-property level shift applied
    /// (`decision_mc_composition`) — the cheap, ~ms-class MC path.
    #[must_use]
    pub fn with_property(self, pipeline: PropertyPipeline) -> Self {
        self.with_property_mode(pipeline, McMode::LevelShift)
    }

    /// Attach a per-property pipeline with an explicit [`McMode`] — `Resimulate`
    /// redraws a fresh SGS pattern per realization (heterogeneity uncertainty, a
    /// simulation per draw); `LevelShift` propagates once and shifts the level.
    #[must_use]
    pub fn with_property_mode(mut self, pipeline: PropertyPipeline, mode: McMode) -> Self {
        self.properties.push(PropertyModel {
            pipeline,
            mode,
            cached: None,
        });
        self
    }

    /// Attach a **zone-scoped** property pipeline (the MC analog of
    /// [`crate::StaticModelBuilder::with_zone_property`]) in the default
    /// [`McMode::LevelShift`] mode: the pipeline is run restricted to `zone_name`'s
    /// `k`-range over the per-zone constant priors, propagated **once** and reused
    /// with only the draw's per-property level shift thereafter. Stack-aware
    /// templates only ([`StaticModelTemplate::from_horizon_stack`]); a zone name
    /// absent from the stack errors at `realize`.
    ///
    /// Without this, a zoned MC over a model whose zones were staged via zone-scoped
    /// pipes realized those zones from the ZONE PRIORS alone, ignoring the
    /// upscale+SGS cubes — so a zero-spread zoned MC did not reproduce the built
    /// model's `in_place_by_zone` (`question_zoned_mc_zone_pipe_parity`).
    #[must_use]
    pub fn with_zone_property(
        self,
        zone_name: impl Into<String>,
        pipeline: PropertyPipeline,
    ) -> Self {
        self.with_zone_property_mode(zone_name, pipeline, McMode::LevelShift)
    }

    /// [`StaticModelTemplate::with_zone_property`] with an explicit [`McMode`] —
    /// `Resimulate` redraws a fresh per-draw SGS pattern in the zone; `LevelShift`
    /// propagates once and shifts the zone's level.
    #[must_use]
    pub fn with_zone_property_mode(
        mut self,
        zone_name: impl Into<String>,
        pipeline: PropertyPipeline,
        mode: McMode,
    ) -> Self {
        self.zone_properties.push(ZonePropertyModel {
            zone: zone_name.into(),
            pipeline,
            mode,
            cached: None,
        });
        self
    }

    /// Label the provenance record with the input-bundle identity.
    #[must_use]
    pub fn with_inputs_ref(mut self, inputs_ref: impl Into<String>) -> Self {
        self.spec.inputs_ref = Some(inputs_ref.into());
        self
    }

    /// Install a whole declarative [`BuildSpec`] in one call — the spec analog of
    /// chaining the individual `with_*` sugar (each of which mutates the same
    /// internal spec the deterministic builder consumes). The extrapolation policy
    /// and any well ties are (re-)applied through the same paths the individual
    /// setters use, so values — and every determinism contract — are identical.
    ///
    /// # Errors
    /// [`StaticError`] as [`StaticModelTemplate::with_extrapolation`] /
    /// [`StaticModelTemplate::with_well_ties`].
    pub fn with_spec(mut self, spec: BuildSpec) -> Result<Self, StaticError> {
        let ties = spec.well_ties.clone();
        self.spec = spec;
        self.spec.well_ties = Vec::new(); // re-stored by with_well_ties as it applies
        let policy = self.spec.extrapolation;
        self = self.with_extrapolation(policy)?;
        if !ties.is_empty() {
            self = self.with_well_ties(ties)?;
        }
        Ok(self)
    }

    /// The effective provenance label: the spec's, or the entry point's default.
    fn inputs_ref_string(&self) -> String {
        self.spec
            .inputs_ref
            .clone()
            .unwrap_or_else(|| self.default_inputs_ref.to_string())
    }

    /// Register the world georeference stamped onto every realized model's view
    /// frames — the MC analog of [`crate::StaticModelBuilder::with_georef`]. Column
    /// `(0, 0)`'s world centroid + world column spacing; non-finite / non-positive
    /// spacing is ignored (the local degenerate frame). Grid geometry is untouched.
    /// See [`Georef`].
    #[must_use]
    pub fn with_georef(
        mut self,
        origin_x: f64,
        origin_y: f64,
        spacing_x: f64,
        spacing_y: f64,
    ) -> Self {
        self.spec.georef = Georef::new(origin_x, origin_y, spacing_x, spacing_y);
        self
    }

    /// Register the world areal boundary ring stamped onto every realized model's map
    /// outline — the MC analog of [`crate::StaticModelBuilder::with_boundary`] (finding
    /// 2). Without it the stack template still emits a world-extent rectangle from the
    /// georef; pass the real closure ring for the true outline shape.
    #[must_use]
    pub fn with_boundary(mut self, ring: Vec<[f64; 2]>) -> Self {
        self.spec.boundary = if ring.is_empty() { None } else { Some(ring) };
        self
    }

    /// Set how well ties are applied ([`crate::TieSettings`]: datum substitution —
    /// the default — or bounded-radius locality). Thin sugar over the spec.
    /// **Set this BEFORE [`StaticModelTemplate::with_well_ties`]** — the template
    /// applies ties at that call (they are draw-invariant), so settings changed
    /// afterwards do not re-tie.
    #[must_use]
    pub fn with_tie_settings(mut self, ties: crate::spec::TieSettings) -> Self {
        self.spec.ties = ties;
        self
    }

    /// Tie the stack surfaces to well markers — the MC analog of
    /// [`crate::StaticModelBuilder::with_well_ties`] (`task_petekstatic_template_ties`).
    /// Ties are **draw-invariant** (they replace map controls with measured tops and
    /// re-solve the surface, which does not depend on the per-draw scalars), so they
    /// are applied **once, at template construction** — exactly like
    /// [`StaticModelTemplate::from_horizon_stack`]'s once-only surface resolution.
    /// Every realization then inherits the tied geometry at **zero per-draw cost**,
    /// and the recorded per-horizon tie residuals ride onto each realization's
    /// provenance ([`StaticModel::well_tie_residuals`]).
    ///
    /// Only a **mapped** horizon is re-tied (a tops-only horizon is already
    /// pick-conditioned), mirroring the builder. Stack-aware templates only.
    ///
    /// # Errors
    /// [`StaticError`] if the template is not stack-aware, a tie node is off the
    /// control lattice, a tie names an unknown horizon, or a re-solve fails.
    pub fn with_well_ties(mut self, ties: Vec<WellTie>) -> Result<Self, StaticError> {
        if ties.is_empty() {
            return Ok(self);
        }
        let (nx, ny) = (self.ni + 1, self.nj + 1);
        let policy = self.spec.extrapolation;
        let settings = self.spec.ties;
        let stack = self.stack.as_mut().ok_or_else(|| {
            StaticError::InvalidInput(
                "with_well_ties requires a stack-aware template (from_horizon_stack)".into(),
            )
        })?;
        // Record the pre-tie residuals against the untied surfaces, then fold each
        // mapped tie top into the source stack's gridded datums per the spec's
        // [`crate::TieSettings`] — byte-for-byte the builder's tie math
        // ([`substitute_tie_datums`], the single tie authority): the tie flows
        // through the same stack resolution as the map data.
        let (records, substituted) = substitute_tie_datums(
            &mut stack.source,
            &ties,
            &stack.surfaces,
            settings,
            stack.tie_dxdy,
            nx,
            ny,
        )?;
        if substituted {
            stack.surfaces = resolve_stack_surfaces(&stack.source, nx, ny, policy)?;
            // Re-run the per-interface order-repair over the tied surfaces (a
            // safety net — the tied resolution is ordered by construction),
            // honouring mapped-over-derived precedence, then rebuild the framework.
            let is_derived: Vec<bool> = stack
                .source
                .horizons
                .iter()
                .map(|h| matches!(h.source, HorizonSource::TopsOnly(_)))
                .collect();
            let mut interface_repairs = Vec::new();
            for i in 0..stack.surfaces.len() - 1 {
                let (upper, lower) = repair_interface(
                    &stack.surfaces[i],
                    &stack.surfaces[i + 1],
                    None,
                    true,
                    is_derived[i],
                    is_derived[i + 1],
                    i,
                    &mut interface_repairs,
                )?;
                stack.surfaces[i] = upper;
                stack.surfaces[i + 1] = lower;
            }
            stack.interface_repairs = interface_repairs;
            let all_contacts: Vec<Contact> =
                stack.zone_contacts.iter().flatten().copied().collect();
            stack.framework = stack_framework_from_names(
                &stack.horizon_names,
                &stack.surfaces,
                all_contacts,
                nx,
                ny,
                stack_boundary_ring(None, None, nx, ny),
            );
        }
        self.well_ties = records;
        self.spec.well_ties = ties;
        Ok(self)
    }

    /// Set the **extrapolation policy** (see
    /// [`crate::StaticModelBuilder::with_extrapolation`]) and — on a stack-aware
    /// template — re-resolve the stored surfaces under it, so the template's
    /// every realization inherits the policy. Default
    /// [`ExtrapolationPolicy::DecayToData`].
    ///
    /// # Errors
    /// [`StaticError`] if the re-resolution fails (stack templates only; a
    /// non-stack template just records the policy).
    pub fn with_extrapolation(mut self, policy: ExtrapolationPolicy) -> Result<Self, StaticError> {
        self.spec.extrapolation = policy;
        let (nx, ny) = (self.ni + 1, self.nj + 1);
        if let Some(stack) = self.stack.as_mut() {
            stack.surfaces = resolve_stack_surfaces(&stack.source, nx, ny, policy)?;
            let is_derived: Vec<bool> = stack
                .source
                .horizons
                .iter()
                .map(|h| matches!(h.source, HorizonSource::TopsOnly(_)))
                .collect();
            let mut interface_repairs = Vec::new();
            for i in 0..stack.surfaces.len() - 1 {
                let (upper, lower) = repair_interface(
                    &stack.surfaces[i],
                    &stack.surfaces[i + 1],
                    None,
                    true,
                    is_derived[i],
                    is_derived[i + 1],
                    i,
                    &mut interface_repairs,
                )?;
                stack.surfaces[i] = upper;
                stack.surfaces[i + 1] = lower;
            }
            stack.interface_repairs = interface_repairs;
            let all_contacts: Vec<Contact> =
                stack.zone_contacts.iter().flatten().copied().collect();
            stack.framework = stack_framework_from_names(
                &stack.horizon_names,
                &stack.surfaces,
                all_contacts,
                nx,
                ny,
                stack_boundary_ring(None, None, nx, ny),
            );
            self.framework = stack.framework.clone();
        }
        Ok(self)
    }

    /// Opt into **sugar-cube** section rendering (flat-box cells) — default `false`
    /// (dip-following trapezoids). The MC analog of
    /// [`crate::StaticModelBuilder::with_sugar_cube`]; flows into each realization's
    /// provenance + section bundle.
    #[must_use]
    pub fn with_sugar_cube(mut self, sugar_cube: bool) -> Self {
        self.spec.sugar_cube = sugar_cube;
        self
    }

    /// Realize one static model from a sampled draw — the cheap per-realization
    /// call. The allocating convenience: it builds an empty model and fills it via
    /// [`StaticModelTemplate::realize_into`] (one code path), so it warm-starts the
    /// surface solve from the previous realization's converged surface and updates
    /// the chain (`&mut self`) identically.
    ///
    /// On the MC hot path prefer [`StaticModelTemplate::realize_into`] with a reused
    /// model — it recycles the ~100 MB/draw of ZCORN + cube buffers in place.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if a draw scalar is non-physical (H2);
    /// [`StaticError::Grid`] if a structural shift is off-lattice or the
    /// warm-started solve/layering fails.
    pub fn realize(&mut self, draw: &RealizationDraw) -> Result<StaticModel, StaticError> {
        let mut model = StaticModel::empty(self.framework.clone());
        self.realize_into(draw, &mut model)?;
        Ok(model)
    }

    /// A fresh, empty [`StaticModel`] to drive a [`StaticModelTemplate::realize_into`]
    /// loop against — one per MC worker. It starts as a trivial single-cell model;
    /// the first `realize_into` grows its buffers and every subsequent draw recycles
    /// them. Cheap (a unit box), so allocate one per shard, not one per draw.
    #[must_use]
    pub fn reusable_model(&self) -> StaticModel {
        StaticModel::empty(self.framework.clone())
    }

    /// Realize into a **reused** [`StaticModel`], recycling its geometry (ZCORN +
    /// COORD) and property-cube allocations in place — the allocation-bound MC hot
    /// path (`task_petekstatic_realize_into`). Every geometry and cube buffer of
    /// `model` is fully overwritten (or grown, on the first draw), so the result is
    /// **bit-identical** to [`StaticModelTemplate::realize`] on the same draw and
    /// carries no stale data across draws. The warm-start chain advances identically,
    /// so determinism is unchanged.
    ///
    /// One reusable model per MC worker + the template's own layering scratch (both
    /// per-worker under one-template-per-shard) drop the steady-state per-draw
    /// allocation from ~100 MB to a few surface buffers.
    ///
    /// # Errors
    /// Same as [`StaticModelTemplate::realize`].
    pub fn realize_into(
        &mut self,
        draw: &RealizationDraw,
        model: &mut StaticModel,
    ) -> Result<(), StaticError> {
        if self.stack.is_some() {
            return self.realize_into_stack(draw, model);
        }
        // H2 at the seam: reject a garbage draw before any geometry work.
        validate_positive("area_m2", draw.area_m2)?;
        validate_positive("gross_height_m", draw.gross_height_m)?;
        validate_fraction("porosity", draw.porosity)?;
        validate_fraction("net_to_gross", draw.net_to_gross)?;
        validate_fraction("water_saturation", draw.water_saturation)?;
        if !draw.contact_depth_m.is_finite() {
            return Err(StaticError::InvalidInput(format!(
                "contact depth must be finite, got {}",
                draw.contact_depth_m
            )));
        }
        if let Some(goc) = draw.goc_depth_m {
            if !goc.is_finite() {
                return Err(StaticError::InvalidInput(format!(
                    "GOC depth must be finite, got {goc}"
                )));
            }
            if goc > draw.contact_depth_m {
                return Err(StaticError::InvalidInput(format!(
                    "GOC ({goc}) must be shallower than the OWC ({})",
                    draw.contact_depth_m
                )));
            }
        }

        // Controls this realization: nominal topology + optional structural shifts.
        // The no-shift branch (the common case) **borrows** the template's control
        // set — the solve only reads it — so no per-realize deep clone; a
        // structural draw takes an owned, shifted copy (copy-on-write).
        let controls: Cow<'_, [Control]> = match &draw.structural {
            None => Cow::Borrowed(&self.top_controls),
            Some(p) => {
                let mut c = self.top_controls.clone();
                for &(ip, jp, dz) in &p.control_shifts {
                    match c.iter_mut().find(|k| k.ip == ip && k.jp == jp) {
                        Some(k) => k.z += dz,
                        None => {
                            return Err(StaticError::Grid(format!(
                                "structural shift at ({ip},{jp}) has no matching control"
                            )))
                        }
                    }
                }
                Cow::Owned(c)
            }
        };

        // Warm-started solve; advance the chain to the just-converged surface. The
        // warm-start chain always advances on the UNPERTURBED converged surface — a
        // structural perturbation is per-draw noise, never structural drift the chain
        // should accumulate.
        let top = solve_surface_seeded(&self.seed, &controls)?;
        self.seed = top.clone();
        let (dx, dy) = spacing(draw.area_m2, self.ni, self.nj);

        // TOP surface for THIS draw: add its correlated depth perturbation field
        // (`decision_structural_uncertainty_isochore`) when the draw carries one —
        // otherwise the converged surface directly. The 2-surface path has no well
        // ties (they are stack-only), so no tie-node pinning here.
        let top_surface: Surface = match &draw.top_structural {
            Some(pf) => {
                let (nx, ny) = (self.ni + 1, self.nj + 1);
                let (ox, oy) = self
                    .spec
                    .georef
                    .map_or((0.0, 0.0), |g| (g.origin_x, g.origin_y));
                let fld = perturbation_field(
                    pf,
                    ox,
                    oy,
                    dx,
                    dy,
                    nx,
                    ny,
                    horizon_seed(draw.seed_index, 0),
                )?;
                top.surface().offset_by_field(&fld)?
            }
            None => top.surface().clone(),
        };

        // Column base: with a Base horizon the template's gross field carries the
        // real relief SHAPE and the draw's `gross_height_m` sets the LEVEL
        // (per-node gross = g × gross_height_m / mean(g)); without one, the
        // constant offset (`decision_template_gross_scaling`).
        let base = match &self.gross_field {
            Some(gf) => {
                let scale = draw.gross_height_m / gf.mean;
                let dz: Vec<f64> = gf.node_dz.iter().map(|g| g * scale).collect();
                top_surface.offset_by_field(&dz)?
            }
            None => top_surface.offset_by(draw.gross_height_m),
        };
        // R1/R-c: guard the scaled base against crossing above the top
        // (thin/crossing columns collapse GRV). Opt-in `min_thickness_m` repairs it
        // (pull the base to a minimum thickness below the top, record a warning);
        // otherwise error by default, or clamp to zero gross if opted in.
        let mut warnings = Vec::new();
        let base = match self.spec.min_thickness_m {
            Some(min_t) => {
                let (repaired, columns, worst_m) =
                    base.repair_min_thickness(&top_surface, min_t)?;
                if columns > 0 {
                    warnings.push(BuildWarning::ThinColumnsRepaired { columns, worst_m });
                }
                repaired
            }
            None => base.guard_below(&top_surface, self.spec.clamp_base_to_top)?,
        };
        // Layering: under a Follow conformity `nk` is dz-derived and thin columns
        // truncate; read the effective nk + report back (`self.nk` only honoured by
        // `Proportional`). The derived nk is geometry-driven, so it is stable across
        // realizations that share a conformity + gross field (MC reproducibility).
        //
        // Recycle the model's geometry buffers: take the (possibly empty, on the
        // first draw) ZCORN + COORD Vecs, refill them through the layering core using
        // the template's reusable scratch, then install them back under the new dims —
        // no fresh 64 MB ZCORN / 1.9 MB COORD alloc on the steady-state path.
        let specs = [ZoneLayerSpec {
            conformity: self.conformity,
            requested_nk: self.nk,
        }];
        let (mut coord, mut zcorn) = model.grid_mut().take_geometry_buffers();
        let layered = layer_grid_stack_into(
            &[&top_surface, &base],
            dx,
            dy,
            &specs,
            self.spec.collapse_below_m,
            &mut self.scratch,
            &mut coord,
            &mut zcorn,
        )?;
        let nk = layered.nk;
        if layered.truncated_cells > 0 {
            warnings.push(BuildWarning::LayersTruncated {
                cells: layered.truncated_cells,
            });
        }
        if layered.collapsed_cells > 0 {
            warnings.push(BuildWarning::CellsCollapsed {
                cells: layered.collapsed_cells,
            });
        }
        if layered.nk_capped {
            warnings.push(BuildWarning::LayerCountCapped { nk });
        }
        model
            .grid_mut()
            .install_geometry(layered.dims, coord, zcorn);
        let grid = model.grid_mut();
        let priors = ConstantPriors {
            porosity: draw.porosity,
            net_to_gross: draw.net_to_gross,
            water_saturation: draw.water_saturation,
        };
        populate(grid, priors, self.logs.as_deref(), self.trend.as_ref())?;

        // P5 per-property geostatistical pipelines, per MC mode
        // (`decision_mc_composition`).
        let property_reports =
            realize_properties(&mut self.properties, grid, draw, self.spec.georef)?;

        // Framework snapshot for this realization: the template's framework with
        // the drawn contact depth(s) (the draw owns the column base). Cloning the
        // `Wireframe` is cheap — its `horizons` are `Arc`-shared (an O(1) refcount
        // bump, not a per-node depth-cube copy); only the two contacts are rebuilt
        // below. A drawn GOC adds a gas cap above the (lower) OWC, making it a
        // two-contact column that `StaticModel::in_place` splits gas/oil.
        let mut framework = self.framework.clone();
        let mut contacts = Vec::with_capacity(2);
        if let Some(goc) = draw.goc_depth_m {
            contacts.push(Contact {
                kind: ContactKind::Goc,
                depth_m: goc,
                hardness: Hardness::Assumed,
            });
        }
        contacts.push(Contact {
            kind: ContactKind::Owc,
            depth_m: draw.contact_depth_m,
            hardness: Hardness::Assumed,
        });
        framework.contacts = contacts;

        let population = if self.logs.is_some() || !self.properties.is_empty() {
            PopulationMode::Logs
        } else {
            PopulationMode::Priors
        };
        let provenance = Provenance {
            inputs_ref: self.inputs_ref_string(),
            solve_opts: self.solve_opts,
            conformity: self.conformity,
            nk,
            population,
            warnings,
            realization: Some(draw.clone()),
            property_reports,
            stack: None,
            well_ties: Vec::new(),
            sugar_cube: self.spec.sugar_cube,
        };
        // Overwrite the reused model's non-geometry state (the grid was recycled in
        // place above). Bit-identical to the fields `StaticModel::new(..)
        // .with_georef_opt(..)` would set for this draw.
        model.reset_state(
            framework,
            ZoneTable::single(nk),
            provenance,
            draw.sw_gas,
            self.spec.georef,
            None,
        );
        Ok(())
    }

    /// The stack-aware realize (`from_horizon_stack`): recycle `model`'s geometry +
    /// cube buffers exactly as the 2-surface path, but over the fixed multi-horizon
    /// framework. Only the areal footprint (spacing), the per-zone contacts, and the
    /// per-zone property levels vary per draw — the surfaces are template-fixed, so
    /// the realization is bit-deterministic and carries no per-draw structural chain.
    fn realize_into_stack(
        &mut self,
        draw: &RealizationDraw,
        model: &mut StaticModel,
    ) -> Result<(), StaticError> {
        // H2 at the seam (base priors); per-zone draws are validated in the inner fn.
        validate_positive("area_m2", draw.area_m2)?;
        validate_fraction("porosity", draw.porosity)?;
        validate_fraction("net_to_gross", draw.net_to_gross)?;
        validate_fraction("water_saturation", draw.water_saturation)?;
        // Take the fixed stack state out so `self.scratch` / `self.properties` are
        // freely borrowable in the body; restore it before returning.
        let stack = self
            .stack
            .take()
            .expect("realize_into_stack requires a stack");
        let result = self.realize_into_stack_inner(&stack, draw, model);
        self.stack = Some(stack);
        result
    }

    fn realize_into_stack_inner(
        &mut self,
        stack: &StackState,
        draw: &RealizationDraw,
        model: &mut StaticModel,
    ) -> Result<(), StaticError> {
        let nzones = stack.zone_names.len();
        // Validate per-zone contact + property draws (H2).
        for zd in &draw.zones {
            if zd.zone >= nzones {
                return Err(StaticError::InvalidInput(format!(
                    "zone draw index {} out of range (nzones={nzones})",
                    zd.zone
                )));
            }
            if let Some(owc) = zd.owc_depth_m {
                if !owc.is_finite() {
                    return Err(StaticError::InvalidInput(format!(
                        "zone {} OWC depth must be finite, got {owc}",
                        zd.zone
                    )));
                }
            }
            if let Some(goc) = zd.goc_depth_m {
                if !goc.is_finite() {
                    return Err(StaticError::InvalidInput(format!(
                        "zone {} GOC depth must be finite, got {goc}",
                        zd.zone
                    )));
                }
                if let Some(owc) = zd.owc_depth_m {
                    if goc > owc {
                        return Err(StaticError::InvalidInput(format!(
                            "zone {} GOC ({goc}) must be shallower than its OWC ({owc})",
                            zd.zone
                        )));
                    }
                }
            }
            if let Some(p) = zd.porosity {
                validate_fraction("porosity", p)?;
            }
            if let Some(n) = zd.net_to_gross {
                validate_fraction("net_to_gross", n)?;
            }
            if let Some(s) = zd.water_saturation {
                validate_fraction("water_saturation", s)?;
            }
        }

        // Layer the surfaces at this draw's spacing, recycling the model's ZCORN +
        // COORD buffers in place (the allocation-bound MC hot path). When the draw
        // carries structural perturbation (a top depth field and/or per-zone isochore
        // fields), reseat the surfaces in isochore space first; otherwise the fixed
        // template surfaces are used directly (zero extra cost).
        let (dx, dy) = spacing(draw.area_m2, self.ni, self.nj);
        let structural_active = draw.top_structural.is_some()
            || draw.zones.iter().any(|z| z.isochore_structural.is_some());
        let perturbed_surfaces: Vec<Surface>;
        let surf_refs: Vec<&Surface> = if structural_active {
            perturbed_surfaces = self.perturb_stack_surfaces(stack, draw, dx, dy)?;
            perturbed_surfaces.iter().collect()
        } else {
            stack.surfaces.iter().collect()
        };
        let (mut coord, mut zcorn) = model.grid_mut().take_geometry_buffers();
        let mut warnings = Vec::new();
        let layered = layer_grid_stack_into(
            &surf_refs,
            dx,
            dy,
            &stack.zone_specs,
            self.spec.collapse_below_m,
            &mut self.scratch,
            &mut coord,
            &mut zcorn,
        )?;
        let nk = layered.nk;
        if layered.truncated_cells > 0 {
            warnings.push(BuildWarning::LayersTruncated {
                cells: layered.truncated_cells,
            });
        }
        if layered.collapsed_cells > 0 {
            warnings.push(BuildWarning::CellsCollapsed {
                cells: layered.collapsed_cells,
            });
        }
        if layered.nk_capped {
            warnings.push(BuildWarning::LayerCountCapped { nk });
        }
        model
            .grid_mut()
            .install_geometry(layered.dims, coord, zcorn);

        // Populate base priors, whole-model pipelines, then per-zone property levels.
        let grid = model.grid_mut();
        let priors = ConstantPriors {
            porosity: draw.porosity,
            net_to_gross: draw.net_to_gross,
            water_saturation: draw.water_saturation,
        };
        populate(grid, priors, self.logs.as_deref(), self.trend.as_ref())?;
        let mut property_reports =
            realize_properties(&mut self.properties, grid, draw, self.spec.georef)?;
        for zd in &draw.zones {
            if zd.porosity.is_some() || zd.net_to_gross.is_some() || zd.water_saturation.is_some() {
                let zp = ConstantPriors {
                    porosity: zd.porosity.unwrap_or(draw.porosity),
                    net_to_gross: zd.net_to_gross.unwrap_or(draw.net_to_gross),
                    water_saturation: zd.water_saturation.unwrap_or(draw.water_saturation),
                };
                override_zone_priors(grid, zp, layered.zones[zd.zone].k_range())?;
            }
        }
        // Per-zone geostatistical pipelines (`with_zone_property`), each restricted to
        // its zone's k-range OVER the per-zone constant priors just written — the MC
        // analog of the builder's per-zone population. Level shifts apply on top of the
        // staged cube (as the whole-model path does), so a zero-spread zoned MC
        // reproduces the built model's in_place_by_zone on every piped zone
        // (question_zoned_mc_zone_pipe_parity).
        if !self.zone_properties.is_empty() {
            let zone_kranges: Vec<core::ops::Range<usize>> =
                layered.zones.iter().map(|z| z.k_range()).collect();
            realize_zone_properties(
                &mut self.zone_properties,
                grid,
                draw,
                self.spec.georef,
                &stack.zone_names,
                &zone_kranges,
                &mut property_reports,
            )?;
        }

        // Per-zone contacts: a draw's ZoneDraw REPLACES the static contacts for that
        // zone (a ZoneDraw with neither GOC nor OWC = an explicitly contactless zone);
        // an unmentioned zone keeps the template's static contacts.
        let mk = |kind: ContactKind, depth_m: f64| Contact {
            kind,
            depth_m,
            hardness: Hardness::Assumed,
        };
        let zone_contacts: Vec<Vec<Contact>> = (0..nzones)
            .map(|z| match draw.zones.iter().find(|zd| zd.zone == z) {
                Some(zd) => {
                    let mut c = Vec::new();
                    if let Some(goc) = zd.goc_depth_m {
                        c.push(mk(ContactKind::Goc, goc));
                    }
                    if let Some(owc) = zd.owc_depth_m {
                        c.push(mk(ContactKind::Owc, owc));
                    }
                    c
                }
                None => stack.zone_contacts[z].clone(),
            })
            .collect();
        let all_contacts: Vec<Contact> = zone_contacts.iter().flatten().copied().collect();
        let mut framework = stack.framework.clone();
        framework.contacts = all_contacts;
        // Resolve the world outline now that the georef / boundary are known (the
        // template built its framework with a placeholder unit square; finding 2).
        framework.boundary.ring = stack_boundary_ring(
            self.spec.boundary.as_ref(),
            self.spec.georef,
            self.ni + 1,
            self.nj + 1,
        );

        // Real zone table + per-zone provenance for this realization.
        let per_zone_nk: Vec<usize> = layered.zones.iter().map(|z| z.nk).collect();
        let zones = ZoneTable::from_stack(
            &stack.horizon_names,
            &stack.zone_names,
            &stack.zone_colors,
            &per_zone_nk,
        );
        let zone_prov: Vec<ZoneProvenance> = layered
            .zones
            .iter()
            .enumerate()
            .map(|(z, sz)| ZoneProvenance {
                name: stack.zone_names[z].clone(),
                top_horizon: stack.horizon_names[z].clone(),
                base_horizon: stack.horizon_names[z + 1].clone(),
                conformity: sz.conformity,
                nk: sz.nk,
                k_start: sz.k_start,
                truncated_cells: sz.truncated_cells,
            })
            .collect();
        let stack_prov = StackProvenance {
            horizons: stack.horizon_names.clone(),
            zones: zone_prov,
            interface_repairs: stack.interface_repairs.clone(),
        };

        let population = if self.logs.is_some() || !self.properties.is_empty() {
            PopulationMode::Logs
        } else {
            PopulationMode::Priors
        };
        let provenance = Provenance {
            inputs_ref: self.inputs_ref_string(),
            solve_opts: self.solve_opts,
            conformity: self.conformity,
            nk,
            population,
            warnings,
            realization: Some(draw.clone()),
            property_reports,
            stack: Some(stack_prov),
            well_ties: self.well_ties.clone(),
            sugar_cube: self.spec.sugar_cube,
        };
        model.reset_state(
            framework,
            zones,
            provenance,
            draw.sw_gas,
            self.spec.georef,
            Some(zone_contacts),
        );
        Ok(())
    }

    /// Build this draw's **perturbed** stack surfaces
    /// (`decision_structural_uncertainty_isochore`): the TOP surface plus its
    /// correlated depth field, then every deeper horizon reseated in ISOCHORE space
    /// — each zone's fixed base thickness plus its (optional) correlated thickness
    /// field, clamped `>= 0` and **zero-masked where the base isochore is exactly 0**
    /// — so a perturbed stack can never invert or resurrect a collapsed zone (the
    /// build-down construction holds by construction). Only called when a draw carries
    /// structural perturbation (the common no-perturbation path realizes the fixed
    /// surfaces at zero extra cost). Perturbation is pinned to zero at well-tie nodes.
    fn perturb_stack_surfaces(
        &self,
        stack: &StackState,
        draw: &RealizationDraw,
        dx: f64,
        dy: f64,
    ) -> Result<Vec<Surface>, StaticError> {
        let (nx, ny) = (self.ni + 1, self.nj + 1);
        let (ox, oy) = self
            .spec
            .georef
            .map_or((0.0, 0.0), |g| (g.origin_x, g.origin_y));
        let nsurf = stack.surfaces.len();
        let mut out: Vec<Surface> = Vec::with_capacity(nsurf);

        // TOP surface (horizon 0): add its correlated depth field.
        let top_surf = match &draw.top_structural {
            Some(pf) => {
                let mut fld = perturbation_field(
                    pf,
                    ox,
                    oy,
                    dx,
                    dy,
                    nx,
                    ny,
                    horizon_seed(draw.seed_index, 0),
                )?;
                pin_ties(&mut fld, &self.well_ties, &stack.horizon_names[0], nx, ny);
                stack.surfaces[0].offset_by_field(&fld)?
            }
            None => stack.surfaces[0].clone(),
        };
        out.push(top_surf);

        // Deeper horizons: perturb each zone's isochore, reseat below the perturbed
        // horizon above. `t >= 0` always → ordering can never invert; a merged
        // (base-isochore == 0) node stays merged in every draw.
        for z in 0..nsurf - 1 {
            let above = &stack.surfaces[z];
            let below = &stack.surfaces[z + 1];
            let mut t = vec![0.0; nx * ny];
            for jp in 0..ny {
                for ip in 0..nx {
                    t[jp * nx + ip] = (below.z(ip, jp) - above.z(ip, jp)).max(0.0);
                }
            }
            if let Some(pf) = draw
                .zones
                .iter()
                .find(|zd| zd.zone == z)
                .and_then(|zd| zd.isochore_structural.as_ref())
            {
                let mut p = perturbation_field(
                    pf,
                    ox,
                    oy,
                    dx,
                    dy,
                    nx,
                    ny,
                    horizon_seed(draw.seed_index, (z + 1) as u64),
                )?;
                pin_ties(&mut p, &self.well_ties, &stack.horizon_names[z + 1], nx, ny);
                for idx in 0..nx * ny {
                    // Zero-masked at an exact merge (base isochore == 0); else clamp
                    // the perturbed thickness to non-negative.
                    t[idx] = if t[idx] == 0.0 {
                        0.0
                    } else {
                        (t[idx] + p[idx]).max(0.0)
                    };
                }
            }
            let prev = &out[z];
            out.push(prev.offset_by_field(&t)?);
        }
        Ok(out)
    }
}

/// Run the attached property pipelines against this realization's `grid` per their
/// MC mode (`decision_mc_composition`), returning the per-property reports.
///
/// - [`McMode::Resimulate`]: re-run the pipeline with a per-draw reseed — a fresh,
///   reproducible pattern each draw.
/// - [`McMode::LevelShift`]: propagate once (cached at the first realization), then
///   set the cube to that pattern plus the draw's per-property additive level shift
///   — same pattern, moved level.
fn realize_properties(
    properties: &mut [PropertyModel],
    grid: &mut Grid,
    draw: &RealizationDraw,
    georef: Option<Georef>,
) -> Result<Vec<PropertyReport>, StaticError> {
    let mut reports = Vec::with_capacity(properties.len());
    for pm in properties.iter_mut() {
        let name = pm.pipeline.name().to_string();
        let report = match pm.mode {
            McMode::Resimulate => pm
                .pipeline
                .reseeded(draw.seed_index)
                .apply_with_georef(grid, georef)?,
            McMode::LevelShift => {
                // Propagate once, then cache the pattern + report; reuse thereafter.
                if pm.cached.is_none() {
                    let r = pm.pipeline.apply_with_georef(grid, georef)?;
                    let pattern = grid
                        .properties()
                        .get(&name)
                        .expect("pipeline set the property cube")
                        .values
                        .clone();
                    pm.cached = Some((pattern, r));
                }
                let (pattern, report) = pm.cached.as_ref().expect("cached above");
                let shift = draw.property_shift(&name);
                // A fraction cube (PORO/NTG/SW) legitimately holds boundary cells
                // (a non-net cell NTG=0; an aquifer cell SW=1). A level shift adds
                // `shift` to every cell, so a boundary cell would escape [0,1] and
                // the per-cell H2 range check (volumetrics `validate_fraction`)
                // would reject the whole draw — making property uncertainty on any
                // log-conditioned real model impossible. Shift **then clamp per
                // cell** so boundary cells *saturate*: an SW=1 aquifer cell stays 1
                // under a positive shift (physically right — it cannot get wetter),
                // and slides into the interior under a negative one. We saturate
                // rather than *skip* saturated cells so a negative shift still moves
                // an SW=1 cell — skipping would freeze it wrongly. NaN (undefined)
                // cells pass through unchanged (`clamp` returns NaN). Non-fraction
                // cubes (e.g. permeability) are shifted without a clamp. H2 still
                // guards the DRAWN INPUTS themselves (`realize` above) — only the
                // per-cell application saturates.
                let is_fraction = matches!(name.as_str(), PORO | NTG | SW);
                // Recycle the cube's value buffer in place (take → refill → reinstall)
                // so the steady-state LevelShift MC path allocates no fresh cube per
                // draw; the shifted pattern fully overwrites every cell.
                let mut values = grid.properties_mut().take_values(&name);
                values.clear();
                values.reserve(pattern.len());
                values.extend(pattern.iter().map(|v| {
                    let shifted = v + shift;
                    if is_fraction {
                        shifted.clamp(0.0, 1.0)
                    } else {
                        shifted
                    }
                }));
                grid.properties_mut().set(Property {
                    name: name.clone(),
                    values,
                })?;
                report.clone()
            }
        };
        reports.push(report);
    }
    Ok(reports)
}

/// Run the attached **zone-scoped** pipelines against this realization's `grid`,
/// each restricted to its named zone's `k`-range (resolved via `zone_names` →
/// `zone_kranges`), appending their reports. The per-zone analog of
/// [`realize_properties`]: it merges each pipeline's cube into the zone's slice
/// over the per-zone constant priors already written, with the draw's per-property
/// level shift applied **only to that zone's cells** on top (so a zero-spread draw
/// reproduces the deterministic per-zone build).
///
/// - [`McMode::Resimulate`]: re-run the pipeline in-zone with a per-draw reseed.
/// - [`McMode::LevelShift`]: propagate once (cached at the first realization), then
///   overwrite only the zone's `k`-range cells with the cached pattern + shift,
///   preserving this draw's out-of-zone cells (base priors / other zones).
fn realize_zone_properties(
    zone_properties: &mut [ZonePropertyModel],
    grid: &mut Grid,
    draw: &RealizationDraw,
    georef: Option<Georef>,
    zone_names: &[String],
    zone_kranges: &[core::ops::Range<usize>],
    reports: &mut Vec<PropertyReport>,
) -> Result<(), StaticError> {
    for zm in zone_properties.iter_mut() {
        let name = zm.pipeline.name().to_string();
        let k_range = zone_names
            .iter()
            .position(|n| n == &zm.zone)
            .map(|i| zone_kranges[i].clone())
            .ok_or_else(|| {
                StaticError::InvalidInput(format!(
                    "zone property '{name}': zone '{}' is not in the stack's zone table",
                    zm.zone
                ))
            })?;
        let zone = zm.zone.clone();
        let report = match zm.mode {
            McMode::Resimulate => zm
                .pipeline
                .reseeded(draw.seed_index)
                .apply_in_zone_with_georef(grid, k_range, georef)
                .map_err(|e| StaticError::InvalidInput(format!("zone '{zone}': {e}")))?,
            McMode::LevelShift => {
                // Propagate once (in-zone, over the current per-zone priors), cache the
                // full cube pattern + report; reuse thereafter. Fall through to apply
                // the draw's level shift on top even on the first realization.
                if zm.cached.is_none() {
                    let r = zm
                        .pipeline
                        .apply_in_zone_with_georef(grid, k_range.clone(), georef)
                        .map_err(|e| StaticError::InvalidInput(format!("zone '{zone}': {e}")))?;
                    let pattern = grid
                        .properties()
                        .get(&name)
                        .expect("zone pipeline set the property cube")
                        .values
                        .clone();
                    zm.cached = Some((pattern, r));
                }
                let (pattern, report) = zm.cached.as_ref().expect("cached above");
                let shift = draw.property_shift(&name);
                let is_fraction = matches!(name.as_str(), PORO | NTG | SW);
                let dims = grid.dims();
                let (ni, nj, nk) = (dims.ni, dims.nj, dims.nk);
                // Start from this draw's current cube (base priors + zone constant
                // priors, or an all-NaN cube if the property is zone-only) so the
                // out-of-zone cells stay THIS draw's values; overwrite only the zone
                // slice with the cached pattern + level shift.
                let mut values = grid.properties_mut().take_values(&name);
                if values.len() != ni * nj * nk {
                    values = vec![f64::NAN; ni * nj * nk];
                }
                for k in k_range.clone() {
                    for j in 0..nj {
                        for i in 0..ni {
                            let idx = (k * nj + j) * ni + i;
                            let shifted = pattern[idx] + shift;
                            values[idx] = if is_fraction {
                                shifted.clamp(0.0, 1.0)
                            } else {
                                shifted
                            };
                        }
                    }
                }
                grid.properties_mut().set(Property {
                    name: name.clone(),
                    values,
                })?;
                report.clone()
            }
        };
        reports.push(report);
    }
    Ok(())
}

#[cfg(test)]
mod stack_tests {
    //! Stack-aware MC template (`from_horizon_stack`, P8
    //! `task_petekstatic_multizone_2`): bit-determinism, the `realize_into`
    //! stale-buffer bit-compare, and per-zone contact/property draws.
    use super::*;
    use crate::draw::ZoneDraw;
    use crate::HorizonStack;
    use srs_gridder::Conformity;
    use srs_wireframe::{GriddedDepth, HorizonRole};

    const N: usize = 9; // 8×8 cells

    fn flat_surf(depth: f64) -> GriddedDepth {
        GriddedDepth {
            ncol: N,
            nrow: N,
            depth_m: vec![depth; N * N],
            is_control: vec![true; N * N],
        }
    }

    /// A 3-horizon / 2-zone flat stack: Top 5000, Mid 5030, Base 5060 → zones
    /// Z0 (5000–5030) and Z1 (5030–5060), both Proportional, contactless by default.
    fn two_zone_stack() -> HorizonStack {
        use crate::builder::{HorizonSource, StackHorizon, StackZone};
        HorizonStack {
            horizons: vec![
                StackHorizon {
                    name: "H0".into(),
                    source: HorizonSource::Mapped(flat_surf(5000.0)),
                },
                StackHorizon {
                    name: "H1".into(),
                    source: HorizonSource::Mapped(flat_surf(5030.0)),
                },
                StackHorizon {
                    name: "H2".into(),
                    source: HorizonSource::Mapped(flat_surf(5060.0)),
                },
            ],
            zone_layers: vec![
                StackZone::new("Z0", Conformity::Proportional, 5, Vec::new()),
                StackZone::new("Z1", Conformity::Proportional, 4, Vec::new()),
            ],
        }
    }

    fn opts() -> BuildOpts {
        BuildOpts {
            area_m2: 1_000_000.0,
            gross_height_m: 0.0, // unused on the stack path
            nk: 0,               // unused
            conformity: Conformity::Proportional,
            solve_opts: SolveOpts::default(),
            priors: ConstantPriors {
                porosity: 0.2,
                net_to_gross: 0.8,
                water_saturation: 0.3,
            },
        }
    }

    /// The full realized state that must be bit-stable across fresh vs stale-buffer
    /// realizes (geometry volume + every cube + per-zone in-place).
    #[derive(PartialEq, Debug)]
    struct Fingerprint {
        bulk_volume: f64,
        poro: Vec<f64>,
        ntg: Vec<f64>,
        sw: Vec<f64>,
        zone_grv: Vec<f64>,
        zone_hcpv: Vec<f64>,
    }

    fn fingerprint(m: &StaticModel) -> Fingerprint {
        let zoned = m.in_place_by_zone().unwrap();
        Fingerprint {
            bulk_volume: m.grid().bulk_volume(),
            poro: m.property(PORO).unwrap().values.clone(),
            ntg: m.property(NTG).unwrap().values.clone(),
            sw: m.property(SW).unwrap().values.clone(),
            zone_grv: zoned.zones.iter().map(|z| z.in_place.grv_m3).collect(),
            zone_hcpv: zoned.zones.iter().map(|z| z.in_place.hcpv_m3).collect(),
        }
    }

    /// A draw with a per-zone oil contact on Z0 and a per-zone porosity level on Z1.
    fn draw_a() -> RealizationDraw {
        RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 1)
            .with_zone_draw(ZoneDraw::new(0).with_owc(5015.0))
            .with_zone_draw(ZoneDraw::new(1).with_priors(0.1, 0.5, 0.4))
    }

    /// A DIFFERENT draw (different area + different per-zone contacts/levels).
    fn draw_b() -> RealizationDraw {
        RealizationDraw::new(900_000.0, 0.0, 0.0, 0.24, 0.85, 0.28, 2)
            .with_zone_draw(ZoneDraw::new(0).with_goc(5010.0).with_owc(5025.0))
            .with_zone_draw(ZoneDraw::new(1).with_owc(5050.0).with_priors(0.3, 0.9, 0.2))
    }

    #[test]
    fn stack_template_realize_is_deterministic_across_fresh_templates() {
        let draw = draw_a();
        let mut t1 = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts()).unwrap();
        let mut t2 = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts()).unwrap();
        let a = t1.realize(&draw).unwrap();
        let b = t2.realize(&draw).unwrap();
        assert_eq!(
            fingerprint(&a),
            fingerprint(&b),
            "stack realize not deterministic"
        );
        // The stack framework is real multi-zone: 3 horizons, 2 zones.
        assert_eq!(a.framework().horizons.len(), 3);
        assert_eq!(a.zones().zones().len(), 2);
        assert_eq!(a.framework().horizons[1].role, HorizonRole::Intermediate);
    }

    #[test]
    fn stack_realize_into_stale_buffer_bit_matches_fresh() {
        // The required stale-buffer pattern: two DIFFERENT stacked draws into ONE
        // reused model must leave the model bit-identical to a fresh realize of the
        // second draw (no stale geometry/cube data survives the recycle).
        let (a, b) = (draw_a(), draw_b());

        let fresh_b = {
            let mut t = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts()).unwrap();
            let mut m = t.reusable_model();
            t.realize_into(&b, &mut m).unwrap();
            fingerprint(&m)
        };

        let reused = {
            let mut t = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts()).unwrap();
            let mut m = t.reusable_model();
            t.realize_into(&a, &mut m).unwrap(); // stale A state in the buffers
            t.realize_into(&b, &mut m).unwrap(); // recycled into B
            fingerprint(&m)
        };

        assert_eq!(
            reused, fresh_b,
            "stale-buffer realize_into != fresh realize"
        );
    }

    #[test]
    fn stack_per_zone_contacts_and_priors_drive_in_place() {
        // Z0 gets an OWC mid-zone (oil leg) at a distinct porosity; Z1 is left
        // contactless -> GRV but zero hydrocarbon. Per-zone levels + contacts both bite.
        let draw = RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 7).with_zone_draw(
            ZoneDraw::new(0)
                .with_owc(5015.0)
                .with_priors(0.30, 0.90, 0.20),
        );
        let mut t = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts()).unwrap();
        let m = t.realize(&draw).unwrap();

        // Per-zone property levels applied: Z0 == its override, Z1 == the base prior.
        let por = m.zone_stats(PORO).unwrap();
        assert!((por[0].mean - 0.30).abs() < 1e-12, "Z0 porosity override");
        assert!((por[1].mean - 0.20).abs() < 1e-12, "Z1 base prior");

        let zoned = m.in_place_by_zone().unwrap();
        assert!(zoned.zones[0].in_place.hcpv_m3 > 0.0, "Z0 oil leg");
        assert_eq!(
            zoned.zones[1].in_place.hcpv_m3, 0.0,
            "Z1 contactless -> zero hydrocarbon"
        );
        // Conservation: total GRV == sum of per-zone GRV.
        let sum: f64 = zoned.zones.iter().map(|z| z.in_place.grv_m3).sum();
        assert!((zoned.total.grv_m3 - sum).abs() <= 1e-9 * sum);

        // An out-of-range zone index is a typed error, not a panic.
        let bad = RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 8)
            .with_zone_draw(ZoneDraw::new(9).with_owc(5015.0));
        assert!(t.realize(&bad).is_err());
    }

    // --- zone-scoped property parity (question_zoned_mc_zone_pipe_parity) ---

    use crate::pipeline::{Gaussian, UpscaleMethod, WellLog};
    use crate::StaticModelBuilder;
    use petektools::{Variogram, VariogramModel};

    /// The two-zone stack with an OWC at each zone's base (full oil column) so HCPV
    /// is non-zero and a zone NTG cube actually moves the per-zone in-place.
    fn two_zone_stack_owc() -> HorizonStack {
        let mut s = two_zone_stack();
        let owc = |d: f64| Contact {
            kind: ContactKind::Owc,
            depth_m: d,
            hardness: Hardness::Assumed,
        };
        s.zone_layers[0].contacts = vec![owc(5030.0)]; // Z0 fully oil
        s.zone_layers[1].contacts = vec![owc(5060.0)]; // Z1 fully oil
        s
    }

    /// An NTG pipeline scoped to Z1: two corner wells reading NTG ≈ 0.5 across Z1's
    /// depth range (5030–5060), well below the base prior NTG = 0.8 — so ignoring the
    /// cube (the bug) vs honouring it gives a clearly different Z1 HCPV.
    fn z1_ntg_pipe() -> PropertyPipeline {
        let samples = || vec![(5033.0, 0.5), (5040.0, 0.5), (5048.0, 0.5), (5055.0, 0.5)];
        let wells = vec![
            WellLog::new(60.0, 60.0, samples()),
            WellLog::new(940.0, 940.0, samples()),
        ];
        let vg = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 400.0).unwrap();
        PropertyPipeline::new(NTG)
            .upscale(wells, UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(vg, 42))
    }

    #[test]
    fn stack_volume_bundle_shell_is_non_empty_and_spill_invariant() {
        // question_volume_bundle_stack_empty: a stack model's volume_bundle used to
        // export an EMPTY shell for large (SPILLED) models — the spilled model's
        // grid() is a 1×1×1 placeholder, so the shell read no geometry/cubes and
        // emitted 0 triangles (surfacing as "N cells - 0 tris" in the viewer). The
        // shell now materializes the backing, so in-core and spilled agree.
        //
        // The fixture is a full 8×8×9 block (no truncation), so the exterior shell is
        // exactly the box surface: 2·(ni·nj + ni·nk + nj·nk) quads · 2 tris/quad.
        let (ni, nj, nk) = (8usize, 8, 9);
        let expected_tris = 2 * 2 * (ni * nj + ni * nk + nj * nk);

        let tris = |m: &StaticModel| m.volume_bundle(PORO).unwrap().indices.len() / 3;

        // In-core deterministic build.
        let built = crate::StaticModelBuilder::from_horizon_stack(two_zone_stack(), opts())
            .unwrap()
            .build()
            .unwrap();
        assert!(!built.is_spilled());
        assert_eq!(tris(&built), expected_tris, "in-core stack shell");

        // Spilled (out-of-core) build — the actual bug scenario. A tiny budget forces
        // the slab-incremental spilled path; its grid() is a placeholder.
        let spilled = crate::StaticModelBuilder::from_horizon_stack(two_zone_stack(), opts())
            .unwrap()
            .with_memory_budget(crate::MemoryBudget::bytes(1))
            .build()
            .unwrap();
        assert!(spilled.is_spilled(), "tiny budget must spill");
        assert_eq!(spilled.grid().dims().cell_count(), 1, "placeholder grid");
        assert_eq!(
            tris(&spilled),
            expected_tris,
            "spilled stack shell must match in-core (not empty)"
        );

        // MC template realize path.
        let mut t = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts()).unwrap();
        let realized = t
            .realize(&RealizationDraw::new(
                1_000_000.0,
                0.0,
                0.0,
                0.2,
                0.8,
                0.3,
                1,
            ))
            .unwrap();
        assert_eq!(tris(&realized), expected_tris, "realized stack shell");
    }

    #[test]
    fn tied_template_zero_spread_matches_tied_builder() {
        // task_petekstatic_template_ties: a well tie is DRAW-INVARIANT, so a tied
        // template applied at construction must, at zero spread, reproduce the tied
        // deterministic build bit-for-bit (geometry + cubes + per-zone in-place) — and
        // carry the same tie residual on provenance.
        use crate::WellTie;
        // Tie H1 (the mid horizon, mapped) at node (4,4) to a measured 5022 (its map
        // depth is 5030 -> a +8 m pull-up that keeps the stack ordered, no repair).
        let ties = || vec![WellTie::new("TIE-1", 500.0, 500.0, 4, 4).with_top("H1", 5022.0)];

        let built = crate::StaticModelBuilder::from_horizon_stack(two_zone_stack(), opts())
            .unwrap()
            .with_well_ties(ties())
            .build()
            .unwrap();

        let draw = RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 1);
        let mut t = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts())
            .unwrap()
            .with_well_ties(ties())
            .unwrap();
        let realized = t.realize(&draw).unwrap();

        // Bit-level geometry + cube + per-zone in-place parity.
        assert_eq!(
            fingerprint(&built),
            fingerprint(&realized),
            "tied template realize must equal tied builder build"
        );
        // The tie actually bit (vs an untied build): zone geometry moved.
        let untied = crate::StaticModelBuilder::from_horizon_stack(two_zone_stack(), opts())
            .unwrap()
            .build()
            .unwrap();
        assert_ne!(
            fingerprint(&untied),
            fingerprint(&built),
            "the tie must change the model"
        );
        // Residual surfaced on the realization's provenance (measured 5022 vs 5030).
        let res = &realized.provenance().well_ties;
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].residuals.len(), 1);
        assert!(
            (res[0].residuals[0].residual_m - (5022.0 - 5030.0)).abs() < 1e-6,
            "tie residual = measured - model"
        );
    }

    #[test]
    fn zoned_mc_zero_spread_matches_deterministic_zone_pipe() {
        // A zone-scoped NTG pipe (with_zone_property on Z1). The deterministic build
        // and a ZERO-SPREAD zoned-MC realization of the same template must agree on
        // in_place_by_zone. Before the fix the template realized Z1 from the ZONE
        // PRIOR NTG (0.8), ignoring the upscale+SGS cube (~0.5), so Z1 HCPV differed
        // by ratio ≈ prior/upscaled.
        let built = StaticModelBuilder::from_horizon_stack(two_zone_stack_owc(), opts())
            .unwrap()
            .with_zone_property("Z1", z1_ntg_pipe())
            .build()
            .unwrap();
        let built_zoned = built.in_place_by_zone().unwrap();

        // Zero spread: base priors == opts priors, no per-zone draws (keep static
        // contacts), no property shifts.
        let draw = RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 1);
        let mut t = StaticModelTemplate::from_horizon_stack(two_zone_stack_owc(), opts())
            .unwrap()
            .with_zone_property("Z1", z1_ntg_pipe());
        let realized = t.realize(&draw).unwrap();
        let realized_zoned = realized.in_place_by_zone().unwrap();

        for (b, r) in built_zoned.zones.iter().zip(realized_zoned.zones.iter()) {
            let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(1.0);
            assert!(
                rel(b.in_place.grv_m3, r.in_place.grv_m3) < 1e-9,
                "zone {} GRV: built {} vs realized {}",
                b.zone,
                b.in_place.grv_m3,
                r.in_place.grv_m3
            );
            assert!(
                rel(b.in_place.hcpv_m3, r.in_place.hcpv_m3) < 1e-9,
                "zone {} HCPV: built {} vs realized {} (zone pipe honoured?)",
                b.zone,
                b.in_place.hcpv_m3,
                r.in_place.hcpv_m3
            );
        }
        // The zone pipe actually moved Z1 off the prior (guards against a trivial
        // pass where the cube was never applied on EITHER side).
        let z1_ntg = realized.zone_stats(NTG).unwrap()[1].mean;
        assert!(
            (z1_ntg - 0.8).abs() > 0.1,
            "Z1 NTG {z1_ntg} should reflect the ~0.5 pipe, not the 0.8 prior"
        );
    }

    // --- structural uncertainty white-box (decision_structural_uncertainty_isochore) ---

    fn field(sd: f64, range: f64) -> PerturbationField {
        PerturbationField::new(
            sd,
            Variogram::new(VariogramModel::Spherical, 0.0, 1.0, range).unwrap(),
        )
    }

    #[test]
    fn perturb_stack_surfaces_is_bit_reproducible() {
        // The perturbed surfaces are a pure function of (draw, fixed surfaces): two
        // calls with the same draw give bit-identical node depths at every horizon.
        let draw = RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 5)
            .with_top_structural(field(8.0, 300.0))
            .with_zone_draw(ZoneDraw::new(0).with_isochore_structural(field(6.0, 300.0)));
        let t = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts()).unwrap();
        let stack = t.stack.as_ref().unwrap();
        let (dx, dy) = spacing(draw.area_m2, t.ni, t.nj);
        let a = t.perturb_stack_surfaces(stack, &draw, dx, dy).unwrap();
        let b = t.perturb_stack_surfaces(stack, &draw, dx, dy).unwrap();
        for (sa, sb) in a.iter().zip(b.iter()) {
            for jp in 0..=t.nj {
                for ip in 0..=t.ni {
                    assert_eq!(
                        sa.z(ip, jp),
                        sb.z(ip, jp),
                        "perturbed surface not reproducible"
                    );
                }
            }
        }
        // The top field actually moved node depths off the fixed surface somewhere.
        assert!(
            (0..=t.nj).any(|jp| (0..=t.ni).any(|ip| a[0].z(ip, jp) != stack.surfaces[0].z(ip, jp))),
            "top perturbation did not move the surface"
        );
    }

    #[test]
    fn perturb_pins_ties_to_zero_at_the_tie_node() {
        // A well tie on the TOP horizon pins the top perturbation to zero at its node
        // (perturbation -> 0 at tie wells) while the field still moves other nodes.
        use crate::WellTie;
        let ties = vec![WellTie::new("T", 500.0, 500.0, 4, 4).with_top("H0", 5001.0)];
        let t = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts())
            .unwrap()
            .with_well_ties(ties)
            .unwrap();
        let draw = RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 3)
            .with_top_structural(field(30.0, 400.0));
        let stack = t.stack.as_ref().unwrap();
        let (dx, dy) = spacing(draw.area_m2, t.ni, t.nj);
        let surfs = t.perturb_stack_surfaces(stack, &draw, dx, dy).unwrap();
        // Pinned exactly at the tie node.
        assert_eq!(
            surfs[0].z(4, 4),
            stack.surfaces[0].z(4, 4),
            "top perturbation must be zero at the tie node"
        );
        // Non-zero somewhere off the tie node (the field is live elsewhere).
        assert!(
            (0..=t.nj).any(|jp| (0..=t.ni).any(|ip| {
                (ip, jp) != (4, 4) && surfs[0].z(ip, jp) != stack.surfaces[0].z(ip, jp)
            })),
            "perturbation should be live away from the tie"
        );
    }

    #[test]
    fn perturb_zero_masks_a_fully_merged_zone() {
        // A fully merged zone (Top == Mid → base isochore 0 everywhere) stays merged in
        // every draw: the isochore perturbation is zero-masked, so the deeper horizon
        // never lifts off the shallower one (no resurrection of a collapsed zone).
        let mut s = two_zone_stack();
        s.horizons[1].source = HorizonSource::Mapped(flat_surf(5000.0)); // Mid == Top
        let t = StaticModelTemplate::from_horizon_stack(s, opts()).unwrap();
        let draw = RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 11)
            .with_zone_draw(ZoneDraw::new(0).with_isochore_structural(field(15.0, 400.0)));
        let stack = t.stack.as_ref().unwrap();
        let (dx, dy) = spacing(draw.area_m2, t.ni, t.nj);
        let surfs = t.perturb_stack_surfaces(stack, &draw, dx, dy).unwrap();
        for jp in 0..=t.nj {
            for ip in 0..=t.ni {
                assert_eq!(
                    surfs[1].z(ip, jp),
                    surfs[0].z(ip, jp),
                    "merged zone must stay merged (zero-masked) at ({ip},{jp})"
                );
            }
        }
    }
}
