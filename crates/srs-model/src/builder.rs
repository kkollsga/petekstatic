//! [`StaticModelBuilder`] — the deterministic single-shot build (the
//! geomodeller's entry point): wireframe (or a flat footprint) + optional logs
//! → one [`StaticModel`]. This is the **relocated structural/population half of
//! petekSim's `RefiningModel`** (`srs-core/src/refine.rs`, moved 2026-07-03 per
//! `task_relocate_refine_orchestration`) — minus nothing: per the layer charter
//! the volumetrics tail now also lives on the produced model
//! ([`StaticModel::in_place`]).
//!
//! The build pipeline (SPEC §3): controls → `warm_surface` (the petekTools
//! warm kernel — the SAME structural-solve path the MC template uses, so the two
//! never diverge by kernel; R2) → base-above-top guard (R1) → `layer_grid` →
//! populate (priors or upscaled logs) → `StaticModel`. The cold `solve_surface`
//! remains srs-gridder's accuracy reference, no longer on the build path.

use crate::model::{Georef, StaticModel};
use crate::pipeline::{PropertyPipeline, PropertyReport};
use crate::population::{override_zone_priors, populate, PetroSample};
use crate::provenance::{
    BuildWarning, HorizonTieResidual, InterfaceRepair, PopulationMode, Provenance, StackProvenance,
    WellTieRecord, ZoneProvenance,
};
use crate::spec::{BuildSpec, TieMethod, TieSettings};
use crate::trend::TrendSurface;
use crate::zones::ZoneTable;
use petekstatic_error::StaticError;
use petektools::{grid_min_curvature_conditioned, Conditioning, Lattice, MinCurvatureOperator};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use srs_gridder::{
    layer_grid_stack, Conformity, Control, ExtrapolationPolicy, KernelSurface, SolveOpts,
    StreamingLayering, Surface, ZoneLayerSpec,
};
use srs_spill::{decide_mode, BuildMode, MemoryBudget, SpillNotice};
use srs_volumetrics::{validate_positive, ConstantPriors, NTG, PORO, SW};
use srs_wireframe::{
    Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
};
use std::path::PathBuf;

/// Build settings: the volumetric scalars a wireframe does not carry (area,
/// gross height), the layering, the gridder settings, and the priors.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BuildOpts {
    /// Areal footprint [m²] (sets cell spacing; square).
    pub area_m2: f64,
    /// Gross column thickness [m] (conformable base = top + this).
    pub gross_height_m: f64,
    /// Number of k-layers.
    pub nk: usize,
    /// Layering scheme.
    pub conformity: Conformity,
    /// Cold-solve gridder settings.
    pub solve_opts: SolveOpts,
    /// Day-1 constant priors (fractions).
    pub priors: ConstantPriors,
}

/// A well pick (a top intersection) on the model node lattice: node `(ip, jp)`
/// at `depth_m` (positive-down). The source of a tops-only horizon — a horizon
/// with no mapped surface, defined solely by picks (`HorizonSource::TopsOnly`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Pick {
    /// Areal node index along i (`0..=ni`).
    pub ip: usize,
    /// Areal node index along j (`0..=nj`).
    pub jp: usize,
    /// Pick depth (m, positive down).
    pub depth_m: f64,
}

/// A raw scatter observation in **world** coordinates: easting/northing + a
/// positive-down depth. The first-class input for a horizon defined by an
/// unstructured point set (the shape petekIO delivers) — the engine grids it
/// itself (`HorizonSource::Scatter`), so no caller pre-grids scatter onto the
/// model lattice before the solve sees it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct WorldPoint {
    /// World easting (same CRS as the stack's [`StackFrame`] georef).
    pub x: f64,
    /// World northing.
    pub y: f64,
    /// Depth (m, positive down).
    pub depth_m: f64,
}

/// How a stack horizon's depth surface is sourced.
///
/// The real regional framework uses every shape: a **scattered** world-coordinate
/// point set the engine grids itself (`Scatter`), a normal **mapped** surface
/// already on the model lattice (`Mapped`), and a **tops-only** internal horizon
/// with *no* mapped surface, defined only by well picks (`TopsOnly`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum HorizonSource {
    /// **Raw scatter** in world coordinates. The engine owns the gridding: the
    /// stack build conditions the points onto the model lattice (via the
    /// [`StackFrame`] georef), leaving genuine data voids **undefined** (`NaN`), so
    /// the single structural solve + [`ExtrapolationPolicy`] + isochore build-down
    /// operate on the actual observations — a data-void margin between merged
    /// horizons collapses to zero instead of carrying independently-extrapolated
    /// fill. Only valid through [`StaticModelBuilder::from_scatter_stack`] /
    /// [`crate::StaticModelTemplate::from_scatter_stack`], which condition it to
    /// `Mapped` before the shared resolution path.
    Scatter(Vec<WorldPoint>),
    /// A mapped depth surface **already gridded** on the model lattice: every
    /// defined node (`!depth_m.is_nan()`) becomes a hard control, honoured exactly.
    /// This **bypasses the engine's solve/conditioning fidelity** — supply it only
    /// for genuinely pre-gridded inputs (a loaded grid surface); raw point sets
    /// belong in [`HorizonSource::Scatter`]. Tied or untied is a well-tie concern,
    /// not a source distinction.
    Mapped(GriddedDepth),
    /// **Tops-only**: no mapped surface. The horizon is constructed at build time as
    /// an internal split **subordinate to its mapped envelope**: the mapped horizon
    /// above it **+ `min(pick isochore, envelope isochore)`** — its pick-thickness
    /// field (min-curvature gridded from the picks; constant-thickness fallback for a
    /// single pick) plainly clamped against the enclosing mapped-zone thickness (see
    /// [`resolve_stack_surfaces`]). It therefore lives strictly inside `[mapped top,
    /// mapped base]`; where the envelope merges to zero the split collapses onto both
    /// bounds. A *trailing* tops-only horizon with no mapped horizon beneath it has no
    /// envelope and falls back to the legacy absolute drape ([`drape_tops_only`]). The
    /// geological statement: an untied internal split follows the mapped horizon above
    /// it at the thickness its well picks record, never breaching the measured
    /// horizons that bound its zone.
    TopsOnly(Vec<Pick>),
}

/// An explicit **well tie** for a horizon-stack build (P8 per-horizon ties,
/// `task_petekstatic_multizone_2`): a well at world `(x, y)` sitting on control node
/// `(ip, jp)`, carrying its measured formation **tops** per named horizon. Attach a
/// set with [`StaticModelBuilder::with_well_ties`].
///
/// At build time each measured top named for a **mapped** horizon is applied as an
/// extra hard control at `(ip, jp)` — the seismic surface is *tied* to the well
/// marker — and the **pre-tie residual** (`measured − untied model surface` at the
/// node) is recorded in [`crate::Provenance::well_ties`]. A **tops-only** horizon is
/// already conditioned by its picks, so its tie is recorded as a QC residual against
/// the draped surface (≈ 0 at the pick wells). The world `(x, y)` positions the well
/// marker in the map bundle; `(ip, jp)` is the control-lattice node the tie pins.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WellTie {
    /// Well identifier (fictional in fixtures, e.g. `99/1-1`).
    pub id: String,
    /// World easting of the well marker.
    pub x: f64,
    /// World northing of the well marker.
    pub y: f64,
    /// Control-lattice node the tie pins (`ip` in `0..=ni`, `jp` in `0..=nj`).
    pub ip: usize,
    /// Control-lattice node the tie pins.
    pub jp: usize,
    /// Measured formation tops: `(horizon_name, measured_depth_m)`, positive-down.
    pub tops: Vec<(String, f64)>,
}

impl WellTie {
    /// A well tie at world `(x, y)` on control node `(ip, jp)`, no tops yet.
    #[must_use]
    pub fn new(id: impl Into<String>, x: f64, y: f64, ip: usize, jp: usize) -> Self {
        Self {
            id: id.into(),
            x,
            y,
            ip,
            jp,
            tops: Vec::new(),
        }
    }

    /// Add a measured top for `horizon` at `depth_m` (positive-down). Replaces any
    /// prior top for the same horizon.
    #[must_use]
    pub fn with_top(mut self, horizon: impl Into<String>, depth_m: f64) -> Self {
        let horizon = horizon.into();
        self.tops.retain(|(h, _)| h != &horizon);
        self.tops.push((horizon, depth_m));
        self
    }
}

/// One horizon in an ordered structural stack (top→down).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StackHorizon {
    /// The horizon's name (used in the zone table + provenance + view bundles).
    pub name: String,
    /// Where its depth surface comes from.
    pub source: HorizonSource,
}

/// One zone between two consecutive stack horizons: its layering and its own
/// fluid contacts.
///
/// ## Per-zone contacts (a domain statement)
/// The real regional framework assigns **separate** contact sets to different
/// zones — distinct accumulations in one stack — and some zones have **no**
/// contact at all. So contacts are scoped per zone here, not globally:
/// `in_place` computes each zone's hydrocarbons against *its* contacts. A zone
/// with **no contacts** contributes its gross bulk volume but **zero hydrocarbon
/// in-place** — no contact means no *known* accumulation, so it is explicitly not
/// treated as a full hydrocarbon column.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StackZone {
    /// This zone's name — identifies it in the zone table, provenance, and the view
    /// bundles (the section bundle's per-zone `zones` list is keyed by it). Folds the
    /// old `HorizonStack::zone_names` parallel array into the zone itself.
    pub name: String,
    /// Optional display colour for this zone (a viewer hint, e.g. `"#ffcc00"` or a
    /// named colour). Carried onto the model's [`crate::Zone`] and surfaced in the
    /// section bundle's `zones` list (`{name, color}`) for colour-by-zone rendering.
    /// `None` = the viewer picks a default. Additive; no effect on geometry/volumes.
    pub color: Option<String>,
    /// This zone's conformity/layering style (per-zone, from the conformity wave).
    pub conformity: Conformity,
    /// Requested layer count — honoured only by [`Conformity::Proportional`]
    /// (a Follow style derives it from the zone's own `dz` and thickness).
    pub nk: usize,
    /// This zone's fluid contacts (OWC/GOC/GWC). Empty = a contactless zone (gross
    /// bulk, zero hydrocarbon in-place). A GOC + a lower OWC/GWC makes it a
    /// two-contact (gas-cap + oil-leg) zone that `in_place_by_zone` splits.
    pub contacts: Vec<Contact>,
}

impl StackZone {
    /// A stack zone named `name` with the given layering + contacts and no colour.
    /// Set the display colour with [`StackZone::with_color`].
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        conformity: Conformity,
        nk: usize,
        contacts: Vec<Contact>,
    ) -> Self {
        Self {
            name: name.into(),
            color: None,
            conformity,
            nk,
            contacts,
        }
    }

    /// Set this zone's display colour (a viewer hint carried into the section
    /// bundle's `zones` list).
    #[must_use]
    pub fn with_color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }
}

/// The canonical regional framework as an **ordered horizon stack**: `N` horizons
/// top→down define `N − 1` named zones (intra-zone splits are simply more
/// horizons), each zone carrying its own layering + fluid contacts. The input to
/// [`StaticModelBuilder::from_horizon_stack`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HorizonStack {
    /// Ordered horizons, top→down (`N >= 2`). The first (top) must be `Mapped`
    /// (a tops-only horizon needs a mapped horizon above to drape from).
    pub horizons: Vec<StackHorizon>,
    /// The `N − 1` zones between consecutive horizons, top→down. `zone_layers[z]`
    /// is the zone between `horizons[z]` and `horizons[z + 1]`, carrying its own
    /// name, colour, layering, and contacts (the old `zone_names` parallel array is
    /// folded into [`StackZone::name`]).
    pub zone_layers: Vec<StackZone>,
}

/// The areal lattice + world frame a **scatter** stack is conditioned onto
/// ([`StaticModelBuilder::from_scatter_stack`]): `ni × nj` areal cells
/// (`(ni+1) × (nj+1)` nodes) registered by `georef`. The caller chooses the
/// lattice resolution (a modelling decision); the engine owns the gridding of the
/// raw points onto it. `ni = √area_m2 / georef.spacing_x` for the build's areal
/// scaling to agree with this frame.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StackFrame {
    /// Areal cells along i (node count along i is `ni + 1`).
    pub ni: usize,
    /// Areal cells along j (node count along j is `nj + 1`).
    pub nj: usize,
    /// The world georeference (column-centroid convention; see [`Georef`]).
    pub georef: Georef,
}

/// Deterministic single-shot builder. Mutable state is the growing control set
/// (`add_top_control`) — the live-refine loop; `build()` converges the current
/// state into a fresh [`StaticModel`].
#[derive(Debug, Clone)]
pub struct StaticModelBuilder {
    opts: BuildOpts,
    ni: usize,
    nj: usize,
    controls: Vec<Control>,
    /// Control set for a supplied `Base` horizon (its defined nodes). `None` =
    /// no Base horizon → the base is the constant `gross_height_m` offset.
    base_controls: Option<Vec<Control>>,
    logs: Option<Vec<PetroSample>>,
    /// Optional areal trend multiplier field applied per-column to NTG (and
    /// optionally φ) after population — lateral shape only (`with_areal_trend`).
    trend: Option<TrendSurface>,
    framework: Wireframe,
    /// Entry-point default provenance label (`"flat-box"` / `"wireframe"` /
    /// `"horizon-stack"`), used when the spec does not name one.
    default_inputs_ref: &'static str,
    /// The declarative build configuration — every `with_*` setter is thin sugar
    /// mutating this; [`StaticModelBuilder::with_spec`] installs a whole one. The
    /// same [`BuildSpec`] drives the MC template, so the two never diverge by
    /// config shape (`task_petekstatic_spec_mirror`).
    spec: BuildSpec,
    /// Build-time advisories (unused supplied horizons); copied into provenance.
    warnings: Vec<BuildWarning>,
    /// Per-property geostatistical pipelines run after the base population (P5):
    /// each upscales positioned logs then SGS-propagates its cube, overwriting the
    /// prior/log cube of that property. Empty = no geostatistical population.
    properties: Vec<PropertyPipeline>,
    /// The ordered horizon stack for a multi-zone build (`from_horizon_stack`).
    /// `Some` routes `build()` down the stack path (N horizons → N−1 zones, each
    /// its own conformity/layering); `None` is the classic 2-surface (Top+Base,
    /// single implicit zone) path. The `controls`/`base_controls`/`opts.conformity`
    /// /`opts.nk`/`opts.gross_height_m` fields are unused on the stack path.
    stack: Option<HorizonStack>,
    /// Per-zone geostatistical pipelines (P8 per-zone population,
    /// `task_petekstatic_multizone_2`): `(zone_name, pipeline)`. Each runs after the
    /// base population + whole-model `properties`, **restricted to its zone's
    /// `k`-range** (`PropertyPipeline::apply_in_zone`) so each zone gets its own
    /// variogram / trend / log-conditioning. Zones not named keep the base field.
    zone_properties: Vec<(String, PropertyPipeline)>,
    /// Per-zone constant-prior overrides (`with_zone_priors`): `(zone_name, priors)`.
    /// Applied after the base population, overwriting `PORO`/`NTG`/`SW` inside that
    /// zone's `k`-range — the per-zone distribution level a stack zone owns. A zone
    /// with a `zone_property` on top gets the pipeline's field over this baseline.
    zone_priors: Vec<(String, ConstantPriors)>,
    /// The declared memory budget (ruling R5). Default: [`MemoryBudget::default`]
    /// (a documented fraction of physical RAM). `build()` compares the live-set
    /// estimate against it — below → today's in-core path (byte-identical); above
    /// → the out-of-core spill (geometry + cubes onto a petekTools store, f32).
    budget: MemoryBudget,
    /// Where a spilled build writes its store (ruling R5). `None` = the platform
    /// temp dir (`std::env::temp_dir`).
    spill_dir: Option<PathBuf>,
    /// Keep the spill store past model drop (`with_spill_persist`); default `false`
    /// = the store is a temp file removed when the last model clone drops.
    spill_persist: bool,
}

impl StaticModelBuilder {
    /// A flat-box start: the four corner controls at `top_depth_m`, so the
    /// initial converged model is the flat box (the model-first walking
    /// skeleton). The framework is synthesized (flat Top horizon + one contact).
    ///
    /// # Errors
    /// [`StaticError::Grid`] if a dimension is zero; [`StaticError::InvalidInput`]
    /// if area/height are non-positive or the contact is not finite.
    pub fn flat(
        ni: usize,
        nj: usize,
        top_depth_m: f64,
        contact_depth_m: f64,
        opts: BuildOpts,
    ) -> Result<Self, StaticError> {
        if ni == 0 || nj == 0 || opts.nk == 0 {
            return Err(StaticError::Grid("dimensions must be non-zero".into()));
        }
        validate_positive("area_m2", opts.area_m2)?;
        validate_positive("gross_height_m", opts.gross_height_m)?;
        if !contact_depth_m.is_finite() {
            return Err(StaticError::InvalidInput(format!(
                "contact depth must be finite, got {contact_depth_m}"
            )));
        }
        let controls: Vec<Control> = [(0, 0), (ni, 0), (0, nj), (ni, nj)]
            .iter()
            .map(|&(ip, jp)| Control {
                ip,
                jp,
                z: top_depth_m,
            })
            .collect();
        let framework = synth_framework(ni + 1, nj + 1, top_depth_m, contact_depth_m);
        Ok(Self {
            opts,
            ni,
            nj,
            controls,
            base_controls: None,
            logs: None,
            trend: None,
            framework,
            default_inputs_ref: "flat-box",
            spec: BuildSpec::default(),
            warnings: Vec::new(),
            properties: Vec::new(),
            stack: None,
            zone_properties: Vec::new(),
            zone_priors: Vec::new(),
            budget: MemoryBudget::default(),
            spill_dir: None,
            spill_persist: false,
        })
    }

    /// Seed from a constraining [`Wireframe`] — the data-layer hand-off. The
    /// wireframe drives the **structure** (its `Top` horizon becomes the control
    /// set; the first fluid contact sets the column base); the scalars it does
    /// not carry come from `opts`. Grid `(ni, nj)` is the top-surface lattice
    /// (`ncol-1, nrow-1`).
    ///
    /// # Errors
    /// [`StaticError`] if the wireframe has no `Top` horizon, the top surface is
    /// degenerate or fully undefined, no fluid contact is present, or `opts`
    /// dimensions/scalars are degenerate.
    pub fn from_wireframe(wf: &Wireframe, opts: BuildOpts) -> Result<Self, StaticError> {
        // Structure extraction also asserts the wireframe carries a contact —
        // the produced model's `in_place` needs one on the framework.
        let s = wireframe_structure(wf, &opts)?;
        Ok(Self {
            opts,
            ni: s.ni,
            nj: s.nj,
            controls: s.top_controls,
            base_controls: s.base_controls,
            logs: None,
            trend: None,
            framework: wf.clone(),
            default_inputs_ref: "wireframe",
            spec: BuildSpec::default(),
            warnings: s.warnings,
            properties: Vec::new(),
            stack: None,
            zone_properties: Vec::new(),
            zone_priors: Vec::new(),
            budget: MemoryBudget::default(),
            spill_dir: None,
            spill_persist: false,
        })
    }

    /// Seed from an ordered **horizon stack** — the canonical regional framework:
    /// `N` horizons top→down define `N − 1` named zones, each with its own
    /// [`Conformity`], layer count, and fluid contacts. This is the multi-zone
    /// generalization of [`Self::from_wireframe`]; the classic Top+Base wireframe
    /// is its 2-horizon / single-zone degenerate case.
    ///
    /// The stack may mix all three horizon source kinds ([`HorizonSource`]): mapped
    /// surfaces, mapped-but-untied surfaces, and tops-only internal horizons draped
    /// conformally from the mapped horizon above at pick-controlled thickness.
    ///
    /// `opts` supplies the areal footprint (`area_m2`), the gridder settings, and
    /// the fallback priors; its `nk` / `conformity` / `gross_height_m` are **unused**
    /// on this path (per-zone layering + the last horizon's surface govern instead).
    /// The builder methods (`with_logs`, `with_property`, `with_min_thickness_m`,
    /// `with_georef`, …) apply as usual; `with_min_thickness_m` becomes the
    /// **per-interface** order-repair floor.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the stack has fewer than 2 horizons, the
    /// zone-count arrays are not `N − 1`, the first (top) horizon is not `Mapped`,
    /// a tops-only horizon has no mapped horizon above it, or a mapped surface is
    /// degenerate / lattice-mismatched / fully undefined.
    pub fn from_horizon_stack(stack: HorizonStack, opts: BuildOpts) -> Result<Self, StaticError> {
        let (ni, nj) = validate_stack(&stack)?;
        Ok(Self {
            opts,
            ni,
            nj,
            controls: Vec::new(),
            base_controls: None,
            logs: None,
            trend: None,
            framework: Wireframe::from_boundary(srs_wireframe::Boundary {
                ring: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
                hardness: Hardness::Interpolated,
            }),
            default_inputs_ref: "horizon-stack",
            spec: BuildSpec::default(),
            warnings: Vec::new(),
            properties: Vec::new(),
            stack: Some(stack),
            zone_properties: Vec::new(),
            zone_priors: Vec::new(),
            budget: MemoryBudget::default(),
            spill_dir: None,
            spill_persist: false,
        })
    }

    /// Build from a stack that may carry **raw scatter** horizons
    /// ([`HorizonSource::Scatter`]) — the engine owns the gridding. Each scatter
    /// horizon is conditioned onto the `frame` lattice (snap-average onto the
    /// nearest node via the frame georef, genuine voids left `NaN`), then the stack
    /// resolves through the **same** merge-safe path as [`Self::from_horizon_stack`]
    /// (converged solve + [`ExtrapolationPolicy`] + isochore build-down over the
    /// actual observations). The `frame` georef is registered as the model's world
    /// frame. `Mapped` / `TopsOnly` horizons pass through untouched.
    ///
    /// This is the single scatter-gridding authority: no caller pre-grids scatter
    /// onto the lattice, so a data-void margin between two merged horizons collapses
    /// to zero instead of carrying independently-extrapolated fill.
    ///
    /// > **Conditioning seam.** Scatter is conditioned onto the frame lattice by
    /// > **bilinear-weighted** minimum curvature (petekTools [`Conditioning::Bilinear`]):
    /// > each off-node datum maps to its *fractional* node position and is honoured
    /// > through the bilinear interpolation of its four surrounding nodes rather than
    /// > snapped to the nearest — holding on-data fidelity to the lattice
    /// > representation floor (`grid_scatter` owns the world→lattice conversion).
    ///
    /// # Errors
    /// [`StaticError`] as [`Self::from_horizon_stack`], plus
    /// [`StaticError::InvalidInput`] if a scatter horizon has no points, or
    /// [`StaticError::Grid`] if it lands no point on the frame (fully undefined).
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

    /// **Condition a raw-scatter stack onto `frame` ONCE** — the expensive
    /// per-horizon cold bilinear solve (`task_suite_scatter_perf`) — returning the
    /// conditioned, all-[`HorizonSource::Mapped`] stack as a **shared handle**.
    ///
    /// This is the dedup seam for callers that build **both** a
    /// [`StaticModel`](crate::StaticModel) and its MC
    /// [`StaticModelTemplate`](crate::StaticModelTemplate) from the same scatter:
    /// conditioning is **draw-invariant** and identical across the builder and the
    /// template, so re-running [`Self::from_scatter_stack`] /
    /// [`StaticModelTemplate::from_scatter_stack`](crate::StaticModelTemplate::from_scatter_stack)
    /// on each re-solves the whole 11-horizon cold gridding redundantly (the
    /// canonical build's dominant cost). Condition once, then feed the returned
    /// stack to both paths **without** re-conditioning:
    ///
    /// ```ignore
    /// let g = frame.georef;
    /// let conditioned = StaticModelBuilder::condition_scatter_stack(stack, &frame)?;
    /// // model:
    /// let model = StaticModelBuilder::from_horizon_stack(conditioned.clone(), opts)?
    ///     .with_georef(g.origin_x, g.origin_y, g.spacing_x, g.spacing_y)
    ///     .build()?;
    /// // MC template (byte-for-byte the same geometry `realize` reproduces):
    /// let template = StaticModelTemplate::from_horizon_stack(conditioned, opts)?
    ///     .with_georef(g.origin_x, g.origin_y, g.spacing_x, g.spacing_y);
    /// ```
    ///
    /// The result is **bit-identical** to conditioning inside each entry point —
    /// this only removes the redundant re-solve. Genuine data voids stay `NaN`.
    ///
    /// # Errors
    /// As [`Self::from_scatter_stack`]'s conditioning step: [`StaticError::InvalidInput`]
    /// if a scatter horizon has no points, [`StaticError::Grid`] if one lands no
    /// point on the frame.
    pub fn condition_scatter_stack(
        mut stack: HorizonStack,
        frame: &StackFrame,
    ) -> Result<HorizonStack, StaticError> {
        condition_scatter(&mut stack, frame)?;
        Ok(stack)
    }

    /// Opt into clamping a crossed base (base above top) to zero gross at the
    /// offending columns instead of erroring. Default is to error
    /// ([`StaticError::CrossedSurfaces`]); this zeroes only those columns and
    /// leaves the rest untouched (R1).
    #[must_use]
    pub fn with_clamp_base_to_top(mut self, clamp: bool) -> Self {
        self.spec.clamp_base_to_top = clamp;
        self
    }

    /// Opt into **post-gridding order-repair**: where the gridded base sits less
    /// than `min_thickness_m` below the top (a thin or crossed column), pull the
    /// base **down** to exactly `top + min_thickness_m`, preserving the top (the
    /// better-constrained seismic pick). Independent gridding of Top and Base can
    /// overshoot at thin margins and undo a pointwise pre-repair; this repairs the
    /// gridded result rather than erroring. Off by default (the crossing guard
    /// stays the default); when enabled, [`Provenance::warnings`] records a
    /// [`BuildWarning::ThinColumnsRepaired`] with the repaired-node count and the
    /// worst violation. Takes precedence over [`Self::with_clamp_base_to_top`] (R-c).
    #[must_use]
    pub fn with_min_thickness_m(mut self, min_thickness_m: f64) -> Self {
        self.spec.min_thickness_m = Some(min_thickness_m);
        self
    }

    /// Opt into the **cell-collapse threshold** (Petrel-style): after layering, any
    /// cell thinner than `collapse_below_m` collapses to zero thickness, its sliver
    /// merged into a thicker **zone-interior** neighbour so rock is conserved (never
    /// deleted, never merged across a zone boundary). Off by default. When enabled,
    /// [`Provenance::warnings`] records a [`BuildWarning::CellsCollapsed`] with the
    /// count. It bites hardest on sub-`dz` slivers under fixed-count proportional
    /// layering in thin columns and on Follow-style truncation partials.
    #[must_use]
    pub fn with_collapse_below_m(mut self, collapse_below_m: f64) -> Self {
        self.spec.collapse_below_m = Some(collapse_below_m);
        self
    }

    /// Set the **extrapolation policy** for the horizon-stack build — how every
    /// solved stack surface/isochore behaves **beyond its data hull**. Default
    /// [`ExtrapolationPolicy::DecayToData`] (dip held for `start_cells`, then a
    /// linear decay to the nearest-data value over `decay_cells`) — conservative,
    /// never silent unbounded natural-dip into a data void. Pass
    /// [`ExtrapolationPolicy::NaturalDip`] to opt back into the legacy unbounded
    /// linear extension (appropriate only when the regional dip is KNOWN to
    /// continue). Stack path only; the classic 2-surface path is unchanged.
    #[must_use]
    pub fn with_extrapolation(mut self, policy: ExtrapolationPolicy) -> Self {
        self.spec.extrapolation = policy;
        self
    }

    /// Set a gas-cap connate-water override applied to gas-zone cells in a
    /// two-contact `in_place` split (R3), so a single shared `SW` cube does not
    /// over-state gas-cap OGIP. No effect on a single-contact (no-GOC) column.
    #[must_use]
    pub fn with_sw_gas(mut self, sw_gas: f64) -> Self {
        self.spec.sw_gas = Some(sw_gas);
        self
    }

    /// Attach positioned petro samples (TVD, φ, Sw) so `build()` populates cells
    /// from upscaled logs instead of constant priors; cells with no samples keep
    /// the priors. Empty = priors everywhere.
    #[must_use]
    pub fn with_logs(mut self, samples: Vec<PetroSample>) -> Self {
        self.logs = if samples.is_empty() {
            None
        } else {
            // Sort by TVD once so population binary-searches each cell (V2).
            Some(crate::population::sort_by_tvd(samples))
        };
        self
    }

    /// Attach an areal trend multiplier field (external-drift-lite): a gridded
    /// lateral shape resampled to the model lattice, mean-normalized, and applied
    /// per-column to NTG (and φ if the trend flags it) after population. The
    /// trend gives lateral *shape*; the prior/log gives the *level*
    /// (`decision_staticmodel_regen_seam`). See [`TrendSurface`].
    ///
    /// **Superseded (deprecation-tracked).** This is the **interim** post-population
    /// multiplier hook; the fuller path is [`PropertyPipeline`] with
    /// [`crate::Gaussian::with_trend`], which steers the *simulation itself* by collocated
    /// (Markov-1) cokriging rather than scaling a populated cube. It is retained for
    /// the simple "lateral NTG/φ shape on a constant-prior fill" case and is still on
    /// the build + MC-template path; a hard `#[deprecated]` waits until the pipeline
    /// subsumes this use (organize wave P10 review, `task_petekstatic_organize`).
    #[must_use]
    pub fn with_areal_trend(mut self, trend: TrendSurface) -> Self {
        self.trend = Some(trend);
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
    /// internal spec). Values are identical either way, so the build — and every
    /// determinism contract over it — is unchanged.
    #[must_use]
    pub fn with_spec(mut self, spec: BuildSpec) -> Self {
        self.spec = spec;
        self
    }

    /// The effective provenance label: the spec's, or the entry point's default.
    fn inputs_ref_string(&self) -> String {
        self.spec
            .inputs_ref
            .clone()
            .unwrap_or_else(|| self.default_inputs_ref.to_string())
    }

    /// Declare the memory budget (ruling R5). Below it `build()` stays on the
    /// in-core path (byte-identical); above it the model spills its geometry +
    /// cubes to a petekTools store and reads through mmap windows. Pass
    /// [`MemoryBudget::unlimited`] to force in-core, or [`MemoryBudget::bytes`] to
    /// force spill at a test scale.
    #[must_use]
    pub fn with_memory_budget(mut self, budget: MemoryBudget) -> Self {
        self.budget = budget;
        self
    }

    /// Where a spilled build writes its store (ruling R5). Default: the platform
    /// temp dir. The store is removed on model drop unless [`Self::with_spill_persist`].
    #[must_use]
    pub fn with_spill_dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.spill_dir = Some(dir.into());
        self
    }

    /// Keep a spilled model's store file past drop (a caller-owned location)
    /// instead of the default temp-file cleanup.
    #[must_use]
    pub fn with_spill_persist(mut self, persist: bool) -> Self {
        self.spill_persist = persist;
        self
    }

    /// Register the built model's **world** georeference: the world `(x, y)` of
    /// column `(0, 0)`'s centroid and the world column spacing. The upstream seam
    /// (which knows the source horizon lattice / CRS) supplies this so the view
    /// bundles emit ONE consistent world frame — the map raster overlays the world
    /// outline + wells, and a world fence / bore trace sections correctly (it maps
    /// world traces through the same `xy↔ij` as well registration).
    ///
    /// Non-finite / non-positive spacing is ignored (the model stays in the local
    /// degenerate frame). This does **not** alter the grid geometry (GRV is still
    /// the area-scaled local square); it only labels the local column lattice with
    /// its world frame. See [`Georef`].
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

    /// Register the **world** areal boundary ring for a horizon-stack build — the
    /// study-area closure the map bundle overlays (`world [x, y]`, first == last point
    /// to close). The `from_horizon_stack` path has no wireframe to source a boundary
    /// from, so without this it would emit the degenerate unit square while the frame
    /// and wells are world coordinates, collapsing the viewer's content extent
    /// (`task_petekstatic_zoned_fixes` finding 2). When omitted, the build still emits a
    /// world-extent **rectangle** derived from the georef (never the unit square against
    /// a world frame); pass the real closure ring here for the true outline shape.
    /// Ignored on the `from_wireframe` / `flat` paths (they carry the wireframe's own
    /// boundary).
    #[must_use]
    pub fn with_boundary(mut self, ring: Vec<[f64; 2]>) -> Self {
        self.spec.boundary = if ring.is_empty() { None } else { Some(ring) };
        self
    }

    /// Attach a per-property geostatistical pipeline (P5): `build()` runs it after
    /// the base population, upscaling its positioned logs and SGS-propagating the
    /// cube — overwriting that property's prior/log cube with the conditioned,
    /// simulated field. Multiple pipelines model properties one at a time (each is
    /// its own [`PropertyPipeline`]). The pipeline's [`crate::PropertyReport`] lands
    /// on the model's [`Provenance::property_reports`].
    #[must_use]
    pub fn with_property(mut self, pipeline: PropertyPipeline) -> Self {
        self.properties.push(pipeline);
        self
    }

    /// Attach a geostatistical pipeline **scoped to one zone** of a horizon-stack
    /// build (P8 per-zone population, `task_petekstatic_multizone_2`): after the base
    /// population and any whole-model [`Self::with_property`] pipelines, `build()`
    /// runs it restricted to `zone_name`'s `k`-range
    /// ([`PropertyPipeline::apply_in_zone`]), so each zone gets its own
    /// variogram/trend/log-conditioning and only that zone's slice is overwritten.
    /// Attach one pipeline per zone; the pipeline's report lands on
    /// [`Provenance::property_reports`]. A `zone_name` absent from the built zone
    /// table is a [`StaticError::InvalidInput`] at build time.
    #[must_use]
    pub fn with_zone_property(
        mut self,
        zone_name: impl Into<String>,
        pipeline: PropertyPipeline,
    ) -> Self {
        self.zone_properties.push((zone_name.into(), pipeline));
        self
    }

    /// Override the constant priors **inside one zone** of a horizon-stack build
    /// (P8 per-zone population): after the base population `build()` overwrites
    /// `PORO`/`NTG`/`SW` across `zone_name`'s `k`-range with `priors` — the per-zone
    /// distribution level a stack zone owns (a sand zone vs a shale zone). A
    /// [`Self::with_zone_property`] pipeline on the same zone then simulates over this
    /// baseline. A `zone_name` absent from the built zone table is a
    /// [`StaticError::InvalidInput`] at build time.
    #[must_use]
    pub fn with_zone_priors(
        mut self,
        zone_name: impl Into<String>,
        priors: ConstantPriors,
    ) -> Self {
        self.zone_priors.push((zone_name.into(), priors));
        self
    }

    /// Set how the well ties are applied ([`TieSettings`]: datum substitution —
    /// the default — or bounded-radius locality). Thin sugar over the spec; both
    /// the ties and their settings are read at `build()` time, so the call order
    /// vs [`StaticModelBuilder::with_well_ties`] does not matter here.
    #[must_use]
    pub fn with_tie_settings(mut self, ties: TieSettings) -> Self {
        self.spec.ties = ties;
        self
    }

    /// Attach explicit **well ties** for a horizon-stack build (P8 per-horizon ties,
    /// `task_petekstatic_multizone_2`): each [`WellTie`]'s measured tops tie their
    /// **mapped** horizons to the well (the top is added as a hard control at the
    /// well's node and the surface re-solved), and every tie's **pre-tie residual**
    /// (`measured − untied model surface`) is recorded in
    /// [`crate::Provenance::well_ties`] and surfaced in the map bundle's
    /// `wells[].ties`. Tops for a **tops-only** horizon are recorded as QC residuals
    /// against the pick-conditioned drape. Honoured on the `from_horizon_stack` path.
    #[must_use]
    pub fn with_well_ties(mut self, ties: Vec<WellTie>) -> Self {
        self.spec.well_ties = ties;
        self
    }

    /// Opt into **sugar-cube** section rendering (flat-box cells) — default `false`,
    /// where the section bundle carries per-edge cell depths and the viewer draws
    /// dip-following trapezoids. The engine geometry is corner-point either way; this
    /// only sets the section view's `sugar_cube` flag (and flattens its edge arrays).
    #[must_use]
    pub fn with_sugar_cube(mut self, sugar_cube: bool) -> Self {
        self.spec.sugar_cube = sugar_cube;
        self
    }

    /// Add a top-surface depth control point on the `(ni+1) x (nj+1)` node
    /// lattice — the new datum the next `build()` re-converges to honour.
    pub fn add_top_control(&mut self, ip: usize, jp: usize, depth_m: f64) {
        self.controls.push(Control { ip, jp, z: depth_m });
    }

    /// Number of control points currently honoured.
    #[must_use]
    pub fn control_count(&self) -> usize {
        self.controls.len()
    }

    /// Converge the model at the current control set: cold minimum-curvature
    /// solve → conformable layering → population → a populated [`StaticModel`].
    ///
    /// # Errors
    /// [`StaticError`] if the surface solve, layering, or population fails.
    pub fn build(&self) -> Result<StaticModel, StaticError> {
        if self.stack.is_some() {
            return self.build_stack();
        }
        let nx = self.ni + 1;
        let ny = self.nj + 1;
        // Structural solve in petekTools kernel space (the SAME path the template
        // uses), so the two builds never diverge by kernel (R2).
        let top_k = warm_surface(nx, ny, &self.controls)?;
        let top = top_k.surface();
        // Base follows a supplied `Base` horizon's real relief if present;
        // otherwise it is the conformable constant `gross_height_m` offset (the
        // backward-compatible fallback).
        let base = match &self.base_controls {
            Some(bc) => warm_surface(nx, ny, bc)?.surface().clone(),
            None => top.offset_by(self.opts.gross_height_m),
        };
        // R1/R-c: a base that crosses above the top collapses GRV. Opt-in
        // `min_thickness_m` repairs it (pull the base to a minimum thickness below
        // the top, record a warning); otherwise error (default) or clamp to zero
        // gross.
        let mut warnings = self.warnings.clone();
        let base = match self.spec.min_thickness_m {
            Some(min_t) => {
                let (repaired, columns, worst_m) = base.repair_min_thickness(top, min_t)?;
                if columns > 0 {
                    warnings.push(BuildWarning::ThinColumnsRepaired { columns, worst_m });
                }
                repaired
            }
            None => base.guard_below(top, self.spec.clamp_base_to_top)?,
        };
        let (dx, dy) = spacing(self.opts.area_m2, self.ni, self.nj);

        // v2 slab-incremental spilled build (true O(slab) PEAK RSS): when the
        // population is the per-slab-safe prior/constant family and collapse is off,
        // stream ZCORN + cubes k-slab-by-k-slab straight into the store — never
        // materialize a whole in-core grid (the O(grid) transient v1 paid). The
        // in-core path (below budget) is untouched, byte-identical.
        if self.streaming_eligible() {
            let streaming = StreamingLayering::prepare(
                &[top, &base],
                dx,
                dy,
                &[ZoneLayerSpec {
                    conformity: self.opts.conformity,
                    requested_nk: self.opts.nk,
                }],
            )?;
            let (mode, estimate) = decide_mode(streaming.dims(), 3, self.budget);
            if mode == BuildMode::Spilled {
                return self.build_spilled_streaming(&streaming, estimate, warnings);
            }
        }

        // Layering: under a Follow conformity `nk` is dz-derived and thin columns
        // truncate against the pinch-out horizon; read the effective nk + report
        // back (`opts.nk` is only honoured by `Proportional`). The single-zone
        // stack path also carries the optional cell-collapse pass.
        let layered = layer_grid_stack(
            &[top, &base],
            dx,
            dy,
            &[ZoneLayerSpec {
                conformity: self.opts.conformity,
                requested_nk: self.opts.nk,
            }],
            self.spec.collapse_below_m,
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
        let mut grid = layered.grid;
        populate(
            &mut grid,
            self.opts.priors,
            self.logs.as_deref(),
            self.trend.as_ref(),
        )?;

        // P5 per-property geostatistical population: each pipeline upscales its logs
        // and SGS-propagates its cube, overwriting the base prior/log cube.
        let mut property_reports = Vec::with_capacity(self.properties.len());
        for pipe in &self.properties {
            property_reports.push(pipe.apply_with_georef(&mut grid, self.spec.georef)?);
        }

        // P8 per-zone population: the single implicit zone is the whole column.
        let zone_kranges = vec![("RESERVOIR".to_string(), 0..nk)];
        self.apply_zone_population(&mut grid, &zone_kranges, &mut property_reports)?;

        let population = self.population_mode();
        self.maybe_spill(
            StaticModel::new(
                self.framework.clone(),
                grid,
                ZoneTable::single(nk),
                Provenance {
                    inputs_ref: self.inputs_ref_string(),
                    solve_opts: self.opts.solve_opts,
                    conformity: self.opts.conformity,
                    nk,
                    population,
                    realization: None,
                    warnings,
                    property_reports,
                    stack: None,
                    well_ties: Vec::new(),
                    sugar_cube: self.spec.sugar_cube,
                },
                self.spec.sw_gas,
            )
            .with_georef_opt(self.spec.georef),
        )
    }

    /// Apply the out-of-core mode decision (ruling R5) to a freshly-built in-core
    /// model: below the budget → return it unchanged (byte-identical); above →
    /// stream its geometry + cubes to a petekTools store (f32 lanes, R4), emit the
    /// loud mode-switch advisory, and return the spilled model (its in-core grid
    /// buffers dropped). The build's transient peak is `O(grid)` in v1 (a
    /// fully-slab-incremental gridder build is the R2 follow-up); the spilled
    /// model's steady-state resident set is `O(slab)`.
    fn maybe_spill(&self, model: StaticModel) -> Result<StaticModel, StaticError> {
        let dims = model.grid().dims();
        let n_cubes = model.property_names().len();
        let (mode, estimate) = decide_mode(dims, n_cubes, self.budget);
        if mode == BuildMode::InCore {
            return Ok(model);
        }
        let dir = self.spill_dir.clone().unwrap_or_else(std::env::temp_dir);
        let backing = srs_spill::spill_grid(model.grid(), &dir, !self.spill_persist)?;
        SpillNotice {
            cells: dims.cell_count(),
            budget_bytes: self.budget.limit_bytes(),
            estimate_bytes: estimate,
            store_path: backing.store_path().to_path_buf(),
        }
        .warn();
        Ok(model.into_spilled(backing))
    }

    /// Whether the single-zone build can take the **slab-incremental** streaming
    /// spill path (true O(slab) build peak). The streaming producer computes ZCORN
    /// on demand from two interface planes and populates each slab with the constant
    /// priors, so it supports exactly the per-slab-safe population: **no** cell
    /// collapse (a whole-zone-column carry), **no** logs / areal trend / SGS property
    /// pipelines / per-zone pipelines (those fill whole cubes and re-read them). Such
    /// builds fall back to build-then-spill (`maybe_spill`, O(grid) transient).
    fn streaming_eligible(&self) -> bool {
        self.spec.collapse_below_m.is_none()
            && self.logs.is_none()
            && self.trend.is_none()
            && self.properties.is_empty()
            && self.zone_properties.is_empty()
    }

    /// The slab-incremental spilled build: drive the [`StreamingLayering`] producer
    /// through the store's `slab_mut_f32` writer, emitting each ZCORN slab + its
    /// constant-prior cube slabs in place. Peak working set is **one slab + two
    /// interface planes** (`O(slab)`), never a whole in-core grid — the R2 promise.
    /// Bit-identical to build-then-spill (same f64 `boundary_depth` narrowed to f32,
    /// same constant cubes).
    fn build_spilled_streaming(
        &self,
        streaming: &StreamingLayering,
        estimate: u64,
        mut warnings: Vec<BuildWarning>,
    ) -> Result<StaticModel, StaticError> {
        let dims = streaming.dims();
        let nk = dims.nk;
        // The constant-prior cube triple, name-sorted → byte-deterministic store
        // layout identical to build-then-spill (`cube_lane_names` sorts likewise).
        let mut cube_names = vec![NTG.to_string(), PORO.to_string(), SW.to_string()];
        cube_names.sort();
        let priors = self.opts.priors;
        let const_for = |name: &str| -> f32 {
            (if name == PORO {
                priors.porosity
            } else if name == NTG {
                priors.net_to_gross
            } else {
                priors.water_saturation
            }) as f32
        };

        let mut coord = Vec::new();
        streaming.fill_coord(&mut coord);
        let plane_len = dims.pillar_count();
        let mut plane_top = vec![0.0f64; plane_len];
        let mut plane_bot = vec![0.0f64; plane_len];
        let mut truncated = 0usize;

        let dir = self.spill_dir.clone().unwrap_or_else(std::env::temp_dir);
        let path = srs_spill::unique_spill_path(&dir);
        let backing = srs_spill::spill_streaming(
            &path,
            dims,
            &coord,
            &cube_names,
            !self.spill_persist,
            |k, out| {
                truncated += streaming.fill_zcorn_slab(k, &mut plane_top, &mut plane_bot, out);
                Ok(())
            },
            |name, _k, out| {
                out.fill(const_for(name));
                Ok(())
            },
        )?;

        if truncated > 0 {
            warnings.push(BuildWarning::LayersTruncated { cells: truncated });
        }
        if streaming.nk_capped() {
            warnings.push(BuildWarning::LayerCountCapped { nk });
        }
        SpillNotice {
            cells: dims.cell_count(),
            budget_bytes: self.budget.limit_bytes(),
            estimate_bytes: estimate,
            store_path: backing.store_path().to_path_buf(),
        }
        .warn();

        Ok(StaticModel::spilled(
            self.framework.clone(),
            ZoneTable::single(nk),
            Provenance {
                inputs_ref: self.inputs_ref_string(),
                solve_opts: self.opts.solve_opts,
                conformity: self.opts.conformity,
                nk,
                population: self.population_mode(),
                realization: None,
                warnings,
                property_reports: Vec::new(),
                stack: None,
                well_ties: Vec::new(),
                sugar_cube: self.spec.sugar_cube,
            },
            self.spec.sw_gas,
            self.spec.georef,
            backing,
        ))
    }

    /// The population mode carried into provenance: `Logs` when any log / whole-model
    /// pipeline / per-zone pipeline informed the cubes, else `Priors` (constant priors,
    /// including a per-zone constant-prior override — still a prior, not a log).
    fn population_mode(&self) -> PopulationMode {
        if self.logs.is_some() || !self.properties.is_empty() || !self.zone_properties.is_empty() {
            PopulationMode::Logs
        } else {
            PopulationMode::Priors
        }
    }

    /// P8 per-zone population (`task_petekstatic_multizone_2`): apply the per-zone
    /// constant-prior overrides, then the per-zone geostatistical pipelines, each
    /// restricted to its zone's `k`-range (resolved from `zone_kranges`). Appends the
    /// per-zone pipeline reports to `reports`.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if a named zone is absent from the model's zone
    /// table, or an override/pipeline fails.
    fn apply_zone_population(
        &self,
        grid: &mut srs_grid::Grid,
        zone_kranges: &[(String, core::ops::Range<usize>)],
        reports: &mut Vec<PropertyReport>,
    ) -> Result<(), StaticError> {
        let find = |name: &str| -> Result<core::ops::Range<usize>, StaticError> {
            zone_kranges
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, r)| r.clone())
                .ok_or_else(|| {
                    StaticError::InvalidInput(format!(
                        "per-zone population: zone '{name}' is not in the model's zone table"
                    ))
                })
        };
        for (name, priors) in &self.zone_priors {
            override_zone_priors(grid, *priors, find(name)?)?;
        }
        for (name, pipe) in &self.zone_properties {
            let report = pipe
                .apply_in_zone_with_georef(grid, find(name)?, self.spec.georef)
                .map_err(|e| StaticError::InvalidInput(format!("zone '{name}': {e}")))?;
            reports.push(report);
        }
        Ok(())
    }

    /// The multi-zone build (`from_horizon_stack`): resolve each horizon surface
    /// top→down (draping tops-only horizons conformally), repair consecutive-horizon
    /// ordering per interface, layer each zone by its own conformity, populate, and
    /// assemble a [`StaticModel`] whose real [`ZoneTable`] + per-zone contacts
    /// describe the stack.
    fn build_stack(&self) -> Result<StaticModel, StaticError> {
        let stack = self.stack.as_ref().expect("build_stack requires a stack");
        let nx = self.ni + 1;
        let ny = self.nj + 1;

        // 1. Resolve every horizon surface by BUILDING DOWN via non-negative
        //    isochores ([`resolve_stack_surfaces`]): the top is gridded once, each
        //    deeper mapped horizon = the horizon above + a clamped gridded zone
        //    isochore, and a tops-only internal split = the mapped horizon above +
        //    `min(pick isochore, envelope isochore)`. Ordering can never invert and
        //    a derived surface can never displace a mapped one — where the inputs
        //    merge the isochore is 0 and the zone collapses to genuine zero.
        let mut surfaces = resolve_stack_surfaces(stack, nx, ny, self.spec.extrapolation)?;
        let is_derived: Vec<bool> = stack
            .horizons
            .iter()
            .map(|h| matches!(h.source, HorizonSource::TopsOnly(_)))
            .collect();

        // 1b. Per-horizon WELL TIES (P8): record the pre-tie residual of every
        //     measured top against the untied surface, then re-solve each mapped
        //     horizon with its tie tops added as hard controls (the surface is tied
        //     to the wells). Applied BEFORE order-repair so the repair still
        //     guarantees ordering after the surfaces move.
        let well_ties = self.apply_well_ties(stack, &mut surfaces, nx, ny)?;

        // 2. Per-interface order-repair, honouring mapped-over-derived precedence
        //    (a derived surface yields to a mapped one). With the isochore
        //    construction a consistent stack never crosses here, so this is a
        //    safety net that only fires under `with_min_thickness_m` or genuinely
        //    inconsistent inputs.
        let mut interface_repairs = Vec::new();
        for i in 0..surfaces.len() - 1 {
            let (upper, lower) = repair_interface(
                &surfaces[i],
                &surfaces[i + 1],
                self.spec.min_thickness_m,
                self.spec.clamp_base_to_top,
                is_derived[i],
                is_derived[i + 1],
                i,
                &mut interface_repairs,
            )?;
            surfaces[i] = upper;
            surfaces[i + 1] = lower;
        }

        // 3. Layer each zone by its own conformity; total nk = sum of per-zone.
        let specs: Vec<ZoneLayerSpec> = stack
            .zone_layers
            .iter()
            .map(|z| ZoneLayerSpec {
                conformity: z.conformity,
                requested_nk: z.nk,
            })
            .collect();
        let surf_refs: Vec<&Surface> = surfaces.iter().collect();
        let (dx, dy) = spacing(self.opts.area_m2, self.ni, self.nj);
        let stacked = layer_grid_stack(&surf_refs, dx, dy, &specs, self.spec.collapse_below_m)?;

        let mut warnings = self.warnings.clone();
        if stacked.truncated_cells > 0 {
            warnings.push(BuildWarning::LayersTruncated {
                cells: stacked.truncated_cells,
            });
        }
        if stacked.collapsed_cells > 0 {
            warnings.push(BuildWarning::CellsCollapsed {
                cells: stacked.collapsed_cells,
            });
        }
        if stacked.nk_capped {
            warnings.push(BuildWarning::LayerCountCapped { nk: stacked.nk });
        }

        let mut grid = stacked.grid;
        populate(
            &mut grid,
            self.opts.priors,
            self.logs.as_deref(),
            self.trend.as_ref(),
        )?;
        let mut property_reports = Vec::with_capacity(self.properties.len());
        for pipe in &self.properties {
            property_reports.push(pipe.apply_with_georef(&mut grid, self.spec.georef)?);
        }

        // Zone names + colours are carried by each `StackZone` now (the old parallel
        // `zone_names` array folded in).
        let zone_names: Vec<String> = stack.zone_layers.iter().map(|z| z.name.clone()).collect();
        let zone_colors: Vec<Option<String>> =
            stack.zone_layers.iter().map(|z| z.color.clone()).collect();

        // P8 per-zone population: per-zone priors then per-zone pipelines, each
        // restricted to its stack zone's k-range (variogram/trend/logs per zone).
        let zone_kranges: Vec<(String, core::ops::Range<usize>)> = zone_names
            .iter()
            .cloned()
            .zip(stacked.zones.iter().map(srs_gridder::StackedZone::k_range))
            .collect();
        self.apply_zone_population(&mut grid, &zone_kranges, &mut property_reports)?;

        // 4. Framework wireframe from the resolved surfaces (Top / Intermediate /
        //    Base roles) + the union of per-zone contacts (for the legacy
        //    whole-model path); per-zone contacts drive `in_place_by_zone`.
        let horizon_names: Vec<String> = stack.horizons.iter().map(|h| h.name.clone()).collect();
        let zone_contacts: Vec<Vec<Contact>> = stack
            .zone_layers
            .iter()
            .map(|z| z.contacts.clone())
            .collect();
        let all_contacts: Vec<Contact> = zone_contacts.iter().flatten().copied().collect();
        let ring = stack_boundary_ring(self.spec.boundary.as_ref(), self.spec.georef, nx, ny);
        let framework = stack_framework(&stack.horizons, &surfaces, all_contacts, nx, ny, ring);

        // 5. Real zone table + per-zone provenance.
        let per_zone_nk: Vec<usize> = stacked.zones.iter().map(|z| z.nk).collect();
        let zones = ZoneTable::from_stack(&horizon_names, &zone_names, &zone_colors, &per_zone_nk);
        let zone_prov: Vec<ZoneProvenance> = stacked
            .zones
            .iter()
            .enumerate()
            .map(|(z, sz)| ZoneProvenance {
                name: zone_names[z].clone(),
                top_horizon: horizon_names[z].clone(),
                base_horizon: horizon_names[z + 1].clone(),
                conformity: sz.conformity,
                nk: sz.nk,
                k_start: sz.k_start,
                truncated_cells: sz.truncated_cells,
            })
            .collect();
        let stack_prov = StackProvenance {
            horizons: horizon_names,
            zones: zone_prov,
            interface_repairs,
        };

        let population = self.population_mode();
        self.maybe_spill(
            StaticModel::new(
                framework,
                grid,
                zones,
                Provenance {
                    inputs_ref: self.inputs_ref_string(),
                    solve_opts: self.opts.solve_opts,
                    conformity: self.opts.conformity,
                    nk: stacked.nk,
                    population,
                    realization: None,
                    warnings,
                    property_reports,
                    stack: Some(stack_prov),
                    well_ties,
                    sugar_cube: self.spec.sugar_cube,
                },
                self.spec.sw_gas,
            )
            .with_georef_opt(self.spec.georef)
            .with_zone_contacts(Some(zone_contacts)),
        )
    }

    /// Apply the supplied [`WellTie`]s (P8 per-horizon ties) per the spec's
    /// [`TieSettings`] — see [`substitute_tie_datums`] for the tie math (datum
    /// substitution / bounded-radius locality) — then re-run the full stack
    /// resolution over the tied data, so the ties flow through the SAME
    /// construction as the map data (cumulative isochores, merge collapse,
    /// extrapolation policy) instead of a separate surface re-solve. Returns the
    /// per-well tie records for provenance.
    ///
    /// # Errors
    /// [`StaticError`] if a tie names an unknown horizon, its node is off-lattice,
    /// the radius setting is degenerate, or a re-solve fails.
    fn apply_well_ties(
        &self,
        stack: &HorizonStack,
        surfaces: &mut Vec<Surface>,
        nx: usize,
        ny: usize,
    ) -> Result<Vec<WellTieRecord>, StaticError> {
        if self.spec.well_ties.is_empty() {
            return Ok(Vec::new());
        }
        let (dx, dy) = spacing(self.opts.area_m2, self.ni, self.nj);
        let mut tied_stack = stack.clone();
        let (records, substituted) = substitute_tie_datums(
            &mut tied_stack,
            &self.spec.well_ties,
            surfaces,
            self.spec.ties,
            (dx, dy),
            nx,
            ny,
        )?;
        // Re-run the full stack resolution over the tied data so the model honours
        // both the map and the wells through ONE construction.
        if substituted {
            *surfaces = resolve_stack_surfaces(&tied_stack, nx, ny, self.spec.extrapolation)?;
        }
        Ok(records)
    }
}

/// Fold each [`WellTie`]'s measured tops into the stack's **mapped** gridded
/// datums per the [`TieSettings`], recording every pre-tie residual (`measured −
/// untied surface` at the well node) — the single tie-math authority shared by
/// the deterministic builder and the MC template.
///
/// - [`TieMethod::Replace`] (default — the historical behaviour): the measured
///   top REPLACES the map datum at the tie node (or defines a previously-
///   undefined node). On a fully-defined lattice every other node is still a
///   hard datum, so the tie moves exactly the tied node (radius of influence 0
///   cells); on a sparse lattice its influence is the solver's interpolation
///   reach, bounded beyond the hull by the extrapolation policy.
/// - [`TieMethod::Radius`]: the tie's residual is blended into every **defined**
///   datum of that horizon within `radius_m` of the tie node with a linear decay
///   (weight `1` at the well → the tie node lands exactly on the measured top;
///   `0` at `radius_m`); datums beyond the radius are bit-untouched, and
///   **undefined** nodes stay undefined except the tie node itself (defined at
///   the measured top), so genuine data voids remain the solve's to taper.
///   Distances are metres on the `(dx, dy)` areal node spacing.
///
/// A tops-only horizon's residual is recorded against its pick-conditioned drape
/// but the drape is not re-tied (the picks already condition it). Returns the
/// per-well records + whether any datum changed (the caller re-resolves then).
///
/// # Errors
/// [`StaticError::Grid`] for an off-lattice tie node; [`StaticError::InvalidInput`]
/// for an unknown horizon name or a non-finite / non-positive `radius_m`.
pub(crate) fn substitute_tie_datums(
    stack: &mut HorizonStack,
    ties: &[WellTie],
    surfaces: &[Surface],
    settings: TieSettings,
    (dx, dy): (f64, f64),
    nx: usize,
    ny: usize,
) -> Result<(Vec<WellTieRecord>, bool), StaticError> {
    if let TieMethod::Radius { radius_m } = settings.method {
        if !(radius_m.is_finite() && radius_m > 0.0) {
            return Err(StaticError::InvalidInput(format!(
                "tie radius must be finite and positive, got {radius_m}"
            )));
        }
    }
    let index_of = |name: &str, horizons: &[StackHorizon]| -> Option<usize> {
        horizons.iter().position(|h| h.name == name)
    };
    let mut substituted = false;
    let mut records = Vec::with_capacity(ties.len());
    for tie in ties {
        if tie.ip >= nx || tie.jp >= ny {
            return Err(StaticError::Grid(format!(
                "well tie '{}' node ({},{}) is off the {nx}x{ny} control lattice",
                tie.id, tie.ip, tie.jp
            )));
        }
        let mut residuals = Vec::with_capacity(tie.tops.len());
        for (horizon, measured) in &tie.tops {
            let h = index_of(horizon, &stack.horizons).ok_or_else(|| {
                StaticError::InvalidInput(format!(
                    "well tie '{}' names unknown horizon '{horizon}'",
                    tie.id
                ))
            })?;
            let model_depth_m = surfaces[h].z(tie.ip, tie.jp);
            let residual_m = *measured - model_depth_m;
            residuals.push(HorizonTieResidual {
                horizon: horizon.clone(),
                measured_depth_m: *measured,
                model_depth_m,
                residual_m,
            });
            // Only a mapped horizon is re-tied (a tops-only horizon is already
            // conditioned by its picks).
            if let HorizonSource::Mapped(gd) = &mut stack.horizons[h].source {
                let tie_idx = tie.jp * gd.ncol + tie.ip;
                match settings.method {
                    TieMethod::Replace => {
                        gd.depth_m[tie_idx] = *measured;
                        gd.is_control[tie_idx] = true;
                        substituted = true;
                    }
                    TieMethod::Radius { radius_m } => {
                        for jp in 0..ny {
                            for ip in 0..nx {
                                let idx = jp * gd.ncol + ip;
                                if gd.depth_m[idx].is_nan() {
                                    continue; // voids stay the solve's to taper
                                }
                                let dist = ((ip as f64 - tie.ip as f64) * dx)
                                    .hypot((jp as f64 - tie.jp as f64) * dy);
                                if dist <= radius_m {
                                    let w = 1.0 - dist / radius_m;
                                    gd.depth_m[idx] += residual_m * w;
                                }
                            }
                        }
                        // The tie node itself always lands on the measured top —
                        // exactly (w = 1 over a defined datum), or by defining a
                        // previously-undefined node.
                        gd.depth_m[tie_idx] = *measured;
                        gd.is_control[tie_idx] = true;
                        substituted = true;
                    }
                }
            }
        }
        records.push(WellTieRecord {
            id: tie.id.clone(),
            x: tie.x,
            y: tie.y,
            ip: tie.ip,
            jp: tie.jp,
            residuals,
        });
    }
    Ok((records, substituted))
}

/// Areal cell spacing `(dx, dy)` [m] from a square area footprint (m²).
pub(crate) fn spacing(area_m2: f64, ni: usize, nj: usize) -> (f64, f64) {
    let side = area_m2.sqrt();
    (side / ni as f64, side / nj as f64)
}

/// Validate a [`HorizonStack`] and return the areal cell dims `(ni, nj)` derived
/// from its (mapped) top horizon's lattice.
/// Condition every [`HorizonSource::Scatter`] horizon in `stack` onto the `frame`
/// lattice, replacing it with the resulting sparse [`HorizonSource::Mapped`]
/// surface (genuine voids `NaN`, so the shared resolution path solves + tapers
/// them). `Mapped` / `TopsOnly` horizons are left untouched. The single
/// scatter-gridding site.
///
/// # Errors
/// [`StaticError::InvalidInput`] if a scatter horizon carries no points;
/// [`StaticError::Grid`] if none of its points land on the frame (fully undefined).
pub(crate) fn condition_scatter(
    stack: &mut HorizonStack,
    frame: &StackFrame,
) -> Result<(), StaticError> {
    // Each scatter horizon conditions **independently** onto the same frame — no
    // cross-horizon coupling — and each conditioning is a cold biharmonic solve
    // (the dominant per-horizon cost). Grid them in parallel across horizons
    // (`task_suite_scatter_perf`): the map is embarrassingly parallel and every
    // solve is deterministic given its own points + frame, so the conditioned
    // result is bit-identical to the serial loop (no floating-point reduction is
    // shared across horizons). Results are collected back in horizon order.
    let prof = std::env::var_os("SRS_PROFILE").is_some();
    stack
        .horizons
        .par_iter_mut()
        .try_for_each(|h| -> Result<(), StaticError> {
            let HorizonSource::Scatter(points) = &h.source else {
                return Ok(());
            };
            if points.is_empty() {
                return Err(StaticError::InvalidInput(format!(
                    "scatter horizon '{}' has no points to grid",
                    h.name
                )));
            }
            let t0 = std::time::Instant::now();
            let gd = grid_scatter(points, frame);
            if prof {
                let ctrl = gd.is_control.iter().filter(|&&b| b).count();
                eprintln!(
                    "[SRS_PROFILE] condition_scatter horizon='{}' points={} support_nodes={ctrl} ms={:.1}",
                    h.name,
                    points.len(),
                    t0.elapsed().as_secs_f64() * 1e3,
                );
            }
            if gd.depth_m.iter().all(|z| z.is_nan()) {
                return Err(StaticError::Grid(format!(
                    "scatter horizon '{}' landed no point on the {}x{} frame lattice \
                     (points outside the georeferenced extent)",
                    h.name,
                    frame.ni + 1,
                    frame.nj + 1
                )));
            }
            h.source = HorizonSource::Mapped(gd);
            Ok(())
        })
}

/// A **factor-once / solve-many** bilinear minimum-curvature conditioner for one
/// scatter horizon (petekTools [`MinCurvatureOperator`]) plus its data-support hull
/// mask — the adopted direct-solve seam (`task_suite_scatter_perf`).
///
/// The sample `(x, y)` geometry — and hence the band-LU factorization — is **fixed**;
/// only the depths (the solve RHS) vary. [`factor`](Self::factor) assembles and
/// factors the conditioning operator once (the dominant cost — the ~0.44 s the direct
/// solve replaced the cap-bound ~60 s SOR with, at the canonical ~40k-sample
/// density); [`resolve`](Self::resolve) back-substitutes a fresh depth vector for
/// ~6 ms. That split is the MC reuse lever: re-seating the **same** sample geometry
/// with new depths across realizations is a resolve, not a re-conditioning. A single
/// build conditions each horizon once — the lever is realised only by a caller that
/// re-solves fixed geometry (a data-perturbation MC mode); today's additive-field
/// structural perturbation re-seats surfaces without re-conditioning, so it never hits
/// this path. See CHANGELOG / SPEC §3.
///
/// [`factor`](Self::factor) returns `None` when no point lands on the frame or the
/// assembled system is degenerate/singular; the caller then falls back to the
/// iterative kernel, matching [`grid_min_curvature_conditioned`]'s own dispatch.
pub(crate) struct ScatterConditioner {
    op: MinCurvatureOperator,
    /// Node-major (`jp * nx + ip`) data-support mask — each kept point's four
    /// bracketing nodes (its bilinear footprint).
    support: Vec<bool>,
    nx: usize,
    ny: usize,
}

impl ScatterConditioner {
    /// Factor the bilinear conditioning operator for `points`' fixed sample geometry
    /// on the `frame` node lattice. Returns the conditioner and the kept-point
    /// **depths** (in the operator's sample order, the order [`resolve`](Self::resolve)
    /// expects), or `None` when nothing lands / the system is degenerate.
    pub(crate) fn factor(points: &[WorldPoint], frame: &StackFrame) -> Option<(Self, Vec<f64>)> {
        let (nx, ny) = (frame.ni + 1, frame.nj + 1);
        let g = &frame.georef;
        // Off-node controls as fractional node positions; support marks each point's
        // four bracketing nodes (its bilinear footprint). Depths are kept in the same
        // order as the sample positions — the RHS `resolve` re-solves.
        let mut sample_xy: Vec<[f64; 2]> = Vec::with_capacity(points.len());
        let mut depths: Vec<f64> = Vec::with_capacity(points.len());
        let mut support = vec![false; nx * ny];
        for p in points {
            if !(p.x.is_finite() && p.y.is_finite() && p.depth_m.is_finite()) {
                continue;
            }
            // World → fractional node: the +0.5 undoes the column-centroid origin so
            // a point at node `ip` maps to exactly `ip`.
            let fi = (p.x - g.origin_x) / g.spacing_x + 0.5;
            let fj = (p.y - g.origin_y) / g.spacing_y + 0.5;
            if fi < -0.5 || fj < -0.5 || fi > nx as f64 - 0.5 || fj > ny as f64 - 0.5 {
                continue; // off the frame — no data support here
            }
            sample_xy.push([fi, fj]);
            depths.push(p.depth_m);
            // Bilinear footprint: the 2×2 nodes bracketing (fi, fj), clamped in-frame.
            let i0 = (fi.floor() as isize).clamp(0, nx as isize - 1) as usize;
            let j0 = (fj.floor() as isize).clamp(0, ny as isize - 1) as usize;
            for jp in j0..=(j0 + 1).min(ny - 1) {
                for ip in i0..=(i0 + 1).min(nx - 1) {
                    support[jp * nx + ip] = true;
                }
            }
        }
        if sample_xy.is_empty() {
            return None;
        }
        let lattice = Lattice::regular(0.0, 0.0, 1.0, 1.0, nx, ny);
        let op = MinCurvatureOperator::factor(&lattice, &sample_xy, Conditioning::Bilinear).ok()?;
        Some((
            ScatterConditioner {
                op,
                support,
                nx,
                ny,
            },
            depths,
        ))
    }

    /// Back-substitute `depths` (in the sample order [`factor`](Self::factor) returned
    /// them) through the factored operator, returning the **support-hull-masked**
    /// node field (`jp * nx + ip` order; voids `NaN`). `None` if the solve fails
    /// (the RHS length must equal the factored sample count).
    pub(crate) fn resolve(&self, depths: &[f64]) -> Option<Vec<f64>> {
        let field = self.op.solve(depths).ok()?;
        let mut depth_m = vec![f64::NAN; self.nx * self.ny];
        for jp in 0..self.ny {
            for ip in 0..self.nx {
                if self.support[jp * self.nx + ip] {
                    depth_m[jp * self.nx + ip] = field[[ip, jp]];
                }
            }
        }
        Some(depth_m)
    }
}

/// Grid raw world scatter onto the frame node lattice by **bilinear-conditioned
/// minimum curvature** (petekTools [`Conditioning::Bilinear`]): each point maps to
/// its fractional node position via the georef (column-centroid convention — node
/// `ip` sits at world `origin_x + (ip − 0.5)·spacing_x`), and an off-node datum is
/// honoured through the bilinear interpolation of its four surrounding nodes
/// rather than snapped to the nearest — eliminating the metres-level snap error
/// dense sub-node scatter otherwise carries. Points off the frame are dropped.
///
/// The solve is the **direct band-LU [`MinCurvatureOperator`]** ([`ScatterConditioner`],
/// factor-once / solve-many), falling back to the iterative kernel
/// ([`grid_min_curvature_conditioned`]) only on a degenerate/singular system — the
/// same dispatch the one-shot kernel makes internally, so the conditioned field is
/// bit-identical to the prior direct path.
///
/// The result is **masked to the data-support hull**: only nodes inside the
/// bilinear footprint of some point (the four nodes bracketing it) carry a value
/// and `is_control = true`; every other node stays `NaN` (a genuine void). So the
/// single surface solve + [`ExtrapolationPolicy`] in [`resolve_stack_surfaces`]
/// re-solve and taper the void instead of inheriting an unbounded fill — the
/// margin-collapse the whole path exists for.
///
/// This is the scatter **conditioning** step (authoring observations onto nodes);
/// the surface **gridding** (converged solve + isochore build-down) stays solely
/// in [`resolve_stack_surfaces`].
fn grid_scatter(points: &[WorldPoint], frame: &StackFrame) -> GriddedDepth {
    let (nx, ny) = (frame.ni + 1, frame.nj + 1);
    // Factor the conditioning operator for this horizon's sample geometry, then
    // solve for its depths (the one-shot use of the factor-once seam). `None` means
    // nothing landed or the system is degenerate → fall back to the iterative kernel.
    let (depth_m, support) = match ScatterConditioner::factor(points, frame) {
        Some((cond, depths)) => match cond.resolve(&depths) {
            Some(field) => (field, cond.support),
            None => grid_scatter_sor_fallback(points, frame, cond.support),
        },
        None => (vec![f64::NAN; nx * ny], vec![false; nx * ny]),
    };
    GriddedDepth {
        ncol: nx,
        nrow: ny,
        depth_m,
        is_control: support,
    }
}

/// Iterative-kernel fallback for [`grid_scatter`]: reached only when the direct
/// factorization or its solve fails on a degenerate system. Rebuilds the world →
/// fractional-node `coords` and runs the SOR path (which honours the same
/// [`Conditioning::Bilinear`] classification), masked to the same support hull.
fn grid_scatter_sor_fallback(
    points: &[WorldPoint],
    frame: &StackFrame,
    support: Vec<bool>,
) -> (Vec<f64>, Vec<bool>) {
    let (nx, ny) = (frame.ni + 1, frame.nj + 1);
    let g = &frame.georef;
    let mut coords: Vec<[f64; 3]> = Vec::with_capacity(points.len());
    for p in points {
        if !(p.x.is_finite() && p.y.is_finite() && p.depth_m.is_finite()) {
            continue;
        }
        let fi = (p.x - g.origin_x) / g.spacing_x + 0.5;
        let fj = (p.y - g.origin_y) / g.spacing_y + 0.5;
        if fi < -0.5 || fj < -0.5 || fi > nx as f64 - 0.5 || fj > ny as f64 - 0.5 {
            continue;
        }
        coords.push([fi, fj, p.depth_m]);
    }
    let mut depth_m = vec![f64::NAN; nx * ny];
    if !coords.is_empty() {
        let lattice = Lattice::regular(0.0, 0.0, 1.0, 1.0, nx, ny);
        if let Ok(field) =
            grid_min_curvature_conditioned(&coords, &lattice, None, Conditioning::Bilinear)
        {
            for jp in 0..ny {
                for ip in 0..nx {
                    if support[jp * nx + ip] {
                        depth_m[jp * nx + ip] = field[[ip, jp]];
                    }
                }
            }
        }
    }
    (depth_m, support)
}

pub(crate) fn validate_stack(stack: &HorizonStack) -> Result<(usize, usize), StaticError> {
    let n = stack.horizons.len();
    if n < 2 {
        return Err(StaticError::InvalidInput(format!(
            "a horizon stack needs at least 2 horizons, got {n}"
        )));
    }
    if stack.zone_layers.len() != n - 1 {
        return Err(StaticError::InvalidInput(format!(
            "a {n}-horizon stack needs {} zones, got {}",
            n - 1,
            stack.zone_layers.len()
        )));
    }
    let (nx, ny) = match &stack.horizons[0].source {
        HorizonSource::Mapped(gd) => (gd.ncol, gd.nrow),
        HorizonSource::TopsOnly(_) => {
            return Err(StaticError::InvalidInput(
                "the first (top) horizon must be a mapped surface".into(),
            ))
        }
        HorizonSource::Scatter(_) => {
            return Err(StaticError::InvalidInput(
                "raw-scatter horizons must be built through `from_scatter_stack` \
                 (they are conditioned to the model lattice before resolution)"
                    .into(),
            ))
        }
    };
    if nx < 2 || ny < 2 {
        return Err(StaticError::Grid(format!(
            "top surface needs at least 2x2 nodes, got {nx}x{ny}"
        )));
    }
    let mut seen_mapped = false;
    for h in &stack.horizons {
        match &h.source {
            HorizonSource::Scatter(_) => {
                return Err(StaticError::InvalidInput(format!(
                    "raw-scatter horizon '{}' must be conditioned via `from_scatter_stack`",
                    h.name
                )))
            }
            HorizonSource::Mapped(gd) => {
                if gd.ncol != nx || gd.nrow != ny {
                    return Err(StaticError::InvalidInput(format!(
                        "mapped horizon '{}' lattice {}x{} does not match the top {}x{}",
                        h.name, gd.ncol, gd.nrow, nx, ny
                    )));
                }
                if gd.depth_m.iter().all(|z| z.is_nan()) {
                    return Err(StaticError::Grid(format!(
                        "mapped horizon '{}' is fully undefined (all NaN)",
                        h.name
                    )));
                }
                seen_mapped = true;
            }
            HorizonSource::TopsOnly(picks) => {
                if !seen_mapped {
                    return Err(StaticError::InvalidInput(format!(
                        "tops-only horizon '{}' has no mapped horizon above to drape from",
                        h.name
                    )));
                }
                if picks.is_empty() {
                    return Err(StaticError::InvalidInput(format!(
                        "tops-only horizon '{}' needs at least one pick",
                        h.name
                    )));
                }
                for p in picks {
                    if p.ip >= nx || p.jp >= ny {
                        return Err(StaticError::Grid(format!(
                            "tops-only horizon '{}' pick ({},{}) is off the {nx}x{ny} lattice",
                            h.name, p.ip, p.jp
                        )));
                    }
                }
            }
        }
    }
    Ok((nx - 1, ny - 1))
}

/// Grid the **non-negative thickness field** (metres, row-major `jp*nx+ip`) a
/// tops-only horizon's picks record below the surface `above`: each pick's
/// thickness `pick.depth − above.z(pick)` is clamped `>= 0` and min-curvature
/// gridded (a single pick → a constant field); the gridded result is clamped
/// `>= 0` per node and tapered beyond the picks per the [`ExtrapolationPolicy`].
/// This is a pure isochore (a thickness), *not* a seated surface — the caller
/// decides how to place it (the trailing absolute drape below, or the envelope
/// min-clamp in [`resolve_stack_surfaces`]).
pub(crate) fn pick_thickness_field(
    picks: &[Pick],
    above: &Surface,
    nx: usize,
    ny: usize,
    policy: ExtrapolationPolicy,
) -> Result<Vec<f64>, StaticError> {
    if picks.len() == 1 {
        let p = picks[0];
        let t = (p.depth_m - above.z(p.ip, p.jp)).max(0.0);
        return Ok(vec![t; nx * ny]);
    }
    let controls: Vec<Control> = picks
        .iter()
        .map(|p| Control {
            ip: p.ip,
            jp: p.jp,
            z: (p.depth_m - above.z(p.ip, p.jp)).max(0.0),
        })
        .collect();
    let field_surf = warm_surface(nx, ny, &controls)?
        .surface()
        .taper_beyond_data(&controls, policy);
    let mut field = vec![0.0; nx * ny];
    for jp in 0..ny {
        for ip in 0..nx {
            field[jp * nx + ip] = field_surf.z(ip, jp).max(0.0);
        }
    }
    Ok(field)
}

/// Legacy **absolute** tops-only drape: seat the pick thickness field directly
/// below `above`. Only used for a *trailing* tops-only horizon with no mapped
/// horizon beneath it (no envelope to be subordinate to). Inside a mapped
/// envelope the isochore construction ([`resolve_stack_surfaces`]) min-clamps the
/// split within `[mapped top, mapped base]` instead, so a derived surface can
/// never push a mapped one down.
pub(crate) fn drape_tops_only(
    picks: &[Pick],
    above: &Surface,
    nx: usize,
    ny: usize,
    policy: ExtrapolationPolicy,
) -> Result<Surface, StaticError> {
    let field = pick_thickness_field(picks, above, nx, ny, policy)?;
    above.offset_by_field(&field)
}

/// Thickness samples between two gridded horizons at their **co-located** defined
/// nodes: `(offset[idx] + (zb − za).max(0))` at every node where both are defined.
/// `offset` (may be `None` = zero) lifts an adjacent-pair thickness into
/// thickness-below-TOP space via the shallower horizon's already-built cumulative
/// field (the fallback chain in [`resolve_stack_surfaces`]).
fn colocated_thickness_samples(
    gd_a: &GriddedDepth,
    gd_b: &GriddedDepth,
    offset: Option<&[f64]>,
    nx: usize,
    ny: usize,
) -> Vec<Control> {
    let mut samples = Vec::new();
    for jp in 0..ny {
        for ip in 0..nx {
            let idx = jp * nx + ip;
            let (za, zb) = (gd_a.depth_m[idx], gd_b.depth_m[idx]);
            if !za.is_nan() && !zb.is_nan() {
                let base = offset.map_or(0.0, |o| o[idx]);
                samples.push(Control {
                    ip,
                    jp,
                    z: base + (zb - za).max(0.0),
                });
            }
        }
    }
    samples
}

/// Resolve every stack surface by **building down via non-negative cumulative
/// isochores anchored at the TOP horizon** — the merge-safe, swap-invariant
/// construction that replaces per-horizon independent gridding (which let
/// independent-regridding overshoot manufacture crossings the inputs never
/// contained, and let a tops-only absolute drape push a mapped base down):
///
/// 1. The TOP mapped horizon is gridded once (kernel space; its defined nodes are
///    hard controls, honoured exactly) and tapered beyond its data per `policy`.
/// 2. Each deeper MAPPED horizon carries a **cumulative isochore vs the TOP
///    anchor**: thickness samples `(z_k − z_top).max(0)` at the co-located
///    defined nodes, min-curvature gridded, clamped `>= 0`, tapered per `policy`,
///    then made **monotone** down the stack (`cum_k = max(cum_{k−1}, iso_k)`).
///    The horizon is seated at `top + cum_k`. Consequences by construction:
///    - ordering can never invert (monotone cums);
///    - an exact merge in the inputs samples exactly `0` → the zone collapses to
///      genuine zero;
///    - **internal-swap invariance**: each mapped horizon's cum depends only on
///      the TOP grid and its OWN grid, so replacing an internal mapped horizon
///      leaves every other mapped horizon bit-unchanged (except where the swapped
///      horizon's own data genuinely crosses a deeper horizon — the monotone
///      clamp then collapses that zone toward the data rather than adding rock).
///
///    Fallback chain when a horizon shares NO defined node with the TOP:
///    co-locate with the nearest shallower mapped horizon (samples lifted by
///    its cum), else (fully disjoint) an independent solve differenced
///    pointwise vs the horizon above, clamped `>= 0`.
/// 3. A tops-only internal horizon inside a mapped envelope = the mapped horizon
///    above + `min(pick isochore, envelope isochore)` — a plain clamp, so
///    envelope zero → split zero → both sub-zones collapse; the derived surface
///    can never displace either mapped bound. A *trailing* tops-only with no
///    mapped horizon beneath falls back to the absolute drape
///    ([`drape_tops_only`]).
///
/// Every solve is tapered beyond its data hull per the [`ExtrapolationPolicy`]
/// (default [`ExtrapolationPolicy::DecayToData`]) — never silent unbounded
/// natural-dip into a data void.
#[allow(clippy::needless_range_loop)] // k indexes both `stack.horizons` and `surfaces`
pub(crate) fn resolve_stack_surfaces(
    stack: &HorizonStack,
    nx: usize,
    ny: usize,
    policy: ExtrapolationPolicy,
) -> Result<Vec<Surface>, StaticError> {
    let n = stack.horizons.len();
    // Mapped-horizon indices (validate_stack guarantees the first horizon is
    // mapped and every tops-only has a mapped horizon above it).
    let mapped_idx: Vec<usize> = (0..n)
        .filter(|&i| matches!(stack.horizons[i].source, HorizonSource::Mapped(_)))
        .collect();
    let top_i = mapped_idx[0];
    let HorizonSource::Mapped(gd_top) = &stack.horizons[top_i].source else {
        unreachable!("mapped_idx only indexes Mapped horizons");
    };
    let top_controls = controls_from_surface(gd_top);
    if top_controls.is_empty() {
        return Err(StaticError::Grid(format!(
            "mapped horizon '{}' is fully undefined (all NaN)",
            stack.horizons[top_i].name
        )));
    }
    let top = warm_surface(nx, ny, &top_controls)?
        .surface()
        .taper_beyond_data(&top_controls, policy);

    let mut surfaces: Vec<Option<Surface>> = vec![None; n];
    let mut cums: Vec<Option<Vec<f64>>> = vec![None; n];
    surfaces[top_i] = Some(top.clone());
    cums[top_i] = Some(vec![0.0; nx * ny]);

    // Pre-solve the PRIMARY-path isochores in PARALLEL. Each mapped envelope's
    // primary solve co-locates vs the TOP anchor (`colocated_thickness_samples(gd_top,
    // gd_b, None)`) and warm-solves that isochore — depending only on the top + the
    // envelope's own gridded datums, never on another envelope. So the (dominant)
    // warm-start solves are independent and run concurrently; the serial build-down
    // below consumes each pre-solved field where the primary path applies and falls
    // back to a serial solve only for the co-location-fallback / degenerate cases,
    // which DO carry a cross-envelope cumulative dependency. Bit-identical to the
    // serial solve — same inputs, same field (`warm_surface` is deterministic).
    let windows: Vec<[usize; 2]> = mapped_idx.windows(2).map(|w| [w[0], w[1]]).collect();
    let primary_iso: Vec<Option<Vec<f64>>> = windows
        .par_iter()
        .map(|&[_, b]| -> Result<Option<Vec<f64>>, StaticError> {
            let HorizonSource::Mapped(gd_b) = &stack.horizons[b].source else {
                unreachable!("mapped_idx only indexes Mapped horizons");
            };
            let samples = colocated_thickness_samples(gd_top, gd_b, None, nx, ny);
            if samples.is_empty() {
                return Ok(None); // primary co-location empty → serial fallback below
            }
            let solved = warm_surface(nx, ny, &samples)?
                .surface()
                .taper_beyond_data(&samples, policy);
            Ok(Some(
                (0..nx * ny)
                    .map(|idx| {
                        let (ip, jp) = (idx % nx, idx / nx);
                        solved.z(ip, jp).max(0.0)
                    })
                    .collect(),
            ))
        })
        .collect::<Result<Vec<_>, StaticError>>()?;

    // Build DOWN one mapped envelope (a → b) at a time, seating any tops-only
    // splits strictly inside it.
    for (wi, w) in mapped_idx.windows(2).enumerate() {
        let (a, b) = (w[0], w[1]);
        let HorizonSource::Mapped(gd_b) = &stack.horizons[b].source else {
            unreachable!("mapped_idx only indexes Mapped horizons");
        };
        let cum_a = cums[a].clone().expect("upper anchor cum built");
        let iso: Vec<f64> = if let Some(pre) = primary_iso[wi].clone() {
            // (1) Primary path: the isochore vs the TOP anchor (the invariance
            //     convention), warm-solved in the parallel pre-pass above.
            pre
        } else {
            // (2) Fallback: co-locate with the nearest shallower mapped horizon,
            //     lifting its samples by that horizon's built cumulative field.
            let mut samples: Vec<Control> = Vec::new();
            for &a2 in mapped_idx.iter().filter(|&&i| i < b).rev() {
                if a2 == top_i {
                    break; // (1) already tried the top anchor (empty → here)
                }
                let HorizonSource::Mapped(gd_a2) = &stack.horizons[a2].source else {
                    unreachable!("mapped_idx only indexes Mapped horizons");
                };
                let cum_a2 = cums[a2].as_deref().expect("shallower cum built");
                samples = colocated_thickness_samples(gd_a2, gd_b, Some(cum_a2), nx, ny);
                if !samples.is_empty() {
                    break;
                }
            }
            if samples.is_empty() {
                // (3) Degenerate (no co-location anywhere above): independent solve,
                //     pointwise vs the built horizon above, clamped >= 0.
                let controls_b = controls_from_surface(gd_b);
                if controls_b.is_empty() {
                    return Err(StaticError::Grid(format!(
                        "mapped horizon '{}' is fully undefined (all NaN)",
                        stack.horizons[b].name
                    )));
                }
                let solved = warm_surface(nx, ny, &controls_b)?
                    .surface()
                    .taper_beyond_data(&controls_b, policy);
                let above = surfaces[a].as_ref().expect("upper anchor built");
                (0..nx * ny)
                    .map(|idx| {
                        let (ip, jp) = (idx % nx, idx / nx);
                        (solved.z(ip, jp) - above.z(ip, jp)).max(0.0)
                    })
                    .collect()
            } else {
                let solved = warm_surface(nx, ny, &samples)?
                    .surface()
                    .taper_beyond_data(&samples, policy);
                (0..nx * ny)
                    .map(|idx| {
                        let (ip, jp) = (idx % nx, idx / nx);
                        solved.z(ip, jp).max(0.0)
                    })
                    .collect()
            }
        };
        // Monotone cumulative: ordering by construction; a genuine data crossing
        // collapses the zone toward the data (never adds rock).
        let cum_b: Vec<f64> = cum_a.iter().zip(&iso).map(|(p, i)| i.max(*p)).collect();
        surfaces[b] = Some(top.offset_by_field(&cum_b)?);

        // Tops-only splits inside envelope (a → b): min(pick iso, envelope iso).
        let above = surfaces[a].clone().expect("upper anchor built");
        for k in (a + 1)..b {
            let HorizonSource::TopsOnly(picks) = &stack.horizons[k].source else {
                unreachable!("a horizon between two mapped anchors must be tops-only");
            };
            let split = pick_thickness_field(picks, &above, nx, ny, policy)?;
            let cum_k: Vec<f64> = split
                .iter()
                .zip(cum_a.iter().zip(&cum_b))
                .map(|(s, (ca, cb))| ca + s.min(cb - ca))
                .collect();
            surfaces[k] = Some(top.offset_by_field(&cum_k)?);
        }
        cums[b] = Some(cum_b);
    }

    // Trailing tops-only horizons below the last mapped anchor: no envelope → the
    // absolute drape from the last mapped surface.
    let last_mapped = *mapped_idx
        .last()
        .expect("validate_stack requires a mapped top");
    for k in (last_mapped + 1)..n {
        let HorizonSource::TopsOnly(picks) = &stack.horizons[k].source else {
            unreachable!("horizons below the last mapped anchor are tops-only");
        };
        let above = surfaces[last_mapped]
            .clone()
            .expect("last mapped anchor built");
        surfaces[k] = Some(drape_tops_only(picks, &above, nx, ny, policy)?);
    }

    let surfaces: Vec<Surface> = surfaces
        .into_iter()
        .map(|s| s.expect("every horizon resolved"))
        .collect();
    Ok(surfaces)
}

/// Repair one stack interface, honouring **mapped-over-derived precedence**: a
/// *mapped* horizon is authoritative, so when a *derived* (tops-only) surface
/// crosses a mapped one the DERIVED surface yields — never the mapped one. When
/// both are mapped, or the lower is itself the derived surface, the lower yields
/// (the historical behaviour: pull `lower` to sit at/below `upper`). With a
/// `min_thickness` floor a thin/crossing node is pushed to the floor (recording an
/// [`InterfaceRepair`]); without it a crossing errors unless `clamp`.
///
/// Returns `(upper, lower)` with exactly one surface moved. With the isochore
/// construction ([`resolve_stack_surfaces`]) a consistent stack never crosses
/// here, so this is a safety net; it still fires under a `min_thickness` floor or
/// genuinely inconsistent inputs.
#[allow(clippy::too_many_arguments)]
pub(crate) fn repair_interface(
    upper: &Surface,
    lower: &Surface,
    min_thickness: Option<f64>,
    clamp: bool,
    upper_derived: bool,
    lower_derived: bool,
    interface: usize,
    repairs: &mut Vec<InterfaceRepair>,
) -> Result<(Surface, Surface), StaticError> {
    // The derived surface yields to the mapped one. Only when the upper is derived
    // AND the lower is mapped does the upper move; every other tie (both mapped,
    // both derived, or a derived lower) keeps the historical lower-yields behaviour.
    if upper_derived && !lower_derived {
        match min_thickness {
            Some(min_t) => {
                let (rep, columns, worst_m) =
                    upper.repair_min_thickness_from_below(lower, min_t)?;
                if columns > 0 {
                    repairs.push(InterfaceRepair {
                        interface,
                        columns,
                        worst_m,
                    });
                }
                Ok((rep, lower.clone()))
            }
            None => Ok((upper.guard_above(lower, clamp)?, lower.clone())),
        }
    } else {
        match min_thickness {
            Some(min_t) => {
                let (rep, columns, worst_m) = lower.repair_min_thickness(upper, min_t)?;
                if columns > 0 {
                    repairs.push(InterfaceRepair {
                        interface,
                        columns,
                        worst_m,
                    });
                }
                Ok((upper.clone(), rep))
            }
            None => Ok((upper.clone(), lower.guard_below(upper, clamp)?)),
        }
    }
}

/// Assemble the framework [`Wireframe`] from the resolved stack surfaces: horizon
/// `0` is `Top`, the last is `Base`, the rest `Intermediate`; every resolved node
/// is a control. `contacts` is the union of the per-zone contacts.
/// The world areal outline ring for a horizon-stack framework (finding 2): an
/// explicit `with_boundary` ring wins; else a **world-extent rectangle** derived from
/// the model georef + areal cell dims (`ni × nj`) — the grid's world cell-edge extent,
/// matching [`crate::view::GridFrame`] so the map outline overlays the raster; else
/// (no georef, a synthetic model) the local unit square. `nx`/`ny` are node dims
/// (`ni + 1`, `nj + 1`).
pub(crate) fn stack_boundary_ring(
    explicit: Option<&Vec<[f64; 2]>>,
    georef: Option<Georef>,
    nx: usize,
    ny: usize,
) -> Vec<[f64; 2]> {
    if let Some(ring) = explicit {
        return ring.clone();
    }
    if let Some(g) = georef {
        let (ni, nj) = ((nx - 1) as f64, (ny - 1) as f64);
        let xmin = g.origin_x - g.spacing_x / 2.0;
        let xmax = g.origin_x + (ni - 0.5) * g.spacing_x;
        let ymin = g.origin_y - g.spacing_y / 2.0;
        let ymax = g.origin_y + (nj - 0.5) * g.spacing_y;
        return vec![
            [xmin, ymin],
            [xmax, ymin],
            [xmax, ymax],
            [xmin, ymax],
            [xmin, ymin],
        ];
    }
    vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]
}

pub(crate) fn stack_framework(
    horizons: &[StackHorizon],
    surfaces: &[Surface],
    contacts: Vec<Contact>,
    nx: usize,
    ny: usize,
    ring: Vec<[f64; 2]>,
) -> Wireframe {
    let names: Vec<String> = horizons.iter().map(|h| h.name.clone()).collect();
    stack_framework_from_names(&names, surfaces, contacts, nx, ny, ring)
}

/// [`stack_framework`] keyed by horizon **names** only (roles derive from position:
/// first `Top`, last `Base`, else `Intermediate`) — so a caller holding just the
/// resolved surfaces + names (the MC template's well-tie re-solve) can rebuild the
/// framework without reconstructing the source [`StackHorizon`]s.
pub(crate) fn stack_framework_from_names(
    names: &[String],
    surfaces: &[Surface],
    contacts: Vec<Contact>,
    nx: usize,
    ny: usize,
    ring: Vec<[f64; 2]>,
) -> Wireframe {
    let n = names.len();
    let hs: Vec<Horizon> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let role = if i == 0 {
                HorizonRole::Top
            } else if i == n - 1 {
                HorizonRole::Base
            } else {
                HorizonRole::Intermediate
            };
            let mut depth = vec![0.0; nx * ny];
            for jp in 0..ny {
                for ip in 0..nx {
                    depth[jp * nx + ip] = surfaces[i].z(ip, jp);
                }
            }
            Horizon {
                name: name.clone(),
                role,
                surface: GriddedDepth {
                    ncol: nx,
                    nrow: ny,
                    depth_m: depth,
                    is_control: vec![true; nx * ny],
                },
            }
        })
        .collect();
    Wireframe {
        boundary: srs_wireframe::Boundary {
            ring,
            hardness: Hardness::Interpolated,
        },
        horizons: std::sync::Arc::new(hs),
        contacts,
    }
}

/// Converge a surface from control points **in petekTools kernel space** — a flat
/// bootstrap at the mean control depth (a fixed point of the kernel) followed by
/// one seeded solve. This is the *single* structural-solve path both the
/// deterministic [`StaticModelBuilder`] and the MC [`crate::StaticModelTemplate`]
/// use, so `template.realize(draw with gross == mean(g))` reproduces the
/// deterministic build on **all** column shapes — not only fully-pinned lattices
/// (graph `decision_gridder_kernel_unification`; R2). The cold
/// [`srs_gridder::solve_surface`] stays the accuracy reference but is no longer on
/// the build path, so the two never diverge by kernel.
///
/// # Errors
/// [`StaticError`] if there are no controls, the lattice is degenerate, or the
/// seeded solve fails.
pub(crate) fn warm_surface(
    nx: usize,
    ny: usize,
    controls: &[Control],
) -> Result<KernelSurface, StaticError> {
    if controls.is_empty() {
        return Err(StaticError::InvalidInput(
            "structural solve needs at least one control point".into(),
        ));
    }
    // Delegate to the converged solve (structure-fidelity audit S2): plane
    // detrending kills the kernel's slow affine mode (a flat bootstrap seed
    // otherwise stalls metres short of the fixed point on sparse-control
    // lattices) and fixed-point restarts verify the result. Kernel space is
    // preserved — the entry composes `solve_surface_seeded` with exact
    // fixed-point (plane) arithmetic.
    srs_gridder::solve_surface_converged(nx, ny, controls)
}

/// The structural skeleton a builder/template extracts off a wireframe.
pub(crate) struct WireframeStructure {
    /// Top-surface lattice dims (`ncol-1`, `nrow-1`).
    pub ni: usize,
    pub nj: usize,
    /// Control set from the `Top` horizon (every defined node = a hard datum).
    pub top_controls: Vec<Control>,
    /// Control set from a supplied `Base` horizon's defined nodes, if usable
    /// (matching lattice + at least one defined node). `None` = no usable Base.
    pub base_controls: Option<Vec<Control>>,
    /// Advisories for supplied horizons the build could not consume.
    pub warnings: Vec<BuildWarning>,
}

/// Turn a `GriddedDepth` surface into the control set of its defined nodes.
pub(crate) fn controls_from_surface(s: &GriddedDepth) -> Vec<Control> {
    let mut controls = Vec::new();
    for r in 0..s.nrow {
        for c in 0..s.ncol {
            let z = s.depth_m[r * s.ncol + c];
            if !z.is_nan() {
                controls.push(Control { ip: c, jp: r, z });
            }
        }
    }
    controls
}

/// Extract the structural skeleton a builder/template needs off a wireframe:
/// lattice dims, the `Top` control set (every defined node is a hard datum), an
/// optional `Base` control set (its real relief), and the column-base contact
/// depth. Supplied horizons that cannot be consumed (Intermediate; a second or
/// lattice-mismatched Base) are reported as [`BuildWarning`]s, not errors.
pub(crate) fn wireframe_structure(
    wf: &Wireframe,
    opts: &BuildOpts,
) -> Result<WireframeStructure, StaticError> {
    let top = wf
        .horizons
        .iter()
        .find(|h| h.role == HorizonRole::Top)
        .ok_or_else(|| StaticError::Grid("wireframe has no Top horizon".into()))?;
    let s = &top.surface;
    if s.ncol < 2 || s.nrow < 2 {
        return Err(StaticError::Grid(
            "top surface needs at least 2x2 nodes".into(),
        ));
    }
    if opts.nk == 0 {
        return Err(StaticError::Grid("nk must be non-zero".into()));
    }
    validate_positive("area_m2", opts.area_m2)?;
    validate_positive("gross_height_m", opts.gross_height_m)?;
    let top_controls = controls_from_surface(s);
    if top_controls.is_empty() {
        return Err(StaticError::Grid(
            "top surface is fully undefined (all NaN)".into(),
        ));
    }

    // The first Base horizon whose lattice matches the Top and carries at least
    // one defined node drives the pillar bases; anything else is an advisory.
    let mut base_controls = None;
    let mut warnings = Vec::new();
    for h in wf.horizons.iter() {
        match h.role {
            HorizonRole::Top => {}
            HorizonRole::Base => {
                let bs = &h.surface;
                let reason = if bs.ncol != s.ncol || bs.nrow != s.nrow {
                    Some(format!(
                        "Base lattice {}x{} does not match Top {}x{}",
                        bs.ncol, bs.nrow, s.ncol, s.nrow
                    ))
                } else if base_controls.is_some() {
                    Some(
                        "a Base horizon is already in use (only the first is honoured)".to_string(),
                    )
                } else {
                    let bc = controls_from_surface(bs);
                    if bc.is_empty() {
                        Some("Base surface is fully undefined (all NaN)".to_string())
                    } else {
                        base_controls = Some(bc);
                        None
                    }
                };
                if let Some(reason) = reason {
                    warnings.push(BuildWarning::UnusedHorizon {
                        name: h.name.clone(),
                        role: h.role,
                        reason,
                    });
                }
            }
            HorizonRole::Intermediate => warnings.push(BuildWarning::UnusedHorizon {
                name: h.name.clone(),
                role: h.role,
                reason: "intermediate horizons are unused until multi-zone layering (P5)"
                    .to_string(),
            }),
        }
    }

    // Assert the wireframe carries a fluid contact — the produced model's
    // `in_place` needs one on the framework to clip the hydrocarbon column.
    if wf.contacts.is_empty() {
        return Err(StaticError::Grid("wireframe has no fluid contact".into()));
    }
    Ok(WireframeStructure {
        ni: s.ncol - 1,
        nj: s.nrow - 1,
        top_controls,
        base_controls,
        warnings,
    })
}

/// A synthesized framework for the flat-box start: one flat `Top` horizon (every
/// node an assumed control) + one assumed OWC.
fn synth_framework(ncol: usize, nrow: usize, top_depth_m: f64, contact_depth_m: f64) -> Wireframe {
    Wireframe {
        boundary: srs_wireframe::Boundary {
            ring: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
            hardness: Hardness::Assumed,
        },
        horizons: std::sync::Arc::new(vec![Horizon {
            name: "Top".to_string(),
            role: HorizonRole::Top,
            surface: srs_wireframe::GriddedDepth {
                ncol,
                nrow,
                depth_m: vec![top_depth_m; ncol * nrow],
                is_control: vec![false; ncol * nrow],
            },
        }]),
        contacts: vec![Contact {
            kind: ContactKind::Owc,
            depth_m: contact_depth_m,
            hardness: Hardness::Assumed,
        }],
    }
}

#[cfg(test)]
mod conformity_tests {
    //! End-to-end layering-conformity tests through the builder + template:
    //! volume/pore conservation across styles, property population on truncated
    //! columns, the two-contact split on truncated columns, MC reproducibility,
    //! and NaN-marking of inactive layers in the section bundle.
    use super::*;
    use crate::draw::RealizationDraw;
    use crate::{
        Gaussian, PropertyPipeline, SectionSpec, StaticModelTemplate, UpscaleMethod, WellLog,
    };
    use petektools::{Variogram, VariogramModel};
    use srs_gridder::Conformity;
    use srs_wireframe::{
        Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
    };

    const DZ: f64 = 20.0; // Follow-style layer thickness → derived nk = 200/20 = 10.

    /// A wedge: flat top at 5000, base dipping only in i so node `ip` sits at
    /// `5000 + 20 + 180·(ip/(n-1))` (thickness 20 updip → 200 downdip). `contacts`
    /// is the fluid-contact stack the framework carries.
    fn wedge_wf(n: usize, contacts: Vec<Contact>) -> Wireframe {
        let top = vec![5000.0; n * n];
        let mut base = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                base[r * n + c] = 5000.0 + 20.0 + 180.0 * (c as f64 / (n as f64 - 1.0));
            }
        }
        let surf = |depth_m: Vec<f64>| GriddedDepth {
            ncol: n,
            nrow: n,
            depth_m,
            is_control: vec![true; n * n],
        };
        Wireframe {
            boundary: Boundary {
                ring: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
                hardness: Hardness::Hard,
            },
            horizons: std::sync::Arc::new(vec![
                Horizon {
                    name: "TopRes".into(),
                    role: HorizonRole::Top,
                    surface: surf(top),
                },
                Horizon {
                    name: "BaseRes".into(),
                    role: HorizonRole::Base,
                    surface: surf(base),
                },
            ]),
            contacts,
        }
    }

    fn owc(depth_m: f64) -> Vec<Contact> {
        vec![Contact {
            kind: ContactKind::Owc,
            depth_m,
            hardness: Hardness::Hard,
        }]
    }

    fn bopts(conformity: Conformity) -> BuildOpts {
        BuildOpts {
            area_m2: 90_000.0, // side 300 → dx = 30 over 10 columns
            gross_height_m: 100.0,
            nk: 8, // honoured only by Proportional
            conformity,
            solve_opts: SolveOpts::default(),
            priors: ConstantPriors {
                porosity: 0.25,
                net_to_gross: 0.8,
                water_saturation: 0.3,
            },
        }
    }

    fn truncated_count(m: &StaticModel) -> usize {
        m.provenance()
            .warnings
            .iter()
            .find_map(|w| match w {
                BuildWarning::LayersTruncated { cells } => Some(*cells),
                _ => None,
            })
            .unwrap_or(0)
    }

    #[test]
    fn follow_styles_conserve_grv_and_hcpv_vs_proportional() {
        // Whole column (contact far below the base) + constant priors → HCPV is
        // bulk × const, so both GRV and HCPV are conformity-invariant to FP.
        let ip = |c| {
            StaticModelBuilder::from_wireframe(&wedge_wf(11, owc(9000.0)), bopts(c))
                .unwrap()
                .build()
                .unwrap()
                .in_place()
                .unwrap()
        };
        let prop = ip(Conformity::Proportional);
        let ftop = ip(Conformity::FollowTop { dz_m: DZ });
        let fbase = ip(Conformity::FollowBase { dz_m: DZ });
        for (name, r) in [("FollowTop", &ftop), ("FollowBase", &fbase)] {
            assert!(
                (r.grv_m3 - prop.grv_m3).abs() / prop.grv_m3 < 1e-9,
                "{name} GRV {} != proportional {}",
                r.grv_m3,
                prop.grv_m3
            );
            assert!(
                (r.hcpv_m3 - prop.hcpv_m3).abs() / prop.hcpv_m3 < 1e-9,
                "{name} HCPV {} != proportional {}",
                r.hcpv_m3,
                prop.hcpv_m3
            );
        }
    }

    #[test]
    fn follow_top_truncates_and_populates_active_cells_without_nan() {
        // priors (NTG) + upscaled logs are irrelevant here; drive an SGS pipeline
        // on PORO so upscale + SGS both run over a truncated grid, then assert no
        // active cell holds NaN in any cube.
        let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 120.0).unwrap();
        let well_lo = WellLog::new(15.0, 15.0, vec![(5010.0, 0.18), (5060.0, 0.20)]);
        let well_hi = WellLog::new(285.0, 285.0, vec![(5010.0, 0.24), (5180.0, 0.26)]);
        let pipe = PropertyPipeline::new(srs_volumetrics::PORO)
            .upscale(vec![well_lo, well_hi], UpscaleMethod::Arithmetic)
            // A truncated wedge legitimately has simulated layers with no
            // conditioning data — opt into the structureless mean-fill rather than
            // the (default) hard error (item 4).
            .propagate(Gaussian::new(vgm, 7).allow_mean_fill());
        let m = StaticModelBuilder::from_wireframe(
            &wedge_wf(11, owc(9000.0)),
            bopts(Conformity::FollowTop { dz_m: DZ }),
        )
        .unwrap()
        .with_property(pipe)
        .build()
        .unwrap();

        assert_eq!(m.provenance().nk, 10, "dz-derived nk = ceil(200/20)");
        assert!(truncated_count(&m) > 0, "thin updip columns truncate");

        let grid = m.grid();
        let dims = grid.dims();
        let (ni, nj) = (dims.ni, dims.nj);
        let cube = |n: &str| m.property(n).unwrap().values.as_slice();
        let (poro, ntg, sw) = (
            cube(srs_volumetrics::PORO),
            cube(srs_volumetrics::NTG),
            cube(srs_volumetrics::SW),
        );
        let mut active = 0usize;
        for c in dims.iter() {
            if grid.cell(c).dz() <= 1e-9 {
                continue; // inactive/truncated — finite values here are irrelevant
            }
            active += 1;
            let idx = (c.k * nj + c.j) * ni + c.i;
            assert!(poro[idx].is_finite(), "PORO NaN in active cell {c:?}");
            assert!(ntg[idx].is_finite(), "NTG NaN in active cell {c:?}");
            assert!(sw[idx].is_finite(), "SW NaN in active cell {c:?}");
        }
        assert!(active > 0);
    }

    #[test]
    fn two_contact_split_on_truncated_columns() {
        // GOC + OWC inside the column depth range, under FollowTop truncation. The
        // split must be clean (finite, gas+oil == total, no NaN).
        let contacts = vec![
            Contact {
                kind: ContactKind::Goc,
                depth_m: 5050.0,
                hardness: Hardness::Hard,
            },
            Contact {
                kind: ContactKind::Owc,
                depth_m: 5120.0,
                hardness: Hardness::Hard,
            },
        ];
        let m = StaticModelBuilder::from_wireframe(
            &wedge_wf(11, contacts),
            bopts(Conformity::FollowTop { dz_m: DZ }),
        )
        .unwrap()
        .build()
        .unwrap();
        assert!(truncated_count(&m) > 0);
        let ip = m.in_place().unwrap();
        let gas = ip.gas.expect("gas cap");
        let oil = ip.oil.expect("oil leg");
        assert!(gas.hcpv_m3.is_finite() && oil.hcpv_m3.is_finite());
        assert!(gas.grv_m3 > 0.0 && oil.grv_m3 > 0.0);
        assert!(
            (ip.hcpv_m3 - (gas.hcpv_m3 + oil.hcpv_m3)).abs() < 1e-6,
            "total HCPV = gas + oil"
        );
        assert!(ip.per_cell_hcpv.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn template_realize_is_deterministic_under_follow_styles() {
        // Two fresh templates + the same draw → identical grids under each Follow
        // style (the dz-derived nk + truncation are geometry-driven, so the MC
        // regeneration seam stays bit-reproducible).
        let draw = RealizationDraw::new(90_000.0, 100.0, 9000.0, 0.25, 0.8, 0.3, 5);
        for style in [
            Conformity::FollowTop { dz_m: DZ },
            Conformity::FollowBase { dz_m: DZ },
        ] {
            let wf = wedge_wf(11, owc(9000.0));
            let mut t1 = StaticModelTemplate::new(&wf, bopts(style)).unwrap();
            let mut t2 = StaticModelTemplate::new(&wf, bopts(style)).unwrap();
            let a = t1.realize(&draw).unwrap();
            let b = t2.realize(&draw).unwrap();
            assert_eq!(a.provenance().nk, 10, "dz-derived nk");
            assert_eq!(a.grid().bulk_volume(), b.grid().bulk_volume());
            assert_eq!(truncated_count(&a), truncated_count(&b));
            assert!(
                truncated_count(&a) > 0,
                "{style:?} truncates the thin columns"
            );
            let (ia, ib) = (a.in_place().unwrap(), b.in_place().unwrap());
            assert_eq!(ia.grv_m3, ib.grv_m3);
            assert_eq!(ia.hcpv_m3, ib.hcpv_m3);
        }
    }

    #[test]
    fn section_bundle_nan_marks_inactive_layers() {
        // The updip thin columns truncate; their deep layers must be NaN in the
        // nk-sized section arrays while active layers stay finite.
        let m = StaticModelBuilder::from_wireframe(
            &wedge_wf(11, owc(9000.0)),
            bopts(Conformity::FollowTop { dz_m: DZ }),
        )
        .unwrap()
        .build()
        .unwrap();
        let nk = m.provenance().nk;
        // Fence along i at the first row's column centres (world = local here).
        let spec = SectionSpec::Polyline(vec![[15.0, 15.0], [285.0, 15.0]]);
        let b = m
            .intersection_bundle(&spec, Some(srs_volumetrics::PORO))
            .unwrap();
        assert!(!b.columns.is_empty());
        let mut saw_nan = false;
        for col in &b.columns {
            assert_eq!(col.layer_tops.len(), nk, "arrays stay nk-sized");
            assert_eq!(col.layer_bases.len(), nk);
            assert_eq!(col.values.len(), nk);
            for k in 0..nk {
                let t = col.layer_tops[k];
                if t.is_nan() {
                    saw_nan = true;
                    assert!(col.layer_bases[k].is_nan(), "base NaN with top");
                    assert!(col.values[k].is_nan(), "value NaN with geometry");
                } else {
                    // An active layer has finite, ordered geometry + property.
                    assert!(col.layer_bases[k].is_finite());
                    assert!(col.layer_bases[k] >= t);
                    assert!(col.values[k].is_finite(), "active PORO finite");
                }
            }
        }
        assert!(saw_nan, "the truncated updip columns must emit NaN layers");
    }

    #[test]
    fn section_bundle_inactive_depths_serialize_to_json_null_not_zero() {
        // The cross-codebase seam contract (SPEC §7e): an inactive/truncated layer's
        // depth is `f64::NAN`, and serde serializes NaN to JSON `null` — NOT `0`.
        // The viewer must frame its depth axis with a null-rejecting guard
        // (`v != null && isFinite(v)`; `isFinite(null) === true` is the footgun this
        // guard defends against). This test pins petekStatic's half of the seam: no
        // inactive depth ever reaches the viewer as a numeric 0 that would drag the
        // frame to depth 0. Companion fix (viewer.js `isFinite` guards + null-safe
        // horizon_traces rendering) routed to peteksim 2026-07-04.
        let m = StaticModelBuilder::from_wireframe(
            &wedge_wf(11, owc(9000.0)),
            bopts(Conformity::FollowTop { dz_m: DZ }),
        )
        .unwrap()
        .build()
        .unwrap();
        let spec = SectionSpec::Polyline(vec![[15.0, 15.0], [285.0, 15.0]]);
        let b = m
            .intersection_bundle(&spec, Some(srs_volumetrics::PORO))
            .unwrap();
        let v = serde_json::to_value(&b).unwrap();
        let cols = v.get("columns").and_then(|c| c.as_array()).unwrap();
        let mut saw_null = false;
        for (ci, col) in b.columns.iter().enumerate() {
            let jtops = cols[ci]
                .get("layer_tops")
                .and_then(|t| t.as_array())
                .unwrap();
            for (k, &t) in col.layer_tops.iter().enumerate() {
                if t.is_nan() {
                    saw_null = true;
                    assert!(
                        jtops[k].is_null(),
                        "inactive layer top must serialize to null"
                    );
                } else {
                    // A finite depth serializes to a JSON number — never null.
                    assert!(
                        jtops[k].is_number(),
                        "active layer top must be a JSON number"
                    );
                }
            }
        }
        assert!(
            saw_null,
            "truncated columns must serialize inactive depths as null"
        );
        // Belt-and-braces: horizon_traces NaN gaps serialize to null too (same pattern).
        if let Some(traces) = v.get("horizon_traces").and_then(|t| t.as_array()) {
            for tr in traces {
                if let Some(depths) = tr.get("depths").and_then(|d| d.as_array()) {
                    for d in depths {
                        assert!(
                            d.is_null() || d.is_number(),
                            "a horizon-trace depth is either a number or null (never 0-for-gap)"
                        );
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod repair_precedence_tests {
    //! Repair precedence: a DERIVED surface yields to a MAPPED one when they cross,
    //! never the reverse; ties (both mapped, both derived) keep the historical
    //! lower-yields behaviour. Exercises the crossing path directly (the isochore
    //! construction keeps it a no-op end-to-end, so it is unit-tested here).
    use super::*;
    use srs_gridder::Surface;

    fn surf(z: Vec<f64>) -> Surface {
        // 2x2 lattice built from its four corner controls.
        warm_surface(
            2,
            2,
            &[
                Control {
                    ip: 0,
                    jp: 0,
                    z: z[0],
                },
                Control {
                    ip: 1,
                    jp: 0,
                    z: z[1],
                },
                Control {
                    ip: 0,
                    jp: 1,
                    z: z[2],
                },
                Control {
                    ip: 1,
                    jp: 1,
                    z: z[3],
                },
            ],
        )
        .unwrap()
        .surface()
        .clone()
    }

    #[test]
    fn derived_upper_yields_to_mapped_lower() {
        // Upper (derived) crosses below the mapped lower at node 0 (12 > 10).
        let upper = surf(vec![12.0, 5.0, 5.0, 5.0]);
        let lower = surf(vec![10.0; 4]);
        let mut reps = Vec::new();
        // upper_derived = true, lower_derived = false → the UPPER moves up.
        let (u, l) =
            repair_interface(&upper, &lower, None, true, true, false, 0, &mut reps).unwrap();
        assert_eq!(l.z(0, 0), 10.0, "mapped lower is authoritative — untouched");
        assert_eq!(
            u.z(0, 0),
            10.0,
            "derived upper clamped up to the mapped lower"
        );
        assert_eq!(u.z(1, 0), 5.0, "non-crossing node untouched");
    }

    #[test]
    fn mapped_lower_yields_when_both_mapped() {
        // Same crossing, but both surfaces are mapped → the historical behaviour:
        // the LOWER yields (pulled down to the upper).
        let upper = surf(vec![12.0, 5.0, 5.0, 5.0]);
        let lower = surf(vec![10.0; 4]);
        let mut reps = Vec::new();
        let (u, l) =
            repair_interface(&upper, &lower, None, true, false, false, 0, &mut reps).unwrap();
        assert_eq!(u.z(0, 0), 12.0, "upper untouched when both mapped");
        assert_eq!(
            l.z(0, 0),
            12.0,
            "lower clamped down to the upper (lower yields)"
        );
    }
}

#[cfg(test)]
mod scatter_conditioner_tests {
    //! The factor-once / solve-many scatter conditioning lever
    //! (`task_suite_scatter_perf`): the sample (x,y) geometry is fixed across depth
    //! draws, so re-seating it with new depths must be a *resolve* on the reused
    //! factorization that is **bit-identical** to conditioning the perturbed scatter
    //! from scratch. Frame-referenced with a fictional world georef (R1).
    use super::*;

    fn world_frame() -> StackFrame {
        // 12x12 cells → 13x13 nodes; fictional world origin (R1 frame rule).
        StackFrame {
            ni: 12,
            nj: 12,
            georef: Georef::new(431_000.0 + 0.5 * 50.0, 6_521_000.0 + 0.5 * 50.0, 50.0, 50.0)
                .unwrap(),
        }
    }

    /// Dense off-node scatter at fixed world (x,y), depths from `depth_at`.
    fn scatter(frame: &StackFrame, depth_at: impl Fn(f64, f64) -> f64) -> Vec<WorldPoint> {
        let g = &frame.georef;
        let mut pts = Vec::new();
        // Off-node fractional positions across the interior (the +0.37/+0.61 keep
        // them off every node so the Bilinear data-fit governs).
        for j in 0..11 {
            for i in 0..11 {
                let fi = i as f64 + 0.37;
                let fj = j as f64 + 0.61;
                let x = g.origin_x + (fi - 0.5) * g.spacing_x;
                let y = g.origin_y + (fj - 0.5) * g.spacing_y;
                pts.push(WorldPoint {
                    x,
                    y,
                    depth_m: depth_at(fi, fj),
                });
            }
        }
        pts
    }

    #[test]
    fn resolve_on_a_reused_factor_is_bit_identical_to_conditioning_from_scratch() {
        let frame = world_frame();
        // Two depth fields over the SAME sample geometry (the MC lever: fixed (x,y),
        // only depths change) — a plane and a curved perturbation of it.
        let pts_a = scatter(&frame, |fi, fj| 2000.0 + 3.0 * fi - 1.5 * fj);
        let pts_b = scatter(&frame, |fi, fj| {
            2000.0 + 3.0 * fi - 1.5 * fj + 7.0 * (0.3 * fi).sin() - 4.0 * (0.2 * fj).cos()
        });

        // Factor ONCE from A's geometry; both depth vectors share that factorization.
        let (cond, depths_a) = ScatterConditioner::factor(&pts_a, &frame).expect("factor A");
        let (_, depths_b) = ScatterConditioner::factor(&pts_b, &frame).expect("factor B");

        let reused_a = cond.resolve(&depths_a).expect("resolve A");
        let reused_b = cond.resolve(&depths_b).expect("resolve B");

        // Independent conditioning of each (a fresh factor+solve per depth field).
        let fresh_a = grid_scatter(&pts_a, &frame);
        let fresh_b = grid_scatter(&pts_b, &frame);

        // Bit-identity (`NaN` voids must match bit-for-bit — `to_bits` so `NaN == NaN`).
        let bit_identical = |x: &[f64], y: &[f64]| {
            x.len() == y.len() && x.iter().zip(y).all(|(a, b)| a.to_bits() == b.to_bits())
        };
        assert!(
            bit_identical(&reused_a, &fresh_a.depth_m),
            "resolve on the reused factor must be bit-identical to conditioning A from scratch"
        );
        assert!(
            bit_identical(&reused_b, &fresh_b.depth_m),
            "resolve on the reused factor with B's depths must equal conditioning B from scratch \
             (the fixed geometry, varying-depth MC lever)"
        );
        // Sanity: the two depth fields actually differ where the support hull covers.
        assert!(
            reused_a
                .iter()
                .zip(&reused_b)
                .any(|(a, b)| a.is_finite() && b.is_finite() && (a - b).abs() > 1e-6),
            "the two draws must produce different conditioned surfaces"
        );
    }
}
