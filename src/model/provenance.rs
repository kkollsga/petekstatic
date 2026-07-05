//! Provenance — the reproducibility record (SPEC §6): what inputs + parameters
//! produced this model. Metadata only — it never affects geometry or cubes, but
//! it lets a consumer label a realization and answer "which inputs gave me this
//! P10 model?".

use crate::gridder::{Conformity, SolveOpts};
use crate::model::draw::RealizationDraw;
use crate::model::pipeline::PropertyReport;
use crate::wireframe::HorizonRole;

/// How the property cubes were populated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PopulationMode {
    /// Constant priors everywhere (the no-data Day-1 case).
    Priors,
    /// Upscaled from positioned petro samples (cells outside log coverage keep
    /// the priors).
    Logs,
}

/// A non-blocking advisory raised during a build: something the caller supplied
/// that the current build could not (yet) honour. Carried on [`Provenance`] so a
/// consumer can surface it without the build failing.
#[derive(Debug, Clone, PartialEq)]
pub enum BuildWarning {
    /// A supplied horizon was present in the wireframe but not consumed by the
    /// build. `Intermediate` horizons are unused until multi-zone layering lands
    /// (P5 `task_petekstatic_zones_faults`); a second `Base`, or a `Base` whose
    /// lattice does not match the `Top`, is likewise skipped (the build falls
    /// back to the constant `gross_height_m` offset for the base).
    UnusedHorizon {
        /// The horizon's `name`.
        name: String,
        /// The horizon's structural role.
        role: HorizonRole,
        /// Why it was not consumed (human-readable).
        reason: String,
    },
    /// Post-gridding order-repair pulled thin/crossing base columns down to a
    /// minimum thickness below the top (`with_min_thickness_m`). The top (the
    /// better-constrained seismic pick) was preserved; only the base moved. Raised
    /// instead of the [`crate::error::StaticError::CrossedSurfaces`] error when the caller
    /// opted into the repair.
    ThinColumnsRepaired {
        /// How many lattice nodes were pushed down to the minimum thickness.
        columns: usize,
        /// The worst (most negative) original `base − top` separation among the
        /// repaired nodes — negative = a true crossing.
        worst_m: f64,
    },
    /// `cells` cells collapsed to zero thickness during layering (zero volume —
    /// excluded from volumetrics, `NaN`-marked in the view bundles). Two causes,
    /// **under any conformity**: a dz-based Follow style (`FollowTop`/`FollowBase`)
    /// truncating the thinner columns against the pinch-out horizon, or a
    /// zero-thickness column where the bounding horizons pinch out / merge (a
    /// Proportional zone over a merged envelope reports here too — the warning is
    /// NOT Follow-specific; see `StackProvenance::zones[*].truncated_cells` +
    /// `conformity` for the per-zone attribution). Informational: total volume is
    /// conformity-invariant, so this never changes the in-place answer.
    LayersTruncated {
        /// Number of cells collapsed to zero thickness by truncation.
        cells: usize,
    },
    /// The sub-threshold **cell-collapse** pass (`with_collapse_below_m`) collapsed
    /// `cells` cells thinner than the threshold to zero thickness, merging each
    /// sliver's rock into a thicker zone-interior neighbour (volume-conserving).
    /// Informational: total volume is preserved, so this never changes the in-place
    /// answer; the collapsed cells are `NaN`-marked in the view bundles.
    CellsCollapsed {
        /// Number of cells collapsed to zero thickness by the threshold.
        cells: usize,
    },
    /// The dz-derived layer count under a Follow conformity hit the `MAX_NK` cap
    /// (the chosen dz is finer than the cap can span over the thickest column), so
    /// the deepest part of the thickest columns is not layered. Coarsen dz or
    /// accept the cap. Carries the capped count.
    LayerCountCapped {
        /// The capped k-layer count (`crate::gridder::MAX_NK`).
        nk: usize,
    },
}

/// One zone's provenance in a multi-zone horizon-stack build: its name, bounding
/// horizons, layering style, effective layer count + first global `k`, and the
/// cells it truncated. Mirrors the gridder's `StackedZone` with the geological
/// names attached.
#[derive(Debug, Clone, PartialEq)]
pub struct ZoneProvenance {
    /// The zone's name.
    pub name: String,
    /// The bounding horizon above (top) and below (base).
    pub top_horizon: String,
    /// The bounding horizon below (base).
    pub base_horizon: String,
    /// This zone's conformity/layering style.
    pub conformity: Conformity,
    /// This zone's effective k-layer count.
    pub nk: usize,
    /// This zone's first global k-layer.
    pub k_start: usize,
    /// Cells in this zone collapsed to zero thickness by truncation.
    pub truncated_cells: usize,
}

/// A per-interface order-repair record from a horizon-stack build: consecutive
/// horizons `interface` and `interface + 1` crossed (or came within the minimum
/// thickness), so the lower horizon was pulled down at `columns` nodes; `worst_m`
/// is the worst (most-negative) original separation among them.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceRepair {
    /// Index of the upper horizon of the repaired interface (repairs pair
    /// `interface` with `interface + 1`).
    pub interface: usize,
    /// How many lattice nodes were pushed down to the minimum thickness.
    pub columns: usize,
    /// The worst (most-negative) original `lower − upper` separation among the
    /// repaired nodes — negative = a true crossing.
    pub worst_m: f64,
}

/// The record of a multi-zone horizon-stack build: the ordered horizon names
/// (top→down), the per-zone layering, and any per-interface order-repairs. `None`
/// on [`Provenance::stack`] for the 2-surface (Top+Base) degenerate path.
#[derive(Debug, Clone, PartialEq)]
pub struct StackProvenance {
    /// Ordered framework horizon names, top→down (`N` of them → `N − 1` zones).
    pub horizons: Vec<String>,
    /// Per-zone layering, top→base; k-ranges partition `[0, nk)`.
    pub zones: Vec<ZoneProvenance>,
    /// Per-interface order-repairs applied to keep consecutive horizons ordered.
    pub interface_repairs: Vec<InterfaceRepair>,
}

/// One horizon's tie residual at a well (P8 per-horizon ties,
/// `task_petekstatic_multizone_2`): the measured formation top, the model surface
/// depth at the well node, and their difference — the **pre-tie** mismatch the tie
/// resolved (`residual_m = measured_depth_m − model_depth_m`, positive = the well is
/// deeper than the untied model surface).
#[derive(Debug, Clone, PartialEq)]
pub struct HorizonTieResidual {
    /// The framework horizon this tie is on.
    pub horizon: String,
    /// The measured formation top at the well \[m, positive-down\].
    pub measured_depth_m: f64,
    /// The **untied** model surface depth at the well node \[m\].
    pub model_depth_m: f64,
    /// `measured_depth_m − model_depth_m` \[m\] — the tie shift / QC residual.
    pub residual_m: f64,
}

/// One well's tie record: its id, world marker position, control node, and the
/// per-horizon tie residuals ([`HorizonTieResidual`]). Carried on
/// [`Provenance::well_ties`] and surfaced in the map bundle's `wells`.
#[derive(Debug, Clone, PartialEq)]
pub struct WellTieRecord {
    /// Well identifier.
    pub id: String,
    /// World easting of the well marker.
    pub x: f64,
    /// World northing of the well marker.
    pub y: f64,
    /// Control-lattice node the ties pinned.
    pub ip: usize,
    /// Control-lattice node the ties pinned.
    pub jp: usize,
    /// Per-horizon tie residuals for this well, in the order supplied.
    pub residuals: Vec<HorizonTieResidual>,
}

/// The record of what produced a [`crate::model::StaticModel`].
#[derive(Debug, Clone)]
pub struct Provenance {
    /// Identity of the input bundle (petekIO `.pproj` / `ModelInputs` reference,
    /// or a caller-supplied label).
    pub inputs_ref: String,
    /// Gridder settings the (cold) surface solve used. For a warm-started
    /// realization the kernel owns its own parameters; these are the template's
    /// nominal settings.
    pub solve_opts: SolveOpts,
    /// The layering scheme.
    pub conformity: Conformity,
    /// Number of k-layers.
    pub nk: usize,
    /// How the cubes were populated.
    pub population: PopulationMode,
    /// The per-realization draw, if this model is a Monte-Carlo realization
    /// (`None` for a deterministic single-shot build). Carries `seed_index`.
    pub realization: Option<RealizationDraw>,
    /// Non-blocking advisories raised during the build (e.g. supplied horizons
    /// the current build could not consume). Empty on a clean build.
    pub warnings: Vec<BuildWarning>,
    /// Per-property geostatistical-pipeline reports (upscale QC + whether SGS ran),
    /// one per attached [`crate::model::PropertyPipeline`]. Empty when no pipeline is
    /// attached (constant-prior / single-well-log population).
    pub property_reports: Vec<PropertyReport>,
    /// The multi-zone horizon-stack record (ordered horizons + per-zone layering +
    /// interface repairs) when the model was built from an ordered stack
    /// (`from_horizon_stack`). `None` for the 2-surface (Top+Base) path — the
    /// single implicit zone is described by [`Provenance::conformity`] +
    /// [`Provenance::nk`].
    pub stack: Option<StackProvenance>,
    /// Explicit per-horizon well ties applied to the build (`with_well_ties`), one
    /// [`WellTieRecord`] per tie well with its per-horizon residuals. Empty when no
    /// ties were supplied. Surfaced in the map bundle's `wells[].ties`.
    pub well_ties: Vec<WellTieRecord>,
    /// **Sugar-cube mode** (`with_sugar_cube`): render section cells as flat boxes
    /// rather than dip-following trapezoids. Default `false` — the engine geometry is
    /// corner-point, so the section bundle carries true per-edge cell depths; this
    /// only flattens the section view's edge arrays to the centroid trace. Flows into
    /// [`crate::model::view::IntersectionBundle::sugar_cube`].
    pub sugar_cube: bool,
}
