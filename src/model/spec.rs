//! [`BuildSpec`] — the ONE declarative build configuration consumed by **both**
//! [`crate::model::StaticModelBuilder`] and [`crate::model::StaticModelTemplate`]
//! (`task_petekstatic_spec_mirror`, suite api-consistency contract).
//!
//! Before this, the builder and the template each carried the same ~12 `with_*`
//! setters over duplicated private fields. Now both hold a `BuildSpec` internally;
//! every `with_*` is **thin sugar** mutating the spec, and `with_spec` installs a
//! whole configuration in one call. The values are identical either way, so the
//! determinism contracts (builder == template at mean gross, sharded == serial,
//! `realize_into` staleness, byte-identity) are untouched — pinned by
//! `tests/spec_conformance.rs`.
//!
//! `BuildSpec` is a **spec value** per the family pattern (modeling-api-v2): it
//! says WHAT, is immutable-valued (`with_*` return new values), serializes
//! (`serde`; a scenario is a savable file), and compares by value. Settings that
//! say HOW ride on it as sub-values ([`TieSettings`]).
//!
//! ## Forward-compat (structural uncertainty, `decision_structural_uncertainty_isochore`)
//! The owner-refined structural-uncertainty design (top-depth + per-zone ISOCHORE
//! perturbation) lands **additively**: `BuildSpec` is `#[non_exhaustive]` +
//! `#[serde(default)]`, and [`crate::model::RealizationDraw`] / [`crate::model::ZoneDraw`] are
//! `#[non_exhaustive]`, so a per-zone thickness-field leg is a new optional field,
//! not a reshape. Nothing here closes that door.

use crate::gridder::ExtrapolationPolicy;
use crate::model::builder::WellTie;
use crate::model::model::Georef;
use serde::{Deserialize, Serialize};

/// How a [`WellTie`]'s measured top is folded into its mapped horizon's gridded
/// datum — the settings mirror over the datum-substitution tie machinery.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub enum TieMethod {
    /// **Datum substitution** (the default — today's behaviour): the measured top
    /// REPLACES the map datum at the tie node (or defines a previously-undefined
    /// node), and the full stack resolution re-runs over the tied data. On a
    /// fully-defined lattice every other node is still a hard datum, so the tie
    /// moves exactly the tied node (radius of influence 0 cells); on a sparse
    /// lattice its influence is the solver's interpolation reach, bounded beyond
    /// the data hull by the [`ExtrapolationPolicy`].
    #[default]
    Replace,
    /// **Bounded locality**: the tie's residual (`measured − untied surface`) is
    /// blended into every *defined* map datum within `radius_m` of the tie node,
    /// with a measured linear decay — weight `1` at the well (the tie node lands
    /// exactly on the measured top) falling to `0` at `radius_m` — then the full
    /// stack resolution re-runs over the adjusted data. Beyond `radius_m` the map
    /// datums are bit-untouched. Distances are metres on the model's areal cell
    /// spacing (the builder's `area_m2`-derived spacing; the template's nominal
    /// `opts.area_m2` footprint — ties are draw-invariant, so a per-draw area does
    /// not re-evaluate the locality). `radius_m` must be finite and positive.
    Radius {
        /// The locality radius \[m\]: the residual decays linearly to zero here.
        radius_m: f64,
    },
}

/// Settings for how well ties are applied ([`TieMethod`]) — attach on the spec
/// via [`BuildSpec::with_tie_settings`]. `Default` = [`TieMethod::Replace`]
/// (today's behaviour).
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct TieSettings {
    /// How a measured top is folded into the mapped datum grid.
    pub method: TieMethod,
}

impl TieSettings {
    /// Datum substitution (the default): the tie replaces the map datum at its node.
    #[must_use]
    pub fn replace() -> Self {
        Self {
            method: TieMethod::Replace,
        }
    }

    /// Bounded locality: the tie's residual decays linearly to zero at `radius_m`.
    #[must_use]
    pub fn radius(radius_m: f64) -> Self {
        Self {
            method: TieMethod::Radius { radius_m },
        }
    }
}

/// The declarative build configuration shared by the deterministic builder and
/// the MC template (see the [module docs](self)). Construct with
/// [`BuildSpec::new`] + the `with_*` sugar (it is `#[non_exhaustive]`), install
/// with `with_spec` on either consumer — or keep using the consumers' own
/// `with_*` methods, which are the same sugar over their internal spec.
///
/// Not in the spec (deliberately): population attachments (`with_logs`,
/// `with_areal_trend`, `with_property*`, `with_zone_*` — pipelines/data, not
/// declarative config), the builder's memory/spill resources (`MemoryBudget`,
/// run resources belong to [`crate::model::McSettings`] per the suite ruling), and the
/// per-draw scalars (those ride [`crate::model::RealizationDraw`]).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
#[serde(default)]
pub struct BuildSpec {
    /// Provenance label for the input bundle. `None` = the entry point's default
    /// (`"flat-box"` / `"wireframe"` / `"horizon-stack"`).
    pub inputs_ref: Option<String>,
    /// The registered world georeference (see `with_georef` on the consumers).
    /// `None` = the local degenerate frame.
    pub georef: Option<Georef>,
    /// World areal boundary ring for a horizon-stack build (map-bundle outline).
    /// `None` → a world-extent rectangle from the georef, or the local unit square.
    pub boundary: Option<Vec<[f64; 2]>>,
    /// How a stack surface behaves beyond its data. Default
    /// [`ExtrapolationPolicy::DecayToData`].
    pub extrapolation: ExtrapolationPolicy,
    /// Clamp a crossed base to zero gross instead of erroring (R1). Default `false`.
    pub clamp_base_to_top: bool,
    /// Post-gridding order-repair floor \[m\] (R-c); `None` = the crossing guard.
    pub min_thickness_m: Option<f64>,
    /// Sub-threshold cell-collapse floor \[m\]; `None` = off.
    pub collapse_below_m: Option<f64>,
    /// Sugar-cube section rendering; default `false` (dip-following trapezoids).
    pub sugar_cube: bool,
    /// Gas-cap connate-water override for a two-contact `in_place` split (R3).
    /// Builder-level; the template's analog rides each draw (`RealizationDraw::sw_gas`).
    pub sw_gas: Option<f64>,
    /// Explicit per-horizon well ties (P8), applied per [`TieSettings`].
    pub well_ties: Vec<WellTie>,
    /// How the well ties are applied. Default [`TieMethod::Replace`].
    pub ties: TieSettings,
}

impl BuildSpec {
    /// The default spec — every knob at the consumers' historical default.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the provenance label for the input bundle.
    #[must_use]
    pub fn with_inputs_ref(mut self, inputs_ref: impl Into<String>) -> Self {
        self.inputs_ref = Some(inputs_ref.into());
        self
    }

    /// Register the world georeference (column `(0,0)` centroid + spacing).
    /// Non-finite / non-positive spacing is ignored (the local degenerate frame).
    #[must_use]
    pub fn with_georef(
        mut self,
        origin_x: f64,
        origin_y: f64,
        spacing_x: f64,
        spacing_y: f64,
    ) -> Self {
        self.georef = Georef::new(origin_x, origin_y, spacing_x, spacing_y);
        self
    }

    /// Register the world areal boundary ring (first == last point to close).
    /// An empty ring clears it.
    #[must_use]
    pub fn with_boundary(mut self, ring: Vec<[f64; 2]>) -> Self {
        self.boundary = if ring.is_empty() { None } else { Some(ring) };
        self
    }

    /// Set the extrapolation policy beyond the data hull.
    #[must_use]
    pub fn with_extrapolation(mut self, policy: ExtrapolationPolicy) -> Self {
        self.extrapolation = policy;
        self
    }

    /// Opt into clamping a crossed base to zero gross instead of erroring (R1).
    #[must_use]
    pub fn with_clamp_base_to_top(mut self, clamp: bool) -> Self {
        self.clamp_base_to_top = clamp;
        self
    }

    /// Opt into the post-gridding order-repair floor (R-c).
    #[must_use]
    pub fn with_min_thickness_m(mut self, min_thickness_m: f64) -> Self {
        self.min_thickness_m = Some(min_thickness_m);
        self
    }

    /// Opt into the sub-threshold cell-collapse floor.
    #[must_use]
    pub fn with_collapse_below_m(mut self, collapse_below_m: f64) -> Self {
        self.collapse_below_m = Some(collapse_below_m);
        self
    }

    /// Opt into sugar-cube section rendering (flat-box cells).
    #[must_use]
    pub fn with_sugar_cube(mut self, sugar_cube: bool) -> Self {
        self.sugar_cube = sugar_cube;
        self
    }

    /// Set the gas-cap connate-water override (R3; builder-level).
    #[must_use]
    pub fn with_sw_gas(mut self, sw_gas: f64) -> Self {
        self.sw_gas = Some(sw_gas);
        self
    }

    /// Attach the explicit per-horizon well ties (replaces any prior set).
    #[must_use]
    pub fn with_well_ties(mut self, ties: Vec<WellTie>) -> Self {
        self.well_ties = ties;
        self
    }

    /// Set how the well ties are applied ([`TieSettings`]).
    #[must_use]
    pub fn with_tie_settings(mut self, ties: TieSettings) -> Self {
        self.ties = ties;
        self
    }
}
