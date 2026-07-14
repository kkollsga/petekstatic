//! The per-property geostatistical population pipeline (P5
//! `task_petekstatic_property_modelling`; coordinator design
//! `python-model-build-api.md`).
//!
//! Properties are modelled **one at a time** through a visible, inspectable
//! pipeline:
//!
//! ```text
//! PropertyPipeline::new("PHIE")
//!     .upscale(wells, UpscaleMethod::Arithmetic)   // logs -> per-cell conditioned values (+ QC)
//!     .propagate(Gaussian::new(variogram, seed))   // SGS fills every cell, conditioned on the upscaled cells
//! ```
//!
//! - **upscale** ([`PropertyPipeline::upscale_cells`]) is a first-class step: each
//!   positioned [`WellLog`] snaps to its areal column, and its samples that fall in a
//!   cell's depth range are upscaled (`srs-petro` power means) into that cell — a
//!   per-cell conditioned field (`NaN` where no log passes) plus an [`UpscaleQc`]
//!   report (upscaled-vs-log statistics).
//! - **propagate** runs petekTools' sequential Gaussian simulation
//!   ([`petektools::geostat::sgs`]) **per k-layer**, conditioned exactly on that
//!   layer's upscaled cells, so every cell is filled with a seeded, reproducible
//!   draw that honours the wells and reproduces the data histogram. A collocated
//!   secondary ([`Gaussian::with_trend`], Phase 3) steers it via cokriging.
//!
//! This supersedes the interim [`crate::model::TrendSurface`] multiplier hook (which only
//! shaped a prior *laterally*): the pipeline produces the property field itself,
//! conditioned on real logs.
//!
//! ## Axis-aligned / regular lattice (documented limitation)
//! The areal SGS lattice is reconstructed from the grid's top-layer cell centroids
//! (origin + node spacing), so — like [`petektools::resample`] — it is designed and
//! tested for **axis-aligned, regular** column layouts (the `layer_grid` box /
//! conformable grids). A rotated or irregularly-warped pillar lattice is future
//! work.

use crate::error::StaticError;
use crate::grid::{Grid, Ijk, Property};
use crate::model::model::Georef;
use crate::model::trend::TrendSurface;
use crate::petro::{arithmetic_mean, geometric_mean, harmonic_mean, WeightedSample};
use petektools::geostat::{sgs_seeded, SgsParams};
use petektools::{Lattice, SpatialVariogram};
use rayon::prelude::*;

/// The golden-ratio odd constant used to derive an independent, reproducible SGS
/// seed per k-layer from the pipeline's base seed.
const SEED_GOLDEN: u64 = 0x9E37_79B9_7F4A_7C15;

/// Minimum fraction of areal nodes a collocated trend secondary must cover (be
/// finite on) after resampling to the model frame, or [`PropertyPipeline::propagate`]
/// errors instead of silently steering nothing.
///
/// The kernel drops a `NaN` secondary node to plain simple kriging **per node**
/// (`petektools::geostat::sgs`), so a trend whose georeference does not overlap the
/// model frame resamples to an all-`NaN` secondary and every node silently falls
/// back — collocated cokriging becomes a no-op indistinguishable from plain SGS
/// (`task_petekstatic_zoned_fixes` finding 1: a world-georeferenced trend resampled
/// onto the model's local lattice covered **0** nodes). Below this fraction that is
/// almost always a frame/georeference mismatch, not intent, so it is a hard error.
/// Above it, genuine partial coverage (a trend defined only over part of the area)
/// keeps the documented per-node SK fallback for the uncovered nodes. Half the grid
/// is the threshold: a trend meant to steer the field covers most of it; a
/// sub-half-coverage secondary cannot meaningfully condition the simulation.
const COLLOCATED_MIN_COVERAGE: f64 = 0.5;

/// A positioned well log for one property: the well's world `(x, y)` and its
/// `(depth_m, value)` samples down-hole (any order; binned per cell at upscale time).
#[derive(Debug, Clone, PartialEq)]
pub struct WellLog {
    /// World easting of the well.
    pub x: f64,
    /// World northing of the well.
    pub y: f64,
    /// `(depth_m, value)` samples — model-internal **positive-down depth in
    /// metres** (the same datum as the grid cells these bin against; any
    /// negative-down ingest flip is already applied upstream) + the property
    /// value at that depth.
    pub samples: Vec<(f64, f64)>,
}

impl WellLog {
    /// A positioned well log from its world `(x, y)` and `(depth_m, value)` samples.
    #[must_use]
    pub fn new(x: f64, y: f64, samples: Vec<(f64, f64)>) -> Self {
        Self { x, y, samples }
    }
}

/// How a property pipeline behaves across Monte-Carlo realizations
/// (`decision_mc_composition`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum McMode {
    /// **Default.** Propagate the field **once** (at the first realization) and
    /// reuse that pattern every draw, shifting only its **level** by the draw's
    /// per-property shift — cheap (~ms realizations), captures the level
    /// uncertainty that dominates this asset class.
    #[default]
    LevelShift,
    /// Re-run SGS with a **fresh per-draw seed** every realization — a new spatial
    /// pattern each draw (captures heterogeneity uncertainty), at the cost of a full
    /// simulation per draw.
    Resimulate,
}

/// How the log samples that fall in a cell's depth range are averaged into the
/// cell's conditioned value — each weighted by what the property conserves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UpscaleMethod {
    /// Length-weighted arithmetic mean (porosity, NTG — additive properties).
    Arithmetic,
    /// Length-weighted harmonic mean (vertical permeability — series flow).
    Harmonic,
    /// Length-weighted geometric mean (isotropic permeability).
    Geometric,
}

impl UpscaleMethod {
    /// Average a set of in-cell samples (unit weight each) by this method.
    fn apply(self, values: &[f64]) -> Result<f64, StaticError> {
        let s: Vec<WeightedSample> = values
            .iter()
            .map(|&v| WeightedSample::new(1.0, v))
            .collect();
        match self {
            UpscaleMethod::Arithmetic => arithmetic_mean(&s),
            UpscaleMethod::Harmonic => harmonic_mean(&s),
            UpscaleMethod::Geometric => geometric_mean(&s),
        }
    }
}

/// The SGS propagation spec (`ps.gaussian(...)`): the spatial-continuity model, the
/// RNG seed, the moving-neighbourhood search, and an optional collocated secondary
/// (trend cokriging, Phase 3).
///
/// The variogram is fitted in **data space**; SGS transforms to normal scores
/// internally (`petektools` `NormalScore`), so its total sill need not be 1.
#[derive(Debug, Clone)]
pub struct Gaussian {
    variogram: SpatialVariogram,
    seed: u64,
    /// Moving-neighbourhood `(max_neighbours, radius_m)`; `None` = derive a
    /// **bounded** neighbourhood at propagate time (`DEFAULT_MAX_NEIGHBOURS` + a
    /// variogram-range/lattice-scaled radius — see [`propagate_sgs_into`]). Use
    /// [`Gaussian::with_unbounded_search`] to restore the old whole-grid window.
    search: Option<(usize, f64)>,
    /// Collocated secondary `(trend, corr)` (Markov-1 cokriging). `None` = plain SGS.
    trend: Option<(TrendSurface, f64)>,
    /// Opt-out: allow a simulated layer with **no conditioning data** to be filled
    /// with the conditioned mean (a flat, structureless layer) instead of erroring.
    /// Default `false` — a data-less layer in the simulated range is a hard error
    /// naming the property (and, via the caller, the zone), because a silent
    /// constant mean-fill loses all spatial structure (and any collocated trend) on
    /// that layer with no warning (`task_petekstatic_canonical_fixes` item 4).
    allow_mean_fill: bool,
    /// Opt-in to the legacy **whole-grid** search window (the old default) instead of
    /// the bounded default, when [`Gaussian::with_search`] is not set. Pathologically
    /// slow on real-scale lattices; kept for small grids / exact pre-bounded results.
    unbounded_search: bool,
}

/// Default moving-neighbourhood node cap when [`Gaussian::with_search`] is not set.
const DEFAULT_MAX_NEIGHBOURS: usize = 16;

impl Gaussian {
    /// A plain (no-secondary) SGS spec from a data-space variogram and a seed. The
    /// search neighbourhood defaults to a **bounded** window (`DEFAULT_MAX_NEIGHBOURS`
    /// nodes within a variogram-range/lattice-scaled radius) — seconds on real-scale
    /// lattices where the old whole-grid default took >15 min/cube. Beyond a
    /// covariance range a conditioning node's kriging weight is ~0, so the bounded
    /// result matches the unbounded one within simulation tolerance.
    #[must_use]
    pub fn new(variogram: impl Into<SpatialVariogram>, seed: u64) -> Self {
        Self {
            variogram: variogram.into(),
            seed,
            search: None,
            trend: None,
            allow_mean_fill: false,
            unbounded_search: false,
        }
    }

    /// Override the moving-neighbourhood search: at most `max_neighbours`
    /// conditioning nodes within `radius_m` (world units).
    #[must_use]
    pub fn with_search(mut self, max_neighbours: usize, radius_m: f64) -> Self {
        self.search = Some((max_neighbours, radius_m));
        self
    }

    /// Restore the legacy **unbounded** (whole-grid) search window — every
    /// previously-simulated node in the lattice is a search candidate. Pathologically
    /// slow on real-scale lattices; kept as an explicit opt-in for small grids or
    /// exact reproduction of pre-bounded-default results.
    #[must_use]
    pub fn with_unbounded_search(mut self) -> Self {
        self.unbounded_search = true;
        self
    }

    /// Opt into filling a simulated layer that carries **no conditioning data** with
    /// the conditioned mean (a flat layer) rather than erroring — for a model where a
    /// data-less layer is expected and a structureless mean-fill is acceptable. The
    /// default errors loudly instead (`task_petekstatic_canonical_fixes` item 4).
    #[must_use]
    pub fn allow_mean_fill(mut self) -> Self {
        self.allow_mean_fill = true;
        self
    }

    /// Steer the simulation with a collocated secondary variable (trend cokriging,
    /// Markov-1): the [`TrendSurface`] is resampled to the model areal lattice and
    /// folded into each node's kriging at correlation `corr`. `corr == 0` recovers
    /// plain SGS bit-for-bit; `corr -> 1` pulls the field toward the trend pattern.
    #[must_use]
    pub fn with_trend(mut self, trend: TrendSurface, corr: f64) -> Self {
        self.trend = Some((trend, corr));
        self
    }
}

/// The upscaled-vs-log QC report for one property's [`PropertyPipeline::upscale_cells`]
/// step — the numbers a modeller inspects before propagating.
#[derive(Debug, Clone, PartialEq)]
pub struct UpscaleQc {
    /// The property name.
    pub property: String,
    /// Cells that received an upscaled log value (the conditioning set).
    pub conditioned_cells: usize,
    /// Total log samples that fell inside a cell's depth range.
    pub log_samples: usize,
    /// Mean of the in-cell log samples (`NaN` if none) — the "before upscaling" level.
    pub log_mean: f64,
    /// Mean of the per-cell upscaled values (`NaN` if none) — the "after" level.
    pub upscaled_mean: f64,
    /// Min / max of the per-cell upscaled values (`NaN` if none).
    pub upscaled_min: f64,
    pub upscaled_max: f64,
}

/// The report from running a full [`PropertyPipeline`] against a grid: the upscale
/// QC plus whether propagation ran.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyReport {
    /// The property name.
    pub property: String,
    /// The upscale-step QC.
    pub upscale: UpscaleQc,
    /// Whether an SGS propagation step filled the whole cube (else only the
    /// conditioned cells carry values).
    pub propagated: bool,
}

/// The per-property pipeline: a fluent spec (property name + optional upscale +
/// optional propagate) that both the deterministic builder and the MC template
/// replay against a grid. Immutable/chainable (house style).
#[derive(Debug, Clone)]
pub struct PropertyPipeline {
    name: String,
    upscale: Option<(Vec<WellLog>, UpscaleMethod)>,
    propagate: Option<Gaussian>,
}

impl PropertyPipeline {
    /// Start a pipeline for the named property (e.g. `"PHIE"`, `"NTG"`).
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            upscale: None,
            propagate: None,
        }
    }

    /// Attach the upscale step: positioned well logs + the averaging method. The
    /// visible first step — [`PropertyPipeline::upscale_cells`] runs it and returns
    /// the conditioned cells + QC.
    #[must_use]
    pub fn upscale(mut self, wells: Vec<WellLog>, method: UpscaleMethod) -> Self {
        self.upscale = Some((wells, method));
        self
    }

    /// Attach the propagation step (SGS). Requires an [`PropertyPipeline::upscale`]
    /// step to condition on.
    #[must_use]
    pub fn propagate(mut self, gaussian: Gaussian) -> Self {
        self.propagate = Some(gaussian);
        self
    }

    /// The property this pipeline populates.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// A copy whose SGS seed is XOR-mixed with `salt` — the per-draw reseed for the
    /// [`McMode::Resimulate`] path (a new pattern per realization, reproducible from
    /// the draw's seed index). A no-op on a pipeline without a propagate step.
    pub(crate) fn reseeded(&self, salt: u64) -> Self {
        let mut p = self.clone();
        if let Some(g) = p.propagate {
            p.propagate = Some(Gaussian {
                seed: g.seed ^ salt.wrapping_mul(SEED_GOLDEN),
                ..g
            });
        }
        p
    }

    /// The upscale step: bin each well's samples into the cells its column passes
    /// through and upscale them, returning the per-cell conditioned field
    /// (`NaN` where no log passes) and the [`UpscaleQc`] report. Row-major
    /// `(k * nj + j) * ni + i` (matches the property cube).
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the pipeline has no upscale step, the grid's
    /// areal lattice is smaller than `2x2` columns or degenerate, or an upscale
    /// average rejects a sample set.
    pub fn upscale_cells(&self, grid: &Grid) -> Result<(Vec<f64>, UpscaleQc), StaticError> {
        self.upscale_cells_with_georef(grid, None)
    }

    /// World-frame variant used by model builders/templates. The public
    /// [`Self::upscale_cells`] remains the exact local-frame compatibility path.
    fn upscale_cells_with_georef(
        &self,
        grid: &Grid,
        georef: Option<Georef>,
    ) -> Result<(Vec<f64>, UpscaleQc), StaticError> {
        let (wells, method) = self.upscale.as_ref().ok_or_else(|| {
            StaticError::InvalidInput(format!(
                "property '{}' has no upscale step (call .upscale(...))",
                self.name
            ))
        })?;
        let dims = grid.dims();
        let (ni, nj, nk) = (dims.ni, dims.nj, dims.nk);
        let lattice = areal_lattice(grid)?;
        let mut cells = vec![f64::NAN; ni * nj * nk];
        let mut log_values: Vec<f64> = Vec::new();

        for well in wells {
            // `WellLog::new` historically accepted local lattice coordinates even
            // when a model carried an unrotated georef. Preserve that path exactly,
            // while oriented models prefer the unambiguous world-frame inverse.
            // The alternate interpretation is only a compatibility fallback when
            // the preferred one lies outside the model.
            let in_grid = |p: Option<(f64, f64)>| {
                p.filter(|(fi, fj)| {
                    let (i, j) = (fi.round(), fj.round());
                    i >= 0.0 && j >= 0.0 && i < ni as f64 && j < nj as f64
                })
            };
            let local = || lattice.xy_to_ij(well.x, well.y);
            let intrinsic = match georef {
                Some(g) if g.rotation_deg != 0.0 || g.yflip => {
                    in_grid(g.world_to_intrinsic(well.x, well.y)).or_else(|| in_grid(local()))
                }
                Some(g) => {
                    in_grid(local()).or_else(|| in_grid(g.world_to_intrinsic(well.x, well.y)))
                }
                None => in_grid(local()),
            };
            let Some((raw_fi, raw_fj)) = intrinsic else {
                continue;
            };
            let local_xy = lattice_intrinsic_to_world(&lattice, raw_fi, raw_fj);
            let (fi, fj) = (raw_fi.round(), raw_fj.round());
            if fi < 0.0 || fj < 0.0 {
                continue;
            }
            let (i, j) = (fi as usize, fj as usize);
            if i >= ni || j >= nj {
                continue;
            }
            for k in 0..nk {
                let cell = grid.cell(Ijk::new(i, j, k));
                // Bin against the cell's depth range **interpolated at the well's true
                // (x, y)** — NOT `cell.top_depth()`/`bottom_depth()`, which are the
                // 4-corner means (the column *centroid* interval). A well rarely sits at
                // its column's centroid; on a dipping horizon the centroid-interpolated
                // zone boundary is offset from the boundary depth *at the well*, so
                // centroid binning mis-assigns near-boundary samples to the neighbouring
                // zone — the per-zone upscale then averages cross-zone samples and
                // compresses each zone's proportion toward the mid-range
                // (`task_petekstatic_zoned_fixes` finding 3). Sampling the cell interval
                // at the well position assigns each sample to the zone its depth truly
                // falls in. A well exactly at the centroid recovers the old means
                // (bilinear at (0.5, 0.5) == the 4-corner mean).
                let (t, b) = cell_depth_at_xy(&cell, local_xy.0, local_xy.1);
                let (lo, hi) = (t.min(b), t.max(b));
                let in_range: Vec<f64> = well
                    .samples
                    .iter()
                    .filter(|(tvd, _)| *tvd >= lo && *tvd <= hi)
                    .map(|(_, v)| *v)
                    .collect();
                if !in_range.is_empty() {
                    log_values.extend_from_slice(&in_range);
                    cells[(k * nj + j) * ni + i] = method.apply(&in_range)?;
                }
            }
        }

        let qc = upscale_qc(&self.name, &cells, &log_values);
        Ok((cells, qc))
    }

    /// Run the full pipeline against a grid: upscale, then SGS-propagate the whole
    /// cube conditioned on the upscaled cells, setting the property cube. Returns
    /// the [`PropertyReport`].
    ///
    /// If the pipeline has no propagate step, only the conditioned cells carry
    /// values (the rest stay `NaN`) and `propagated` is `false`.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] as [`PropertyPipeline::upscale_cells`], if there
    /// is no conditioning data to propagate from, or if a variogram/SGS solve fails.
    pub fn apply(&self, grid: &mut Grid) -> Result<PropertyReport, StaticError> {
        self.apply_with_georef(grid, None)
    }

    /// [`PropertyPipeline::apply`] with the model's world [`Georef`] (crate-internal;
    /// the builder / template thread their own `with_georef`). The georef maps the
    /// grid's local column lattice to world coordinates so a **world-georeferenced
    /// collocated trend** ([`Gaussian::with_trend`]) is resampled at each column's
    /// true world position rather than onto the local lattice — without it a loaded
    /// world trend never overlaps the local frame and the cokriging silently no-ops
    /// (`task_petekstatic_zoned_fixes` finding 1). `None` keeps the local/index-space
    /// frame (a synthetic model with no world georeference).
    pub(crate) fn apply_with_georef(
        &self,
        grid: &mut Grid,
        georef: Option<Georef>,
    ) -> Result<PropertyReport, StaticError> {
        let (cells, qc) = self.upscale_cells_with_georef(grid, georef)?;
        let propagated = match &self.propagate {
            None => {
                grid.properties_mut().set(Property {
                    name: self.name.clone(),
                    values: cells,
                })?;
                false
            }
            Some(g) => {
                let values = propagate_sgs(grid, &self.name, &cells, g, georef)?;
                grid.properties_mut().set(Property {
                    name: self.name.clone(),
                    values,
                })?;
                true
            }
        };
        Ok(PropertyReport {
            property: self.name.clone(),
            upscale: qc,
            propagated,
        })
    }

    /// Run the pipeline **restricted to one zone's `k`-range** (P8 per-zone
    /// population, `task_petekstatic_multizone_2`): condition + propagate only the
    /// cells in `k_range`, **merging** into the property cube already on the grid so
    /// the other zones' slices are untouched. This is the per-zone analog of
    /// [`PropertyPipeline::apply`] — each zone in a stack gets its own
    /// variogram/trend/log-conditioning by attaching one pipeline per zone
    /// ([`crate::model::StaticModelBuilder::with_zone_property`]).
    ///
    /// - Conditioning is scoped to `k_range` (upscaled cells outside it are ignored),
    ///   so a zone's SGS sees only its own logs and its own normal-score transform.
    /// - With a propagate step, every cell in `k_range` is filled (a data-less layer
    ///   in-range falls back to the zone's conditioned mean); cells outside `k_range`
    ///   keep the grid's existing value (the base priors / another zone's field).
    /// - Without a propagate step, only the in-range conditioned (well-column) cells
    ///   are overwritten; the rest of the zone keeps its base value.
    ///
    /// # Errors
    /// As [`PropertyPipeline::apply`]; additionally the zone must carry at least one
    /// conditioning datum when a propagate step is present.
    pub fn apply_in_zone(
        &self,
        grid: &mut Grid,
        k_range: core::ops::Range<usize>,
    ) -> Result<PropertyReport, StaticError> {
        self.apply_in_zone_with_georef(grid, k_range, None)
    }

    /// [`PropertyPipeline::apply_in_zone`] with the model's world [`Georef`]
    /// (crate-internal) — the per-zone analog of [`PropertyPipeline::apply_with_georef`],
    /// so a per-zone collocated trend resamples in the correct world frame (finding 1).
    pub(crate) fn apply_in_zone_with_georef(
        &self,
        grid: &mut Grid,
        k_range: core::ops::Range<usize>,
        georef: Option<Georef>,
    ) -> Result<PropertyReport, StaticError> {
        let dims = grid.dims();
        let (ni, nj, nk) = (dims.ni, dims.nj, dims.nk);
        let (mut cells, _) = self.upscale_cells_with_georef(grid, georef)?;
        // Scope conditioning to the zone: drop any upscaled cell outside k_range so
        // the QC + the SGS normal-score transform reflect only this zone's data.
        for k in 0..nk {
            if !k_range.contains(&k) {
                for j in 0..nj {
                    for i in 0..ni {
                        cells[(k * nj + j) * ni + i] = f64::NAN;
                    }
                }
            }
        }
        let log_values: Vec<f64> = cells.iter().copied().filter(|v| v.is_finite()).collect();
        let qc = upscale_qc(&self.name, &cells, &log_values);
        // Base = the cube already on the grid (priors / log population), or NaN if the
        // property is not yet present — the zone's slice merges into it.
        let mut base = grid
            .properties()
            .get(&self.name)
            .map(|p| p.values.clone())
            .unwrap_or_else(|| vec![f64::NAN; ni * nj * nk]);
        let propagated = match &self.propagate {
            None => {
                for k in k_range.clone() {
                    for j in 0..nj {
                        for i in 0..ni {
                            let idx = (k * nj + j) * ni + i;
                            if cells[idx].is_finite() {
                                base[idx] = cells[idx];
                            }
                        }
                    }
                }
                false
            }
            Some(g) => {
                propagate_sgs_into(grid, &self.name, &cells, g, k_range, &mut base, georef)?;
                true
            }
        };
        grid.properties_mut().set(Property {
            name: self.name.clone(),
            values: base,
        })?;
        Ok(PropertyReport {
            property: self.name.clone(),
            upscale: qc,
            propagated,
        })
    }
}

/// The cell's `(top_depth, bottom_depth)` **bilinearly interpolated at world
/// `(wx, wy)`** — the depth range the cell presents *at the well*, rather than the
/// 4-corner means ([`crate::grid::Cell::top_depth`]/`bottom_depth`, the column
/// centroid). Used by the log-upscale binning so a well that does not sit at its
/// column centroid assigns each sample to the zone its depth truly falls in on a
/// dipping horizon (finding 3). A query at the cell centre recovers the 4-corner
/// mean, so centred wells are unchanged.
///
/// Corners are ordered `di + 2*dj + 4*dk`: top face `0..4` (nodes (i,j), (i+1,j),
/// (i,j+1), (i+1,j+1)), bottom face `4..8` in the same areal order. The lattice is
/// axis-aligned/regular (the module limitation), so the areal quad is a rectangle;
/// the interpolation weights clamp to the cell so an off-lattice query stays inside.
fn cell_depth_at_xy(cell: &crate::grid::Cell, wx: f64, wy: f64) -> (f64, f64) {
    let c = &cell.corners;
    let (x0, x1) = (c[0].x, c[1].x);
    let (y0, y1) = (c[0].y, c[2].y);
    let frac = |q: f64, a: f64, b: f64| {
        let d = b - a;
        if d.abs() > f64::EPSILON {
            ((q - a) / d).clamp(0.0, 1.0)
        } else {
            0.5
        }
    };
    let (tx, ty) = (frac(wx, x0, x1), frac(wy, y0, y1));
    let bilinear = |z00: f64, z10: f64, z01: f64, z11: f64| {
        (1.0 - tx) * (1.0 - ty) * z00
            + tx * (1.0 - ty) * z10
            + (1.0 - tx) * ty * z01
            + tx * ty * z11
    };
    let top = bilinear(c[0].z, c[1].z, c[2].z, c[3].z);
    let bottom = bilinear(c[4].z, c[5].z, c[6].z, c[7].z);
    (top, bottom)
}

/// Place fractional intrinsic coordinates through petekTools' canonical lattice
/// transform without reimplementing rotation math in petekStatic.
fn lattice_intrinsic_to_world(lat: &Lattice, fi: f64, fj: f64) -> (f64, f64) {
    Lattice {
        xori: lat.xori,
        yori: lat.yori,
        xinc: lat.xinc * fi,
        yinc: lat.yinc * fj,
        ncol: 2,
        nrow: 2,
        rotation_deg: lat.rotation_deg,
        yflip: lat.yflip,
    }
    .node_xy(1, 1)
}

/// Reconstruct the areal SGS lattice from the grid's top-layer cell centroids
/// (origin + node spacing). Requires an at-least-`2x2`, regular, axis-aligned
/// column layout (see the module limitation).
pub(crate) fn areal_lattice(grid: &Grid) -> Result<Lattice, StaticError> {
    let dims = grid.dims();
    let (ni, nj) = (dims.ni, dims.nj);
    if ni < 2 || nj < 2 {
        return Err(StaticError::InvalidInput(format!(
            "property pipeline needs an areal lattice of at least 2x2 columns, got {ni}x{nj}"
        )));
    }
    let c00 = grid.cell(Ijk::new(0, 0, 0)).centroid();
    let c10 = grid.cell(Ijk::new(1, 0, 0)).centroid();
    let c01 = grid.cell(Ijk::new(0, 1, 0)).centroid();
    let dx = c10.x - c00.x;
    let dy = c01.y - c00.y;
    if !(dx.is_finite() && dx > 0.0 && dy.is_finite() && dy > 0.0) {
        return Err(StaticError::InvalidInput(format!(
            "property pipeline needs a regular axis-aligned column lattice (node spacing dx={dx}, dy={dy})"
        )));
    }
    Ok(Lattice::regular(c00.x, c00.y, dx, dy, ni, nj))
}

/// Per-layer SGS conditioned on the upscaled cells over the whole cube. Returns the
/// full cube (row-major `(k * nj + j) * ni + i`); a layer with no conditioning data
/// is filled with the global conditioned mean (never left `NaN`). Thin wrapper over
/// [`propagate_sgs_into`] with the full `0..nk` range and a fresh `NaN` base.
fn propagate_sgs(
    grid: &Grid,
    property: &str,
    cells: &[f64],
    g: &Gaussian,
    georef: Option<Georef>,
) -> Result<Vec<f64>, StaticError> {
    let dims = grid.dims();
    let mut base = vec![f64::NAN; dims.ni * dims.nj * dims.nk];
    propagate_sgs_into(grid, property, cells, g, 0..dims.nk, &mut base, georef)?;
    Ok(base)
}

fn gaussian_search_range(g: &Gaussian) -> f64 {
    match g.variogram {
        SpatialVariogram::Isotropic(v) => v.range,
        SpatialVariogram::Anisotropic(v) => v.major.max(v.minor).max(v.vertical),
    }
}

/// Per-layer SGS over `k_range` only, writing into `base` (the cube already on the
/// grid) — the per-zone population core. Only layers in `k_range` are simulated and
/// overwritten; layers outside it keep their `base` value. The seed derivation is
/// **per-global-k** (`g.seed ^ (k+1)·GOLDEN`), independent of the range, so a zone's
/// field is identical whether it is populated alone or as part of a wider run
/// (reproducibility across the whole-grid vs per-zone paths).
fn propagate_sgs_into(
    grid: &Grid,
    property: &str,
    cells: &[f64],
    g: &Gaussian,
    k_range: core::ops::Range<usize>,
    base: &mut [f64],
    georef: Option<Georef>,
) -> Result<(), StaticError> {
    let dims = grid.dims();
    let (ni, nj) = (dims.ni, dims.nj);
    let lattice = areal_lattice(grid)?;

    // Conditioning stats over the range being populated (fallback + "must have data"
    // gate). For the full-grid path this is the global set; for a zone it is the
    // zone's own conditioning (the caller nulls out-of-range cells).
    let conditioned: Vec<f64> = k_range
        .clone()
        .flat_map(|k| (0..ni * nj).map(move |c| cells[k * ni * nj + c]))
        .filter(|v| v.is_finite())
        .collect();
    if conditioned.is_empty() {
        return Err(StaticError::InvalidInput(format!(
            "property '{property}': no conditioning data in the simulated range \
             (upscale produced no informed cells) — the wells' samples do not fall in \
             any cell of this range. Check the well positions/depths against the grid \
             (and, for a zone-scoped pipe, that the wells penetrate the zone)."
        )));
    }
    let global_mean = conditioned.iter().sum::<f64>() / conditioned.len() as f64;

    // Column world positions (top-layer centroids) — the areal position of every
    // conditioning datum in a column. Regular lattice => these snap back to (i, j).
    let mut column_xy = vec![(0.0, 0.0); ni * nj];
    for j in 0..nj {
        for i in 0..ni {
            let c = grid.cell(Ijk::new(i, j, 0)).centroid();
            column_xy[j * ni + i] = (c.x, c.y);
        }
    }

    // Search neighbourhood: caller override, else a BOUNDED default —
    // `DEFAULT_MAX_NEIGHBOURS` nodes within a radius scaled to the covariance range
    // (with a lattice-spacing floor so a tiny range still finds neighbours). The old
    // whole-grid default made every simulated node a candidate → O(N²) scans, >15
    // min/cube on real-scale lattices. Beyond a covariance range a node's kriging
    // weight is ~0, so 1.5× the range captures every materially-conditioning node;
    // results match the unbounded window within simulation tolerance. Explicit
    // opt-in to the old behaviour via `Gaussian::with_unbounded_search`.
    let (max_neighbours, radius) = g.search.unwrap_or_else(|| {
        let spacing = lattice.xinc.max(lattice.yinc);
        if g.unbounded_search {
            // Legacy whole-grid window (opt-in): the diagonal extent, same node cap.
            let extent = ((ni as f64) * lattice.xinc).hypot((nj as f64) * lattice.yinc);
            (DEFAULT_MAX_NEIGHBOURS, extent.max(spacing))
        } else {
            let radius = (gaussian_search_range(g) * 1.5).max(spacing * 4.0);
            (DEFAULT_MAX_NEIGHBOURS, radius)
        }
    });

    // Collocated secondary resampled to the model lattice (Phase 3). A
    // world-georeferenced trend must be sampled at each column's **world** position,
    // not on the grid's local lattice — the grid geometry stays local (area-scaled at
    // a local origin) while a loaded trend carries its world georeference, so
    // resampling the world trend onto the local lattice yields an all-NaN secondary
    // and the kernel silently drops every node to plain SK (`task_petekstatic_zoned_fixes`
    // finding 1). Map the local lattice to world through the model [`Georef`] before
    // handing the secondary to the frame-agnostic kernel. When there is no world georef
    // (a synthetic model) or the trend is not georeferenced, the local lattice IS the
    // sampling frame (index-space / matching-local trends).
    let secondary = match &g.trend {
        None => None,
        Some((trend, corr)) => {
            let sample_lattice = match georef {
                Some(gr) if trend.is_georeferenced() => Lattice {
                    xori: gr.origin_x,
                    yori: gr.origin_y,
                    xinc: gr.spacing_x,
                    yinc: gr.spacing_y,
                    ncol: ni,
                    nrow: nj,
                    rotation_deg: gr.rotation_deg,
                    yflip: gr.yflip,
                },
                _ => lattice.clone(),
            };
            let field = trend.resample_to(&sample_lattice)?;
            // Loud fallback: a secondary that covers too little of the frame is almost
            // always a georef mismatch (the world-vs-local no-op covers 0 nodes), never
            // intent — error rather than silently steering nothing.
            let finite = field.iter().filter(|v| v.is_finite()).count();
            let total = field.len().max(1);
            if (finite as f64) < COLLOCATED_MIN_COVERAGE * total as f64 {
                return Err(StaticError::InvalidInput(format!(
                    "collocated trend covers only {finite}/{total} areal nodes (< {:.0}% of the \
                     model frame) — the trend's georeference does not overlap the grid. A \
                     world-georeferenced trend needs the model's world georef (build with \
                     `with_georef`); a local trend must match the grid's local lattice.",
                    COLLOCATED_MIN_COVERAGE * 100.0
                )));
            }
            Some((field, *corr))
        }
    };

    // Built once and shared **immutably** across layers: the collocated secondary
    // (a full lattice `Array2`) is moved in ONCE, and `params.seed` is overridden
    // per layer at the call site ([`sgs_seeded`]) rather than mutated — so a single
    // `&params` can be borrowed across the parallel layers below.
    let params = SgsParams {
        variogram: g.variogram,
        max_neighbours,
        radius,
        seed: g.seed, // per-layer seed overridden via `sgs_seeded`
        collocated: secondary,
    };

    // Each k-layer is an INDEPENDENT 2D SGS: its own per-global-`k` seed
    // (`g.seed ^ (k+1)·GOLDEN`, range-independent), its own conditioning, its own
    // disjoint output slice, and only an immutable `&params`/`&lattice` read. So
    // simulate the layers in **parallel** and scatter each field back in `k`-order
    // — bit-for-bit identical to the serial sweep (the seed is a pure function of
    // `k`, never of execution order). Real-scale zones carry hundreds of 1 m
    // layers, so this is the dominant per-build property-population lever.
    let layers: Vec<usize> = k_range.collect();
    let fields: Vec<(usize, Vec<f64>)> = layers
        .par_iter()
        .map(|&k| -> Result<(usize, Vec<f64>), StaticError> {
            // This layer's conditioning data as [x, y, value] rows.
            let mut coords: Vec<[f64; 3]> = Vec::new();
            for j in 0..nj {
                for i in 0..ni {
                    let v = cells[(k * nj + j) * ni + i];
                    if v.is_finite() {
                        let (x, y) = column_xy[j * ni + i];
                        coords.push([x, y, v]);
                    }
                }
            }
            let mut layer = vec![0.0_f64; ni * nj];
            if coords.is_empty() {
                // No conditioning data in this layer. A silent constant mean-fill loses
                // all spatial structure (and any collocated trend) on the layer, so by
                // default this is a hard error naming the property; opt in to the fill
                // with `Gaussian::allow_mean_fill` (`task_petekstatic_canonical_fixes`
                // item 4). The zone-scoped callers wrap this to name the zone too.
                if !g.allow_mean_fill {
                    return Err(StaticError::InvalidInput(format!(
                        "property '{property}': simulated layer {k} has no conditioning data \
                         — a silent constant mean-fill would erase its spatial structure. \
                         Ensure the wells condition every simulated layer, or opt into the \
                         structureless mean-fill with `Gaussian::allow_mean_fill`."
                    )));
                }
                // Opted in: fill with the (range) conditioned mean.
                layer.iter_mut().for_each(|v| *v = global_mean);
                return Ok((k, layer));
            }
            // Independent yet reproducible per-layer seed (a pure function of `k`).
            let seed = g.seed ^ (k as u64).wrapping_add(1).wrapping_mul(SEED_GOLDEN);
            let field = sgs_seeded(&coords, &lattice, &params, seed).map_err(|e| {
                StaticError::Grid(format!("SGS propagate failed on layer {k}: {e}"))
            })?;
            for j in 0..nj {
                for i in 0..ni {
                    layer[j * ni + i] = field[[i, j]];
                }
            }
            Ok((k, layer))
        })
        .collect::<Result<Vec<_>, StaticError>>()?;
    for (k, layer) in fields {
        base[k * ni * nj..(k + 1) * ni * nj].copy_from_slice(&layer);
    }
    Ok(())
}

/// Assemble the [`UpscaleQc`] from the conditioned cells and the in-cell log values.
fn upscale_qc(property: &str, cells: &[f64], log_values: &[f64]) -> UpscaleQc {
    let conditioned: Vec<f64> = cells.iter().copied().filter(|v| v.is_finite()).collect();
    let mean = |v: &[f64]| {
        if v.is_empty() {
            f64::NAN
        } else {
            v.iter().sum::<f64>() / v.len() as f64
        }
    };
    let (mut min, mut max) = (f64::NAN, f64::NAN);
    if !conditioned.is_empty() {
        min = conditioned.iter().copied().fold(f64::INFINITY, f64::min);
        max = conditioned
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
    }
    UpscaleQc {
        property: property.to_string(),
        conditioned_cells: conditioned.len(),
        log_samples: log_values.len(),
        log_mean: mean(log_values),
        upscaled_mean: mean(&conditioned),
        upscaled_min: min,
        upscaled_max: max,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::{build_box, BoxSpec, Dims};
    use petektools::{AnisotropicVariogram, Variogram, VariogramModel};

    const NI: usize = 5;
    const NJ: usize = 5;
    const NK: usize = 4;
    // 100 m square, dx = dy = 20 m; gross 40 m over 4 layers => dz = 10 m, top 0.
    // Column (i,j) centroid = ((i+0.5)*20, (j+0.5)*20); layer k spans [10k, 10k+10].

    fn test_grid() -> Grid {
        build_box(BoxSpec::square(
            10_000.0,
            40.0,
            Dims::new(NI, NJ, NK).unwrap(),
        ))
        .unwrap()
    }

    /// Two wells at opposite corner columns (0,0) and (4,4), each carrying one
    /// sample in every layer so both condition all four layers. Column (0,0) reads
    /// low, column (4,4) high — a strong areal contrast.
    fn two_wells() -> Vec<WellLog> {
        // tvd 5,15,25,35 -> layers 0,1,2,3.
        let low = WellLog::new(
            10.0,
            10.0,
            vec![(5.0, 0.10), (15.0, 0.12), (25.0, 0.14), (35.0, 0.16)],
        );
        let high = WellLog::new(
            90.0,
            90.0,
            vec![(5.0, 0.28), (15.0, 0.26), (25.0, 0.24), (35.0, 0.22)],
        );
        vec![low, high]
    }

    fn nscore_variogram() -> Variogram {
        // Fitted on normal scores: unit sill, range ~half the field, no nugget.
        Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 50.0).unwrap()
    }

    fn idx(i: usize, j: usize, k: usize) -> usize {
        (k * NJ + j) * NI + i
    }

    #[test]
    fn upscale_honours_logs_at_well_columns() {
        let grid = test_grid();
        let pipe = PropertyPipeline::new("PHIE").upscale(two_wells(), UpscaleMethod::Arithmetic);
        let (cells, qc) = pipe.upscale_cells(&grid).unwrap();

        // Column (0,0): one sample per layer, so the upscaled cell == that sample.
        for (k, v) in [(0, 0.10), (1, 0.12), (2, 0.14), (3, 0.16)] {
            assert!(
                (cells[idx(0, 0, k)] - v).abs() < 1e-12,
                "cell (0,0,{k}) = {} != {v}",
                cells[idx(0, 0, k)]
            );
        }
        // Column (4,4).
        for (k, v) in [(0, 0.28), (1, 0.26), (2, 0.24), (3, 0.22)] {
            assert!((cells[idx(4, 4, k)] - v).abs() < 1e-12);
        }
        // Every other column is unconditioned (NaN).
        assert!(cells[idx(2, 2, 0)].is_nan());
        // QC: 8 conditioned cells (2 wells x 4 layers), 8 log samples.
        assert_eq!(qc.conditioned_cells, 8);
        assert_eq!(qc.log_samples, 8);
        assert!((qc.log_mean - qc.upscaled_mean).abs() < 1e-12);
        assert!((qc.upscaled_min - 0.10).abs() < 1e-12);
        assert!((qc.upscaled_max - 0.28).abs() < 1e-12);
    }

    #[test]
    fn upscale_averages_multiple_samples_in_a_cell() {
        let grid = test_grid();
        // Three samples all in layer 0 ([0,10]) of column (0,0): mean 0.20.
        let well = WellLog::new(10.0, 10.0, vec![(2.0, 0.18), (5.0, 0.20), (8.0, 0.22)]);
        let pipe = PropertyPipeline::new("PHIE").upscale(vec![well], UpscaleMethod::Arithmetic);
        let (cells, qc) = pipe.upscale_cells(&grid).unwrap();
        assert!((cells[idx(0, 0, 0)] - 0.20).abs() < 1e-12);
        assert_eq!(qc.conditioned_cells, 1);
        assert_eq!(qc.log_samples, 3);
    }

    #[test]
    fn sgs_propagation_honours_conditioning_and_fills_every_cell() {
        let mut grid = test_grid();
        let pipe = PropertyPipeline::new("PHIE")
            .upscale(two_wells(), UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(nscore_variogram(), 1));
        let report = pipe.apply(&mut grid).unwrap();
        assert!(report.propagated);

        let prop = grid.properties().get("PHIE").unwrap();
        // Every cell filled (no NaN).
        assert!(prop.values.iter().all(|v| v.is_finite()), "cube has NaN");
        // SGS honours conditioning exactly at the data nodes (well columns).
        for k in 0..NK {
            let lo = [0.10, 0.12, 0.14, 0.16][k];
            let hi = [0.28, 0.26, 0.24, 0.22][k];
            assert!(
                (prop.values[idx(0, 0, k)] - lo).abs() < 1e-6,
                "well (0,0,{k}) not honoured: {}",
                prop.values[idx(0, 0, k)]
            );
            assert!((prop.values[idx(4, 4, k)] - hi).abs() < 1e-6);
        }
    }

    #[test]
    fn sgs_propagation_accepts_anisotropic_variogram() {
        let mut grid = test_grid();
        let vgm =
            AnisotropicVariogram::new(VariogramModel::Spherical, 0.0, 1.0, 80.0, 30.0, 10.0, 35.0)
                .unwrap();
        let pipe = PropertyPipeline::new("PHIE")
            .upscale(two_wells(), UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(vgm, 3));

        let report = pipe.apply(&mut grid).unwrap();
        let prop = grid.properties().get("PHIE").unwrap();

        assert!(report.propagated);
        assert!(prop.values.iter().all(|v| v.is_finite()));
        assert!((prop.values[idx(0, 0, 0)] - 0.10).abs() < 1e-6);
        assert!((prop.values[idx(4, 4, 0)] - 0.28).abs() < 1e-6);
    }

    #[test]
    fn sgs_reproduces_the_data_histogram_loosely() {
        let mut grid = test_grid();
        let pipe = PropertyPipeline::new("PHIE")
            .upscale(two_wells(), UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(nscore_variogram(), 7));
        pipe.apply(&mut grid).unwrap();
        let prop = grid.properties().get("PHIE").unwrap();
        let field_mean = prop.values.iter().sum::<f64>() / prop.values.len() as f64;
        // Conditioning mean is 0.19; the simulated field tracks it within a few
        // hundredths (loose, seeded), and stays inside the data range [0.10, 0.28].
        assert!((field_mean - 0.19).abs() < 0.05, "field mean {field_mean}");
        assert!(
            prop.values.iter().all(|&v| (0.09..=0.29).contains(&v)),
            "field escaped the data range"
        );
    }

    #[test]
    fn sgs_is_seeded_reproducible_and_seed_sensitive() {
        let run = |seed: u64| {
            let mut grid = test_grid();
            PropertyPipeline::new("PHIE")
                .upscale(two_wells(), UpscaleMethod::Arithmetic)
                .propagate(Gaussian::new(nscore_variogram(), seed))
                .apply(&mut grid)
                .unwrap();
            grid.properties().get("PHIE").unwrap().values.clone()
        };
        assert_eq!(
            run(1),
            run(1),
            "same seed must reproduce the cube bit-for-bit"
        );
        assert_ne!(
            run(1),
            run(2),
            "a different seed must give a different field"
        );
    }

    #[test]
    fn propagate_needs_conditioning_data() {
        let mut grid = test_grid();
        // A well that misses the grid extent entirely -> no conditioned cells.
        let off = WellLog::new(1e6, 1e6, vec![(5.0, 0.2)]);
        let pipe = PropertyPipeline::new("PHIE")
            .upscale(vec![off], UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(nscore_variogram(), 1));
        assert!(pipe.apply(&mut grid).is_err());
    }

    #[test]
    fn apply_without_propagate_sets_only_conditioned_cells() {
        let mut grid = test_grid();
        let report = PropertyPipeline::new("PHIE")
            .upscale(two_wells(), UpscaleMethod::Arithmetic)
            .apply(&mut grid)
            .unwrap();
        assert!(!report.propagated);
        let prop = grid.properties().get("PHIE").unwrap();
        assert!((prop.values[idx(0, 0, 0)] - 0.10).abs() < 1e-12);
        assert!(prop.values[idx(2, 2, 0)].is_nan());
    }

    #[test]
    fn apply_in_zone_fills_only_its_krange_and_merges_into_the_base() {
        // Per-zone population (P8): pre-fill the whole PHIE cube with a constant
        // baseline, then run the pipeline restricted to layers 1..3. Only those
        // layers are overwritten (SGS honouring the wells); layers 0 and 3 keep the
        // baseline exactly, and the well columns in-range are honoured.
        let mut grid = test_grid();
        grid.properties_mut()
            .set(Property::constant("PHIE", 0.5, NI * NJ * NK))
            .unwrap();
        let report = PropertyPipeline::new("PHIE")
            .upscale(two_wells(), UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(nscore_variogram(), 3))
            .apply_in_zone(&mut grid, 1..3)
            .unwrap();
        assert!(report.propagated);
        let v = &grid.properties().get("PHIE").unwrap().values;
        // Out-of-range layers untouched (still the 0.5 baseline).
        for (i, j, k) in [(0, 0, 0), (2, 2, 0), (4, 4, 3), (2, 2, 3)] {
            assert!(
                (v[idx(i, j, k)] - 0.5).abs() < 1e-12,
                "layer {k} outside 1..3 was overwritten"
            );
        }
        // In-range layers filled everywhere (no NaN) and honour the wells.
        for k in 1..3 {
            for j in 0..NJ {
                for i in 0..NI {
                    assert!(v[idx(i, j, k)].is_finite());
                }
            }
            let lo = [0.10, 0.12, 0.14, 0.16][k];
            let hi = [0.28, 0.26, 0.24, 0.22][k];
            assert!(
                (v[idx(0, 0, k)] - lo).abs() < 1e-6,
                "well (0,0,{k}) honoured"
            );
            assert!((v[idx(4, 4, k)] - hi).abs() < 1e-6);
        }
    }

    #[test]
    fn apply_in_zone_matches_full_apply_over_the_whole_range() {
        // apply_in_zone(0..nk) into a NaN-seeded cube reproduces apply() bit-for-bit
        // (the per-layer seed is per-global-k, range-independent) — so a zone's field
        // is identical whether populated alone or as part of the whole grid.
        let full = {
            let mut g = test_grid();
            PropertyPipeline::new("PHIE")
                .upscale(two_wells(), UpscaleMethod::Arithmetic)
                .propagate(Gaussian::new(nscore_variogram(), 5))
                .apply(&mut g)
                .unwrap();
            g.properties().get("PHIE").unwrap().values.clone()
        };
        let zoned = {
            let mut g = test_grid();
            PropertyPipeline::new("PHIE")
                .upscale(two_wells(), UpscaleMethod::Arithmetic)
                .propagate(Gaussian::new(nscore_variogram(), 5))
                .apply_in_zone(&mut g, 0..NK)
                .unwrap();
            g.properties().get("PHIE").unwrap().values.clone()
        };
        assert_eq!(full, zoned, "apply_in_zone(0..nk) must equal apply()");
    }

    #[test]
    fn tiny_lattice_is_rejected() {
        let grid = build_box(BoxSpec::square(10_000.0, 40.0, Dims::new(1, 1, 4).unwrap())).unwrap();
        let pipe = PropertyPipeline::new("PHIE").upscale(two_wells(), UpscaleMethod::Arithmetic);
        assert!(pipe.upscale_cells(&grid).is_err());
    }

    // --- collocated cokriging (with_trend) ---

    /// A trend georeferenced onto the exact column-centroid lattice, increasing with
    /// the column index i (nodes at world x = 10,30,50,70,90; value = i).
    fn i_increasing_trend() -> TrendSurface {
        let values: Vec<f64> = (0..NI * NJ).map(|k| (k % NI) as f64).collect();
        TrendSurface::new(NI, NJ, values)
            .unwrap()
            .with_georef(10.0, 10.0, 20.0, 20.0)
    }

    /// Two wells at the SAME column i=2, differing in j — so the conditioning data
    /// carries a j-gradient but NO i-gradient (isolates the collocated trend's areal
    /// effect along i).
    fn wells_j_varying() -> Vec<WellLog> {
        let south = WellLog::new(
            50.0,
            10.0,
            vec![(5.0, 0.15), (15.0, 0.16), (25.0, 0.14), (35.0, 0.15)],
        );
        let north = WellLog::new(
            50.0,
            90.0,
            vec![(5.0, 0.23), (15.0, 0.24), (25.0, 0.22), (35.0, 0.23)],
        );
        vec![south, north]
    }

    /// Least-squares slope of the cube's values against the column index i.
    fn slope_vs_i(values: &[f64]) -> f64 {
        let (mut sx, mut sy, mut sxy, mut sxx, mut n) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for k in 0..NK {
            for j in 0..NJ {
                for i in 0..NI {
                    let (x, y) = (i as f64, values[idx(i, j, k)]);
                    sx += x;
                    sy += y;
                    sxy += x * y;
                    sxx += x * x;
                    n += 1.0;
                }
            }
        }
        (sxy / n - (sx / n) * (sy / n)) / (sxx / n - (sx / n).powi(2))
    }

    #[test]
    fn collocated_corr_zero_is_a_bitwise_noop() {
        // corr = 0 collocated cokriging must reduce to plain SGS bit-for-bit (the
        // secondary decouples, the RNG stream is unchanged).
        let plain = {
            let mut g = test_grid();
            PropertyPipeline::new("PHIE")
                .upscale(two_wells(), UpscaleMethod::Arithmetic)
                .propagate(Gaussian::new(nscore_variogram(), 9))
                .apply(&mut g)
                .unwrap();
            g.properties().get("PHIE").unwrap().values.clone()
        };
        let with_zero = {
            let mut g = test_grid();
            PropertyPipeline::new("PHIE")
                .upscale(two_wells(), UpscaleMethod::Arithmetic)
                .propagate(
                    Gaussian::new(nscore_variogram(), 9).with_trend(i_increasing_trend(), 0.0),
                )
                .apply(&mut g)
                .unwrap();
            g.properties().get("PHIE").unwrap().values.clone()
        };
        assert_eq!(plain, with_zero, "corr=0 must equal plain SGS bit-for-bit");
    }

    #[test]
    fn collocated_world_trend_without_georef_is_a_loud_error() {
        // Finding 1: a WORLD-georeferenced trend on a LOCAL model (no model georef)
        // resamples to an all-NaN secondary (the world extent never overlaps the local
        // lattice). Before the fix the kernel silently dropped every node to plain SGS
        // (collocated cokriging a no-op); now sub-threshold coverage is a hard error, so
        // the misalignment surfaces instead of hiding.
        let mut g = test_grid();
        let world_trend =
            TrendSurface::new(NI, NJ, (0..NI * NJ).map(|k| (k % NI) as f64).collect())
                .unwrap()
                .with_georef(431_000.0, 6_521_000.0, 20.0, 20.0);
        let res = PropertyPipeline::new("PHIE")
            .upscale(two_wells(), UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(nscore_variogram(), 9).with_trend(world_trend, 0.6))
            .apply(&mut g);
        assert!(
            matches!(res, Err(StaticError::InvalidInput(_))),
            "world trend on a local model must error, not silently no-op: {res:?}"
        );
    }

    #[test]
    fn collocated_positive_corr_shifts_the_lateral_pattern() {
        // Conditioning has no i-gradient; a trend increasing with i at corr>0 must
        // steepen the field's i-slope relative to plain SGS.
        let run = |trend: Option<(TrendSurface, f64)>| {
            let mut g = test_grid();
            let mut gauss = Gaussian::new(nscore_variogram(), 4);
            if let Some((t, c)) = trend {
                gauss = gauss.with_trend(t, c);
            }
            PropertyPipeline::new("PHIE")
                .upscale(wells_j_varying(), UpscaleMethod::Arithmetic)
                .propagate(gauss)
                .apply(&mut g)
                .unwrap();
            g.properties().get("PHIE").unwrap().values.clone()
        };
        let plain_slope = slope_vs_i(&run(None));
        let co_slope = slope_vs_i(&run(Some((i_increasing_trend(), 0.9))));
        assert!(
            co_slope > plain_slope,
            "collocated trend should steepen the i-slope: plain {plain_slope} vs co {co_slope}"
        );
        assert!(co_slope > 0.0, "co-slope should track the increasing trend");
    }

    // --- bounded SGS default + loud conditioning errors (items 4/5/6) ---

    /// A single well conditioning only layer 0 (tvd 5). Layers 1..3 of a 0..4 pipe
    /// then have no data.
    fn one_layer_well() -> Vec<WellLog> {
        vec![WellLog::new(10.0, 10.0, vec![(5.0, 0.2)])]
    }

    #[test]
    fn data_less_layer_is_a_loud_named_error_by_default() {
        // Item 4: a simulated layer with no conditioning data used to be silently
        // filled with the conditioned mean (structureless). It is now a hard error
        // naming the property; `allow_mean_fill` opts back into the fill.
        let mut grid = test_grid();
        let res = PropertyPipeline::new("PHIE")
            .upscale(one_layer_well(), UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(nscore_variogram(), 1))
            .apply(&mut grid);
        match res {
            Err(StaticError::InvalidInput(m)) => {
                assert!(m.contains("PHIE"), "error names the property: {m}");
                assert!(m.contains("no conditioning data"), "error explains: {m}");
            }
            other => panic!("expected a named InvalidInput, got {other:?}"),
        }
        // Opt-in succeeds and mean-fills the data-less layers.
        let mut grid = test_grid();
        let r = PropertyPipeline::new("PHIE")
            .upscale(one_layer_well(), UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(nscore_variogram(), 1).allow_mean_fill())
            .apply(&mut grid);
        assert!(r.is_ok(), "allow_mean_fill opts into the fill: {r:?}");
        assert!(grid
            .properties()
            .get("PHIE")
            .unwrap()
            .values
            .iter()
            .all(|v| v.is_finite()));
    }

    #[test]
    fn no_conditioning_error_names_the_property() {
        // Item 6: a well that misses the grid entirely errors with the property name
        // (previously an anonymous "propagate: no conditioning data").
        let mut grid = test_grid();
        let off = WellLog::new(1e6, 1e6, vec![(5.0, 0.2)]);
        let res = PropertyPipeline::new("KLOGH")
            .upscale(vec![off], UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(nscore_variogram(), 1))
            .apply(&mut grid);
        match res {
            Err(StaticError::InvalidInput(m)) => {
                assert!(m.contains("KLOGH"), "names property: {m}")
            }
            other => panic!("expected named error, got {other:?}"),
        }
    }

    #[test]
    fn bounded_search_is_the_default_and_matches_unbounded() {
        // Item 5: the default search is now BOUNDED (a variogram-range/lattice-scaled
        // radius). On this fixture the variogram range (50) is within the default
        // radius, so beyond-range nodes carry ~zero kriging weight and the bounded
        // field matches the (opt-in) unbounded field to tight tolerance, while the
        // wells stay honoured exactly.
        let run = |gauss: Gaussian| {
            let mut g = test_grid();
            PropertyPipeline::new("PHIE")
                .upscale(two_wells(), UpscaleMethod::Arithmetic)
                .propagate(gauss)
                .apply(&mut g)
                .unwrap();
            g.properties().get("PHIE").unwrap().values.clone()
        };
        let bounded = run(Gaussian::new(nscore_variogram(), 11));
        let unbounded = run(Gaussian::new(nscore_variogram(), 11).with_unbounded_search());
        // Wells honoured in both (SGS reproduces data nodes regardless of search).
        for k in 0..NK {
            assert!((bounded[idx(0, 0, k)] - [0.10, 0.12, 0.14, 0.16][k]).abs() < 1e-6);
        }
        // Fields agree to a tight tolerance (max abs diff over the whole cube).
        let maxdiff = bounded
            .iter()
            .zip(&unbounded)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(
            maxdiff < 5e-3,
            "bounded default should match unbounded within tolerance: maxdiff {maxdiff}"
        );
    }
}
