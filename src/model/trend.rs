//! [`TrendSurface`] — the minimal areal-trend hook (external-drift-lite).
//!
//! A gridded areal multiplier field that shapes property population *laterally*
//! without touching the *level*: the field is resampled to the model column
//! lattice, normalized around its own mean, and applied per-column to the NTG
//! cube (and, optionally, porosity). This honours the ratified
//! logs-give-shape / draw-gives-level principle
//! (`decision_staticmodel_regen_seam`): the trend supplies lateral shape, the
//! prior/draw supplies the level. Because the multipliers are mean-normalized,
//! the field-mean of the populated property is preserved.
//!
//! The full geostatistical treatment (variograms, kriging-with-trend, SGS)
//! stays P5 (`task_petekstatic_property_modelling`); this is the deterministic
//! MVP a real appraisal project needs today.
//!
//! **Resampling (2026-07-04): the shared kernel.** The hand-rolled nearest-node
//! sampler was retired onto [`petektools::resample`] (bilinear, null-aware) — the
//! one-resampler-in-the-family rule. A [`TrendSurface::with_georef`] trend resamples
//! by the model column's **world** `(x, y)` (from the grid cell centroids); a bare
//! trend is stretched across the model extent (index-space alignment, correct only
//! when the frames coincide). Nodes outside the source extent — and `NaN` trend
//! nodes — fall back to a multiplier of `1.0`. The same `resample_to` also serves
//! the trend as a **collocated cokriging secondary** in the P5 property pipeline
//! ([`crate::model::Gaussian::with_trend`]).

use crate::error::StaticError;
use ndarray::Array2;
use petektools::{resample, Lattice, ResampleMethod};

/// World georeference for a trend field: the world `(x, y)` of node `(0, 0)`'s
/// centre and the node spacing, so a model column's world coordinate maps to the
/// nearest trend node (R4).
#[derive(Debug, Clone, Copy, PartialEq)]
struct Georef {
    origin_x: f64,
    origin_y: f64,
    node_dx: f64,
    node_dy: f64,
}

/// A gridded areal multiplier field on its own `ncol × nrow` lattice (row-major,
/// `NaN` = undefined → multiplier `1.0`). Resampled + mean-normalized to the
/// model column lattice at population time.
#[derive(Debug, Clone, PartialEq)]
pub struct TrendSurface {
    ncol: usize,
    nrow: usize,
    /// Row-major multiplier values; `NaN` = undefined.
    values: Vec<f64>,
    /// Apply the trend to porosity in addition to NTG (NTG is always targeted).
    apply_porosity: bool,
    /// Optional world georeference (`None` = index-space fraction resampling).
    georef: Option<Georef>,
}

impl TrendSurface {
    /// Build an areal trend field of `ncol × nrow` row-major multiplier values.
    /// `NaN` marks an undefined node (falls back to `1.0`); every finite value
    /// must be non-negative (a multiplier).
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the lattice is empty, `values.len()` !=
    /// `ncol*nrow`, or a finite value is negative.
    pub fn new(ncol: usize, nrow: usize, values: Vec<f64>) -> Result<Self, StaticError> {
        if ncol == 0 || nrow == 0 {
            return Err(StaticError::InvalidInput(format!(
                "trend lattice must be non-empty, got {ncol}x{nrow}"
            )));
        }
        if values.len() != ncol * nrow {
            return Err(StaticError::InvalidInput(format!(
                "trend has {} values, expected {}",
                values.len(),
                ncol * nrow
            )));
        }
        if let Some(bad) = values.iter().find(|v| v.is_finite() && **v < 0.0) {
            return Err(StaticError::InvalidInput(format!(
                "trend multiplier must be non-negative, got {bad}"
            )));
        }
        Ok(Self {
            ncol,
            nrow,
            values,
            apply_porosity: false,
            georef: None,
        })
    }

    /// Also apply this trend to the porosity cube (default: NTG only).
    #[must_use]
    pub fn with_porosity(mut self) -> Self {
        self.apply_porosity = true;
        self
    }

    /// Georeference the field to world coordinates (R4): `origin_x`/`origin_y` are
    /// the world `(x, y)` of node `(0, 0)`'s centre, `node_dx`/`node_dy` its node
    /// spacing. Population then resamples by the model column's **world**
    /// coordinate (nearest node) rather than by index fraction, so the trend is
    /// not aligned to the model lattice by luck. Non-positive spacing is ignored
    /// (stays index-space).
    #[must_use]
    pub fn with_georef(mut self, origin_x: f64, origin_y: f64, node_dx: f64, node_dy: f64) -> Self {
        if node_dx.is_finite() && node_dx > 0.0 && node_dy.is_finite() && node_dy > 0.0 {
            self.georef = Some(Georef {
                origin_x,
                origin_y,
                node_dx,
                node_dy,
            });
        }
        self
    }

    /// Whether this trend carries a world georeference (R4).
    #[must_use]
    pub fn is_georeferenced(&self) -> bool {
        self.georef.is_some()
    }

    /// Whether porosity is also modulated by this trend.
    #[must_use]
    pub fn applies_to_porosity(&self) -> bool {
        self.apply_porosity
    }

    /// Resample the raw field onto a target [`Lattice`] via the shared
    /// [`petektools::resample`] kernel (bilinear, null-aware) — returning the
    /// `(target.ncol × target.nrow)` field indexed `[[i, j]]`, `NaN` outside the
    /// source extent. This is the trend as a **collocated secondary** for
    /// [`crate::model::Gaussian::with_trend`] (SGS cokriging standardizes it internally, so
    /// no mean-normalization here — unlike the multiplier path).
    ///
    /// A georeferenced trend ([`TrendSurface::with_georef`]) resamples by world
    /// coordinate; otherwise the field is stretched across the target's extent
    /// (index-space alignment — correct only when the frames coincide).
    ///
    /// # Errors
    /// [`StaticError::Grid`] if the shared resampler rejects the geometry.
    pub(crate) fn resample_to(&self, target: &Lattice) -> Result<Array2<f64>, StaticError> {
        // Row-major `values[r*ncol+c]` -> Array2 shape (ncol, nrow) indexed [[c, r]].
        let src = Array2::from_shape_fn((self.ncol, self.nrow), |(c, r)| {
            self.values[r * self.ncol + c]
        });
        let src_georef = match self.georef {
            Some(g) => Lattice::regular(
                g.origin_x, g.origin_y, g.node_dx, g.node_dy, self.ncol, self.nrow,
            ),
            None => {
                // Index-space: stretch the source across the target's extent so node
                // (0,0) lands on the target origin and (ncol-1, nrow-1) on its far node.
                let o = target.node_xy(0, 0);
                let far = target.node_xy(target.ncol - 1, target.nrow - 1);
                let dx = span_inc(o.0, far.0, self.ncol);
                let dy = span_inc(o.1, far.1, self.nrow);
                Lattice::regular(o.0, o.1, dx, dy, self.ncol, self.nrow)
            }
        };
        resample(&src, &src_georef, target, ResampleMethod::Bilinear)
            .map_err(|e| StaticError::Grid(format!("trend resample failed: {e}")))
    }

    /// Mean-normalize a raw multiplier field around the mean of its *defined*
    /// nodes (so the property field-mean is preserved); undefined → `1.0`.
    fn normalize(raw: Vec<f64>) -> Vec<f64> {
        let defined: Vec<f64> = raw.iter().copied().filter(|v| v.is_finite()).collect();
        let mean = if defined.is_empty() {
            1.0
        } else {
            defined.iter().sum::<f64>() / defined.len() as f64
        };
        raw.into_iter()
            .map(|v| {
                if v.is_finite() && mean > 0.0 {
                    v / mean
                } else {
                    1.0
                }
            })
            .collect()
    }

    /// The per-column multiplier field resampled to a model areal `lattice` and
    /// mean-normalized (undefined columns → `1.0`), indexed row-major `j * ni + i`.
    ///
    /// The resample rides the **shared** [`petektools::resample`] kernel via
    /// [`TrendSurface::resample_to`] (the one-resampler-in-the-family rule — the
    /// hand-rolled nearest-node sampler was retired 2026-07-04): a georeferenced
    /// trend samples by the lattice nodes' world positions, a non-georeferenced one
    /// is stretched across the lattice extent.
    ///
    /// # Errors
    /// [`StaticError::Grid`] if the shared resampler rejects the geometry.
    pub(crate) fn column_multipliers_on(&self, lattice: &Lattice) -> Result<Vec<f64>, StaticError> {
        let field = self.resample_to(lattice)?;
        let (ni, nj) = (lattice.ncol, lattice.nrow);
        let mut raw = vec![f64::NAN; ni * nj];
        for j in 0..nj {
            for i in 0..ni {
                raw[j * ni + i] = field[[i, j]];
            }
        }
        Ok(Self::normalize(raw))
    }
}

/// Node spacing that spreads `n` nodes across `[start, end]`. A single-node or
/// degenerate span falls back to unit spacing (keeps [`Lattice::regular`] valid).
fn span_inc(start: f64, end: f64, n: usize) -> f64 {
    if n < 2 {
        return 1.0;
    }
    let inc = (end - start) / (n - 1) as f64;
    if inc.is_finite() && inc > 0.0 {
        inc
    } else {
        1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_shape_and_negatives() {
        assert!(TrendSurface::new(0, 3, vec![]).is_err());
        assert!(TrendSurface::new(2, 2, vec![1.0; 3]).is_err()); // wrong len
        assert!(TrendSurface::new(2, 1, vec![1.0, -0.5]).is_err()); // negative
        assert!(TrendSurface::new(2, 1, vec![1.0, f64::NAN]).is_ok()); // NaN allowed
    }

    #[test]
    fn multipliers_are_mean_normalized() {
        // A 2x1 field {1, 3} resampled (via the shared kernel) 1:1 to a 2x1 column
        // lattice -> mean 2 -> multipliers {0.5, 1.5}, averaging to exactly 1.0.
        let t = TrendSurface::new(2, 1, vec![1.0, 3.0]).unwrap();
        let m = t
            .column_multipliers_on(&Lattice::regular(0.0, 0.0, 1.0, 1.0, 2, 1))
            .unwrap();
        assert!((m[0] - 0.5).abs() < 1e-12, "{m:?}");
        assert!((m[1] - 1.5).abs() < 1e-12, "{m:?}");
        assert!(((m[0] + m[1]) / 2.0 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn georef_resamples_by_world_coordinate() {
        // 3x1 nodes at world x = 0, 100, 200 (node_dx 100), values {1,2,3}. The
        // shared resampler samples the target lattice's node world positions.
        let t = TrendSurface::new(3, 1, vec![1.0, 2.0, 3.0])
            .unwrap()
            .with_georef(0.0, 0.0, 100.0, 100.0);
        assert!(t.is_georeferenced());
        // A target lattice whose nodes coincide with the trend nodes.
        let m = t
            .column_multipliers_on(&Lattice::regular(0.0, 0.0, 100.0, 100.0, 3, 1))
            .unwrap();
        // Raw {1,2,3}, mean 2 -> normalized {0.5, 1.0, 1.5}.
        assert!((m[0] - 0.5).abs() < 1e-12, "{m:?}");
        assert!((m[1] - 1.0).abs() < 1e-12, "{m:?}");
        assert!((m[2] - 1.5).abs() < 1e-12, "{m:?}");
    }

    #[test]
    fn georef_out_of_extent_falls_back_to_unity() {
        // Trend covers world x in [0, 100] (nodes {1, 3}); a target node beyond that
        // resamples to NaN (no extrapolation) -> multiplier 1.0. The two in-extent
        // nodes define the mean (2) -> {0.5, 1.5}; the third (x=200) is unity.
        let t = TrendSurface::new(2, 1, vec![1.0, 3.0])
            .unwrap()
            .with_georef(0.0, 0.0, 100.0, 100.0);
        let m = t
            .column_multipliers_on(&Lattice::regular(0.0, 0.0, 100.0, 100.0, 3, 1))
            .unwrap();
        assert!((m[0] - 0.5).abs() < 1e-12, "in-extent {m:?}");
        assert!((m[1] - 1.5).abs() < 1e-12, "in-extent {m:?}");
        assert!((m[2] - 1.0).abs() < 1e-12, "out-of-extent -> unity {m:?}");
    }

    #[test]
    fn georef_off_by_default_and_ignores_bad_spacing() {
        let t = TrendSurface::new(2, 1, vec![1.0, 2.0]).unwrap();
        assert!(!t.is_georeferenced());
        // Non-positive spacing does not georeference.
        assert!(!t
            .clone()
            .with_georef(0.0, 0.0, 0.0, 100.0)
            .is_georeferenced());
    }

    #[test]
    fn undefined_nodes_fall_back_to_unity() {
        let t = TrendSurface::new(2, 1, vec![2.0, f64::NAN]).unwrap();
        let m = t
            .column_multipliers_on(&Lattice::regular(0.0, 0.0, 1.0, 1.0, 2, 1))
            .unwrap();
        // Only node 0 defines the mean (2.0) -> its multiplier is 1.0; the NaN
        // column falls back to 1.0 too.
        assert!((m[0] - 1.0).abs() < 1e-12, "{m:?}");
        assert!((m[1] - 1.0).abs() < 1e-12, "{m:?}");
    }
}
