//! Shared view primitives: the areal [`GridFrame`] georeference, a georeferenced
//! [`ScalarLayer`], and the [`ValueRange`] legend metadata.

use crate::model::Georef;
use crate::pipeline::areal_lattice;
use petekstatic_error::StaticError;
use serde::{Deserialize, Serialize};
use srs_grid::Grid;

/// The value span of a layer over its **finite** entries — the metadata a viewer
/// needs to build a colour-ramp legend. Both fields are `NaN` when the layer has
/// no finite value.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ValueRange {
    /// Minimum finite value (`NaN` if none).
    pub min: f64,
    /// Maximum finite value (`NaN` if none).
    pub max: f64,
}

impl ValueRange {
    /// The `[min, max]` span over the finite values in `it` (`NaN`/`inf` skipped).
    pub(crate) fn of(it: impl Iterator<Item = f64>) -> Self {
        let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
        for v in it {
            if v.is_finite() {
                min = min.min(v);
                max = max.max(v);
            }
        }
        if min <= max {
            Self { min, max }
        } else {
            Self {
                min: f64::NAN,
                max: f64::NAN,
            }
        }
    }
}

/// The areal georeference all map layers share: a regular, axis-aligned lattice
/// of the grid's **column centroids** (the same `xy↔ij` frame the property
/// pipeline + well registration use). `origin_*` is the world `(x, y)` of node
/// `(0, 0)` (column `(0, 0)`'s centroid); `spacing_*` the node spacing;
/// `ncol == ni`, `nrow == nj`. World `(x, y)` of node `(i, j)` is
/// `(origin_x + i * spacing_x, origin_y + j * spacing_y)`.
///
/// When the model carries a registered world [`Georef`] the frame is that world
/// lattice (so the raster overlays the world outline / wells and a world
/// fence/bore section traces through it); with no georeference it degenerates to
/// the grid's own local column-centroid lattice (the synthetic square / box).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GridFrame {
    pub origin_x: f64,
    pub origin_y: f64,
    pub spacing_x: f64,
    pub spacing_y: f64,
    /// Node count along x (== grid `ni`).
    pub ncol: usize,
    /// Node count along y (== grid `nj`).
    pub nrow: usize,
}

impl GridFrame {
    /// The areal frame for a model's grid. With a registered world [`Georef`] the
    /// frame is that world column lattice (world origin + spacing, `ncol`/`nrow`
    /// from the grid dims) — so raster, outline, wells and world fence/bore
    /// sections all speak ONE frame. Without one it degenerates to the grid's own
    /// local column-centroid lattice (requires an at-least-`2x2`, regular,
    /// axis-aligned column layout — the box / conformable case).
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the column lattice is smaller than `2x2`
    /// or (local case) not regular/axis-aligned.
    pub(crate) fn of_grid(grid: &Grid, georef: Option<Georef>) -> Result<Self, StaticError> {
        if let Some(g) = georef {
            let dims = grid.dims();
            let (ncol, nrow) = (dims.ni, dims.nj);
            if ncol < 2 || nrow < 2 {
                return Err(StaticError::InvalidInput(format!(
                    "view frame needs an areal lattice of at least 2x2 columns, got {ncol}x{nrow}"
                )));
            }
            return Ok(Self {
                origin_x: g.origin_x,
                origin_y: g.origin_y,
                spacing_x: g.spacing_x,
                spacing_y: g.spacing_y,
                ncol,
                nrow,
            });
        }
        let lat = areal_lattice(grid)?;
        Ok(Self {
            origin_x: lat.xori,
            origin_y: lat.yori,
            spacing_x: lat.xinc,
            spacing_y: lat.yinc,
            ncol: lat.ncol,
            nrow: lat.nrow,
        })
    }
}

/// A named, georeferenced scalar field on the shared [`GridFrame`], row-major
/// `values[j * ncol + i]`, `NaN` where undefined. Carries its unit label and the
/// finite value span for a legend.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScalarLayer {
    /// Layer name (property name, or a horizon/surface label).
    pub name: String,
    /// SI unit label for the values (e.g. `"m"`, `"fraction"`).
    pub units: String,
    /// Row-major `j * ncol + i` values (`NaN` = undefined / outside boundary).
    pub values: Vec<f64>,
    /// Value span over the finite entries (legend metadata).
    pub range: ValueRange,
}

impl ScalarLayer {
    /// Assemble a layer, computing its [`ValueRange`] from the values.
    pub(crate) fn new(name: impl Into<String>, units: impl Into<String>, values: Vec<f64>) -> Self {
        let range = ValueRange::of(values.iter().copied());
        Self {
            name: name.into(),
            units: units.into(),
            values,
            range,
        }
    }
}
