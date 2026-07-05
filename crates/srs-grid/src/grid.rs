//! The `Grid`: corner-point geometry + cell properties + layering. The single
//! spine type every downstream crate (gridder, volumetrics, renderer) operates
//! on.

use crate::cell::Cell;
use crate::geometry::{CornerPointGeom, Pillar};
use crate::index::{Dims, Ijk};
use crate::klayer::KLayer;
use crate::property::Properties;

/// A reservoir grid: geometry, per-cell properties, and stratigraphic layers.
#[derive(Debug, Clone)]
pub struct Grid {
    geom: CornerPointGeom,
    properties: Properties,
    layers: Vec<KLayer>,
}

impl Grid {
    /// Wrap a corner-point geometry in a grid with empty properties and one
    /// `KLayer` per `k`.
    #[must_use]
    pub fn new(geom: CornerPointGeom) -> Self {
        let dims = geom.dims();
        let layers = (0..dims.nk).map(KLayer::new).collect();
        let properties = Properties::new(dims.cell_count());
        Self {
            geom,
            properties,
            layers,
        }
    }

    /// Take this grid's geometry buffers **out** for in-place refill — the
    /// allocation-recycling entry for `StaticModelTemplate::realize_into`. The
    /// property cube allocations are left untouched (recycle them separately via
    /// [`Properties::take_values`]); reinstall the refilled geometry with
    /// [`Grid::install_geometry`]. The grid is invalid until then.
    #[must_use]
    pub fn take_geometry_buffers(&mut self) -> (Vec<Pillar>, Vec<f64>) {
        self.geom.take_buffers()
    }

    /// Reinstall refilled geometry under `dims`: rebuild the k-layers **only if**
    /// the layer count changed (else the layer Strings are kept, saving `nk`
    /// allocations per draw), and update the property cell count so the following
    /// populate refills the recycled cubes to the new length. Bit-identical to a
    /// grid freshly built from the same `(dims, coord, zcorn)`.
    pub fn install_geometry(&mut self, dims: Dims, coord: Vec<Pillar>, zcorn: Vec<f64>) {
        let nk = dims.nk;
        let cell_count = dims.cell_count();
        self.geom.install(dims, coord, zcorn);
        if self.layers.len() != nk {
            self.layers = (0..nk).map(KLayer::new).collect();
        }
        self.properties.set_cell_count(cell_count);
    }

    /// Grid dimensions.
    #[must_use]
    pub fn dims(&self) -> Dims {
        self.geom.dims()
    }

    /// Borrow the corner-point geometry (COORD + ZCORN + the vertical flag) — the
    /// read side the out-of-core spill path needs to persist a built grid to a
    /// petekTools store lane family.
    #[must_use]
    pub fn geometry(&self) -> &CornerPointGeom {
        &self.geom
    }

    /// Total cell count.
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.geom.dims().cell_count()
    }

    /// The cell view at `(i,j,k)`.
    ///
    /// # Panics
    /// Panics if `c` is out of bounds.
    #[must_use]
    #[inline]
    pub fn cell(&self, c: Ijk) -> Cell {
        Cell {
            ijk: c,
            corners: self.geom.cell_corners(c),
        }
    }

    /// Iterate every cell in linear order.
    pub fn cells(&self) -> impl Iterator<Item = Cell> + '_ {
        self.geom.dims().iter().map(|c| self.cell(c))
    }

    /// Centroid depth (`z`) of the cell at **linear** index `lin`, read straight
    /// from the corner depths — the cheap contact test for volumetrics, avoiding
    /// the full corner+`xy` build of [`Grid::cell`]. `lin` must be `< cell_count`.
    /// Equals `self.cell(c).centroid().z` for the matching `(i,j,k)`.
    #[must_use]
    #[inline]
    pub fn cell_centroid_z_at(&self, lin: usize) -> f64 {
        self.geom.centroid_z_at(lin)
    }

    /// Bulk volume \[m³\] of the cell at **linear** index `lin` — the cheap
    /// `area × mean_dz` path on a vertical lattice (every grid today), else the
    /// general hexahedron split. Equals `self.cell(c).volume()` for the matching
    /// `(i,j,k)`. `lin` must be `< cell_count`.
    #[must_use]
    #[inline]
    pub fn cell_volume_at(&self, lin: usize) -> f64 {
        self.geom.cell_volume_at(lin)
    }

    /// Total bulk (gross rock) volume in m³.
    #[must_use]
    pub fn bulk_volume(&self) -> f64 {
        (0..self.cell_count())
            .map(|i| self.geom.cell_volume_at(i))
            .sum()
    }

    /// Shared cell properties.
    #[must_use]
    pub fn properties(&self) -> &Properties {
        &self.properties
    }

    /// Mutable cell properties.
    pub fn properties_mut(&mut self) -> &mut Properties {
        &mut self.properties
    }

    /// The stratigraphic layers (one per `k`).
    #[must_use]
    pub fn layers(&self) -> &[KLayer] {
        &self.layers
    }
}
