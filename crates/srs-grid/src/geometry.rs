//! Corner-point geometry: a pillar lattice (COORD) plus per-cell corner depths
//! (ZCORN). The box grid is the degenerate case (vertical, uniform pillars);
//! the same structure carries a faulted/tilted grid later, so the box upgrades
//! in place (Ponting 1989).

use crate::index::{Dims, Ijk};
use crate::point::Point3;

/// A coordinate pillar: a straight line between a top and bottom anchor. Cell
/// corners ride on the pillar at their ZCORN depth (vertical for a box).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pillar {
    pub top: Point3,
    pub bottom: Point3,
}

impl Pillar {
    /// `(x,y)` on the pillar at depth `z`, interpolated between the anchors.
    /// Degenerate (top.z == bottom.z) pillars return the top `(x,y)`.
    #[inline]
    fn xy_at(self, z: f64) -> (f64, f64) {
        let dz = self.bottom.z - self.top.z;
        if dz.abs() < f64::EPSILON {
            return (self.top.x, self.top.y);
        }
        let t = (z - self.top.z) / dz;
        (
            self.top.x + t * (self.bottom.x - self.top.x),
            self.top.y + t * (self.bottom.y - self.top.y),
        )
    }
}

/// Pillar lattice + corner depths for every cell.
///
/// `zcorn` is `cell_count * 8`, the 8 depths per cell ordered by local
/// `(di,dj,dk)` bits (`di + 2*dj + 4*dk`), matching the `volume` corner order.
#[derive(Debug, Clone)]
pub struct CornerPointGeom {
    dims: Dims,
    coord: Vec<Pillar>,
    zcorn: Vec<f64>,
    /// `true` **iff every pillar is vertical** (`top.xy == bottom.xy`) — the case
    /// for every box / layered grid today. On a vertical pillar `xy_at` is
    /// z-independent (the interpolation term is `t * 0.0 == 0.0`), so a corner's
    /// `xy` is exactly its pillar's `top.xy`; the volume hot path then skips the
    /// per-corner division while staying **bit-identical** to the general
    /// [`CornerPointGeom::cell_corners`] build. Computed once at construction (the
    /// geometry is immutable thereafter). `false` → the general per-corner path.
    vertical: bool,
}

impl CornerPointGeom {
    /// Build from a pillar lattice and a flat ZCORN array.
    ///
    /// # Panics
    /// Panics if `coord`/`zcorn` lengths disagree with `dims` (internal
    /// invariant; constructors in this crate always satisfy it).
    #[must_use]
    pub fn new(dims: Dims, coord: Vec<Pillar>, zcorn: Vec<f64>) -> Self {
        assert_eq!(coord.len(), dims.pillar_count(), "coord length mismatch");
        assert_eq!(zcorn.len(), dims.cell_count() * 8, "zcorn length mismatch");
        // A vertical lattice (every grid today) lets the volume hot path build
        // corners without the per-corner `xy` interpolation division — bit-for-bit
        // the same corners, since a vertical pillar's `xy_at` returns `top.xy`.
        let vertical = coord
            .iter()
            .all(|p| p.top.x == p.bottom.x && p.top.y == p.bottom.y);
        Self {
            dims,
            coord,
            zcorn,
            vertical,
        }
    }

    /// Bulk volume of the cell at **linear** index `lin`. On a vertical lattice
    /// (every grid today) this builds the 8 corners straight from the pillar tops +
    /// ZCORN — skipping the per-corner `xy` interpolation division — then runs the
    /// same [`hexahedron_volume`] split, so the result is **bit-identical** to
    /// `Grid::cell(c).volume()`. On a non-vertical lattice it takes the general
    /// [`CornerPointGeom::cell_corners`] path. `lin` must be `< cell_count`.
    #[must_use]
    #[inline]
    pub fn cell_volume_at(&self, lin: usize) -> f64 {
        crate::volume::hexahedron_volume(&self.cell_corners_at(lin))
    }

    /// The 8 corners of the cell at **linear** index `lin`. On a vertical lattice
    /// the `xy` come free from the pillar tops (no interpolation); otherwise this
    /// delegates to the general [`CornerPointGeom::cell_corners`].
    #[inline]
    fn cell_corners_at(&self, lin: usize) -> [Point3; 8] {
        let per_layer = self.dims.ni * self.dims.nj;
        let layer = lin / per_layer;
        let rem = lin % per_layer;
        let j = rem / self.dims.ni;
        let i = rem % self.dims.ni;
        if !self.vertical {
            return self.cell_corners(Ijk::new(i, j, layer));
        }
        let z = &self.zcorn[lin * 8..lin * 8 + 8];
        let mut corners = [Point3::new(0.0, 0.0, 0.0); 8];
        for (idx, slot) in corners.iter_mut().enumerate() {
            let di = idx & 1;
            let dj = (idx >> 1) & 1;
            let top = self.coord[self.dims.pillar_linear(i + di, j + dj)].top;
            *slot = Point3::new(top.x, top.y, z[idx]);
        }
        corners
    }

    /// Take the geometry buffers **out** for in-place refill, leaving the geometry
    /// empty (`coord`/`zcorn` emptied, capacity retained). The allocation-recycling
    /// entry for `StaticModelTemplate::realize_into`: a gridder refills the returned
    /// `(coord, zcorn)` and [`CornerPointGeom::install`] puts them back. The geometry
    /// is invalid (zero-length buffers) between the two calls — install before any
    /// read.
    #[must_use]
    pub fn take_buffers(&mut self) -> (Vec<Pillar>, Vec<f64>) {
        (
            std::mem::take(&mut self.coord),
            std::mem::take(&mut self.zcorn),
        )
    }

    /// Reinstall refilled geometry buffers under `dims`, recomputing the `vertical`
    /// fast-path flag. Same length invariants as [`CornerPointGeom::new`], so the
    /// recycled geometry is **bit-identical** to a freshly built one.
    ///
    /// # Panics
    /// Panics if `coord`/`zcorn` lengths disagree with `dims` (the constructor
    /// invariant).
    pub fn install(&mut self, dims: Dims, coord: Vec<Pillar>, zcorn: Vec<f64>) {
        assert_eq!(coord.len(), dims.pillar_count(), "coord length mismatch");
        assert_eq!(zcorn.len(), dims.cell_count() * 8, "zcorn length mismatch");
        let vertical = coord
            .iter()
            .all(|p| p.top.x == p.bottom.x && p.top.y == p.bottom.y);
        self.dims = dims;
        self.coord = coord;
        self.zcorn = zcorn;
        self.vertical = vertical;
    }

    /// The grid dimensions.
    #[must_use]
    pub fn dims(&self) -> Dims {
        self.dims
    }

    /// Borrow the pillar lattice (COORD) — the small, k-invariant areal frame.
    /// Used by the out-of-core spill path to persist COORD as a flat store lane
    /// (it stays resident even for a spilled model; only ZCORN/cubes spill).
    #[must_use]
    pub fn coord(&self) -> &[Pillar] {
        &self.coord
    }

    /// Borrow the flat ZCORN array (`cell_count * 8`, k-slab-major) — the big
    /// per-cell corner-depth lane the out-of-core spill path streams slab-by-slab
    /// into a petekTools store (as f32 at spill scale, ruling R4).
    #[must_use]
    pub fn zcorn(&self) -> &[f64] {
        &self.zcorn
    }

    /// Whether every pillar is vertical (the case for every box / layered grid
    /// today). The spill path requires it: corner `xy` come straight from the
    /// pillar tops, so a spilled cell rebuilds bit-parallel to the in-core
    /// [`CornerPointGeom::cell_corners_at`] fast path (only ZCORN narrows to f32).
    #[must_use]
    pub fn is_vertical(&self) -> bool {
        self.vertical
    }

    /// Mean corner depth (centroid `z`) of the cell at **linear** index `lin`,
    /// read straight from ZCORN — no pillar `xy` interpolation, no corner
    /// materialization. The cheap contact test on the volumetrics hot path: a
    /// below-contact cell is rejected on 8 depth reads without ever building its
    /// corners or volume. `lin` must be a valid cell index (`< cell_count`).
    ///
    /// The plain mean of the 8 corner depths equals `centroid().z` because
    /// [`Cell::centroid`] averages the same 8 corners; the `xy` interpolation the
    /// full corner build performs never moves a corner's `z` (a corner rides its
    /// pillar *at* its ZCORN depth).
    #[must_use]
    #[inline]
    pub fn centroid_z_at(&self, lin: usize) -> f64 {
        let z = &self.zcorn[lin * 8..lin * 8 + 8];
        (z[0] + z[1] + z[2] + z[3] + z[4] + z[5] + z[6] + z[7]) / 8.0
    }

    /// The 8 corners of cell `(i,j,k)`, ordered `di + 2*dj + 4*dk`.
    ///
    /// # Panics
    /// Panics if `c` is out of bounds.
    #[must_use]
    #[inline]
    pub fn cell_corners(&self, c: Ijk) -> [Point3; 8] {
        let cell = self.dims.linear(c).expect("cell index out of bounds");
        let z = &self.zcorn[cell * 8..cell * 8 + 8];
        let mut corners = [Point3::new(0.0, 0.0, 0.0); 8];
        for (idx, slot) in corners.iter_mut().enumerate() {
            // Areal pillar from the di/dj bits; dk (bit 2) selects top/bottom
            // depth, which is already encoded in zcorn[idx].
            let di = idx & 1;
            let dj = (idx >> 1) & 1;
            let pillar = self.coord[self.dims.pillar_linear(c.i + di, c.j + dj)];
            let zc = z[idx];
            let (x, y) = pillar.xy_at(zc);
            *slot = Point3::new(x, y, zc);
        }
        corners
    }
}
