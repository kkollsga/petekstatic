//! A single cell: its 8 corners and the geometric quantities volumetrics needs
//! (volume, centroid, top/bottom depth). Built on demand from the grid geometry.
//! All lengths in metres; depth is positive-down.

use crate::grid::index::Ijk;
use crate::grid::point::Point3;
use crate::grid::volume::hexahedron_volume;

/// A cell view: its index and 8 corners (order `di + 2*dj + 4*dk`).
#[derive(Debug, Clone, Copy)]
pub struct Cell {
    pub ijk: Ijk,
    pub corners: [Point3; 8],
}

impl Cell {
    /// Bulk (gross) volume in m³.
    #[must_use]
    #[inline]
    pub fn volume(&self) -> f64 {
        hexahedron_volume(&self.corners)
    }

    /// Geometric centroid (mean of the 8 corners).
    #[must_use]
    #[inline]
    pub fn centroid(&self) -> Point3 {
        let s = self
            .corners
            .iter()
            .fold(Point3::new(0.0, 0.0, 0.0), |a, &c| a + c);
        Point3::new(s.x / 8.0, s.y / 8.0, s.z / 8.0)
    }

    /// Mean depth of the four top corners (k face, smaller z).
    #[must_use]
    #[inline]
    pub fn top_depth(&self) -> f64 {
        (self.corners[0].z + self.corners[1].z + self.corners[2].z + self.corners[3].z) / 4.0
    }

    /// Mean depth of the four bottom corners (k+1 face, larger z).
    #[must_use]
    #[inline]
    pub fn bottom_depth(&self) -> f64 {
        (self.corners[4].z + self.corners[5].z + self.corners[6].z + self.corners[7].z) / 4.0
    }

    /// Mean vertical thickness (bottom - top depth).
    #[must_use]
    #[inline]
    pub fn dz(&self) -> f64 {
        self.bottom_depth() - self.top_depth()
    }
}
