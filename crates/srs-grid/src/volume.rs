//! Hexahedral cell volume.
//!
//! Corners are ordered by local `(di,dj,dk)` bits as `di + 2*dj + 4*dk`, i.e.
//! `corners[0]` = `(i,  j,  k)`, `corners[7]` = `(i+1, j+1, k+1)`. The cell is
//! split into six tetrahedra fanned around the body diagonal `0->7`; their
//! signed volumes sum to the hexahedron volume. Exact for an affine
//! (box/sheared) cell; the conventional approximation for a trilinear
//! corner-point cell.

use crate::point::Point3;

/// Signed volume of the tetrahedron `(a, b, c, d)`.
#[inline]
fn tet_volume(a: Point3, b: Point3, c: Point3, d: Point3) -> f64 {
    (a - d).dot((b - d).cross(c - d)) / 6.0
}

/// Volume of a hexahedral cell from its 8 corners (see module docs for order).
#[must_use]
#[inline]
pub fn hexahedron_volume(c: &[Point3; 8]) -> f64 {
    // Six tets sharing the 0->7 diagonal, walking the three faces not on it.
    const TETS: [(usize, usize); 6] = [(1, 3), (3, 2), (2, 6), (6, 4), (4, 5), (5, 1)];
    let v: f64 = TETS
        .iter()
        .map(|&(a, b)| tet_volume(c[0], c[a], c[b], c[7]))
        .sum();
    v.abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build the 8 corners of an axis-aligned box of size `(lx,ly,lz)` at origin.
    fn box_corners(lx: f64, ly: f64, lz: f64) -> [Point3; 8] {
        let mut c = [Point3::new(0.0, 0.0, 0.0); 8];
        for (idx, slot) in c.iter_mut().enumerate() {
            let di = (idx & 1) as f64;
            let dj = ((idx >> 1) & 1) as f64;
            let dk = ((idx >> 2) & 1) as f64;
            *slot = Point3::new(di * lx, dj * ly, dk * lz);
        }
        c
    }

    #[test]
    fn unit_cube_is_one() {
        assert!((hexahedron_volume(&box_corners(1.0, 1.0, 1.0)) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn box_volume_is_lx_ly_lz() {
        let v = hexahedron_volume(&box_corners(10.0, 20.0, 3.0));
        assert!((v - 600.0).abs() < 1e-9, "got {v}");
    }

    #[test]
    fn sheared_cell_preserves_volume() {
        // Shear the top face in x by 5 ft: a parallelepiped has base-area * height.
        let mut c = box_corners(10.0, 20.0, 3.0);
        for corner in &mut c[4..8] {
            corner.x += 5.0;
        }
        let v = hexahedron_volume(&c);
        assert!(
            (v - 600.0).abs() < 1e-9,
            "shear should preserve volume, got {v}"
        );
    }
}
