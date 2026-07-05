//! Box-grid construction — the degenerate corner-point grid (`AlgorithmSpec`
//! `box_grid_construction_spec`, Ponting 1989). Vertical, evenly-spaced,
//! untilted, unfaulted pillars: the Day-1 model that upgrades in place.
//!
//! Footprint from area `A` (m²): with aspect `r = Lx/Ly` (default square),
//! `Lx = sqrt(A*r)`, `Ly = sqrt(A/r)`, `Lz = h`. Cell sizes `dx=Lx/ni`,
//! `dy=Ly/nj`, `dz=h/nk`; cell origin `x=i*dx, y=j*dy, z=z_top + k*dz`. All
//! lengths in metres (family SI standard); depth is positive-down.

use crate::error::StaticError;
use crate::grid::geometry::{CornerPointGeom, Pillar};
use crate::grid::grid::Grid;
use crate::grid::index::Dims;
use crate::grid::point::Point3;

/// Inputs to box-grid construction (the spec's `inputs_contract`).
#[derive(Debug, Clone, Copy)]
pub struct BoxSpec {
    /// Reservoir area \[m²\].
    pub area_m2: f64,
    /// Gross height \[m\].
    pub gross_height_m: f64,
    /// Grid resolution.
    pub dims: Dims,
    /// Top depth \[m\] (positive-down; z increases downward).
    pub top_depth_m: f64,
    /// Aspect ratio `Lx/Ly` (1.0 = square).
    pub aspect_ratio: f64,
}

impl BoxSpec {
    /// A square box at the surface with the given area, height and resolution.
    #[must_use]
    pub fn square(area_m2: f64, gross_height_m: f64, dims: Dims) -> Self {
        Self {
            area_m2,
            gross_height_m,
            dims,
            top_depth_m: 0.0,
            aspect_ratio: 1.0,
        }
    }
}

/// Build the box grid from a [`BoxSpec`].
///
/// # Errors
/// Returns [`StaticError::InvalidInput`] if area, height, or aspect ratio is not
/// strictly positive.
pub fn build_box(spec: BoxSpec) -> Result<Grid, StaticError> {
    // Strictly positive and finite (rejects 0, negatives, NaN, inf).
    fn positive(x: f64) -> bool {
        x.is_finite() && x > 0.0
    }
    let require = |ok: bool, what: &str, value: f64| -> Result<(), StaticError> {
        if ok {
            Ok(())
        } else {
            Err(StaticError::InvalidInput(format!(
                "{what} must be a finite value > 0, got {value}"
            )))
        }
    };
    require(positive(spec.area_m2), "area", spec.area_m2)?;
    require(
        positive(spec.gross_height_m),
        "gross height",
        spec.gross_height_m,
    )?;
    require(
        positive(spec.aspect_ratio),
        "aspect ratio",
        spec.aspect_ratio,
    )?;

    let dims = spec.dims;
    let area_m2 = spec.area_m2;
    let lx = (area_m2 * spec.aspect_ratio).sqrt();
    let ly = (area_m2 / spec.aspect_ratio).sqrt();
    let dx = lx / dims.ni as f64;
    let dy = ly / dims.nj as f64;
    let dz = spec.gross_height_m / dims.nk as f64;
    let z_top = spec.top_depth_m;
    let z_bot = z_top + spec.gross_height_m;

    // COORD: vertical pillars on a uniform (ni+1)x(nj+1) lattice.
    let mut coord = Vec::with_capacity(dims.pillar_count());
    for jp in 0..=dims.nj {
        for ip in 0..=dims.ni {
            let x = ip as f64 * dx;
            let y = jp as f64 * dy;
            coord.push(Pillar {
                top: Point3::new(x, y, z_top),
                bottom: Point3::new(x, y, z_bot),
            });
        }
    }

    // ZCORN: 8 depths per cell; top four = layer top, bottom four = layer base.
    let mut zcorn = Vec::with_capacity(dims.cell_count() * 8);
    for c in dims.iter() {
        let zt = z_top + c.k as f64 * dz;
        let zb = zt + dz;
        // corners 0..4 are the top face (dk=0), 4..8 the bottom (dk=1).
        zcorn.extend_from_slice(&[zt, zt, zt, zt, zb, zb, zb, zb]);
    }

    Ok(Grid::new(CornerPointGeom::new(dims, coord, zcorn)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::index::Ijk;

    #[test]
    fn bulk_volume_equals_area_times_height() {
        // 400_000 m² * 50 m = 20_000_000 m³ (bulk = area × height, SI-native).
        let spec = BoxSpec::square(400_000.0, 50.0, Dims::new(10, 8, 5).unwrap());
        let grid = build_box(spec).unwrap();
        let expected = 400_000.0 * 50.0;
        let v = grid.bulk_volume();
        assert!(
            (v - expected).abs() / expected < 1e-10,
            "bulk volume {v} != {expected}"
        );
    }

    #[test]
    fn cell_count_matches_dims() {
        let grid = build_box(BoxSpec::square(40.0, 30.0, Dims::new(4, 5, 6).unwrap())).unwrap();
        assert_eq!(grid.cell_count(), 120);
    }

    #[test]
    fn cells_are_uniform_cuboids() {
        let dims = Dims::new(3, 3, 3).unwrap();
        let grid = build_box(BoxSpec::square(90.0, 60.0, dims)).unwrap();
        let v0 = grid.cell(Ijk::new(0, 0, 0)).volume();
        for c in grid.cells() {
            assert!(
                (c.volume() - v0).abs() / v0 < 1e-9,
                "non-uniform cell {c:?}"
            );
        }
    }

    #[test]
    fn depth_increases_with_k() {
        let dims = Dims::new(2, 2, 4).unwrap();
        let grid = build_box(BoxSpec {
            area_m2: 400_000.0,
            gross_height_m: 80.0,
            dims,
            top_depth_m: 5000.0,
            aspect_ratio: 1.0,
        })
        .unwrap();
        let top = grid.cell(Ijk::new(0, 0, 0)).top_depth();
        let deep = grid.cell(Ijk::new(0, 0, 3)).top_depth();
        assert!((top - 5000.0).abs() < 1e-9);
        assert!(deep > top);
        assert!((grid.cell(Ijk::new(0, 0, 0)).dz() - 20.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_nonpositive_inputs() {
        let d = Dims::new(2, 2, 2).unwrap();
        assert!(build_box(BoxSpec::square(0.0, 10.0, d)).is_err());
        assert!(build_box(BoxSpec::square(10.0, 0.0, d)).is_err());
    }
}
