//! Adoption-readiness probe (`decision_kernel_home`): cross-check our `solve_surface`
//! against `petektools::gridding::grid(MinimumCurvature)` before we adopt the toolkit kernel.
//!
//! The production swap is deferred until we adopt petekTools' warm-start
//! `ConvergentGridder` (the actual perf win for the refine loop). These
//! tests de-risk that swap: the critical property is that the kernel **reproduces a
//! plane** (honours regional dip) the way our accuracy reference (`reproduces_a_plane`)
//! requires — a kernel that flattened edges would regress us. Aggregates are
//! orientation-independent (no assumption about the returned Array2's axis order).

use petektools::{gridding::grid, GridMethod, Lattice};
use srs_gridder::{solve_surface, Control, SolveOpts};

const NX: usize = 10;
const NY: usize = 10;

fn plane(ip: usize, jp: usize) -> f64 {
    5000.0 + 2.0 * ip as f64 + 3.0 * jp as f64
}

/// Sparse control points (corners + edges + centre) lying on the plane.
fn control_nodes() -> Vec<(usize, usize)> {
    vec![
        (0, 0),
        (9, 0),
        (0, 9),
        (9, 9), // corners
        (5, 5), // centre
        (0, 5),
        (9, 5),
        (5, 0),
        (5, 9), // edge midpoints
    ]
}

fn stats(vals: impl Iterator<Item = f64>) -> (f64, f64, f64) {
    let defined: Vec<f64> = vals.filter(|v| !v.is_nan()).collect();
    let min = defined.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = defined.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean = defined.iter().sum::<f64>() / defined.len() as f64;
    (min, max, mean)
}

#[test]
fn our_solver_reproduces_the_plane() {
    let controls: Vec<Control> = control_nodes()
        .iter()
        .map(|&(ip, jp)| Control {
            ip,
            jp,
            z: plane(ip, jp),
        })
        .collect();
    let s = solve_surface(
        NX,
        NY,
        &controls,
        SolveOpts {
            tol: 1e-9,
            max_iter: 60_000,
            ..Default::default()
        },
    )
    .unwrap();
    let vals = (0..NY)
        .flat_map(|jp| (0..NX).map(move |ip| (ip, jp)))
        .map(|(ip, jp)| s.z(ip, jp));
    let (min, max, mean) = stats(vals);
    // plane: 5000 .. 5045, mean 5022.5
    assert!((min - 5000.0).abs() < 0.1, "min {min}");
    assert!((max - 5045.0).abs() < 0.1, "max {max}");
    assert!((mean - 5022.5).abs() < 0.1, "mean {mean}");
}

#[test]
fn petektools_min_curvature_reproduces_the_plane() {
    let coords: Vec<[f64; 3]> = control_nodes()
        .iter()
        .map(|&(ip, jp)| [ip as f64, jp as f64, plane(ip, jp)])
        .collect();
    let lattice = Lattice::regular(0.0, 0.0, 1.0, 1.0, NX, NY);
    let field = grid(&coords, &lattice, GridMethod::MinimumCurvature).expect("grid()");

    let (min, max, mean) = stats(field.iter().copied());
    // If the external kernel honours regional dip (as our reference requires), a plane is
    // reproduced: min~5000, max~5045, mean~5022.5. A flattening kernel would pull these in.
    assert!(
        (min - 5000.0).abs() < 1.0,
        "petektools plane min {min} (edge flattening?)"
    );
    assert!(
        (max - 5045.0).abs() < 1.0,
        "petektools plane max {max} (edge flattening?)"
    );
    assert!((mean - 5022.5).abs() < 1.0, "petektools plane mean {mean}");

    // PER-NODE gate (`decision_gridder_kernel_unification`: aggregates hid the old
    // 12.48-ft interior sag while min/max/mean passed). After petekTools f81b6a6
    // (natural-dip boundary) the interior no longer sags — every node reproduces the
    // plane to well under a foot. `grid` returns shape (ncol, nrow) = (NX, NY),
    // indexed [[ip, jp]].
    let mut worst = 0.0_f64;
    for jp in 0..NY {
        for ip in 0..NX {
            worst = worst.max((field[[ip, jp]] - plane(ip, jp)).abs());
        }
    }
    assert!(
        worst < 0.1,
        "petektools per-node plane drift {worst} (interior sag not eliminated?)"
    );
}
