//! Scatter-conditioning micro-bench (`task_suite_scatter_perf`, `dev-docs/bench`):
//! the two halves of the adopted petekTools direct `MinCurvatureOperator` that
//! `grid_scatter`'s `ScatterConditioner` splits —
//!
//! * **factor** — assemble + band-LU-factor the bilinear conditioning operator for a
//!   horizon's fixed off-node sample geometry. The per-horizon conditioning cost
//!   (the direct solve that replaced the cap-bound ~60 s SOR); done once per surface.
//! * **resolve** — back-substitute a fresh depth vector through the reused factor.
//!   The MC lever: the sample (x,y) are fixed across draws, so re-seating with new
//!   depths is this cheap solve, not a re-factor.
//!
//! Scale is representative dense off-node scatter kept under the owner 60 s cap
//! (`sample_size(10)`): ~100×100 nodes, ~20k off-node samples. The operator is the
//! petekTools primitive `ScatterConditioner` wraps, benched directly here (the
//! crate-private conditioner is not reachable from a bench crate).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use petektools::{Conditioning, Lattice, MinCurvatureOperator};

const NODES: usize = 100; // 100×100 node solve lattice
const PTS: usize = 140; // 140×140 off-node samples → 19_600 controls

/// A dense off-node scatter geometry on the unit-spaced solve lattice: fractional
/// node positions (kept off every node so the Bilinear data-fit governs) plus a
/// dome-shaped depth field.
fn scatter() -> (Vec<[f64; 2]>, Vec<f64>) {
    let mut xy = Vec::with_capacity(PTS * PTS);
    let mut z = Vec::with_capacity(PTS * PTS);
    let span = (NODES - 1) as f64;
    for j in 0..PTS {
        for i in 0..PTS {
            // Spread across the interior, offset off-node by (0.37, 0.61).
            let fi = 0.5 + (span - 1.0) * (i as f64 / (PTS - 1) as f64) + 0.37;
            let fj = 0.5 + (span - 1.0) * (j as f64 / (PTS - 1) as f64) + 0.61;
            xy.push([fi, fj]);
            z.push(2000.0 + 0.02 * (fi - span / 2.0).powi(2) + 0.015 * (fj - span / 2.0).powi(2));
        }
    }
    (xy, z)
}

fn bench(c: &mut Criterion) {
    let lattice = Lattice::regular(0.0, 0.0, 1.0, 1.0, NODES, NODES);
    let (xy, z) = scatter();

    let mut g = c.benchmark_group("scatter_condition");
    g.sample_size(10);

    // factor: per-horizon conditioning (assemble + band-LU factor), done once.
    g.bench_function("factor_100x100_19600pts", |b| {
        b.iter(|| {
            MinCurvatureOperator::factor(
                black_box(&lattice),
                black_box(&xy),
                Conditioning::Bilinear,
            )
            .unwrap()
        })
    });

    // resolve: the MC lever — back-substitute a fresh depth vector through the
    // already-factored operator (fixed geometry, varying depths).
    let op = MinCurvatureOperator::factor(&lattice, &xy, Conditioning::Bilinear).unwrap();
    g.bench_function("resolve_100x100_19600pts", |b| {
        b.iter(|| op.solve(black_box(&z)).unwrap())
    });

    g.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
