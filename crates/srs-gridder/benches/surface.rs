use criterion::{criterion_group, criterion_main, Criterion};
use srs_gridder::{solve_surface, solve_surface_seeded, Control, KernelSurface, SolveOpts};

fn controls(n: usize) -> Vec<Control> {
    vec![
        Control {
            ip: 0,
            jp: 0,
            z: 5000.0,
        },
        Control {
            ip: n - 1,
            jp: 0,
            z: 5000.0,
        },
        Control {
            ip: 0,
            jp: n - 1,
            z: 5000.0,
        },
        Control {
            ip: n - 1,
            jp: n - 1,
            z: 5000.0,
        },
        Control {
            ip: n / 2,
            jp: n / 2,
            z: 4960.0,
        },
    ]
}

fn bench_surface(c: &mut Criterion) {
    let n = 51; // a 50x50-cell lattice
    let controls = controls(n);

    // Cold solve (the accuracy reference; `50x50` is the hot lifecycle flag).
    c.bench_function("solve_surface 50x50", |b| {
        b.iter(|| solve_surface(n, n, &controls, SolveOpts::default()).unwrap())
    });

    // Warm-start refine (petekTools ConvergentGridder): re-solve from a nearly-
    // converged seed — the per-realization regeneration cost. The seed is a prior
    // petekTools-kernel field, perturbed by nudging one control, so the warm
    // re-solve is a few sweeps rather than a cold relaxation.
    // Kernel-space bootstrap (the KernelSurface newtype forbids seeding from the
    // cold solver's output — warm==cold holds only within one kernel).
    let flat = solve_surface_seeded(&KernelSurface::flat(n, n, 5000.0), &controls).unwrap();
    let mut nudged = controls.clone();
    nudged[4].z = 4955.0; // one control moves between realizations
    c.bench_function("solve_surface_seeded 50x50 (warm refine)", |b| {
        b.iter(|| solve_surface_seeded(&flat, &nudged).unwrap())
    });
}

criterion_group!(benches, bench_surface);
criterion_main!(benches);
