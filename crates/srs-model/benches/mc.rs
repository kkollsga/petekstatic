//! Monte-Carlo hot-path benches (`dev-docs/bench`): the per-realization cost the
//! MC driver pays. Two comparisons:
//!
//! * **V2 log-population presort** — the old per-cell linear scan over all log
//!   samples vs the new TVD-sorted binary-search window. Isolated micro-bench
//!   over representative cell depth ranges + a fine log so the delta is the
//!   search cost, not the surrounding build.
//! * **V7 realize / in-place** — a realistic `realize` + full `in_place` (which
//!   materializes a per-cell HCPV cube) vs `realize` + the summary-only
//!   `in_place_summary` (aggregates only, no per-cell Vec).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use srs_gridder::{Conformity, SolveOpts};
use srs_model::{BuildOpts, ConstantPriors, RealizationDraw, StaticModelTemplate};
use srs_wireframe::{
    Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
};

// --- V2: log population, linear scan vs sorted binary search ---

/// Cell depth ranges for an `nk`-layer column over `[top, top+gross]`, repeated
/// over `ncol*ncol` areal columns (every column sees the same layer-cake — the
/// current single-well placement).
fn cell_ranges(ncol: usize, nk: usize, top: f64, gross: f64) -> Vec<(f64, f64)> {
    let dz = gross / nk as f64;
    let mut v = Vec::with_capacity(ncol * ncol * nk);
    for _ in 0..(ncol * ncol) {
        for k in 0..nk {
            let lo = top + k as f64 * dz;
            v.push((lo, lo + dz));
        }
    }
    v
}

/// A fine, TVD-sorted synthetic log over `[top, top+gross]`.
fn samples(n: usize, top: f64, gross: f64) -> Vec<(f64, f64, f64)> {
    (0..n)
        .map(|i| {
            let tvd = top + gross * (i as f64 / n as f64);
            (tvd, 0.25, 0.3)
        })
        .collect()
}

fn population(c: &mut Criterion) {
    let ranges = cell_ranges(40, 25, 5000.0, 300.0); // 40x40x25 = 40k cells
    let s = samples(1500, 5000.0, 300.0);
    let mut g = c.benchmark_group("v2_log_population");

    // Old path: per cell, filter the whole sample array.
    g.bench_function("linear_scan", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for &(lo, hi) in &ranges {
                let mut sum = 0.0;
                let mut n = 0u32;
                for &(tvd, phi, _) in &s {
                    if tvd >= lo && tvd <= hi {
                        sum += phi;
                        n += 1;
                    }
                }
                acc += if n > 0 { sum / f64::from(n) } else { 0.25 };
            }
            black_box(acc)
        })
    });

    // New path: two binary searches bound a contiguous window.
    g.bench_function("binary_search", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for &(lo, hi) in &ranges {
                let start = s.partition_point(|(tvd, _, _)| *tvd < lo);
                let end = s.partition_point(|(tvd, _, _)| *tvd <= hi);
                let win = &s[start..end];
                acc += if win.is_empty() {
                    0.25
                } else {
                    win.iter().map(|(_, p, _)| *p).sum::<f64>() / win.len() as f64
                };
            }
            black_box(acc)
        })
    });
    g.finish();
}

// --- V7: realize + in-place, full vs summary ---

fn flat_wireframe(n: usize, depth_m: f64, owc_m: f64) -> Wireframe {
    Wireframe {
        boundary: Boundary {
            ring: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
            hardness: Hardness::Interpolated,
        },
        horizons: std::sync::Arc::new(vec![Horizon {
            name: "top".into(),
            role: HorizonRole::Top,
            surface: GriddedDepth {
                ncol: n,
                nrow: n,
                depth_m: vec![depth_m; n * n],
                is_control: vec![true; n * n],
            },
        }]),
        contacts: vec![Contact {
            kind: ContactKind::Owc,
            depth_m: owc_m,
            hardness: Hardness::Hard,
        }],
    }
}

fn opts(nk: usize) -> BuildOpts {
    BuildOpts {
        area_m2: 100.0,
        gross_height_m: 50.0,
        nk,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.25,
            net_to_gross: 0.8,
            water_saturation: 0.3,
        },
    }
}

fn realize(c: &mut Criterion) {
    let wf = flat_wireframe(51, 5000.0, 5025.0); // 50x50 lattice
    let o = opts(20);
    let draw = RealizationDraw::new(100.0, 50.0, 5025.0, 0.25, 0.8, 0.3, 0);
    let mut g = c.benchmark_group("v7_realize");

    g.bench_function("realize_then_full_in_place", |b| {
        let mut t = StaticModelTemplate::new(&wf, o).unwrap();
        b.iter(|| {
            let m = t.realize(black_box(&draw)).unwrap();
            black_box(m.in_place().unwrap().hcpv_m3)
        })
    });

    g.bench_function("realize_then_summary_in_place", |b| {
        let mut t = StaticModelTemplate::new(&wf, o).unwrap();
        b.iter(|| {
            let m = t.realize(black_box(&draw)).unwrap();
            black_box(m.in_place_summary().unwrap().hcpv_m3)
        })
    });

    g.finish();
}

// --- P5: property-pipeline realize cost, LevelShift vs Resimulate ---

fn realize_property_modes(c: &mut Criterion) {
    use petektools::{Variogram, VariogramModel};
    use srs_model::{Gaussian, McMode, PropertyPipeline, UpscaleMethod, WellLog};

    let wf = flat_wireframe(51, 5000.0, 5025.0); // 50x50 lattice, side 10 m
    let o = opts(20); // nk = 20, gross 50 -> dz 2.5
    let draw = RealizationDraw::new(100.0, 50.0, 5025.0, 0.25, 0.8, 0.3, 0);

    // Four positioned wells spanning all 20 layers.
    let col = |x: f64, y: f64, v: f64| {
        WellLog::new(
            x,
            y,
            (0..20)
                .map(|k| (5000.0 + 2.5 * k as f64 + 1.25, v))
                .collect(),
        )
    };
    let wells = vec![
        col(1.0, 1.0, 0.20),
        col(9.0, 9.0, 0.28),
        col(1.0, 9.0, 0.24),
        col(9.0, 1.0, 0.22),
    ];
    let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 3.0).unwrap();
    let pipe = PropertyPipeline::new("PHIE")
        .upscale(wells, UpscaleMethod::Arithmetic)
        .propagate(Gaussian::new(vgm, 1));

    let mut g = c.benchmark_group("mc_property_modes");

    // LevelShift: propagate once (warmed before the loop), then each realize reuses
    // the pattern with only a level shift — must stay ~ms-class.
    g.bench_function("realize_level_shift", |b| {
        let mut t = StaticModelTemplate::new(&wf, o)
            .unwrap()
            .with_property(pipe.clone());
        let _ = t.realize(&draw).unwrap(); // warm the cached pattern
        b.iter(|| {
            let m = t.realize(black_box(&draw)).unwrap();
            black_box(m.property("PHIE").unwrap().values.len())
        })
    });

    // Resimulate: a full per-layer SGS every realize — its true (heavier) cost.
    g.bench_function("realize_resimulate", |b| {
        let mut t = StaticModelTemplate::new(&wf, o)
            .unwrap()
            .with_property_mode(pipe.clone(), McMode::Resimulate);
        b.iter(|| {
            let m = t.realize(black_box(&draw)).unwrap();
            black_box(m.property("PHIE").unwrap().values.len())
        })
    });

    g.finish();
}

// --- Layering conformity: Proportional vs FollowTop at ~1 m dz ---

/// A wedge wireframe: flat top at 5000, a Base horizon dipping in i from ~5 m
/// (updip) to ~50 m (downdip) thick — so a fine-dz FollowTop truncates the thin
/// updip columns while fully layering the thick downdip ones.
fn wedge_wireframe(n: usize) -> Wireframe {
    let top = vec![5000.0; n * n];
    let mut base = vec![0.0; n * n];
    for r in 0..n {
        for c in 0..n {
            base[r * n + c] = 5000.0 + 5.0 + 45.0 * (c as f64 / (n as f64 - 1.0));
        }
    }
    let surf = |depth_m: Vec<f64>| GriddedDepth {
        ncol: n,
        nrow: n,
        depth_m,
        is_control: vec![true; n * n],
    };
    Wireframe {
        boundary: Boundary {
            ring: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
            hardness: Hardness::Hard,
        },
        horizons: std::sync::Arc::new(vec![
            Horizon {
                name: "top".into(),
                role: HorizonRole::Top,
                surface: surf(top),
            },
            Horizon {
                name: "base".into(),
                role: HorizonRole::Base,
                surface: surf(base),
            },
        ]),
        contacts: vec![Contact {
            kind: ContactKind::Owc,
            depth_m: 9000.0,
            hardness: Hardness::Hard,
        }],
    }
}

/// Build (template) + one realize + full in-place under Proportional vs FollowTop
/// at ~1 m dz — the fine-layering scale the polish work makes default. FollowTop's
/// derived nk (= ceil(50/1) = 50) matches the Proportional nk so the two do equal
/// layer work; the delta is the conformity's per-column truncation bookkeeping.
fn realize_conformity(c: &mut Criterion) {
    let wf = wedge_wireframe(51); // 50x50 lattice
    let draw = RealizationDraw::new(100.0, 50.0, 9000.0, 0.25, 0.8, 0.3, 0);
    let prop = {
        let mut o = opts(50); // nk = 50 to match FollowTop's dz-derived count
        o.gross_height_m = 50.0;
        o
    };
    let ftop = {
        let mut o = prop;
        o.conformity = Conformity::FollowTop { dz_m: 1.0 };
        o
    };

    let mut g = c.benchmark_group("conformity_realize_1m_dz");
    g.bench_function("proportional", |b| {
        let mut t = StaticModelTemplate::new(&wf, prop).unwrap();
        b.iter(|| {
            let m = t.realize(black_box(&draw)).unwrap();
            black_box(m.in_place().unwrap().hcpv_m3)
        })
    });
    g.bench_function("follow_top", |b| {
        let mut t = StaticModelTemplate::new(&wf, ftop).unwrap();
        b.iter(|| {
            let m = t.realize(black_box(&draw)).unwrap();
            black_box(m.in_place().unwrap().hcpv_m3)
        })
    });
    g.finish();
}

// --- structural uncertainty: per-draw perturbation-field cost ---
//
// The per-draw cost of the correlated structural perturbation fields
// (`decision_structural_uncertainty_isochore`): one unconditional SGS on the areal
// node lattice per uncertain horizon per draw. Benched on a stack template against
// the fixed-surface (no-perturbation) realize, so the delta IS the field cost.
fn realize_structural(c: &mut Criterion) {
    use petektools::{Variogram, VariogramModel};
    use srs_gridder::Conformity;
    use srs_model::{
        HorizonSource, HorizonStack, PerturbationField, StackHorizon, StackZone, ZoneDraw,
    };
    use srs_wireframe::GriddedDepth;

    // A 50×50 node lattice (2 500 nodes/field), 3-horizon / 2-zone flat stack.
    let n = 50usize;
    let surf = |d: f64| GriddedDepth {
        ncol: n,
        nrow: n,
        depth_m: vec![d; n * n],
        is_control: vec![true; n * n],
    };
    let mapped = |name: &str, d: f64| StackHorizon {
        name: name.into(),
        source: HorizonSource::Mapped(surf(d)),
    };
    let stack = || HorizonStack {
        horizons: vec![
            mapped("H0", 5000.0),
            mapped("H1", 5030.0),
            mapped("H2", 5060.0),
        ],
        zone_layers: vec![
            StackZone::new("Z0", Conformity::Proportional, 8, Vec::new()),
            StackZone::new("Z1", Conformity::Proportional, 8, Vec::new()),
        ],
    };
    let o = opts(0);
    let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 800.0).unwrap();
    let field = || PerturbationField::new(8.0, vgm);
    let plain = RealizationDraw::new(1_000_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 1);
    let structural = plain
        .clone()
        .with_top_structural(field())
        .with_zone_draw(ZoneDraw::new(0).with_isochore_structural(field()))
        .with_zone_draw(ZoneDraw::new(1).with_isochore_structural(field()));

    let mut g = c.benchmark_group("mc_structural");
    g.bench_function("realize_fixed_surfaces", |b| {
        let mut t = StaticModelTemplate::from_horizon_stack(stack(), o).unwrap();
        b.iter(|| {
            let m = t.realize(black_box(&plain)).unwrap();
            black_box(m.grid().bulk_volume())
        })
    });
    // Three perturbation fields per draw (top + two isochores) — the delta vs the
    // fixed-surface bench above is the structural field cost.
    g.bench_function("realize_three_fields", |b| {
        let mut t = StaticModelTemplate::from_horizon_stack(stack(), o).unwrap();
        b.iter(|| {
            let m = t.realize(black_box(&structural)).unwrap();
            black_box(m.grid().bulk_volume())
        })
    });
    g.finish();
}

criterion_group!(
    benches,
    population,
    realize,
    realize_property_modes,
    realize_conformity,
    realize_structural
);
criterion_main!(benches);
