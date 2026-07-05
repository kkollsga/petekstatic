//! Canonical-scale scatter-build **profiler** (`task_suite_scatter_perf`).
//! `#[ignore]` — release-only, wall-clock, prints a stage breakdown; not a CI
//! assertion. Run with the profile eprintlns on:
//!
//! ```text
//! SRS_PROFILE=1 cargo test -p srs-model --release --test scatter_profile -- --ignored --nocapture --test-threads=1
//! ```
//!
//! Shape mirrors the canonical real model: an ~122×116 node lattice, 11 horizons
//! resolved from dense off-node scatter (~39k points total), 10 proportional
//! zones. All-synthetic (a fictional dipping dome) — no dataset content.

use std::time::Instant;

use petekstatic::gridder::{solve_surface_converged, Conformity, Control, SolveOpts};
use petekstatic::model::{
    BuildOpts, ConstantPriors, Georef, HorizonSource, HorizonStack, StackFrame, StackHorizon,
    StackZone, StaticModelBuilder, StaticModelTemplate, WorldPoint,
};

const NX: usize = 122; // nodes along i
const NY: usize = 116; // nodes along j
const SPACING: f64 = 100.0;
const NH: usize = 11; // horizons
const PTS_PER_AXIS: usize = 60; // 60×60 = 3600 pts/horizon → ~39.6k total
const ORIGIN_X: f64 = 431_000.0;
const ORIGIN_Y: f64 = 6_521_000.0;

/// Smooth truth depth (positive-down) for horizon `h`: a regional dip + a broad
/// dome, offset down per horizon so the stack never crosses.
fn truth(h: usize, x: f64, y: f64) -> f64 {
    let regional = 2000.0 + 0.02 * x + 0.015 * y;
    let (cx, cy) = (NX as f64 * SPACING * 0.55, NY as f64 * SPACING * 0.45);
    let r2 = (x - cx).powi(2) + (y - cy).powi(2);
    let dome = 80.0 * (-r2 / (2.0 * 3000.0_f64.powi(2))).exp();
    regional - dome + h as f64 * 30.0
}

/// Dense off-node scatter for horizon `h` over the frame extent (points sit at
/// non-integer node fractions by construction).
fn scatter(h: usize) -> Vec<WorldPoint> {
    let ext_x = (NX - 1) as f64 * SPACING;
    let ext_y = (NY - 1) as f64 * SPACING;
    let step_x = ext_x / PTS_PER_AXIS as f64;
    let step_y = ext_y / PTS_PER_AXIS as f64;
    let mut pts = Vec::with_capacity(PTS_PER_AXIS * PTS_PER_AXIS);
    for j in 0..PTS_PER_AXIS {
        for i in 0..PTS_PER_AXIS {
            let x = 13.0 + i as f64 * step_x; // 13 m offset → off-node
            let y = 17.0 + j as f64 * step_y;
            pts.push(WorldPoint {
                x: ORIGIN_X + x,
                y: ORIGIN_Y + y,
                depth_m: truth(h, x, y),
            });
        }
    }
    pts
}

fn frame() -> StackFrame {
    StackFrame {
        ni: NX - 1,
        nj: NY - 1,
        georef: Georef::new(
            ORIGIN_X + 0.5 * SPACING,
            ORIGIN_Y + 0.5 * SPACING,
            SPACING,
            SPACING,
        )
        .unwrap(),
    }
}

fn scatter_stack() -> HorizonStack {
    let horizons: Vec<StackHorizon> = (0..NH)
        .map(|h| StackHorizon {
            name: format!("H{h}"),
            source: HorizonSource::Scatter(scatter(h)),
        })
        .collect();
    let zone_layers: Vec<StackZone> = (0..NH - 1)
        .map(|z| StackZone::new(format!("Z{z}"), Conformity::Proportional, 4, vec![]))
        .collect();
    HorizonStack {
        horizons,
        zone_layers,
    }
}

fn opts() -> BuildOpts {
    let ext = (NX - 1) as f64 * SPACING;
    BuildOpts {
        area_m2: ext * ext,
        gross_height_m: 30.0,
        nk: 4,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.2,
            net_to_gross: 1.0,
            water_saturation: 0.25,
        },
    }
}

/// Seed-sensitivity of the cold Bilinear conditioning solve (the ~4.7s/horizon
/// hotspot). If a near-answer seed collapses the time, the cold solve is
/// convergence/cap-bound → cross-horizon warm-start is a lever. If not, it is
/// per-sweep-cost-bound → parallelism is the only in-repo lever.
#[test]
#[ignore = "profiler: SRS_PROFILE=1 cargo test --release --test scatter_profile -- --ignored --nocapture"]
fn profile_conditioning_seed_sensitivity() {
    use petektools::{grid_min_curvature_conditioned, Conditioning, Lattice};
    let g = frame().georef;
    let lattice = Lattice::regular(0.0, 0.0, 1.0, 1.0, NX, NY);
    let pts = scatter(0);
    let coords: Vec<[f64; 3]> = pts
        .iter()
        .map(|p| {
            let fi = (p.x - g.origin_x) / g.spacing_x + 0.5;
            let fj = (p.y - g.origin_y) / g.spacing_y + 0.5;
            [fi, fj, p.depth_m]
        })
        .collect();
    // Cold (None → IDW seed): the current path.
    let t0 = Instant::now();
    let cold =
        grid_min_curvature_conditioned(&coords, &lattice, None, Conditioning::Bilinear).unwrap();
    eprintln!(
        "[profile] conditioning COLD (None seed) = {:.1} ms",
        t0.elapsed().as_secs_f64() * 1e3
    );
    // Warm: seed with the near-answer field (the previous cold result perturbed
    // like a neighbouring stacked horizon would be).
    let mut seed = cold.clone();
    seed.iter_mut().for_each(|z| *z += 30.0);
    let t1 = Instant::now();
    let _warm =
        grid_min_curvature_conditioned(&coords, &lattice, Some(&seed), Conditioning::Bilinear)
            .unwrap();
    eprintln!(
        "[profile] conditioning WARM (near-answer seed) = {:.1} ms",
        t1.elapsed().as_secs_f64() * 1e3
    );
}

/// ONE 122×116 converged solve with a realistic sparse-ish control set (the
/// per-horizon inner cost). The smallest thing that reproduces the solve cost.
#[test]
#[ignore = "profiler: SRS_PROFILE=1 cargo test --release --test scatter_profile -- --ignored --nocapture"]
fn profile_single_solve() {
    // ~70% of nodes controlled (mirrors a dense scatter's support hull).
    let mut controls = Vec::new();
    for jp in 0..NY {
        for ip in 0..NX {
            if (ip * 7 + jp * 13) % 10 < 7 {
                let (x, y) = (ip as f64 * SPACING, jp as f64 * SPACING);
                controls.push(Control {
                    ip,
                    jp,
                    z: truth(0, x, y),
                });
            }
        }
    }
    eprintln!(
        "[profile] single solve: {}x{} lattice, {} controls",
        NX,
        NY,
        controls.len()
    );
    let t0 = Instant::now();
    let _s = solve_surface_converged(NX, NY, &controls).unwrap();
    eprintln!(
        "[profile] single solve total: {:.1} ms",
        t0.elapsed().as_secs_f64() * 1e3
    );
}

/// The full 11-horizon scatter build, staged: builder path only (no MC template,
/// no parity) — the base cost one build pays.
#[test]
#[ignore = "profiler: SRS_PROFILE=1 cargo test --release --test scatter_profile -- --ignored --nocapture"]
fn profile_full_scatter_stack() {
    let t_all = Instant::now();
    let t0 = Instant::now();
    let builder = StaticModelBuilder::from_scatter_stack(scatter_stack(), opts(), frame()).unwrap();
    let t_cond_resolve = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let _m = builder.build().unwrap();
    let t_build = t1.elapsed().as_secs_f64();
    eprintln!(
        "[profile] from_scatter_stack (condition+resolve) = {:.2} s | build() = {:.2} s | TOTAL = {:.2} s",
        t_cond_resolve,
        t_build,
        t_all.elapsed().as_secs_f64()
    );
}

/// The template path (MC template re-solve) — the second of the 3× redundant
/// full solves the task note calls out.
#[test]
#[ignore = "profiler: SRS_PROFILE=1 cargo test --release --test scatter_profile -- --ignored --nocapture"]
fn profile_template_scatter_stack() {
    let t0 = Instant::now();
    let _t = StaticModelTemplate::from_scatter_stack(scatter_stack(), opts(), frame()).unwrap();
    eprintln!(
        "[profile] StaticModelTemplate::from_scatter_stack = {:.2} s",
        t0.elapsed().as_secs_f64()
    );
}
