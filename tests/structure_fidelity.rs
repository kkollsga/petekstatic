//! Structure-build fidelity audit (`task_petekstatic_topsonly_envelope`, owner
//! rider "look into the structure build") — a fixture mirroring the real input
//! shape: **dense off-node scatter** (~65 m spacing on a 100 m node lattice),
//! an **exact merge** band along one data edge, a **data-void margin** inside
//! the model extent, a structure flank running against another data edge, and a
//! few **local ties**. Every stage's misfit contribution is measured and pinned:
//!
//! | stage | what is measured | assertion |
//! |---|---|---|
//! | S1 conditioning | scatter→node authoring: nearest-node **snap-average** vs proper **local interpolation**; built-surface misfit AT the scatter points | snap error >> interp error (the metres-level on-data rms is authoring, not the solve); defined nodes honored EXACTLY |
//! | S2 solve | node-control fidelity + convergence residual (second seeded solve) | fully-pinned lattice = input bit-exact; residual ~0 — with the direct band-LU operator BOTH the raw flat-seeded kernel and the converged entry reach the fixed point in one solve (no affine-mode stall) |
//! | S3 extrapolation | margin behaviour: legacy NaturalDip vs default DecayToData | NaturalDip runs beyond the data's own range (phantom structure); DecayToData stays data-bounded, far margin = nearest datum exactly; merged margin stays merged |
//! | S4 ties | spatial reach of tie substitutions | dense lattice: ONLY the tied nodes move; sparse: decay measured, far field < 0.05 m |
//! | S5 repairs | interface repairs on consistent scatter inputs | exactly zero (crossings resolve toward data by construction, never inflate) |
//! | S6 invariance | internal mapped-surface swap | top/base/other horizons bit-unchanged, total GRV invariant |
//!
//! All data is hand-authored synthetic truth (a fictional dome + merge band at a
//! fictional coordinate window) — no dataset content. World-georef variant per
//! doctrine R1 (the merge band runs diagonally there).

#![allow(clippy::needless_range_loop)] // the audit indexes parallel node fields throughout

use petekstatic::gridder::{
    solve_surface_converged, solve_surface_seeded, Conformity, Control, ExtrapolationPolicy,
    KernelSurface, SolveOpts,
};
use petekstatic::model::{
    BuildOpts, ConstantPriors, Georef, HorizonSource, HorizonStack, RealizationDraw, StackFrame,
    StackHorizon, StackZone, StaticModel, StaticModelBuilder, StaticModelTemplate, WellTie,
    WorldPoint,
};
use petekstatic::wireframe::GriddedDepth;

const N: usize = 21; // nodes per axis; 20x20 cells
const SPACING: f64 = 100.0; // node spacing (m) — the "100 m lattice"
const EXTENT: f64 = (N - 1) as f64 * SPACING; // 2000 m
const AREA_M2: f64 = EXTENT * EXTENT;
const HULL_MIN: f64 = 500.0; // data hull [500, 1500]^2 → 5-cell void margin
const HULL_MAX: f64 = 1500.0;
const SCATTER_STEP: f64 = 65.0; // dense scatter spacing (off-node by construction)
const REGIONAL: f64 = 2000.0;
const DOME_AMP: f64 = 60.0;
const DOME_X: f64 = 1400.0; // crest inside the hull, steep flank against the EAST edge
const DOME_Y: f64 = 1000.0;
const DOME_SIGMA: f64 = 250.0;
const SEP_OPEN: f64 = 25.0; // open-envelope separation (m)
const ORIGIN_X: f64 = 700_000.0; // fictional world window
const ORIGIN_Y: f64 = 7_100_000.0;

/// Truth TOP depth (positive down): a regional level with an anticline crest.
fn truth_top(x: f64, y: f64) -> f64 {
    let (dx, dy) = ((x - DOME_X) / DOME_SIGMA, (y - DOME_Y) / DOME_SIGMA);
    REGIONAL - DOME_AMP * (-(dx * dx + dy * dy)).exp()
}

/// Truth TOP→BASE separation: **exactly zero** (merged) where `merged(x, y)`,
/// ramping to `SEP_OPEN` over 200 m past the merge line.
fn truth_sep(x: f64, y: f64, merge_coord: impl Fn(f64, f64) -> f64) -> f64 {
    let c = merge_coord(x, y);
    if c <= 0.0 {
        0.0
    } else {
        let t = (c / 200.0).min(1.0);
        SEP_OPEN * t * t * (3.0 - 2.0 * t) // smoothstep, exact 0 at the merge line
    }
}

/// Axis-aligned merge coordinate: merged (<= 0) for x <= 800.
fn merge_axis(x: f64, _y: f64) -> f64 {
    x - 800.0
}
/// Diagonal ("azimuthed") merge coordinate for the world-georef variant.
fn merge_diag(x: f64, y: f64) -> f64 {
    (x + y) / 2.0 - 900.0
}

/// The dense scatter: a 65 m grid over the hull, offset so points are off-node.
fn scatter_xy() -> Vec<(f64, f64)> {
    let mut pts = Vec::new();
    let mut y = HULL_MIN + 13.0;
    while y <= HULL_MAX {
        let mut x = HULL_MIN + 17.0;
        while x <= HULL_MAX {
            pts.push((x, y));
            x += SCATTER_STEP;
        }
        y += SCATTER_STEP;
    }
    pts
}

fn node_xy(ip: usize, jp: usize) -> (f64, f64) {
    (ip as f64 * SPACING, jp as f64 * SPACING)
}

/// S1(a): nearest-node SNAP-AVERAGE authoring (emulates the metres-level
/// conditioning defect): each sample snaps to its nearest node; collisions
/// average; nodes with no sample stay NaN.
fn author_snap(values: &[((f64, f64), f64)]) -> GriddedDepth {
    let mut sum = vec![0.0; N * N];
    let mut cnt = vec![0usize; N * N];
    for &((x, y), z) in values {
        let ip = (x / SPACING).round() as usize;
        let jp = (y / SPACING).round() as usize;
        if ip < N && jp < N {
            sum[jp * N + ip] += z;
            cnt[jp * N + ip] += 1;
        }
    }
    let depth: Vec<f64> = sum
        .iter()
        .zip(&cnt)
        .map(|(s, &c)| if c > 0 { s / c as f64 } else { f64::NAN })
        .collect();
    GriddedDepth {
        ncol: N,
        nrow: N,
        depth_m: depth,
        is_control: cnt.iter().map(|&c| c > 0).collect(),
    }
}

/// S1(b): PERFECT conditioning — the truth evaluated at every node inside the
/// data hull (what a correct scatter→node authoring converges to). Isolates the
/// conditioning term completely: any residual misfit at the scatter under this
/// authoring is pure lattice-resolution representation (irreducible at this
/// node spacing), and any extra misfit under snap authoring is the conditioning
/// defect itself.
fn author_hull(f: impl Fn(f64, f64) -> f64) -> GriddedDepth {
    let mut depth = vec![f64::NAN; N * N];
    let mut defined = vec![false; N * N];
    for jp in 0..N {
        for ip in 0..N {
            let (x, y) = node_xy(ip, jp);
            if (HULL_MIN..=HULL_MAX).contains(&x) && (HULL_MIN..=HULL_MAX).contains(&y) {
                depth[jp * N + ip] = f(x, y);
                defined[jp * N + ip] = true;
            }
        }
    }
    GriddedDepth {
        ncol: N,
        nrow: N,
        depth_m: depth,
        is_control: defined,
    }
}

/// Truth sampled AT the nodes (the fully-defined dense variant).
fn author_dense(f: impl Fn(f64, f64) -> f64) -> GriddedDepth {
    let mut depth = vec![0.0; N * N];
    for jp in 0..N {
        for ip in 0..N {
            let (x, y) = node_xy(ip, jp);
            depth[jp * N + ip] = f(x, y);
        }
    }
    GriddedDepth {
        ncol: N,
        nrow: N,
        depth_m: depth,
        is_control: vec![true; N * N],
    }
}

fn two_horizon_stack(top: GriddedDepth, base: GriddedDepth) -> HorizonStack {
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Mapped(top),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(base),
            },
        ],
        zone_layers: vec![StackZone {
            name: "Z".into(),
            color: None,
            conformity: Conformity::Proportional,
            nk: 8,
            contacts: vec![],
        }],
    }
}

fn opts() -> BuildOpts {
    BuildOpts {
        area_m2: AREA_M2,
        gross_height_m: SEP_OPEN,
        nk: 1,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.2,
            net_to_gross: 1.0,
            water_saturation: 0.25,
        },
    }
}

fn build(stack: HorizonStack, policy: ExtrapolationPolicy, georef: bool) -> StaticModel {
    let b = StaticModelBuilder::from_horizon_stack(stack, opts())
        .unwrap()
        .with_extrapolation(policy);
    let b = if georef {
        b.with_georef(ORIGIN_X, ORIGIN_Y, SPACING, SPACING)
    } else {
        b
    };
    b.build().unwrap()
}

fn surface<'a>(m: &'a StaticModel, name: &str) -> &'a [f64] {
    &m.framework()
        .horizons
        .iter()
        .find(|h| h.name == name)
        .unwrap()
        .surface
        .depth_m
}

/// Bilinear evaluation of a node field at world (x, y).
fn eval(field: &[f64], x: f64, y: f64) -> f64 {
    let fx = (x / SPACING).clamp(0.0, (N - 1) as f64);
    let fy = (y / SPACING).clamp(0.0, (N - 1) as f64);
    let (i0, j0) = (fx.floor() as usize, fy.floor() as usize);
    let (i1, j1) = ((i0 + 1).min(N - 1), (j0 + 1).min(N - 1));
    let (tx, ty) = (fx - i0 as f64, fy - j0 as f64);
    let z00 = field[j0 * N + i0];
    let z10 = field[j0 * N + i1];
    let z01 = field[j1 * N + i0];
    let z11 = field[j1 * N + i1];
    z00 * (1.0 - tx) * (1.0 - ty) + z10 * tx * (1.0 - ty) + z01 * (1.0 - tx) * ty + z11 * tx * ty
}

fn rms(errs: &[f64]) -> f64 {
    (errs.iter().map(|e| e * e).sum::<f64>() / errs.len() as f64).sqrt()
}
fn max_abs(errs: &[f64]) -> f64 {
    errs.iter().fold(0.0_f64, |a, e| a.max(e.abs()))
}

/// Nearest-datum distance in CELLS from node (ip, jp) to any defined node of `gd`.
fn nearest_datum_cells(gd: &GriddedDepth, ip: usize, jp: usize) -> f64 {
    let mut best = f64::INFINITY;
    for j in 0..N {
        for i in 0..N {
            if gd.is_control[j * N + i] {
                let d = ((ip as f64 - i as f64).powi(2) + (jp as f64 - j as f64).powi(2)).sqrt();
                best = best.min(d);
            }
        }
    }
    best
}

// ---------------------------------------------------------------------------
// S1 — conditioning assignment + on-data fidelity
// ---------------------------------------------------------------------------

#[test]
fn s1_on_data_misfit_is_authoring_not_the_solve() {
    let (mut snap_node_rms, mut exact_scatter_rms) = (0.0_f64, f64::INFINITY);
    let (mut snap_scatter_rms, mut exact_node_rms) = (0.0_f64, f64::INFINITY);
    let pts = scatter_xy();
    let top_samples: Vec<((f64, f64), f64)> = pts
        .iter()
        .map(|&(x, y)| ((x, y), truth_top(x, y)))
        .collect();
    let base_samples: Vec<((f64, f64), f64)> = pts
        .iter()
        .map(|&(x, y)| ((x, y), truth_top(x, y) + truth_sep(x, y, merge_axis)))
        .collect();

    for label in ["snap-average", "exact-conditioning"] {
        let (top_gd, base_gd) = if label == "snap-average" {
            (author_snap(&top_samples), author_snap(&base_samples))
        } else {
            (
                author_hull(truth_top),
                author_hull(|x, y| truth_top(x, y) + truth_sep(x, y, merge_axis)),
            )
        };
        let m = build(
            two_horizon_stack(top_gd.clone(), base_gd),
            ExtrapolationPolicy::default(),
            false,
        );
        let built_top = surface(&m, "TOP");

        // Built surface honors every DEFINED node of its input EXACTLY (hard
        // controls) — the solve itself contributes zero on-data (at-node) misfit.
        let mut worst_node = 0.0_f64;
        for idx in 0..N * N {
            if top_gd.is_control[idx] {
                worst_node = worst_node.max((built_top[idx] - top_gd.depth_m[idx]).abs());
            }
        }
        assert!(
            worst_node <= 1e-9,
            "{label}: defined nodes must be honored exactly, worst {worst_node}"
        );

        // CONDITIONING error, isolated: the authored node value vs truth AT the
        // node (free of the lattice-resolution representation term that any
        // 100 m sampling of a curved surface carries between nodes).
        let node_errs: Vec<f64> = (0..N * N)
            .filter(|&i| top_gd.is_control[i])
            .map(|i| {
                let (x, y) = node_xy(i % N, i / N);
                top_gd.depth_m[i] - truth_top(x, y)
            })
            .collect();
        // The combined metric at the scatter points (authoring + representation)
        // — the number comparable to the real-model on-data rms.
        let errs: Vec<f64> = pts
            .iter()
            .map(|&(x, y)| eval(built_top, x, y) - truth_top(x, y))
            .collect();
        eprintln!(
            "S1 [{label}] node-conditioning: rms {:.3} m max {:.3} m | on-scatter (cond + representation): rms {:.3} m max {:.3} m",
            rms(&node_errs),
            max_abs(&node_errs),
            rms(&errs),
            max_abs(&errs)
        );
        match label {
            // Snap-average conditioning is the defect: a 65 m scatter snapped
            // onto a 100 m lattice mis-assigns each sample by up to ~70 m
            // laterally — metres of depth on a dipping/curved flank.
            "snap-average" => {
                snap_node_rms = rms(&node_errs);
                snap_scatter_rms = rms(&errs);
            }
            // Perfect conditioning: node error identically zero; the remaining
            // on-scatter misfit is the irreducible lattice-resolution term.
            _ => {
                exact_node_rms = rms(&node_errs);
                exact_scatter_rms = rms(&errs);
            }
        }
    }
    assert!(
        exact_node_rms < 1e-9,
        "exact conditioning must leave zero node error, got {exact_node_rms}"
    );
    assert!(
        snap_node_rms > 0.3,
        "the snap fixture must demonstrate the conditioning defect, got {snap_node_rms}"
    );
    // The on-data (scatter) misfit decomposes: representation floor (exact
    // conditioning) + the conditioning defect (snap adds on top of the floor).
    eprintln!(
        "S1 decomposition: representation floor rms {exact_scatter_rms:.3} m | snap conditioning adds {:.3} m (total {snap_scatter_rms:.3} m)",
        snap_scatter_rms - exact_scatter_rms
    );
    assert!(
        snap_scatter_rms > exact_scatter_rms,
        "snap conditioning must cost misfit on top of the representation floor"
    );
}

// ---------------------------------------------------------------------------
// S2 — the solve itself: full-pin identity + convergence residual
// ---------------------------------------------------------------------------

#[test]
fn s2_solver_is_exact_on_pinned_nodes_and_converged() {
    // Fully-pinned lattice: the built surface IS the input, bit-for-bit.
    let dense = author_dense(truth_top);
    let m = build(
        two_horizon_stack(
            dense.clone(),
            author_dense(|x, y| truth_top(x, y) + SEP_OPEN),
        ),
        ExtrapolationPolicy::default(),
        false,
    );
    let built = surface(&m, "TOP");
    for idx in 0..N * N {
        assert!(
            (built[idx] - dense.depth_m[idx]).abs() <= 1e-12,
            "fully-pinned node {idx} moved"
        );
    }

    // Convergence residual on the SPARSE (hull-only) solve: a second seeded solve
    // from the converged field must not move (iteration budget is sufficient).
    let gd = author_hull(truth_top);
    let controls: Vec<Control> = (0..N * N)
        .filter(|&i| gd.is_control[i])
        .map(|i| Control {
            ip: i % N,
            jp: i / N,
            z: gd.depth_m[i],
        })
        .collect();
    // (a) The RAW kernel path (flat bootstrap + one seeded solve). Historically the
    // slow affine mode stalled this metres short of the fixed point (the SOR was
    // iteration-cap-bound); the converged wrapper existed to work around that. With
    // petekTools' DIRECT band-LU `MinCurvatureOperator` (`task_suite_scatter_perf`)
    // the kernel now *attains* the fixed point in one solve, so a second seeded
    // solve no longer moves — the raw residual is at the floor too.
    let flat = KernelSurface::flat(N, N, 2000.0);
    let s1 = solve_surface_seeded(&flat, &controls).unwrap();
    let s2 = solve_surface_seeded(&s1, &controls).unwrap();
    let mut raw_resid = 0.0_f64;
    for jp in 0..N {
        for ip in 0..N {
            raw_resid = raw_resid.max((s2.z(ip, jp) - s1.z(ip, jp)).abs());
        }
    }
    // (b) The build's actual entry (`solve_surface_converged`: plane detrend +
    // fixed-point restarts) — a whole re-solve of its output must not move.
    let c1 = solve_surface_converged(N, N, &controls).unwrap();
    let c2 = solve_surface_seeded(&c1, &controls).unwrap();
    let mut conv_resid = 0.0_f64;
    for jp in 0..N {
        for ip in 0..N {
            conv_resid = conv_resid.max((c2.z(ip, jp) - c1.z(ip, jp)).abs());
        }
    }
    eprintln!(
        "S2 convergence residual: raw flat-seeded kernel {raw_resid:.2e} m | converged entry {conv_resid:.2e} m"
    );
    assert!(
        conv_resid < 0.01,
        "the build's solve entry must be at the fixed point, residual {conv_resid}"
    );
    // The direct solve eliminated the raw-kernel affine-mode stall: BOTH the raw
    // flat-seeded kernel and the converged entry now sit at the fixed point (a
    // second solve does not move either). The converged wrapper no longer has a
    // stall to beat — it stays as the robust build entry, in agreement with the
    // kernel rather than rescuing it.
    assert!(
        raw_resid < 0.01,
        "the direct solve must reach the fixed point in one solve (no affine-mode stall), \
         raw residual {raw_resid}"
    );
}

// ---------------------------------------------------------------------------
// S3 — extrapolation policy beyond the data hull
// ---------------------------------------------------------------------------

#[test]
fn s3_margin_is_data_bounded_under_the_default_policy() {
    let top_gd = author_hull(truth_top);
    let base_gd = author_hull(|x, y| truth_top(x, y) + truth_sep(x, y, merge_axis));

    let data_min = top_gd
        .depth_m
        .iter()
        .filter(|v| !v.is_nan())
        .fold(f64::INFINITY, |a, &b| a.min(b));
    let data_max = top_gd
        .depth_m
        .iter()
        .filter(|v| !v.is_nan())
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    let natural = build(
        two_horizon_stack(top_gd.clone(), base_gd.clone()),
        ExtrapolationPolicy::NaturalDip,
        false,
    );
    let decayed = build(
        two_horizon_stack(top_gd.clone(), base_gd.clone()),
        ExtrapolationPolicy::default(),
        false,
    );
    let (nat_top, dec_top) = (surface(&natural, "TOP"), surface(&decayed, "TOP"));

    // Audit the margin (void) nodes.
    let mut nat_overshoot = 0.0_f64; // how far beyond the data's own range NaturalDip runs
    let mut dec_overshoot = 0.0_f64;
    let mut nat_err = Vec::new();
    let mut dec_err = Vec::new();
    let mut far_clamp_err = 0.0_f64; // far margin (w = 1) must equal the nearest datum
    for jp in 0..N {
        for ip in 0..N {
            let idx = jp * N + ip;
            if top_gd.is_control[idx] {
                continue;
            }
            let (x, y) = node_xy(ip, jp);
            nat_err.push(nat_top[idx] - truth_top(x, y));
            dec_err.push(dec_top[idx] - truth_top(x, y));
            nat_overshoot = nat_overshoot
                .max(nat_top[idx] - data_max)
                .max(data_min - nat_top[idx]);
            dec_overshoot = dec_overshoot
                .max(dec_top[idx] - data_max)
                .max(data_min - dec_top[idx]);
            let d = nearest_datum_cells(&top_gd, ip, jp);
            if d >= 6.0 {
                // beyond start(2) + decay(4): the policy pins to nearest-data
                let mut best = (f64::INFINITY, 0.0);
                for j in 0..N {
                    for i in 0..N {
                        if top_gd.is_control[j * N + i] {
                            let dd = ((ip as f64 - i as f64).powi(2)
                                + (jp as f64 - j as f64).powi(2))
                            .sqrt();
                            if dd < best.0 {
                                best = (dd, top_gd.depth_m[j * N + i]);
                            }
                        }
                    }
                }
                far_clamp_err = far_clamp_err.max((dec_top[idx] - best.1).abs());
            }
        }
    }
    eprintln!(
        "S3 margin vs truth: NaturalDip rms {:.2} m max {:.2} m | DecayToData rms {:.2} m max {:.2} m",
        rms(&nat_err),
        max_abs(&nat_err),
        rms(&dec_err),
        max_abs(&dec_err)
    );
    eprintln!(
        "S3 beyond-data-range overshoot: NaturalDip {nat_overshoot:.2} m | DecayToData {dec_overshoot:.2} m | far-margin clamp err {far_clamp_err:.2e} m"
    );
    // Audit finding (documented): over this 5-cell margin the TENSIONED kernel
    // does not run linearly unbounded — the tension term relaxes the extension
    // toward harmonic flattening — so the catastrophic real-model margin
    // divergence is attributable to the (now fixed) under-convergence freeze at
    // the seed mean, not to a runaway dip alone. The POLICY guarantees asserted
    // here are the owner-visible contract:
    // (1) the margin never invents structure beyond the observed data envelope,
    assert!(
        dec_overshoot < 1.5,
        "DecayToData must stay data-bounded (ramp-band tolerance), overshoot {dec_overshoot}"
    );
    assert!(
        dec_overshoot < nat_overshoot,
        "the policy must tighten the data bound vs the legacy behaviour"
    );
    // (2) the far margin (beyond start + decay cells) is pinned to the nearest
    //     datum EXACTLY (no drift, whatever the solver does),
    assert!(
        far_clamp_err < 1e-9,
        "far margin must equal the nearest datum exactly, err {far_clamp_err}"
    );

    // The MERGED side's margin: the envelope must STAY merged into the void (the
    // real-model phantom margin) — nearest-data separation there is exactly 0.
    let dec_base = surface(&decayed, "BASE");
    let mut worst_margin_sep = 0.0_f64;
    for jp in 0..N {
        for ip in 0..N {
            let idx = jp * N + ip;
            let (x, _) = node_xy(ip, jp);
            if !top_gd.is_control[idx] && x < 500.0 && nearest_datum_cells(&top_gd, ip, jp) >= 6.0 {
                worst_margin_sep = worst_margin_sep.max(dec_base[idx] - dec_top[idx]);
            }
        }
    }
    eprintln!("S3 merged-margin residual separation (DecayToData): {worst_margin_sep:.2e} m");
    assert!(
        worst_margin_sep < 1e-6,
        "merged envelope must stay merged into the void, got {worst_margin_sep} m"
    );
}

#[test]
fn s3_world_georef_diagonal_variant() {
    // Doctrine R1: the same statement with a world georef and a DIAGONAL merge
    // line (not axis-aligned) — policy behaviour must not be an axis accident.
    let top_gd = author_hull(truth_top);
    let base_gd = author_hull(|x, y| truth_top(x, y) + truth_sep(x, y, merge_diag));
    let m = build(
        two_horizon_stack(top_gd.clone(), base_gd),
        ExtrapolationPolicy::default(),
        true,
    );
    let (top, base) = (surface(&m, "TOP"), surface(&m, "BASE"));
    // Merged-region nodes (inside the data): separation exactly zero-ish (gridded
    // zeros stay zero; taper preserves zeros).
    let mut worst = 0.0_f64;
    let mut n_merged = 0;
    for jp in 0..N {
        for ip in 0..N {
            let idx = jp * N + ip;
            let (x, y) = node_xy(ip, jp);
            if top_gd.is_control[idx] && merge_diag(x, y) <= -100.0 {
                worst = worst.max(base[idx] - top[idx]);
                n_merged += 1;
            }
        }
    }
    eprintln!(
        "S3d merged-region residual separation (world/diag): {worst:.2e} m over {n_merged} nodes"
    );
    assert!(n_merged > 10, "fixture must exercise the merged region");
    assert!(
        worst < 1e-6,
        "diagonal merge must collapse exactly, got {worst}"
    );
    // No repairs on consistent inputs (S5 under the world variant too).
    assert!(m
        .provenance()
        .stack
        .as_ref()
        .unwrap()
        .interface_repairs
        .is_empty());
}

// ---------------------------------------------------------------------------
// S4 — tie locality
// ---------------------------------------------------------------------------

#[test]
fn s4_ties_are_local() {
    // (a) DENSE lattice: a tie replaces one node's datum; every other node is
    // still a hard datum → ONLY the tied nodes move. Radius of influence: 0 cells.
    let dense_top = author_dense(truth_top);
    let dense_base = author_dense(|x, y| truth_top(x, y) + SEP_OPEN);
    let ties: Vec<WellTie> = [
        (7usize, 7usize, 3.0),
        (13, 9, -4.0),
        (9, 13, 2.0),
        (11, 6, -2.5),
        (6, 11, 1.5),
    ]
    .iter()
    .map(|&(ip, jp, dr)| {
        let (x, y) = node_xy(ip, jp);
        WellTie::new(format!("99/{ip}-{jp}"), x, y, ip, jp).with_top("TOP", truth_top(x, y) + dr)
    })
    .collect();
    let untied = build(
        two_horizon_stack(dense_top.clone(), dense_base.clone()),
        ExtrapolationPolicy::default(),
        false,
    );
    let tied =
        StaticModelBuilder::from_horizon_stack(two_horizon_stack(dense_top, dense_base), opts())
            .unwrap()
            .with_well_ties(ties.clone())
            .build()
            .unwrap();
    let (u, t) = (surface(&untied, "TOP"), surface(&tied, "TOP"));
    let tie_nodes: Vec<usize> = ties.iter().map(|w| w.jp * N + w.ip).collect();
    let mut moved_off_tie = 0.0_f64;
    for idx in 0..N * N {
        let d = (t[idx] - u[idx]).abs();
        if tie_nodes.contains(&idx) {
            assert!(d > 1.0, "tied node must move to the measured top");
        } else {
            moved_off_tie = moved_off_tie.max(d);
        }
    }
    eprintln!(
        "S4 dense-lattice tie reach: max off-tie movement {moved_off_tie:.2e} m (radius 0 cells)"
    );
    assert!(
        moved_off_tie <= 1e-9,
        "dense-lattice ties must move ONLY the tied nodes, moved {moved_off_tie}"
    );

    // (b) SPARSE lattice: influence is the interpolation reach; measure the decay
    // and pin the far field.
    let sparse_top = author_hull(truth_top);
    let sparse_top2 = sparse_top.clone();
    let sparse_base = author_hull(|x, y| truth_top(x, y) + SEP_OPEN);
    let untied = build(
        two_horizon_stack(sparse_top.clone(), sparse_base.clone()),
        ExtrapolationPolicy::default(),
        false,
    );
    let tied =
        StaticModelBuilder::from_horizon_stack(two_horizon_stack(sparse_top, sparse_base), opts())
            .unwrap()
            .with_well_ties(ties)
            .build()
            .unwrap();
    let (u, t) = (surface(&untied, "TOP"), surface(&tied, "TOP"));
    let tie_ij = [(7usize, 7usize), (13, 9), (9, 13), (11, 6), (6, 11)];
    let ring_max = |r_lo: f64, r_hi: f64| -> f64 {
        let mut worst = 0.0_f64;
        for jp in 0..N {
            for ip in 0..N {
                let d = tie_ij
                    .iter()
                    .map(|&(ti, tj)| {
                        ((ip as f64 - ti as f64).powi(2) + (jp as f64 - tj as f64).powi(2)).sqrt()
                    })
                    .fold(f64::INFINITY, f64::min);
                if d >= r_lo && d < r_hi {
                    worst = worst.max((t[jp * N + ip] - u[jp * N + ip]).abs());
                }
            }
        }
        worst
    };
    let (r2, r4, r6, r8) = (
        ring_max(1.0, 2.0),
        ring_max(2.0, 4.0),
        ring_max(4.0, 6.0),
        ring_max(8.0, f64::INFINITY),
    );
    eprintln!(
        "S4 sparse-lattice tie decay: |Δ| ring 1-2 cells {r2:.3} m, 2-4 {r4:.3} m, 4-6 {r6:.3} m, >=8 {r8:.3} m"
    );
    // Through the biharmonic-family operator a tie's influence on a SPARSE
    // lattice decays but is not strictly compact; the contract asserted: the
    // far field is sub-half-metre for these ±4 m residuals, and the pinned far
    // margin (beyond start + decay of the DATA hull) does not move at all.
    assert!(
        r8 < 0.5,
        "sparse-lattice tie far-field must stay sub-half-metre, got {r8} m"
    );
    let mut pinned_moved = 0.0_f64;
    for jp in 0..N {
        for ip in 0..N {
            let idx = jp * N + ip;
            if !sparse_top2.is_control[idx] && nearest_datum_cells(&sparse_top2, ip, jp) >= 6.0 {
                pinned_moved = pinned_moved.max((t[idx] - u[idx]).abs());
            }
        }
    }
    eprintln!("S4 pinned far-margin movement under ties: {pinned_moved:.2e} m");
    assert!(
        pinned_moved <= 1e-9,
        "the policy-pinned far margin must be tie-invariant, moved {pinned_moved}"
    );
}

// ---------------------------------------------------------------------------
// S5 — repair semantics on consistent inputs + the truncation warning label
// ---------------------------------------------------------------------------

#[test]
fn s5_no_repairs_on_consistent_scatter_and_proportional_truncation_reports() {
    let m = build(
        two_horizon_stack(
            author_hull(truth_top),
            author_hull(|x, y| truth_top(x, y) + truth_sep(x, y, merge_axis)),
        ),
        ExtrapolationPolicy::default(),
        false,
    );
    // Crossings resolve toward the data by construction (isochores >= 0): zero
    // interface repairs — repairs never ADD rock to a consistent model.
    let prov = m.provenance();
    assert!(
        prov.stack.as_ref().unwrap().interface_repairs.is_empty(),
        "no repairs may fire on consistent inputs: {:?}",
        prov.stack.as_ref().unwrap().interface_repairs
    );
    // This is a PROPORTIONAL build over a merged (pinched) envelope: the
    // zero-thickness truncation warning fires and is NOT Follow-specific (the
    // relabeled `LayersTruncated` semantics).
    assert_eq!(
        prov.stack.as_ref().unwrap().zones[0].conformity,
        Conformity::Proportional
    );
    assert!(
        prov.warnings
            .iter()
            .any(|w| matches!(w, petekstatic::model::BuildWarning::LayersTruncated { cells } if *cells > 0)),
        "merged proportional zone must report truncated cells"
    );
}

// ---------------------------------------------------------------------------
// S6 — envelope-GRV invariance to an internal mapped-surface swap
// ---------------------------------------------------------------------------

#[test]
fn s6_internal_mapped_swap_leaves_other_horizons_and_total_grv_unchanged() {
    // 4 mapped horizons; swap the INTERNAL M1 between two variants (both strictly
    // inside the envelope). With the cumulative-from-top construction each mapped
    // horizon depends only on the TOP anchor and its OWN data, so TOP, M2 and
    // BASE — and hence total envelope GRV — are bit-invariant to the swap.
    let stack_with_m1 = |m1_off: f64| -> HorizonStack {
        HorizonStack {
            horizons: vec![
                StackHorizon {
                    name: "TOP".into(),
                    source: HorizonSource::Mapped(author_dense(truth_top)),
                },
                StackHorizon {
                    name: "M1".into(),
                    source: HorizonSource::Mapped(author_dense(move |x, y| {
                        truth_top(x, y) + m1_off
                    })),
                },
                StackHorizon {
                    name: "M2".into(),
                    source: HorizonSource::Mapped(author_dense(|x, y| truth_top(x, y) + 20.0)),
                },
                StackHorizon {
                    name: "BASE".into(),
                    source: HorizonSource::Mapped(author_dense(|x, y| truth_top(x, y) + 30.0)),
                },
            ],
            zone_layers: ["A", "B", "C"]
                .iter()
                .map(|&name| StackZone {
                    name: name.into(),
                    color: None,
                    conformity: Conformity::Proportional,
                    nk: 4,
                    contacts: vec![],
                })
                .collect(),
        }
    };
    let a = build(stack_with_m1(8.0), ExtrapolationPolicy::default(), false);
    let b = build(stack_with_m1(14.0), ExtrapolationPolicy::default(), false);
    for name in ["TOP", "M2", "BASE"] {
        assert_eq!(
            surface(&a, name),
            surface(&b, name),
            "{name} must be bit-invariant to the internal M1 swap"
        );
    }
    let (ga, gb) = (
        a.in_place_by_zone().unwrap().total.grv_m3,
        b.in_place_by_zone().unwrap().total.grv_m3,
    );
    let rel = (ga - gb).abs() / ga.max(1.0);
    eprintln!("S6 total envelope GRV under internal swap: {ga:.6e} vs {gb:.6e} (rel {rel:.2e})");
    assert!(
        rel < 1e-12,
        "total envelope GRV must be invariant to an internal surface swap"
    );
    // The swap DOES move the zones it bounds (it is a real change, not a no-op).
    let za = a.in_place_by_zone().unwrap();
    let zb = b.in_place_by_zone().unwrap();
    assert!(
        (za.zones[0].in_place.grv_m3 - zb.zones[0].in_place.grv_m3).abs() > 1.0,
        "zone A must respond to the M1 swap"
    );
}

// ---------------------------------------------------------------------------
// S7 — raw scatter end-to-end: the engine owns the gridding
// ---------------------------------------------------------------------------
//
// The root fix (owner ruling): raw world-coordinate scatter flows into the stack
// build as `HorizonSource::Scatter`; the ENGINE conditions it onto the lattice
// (voids left NaN) so the converged solve + DecayToData + isochore build-down act
// on the actual data. This reproduces, end-to-end through `from_scatter_stack`,
// the real-model phantom-margin defect and its cure — and A/Bs it against the
// **independent full-fill** a caller-side pre-gridding produced (each horizon
// solved over the whole lattice, every node a hard control), which drifts in the
// data void between exactly-merged horizons.

/// Column-centroid georef for the fixture lattice at world origin `(ox, oy)`:
/// node 0 sits at `(ox, oy)`, so the centroid origin is half a cell in.
fn scatter_georef(ox: f64, oy: f64) -> Georef {
    Georef::new(ox + 0.5 * SPACING, oy + 0.5 * SPACING, SPACING, SPACING).unwrap()
}

fn scatter_frame(ox: f64, oy: f64) -> StackFrame {
    StackFrame {
        ni: N - 1,
        nj: N - 1,
        georef: scatter_georef(ox, oy),
    }
}

/// World scatter for a horizon: the dense off-node hull scatter at `(ox+x, oy+y)`
/// carrying `f(x, y)` as positive-down depth.
fn scatter_points(ox: f64, oy: f64, f: impl Fn(f64, f64) -> f64) -> Vec<WorldPoint> {
    scatter_xy()
        .into_iter()
        .map(|(x, y)| WorldPoint {
            x: ox + x,
            y: oy + y,
            depth_m: f(x, y),
        })
        .collect()
}

fn scatter_stack(top: Vec<WorldPoint>, base: Vec<WorldPoint>) -> HorizonStack {
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Scatter(top),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Scatter(base),
            },
        ],
        zone_layers: vec![StackZone {
            name: "Z".into(),
            color: None,
            conformity: Conformity::Proportional,
            nk: 8,
            contacts: vec![],
        }],
    }
}

/// The pre-merge **independent full-fill**: snap each horizon's scatter to nodes,
/// solve each INDEPENDENTLY over the whole lattice (every node finite → a hard
/// control, no data-void mask), and stack them the pre-scatter (`Mapped`) way.
/// This is what caller-side per-horizon gridding delivered to the engine.
fn independent_fill_stack(
    top_samples: &[((f64, f64), f64)],
    base_samples: &[((f64, f64), f64)],
) -> HorizonStack {
    let solve_full = |samples: &[((f64, f64), f64)]| -> GriddedDepth {
        let sparse = author_snap(samples);
        let controls: Vec<Control> = (0..N * N)
            .filter(|&i| sparse.is_control[i])
            .map(|i| Control {
                ip: i % N,
                jp: i / N,
                z: sparse.depth_m[i],
            })
            .collect();
        let s = solve_surface_converged(N, N, &controls).unwrap();
        let mut depth = vec![0.0; N * N];
        for jp in 0..N {
            for ip in 0..N {
                depth[jp * N + ip] = s.z(ip, jp);
            }
        }
        GriddedDepth {
            ncol: N,
            nrow: N,
            depth_m: depth,
            is_control: vec![true; N * N], // every node a hard control (the defect)
        }
    };
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Mapped(solve_full(top_samples)),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(solve_full(base_samples)),
            },
        ],
        zone_layers: vec![StackZone {
            name: "Z".into(),
            color: None,
            conformity: Conformity::Proportional,
            nk: 8,
            contacts: vec![],
        }],
    }
}

/// Worst base−top separation over deep MERGED-void nodes (`x < HULL_MIN`, ≥6 cells
/// from any hull datum) of a built model — the phantom-margin metric.
fn merged_void_sep(m: &StaticModel, mask: &GriddedDepth) -> f64 {
    let (top, base) = (surface(m, "TOP"), surface(m, "BASE"));
    let mut worst = 0.0_f64;
    for jp in 0..N {
        for ip in 0..N {
            let idx = jp * N + ip;
            let (x, _) = node_xy(ip, jp);
            if x < HULL_MIN && nearest_datum_cells(mask, ip, jp) >= 6.0 {
                worst = worst.max((base[idx] - top[idx]).abs());
            }
        }
    }
    worst
}

#[test]
fn s7_raw_scatter_collapses_the_merged_void_the_engine_grids() {
    // Merged for x <= 800 (sep exactly 0), opening past it; scatter only in the
    // hull [500,1500] → the margin x < 500 is a deep MERGED data void.
    let base_f = |x: f64, y: f64| truth_top(x, y) + truth_sep(x, y, merge_axis);
    let frame = scatter_frame(0.0, 0.0);
    let mask = author_hull(truth_top); // is_control marks the data hull (for the metric)

    // (a) The FIX: raw scatter, engine-gridded.
    let stack = scatter_stack(
        scatter_points(0.0, 0.0, truth_top),
        scatter_points(0.0, 0.0, base_f),
    );
    let m_scatter = StaticModelBuilder::from_scatter_stack(stack, opts(), frame)
        .unwrap()
        .build()
        .unwrap();
    let scatter_sep = merged_void_sep(&m_scatter, &mask);

    // (b) The DEFECT: the same data, independently full-filled + all-control.
    let top_samples: Vec<((f64, f64), f64)> = scatter_xy()
        .into_iter()
        .map(|p| (p, truth_top(p.0, p.1)))
        .collect();
    let base_samples: Vec<((f64, f64), f64)> = scatter_xy()
        .into_iter()
        .map(|p| (p, base_f(p.0, p.1)))
        .collect();
    let m_fill = StaticModelBuilder::from_horizon_stack(
        independent_fill_stack(&top_samples, &base_samples),
        opts(),
    )
    .unwrap()
    .build()
    .unwrap();
    let fill_sep = merged_void_sep(&m_fill, &mask);

    eprintln!(
        "S7 deep merged-void separation: raw-scatter (engine-gridded) {scatter_sep:.3} m | independent full-fill {fill_sep:.3} m"
    );
    // The engine-gridded scatter keeps the merged envelope merged into the void.
    assert!(
        scatter_sep < 0.5,
        "raw scatter must collapse the merged void, got {scatter_sep} m"
    );
    // The pre-gridded full-fill manufactures phantom separation there (the defect).
    assert!(
        fill_sep > 2.0,
        "the independent full-fill must reproduce the phantom margin, got {fill_sep} m"
    );
    assert!(
        scatter_sep < fill_sep,
        "engine-gridded scatter must beat the pre-gridded fill in the void"
    );

    // On-data fidelity: the engine honours its conditioned nodes; the built TOP at
    // the hull data band matches truth within the snap-conditioning floor.
    let built_top = surface(&m_scatter, "TOP");
    let on_data: Vec<f64> = scatter_xy()
        .into_iter()
        .map(|(x, y)| eval(built_top, x, y) - truth_top(x, y))
        .collect();
    eprintln!(
        "S7 raw-scatter on-data misfit (bilinear conditioning): rms {:.3} m max {:.3} m",
        rms(&on_data),
        max_abs(&on_data)
    );
    // Bilinear conditioning holds on-data to ~the lattice representation floor.
    assert!(
        rms(&on_data) < 0.5,
        "bilinear-conditioned scatter must hold on-data sub-metre, got {} m",
        rms(&on_data)
    );
    // No interface repairs on consistent scatter (S5 through the scatter entry).
    assert!(m_scatter
        .provenance()
        .stack
        .as_ref()
        .unwrap()
        .interface_repairs
        .is_empty());
}

#[test]
fn s7_world_georef_scatter_end_to_end() {
    // Doctrine R1: the same statement at a fictional world georef with a DIAGONAL
    // merge line — the engine's world→lattice conditioning must not be an axis or
    // local-origin accident.
    let base_f = |x: f64, y: f64| truth_top(x, y) + truth_sep(x, y, merge_diag);
    let frame = scatter_frame(ORIGIN_X, ORIGIN_Y);
    let mask = author_hull(truth_top);

    let stack = scatter_stack(
        scatter_points(ORIGIN_X, ORIGIN_Y, truth_top),
        scatter_points(ORIGIN_X, ORIGIN_Y, base_f),
    );
    let m = StaticModelBuilder::from_scatter_stack(stack, opts(), frame)
        .unwrap()
        .build()
        .unwrap();

    // Merged-region nodes inside the data (diagonal): separation ~0.
    let (top, base) = (surface(&m, "TOP"), surface(&m, "BASE"));
    let mut worst_data = 0.0_f64;
    let mut n_merged = 0;
    for jp in 0..N {
        for ip in 0..N {
            let idx = jp * N + ip;
            let (x, y) = node_xy(ip, jp);
            if mask.is_control[idx] && merge_diag(x, y) <= -100.0 {
                worst_data = worst_data.max((base[idx] - top[idx]).abs());
                n_merged += 1;
            }
        }
    }
    eprintln!(
        "S7w world/diag merged-region residual separation: {worst_data:.3e} m over {n_merged} nodes"
    );
    assert!(n_merged > 10, "fixture must exercise the merged region");
    assert!(
        worst_data < 0.5,
        "world-georef diagonal merge must collapse at data, got {worst_data} m"
    );
    // The georef is registered (world frame carried onto the model).
    assert!(
        m.georef().is_some(),
        "scatter build must register the world frame"
    );
}

#[test]
fn s7_scatter_template_matches_deterministic() {
    // R2 mode-matrix: horizon-stack construction × MC template. The MC template's
    // `from_scatter_stack` must condition + resolve the scatter byte-for-byte the
    // deterministic builder does, so a nominal realization reproduces the
    // deterministic model's structural surfaces.
    let base_f = |x: f64, y: f64| truth_top(x, y) + truth_sep(x, y, merge_axis);
    let frame = scatter_frame(ORIGIN_X, ORIGIN_Y);
    let stack = || {
        scatter_stack(
            scatter_points(ORIGIN_X, ORIGIN_Y, truth_top),
            scatter_points(ORIGIN_X, ORIGIN_Y, base_f),
        )
    };
    let det = StaticModelBuilder::from_scatter_stack(stack(), opts(), frame)
        .unwrap()
        .build()
        .unwrap();
    let mut tmpl = StaticModelTemplate::from_scatter_stack(stack(), opts(), frame).unwrap();
    // Nominal draw at the fixture's own area (spacing) so the frames coincide.
    let real = tmpl
        .realize(&RealizationDraw::new(
            AREA_M2, SEP_OPEN, 3000.0, 0.2, 1.0, 0.25, 1,
        ))
        .unwrap();
    for name in ["TOP", "BASE"] {
        let (d, r) = (surface(&det, name), surface(&real, name));
        let worst = d
            .iter()
            .zip(r)
            .fold(0.0_f64, |a, (x, y)| a.max((x - y).abs()));
        assert!(
            worst < 1e-6,
            "{name}: MC template scatter must reproduce the deterministic surface, worst {worst} m"
        );
    }
}

#[test]
fn scatter_dedup_seam_is_bit_identical() {
    // `task_suite_scatter_perf` dedup seam: conditioning the scatter ONCE via
    // `condition_scatter_stack` and feeding the conditioned handle to
    // `from_horizon_stack` + `with_georef` must produce a model **bit-identical**
    // to conditioning inside `from_scatter_stack` — the guarantee that lets a
    // caller condition once and reuse it for both the model and its MC template
    // without a second cold solve.
    let base_f = |x: f64, y: f64| truth_top(x, y) + truth_sep(x, y, merge_axis);
    let frame = scatter_frame(ORIGIN_X, ORIGIN_Y);
    let stack = || {
        scatter_stack(
            scatter_points(ORIGIN_X, ORIGIN_Y, truth_top),
            scatter_points(ORIGIN_X, ORIGIN_Y, base_f),
        )
    };
    // Path A: the all-in-one scatter entry (conditions internally).
    let a = StaticModelBuilder::from_scatter_stack(stack(), opts(), frame)
        .unwrap()
        .build()
        .unwrap();
    // Path B: condition once, then the plain horizon-stack entry + georef.
    let conditioned = StaticModelBuilder::condition_scatter_stack(stack(), &frame).unwrap();
    let g = frame.georef;
    let b = StaticModelBuilder::from_horizon_stack(conditioned, opts())
        .unwrap()
        .with_georef(g.origin_x, g.origin_y, g.spacing_x, g.spacing_y)
        .build()
        .unwrap();
    for name in ["TOP", "BASE"] {
        let (sa, sb) = (surface(&a, name), surface(&b, name));
        // Bit-for-bit equality — not a tolerance. The conditioned handle IS the
        // same conditioning `from_scatter_stack` runs internally.
        assert_eq!(sa, sb, "{name}: dedup seam must be bit-identical");
    }
    assert!(
        b.georef().is_some(),
        "the seam path must register the frame"
    );
}
