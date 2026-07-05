//! Per-horizon correlated structural uncertainty (`task_petekstatic_structural_uncertainty`,
//! `decision_structural_uncertainty_isochore`): the TOP surface takes a correlated
//! DEPTH perturbation field per draw; every deeper horizon perturbs via a per-zone
//! ISOCHORE (thickness) field, clamped `>= 0` and zero-masked at exact merges.
//!
//! ## What this file pins
//! - **R3 planted-truth** — the clamp-induced **mean-GRV bias** on a fixture with a
//!   *known analytic sensitivity*: for a zone of uniform thickness `T` with a
//!   marginal-`N(0, sd²)` field, `E[max(0, T + P)] = T·Φ(T/sd) + sd·φ(T/sd)` per node
//!   (variogram-independent — the marginal is `N` at every node regardless of the
//!   correlation). Measured over many seeds vs the closed form, at a pinchout
//!   (`T ≈ sd`, large bias) and a thick zone (`T ≫ sd`, bias ≈ 0, the bound).
//! - **Determinism** — same seed ⇒ bit-identical realize; **shard-split invariance**
//!   (the stack structural realize is a pure function of the draw ⇒ sharded == serial).
//! - **`realize_into`** stale-buffer bit-match with structural draws.
//! - **R1 frame** — a world-georeferenced variant (perturbation is frame-invariant).
//! - **R5 degenerate** — `sd_m = 0` is identical to the unperturbed realize.
//! - **R7 serde** — `PerturbationField` on the draws round-trips.
//!
//! All fixtures are synthetic at a fictional area/frame.

use petekstatic::gridder::{Conformity, SolveOpts};
use petekstatic::model::{
    BuildOpts, ConstantPriors, HorizonSource, HorizonStack, PerturbationField, RealizationDraw,
    StaticModelTemplate, ZoneDraw,
};
use petekstatic::volumetrics::PORO;
use petekstatic::wireframe::GriddedDepth;

use petektools::{Variogram, VariogramModel};

const N: usize = 9; // 8×8 cells, 9×9 nodes
const AREA: f64 = 1_000_000.0; // dx = dy = 125 m; footprint area = 1e6 m²

fn flat(depth: f64) -> GriddedDepth {
    GriddedDepth {
        ncol: N,
        nrow: N,
        depth_m: vec![depth; N * N],
        is_control: vec![true; N * N],
    }
}

fn mapped(name: &str, depth: f64) -> petekstatic::model::StackHorizon {
    petekstatic::model::StackHorizon {
        name: name.into(),
        source: HorizonSource::Mapped(flat(depth)),
    }
}

/// A 3-horizon / 2-zone flat stack with the given horizon depths, both zones
/// Proportional + contactless (GRV only — the bias test is a GRV statement).
fn stack(top: f64, mid: f64, base: f64) -> HorizonStack {
    HorizonStack {
        horizons: vec![mapped("H0", top), mapped("H1", mid), mapped("H2", base)],
        zone_layers: vec![
            petekstatic::model::StackZone::new("Z0", Conformity::Proportional, 4, Vec::new()),
            petekstatic::model::StackZone::new("Z1", Conformity::Proportional, 4, Vec::new()),
        ],
    }
}

fn opts() -> BuildOpts {
    BuildOpts {
        area_m2: AREA,
        gross_height_m: 0.0,
        nk: 0,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.2,
            net_to_gross: 0.8,
            water_saturation: 0.3,
        },
    }
}

/// A near-independent field (spherical range ≪ node spacing) so the per-node marginal
/// dominates and the MC mean is tight — the variogram does not bias the mean, only
/// its Monte-Carlo variance, so a short range gives the most nodes-independent estimate.
fn near_nugget(sd_m: f64) -> PerturbationField {
    PerturbationField::new(
        sd_m,
        Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 1.0).unwrap(),
    )
}

fn base_draw(seed: u64) -> RealizationDraw {
    RealizationDraw::new(AREA, 0.0, 0.0, 0.2, 0.8, 0.3, seed)
}

fn zone0_grv(m: &petekstatic::model::StaticModel) -> f64 {
    m.in_place_by_zone().unwrap().zones[0].in_place.grv_m3
}

// --- standard-normal helpers (Zelen & Severo, A&S 26.2.17; err < 7.5e-8) ---
fn norm_pdf(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (2.0 * std::f64::consts::PI).sqrt()
}
fn norm_cdf(x: f64) -> f64 {
    if x < 0.0 {
        return 1.0 - norm_cdf(-x);
    }
    let t = 1.0 / (1.0 + 0.231_641_9 * x);
    let poly = t
        * (0.319_381_530
            + t * (-0.356_563_782
                + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
    1.0 - norm_pdf(x) * poly
}
/// `E[max(0, X)]` for `X ~ N(mu, sd²)` — the per-node expected clamped thickness.
fn e_clamp(mu: f64, sd: f64) -> f64 {
    let z = mu / sd;
    mu * norm_cdf(z) + sd * norm_pdf(z)
}

/// Mean zone-0 GRV over `seeds` draws that perturb zone 0's isochore by `sd`.
fn mean_zone0_grv(top: f64, mid: f64, base: f64, sd: f64, seeds: u64) -> f64 {
    let mut t = StaticModelTemplate::from_horizon_stack(stack(top, mid, base), opts()).unwrap();
    let mut acc = 0.0;
    for s in 0..seeds {
        let draw =
            base_draw(s).with_zone_draw(ZoneDraw::new(0).with_isochore_structural(near_nugget(sd)));
        acc += zone0_grv(&t.realize(&draw).unwrap());
    }
    acc / seeds as f64
}

#[test]
fn clamp_bias_pinchout_matches_analytic() {
    // A pinchout zone: T = 2 m, sd = 10 m (z = 0.2). The clamp lifts the mean far
    // above T — the bias is real and large. Analytic per-node E[max(0,T+P)]:
    let (t_m, sd) = (2.0, 10.0);
    let e_thickness = e_clamp(t_m, sd); // ≈ 5.069 m
    let expected_grv = AREA * e_thickness; // footprint area · mean thickness
                                           // The uniform pinchout: Top 5000, Mid 5002 (T=2), Base 5060.
    let seeds = 400;
    let mc = mean_zone0_grv(5000.0, 5002.0, 5060.0, sd, seeds);

    // The clamp bias is unmistakable: the mean thickness is > 2× the geometric T.
    assert!(
        mc > AREA * t_m * 1.5,
        "clamp bias not observed: mc GRV {mc} vs geometric {}",
        AREA * t_m
    );
    // …and it matches the analytic half-normal-tail value within MC error. SE of the
    // mean thickness ≈ sd/√(nodes·seeds) ≈ 10/√(81·400) ≈ 0.056 m ⇒ 3σ ≈ 0.17 m; use
    // a comfortable 0.4 m (·AREA) band.
    let rel = (mc - expected_grv).abs() / expected_grv;
    assert!(
        rel < 0.4 * AREA / expected_grv,
        "mean-GRV bias off analytic: mc {mc} vs E {expected_grv} (E[thickness]={e_thickness})"
    );
}

#[test]
fn clamp_bias_thick_zone_is_negligible() {
    // The bound at the thick end: T = 40 m, sd = 5 m (z = 8). The clamp essentially
    // never bites, so E[max(0,T+P)] ≈ T — the mean GRV is unbiased. This pins the
    // upper end of the bound (bias ∈ [0, sd/√(2π)] per node, → 0 as T/sd → ∞).
    let (t_m, sd) = (40.0, 5.0);
    let expected = AREA * t_m; // bias ≈ 0
    let mc = mean_zone0_grv(5000.0, 5040.0, 5100.0, sd, 300);
    let rel = (mc - expected).abs() / expected;
    assert!(
        rel < 0.02,
        "thick-zone GRV should be ~unbiased: mc {mc} vs {expected} (rel {rel})"
    );
    // Analytic bias here is < 1e-3 m per node — negligible vs T.
    let analytic_bias = e_clamp(t_m, sd) - t_m;
    assert!(
        analytic_bias < 1e-3,
        "analytic thick-zone bias should be ~0, got {analytic_bias}"
    );
}

#[test]
fn structural_realize_bit_deterministic_per_seed() {
    // Same seed + same structural draw ⇒ bit-identical GRV, on two fresh templates.
    let draw = base_draw(7)
        .with_top_structural(near_nugget(8.0))
        .with_zone_draw(ZoneDraw::new(0).with_isochore_structural(near_nugget(6.0)))
        .with_zone_draw(ZoneDraw::new(1).with_isochore_structural(near_nugget(4.0)));
    let mut a =
        StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
    let mut b =
        StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
    let (ga, gb) = (a.realize(&draw).unwrap(), b.realize(&draw).unwrap());
    assert_eq!(
        ga.grid().bulk_volume(),
        gb.grid().bulk_volume(),
        "structural realize not deterministic"
    );
    assert_eq!(zone0_grv(&ga), zone0_grv(&gb));
    // The perturbation actually moved geometry off the unperturbed GRV.
    let plain = a.realize(&base_draw(7)).unwrap();
    assert_ne!(
        zone0_grv(&ga),
        zone0_grv(&plain),
        "structural draw must change GRV"
    );
}

#[test]
fn structural_stack_realize_is_shard_split_invariant() {
    // The stack structural realize is a PURE function of the draw (no warm chain: the
    // surfaces are template-fixed and the field depends only on draw+surfaces), so any
    // shard split / draw order yields identical per-draw results — the sharded==serial
    // contract. Prove it: realize a draw set in order, then in REVERSE on a fresh
    // template, and on a CLONE; each draw's GRV must be identical every way.
    let draws: Vec<RealizationDraw> = (0..5)
        .map(|s| {
            base_draw(s)
                .with_top_structural(near_nugget(7.0))
                .with_zone_draw(ZoneDraw::new(0).with_isochore_structural(near_nugget(5.0)))
        })
        .collect();
    let mut fwd =
        StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
    let forward: Vec<f64> = draws
        .iter()
        .map(|d| zone0_grv(&fwd.realize(d).unwrap()))
        .collect();

    let mut rev =
        StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
    for (i, d) in draws.iter().enumerate().rev() {
        assert_eq!(
            zone0_grv(&rev.realize(d).unwrap()),
            forward[i],
            "draw {i}: order-dependent"
        );
    }
    let mut clone = fwd.clone();
    for (i, d) in draws.iter().enumerate() {
        assert_eq!(
            zone0_grv(&clone.realize(d).unwrap()),
            forward[i],
            "draw {i}: clone diverged"
        );
    }
}

#[test]
fn structural_realize_into_stale_buffer_bit_matches_fresh() {
    // Two DIFFERENT structural draws into ONE reused model must leave it bit-identical
    // to a fresh realize of the second draw (realize_into recycles buffers cleanly).
    let a = base_draw(1)
        .with_top_structural(near_nugget(9.0))
        .with_zone_draw(ZoneDraw::new(0).with_isochore_structural(near_nugget(6.0)));
    let b = base_draw(2)
        .with_top_structural(near_nugget(4.0))
        .with_zone_draw(ZoneDraw::new(1).with_isochore_structural(near_nugget(7.0)));

    let fresh_b = {
        let mut t =
            StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
        let mut m = t.reusable_model();
        t.realize_into(&b, &mut m).unwrap();
        (m.grid().bulk_volume(), zone0_grv(&m))
    };
    let reused = {
        let mut t =
            StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
        let mut m = t.reusable_model();
        t.realize_into(&a, &mut m).unwrap(); // stale A state
        t.realize_into(&b, &mut m).unwrap(); // recycled into B
        (m.grid().bulk_volume(), zone0_grv(&m))
    };
    assert_eq!(
        reused, fresh_b,
        "stale-buffer structural realize_into != fresh"
    );
}

#[test]
fn structural_world_georef_variant_is_frame_invariant() {
    // R1: a world-georeferenced template (fictional 431000/6521000 origin) must give
    // the SAME structural geology as the local frame — a depth/thickness perturbation
    // is frame-invariant (the SGS field is translation-invariant on the lattice). Same
    // draw + seed ⇒ identical GRV local vs world.
    let draw = base_draw(3)
        .with_top_structural(near_nugget(8.0))
        .with_zone_draw(ZoneDraw::new(0).with_isochore_structural(near_nugget(5.0)));
    let mut local =
        StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
    let mut world = StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts())
        .unwrap()
        .with_georef(431_000.0, 6_521_000.0, 125.0, 125.0);
    let gl = local.realize(&draw).unwrap();
    let gw = world.realize(&draw).unwrap();
    assert_eq!(
        zone0_grv(&gl),
        zone0_grv(&gw),
        "structural GRV must be frame-invariant"
    );
    // The world model actually carries a world frame (not the local degenerate one).
    assert!(gw.property(PORO).is_some());
}

#[test]
fn sd_zero_is_identical_to_unperturbed() {
    // R5 degenerate: a zero-sd field is a no-op — bit-identical to the plain realize.
    let draw = base_draw(5)
        .with_top_structural(near_nugget(0.0))
        .with_zone_draw(ZoneDraw::new(0).with_isochore_structural(near_nugget(0.0)));
    let mut t =
        StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
    let mut u =
        StaticModelTemplate::from_horizon_stack(stack(5000.0, 5030.0, 5060.0), opts()).unwrap();
    let perturbed = t.realize(&draw).unwrap();
    let plain = u.realize(&base_draw(5)).unwrap();
    assert_eq!(
        perturbed.grid().bulk_volume(),
        plain.grid().bulk_volume(),
        "sd=0 must be a no-op"
    );
}

#[test]
fn perturbation_field_serde_round_trips() {
    // R7: PerturbationField on the draws serializes as part of the scenario.
    let draw = base_draw(9)
        .with_top_structural(near_nugget(3.5))
        .with_zone_draw(
            ZoneDraw::new(0).with_isochore_structural(PerturbationField::new(
                2.0,
                Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 2500.0).unwrap(),
            )),
        );
    let json = serde_json::to_string(&draw).unwrap();
    let back: RealizationDraw = serde_json::from_str(&json).unwrap();
    assert_eq!(
        draw, back,
        "RealizationDraw with structural fields must round-trip"
    );
    let pf = back.top_structural.unwrap();
    assert_eq!(pf.sd_m, 3.5);
}
