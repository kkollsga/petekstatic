//! The R7 **conformance battery** — the engine half (testing doctrine R7,
//! `task_petekstatic_spec_mirror`): every config-layer value type round-trips
//! through serde (`from_json(to_json(x)) == x`, value equality) and behaves as a
//! **value** (`with_*` derives a new value; the original is unchanged), and the
//! ONE declarative [`BuildSpec`] drives the builder and the template to
//! **bit-identical** results vs the equivalent `with_*` chains — so the spec
//! consolidation cannot drift the determinism contracts.
//!
//! ## The McInputs slot (pending upstream)
//! `McInputs` is deliberately **absent** from the round-trip registry: it wraps
//! `petektools::sampling::Sampler` / `Clamped`, which do not yet derive serde
//! (petekTools is an upstream DAG dep — the serde derives are queued there by the
//! coordinator, `task_petektools_specs`). When that lands, add `McInputs` to
//! `round_trips_value_equal` exactly like the other entries.
//!
//! Mode-matrix note (R2): these are config-layer value tests — mode-independent
//! by construction; the mode matrix lives in `mode_matrix.rs`.

use petekstatic::gridder::{Conformity, ExtrapolationPolicy, SolveOpts};
use petekstatic::model::{
    BuildOpts, BuildSpec, ConstantPriors, Georef, HorizonSource, HorizonStack, McSettings, Pick,
    RealizationDraw, StackFrame, StackHorizon, StackZone, StaticModel, StaticModelBuilder,
    StaticModelTemplate, StructuralPerturbation, TieMethod, TieSettings, WellTie, WorldPoint,
    ZoneDraw,
};
use petekstatic::wireframe::{Contact, ContactKind, GriddedDepth, Hardness};
use serde::{de::DeserializeOwned, Serialize};

/// Round-trip `v` through JSON and assert **value equality** (R7 rule 1).
fn round_trip<T: Serialize + DeserializeOwned + PartialEq + std::fmt::Debug>(v: &T) {
    let json = serde_json::to_string(v).expect("serializes");
    let back: T = serde_json::from_str(&json).expect("deserializes");
    assert_eq!(*v, back, "round-trip must be value-equal:\n{json}");
}

const N: usize = 7;

fn flat_surf(depth: f64) -> GriddedDepth {
    GriddedDepth {
        ncol: N,
        nrow: N,
        depth_m: vec![depth; N * N],
        is_control: vec![true; N * N],
    }
}

fn owc(depth_m: f64) -> Contact {
    Contact {
        kind: ContactKind::Owc,
        depth_m,
        hardness: Hardness::Hard,
    }
}

/// A fully-defined (no-`NaN` — JSON has no NaN) 3-horizon / 2-zone stack that
/// exercises every [`HorizonSource`] variant.
fn fixture_stack() -> HorizonStack {
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Mapped(flat_surf(2500.0)),
            },
            StackHorizon {
                name: "MID".into(),
                source: HorizonSource::TopsOnly(vec![Pick {
                    ip: 3,
                    jp: 3,
                    depth_m: 2520.0,
                }]),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(flat_surf(2540.0)),
            },
        ],
        zone_layers: vec![
            StackZone::new("UPPER", Conformity::Proportional, 4, vec![owc(2535.0)])
                .with_color("#ffcc00"),
            StackZone::new("LOWER", Conformity::FollowTop { dz_m: 2.0 }, 4, Vec::new()),
        ],
    }
}

fn fixture_ties() -> Vec<WellTie> {
    vec![WellTie::new("99/1-1", 300.0, 300.0, 3, 3).with_top("TOP", 2504.0)]
}

fn fixture_opts() -> BuildOpts {
    BuildOpts {
        area_m2: 360_000.0, // side 600 → dx = dy = 100 over 6 cells
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

fn fixture_spec() -> BuildSpec {
    BuildSpec::new()
        .with_inputs_ref("conformance-fixture")
        .with_georef(431_000.0, 6_521_000.0, 100.0, 100.0)
        .with_boundary(vec![
            [431_000.0, 6_521_000.0],
            [431_600.0, 6_521_000.0],
            [431_600.0, 6_521_600.0],
            [431_000.0, 6_521_600.0],
            [431_000.0, 6_521_000.0],
        ])
        .with_extrapolation(ExtrapolationPolicy::DecayToData {
            start_cells: 1.0,
            decay_cells: 3.0,
        })
        .with_clamp_base_to_top(true)
        .with_min_thickness_m(0.5)
        .with_collapse_below_m(0.05)
        .with_sugar_cube(true)
        .with_sw_gas(0.12)
        .with_well_ties(fixture_ties())
        .with_tie_settings(TieSettings::radius(250.0))
}

// ---------------------------------------------------------------------------
// R7 rule 1 — round-trip == value equality, on every config type
// ---------------------------------------------------------------------------

#[test]
fn round_trips_value_equal() {
    // Scalars / small values.
    round_trip(&Pick {
        ip: 2,
        jp: 5,
        depth_m: 2513.25,
    });
    round_trip(&WorldPoint {
        x: 431_250.0,
        y: 6_521_125.5,
        depth_m: 2507.75,
    });
    round_trip(&Georef::new(431_000.0, 6_521_000.0, 50.0, 62.5).unwrap());
    round_trip(&fixture_ties()[0]);
    round_trip(&TieSettings::replace());
    round_trip(&TieSettings::radius(1500.0));
    assert_eq!(
        TieSettings::default().method,
        TieMethod::Replace,
        "the default tie method is datum substitution (today's behaviour)"
    );
    round_trip(&fixture_opts());

    // The horizon-stack family, every HorizonSource variant (incl. raw Scatter).
    round_trip(&HorizonSource::Scatter(vec![
        WorldPoint {
            x: 431_100.0,
            y: 6_521_050.0,
            depth_m: 2502.0,
        },
        WorldPoint {
            x: 431_400.0,
            y: 6_521_450.0,
            depth_m: 2511.0,
        },
    ]));
    round_trip(&fixture_stack());
    round_trip(&StackFrame {
        ni: 6,
        nj: 6,
        georef: Georef::new(431_000.0, 6_521_000.0, 100.0, 100.0).unwrap(),
    });

    // The draws.
    round_trip(&StructuralPerturbation {
        control_shifts: vec![(3, 3, -12.5), (0, 6, 4.0)],
    });
    round_trip(&ZoneDraw::new(1).with_owc(2535.0).with_goc(2510.0));
    round_trip(
        &RealizationDraw::new(360_000.0, 40.0, 2535.0, 0.22, 0.85, 0.28, 7)
            .with_goc(2510.0)
            .with_sw_gas(0.1)
            .with_property_shift("PORO", 0.015)
            .with_zone_draw(ZoneDraw::new(0).with_priors(0.21, 0.8, 0.3))
            .with_structural(StructuralPerturbation {
                control_shifts: vec![(1, 1, 2.0)],
            }),
    );

    // The one declarative build spec + the one MC run-settings value.
    round_trip(&fixture_spec());
    round_trip(
        &McSettings::new(500, 42)
            .with_workers(4)
            .with_spill_dir("/tmp/mc-conformance"),
    );

    // McInputs: deliberately absent — blocked on petekTools Sampler/Clamped serde
    // (see the module docs). Slot it in here when the upstream derives land.
}

// ---------------------------------------------------------------------------
// R7 rule 2 — value semantics: `with_*` derives a NEW value, original unchanged
// ---------------------------------------------------------------------------

#[test]
fn with_sugar_derives_new_values_original_unchanged() {
    let base = BuildSpec::new();
    let derived = base.clone().with_min_thickness_m(1.0).with_sugar_cube(true);
    assert_eq!(base, BuildSpec::new(), "original spec unchanged");
    assert_ne!(base, derived, "derived spec compares unequal");

    let ms = McSettings::new(100, 1);
    let ms2 = ms.clone().with_workers(4);
    assert_eq!(ms, McSettings::new(100, 1));
    assert_ne!(ms, ms2);

    let zd = ZoneDraw::new(0);
    let zd2 = zd.clone().with_owc(2535.0);
    assert_eq!(zd, ZoneDraw::new(0));
    assert_ne!(zd, zd2);

    let sz = StackZone::new("Z", Conformity::Proportional, 4, Vec::new());
    let sz2 = sz.clone().with_color("#0af");
    assert_eq!(
        sz,
        StackZone::new("Z", Conformity::Proportional, 4, Vec::new())
    );
    assert_ne!(sz, sz2);
}

// ---------------------------------------------------------------------------
// The ONE BuildSpec: with_spec == the with_* chain, bit-identically, on BOTH
// consumers (the determinism guard for the spec consolidation)
// ---------------------------------------------------------------------------

/// Bit-compare two models: per-cell geometry + every cube + provenance label.
fn assert_bit_identical(a: &StaticModel, b: &StaticModel, ctx: &str) {
    assert_eq!(a.grid().cell_count(), b.grid().cell_count(), "{ctx}: cells");
    for lin in 0..a.grid().cell_count() {
        assert_eq!(
            a.grid().cell_centroid_z_at(lin).to_bits(),
            b.grid().cell_centroid_z_at(lin).to_bits(),
            "{ctx}: ZCORN centroid-z diverged at cell {lin}"
        );
    }
    for name in a.property_names() {
        assert_eq!(
            a.property(name).map(|p| &p.values),
            b.property(name).map(|p| &p.values),
            "{ctx}: cube {name}"
        );
    }
    assert_eq!(
        a.provenance().inputs_ref,
        b.provenance().inputs_ref,
        "{ctx}: inputs_ref"
    );
    assert_eq!(
        a.provenance().sugar_cube,
        b.provenance().sugar_cube,
        "{ctx}: sugar_cube"
    );
    assert_eq!(
        a.provenance().well_ties.len(),
        b.provenance().well_ties.len(),
        "{ctx}: tie records"
    );
}

/// The `with_*` chain equivalent of [`fixture_spec`], applied to a builder.
fn chained_builder() -> StaticModelBuilder {
    StaticModelBuilder::from_horizon_stack(fixture_stack(), fixture_opts())
        .unwrap()
        .with_inputs_ref("conformance-fixture")
        .with_georef(431_000.0, 6_521_000.0, 100.0, 100.0)
        .with_boundary(vec![
            [431_000.0, 6_521_000.0],
            [431_600.0, 6_521_000.0],
            [431_600.0, 6_521_600.0],
            [431_000.0, 6_521_600.0],
            [431_000.0, 6_521_000.0],
        ])
        .with_extrapolation(ExtrapolationPolicy::DecayToData {
            start_cells: 1.0,
            decay_cells: 3.0,
        })
        .with_clamp_base_to_top(true)
        .with_min_thickness_m(0.5)
        .with_collapse_below_m(0.05)
        .with_sugar_cube(true)
        .with_sw_gas(0.12)
        .with_tie_settings(TieSettings::radius(250.0))
        .with_well_ties(fixture_ties())
}

#[test]
fn builder_with_spec_is_bit_identical_to_the_chain() {
    let via_chain = chained_builder().build().unwrap();
    let via_spec = StaticModelBuilder::from_horizon_stack(fixture_stack(), fixture_opts())
        .unwrap()
        .with_spec(fixture_spec())
        .build()
        .unwrap();
    assert_bit_identical(&via_chain, &via_spec, "builder chain vs spec");
    assert_eq!(via_spec.provenance().inputs_ref, "conformance-fixture");
}

#[test]
fn template_with_spec_is_bit_identical_to_the_chain() {
    let draw = RealizationDraw::new(360_000.0, 0.0, 0.0, 0.2, 0.8, 0.3, 3)
        .with_zone_draw(ZoneDraw::new(0).with_owc(2535.0));

    let mut via_chain = StaticModelTemplate::from_horizon_stack(fixture_stack(), fixture_opts())
        .unwrap()
        .with_inputs_ref("conformance-fixture")
        .with_georef(431_000.0, 6_521_000.0, 100.0, 100.0)
        .with_boundary(vec![
            [431_000.0, 6_521_000.0],
            [431_600.0, 6_521_000.0],
            [431_600.0, 6_521_600.0],
            [431_000.0, 6_521_600.0],
            [431_000.0, 6_521_000.0],
        ])
        .with_clamp_base_to_top(true)
        .with_min_thickness_m(0.5)
        .with_collapse_below_m(0.05)
        .with_sugar_cube(true)
        .with_tie_settings(TieSettings::radius(250.0))
        .with_extrapolation(ExtrapolationPolicy::DecayToData {
            start_cells: 1.0,
            decay_cells: 3.0,
        })
        .unwrap()
        .with_well_ties(fixture_ties())
        .unwrap();

    let mut via_spec = StaticModelTemplate::from_horizon_stack(fixture_stack(), fixture_opts())
        .unwrap()
        .with_spec(fixture_spec())
        .unwrap();

    let a = via_chain.realize(&draw).unwrap();
    let b = via_spec.realize(&draw).unwrap();
    assert_bit_identical(&a, &b, "template chain vs spec");

    // And the template realization matches the deterministic builder's tie
    // geometry: the spec is ONE value consumed by both (the mirror's point).
    let built = chained_builder().build().unwrap();
    let hz = |m: &StaticModel| m.framework().horizons[0].surface.depth_m.clone();
    let (hb, ht) = (hz(&built), hz(&a));
    for (i, (x, y)) in hb.iter().zip(&ht).enumerate() {
        assert!(
            (x - y).abs() < 1e-9,
            "tied TOP diverged between builder and template at node {i}: {x} vs {y}"
        );
    }
}

// ---------------------------------------------------------------------------
// TieSettings — Replace (default) vs Radius (bounded locality)
// ---------------------------------------------------------------------------

/// A dense flat 2-horizon stack for the tie-locality tests: 11×11 nodes,
/// dx = dy = 100 m (area 1 km²), TOP 2500 / BASE 2540.
const TN: usize = 11;

fn tie_stack() -> HorizonStack {
    let flat = |d: f64| GriddedDepth {
        ncol: TN,
        nrow: TN,
        depth_m: vec![d; TN * TN],
        is_control: vec![true; TN * TN],
    };
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Mapped(flat(2500.0)),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(flat(2540.0)),
            },
        ],
        zone_layers: vec![StackZone::new("Z", Conformity::Proportional, 4, Vec::new())],
    }
}

fn tie_opts() -> BuildOpts {
    let mut o = fixture_opts();
    o.area_m2 = 1_000_000.0; // side 1000 → dx = dy = 100 over 10 cells
    o
}

fn tie_build(settings: TieSettings) -> StaticModel {
    StaticModelBuilder::from_horizon_stack(tie_stack(), tie_opts())
        .unwrap()
        .with_tie_settings(settings)
        .with_well_ties(vec![
            WellTie::new("99/2-1", 500.0, 500.0, 5, 5).with_top("TOP", 2510.0)
        ])
        .build()
        .unwrap()
}

fn top_depth(m: &StaticModel, ip: usize, jp: usize) -> f64 {
    m.framework().horizons[0].surface.depth_m[jp * TN + ip]
}

#[test]
fn tie_replace_is_the_default_and_moves_exactly_the_tied_node() {
    // Replace on a fully-defined lattice: every other node is a hard datum, so
    // the tie moves exactly the tied node (radius of influence 0 cells).
    let untied = StaticModelBuilder::from_horizon_stack(tie_stack(), tie_opts())
        .unwrap()
        .build()
        .unwrap();
    let replaced = tie_build(TieSettings::replace());
    // The default IS Replace: an explicit setting changes nothing.
    let defaulted = StaticModelBuilder::from_horizon_stack(tie_stack(), tie_opts())
        .unwrap()
        .with_well_ties(vec![
            WellTie::new("99/2-1", 500.0, 500.0, 5, 5).with_top("TOP", 2510.0)
        ])
        .build()
        .unwrap();
    assert_bit_identical(&replaced, &defaulted, "explicit Replace vs default");

    assert!(
        (top_depth(&replaced, 5, 5) - 2510.0).abs() < 1e-9,
        "tie node on measured"
    );
    for &(ip, jp) in &[(4usize, 5usize), (5, 4), (6, 5), (5, 6), (0, 0), (10, 10)] {
        assert!(
            (top_depth(&replaced, ip, jp) - top_depth(&untied, ip, jp)).abs() < 1e-12,
            "Replace must leave every other node untouched (node {ip},{jp})"
        );
    }
}

#[test]
fn tie_radius_decays_the_residual_and_is_bounded() {
    // Radius 250 m over dx = 100 m: the +10 residual decays linearly — full at
    // the well, 10·(1 − d/250) inside, ZERO beyond 250 m (bit-untouched).
    let untied = StaticModelBuilder::from_horizon_stack(tie_stack(), tie_opts())
        .unwrap()
        .build()
        .unwrap();
    let tied = tie_build(TieSettings::radius(250.0));

    assert!(
        (top_depth(&tied, 5, 5) - 2510.0).abs() < 1e-9,
        "tie node on measured"
    );
    // One cell away (100 m): w = 0.6 → 2506.
    for &(ip, jp) in &[(4usize, 5usize), (6, 5), (5, 4), (5, 6)] {
        assert!(
            (top_depth(&tied, ip, jp) - 2506.0).abs() < 1e-9,
            "100 m from the well: expected 2506, got {}",
            top_depth(&tied, ip, jp)
        );
    }
    // Two cells (200 m): w = 0.2 → 2502.
    assert!((top_depth(&tied, 7, 5) - 2502.0).abs() < 1e-9);
    // Diagonal neighbour (141.4 m): w = 1 − √2·100/250.
    let w_diag = 1.0 - (2.0f64).sqrt() * 100.0 / 250.0;
    assert!((top_depth(&tied, 6, 6) - (2500.0 + 10.0 * w_diag)).abs() < 1e-9);
    // Beyond the radius (300 m +): bit-untouched vs the untied build.
    for &(ip, jp) in &[(8usize, 5usize), (5, 8), (2, 5), (0, 0), (10, 10)] {
        assert!(
            (top_depth(&tied, ip, jp) - top_depth(&untied, ip, jp)).abs() < 1e-12,
            "beyond radius must be untouched (node {ip},{jp})"
        );
    }
    // The BASE keeps honouring its own datums (the isochore construction).
    let base = |m: &StaticModel, ip: usize, jp: usize| {
        m.framework().horizons[1].surface.depth_m[jp * TN + ip]
    };
    assert!(
        (base(&tied, 5, 5) - 2540.0).abs() < 1e-9,
        "BASE stays on its datum"
    );

    // A degenerate radius is a typed error, not silence.
    let err = StaticModelBuilder::from_horizon_stack(tie_stack(), tie_opts())
        .unwrap()
        .with_tie_settings(TieSettings::radius(0.0))
        .with_well_ties(vec![
            WellTie::new("99/2-1", 500.0, 500.0, 5, 5).with_top("TOP", 2510.0)
        ])
        .build()
        .unwrap_err();
    assert!(
        err.to_string().contains("radius"),
        "loud degenerate-radius error: {err}"
    );
}

#[test]
fn tie_radius_template_matches_builder() {
    // The template applies the SAME tie math (one authority): its tied TOP
    // equals the builder's, node for node.
    let built = tie_build(TieSettings::radius(250.0));
    let mut tmpl = StaticModelTemplate::from_horizon_stack(tie_stack(), tie_opts())
        .unwrap()
        .with_tie_settings(TieSettings::radius(250.0))
        .with_well_ties(vec![
            WellTie::new("99/2-1", 500.0, 500.0, 5, 5).with_top("TOP", 2510.0)
        ])
        .unwrap();
    let real = tmpl
        .realize(&RealizationDraw::new(
            1_000_000.0,
            0.0,
            0.0,
            0.2,
            0.8,
            0.3,
            1,
        ))
        .unwrap();
    for jp in 0..TN {
        for ip in 0..TN {
            assert!(
                (top_depth(&built, ip, jp) - top_depth(&real, ip, jp)).abs() < 1e-9,
                "template tie diverged from builder at ({ip},{jp})"
            );
        }
    }
}
