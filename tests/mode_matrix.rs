//! Mode-matrix acceptance (testing doctrine **R2**, `task_petekstatic_test_matrix`,
//! `petekSuite/dev-docs/designs/testing-doctrine.md`). A model-level feature is
//! tested across the mode matrix it claims to support; an unsupported cell is a
//! **documented typed error**, not an untested hole.
//!
//! ## The support matrix (the family convention — declared here, honoured everywhere)
//!
//! The axes a `StaticModel` feature crosses:
//!   - **in-core × spilled** — a resident f64 grid vs the out-of-core, mmap-backed
//!     f32 spill store (above the memory budget).
//!   - **serial × sharded** — single-threaded MC vs the rayon-sharded driver.
//!   - **single-zone × horizon-stack** — one zone vs an ordered multizone stack.
//!   - **wireframe × horizon-stack construction** — `StaticModelTemplate::new`
//!     (single-surface wireframe) vs `from_horizon_stack`.
//!
//! ## Cell coverage (feature → cell → the test that pins it)
//!
//! | feature | in-core | spilled | serial | sharded | single-zone | stack |
//! |---|---|---|---|---|---|---|
//! | volumetrics `in_place`/`summary` | grv.rs units | out_of_core.rs:99,168 | — | — | grv.rs units | multizone_acceptance.rs:924 |
//! | volumetrics `by_zone` | multizone_acceptance.rs:924 | multizone_acceptance.rs:1000 | — | — | n/a | ✓ |
//! | map bundle | map.rs units | **`spilled_map_bundle_matches_in_core`** (this file) | — | — | ✓ | multizone_acceptance.rs:1067 |
//! | section bundle | section.rs units | **`spilled_section_bundle_matches_in_core`** (this file) | — | — | ✓ | multizone_acceptance.rs:1067 |
//! | volume bundle | volume.rs units | template.rs:1573 (spill-invariant) | — | — | ✓ | multizone_acceptance.rs:1067 |
//! | `zone_stats` | model.rs units | **UNSUPPORTED → typed error** (`spilled_zone_stats_is_a_typed_error`, this file) | — | — | ✓ | ✓ |
//! | structured MC | mc.rs:777 | out_of_core.rs:284 | mc.rs:777 | mc.rs:793 | ✓ | — |
//! | zoned MC | template.rs:1677 | out_of_core.rs:284 (flat)† | template.rs:1677 | **`zoned_mc_sharded_matches_serial`** (this file) | — | ✓ |
//! | `realize_into` | lib.rs:1330 | out_of_core.rs:331 | ✓ | ✓ | ✓ | template.rs:1477 |
//!
//! † The spilled MC driver is template-generic — the flat-template spilled coverage
//! (out_of_core.rs) plus the zoned-template serial/sharded coverage here exercise the
//! same `run_structured_mc*` code path over both template kinds; there is no
//! zoned-specific spilled branch. All fixtures are synthetic at a fictional area.

use petekstatic::gridder::{Conformity, SolveOpts};
use petekstatic::model::{
    run_mc, BuildOpts, ConstantPriors, HorizonSource, HorizonStack, Input, MapSpec, McInputs,
    McSettings, MemoryBudget, SectionSpec, StackHorizon, StackZone, StaticModel,
    StaticModelBuilder, StaticModelTemplate,
};
use petekstatic::wireframe::{Contact, ContactKind, GriddedDepth, Hardness};

fn priors() -> ConstantPriors {
    ConstantPriors {
        porosity: 0.25,
        net_to_gross: 0.8,
        water_saturation: 0.3,
    }
}

fn opts(nk: usize) -> BuildOpts {
    BuildOpts {
        area_m2: 1.0e6,
        gross_height_m: 40.0,
        nk,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: priors(),
    }
}

/// A single-contact flat model, built in-core or force-spilled by the budget.
fn flat_model(spilled: bool) -> StaticModel {
    let budget = if spilled {
        MemoryBudget::bytes(1024)
    } else {
        MemoryBudget::unlimited()
    };
    let m = StaticModelBuilder::flat(24, 24, 2500.0, 9000.0, opts(10))
        .unwrap()
        .with_memory_budget(budget)
        .build()
        .unwrap();
    assert_eq!(
        m.is_spilled(),
        spilled,
        "budget did not select the wanted mode"
    );
    m
}

fn rel(a: f64, b: f64) -> f64 {
    if b == 0.0 {
        a.abs()
    } else {
        (a - b).abs() / b.abs()
    }
}

// --- spilled × view bundles: a spilled model materializes its backing for the
// (non-hot-path) exports, so every bundle kind is a SUPPORTED cell, not a
// degenerate 1×1 placeholder read. ---

#[test]
fn spilled_map_bundle_matches_in_core() {
    // Regression: `map_bundle` used to read the spilled model's 1×1×1 placeholder
    // grid, failing with a misleading "1x1 lattice" error on EVERY spilled model
    // (the map sibling of `question_volume_bundle_stack_empty`). It now materializes
    // the backing like the volume/section exports.
    let core = flat_model(false);
    let spill = flat_model(true);
    let spec = MapSpec::new().property("PORO");

    let mc = core.map_bundle(&spec).unwrap();
    let ms = spill
        .map_bundle(&spec)
        .expect("spilled map_bundle must materialize the backing, not error on 1x1");

    // Not a degenerate placeholder: the real areal frame comes through.
    assert!(
        ms.frame.ncol > 1 && ms.frame.nrow > 1,
        "spilled map is degenerate"
    );
    assert_eq!(
        ms.frame.ncol, mc.frame.ncol,
        "areal frame differs from in-core"
    );
    assert_eq!(ms.frame.nrow, mc.frame.nrow);
    assert_eq!(ms.zone_averages.len(), mc.zone_averages.len());

    // The zone-average PORO maps agree within the f32 spill-lane tolerance.
    let za_core = &mc.zone_averages[0].values;
    let za_spill = &ms.zone_averages[0].values;
    assert_eq!(za_core.len(), za_spill.len());
    let worst = za_core
        .iter()
        .zip(za_spill)
        .filter(|(a, b)| a.is_finite() && b.is_finite())
        .map(|(a, b)| rel(*b, *a))
        .fold(0.0f64, f64::max);
    assert!(
        worst <= 1e-5,
        "spilled map PORO diverged from in-core by {worst:.2e}"
    );
}

#[test]
fn spilled_section_bundle_matches_in_core() {
    let core = flat_model(false);
    let spill = flat_model(true);
    // A diagonal fence across the model, authored in the model's own areal frame
    // (real metres for this fixture) so it actually crosses many columns.
    let f = core.map_bundle(&MapSpec::new()).unwrap().frame;
    let (x0, y0) = (
        f.origin_x + 0.5 * f.spacing_x,
        f.origin_y + 0.5 * f.spacing_y,
    );
    let (x1, y1) = (
        f.origin_x + (f.ncol as f64 - 1.5) * f.spacing_x,
        f.origin_y + (f.nrow as f64 - 1.5) * f.spacing_y,
    );
    let spec = SectionSpec::Polyline(vec![[x0, y0], [x1, y1]]);

    let bc = core.intersection_bundle(&spec, Some("PORO")).unwrap();
    let bs = spill
        .intersection_bundle(&spec, Some("PORO"))
        .expect("spilled section must materialize the backing");

    assert!(bc.columns.len() > 1, "in-core section is degenerate");
    assert_eq!(
        bs.columns.len(),
        bc.columns.len(),
        "spilled section column count differs"
    );
    // Column depths are finite and match in-core within the f32 tolerance.
    let mut worst = 0.0f64;
    for (cc, cs) in bc.columns.iter().zip(&bs.columns) {
        for (a, b) in cc.layer_tops.iter().zip(&cs.layer_tops) {
            if a.is_finite() && b.is_finite() {
                worst = worst.max(rel(*b, *a));
            }
        }
    }
    assert!(
        worst <= 1e-4,
        "spilled section tops diverged by {worst:.2e}"
    );
}

// --- spilled × zone_stats: an UNSUPPORTED cell. Per-zone cube statistics window the
// cubes randomly across all k-slabs; the spilled read surface is v1-scoped to the
// volumetric re-cuts. The cell is a documented, message-pinned typed error — never a
// silent wrong answer. (This is the doctrine's template for an unsupported cell.) ---

#[test]
fn spilled_zone_stats_is_a_typed_error() {
    let spill = flat_model(true);
    let err = spill
        .zone_stats("PORO")
        .expect_err("zone_stats on a spilled model must be a typed error in v1");
    let msg = format!("{err}");
    assert!(
        msg.contains("spilled") || msg.contains("out-of-core"),
        "the error must name the unsupported spilled mode, got: {msg}"
    );
    // The in-core cell is supported (contrast): same call succeeds.
    assert!(flat_model(false).zone_stats("PORO").is_ok());
}

// --- zoned (horizon-stack) × sharded MC: the missing MC cell. The sharded driver
// must match the serial driver bit-for-bit at every worker count when driving a
// MULTIZONE template (not just the single-zone flat one). ---

/// A compact flat 2-zone stack with an OWC in the lower zone (so `in_place` resolves).
fn two_zone_stack() -> HorizonStack {
    const N: usize = 11;
    let flat = |d: f64| GriddedDepth {
        ncol: N,
        nrow: N,
        depth_m: vec![d; N * N],
        is_control: vec![true; N * N],
    };
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Mapped(flat(2500.0)),
            },
            StackHorizon {
                name: "MID".into(),
                source: HorizonSource::Mapped(flat(2540.0)),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(flat(2580.0)),
            },
        ],
        zone_layers: vec![
            StackZone {
                name: "UPPER".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 6,
                contacts: vec![],
            },
            StackZone {
                name: "LOWER".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 6,
                contacts: vec![Contact {
                    kind: ContactKind::Owc,
                    depth_m: 2565.0,
                    hardness: Hardness::Hard,
                }],
            },
        ],
    }
}

fn mc_inputs() -> McInputs {
    let tri = |min, mode, max| {
        Input::plain(petektools::sampling::Sampler::new_triangular(min, mode, max).unwrap())
    };
    McInputs::new(
        tri(0.9e6, 1.0e6, 1.1e6),    // area
        tri(35.0, 40.0, 45.0),       // gross
        tri(2560.0, 2565.0, 2570.0), // contact
        tri(0.18, 0.22, 0.26),       // porosity
        tri(0.70, 0.80, 0.90),       // ntg
        tri(0.25, 0.30, 0.35),       // sw
        tri(1.20, 1.30, 1.45),       // boi
    )
}

#[test]
fn zoned_mc_sharded_matches_serial() {
    let inputs = mc_inputs();
    let (n, seed) = (48usize, 11u64);

    let mut ts = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts(12)).unwrap();
    let serial = run_mc(&mut ts, &inputs, &McSettings::new(n, seed)).unwrap();

    for workers in [1usize, 2, 3, 5] {
        let mut tp = StaticModelTemplate::from_horizon_stack(two_zone_stack(), opts(12)).unwrap();
        let par = run_mc(
            &mut tp,
            &inputs,
            &McSettings::new(n, seed).with_workers(workers),
        )
        .unwrap();
        assert_eq!(
            par.oil_sm3, serial.oil_sm3,
            "zoned sharded workers={workers}: oil diverged from serial (non-deterministic)"
        );
        assert_eq!(
            par.grv_m3, serial.grv_m3,
            "zoned sharded workers={workers}: grv diverged from serial"
        );
    }
}
