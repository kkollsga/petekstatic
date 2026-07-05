//! Out-of-core pipeline acceptance (rulings R2/R3/R4/R5,
//! `petekSuite/dev-docs/designs/out-of-core-strategy.md`): the memory-budget mode
//! switch, in-core↔spilled tolerance parity, the loud switch + cleanup, and the
//! spilled-mode structured-MC determinism + stale-buffer contracts.

use petekstatic::gridder::{Conformity, SolveOpts};
use petekstatic::model::{
    run_mc, spill_grid_to, BuildOpts, ConstantPriors, Input, McInputs, McSettings, MemoryBudget,
    StaticModelBuilder, StaticModelTemplate,
};
use petekstatic::volumetrics::{compute_clipped, Clip};
use petekstatic::wireframe::{
    Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
};

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

/// A single-contact flat model, `ni×nj×nk`, top at 2500 m, OWC below base — the
/// whole column is hydrocarbon (a clean parity fixture, moderate absolute depth so
/// the f32 quantization stays small).
fn build(ni: usize, nj: usize, nk: usize, budget: MemoryBudget) -> petekstatic::model::StaticModel {
    StaticModelBuilder::flat(ni, nj, 2500.0, 9000.0, opts(nk))
        .unwrap()
        .with_memory_budget(budget)
        .build()
        .unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    if b == 0.0 {
        a.abs()
    } else {
        (a - b).abs() / b.abs()
    }
}

#[test]
fn below_budget_stays_in_core_and_is_byte_identical() {
    // A generous budget keeps the in-core path; two builds are bit-identical, and
    // the model is NOT spilled (zero behaviour change, R5).
    let a = build(30, 30, 10, MemoryBudget::unlimited());
    let b = build(30, 30, 10, MemoryBudget::bytes(4 * 1024 * 1024 * 1024));
    assert!(!a.is_spilled(), "generous budget must stay in-core");
    assert!(!b.is_spilled());
    let ia = a.in_place().unwrap();
    let ib = b.in_place().unwrap();
    // Same in-core code path → bit-identical volumes.
    assert_eq!(
        ia.grv_m3.to_bits(),
        ib.grv_m3.to_bits(),
        "GRV not byte-identical"
    );
    assert_eq!(
        ia.hcpv_m3.to_bits(),
        ib.hcpv_m3.to_bits(),
        "HCPV not byte-identical"
    );
    assert_eq!(a.bulk_volume().to_bits(), b.bulk_volume().to_bits());
}

#[test]
fn above_budget_spills_loudly_with_a_store_path() {
    let path;
    {
        let m = build(30, 30, 10, MemoryBudget::bytes(1024));
        assert!(m.is_spilled(), "tiny budget must spill");
        let p = m
            .spill_store_path()
            .expect("spilled model exposes its store path");
        assert!(
            p.exists(),
            "spill store file must exist while the model is alive"
        );
        path = p.to_path_buf();
    }
    // Drop semantics (R5): the temp store is removed when the last clone drops.
    assert!(!path.exists(), "spill store must be cleaned up on drop");
}

#[test]
fn in_core_vs_spilled_parity_within_measured_bound() {
    // The R4 tolerance-parity assertion: f32 storage lanes change volumes at a
    // small relative bound; accumulations stay f64. Measure it, print it, assert it.
    let core = build(40, 40, 20, MemoryBudget::unlimited());
    let spill = build(40, 40, 20, MemoryBudget::bytes(1024));
    assert!(!core.is_spilled() && spill.is_spilled());

    let ic = core.in_place().unwrap();
    let is = spill.in_place().unwrap();
    let g = rel(is.grv_m3, ic.grv_m3);
    let h = rel(is.hcpv_m3, ic.hcpv_m3);
    let bulk = rel(spill.bulk_volume(), core.bulk_volume());
    eprintln!(
        "in-core↔spilled parity (40×40×20, top 2500 m): GRV {g:.2e}, HCPV {h:.2e}, bulk {bulk:.2e}"
    );
    // Documented bound (measured): f32 ZCORN at ~2500 m over ~2 m layers keeps the
    // relative volume error well under 1e-5. Asserted with margin.
    let bound = 1e-5;
    assert!(g <= bound, "GRV parity {g:.2e} exceeds {bound:.0e}");
    assert!(h <= bound, "HCPV parity {h:.2e} exceeds {bound:.0e}");
    assert!(bulk <= bound, "bulk parity {bulk:.2e} exceeds {bound:.0e}");
    // Cell counts (in the column) are integer-identical — the geometry topology
    // is unchanged, only the stored precision narrows.
    assert_eq!(is.cells_in_column, ic.cells_in_column);
}

#[test]
fn slab_incremental_build_is_bit_identical_to_build_then_spill() {
    // v2 item 1: the forced-spill build now streams ZCORN + cubes k-slab-by-k-slab
    // straight into the store (never a whole in-core grid). Prove it produces the
    // EXACT same store as the v1 build-then-spill would: the streaming ZCORN uses the
    // same f64 `boundary_depth` narrowed to f32, and the constant cubes are the same
    // f32, so a streamed spilled model's volumetrics are BIT-IDENTICAL to spilling a
    // fully-built in-core grid. (Above the tolerance parity vs in-core — this is the
    // exactness of the new producer against the old one.)
    let (ni, nj, nk) = (24, 24, 12);
    let contact = 9000.0; // OWC below base — whole column, single contact
    let core = build(ni, nj, nk, MemoryBudget::unlimited());
    assert!(!core.is_spilled());
    // v1 reference: spill the fully-built in-core grid (build-then-spill).
    let dir = std::env::temp_dir();
    let refpath = dir.join(format!("petekstatic-bitref-{}.pts", std::process::id()));
    let backing = spill_grid_to(core.grid(), &refpath, true).unwrap();
    let reference =
        compute_clipped(&backing.source(), Clip::Single(contact), 0..nk, false).unwrap();

    // v2: the streaming spilled build.
    let streamed = build(ni, nj, nk, MemoryBudget::bytes(1024));
    assert!(streamed.is_spilled());
    let got = streamed.in_place_summary().unwrap();

    assert_eq!(
        got.grv_m3.to_bits(),
        reference.grv_m3.to_bits(),
        "GRV: streamed build not bit-identical to build-then-spill"
    );
    assert_eq!(
        got.hcpv_m3.to_bits(),
        reference.hcpv_m3.to_bits(),
        "HCPV: streamed build not bit-identical to build-then-spill"
    );
    assert_eq!(got.cells_in_column, reference.cells_in_column);
    assert_eq!(
        streamed.bulk_volume().to_bits(),
        backing.bulk_volume().unwrap().to_bits()
    );
}

#[test]
fn spilled_in_place_full_materializes_per_cell() {
    let spill = build(20, 20, 8, MemoryBudget::bytes(1024));
    let ip = spill.in_place().unwrap();
    assert_eq!(ip.per_cell_hcpv.len(), spill.dims().cell_count());
    let summ = spill.in_place_summary().unwrap();
    assert!(
        summ.per_cell_hcpv.is_empty(),
        "summary path leaves per-cell empty"
    );
    assert_eq!(
        ip.grv_m3.to_bits(),
        summ.grv_m3.to_bits(),
        "full == summary aggregates"
    );
}

// --- two-contact spilled parity (gas cap + oil rim) ---

fn two_contact_wf(n: usize, top_m: f64) -> Wireframe {
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
                depth_m: vec![top_m; n * n],
                is_control: vec![true; n * n],
            },
        }]),
        contacts: vec![
            Contact {
                kind: ContactKind::Goc,
                depth_m: top_m + 12.0,
                hardness: Hardness::Hard,
            },
            Contact {
                kind: ContactKind::Owc,
                depth_m: top_m + 28.0,
                hardness: Hardness::Hard,
            },
        ],
    }
}

#[test]
fn spilled_two_contact_parity() {
    let wf = two_contact_wf(21, 2500.0);
    let core = StaticModelBuilder::from_wireframe(&wf, opts(20))
        .unwrap()
        .with_memory_budget(MemoryBudget::unlimited())
        .build()
        .unwrap();
    let spill = StaticModelBuilder::from_wireframe(&wf, opts(20))
        .unwrap()
        .with_memory_budget(MemoryBudget::bytes(1024))
        .build()
        .unwrap();
    assert!(spill.is_spilled());
    let ic = core.in_place().unwrap();
    let is = spill.in_place().unwrap();
    let (gc, gs) = (ic.gas.unwrap(), is.gas.unwrap());
    let (oc, os) = (ic.oil.unwrap(), is.oil.unwrap());
    assert_eq!(gs.cells, gc.cells, "gas-cap cell count identical");
    assert_eq!(os.cells, oc.cells, "oil-leg cell count identical");
    assert!(rel(gs.hcpv_m3, gc.hcpv_m3) <= 1e-5, "gas HCPV parity");
    assert!(rel(os.hcpv_m3, oc.hcpv_m3) <= 1e-5, "oil HCPV parity");
}

// --- spilled structured MC (R3/R4) ---

fn mc_wf(n: usize, top_m: f64, owc_m: f64) -> Wireframe {
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
                depth_m: vec![top_m; n * n],
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

fn tri(min: f64, mode: f64, max: f64) -> Input {
    Input::plain(petektools::sampling::Sampler::new_triangular(min, mode, max).unwrap())
}

fn mc_inputs() -> McInputs {
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
fn spilled_mc_is_deterministic_and_matches_in_core_within_tolerance() {
    let wf = mc_wf(11, 2500.0, 2565.0);
    let inputs = mc_inputs();
    let (n, seed) = (60usize, 7u64);

    // In-core reference.
    let mut t0 = StaticModelTemplate::new(&wf, opts(20)).unwrap();
    let core = run_mc(&mut t0, &inputs, &McSettings::new(n, seed)).unwrap();

    // Spilled serial: bit-reproducible within the spilled mode.
    let mut ta = StaticModelTemplate::new(&wf, opts(20)).unwrap();
    let mut tb = StaticModelTemplate::new(&wf, opts(20)).unwrap();
    let spill_settings = McSettings::new(n, seed).with_spill_dir(std::env::temp_dir());
    let sa = run_mc(&mut ta, &inputs, &spill_settings).unwrap();
    let sb = run_mc(&mut tb, &inputs, &spill_settings).unwrap();
    assert_eq!(sa.oil_sm3, sb.oil_sm3, "spilled MC not bit-reproducible");

    // Sharded spilled == serial spilled at every worker count (the R3 determinism
    // contract, in the spilled mode).
    for workers in [1usize, 2, 3, 5] {
        let mut t = StaticModelTemplate::new(&wf, opts(20)).unwrap();
        let par = run_mc(
            &mut t,
            &inputs,
            &spill_settings.clone().with_workers(workers),
        )
        .unwrap();
        assert_eq!(
            par.oil_sm3, sa.oil_sm3,
            "spilled workers={workers}: oil diverged from serial"
        );
        assert_eq!(
            par.grv_m3, sa.grv_m3,
            "spilled workers={workers}: grv diverged"
        );
    }

    // Spilled ≈ in-core within the f32 tolerance (per-draw, worst relative).
    let worst = core
        .oil_sm3
        .iter()
        .zip(&sa.oil_sm3)
        .map(|(c, s)| rel(*s, *c))
        .fold(0.0f64, f64::max);
    eprintln!("spilled↔in-core MC oil worst per-draw relative error: {worst:.2e}");
    assert!(
        worst <= 1e-5,
        "spilled MC oil diverged from in-core by {worst:.2e}"
    );
}

#[test]
fn spilled_realize_into_has_no_stale_buffer() {
    // WIDE-range draws stress the reused-model + reused-store loop: consecutive
    // draws differ a lot, so any stale ZCORN/cube slab left over from the prior
    // draw would corrupt the next. Each spilled draw must still equal its in-core
    // counterpart within the f32 tolerance — the R3 realize_into stale-buffer proof
    // in the spilled mode.
    let wf = mc_wf(11, 2500.0, 2565.0);
    let wide = McInputs::new(
        tri(0.5e6, 1.0e6, 1.5e6),    // area (wide)
        tri(20.0, 40.0, 60.0),       // gross (wide)
        tri(2555.0, 2565.0, 2575.0), // contact (wide)
        tri(0.15, 0.22, 0.30),       // porosity
        tri(0.55, 0.75, 0.95),       // ntg
        tri(0.15, 0.30, 0.45),       // sw
        tri(1.20, 1.30, 1.45),       // boi
    );
    let (n, seed) = (40usize, 3u64);
    let mut tc = StaticModelTemplate::new(&wf, opts(20)).unwrap();
    let core = run_mc(&mut tc, &wide, &McSettings::new(n, seed)).unwrap();
    let mut ts = StaticModelTemplate::new(&wf, opts(20)).unwrap();
    let spill = run_mc(
        &mut ts,
        &wide,
        &McSettings::new(n, seed).with_spill_dir(std::env::temp_dir()),
    )
    .unwrap();

    // Draws genuinely vary (the loop is stressed, not a constant).
    let (lo, hi) = core
        .oil_sm3
        .iter()
        .fold((f64::MAX, f64::MIN), |(l, h), &v| (l.min(v), h.max(v)));
    assert!(
        hi > lo * 2.0,
        "draws must vary widely to stress the reuse loop"
    );

    // Every spilled draw matches its in-core counterpart → no stale carryover.
    for (i, (c, s)) in core.oil_sm3.iter().zip(&spill.oil_sm3).enumerate() {
        assert!(
            rel(*s, *c) <= 1e-5,
            "draw {i}: spilled {s} != in-core {c} (stale buffer?)"
        );
    }
}
