//! Envelope-subordination acceptance for a tops-only internal split inside a
//! mapped envelope that MERGES (`task_petekstatic_topsonly_envelope`).
//!
//! ## The defect this pins (owner-reported, coordinator-verified on the real model)
//! A zone bounded by two MAPPED horizons that merge — separation exactly 0 in the
//! input, the zone geologically absent there — used to come out of the build with a
//! near-constant PHANTOM thickness ≈ the tops-only internal horizon's pick
//! thickness, across the whole merged region. Mechanism: the tops-only conformal
//! drape applied its pick thickness as an ABSOLUTE offset from the zone top; where
//! that crossed the mapped BASE, the per-interface order-repair pushed the MAPPED
//! base DOWN to sit below the derived internal horizon — derived geometry
//! overriding measured geometry, propagating phantom thickness down the stack.
//!
//! ## The construction that removes it (no special rules — it falls out of the build)
//! `from_horizon_stack` now BUILDS DOWN via non-negative isochores: the top is
//! gridded once, each deeper mapped horizon = the horizon above + a clamped gridded
//! zone isochore (exact-0 where the inputs merge → the zone collapses to genuine
//! zero), and a tops-only internal split = the mapped horizon above +
//! `min(pick isochore, envelope isochore)` — a plain clamp, so envelope 0 ⇒ split 0
//! ⇒ both sub-zones collapse. Mapped horizons stay bit-authoritative; a derived
//! surface can never displace them.
//!
//! Everything here is a hand-authored FICTIONAL fixture (flat/linear surfaces, made-
//! up depths + well nodes) — no dataset content of any kind. Per doctrine R1 the
//! assertions run on both a local build and a world-georeferenced, azimuthed variant.

use petekstatic::grid::Ijk;
use petekstatic::gridder::{Conformity, SolveOpts};
use petekstatic::model::{
    BuildOpts, BuildWarning, ConstantPriors, HorizonSource, HorizonStack, Pick, StackHorizon,
    StackZone, StaticModel, StaticModelBuilder,
};
use petekstatic::wireframe::GriddedDepth;

const N: usize = 9; // 8×8 cells
const AREA_M2: f64 = 1_000_000.0; // 1000 m side
const DX: f64 = 1000.0 / (N - 1) as f64;
const TOP_DEPTH: f64 = 2000.0;
const GROSS: f64 = 30.0; // open-half envelope thickness (m)
const PICK_THICK: f64 = 15.0; // internal split thickness at the picks (m)
const ORIGIN_X: f64 = 700_000.0; // fictional study-area window
const ORIGIN_Y: f64 = 7_100_000.0;

fn opts() -> BuildOpts {
    BuildOpts {
        area_m2: AREA_M2,
        gross_height_m: GROSS,
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

fn gridded(field: Vec<f64>) -> GriddedDepth {
    GriddedDepth {
        ncol: N,
        nrow: N,
        depth_m: field,
        is_control: vec![true; N * N],
    }
}

/// A gently dipping TOP (so the isochore build-down is exercised on a non-constant
/// surface, not just a flat one).
fn top_field() -> Vec<f64> {
    let mut f = vec![0.0; N * N];
    for jp in 0..N {
        for ip in 0..N {
            f[jp * N + ip] = TOP_DEPTH + 2.0 * ip as f64 + 1.0 * jp as f64;
        }
    }
    f
}

/// The zone thickness (isochore) the BASE encodes: exactly `0` where `merged(ip,jp)`
/// (the two mapped horizons coincide — the zone absent), `GROSS` where open.
fn base_field(merged: impl Fn(usize, usize) -> bool) -> Vec<f64> {
    let top = top_field();
    let mut f = vec![0.0; N * N];
    for jp in 0..N {
        for ip in 0..N {
            let t = if merged(ip, jp) { 0.0 } else { GROSS };
            f[jp * N + ip] = top[jp * N + ip] + t;
        }
    }
    f
}

/// The stack: TOP (mapped) → H_MID (tops-only split, picks only in the open half) →
/// BASE (mapped, merging with TOP where `merged`). `pick_nodes` are open-half nodes.
fn split_stack(
    merged: impl Fn(usize, usize) -> bool + Copy,
    pick_nodes: &[(usize, usize)],
) -> HorizonStack {
    let top = top_field();
    let picks = pick_nodes
        .iter()
        .map(|&(ip, jp)| Pick {
            ip,
            jp,
            depth_m: top[jp * N + ip] + PICK_THICK,
        })
        .collect();
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Mapped(gridded(top)),
            },
            StackHorizon {
                name: "MID".into(),
                source: HorizonSource::TopsOnly(picks),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(gridded(base_field(merged))),
            },
        ],
        zone_layers: vec![
            StackZone {
                name: "UPPER".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 8,
                contacts: vec![],
            },
            StackZone {
                name: "LOWER".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 8,
                contacts: vec![],
            },
        ],
    }
}

/// The single-zone envelope (TOP → BASE, no internal split) — the reference the
/// mapped BASE and the total GRV must match bit-for-bit / to the metre.
fn envelope_stack(merged: impl Fn(usize, usize) -> bool) -> HorizonStack {
    let top = top_field();
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Mapped(gridded(top)),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(gridded(base_field(merged))),
            },
        ],
        zone_layers: vec![StackZone {
            name: "ENVELOPE".into(),
            color: None,
            conformity: Conformity::Proportional,
            nk: 16,
            contacts: vec![],
        }],
    }
}

fn build(stack: HorizonStack, georef: bool) -> StaticModel {
    let b = StaticModelBuilder::from_horizon_stack(stack, opts()).unwrap();
    // No min-thickness / clamp opts: the pure default path, so any interface repair
    // that fires is a real crossing — the construction must produce ordered surfaces
    // with NO repair on these consistent inputs.
    let b = if georef {
        b.with_georef(ORIGIN_X, ORIGIN_Y, DX, DX)
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

/// Axis-aligned merge: the left half (`ip < 4`) is merged, the right half open.
fn merged_axis(ip: usize, _jp: usize) -> bool {
    ip < 4
}

/// Diagonal ("azimuthed") merge: the merge line runs NE, so it is not aligned with
/// either lattice axis — `ip + jp < 5` is merged.
fn merged_diag(ip: usize, jp: usize) -> bool {
    ip + jp < 5
}

fn run_case(
    merged: impl Fn(usize, usize) -> bool + Copy,
    pick_nodes: &[(usize, usize)],
    georef: bool,
) {
    let m = build(split_stack(merged, pick_nodes), georef);
    let env = build(envelope_stack(merged), georef);

    let top = surface(&m, "TOP");
    let mid = surface(&m, "MID");
    let base = surface(&m, "BASE");
    let base_env = surface(&env, "BASE");

    // (1) The mapped BASE is BIT-unchanged by the internal horizon's presence: the
    //     split stack's BASE equals the envelope-only stack's BASE exactly. A derived
    //     surface never displaces a mapped one.
    assert_eq!(
        base, base_env,
        "mapped BASE must be bit-identical with and without the internal split"
    );

    // (2) Where the envelope MERGES: both sub-zones collapse to EXACTLY zero — the
    //     internal horizon coincides with both mapped bounds (no phantom thickness).
    let mut saw_merged = false;
    for jp in 0..N {
        for ip in 0..N {
            if !merged(ip, jp) {
                continue;
            }
            saw_merged = true;
            let n = jp * N + ip;
            assert_eq!(base[n], top[n], "merged node ({ip},{jp}): BASE != TOP");
            assert_eq!(
                mid[n], top[n],
                "merged node ({ip},{jp}): internal split not collapsed onto TOP (phantom)"
            );
            assert_eq!(mid[n], base[n], "merged node ({ip},{jp}): split != BASE");
        }
    }
    assert!(saw_merged, "fixture must exercise a merged region");

    // (3) Where the envelope is OPEN: the split sits strictly inside [TOP, BASE] at
    //     the pick-derived thickness, so both sub-zones are genuinely positive.
    let mut saw_open = false;
    for jp in 0..N {
        for ip in 0..N {
            if merged(ip, jp) {
                continue;
            }
            saw_open = true;
            let n = jp * N + ip;
            assert!(
                base[n] - top[n] > GROSS - 1e-9,
                "open node ({ip},{jp}): envelope should be ~{GROSS} m"
            );
            assert!(
                mid[n] > top[n] + 1e-6,
                "open node ({ip},{jp}): upper sub-zone zero"
            );
            assert!(
                mid[n] < base[n] - 1e-6,
                "open node ({ip},{jp}): lower sub-zone zero"
            );
        }
    }
    assert!(saw_open, "fixture must exercise an open region");

    // (4) At the picks the split lands on the measured pick thickness exactly.
    let tf = top_field();
    for &(ip, jp) in pick_nodes {
        let n = jp * N + ip;
        assert!(
            (mid[n] - (tf[n] + PICK_THICK)).abs() < 1e-6,
            "pick ({ip},{jp}): split {} not at pick depth {}",
            mid[n],
            tf[n] + PICK_THICK
        );
    }

    // (5) NO interface repair fires on these consistent inputs — mapped-vs-mapped
    //     never crosses by construction, and the derived split is clamped inside the
    //     envelope, so the order-repair is a no-op. (This is the phantom's old source.)
    let prov = m.provenance();
    let stack_prov = prov.stack.as_ref().unwrap();
    assert!(
        stack_prov.interface_repairs.is_empty(),
        "expected zero interface repairs, got {:?}",
        stack_prov.interface_repairs
    );
    assert!(
        !prov
            .warnings
            .iter()
            .any(|w| matches!(w, BuildWarning::ThinColumnsRepaired { .. })),
        "no thin-column repair should be needed"
    );

    // (6) GRV conservation: the two sub-zones partition the envelope exactly — their
    //     GRV sums to the single-zone envelope GRV (the merged region contributes
    //     ZERO, not a phantom slab).
    let zoned = m.in_place_by_zone().unwrap();
    let split_grv: f64 = zoned.zones.iter().map(|z| z.in_place.grv_m3).sum();
    let env_grv = env.in_place_by_zone().unwrap().total.grv_m3;
    let rel = (split_grv - env_grv).abs() / env_grv.max(1.0);
    assert!(
        rel < 1e-9,
        "zone GRV {split_grv} != envelope GRV {env_grv} (rel {rel:.2e})"
    );

    // Each sub-zone's GRV is a genuine, positive share of the envelope (the split is
    // real in the open half, not degenerate).
    for z in &zoned.zones {
        assert!(
            z.in_place.grv_m3 > 0.0,
            "sub-zone {} GRV must be positive",
            z.zone
        );
    }
}

#[test]
fn merged_envelope_collapses_both_subzones_no_phantom_axis_local() {
    run_case(merged_axis, &[(6, 2), (7, 5), (6, 6)], false);
}

#[test]
fn merged_envelope_collapses_both_subzones_no_phantom_diag_world_georef() {
    // Doctrine R1: the world-georeferenced, azimuthed (diagonal merge line) variant —
    // the fix must not be an axis-aligned accident, and the world frame must not
    // change the geometry (GRV is the local area-scaled volume regardless).
    run_case(merged_diag, &[(6, 3), (5, 5), (7, 6)], true);
}

#[test]
fn subzone_thickness_is_exactly_zero_over_the_merged_cells() {
    // Cell-level statement of the phantom's absence: every CELL whose four areal
    // corners are all merged has exactly zero thickness in BOTH sub-zones — no slab
    // of phantom volume anywhere in the merged region.
    let m = build(split_stack(merged_axis, &[(6, 2), (7, 5)]), false);
    let dims = m.grid().dims();
    let merged_cell = |i: usize, j: usize| {
        merged_axis(i, j)
            && merged_axis(i + 1, j)
            && merged_axis(i, j + 1)
            && merged_axis(i + 1, j + 1)
    };
    let mut checked = 0;
    for k in 0..dims.nk {
        for j in 0..dims.nj {
            for i in 0..dims.ni {
                if merged_cell(i, j) {
                    let dz = m.grid().cell(Ijk::new(i, j, k)).dz();
                    assert!(
                        dz.abs() < 1e-9,
                        "merged cell ({i},{j},{k}) dz {dz} (phantom)"
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 0, "must exercise fully-merged cells");
}
