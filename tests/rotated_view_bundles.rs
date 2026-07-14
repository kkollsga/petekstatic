//! Rotated world-frame acceptance for map/section bundles (suite Phase 7c).
//!
//! **Mode matrix:** this file owns the in-core, serial, two-zone horizon-stack
//! rotated-frame cell. Existing `mode_matrix.rs` retains the in-core/spilled and
//! single-zone parity cells; rotation is a world-seam transform over those same
//! local kernels.

use petekstatic::gridder::{Conformity, SolveOpts};
use petekstatic::model::{
    BuildOpts, ConstantPriors, Georef, HorizonSource, HorizonStack, MapBundle, MapSpec,
    PropertyPipeline, SectionSpec, StackFrame, StackHorizon, StackZone, StaticModel,
    StaticModelBuilder, TrendSurface, UpscaleMethod, WellLog, WellTie, WorldPoint,
};
use petekstatic::wireframe::GriddedDepth;
use petektools::Lattice;

const N: usize = 5;
const NI: usize = N - 1;
const SPACING: f64 = 25.0;
const ORIGIN_X: f64 = 510_000.0;
const ORIGIN_Y: f64 = 6_710_000.0;
const ROTATION_DEG: f64 = 30.0;

fn surface(base: f64, dip_i: f64, dip_j: f64) -> GriddedDepth {
    let mut depth_m = Vec::with_capacity(N * N);
    for j in 0..N {
        for i in 0..N {
            depth_m.push(base + dip_i * i as f64 + dip_j * j as f64);
        }
    }
    GriddedDepth {
        ncol: N,
        nrow: N,
        depth_m,
        is_control: vec![true; N * N],
    }
}

fn model() -> StaticModel {
    let stack = HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Mapped(surface(2_000.0, 2.0, 1.0)),
            },
            StackHorizon {
                name: "MID".into(),
                source: HorizonSource::Mapped(surface(2_025.0, 2.0, 1.0)),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(surface(2_055.0, 2.0, 1.0)),
            },
        ],
        zone_layers: vec![
            StackZone::new("UPPER", Conformity::Proportional, 2, Vec::new()).with_color("#4477aa"),
            StackZone::new("LOWER", Conformity::Proportional, 2, Vec::new()).with_color("#cc6677"),
        ],
    };
    let georef =
        Georef::oriented(ORIGIN_X, ORIGIN_Y, SPACING, SPACING, ROTATION_DEG, false).unwrap();
    let well_xy = georef.intrinsic_to_world(1.0, 2.0);
    let opts = BuildOpts {
        area_m2: (SPACING * NI as f64).powi(2),
        gross_height_m: 55.0,
        nk: 4,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.22,
            net_to_gross: 0.71,
            water_saturation: 0.31,
        },
    };
    let trend = TrendSurface::new(
        NI,
        NI,
        (0..NI)
            .flat_map(|j| (0..NI).map(move |i| 0.75 + 0.10 * i as f64 + 0.02 * j as f64))
            .collect(),
    )
    .unwrap()
    .with_oriented_georef(ORIGIN_X, ORIGIN_Y, SPACING, SPACING, ROTATION_DEG, false);
    StaticModelBuilder::from_horizon_stack(stack, opts)
        .unwrap()
        .with_oriented_georef(ORIGIN_X, ORIGIN_Y, SPACING, SPACING, ROTATION_DEG, false)
        .with_areal_trend(trend)
        .with_property(PropertyPipeline::new("ROT_LOG").upscale(
            vec![WellLog::new(
                well_xy.0,
                well_xy.1,
                vec![(2_010.0, 0.10), (2_045.0, 0.80)],
            )],
            UpscaleMethod::Arithmetic,
        ))
        .with_well_ties(vec![WellTie::new("SYNTH-WELL", well_xy.0, well_xy.1, 1, 2)])
        .build()
        .unwrap()
}

fn lattice(bundle: &MapBundle) -> Lattice {
    let f = bundle.frame;
    Lattice {
        xori: f.origin_x,
        yori: f.origin_y,
        xinc: f.spacing_x,
        yinc: f.spacing_y,
        ncol: f.ncol,
        nrow: f.nrow,
        rotation_deg: f.rotation_deg,
        yflip: f.yflip,
    }
}

#[test]
fn rotated_map_section_and_cursor_share_one_world_frame() {
    let model = model();
    let rotated_log = model.property("ROT_LOG").unwrap();
    let informed: Vec<f64> = rotated_log
        .values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect();
    assert_eq!(informed.len(), 2, "the rotated world well hits one column");
    assert!(informed.iter().any(|value| (*value - 0.10).abs() < 1e-12));
    assert!(informed.iter().any(|value| (*value - 0.80).abs() < 1e-12));
    let map = model
        .map_bundle(&MapSpec::new().property("PORO").property("NTG").k_slice(0))
        .unwrap();
    assert_eq!(map.schema_version, 6);
    assert_eq!(map.frame.rotation_deg, ROTATION_DEG);
    assert!(!map.frame.yflip);
    assert_eq!(map.zone_averages.len(), 4);
    assert_eq!(map.k_slices.len(), 2);
    let ntg = map
        .zone_averages
        .iter()
        .find(|layer| layer.name == "NTG::UPPER")
        .unwrap();
    assert_ne!(ntg.values[0].to_bits(), ntg.values[NI - 1].to_bits());

    // The derived stack outline is the oriented half-cell footprint, not its
    // axis-aligned bounding box.
    let expected = [(-0.5, -0.5), (3.5, -0.5), (3.5, 3.5), (-0.5, 3.5)];
    for (point, (fi, fj)) in map.outline[0].iter().take(4).zip(expected) {
        let xy = map.frame.intrinsic_to_world(fi, fj);
        assert!((point[0] - xy.0).abs() < 1e-8);
        assert!((point[1] - xy.1).abs() < 1e-8);
    }

    // Well marker, map cursor and section march all invert through the exact
    // petekTools lattice convention.
    let lat = lattice(&map);
    let well = &map.wells[0];
    let well_ij = lat.xy_to_ij(well.x, well.y).unwrap();
    assert!((well_ij.0 - 1.0).abs() < 1e-9);
    assert!((well_ij.1 - 2.0).abs() < 1e-9);

    let a = map.frame.intrinsic_to_world(-0.25, 1.0);
    let b = map.frame.intrinsic_to_world(3.25, 1.0);
    let section = model
        .intersection_bundle(
            &SectionSpec::Polyline(vec![[a.0, a.1], [b.0, b.1]]),
            Some("PORO"),
        )
        .unwrap();
    assert_eq!(section.frame, Some(map.frame));
    assert!(!section.columns.is_empty());
    assert_eq!(section.zones.len(), 2);
    for column in &section.columns {
        let ij = lat.xy_to_ij(column.x, column.y).unwrap();
        assert!((ij.0 - column.i as f64).abs() <= 0.5 + 1e-9);
        assert!((ij.1 - column.j as f64).abs() <= 0.5 + 1e-9);
        assert_eq!(column.values.len(), 4);
        assert!(column.zone_ids.contains(&0));
        assert!(column.zone_ids.contains(&1));
    }

    let map_back: MapBundle = serde_json::from_str(&serde_json::to_string(&map).unwrap()).unwrap();
    assert_eq!(map_back, map);
    let section_back = serde_json::from_str::<petekstatic::model::IntersectionBundle>(
        &serde_json::to_string(&section).unwrap(),
    )
    .unwrap();
    assert_eq!(section_back, section);
}

#[test]
fn oriented_georef_normalizes_rotation_and_round_trips_yflip() {
    let georef = Georef::oriented(10.0, 20.0, 2.0, 3.0, 390.0, true).unwrap();
    assert_eq!(georef.rotation_deg, 30.0);
    assert!(georef.yflip);
    let xy = georef.intrinsic_to_world(1.25, 2.5);
    let ij = georef.world_to_intrinsic(xy.0, xy.1).unwrap();
    assert!((ij.0 - 1.25).abs() < 1e-9);
    assert!((ij.1 - 2.5).abs() < 1e-9);
    let back: Georef = serde_json::from_str(&serde_json::to_string(&georef).unwrap()).unwrap();
    assert_eq!(back, georef);
    assert!(Georef::oriented(0.0, 0.0, 1.0, 1.0, f64::NAN, false).is_none());
}

#[test]
fn rotated_scatter_conditions_on_the_intrinsic_node_lattice() {
    let georef =
        Georef::oriented(ORIGIN_X, ORIGIN_Y, SPACING, SPACING, ROTATION_DEG, true).unwrap();
    let mut points = Vec::with_capacity(N * N);
    for jp in 0..N {
        for ip in 0..N {
            let xy = georef.intrinsic_to_world(ip as f64 - 0.5, jp as f64 - 0.5);
            points.push(WorldPoint {
                x: xy.0,
                y: xy.1,
                depth_m: 1_900.0 + 2.0 * ip as f64 + jp as f64,
            });
        }
    }
    let stack = HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "TOP".into(),
                source: HorizonSource::Scatter(points),
            },
            StackHorizon {
                name: "BASE".into(),
                source: HorizonSource::Mapped(surface(1_960.0, 2.0, 1.0)),
            },
        ],
        zone_layers: vec![StackZone::new(
            "RESERVOIR",
            Conformity::Proportional,
            3,
            Vec::new(),
        )],
    };
    let opts = BuildOpts {
        area_m2: (SPACING * NI as f64).powi(2),
        gross_height_m: 60.0,
        nk: 3,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.2,
            net_to_gross: 0.7,
            water_saturation: 0.3,
        },
    };
    let model = StaticModelBuilder::from_scatter_stack(
        stack,
        opts,
        StackFrame {
            ni: NI,
            nj: NI,
            georef,
        },
    )
    .unwrap()
    .build()
    .unwrap();
    let map = model.map_bundle(&MapSpec::new()).unwrap();
    assert_eq!(map.frame.rotation_deg, ROTATION_DEG);
    assert!(map.frame.yflip);
    assert!(map.horizons[0].values.iter().all(|value| value.is_finite()));
}
