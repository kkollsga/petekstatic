//! Integration tests against the REAL `petekio = "0.2.1"` crate — petekio's GATE-0
//! verification (their "integration-tests-required" note). Each test exercises one
//! requirement of the locked `ModelInputs` contract end-to-end, from a real `GeoData`
//! built off shipped fixtures through srs-data's adapter + `data_to_wireframe`.
//!
//! Fixtures (copied from petekio's repo) live in `tests/fixtures/`.
//!
//! Scope: the standalone petekio-contract checks that srs-data can verify in
//! isolation. The sampler-path + refine-loop end-to-end tests (which need petekSim's
//! `srs-uncertainty` sampler and `srs-core` refine loop) live downstream in petekSim
//! (`srs-core/tests/integration_data.rs`) — srs-data no longer reaches up to them.
//!
//! Findings reported back to petekio are noted inline as `GAP:`.

use srs_data::adapter::ModelScalars;
use srs_data::petekio::{GeoData, GridGeometry, Surface, Unit};
use srs_petro::{upscale_porosity, WeightedSample};

const SURFACE: &str = "tests/fixtures/simple.irap";
const WELL_DIR: &str = "tests/fixtures/wells/15_9-A1";

/// A real GeoData loaded from the shipped fixtures (a surface + a well).
fn loaded_geo() -> GeoData {
    let mut geo = GeoData::new(Unit::Feet);
    geo.load_surface("top", SURFACE).expect("load_surface");
    geo.load_well("15/9-A1", (1200.0, 1500.0), 0.0, WELL_DIR)
        .expect("load_well");
    geo
}

// 1. Compiles + runs against the real crate.
#[test]
fn builds_real_model_inputs() {
    let geo = loaded_geo();
    let mi = geo.model_inputs().expect("model_inputs() must succeed");
    // it produced something usable
    assert!(mi.summary.area_m2.value.is_finite());
    assert!(
        !mi.spatial.horizons.is_empty(),
        "the loaded surface yields a horizon"
    );
}

// 2. Every summary field consumed, incl. the Option contacts == None.
#[test]
fn every_summary_field_consumed() {
    let mi = loaded_geo().model_inputs().unwrap();
    let s = &mi.summary;
    for (name, v) in [
        ("area", s.area_m2.value),
        ("net_pay", s.net_pay_m.value),
        ("porosity", s.porosity_frac.value),
        ("sw", s.water_saturation_frac.value),
        ("ntg", s.net_to_gross_frac.value),
    ] {
        assert!(v.is_finite(), "{name} should be finite, got {v}");
    }
    // Contacts are Option and always None today — our code must handle None.
    assert!(
        s.owc_depth_m.is_none(),
        "owc_depth_m is None in the current fixtures"
    );
    assert!(
        s.goc_depth_m.is_none(),
        "goc_depth_m is None in the current fixtures"
    );
    // The adapter consumes the whole summary (incl. the None options) without panic.
    let scalars = ModelScalars::from_summary(s);
    assert!(scalars.owc_depth_m.is_none());
    assert!(scalars.goc_depth_m.is_none());
}

// 3. Every spatial field consumed — horizons, well_curves (incl. xyz), boundary Option.
#[test]
fn every_spatial_field_consumed() {
    let mi = loaded_geo().model_inputs().unwrap();
    let sp = &mi.spatial;

    // horizons: each surface is readable via geom + values()
    for h in &sp.horizons {
        assert!(h.surface.geom.ncol >= 1 && h.surface.geom.nrow >= 1);
        assert_eq!(
            h.surface.values().len(),
            h.surface.geom.ncol * h.surface.geom.nrow
        );
    }
    // well_curves: md/values/xyz aligned 1:1 (the positioning field is present)
    assert!(!sp.well_curves.is_empty(), "the loaded well yields curves");
    for c in &sp.well_curves {
        assert_eq!(c.md.len(), c.values.len(), "{} md/values", c.mnemonic);
        assert_eq!(c.md.len(), c.xyz.len(), "{} md/xyz", c.mnemonic);
    }
    // boundary is Option<PolygonSet>; handle both arms. model_inputs() derives it from
    // the surface's convex hull, so it is Some here.
    match &sp.boundary {
        Some(p) => {
            let bb = p.bbox();
            assert!(bb.xmax >= bb.xmin && bb.ymax >= bb.ymin);
            // GAP (q_petekio_polygon_rings): PolygonSet exposes only bbox/area/contains/
            // clip — no ring-vertex accessor — so we cannot read the true outline.
        }
        None => { /* also valid; assemble_wireframe falls back to the area square */ }
    }
}

// 5. Positioning seam — upscale a curve onto cells via WellCurveInput.xyz, incl. [NaN;3].
#[test]
fn positioning_seam_upscales_from_xyz() {
    let mi = loaded_geo().model_inputs().unwrap();
    let curve = mi
        .spatial
        .well_curves
        .iter()
        .find(|c| c.values.iter().any(|v| (0.0..=1.0).contains(v)))
        .expect("a fraction-valued curve (e.g. PHIE/SW)");

    // build positioned samples, skipping unpositioned [NaN;3] and NaN values
    let samples: Vec<WeightedSample> = curve
        .xyz
        .iter()
        .zip(&curve.values)
        .filter(|([x, y, z], v)| {
            !x.is_nan() && !y.is_nan() && !z.is_nan() && !v.is_nan() && (0.0..=1.0).contains(*v)
        })
        .map(|(_, &v)| WeightedSample::new(1.0, v)) // unit thickness per sample
        .collect();
    assert!(
        !samples.is_empty(),
        "at least one positioned in-range sample"
    );
    let upscaled = upscale_porosity(&samples).expect("upscale onto a cell");
    assert!(
        (0.0..=1.0).contains(&upscaled),
        "upscaled fraction {upscaled}"
    );
}

// 6. Lattice hand-off — Surface::resample onto our (denser) GridGeometry.
// petekio ≥0.2.8 (io-centralize): `resample` is FALLIBLE — a rotated geometry is
// a typed `GeoError::Unsupported` (until suite task_suite_grid_rotation), never a
// silent wrong answer. The fixture surface is rotated 30°, so it pins the loud
// arm; the supported axis-aligned hand-off is exercised on a derotated copy.
#[test]
fn lattice_handoff_resamples_onto_our_geometry() {
    let geo = loaded_geo();
    let surf = geo.surface("top").expect("surface present");
    // our target lattice: same extent, doubled resolution
    let target_for = |g: &GridGeometry| GridGeometry {
        xori: g.xori,
        yori: g.yori,
        xinc: g.xinc / 2.0,
        yinc: g.yinc / 2.0,
        ncol: g.ncol * 2,
        nrow: g.nrow * 2,
        rotation_deg: 0.0,
        yflip: false,
    };

    // Rotated source (the fixture is rotated): a LOUD typed error, not a wrong grid.
    let Err(err) = surf.resample(&target_for(&surf.geom)) else {
        panic!("rotated resample must be a typed error (Surface is not Debug)");
    };
    assert!(
        err.to_string().contains("rotated"),
        "rotated resample must error loudly: {err}"
    );

    // Axis-aligned source: the supported hand-off resamples onto our lattice.
    let mut g0 = surf.geom.clone();
    g0.rotation_deg = 0.0;
    let axis = Surface::new(g0.clone(), surf.values().clone()).expect("axis-aligned copy");
    let target = target_for(&g0);
    let resampled = axis
        .resample(&target)
        .expect("axis-aligned resample is supported");
    assert_eq!(resampled.geom.ncol, target.ncol);
    assert_eq!(resampled.geom.nrow, target.nrow);
    let defined = resampled.values().iter().filter(|v| !v.is_nan()).count();
    assert!(defined > 0, "resampled surface has defined nodes");
}
