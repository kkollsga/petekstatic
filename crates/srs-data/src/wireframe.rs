//! `data_to_wireframe` assembly (spec: `data_to_wireframe_spec`, confidence medium).
//!
//! Assemble the constraining [`Wireframe`] `{boundary, horizons, contacts}` from
//! petekio's [`ModelInputs`] per the spec. Still a thin consumer: we *classify* and
//! *arrange* constraints (the wireframe's own job — its types are the output) but do
//! no input-data processing. Faults are deferred (spec MVP boundary).
//!
//! ## Real boundary footprints (petekio `PolygonSet::rings()`, ≥0.2.2)
//! petekio's `PolygonSet` exposes `rings()` — the exterior ring vertices of each
//! polygon (`[x, y, z]`, closed, Z dropped) — so we read the **true** field
//! outline into [`Boundary::ring`]. The bounding-box rectangle is now only a
//! **documented fallback** for a ringless set (all rings dropped as degenerate).
//! Resolves `q_petekio_polygon_rings` (`task_petekstatic_boundary_rings`).
//!
//! ## z-datum: one positive-down convention inside the [`Wireframe`]
//! petekio delivers surfaces as **negative-down subsea elevation** and contacts
//! as **positive-down `depth_m`**. The whole GEOMODEL layer (srs-model's
//! `BuildOpts`/gridder, `Contact.depth_m`) works in positive-down `depth_m`, so
//! surfaces are **negated at ingest** ([`surface_depths`]) to land on that same
//! datum — horizons and contacts then share one convention, and structural role
//! assignment (`Top` = shallowest) reads correctly. Contacts pass straight
//! through (already positive-down metres).
//!
//! ## Stub-stage simplification (resolve at integration — G2)
//! A petekio [`Surface`] is copied node-for-node (`surface.values()` over
//! `surface.geom`) into a [`GriddedDepth`]; resampling onto *our* model lattice
//! (`Surface::resample(&GridGeometry)`) lands when our 3D `GridGeometry` exists.

use crate::adapter::hardness_of;
use crate::petekio::{HorizonInput, ModelInputs, PolygonSet, Surface};
use petekstatic_error::StaticError;
use srs_wireframe::{
    Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
};

/// Mean of the defined (non-NaN) depths, or `NaN` if none are defined.
fn mean_defined(depths: &[f64]) -> f64 {
    let (sum, n) = depths
        .iter()
        .filter(|d| !d.is_nan())
        .fold((0.0, 0usize), |(s, n), d| (s + d, n + 1));
    if n == 0 {
        f64::NAN
    } else {
        sum / n as f64
    }
}

/// Areal boundary. With a supplied polygon we read its **real exterior ring**
/// (petekio `PolygonSet::rings()`) into the boundary outline (`Interpolated`),
/// falling back to the bounding-box rectangle only for a ringless set. With no
/// polygon we seed a square of side `sqrt(area_m2)` (`Assumed`).
fn boundary_from(boundary: Option<&PolygonSet>, area_m2: f64) -> Result<Boundary, StaticError> {
    match boundary {
        Some(p) => {
            // The true field outline: the first usable exterior ring (closed,
            // ≥3 vertices), projected to 2-D. Non-rectangular footprints are
            // preserved instead of being flattened to their bbox.
            if let Some(exterior) = p.rings().into_iter().find(|r| r.len() >= 3) {
                return Ok(Boundary {
                    ring: exterior.iter().map(|c| [c[0], c[1]]).collect(),
                    hardness: Hardness::Interpolated,
                });
            }
            // Documented fallback: a ringless PolygonSet (all rings degenerate) has
            // no outline to read — use its bounding-box rectangle.
            let bb = p.bbox();
            Ok(Boundary {
                ring: vec![
                    [bb.xmin, bb.ymin],
                    [bb.xmax, bb.ymin],
                    [bb.xmax, bb.ymax],
                    [bb.xmin, bb.ymax],
                    [bb.xmin, bb.ymin],
                ],
                hardness: Hardness::Interpolated,
            })
        }
        None => {
            if !(area_m2.is_finite() && area_m2 > 0.0) {
                return Err(StaticError::InvalidInput(format!(
                    "no boundary polygon and area not positive ({area_m2} m²)"
                )));
            }
            // petekio ≥0.3.0 delivers SI area (m²) — no conversion; the metric
            // boundary ring is a square of the given footprint.
            let side = area_m2.sqrt();
            Ok(Boundary {
                ring: vec![
                    [0.0, 0.0],
                    [side, 0.0],
                    [side, side],
                    [0.0, side],
                    [0.0, 0.0],
                ],
                hardness: Hardness::Assumed,
            })
        }
    }
}

/// petekio [`Surface`] depths as **positive-down depth in metres** (row-major over
/// `geom`), read via the public `values()`.
///
/// ## z-datum flip (the one Wireframe convention)
/// petekio delivers surface values as **negative-down subsea elevation**
/// (`xyz()` = subsea elevation, ≥0.3.0); the whole GEOMODEL layer works in
/// **positive-down `depth_m`** — the datum srs-model's `BuildOpts`/gridder and
/// our `Contact.depth_m` already use. So we **negate at this ingest boundary**,
/// unifying surfaces and contacts onto one positive-down convention inside the
/// [`Wireframe`]. This retires the deferral once tracked in
/// `q_petekio_modelinputs_si` (its blocker — imperial `SummaryInputs` contacts —
/// resolved when petekio went SI, 2026-07-04). `NaN` (undefined) negates to
/// `NaN`, so the defined-node mask is unaffected.
fn surface_depths(surface: &Surface) -> Vec<f64> {
    surface.values().iter().map(|z| -z).collect()
}

/// petekio [`Surface`] → wireframe [`GriddedDepth`] (node-for-node; resample deferred).
///
/// Surface node values are flipped from petekio's negative-down elevation to our
/// positive-down `depth_m` at the [`surface_depths`] ingest boundary (petekio is
/// metric-native, so no unit conversion — only the sign). See [`surface_depths`]
/// for the datum rationale.
fn gridded_depth_of(surface: &Surface, hardness: Hardness) -> GriddedDepth {
    let is_hard = hardness == Hardness::Hard;
    let depth_m = surface_depths(surface);
    let is_control = depth_m.iter().map(|d| is_hard && !d.is_nan()).collect();
    GriddedDepth {
        ncol: surface.geom.ncol,
        nrow: surface.geom.nrow,
        depth_m,
        is_control,
    }
}

/// Assign structural roles by depth: with multiple horizons the shallowest is
/// `Top`, the deepest is `Base`, the rest `Intermediate`; a single horizon is the
/// `Top`. Depth ordering is a geometric/structural determination (the gridder needs
/// stacking order), not input processing.
///
/// `mean_depths` are **positive-down `depth_m`** (flipped from petekio's
/// negative-down elevation at [`surface_depths`]), so the numerically-smallest
/// depth is the structurally shallowest = `Top`. Feeding raw negative-down
/// elevations here would invert the pick (deepest labelled `Top`).
fn role_for(idx: usize, mean_depths: &[f64]) -> HorizonRole {
    if mean_depths.len() == 1 {
        return HorizonRole::Top;
    }
    let shallowest = mean_depths
        .iter()
        .enumerate()
        .filter(|(_, d)| !d.is_nan())
        .min_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i);
    let deepest = mean_depths
        .iter()
        .enumerate()
        .filter(|(_, d)| !d.is_nan())
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i);
    if Some(idx) == shallowest {
        HorizonRole::Top
    } else if Some(idx) == deepest {
        HorizonRole::Base
    } else {
        HorizonRole::Intermediate
    }
}

fn horizons_from(inputs: &[HorizonInput]) -> Vec<Horizon> {
    let mean_depths: Vec<f64> = inputs
        .iter()
        .map(|h| mean_defined(&surface_depths(&h.surface)))
        .collect();
    inputs
        .iter()
        .enumerate()
        .map(|(i, h)| {
            let hardness = hardness_of(h.provenance);
            Horizon {
                name: h.name.clone(),
                role: role_for(i, &mean_depths),
                surface: gridded_depth_of(&h.surface, hardness),
            }
        })
        .collect()
}

/// Assemble the constraining [`Wireframe`] from petekio's model-ready inputs.
///
/// # Errors
/// [`StaticError::InvalidInput`] if no areal extent can be derived (neither a boundary
/// polygon nor a positive `area_m2`).
pub fn assemble_wireframe(inputs: &ModelInputs) -> Result<Wireframe, StaticError> {
    let boundary = boundary_from(
        inputs.spatial.boundary.as_ref(),
        inputs.summary.area_m2.value,
    )?;
    let horizons = horizons_from(&inputs.spatial.horizons);

    // petekio ≥0.3.0 delivers contacts as positive-down depth in METRES
    // (`owc_depth_m`/`goc_depth_m`) — already the Wireframe's positive-down
    // `depth_m` datum (the one surfaces are flipped onto at `surface_depths`), so
    // they pass straight through (the ft→m seam shim retired 2026-07-04).
    let mut contacts = Vec::new();
    if let Some(owc) = inputs.summary.owc_depth_m {
        contacts.push(Contact {
            kind: ContactKind::Owc,
            depth_m: owc.value,
            hardness: hardness_of(owc.provenance),
        });
    }
    if let Some(goc) = inputs.summary.goc_depth_m {
        contacts.push(Contact {
            kind: ContactKind::Goc,
            depth_m: goc.value,
            hardness: hardness_of(goc.provenance),
        });
    }

    Ok(Wireframe {
        boundary,
        horizons: std::sync::Arc::new(horizons),
        contacts,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::petekio::{Distribution, Provenance, SpatialInputs, SummaryInputs, Uncertain};

    /// Build a ModelInputs with no spatial geometry (boundary None, no horizons/curves)
    /// — exercises the scalar-driven paths (boundary fallback, contacts, errors) without
    /// the unconstructable real `Surface`/`PolygonSet`. The horizon + polygon paths are
    /// covered by the real-petekio integration tests.
    fn scalar_inputs(area_m2: f64, owc_depth_m: Option<f64>) -> ModelInputs {
        let det = |value, provenance| Uncertain {
            value,
            distribution: Distribution::Deterministic,
            provenance,
        };
        ModelInputs {
            summary: SummaryInputs {
                area_m2: det(area_m2, Provenance::Assumed),
                net_pay_m: det(25.0, Provenance::HardData),
                porosity_frac: det(0.22, Provenance::HardData),
                water_saturation_frac: det(0.30, Provenance::HardData),
                net_to_gross_frac: det(0.80, Provenance::Interpolated),
                owc_depth_m: owc_depth_m.map(|d| det(d, Provenance::HardData)),
                goc_depth_m: None,
            },
            spatial: SpatialInputs {
                boundary: None,
                horizons: vec![],
                well_curves: vec![],
            },
        }
    }

    #[test]
    fn boundary_falls_back_to_area_square_when_no_polygon() {
        let wf = assemble_wireframe(&scalar_inputs(2_509_000.0, None)).unwrap();
        assert_eq!(wf.boundary.hardness, Hardness::Assumed);
        // Metric side: sqrt of the SI area, straight through (no conversion).
        let side = 2_509_000.0_f64.sqrt();
        assert!((wf.boundary.ring[2][0] - side).abs() < 1e-6);
        assert!(wf.horizons.is_empty());
    }

    #[test]
    fn contacts_from_owc_with_provenance_hardness() {
        let wf = assemble_wireframe(&scalar_inputs(2_509_000.0, Some(2511.6))).unwrap();
        assert_eq!(wf.contacts.len(), 1);
        assert_eq!(wf.contacts[0].kind, ContactKind::Owc);
        // petekio delivers owc_depth_m (metres) — stored straight through.
        assert!((wf.contacts[0].depth_m - 2511.6).abs() < 1e-9);
        assert_eq!(wf.contacts[0].hardness, Hardness::Hard);
    }

    #[test]
    fn no_contact_when_owc_goc_absent() {
        let wf = assemble_wireframe(&scalar_inputs(2_509_000.0, None)).unwrap();
        assert!(wf.contacts.is_empty());
    }

    #[test]
    fn no_extent_at_all_is_an_error() {
        assert!(assemble_wireframe(&scalar_inputs(0.0, None)).is_err());
    }

    /// A single-node-uniform horizon surface at `elevation_m` (negative-down
    /// subsea elevation, petekio's convention). A 2×2 constant grid is enough to
    /// drive role assignment; `Surface::constant` avoids an ndarray dep here.
    fn horizon_at(name: &str, elevation_m: f64) -> HorizonInput {
        use crate::petekio::GridGeometry;
        let geom = GridGeometry {
            xori: 0.0,
            yori: 0.0,
            xinc: 100.0,
            yinc: 100.0,
            ncol: 2,
            nrow: 2,
            rotation_deg: 0.0,
            yflip: false,
        };
        HorizonInput {
            name: name.to_string(),
            surface: Surface::constant(geom, elevation_m),
            provenance: Provenance::Interpolated,
        }
    }

    #[test]
    fn top_is_structurally_shallowest_under_negative_down_elevation() {
        // petekio delivers surfaces as NEGATIVE-DOWN subsea elevation: the
        // structurally SHALLOWER surface carries the LARGER (less negative) z,
        // the deeper one the more negative z. So the shallow surface (2000 m
        // deep, elevation −2000) is Top; the deep one (2100 m, elevation −2100)
        // is Base. Under the old positive-down assumption `role_for` picked Top
        // by numeric MIN — which is the MOST-negative (deepest) elevation — and
        // labelled the deepest surface Top. This asserts the corrected ordering.
        let horizons = horizons_from(&[
            horizon_at("ShallowTop", -2000.0),
            horizon_at("DeepBase", -2100.0),
        ]);
        assert_eq!(
            horizons[0].role,
            HorizonRole::Top,
            "the structurally shallowest surface must be Top"
        );
        assert_eq!(
            horizons[1].role,
            HorizonRole::Base,
            "the structurally deepest surface must be Base"
        );
    }

    #[test]
    fn roles_order_by_depth() {
        // single horizon -> Top
        assert_eq!(role_for(0, &[8000.0]), HorizonRole::Top);
        // multiple -> shallowest Top, deepest Base, middle Intermediate
        let depths = [8200.0, 8000.0, 8100.0];
        assert_eq!(role_for(1, &depths), HorizonRole::Top); // 8000 shallowest
        assert_eq!(role_for(0, &depths), HorizonRole::Base); // 8200 deepest
        assert_eq!(role_for(2, &depths), HorizonRole::Intermediate);
    }

    #[test]
    fn mean_defined_skips_nan() {
        assert_eq!(mean_defined(&[8000.0, f64::NAN, 8200.0]), 8100.0);
        assert!(mean_defined(&[f64::NAN, f64::NAN]).is_nan());
    }

    /// Shoelace area of a closed 2-D ring.
    fn shoelace(ring: &[[f64; 2]]) -> f64 {
        let mut s = 0.0;
        for w in ring.windows(2) {
            s += w[0][0] * w[1][1] - w[1][0] * w[0][1];
        }
        (s / 2.0).abs()
    }

    #[test]
    fn boundary_reads_the_real_ring_not_the_bbox() {
        use crate::petekio::PolygonSet;
        use std::io::Write;
        // A synthetic non-rectangular outline: a right triangle (0,0)-(100,0)-
        // (0,100). Its true area is 5000; its bbox is 100x100 = 10000. The
        // boundary must carry the triangle, not the bounding rectangle.
        let path = std::env::temp_dir().join(format!("srs_data_ring_{}.irap", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "0 0 0\n100 0 0\n0 100 0").unwrap();
        }
        let poly = PolygonSet::load_irap_polygons(&path).unwrap();
        let b = boundary_from(Some(&poly), 999.0).unwrap();
        let _ = std::fs::remove_file(&path);

        let area = shoelace(&b.ring);
        // Matches the true triangle (== petekio's own polygon area), not the bbox.
        assert!(
            (area - 5000.0).abs() < 1.0,
            "ring area {area} should be the triangle (5000), not the bbox (10000)"
        );
        assert!(
            (area - poly.area()).abs() < 1.0,
            "ring area tracks polygon area"
        );
        assert!(
            area < 9000.0,
            "ring must not collapse to the 10000 bbox rectangle"
        );
        assert_eq!(b.hardness, Hardness::Interpolated);
    }
}
