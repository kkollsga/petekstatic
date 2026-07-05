//! [`MapBundle`] — the areal (plan-view) inspection bundle
//! ([`StaticModel::map_bundle`]).
//!
//! Everything a map view draws, pre-computed on the shared [`GridFrame`]:
//! structural depth surfaces, property maps (a single k-slice and the
//! zone/interval average — the useful default), the areal outline ring(s), well
//! surface markers, and per-contact subcrop masks. The viewer renders; it does
//! not compute.

use super::frame::{GridFrame, ScalarLayer};
use super::SCHEMA_VERSION;
use crate::error::StaticError;
use crate::grid::Ijk;
use crate::model::model::StaticModel;
use crate::wireframe::HorizonRole;
use serde::{Deserialize, Serialize};

/// What to put in a [`MapBundle`]: which property cubes to map (each yields a
/// per-zone average map — the useful default) and, optionally, a single k-slice
/// of each. The structural surfaces, outline, wells and contact masks are always
/// included. Immutable/chainable (house style).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MapSpec {
    properties: Vec<String>,
    k_slice: Option<usize>,
}

impl MapSpec {
    /// An empty spec — structural surfaces, outline, wells and contact masks only
    /// (no property maps).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a property cube to map (a zone-average map per zone, plus a k-slice if
    /// [`MapSpec::k_slice`] is set). The named cube must exist on the model.
    #[must_use]
    pub fn property(mut self, name: impl Into<String>) -> Self {
        self.properties.push(name.into());
        self
    }

    /// Also emit a single-layer k-slice map (at layer `k`) for every requested
    /// property, alongside the zone-average maps.
    #[must_use]
    pub fn k_slice(mut self, k: usize) -> Self {
        self.k_slice = Some(k);
        self
    }
}

/// One horizon's tie residual for a well marker (SCHEMA_VERSION 4): the framework
/// `horizon` and the `residual_m` (measured formation top − untied model surface at
/// the well node, metres, positive = well deeper than the model).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WellTieResidual {
    pub horizon: String,
    pub residual_m: f64,
}

/// A well/bore surface marker: its id, world surface `(x, y)`, a summary
/// `tie_residual_m` (the mean per-horizon residual, `None` when the well carries no
/// ties), and — SCHEMA_VERSION 4 — the full per-horizon `ties`
/// (`task_petekstatic_multizone_2`). `ties` is populated from
/// [`crate::model::Provenance::well_ties`] when the model was built with
/// [`crate::model::StaticModelBuilder::with_well_ties`]; it is empty otherwise.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WellMarker {
    pub id: String,
    pub x: f64,
    pub y: f64,
    pub tie_residual_m: Option<f64>,
    /// Per-horizon tie residuals (SCHEMA_VERSION 4); empty when the well has no ties.
    pub ties: Vec<WellTieResidual>,
}

/// A fluid contact's areal subcrop: the contact `kind` (`"OWC"` / `"GOC"` /
/// `"GWC"`), its `depth_m`, and a row-major `j * ncol + i` mask marking the
/// columns whose vertical extent the contact plane crosses (`true` where
/// `top_depth <= depth_m <= base_depth` for that column — the subcrop band the
/// viewer draws as the contact's areal trace). The simpler, honest form of the
/// "contact contour" (a marked mask, not a polyline extraction).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ContactMask {
    pub kind: String,
    pub depth_m: f64,
    /// Row-major `j * ncol + i`; `true` where the contact crosses the column.
    pub crossing: Vec<bool>,
}

/// The areal inspection bundle. All areal layers share [`MapBundle::frame`];
/// coordinates are the model world frame, depths metres positive-down.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapBundle {
    /// View-bundle schema version ([`SCHEMA_VERSION`]).
    pub schema_version: u32,
    /// Input-bundle identity (from provenance).
    pub inputs_ref: String,
    /// The shared areal georeference.
    pub frame: GridFrame,
    /// Areal outline ring(s), world `[x, y]` (the framework boundary).
    pub outline: Vec<Vec<[f64; 2]>>,
    /// Realized structural depth surfaces (model top + base), grid-georeferenced,
    /// named after the framework horizons. Units metres, positive-down.
    pub horizons: Vec<ScalarLayer>,
    /// Per-property, per-zone interval-average maps (name `"{property}::{zone}"`).
    pub zone_averages: Vec<ScalarLayer>,
    /// Per-property single-k-slice maps (name `"{property}::k{n}"`), empty unless
    /// [`MapSpec::k_slice`] was set.
    pub k_slices: Vec<ScalarLayer>,
    /// Well/bore surface markers with per-horizon tie residuals — populated from
    /// [`crate::model::Provenance::well_ties`] (built via
    /// [`crate::model::StaticModelBuilder::with_well_ties`]); empty when the model carries
    /// no ties.
    pub wells: Vec<WellMarker>,
    /// Per-contact subcrop masks.
    pub contacts: Vec<ContactMask>,
}

impl MapBundle {
    /// Stream this bundle to `w` as JSON (no intermediate `Value` tree) — the same
    /// streaming path the volume bundle's envelope uses (`view::wire`). Areal
    /// bundles are small, so they stay plain JSON.
    pub fn write_json<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()> {
        super::wire::write_json(self, w)
    }
}

impl StaticModel {
    /// Export the areal ([`MapBundle`]) inspection bundle: structural depth
    /// surfaces, the requested property maps (zone-average + optional k-slice),
    /// the outline, well markers, and contact subcrop masks — all on one shared
    /// [`GridFrame`].
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the column lattice is smaller than `2x2`
    /// or not axis-aligned/regular, a requested property is absent, or the
    /// `k_slice` index is out of range.
    pub fn map_bundle(&self, spec: &MapSpec) -> Result<MapBundle, StaticError> {
        // The real grid whether in-core or spilled: a spilled (out-of-core) model's
        // `grid()` is a 1×1×1 placeholder with no cubes, so — like the volume and
        // section bundles — the map export must materialize the backing. Reading the
        // placeholder instead produced a misleading "1x1 lattice" error on every
        // spilled model (the map-bundle sibling of `question_volume_bundle_stack_empty`).
        let grid = self.view_grid()?;
        let frame = GridFrame::of_grid(&grid, self.georef())?;
        let dims = grid.dims();
        let (ni, nj, nk) = (dims.ni, dims.nj, dims.nk);

        // Realized structural surfaces from the grid geometry (georef-exact),
        // named after the framework Top/Base horizons where present.
        let name_for = |role: HorizonRole, fallback: &str| -> String {
            self.framework()
                .horizons
                .iter()
                .find(|h| h.role == role)
                .map_or_else(|| fallback.to_string(), |h| h.name.clone())
        };
        let mut top = vec![f64::NAN; ni * nj];
        let mut base = vec![f64::NAN; ni * nj];
        for j in 0..nj {
            for i in 0..ni {
                top[j * ni + i] = grid.cell(Ijk::new(i, j, 0)).top_depth();
                base[j * ni + i] = grid.cell(Ijk::new(i, j, nk - 1)).bottom_depth();
            }
        }
        let horizons = vec![
            ScalarLayer::new(name_for(HorizonRole::Top, "TOP"), "m", top),
            ScalarLayer::new(name_for(HorizonRole::Base, "BASE"), "m", base),
        ];

        // Property maps: a zone-average per zone (the useful default) + an
        // optional single k-slice, for every requested cube.
        let mut zone_averages = Vec::new();
        let mut k_slices = Vec::new();
        for prop_name in &spec.properties {
            let prop = grid.properties().get(prop_name).ok_or_else(|| {
                StaticError::InvalidInput(format!("map_bundle: no property '{prop_name}'"))
            })?;
            for zone in self.zones().zones() {
                let mut m = vec![f64::NAN; ni * nj];
                for j in 0..nj {
                    for i in 0..ni {
                        let (mut sum, mut n) = (0.0f64, 0usize);
                        for k in zone.k_range.clone() {
                            let v = prop.values[(k * nj + j) * ni + i];
                            if v.is_finite() {
                                sum += v;
                                n += 1;
                            }
                        }
                        m[j * ni + i] = if n > 0 { sum / n as f64 } else { f64::NAN };
                    }
                }
                zone_averages.push(ScalarLayer::new(
                    format!("{prop_name}::{}", zone.name),
                    "fraction",
                    m,
                ));
            }
            if let Some(k) = spec.k_slice {
                if k >= nk {
                    return Err(StaticError::InvalidInput(format!(
                        "map_bundle: k_slice {k} out of range (nk={nk})"
                    )));
                }
                let mut m = vec![f64::NAN; ni * nj];
                for j in 0..nj {
                    for i in 0..ni {
                        m[j * ni + i] = prop.values[(k * nj + j) * ni + i];
                    }
                }
                k_slices.push(ScalarLayer::new(
                    format!("{prop_name}::k{k}"),
                    "fraction",
                    m,
                ));
            }
        }

        // Contact subcrop masks: a column is crossed where the contact plane lies
        // between the column's realized top and base depth.
        let contacts = self
            .contacts()
            .iter()
            .map(|c| {
                let mut crossing = vec![false; ni * nj];
                for j in 0..nj {
                    for i in 0..ni {
                        let t = grid.cell(Ijk::new(i, j, 0)).top_depth();
                        let b = grid.cell(Ijk::new(i, j, nk - 1)).bottom_depth();
                        crossing[j * ni + i] = c.depth_m >= t && c.depth_m <= b;
                    }
                }
                ContactMask {
                    kind: format!("{:?}", c.kind).to_uppercase(),
                    depth_m: c.depth_m,
                    crossing,
                }
            })
            .collect();

        let outline = vec![self.framework().boundary.ring.clone()];

        // Well markers with per-horizon tie residuals, from provenance (built via
        // `with_well_ties`). The summary `tie_residual_m` is the mean of the well's
        // per-horizon residuals; `ties` carries the full per-horizon breakdown.
        let wells = self
            .provenance()
            .well_ties
            .iter()
            .map(|w| {
                let ties: Vec<WellTieResidual> = w
                    .residuals
                    .iter()
                    .map(|r| WellTieResidual {
                        horizon: r.horizon.clone(),
                        residual_m: r.residual_m,
                    })
                    .collect();
                let tie_residual_m = if ties.is_empty() {
                    None
                } else {
                    Some(ties.iter().map(|t| t.residual_m).sum::<f64>() / ties.len() as f64)
                };
                WellMarker {
                    id: w.id.clone(),
                    x: w.x,
                    y: w.y,
                    tie_residual_m,
                    ties,
                }
            })
            .collect();

        Ok(MapBundle {
            schema_version: SCHEMA_VERSION,
            inputs_ref: self.provenance().inputs_ref.clone(),
            frame,
            outline,
            horizons,
            zone_averages,
            k_slices,
            wells,
            contacts,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Ijk;
    use crate::gridder::{Conformity, SolveOpts};
    use crate::model::{BuildOpts, ConstantPriors, StaticModelBuilder};

    // 4x4x4 flat box: area 10_000 (side 100 -> dx=dy=25), gross 40 (dz=10), top
    // 5000. Logs put phi=0.30 in the upper two layers (tvd 5000..5020), the 0.25
    // prior in the lower two. Contact 5025 lies inside every column.
    fn model() -> StaticModel {
        let opts = BuildOpts {
            area_m2: 10_000.0,
            gross_height_m: 40.0,
            nk: 4,
            conformity: Conformity::Proportional,
            solve_opts: SolveOpts::default(),
            priors: ConstantPriors {
                porosity: 0.25,
                net_to_gross: 0.8,
                water_saturation: 0.3,
            },
        };
        let logs = vec![
            (5001.0, 0.30, 0.20),
            (5005.0, 0.30, 0.20),
            (5015.0, 0.30, 0.20),
        ];
        StaticModelBuilder::flat(4, 4, 5000.0, 5025.0, opts)
            .unwrap()
            .with_logs(logs)
            .build()
            .unwrap()
    }

    // A UTM-origin (world-frame) synthetic model: a flat 11x11 Top at 5000 on a
    // WORLD boundary ring (UTM31N-magnitude), built through the same
    // `from_wireframe` path real data takes, then given its registered world
    // georeference. The grid stays a local area-scaled square (side 300, dx=30);
    // the georef labels column (0,0)'s world centroid (431015, 6521015) + spacing.
    const UTM_X0: f64 = 431_000.0; // grid corner / cell edge (UTM easting)
    const UTM_Y0: f64 = 6_521_000.0; // grid corner / cell edge (UTM northing)
    const UTM_INC: f64 = 30.0; // column spacing (side 300 / 10 columns)

    fn utm_model() -> StaticModel {
        use crate::wireframe::{
            Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
        };
        let (nc, nr) = (11usize, 11usize);
        let wf = Wireframe {
            boundary: Boundary {
                ring: vec![
                    [UTM_X0, UTM_Y0],
                    [UTM_X0 + 300.0, UTM_Y0],
                    [UTM_X0 + 300.0, UTM_Y0 + 300.0],
                    [UTM_X0, UTM_Y0 + 300.0],
                    [UTM_X0, UTM_Y0],
                ],
                hardness: Hardness::Hard,
            },
            horizons: std::sync::Arc::new(vec![Horizon {
                name: "TopRes".into(),
                role: HorizonRole::Top,
                surface: GriddedDepth {
                    ncol: nc,
                    nrow: nr,
                    depth_m: vec![5000.0; nc * nr],
                    is_control: vec![true; nc * nr],
                },
            }]),
            contacts: vec![Contact {
                kind: ContactKind::Owc,
                depth_m: 5025.0,
                hardness: Hardness::Hard,
            }],
        };
        let opts = BuildOpts {
            area_m2: 90_000.0, // side 300 -> dx = 30 over 10 columns
            gross_height_m: 40.0,
            nk: 4,
            conformity: Conformity::Proportional,
            solve_opts: SolveOpts::default(),
            priors: ConstantPriors {
                porosity: 0.25,
                net_to_gross: 0.8,
                water_saturation: 0.3,
            },
        };
        StaticModelBuilder::from_wireframe(&wf, opts)
            .unwrap()
            // Column (0,0)'s world centroid is half a cell in from the corner.
            .with_georef(
                UTM_X0 + UTM_INC / 2.0,
                UTM_Y0 + UTM_INC / 2.0,
                UTM_INC,
                UTM_INC,
            )
            .build()
            .unwrap()
    }

    #[test]
    fn utm_map_frame_is_world_not_local() {
        // F5-class fix: a world-georeferenced model emits its map frame in WORLD
        // coordinates (≈ the UTM origin), NOT the grid's local lattice (≈ 15).
        let m = utm_model();
        let b = m.map_bundle(&MapSpec::new()).unwrap();
        let f = &b.frame;
        assert_eq!((f.ncol, f.nrow), (10, 10));
        // Frame origin ≈ the world column-(0,0) centroid, worlds away from local 15.
        assert!(
            (f.origin_x - (UTM_X0 + UTM_INC / 2.0)).abs() < 1e-6,
            "{f:?}"
        );
        assert!(
            (f.origin_y - (UTM_Y0 + UTM_INC / 2.0)).abs() < 1e-6,
            "{f:?}"
        );
        assert!((f.spacing_x - UTM_INC).abs() < 1e-9 && (f.spacing_y - UTM_INC).abs() < 1e-9);
        assert!(f.origin_x > 400_000.0, "frame is world, not local: {f:?}");

        // The world outline now overlays the frame: every outline point falls
        // inside the frame's world extent (centroid lattice ± half a cell = the
        // grid's world cell-edge extent). Under the old local frame these UTM
        // points sat ~431_000 away from a 0..300 lattice — no overlay at all.
        let xmin = f.origin_x - f.spacing_x / 2.0;
        let xmax = f.origin_x + (f.ncol as f64 - 0.5) * f.spacing_x;
        let ymin = f.origin_y - f.spacing_y / 2.0;
        let ymax = f.origin_y + (f.nrow as f64 - 0.5) * f.spacing_y;
        assert!(!b.outline.is_empty());
        for ring in &b.outline {
            for p in ring {
                assert!(
                    p[0] >= xmin - 1e-6 && p[0] <= xmax + 1e-6,
                    "outline x {} outside [{xmin}, {xmax}]",
                    p[0]
                );
                assert!(
                    p[1] >= ymin - 1e-6 && p[1] <= ymax + 1e-6,
                    "outline y {} outside [{ymin}, {ymax}]",
                    p[1]
                );
            }
        }
    }

    #[test]
    fn map_frame_matches_grid_world_frame() {
        // Georef round-trip: the bundle frame == the grid's own centroid lattice.
        let m = model();
        let b = m.map_bundle(&MapSpec::new()).unwrap();
        let g = m.grid();
        let c00 = g.cell(Ijk::new(0, 0, 0)).centroid();
        let c10 = g.cell(Ijk::new(1, 0, 0)).centroid();
        let c01 = g.cell(Ijk::new(0, 1, 0)).centroid();
        assert_eq!(b.frame.ncol, 4);
        assert_eq!(b.frame.nrow, 4);
        assert!((b.frame.origin_x - c00.x).abs() < 1e-9);
        assert!((b.frame.origin_y - c00.y).abs() < 1e-9);
        assert!((b.frame.spacing_x - (c10.x - c00.x)).abs() < 1e-9);
        assert!((b.frame.spacing_y - (c01.y - c00.y)).abs() < 1e-9);
        assert!((b.frame.spacing_x - 25.0).abs() < 1e-9);
    }

    #[test]
    fn structural_surfaces_read_the_grid_envelope() {
        let m = model();
        let b = m.map_bundle(&MapSpec::new()).unwrap();
        assert_eq!(b.horizons.len(), 2);
        // Top surface == 5000 everywhere, base == 5040 (top + gross).
        assert!(b.horizons[0]
            .values
            .iter()
            .all(|d| (d - 5000.0).abs() < 1e-6));
        assert!(b.horizons[1]
            .values
            .iter()
            .all(|d| (d - 5040.0).abs() < 1e-6));
        assert!((b.horizons[0].range.min - 5000.0).abs() < 1e-6);
    }

    #[test]
    fn k_slice_map_equals_the_cube_layer() {
        let m = model();
        let b = m
            .map_bundle(&MapSpec::new().property("PORO").k_slice(0))
            .unwrap();
        let cube = &m.property("PORO").unwrap().values;
        let dims = m.grid().dims();
        let (ni, nj) = (dims.ni, dims.nj);
        assert_eq!(b.k_slices.len(), 1);
        let slice = &b.k_slices[0];
        assert_eq!(slice.name, "PORO::k0");
        for j in 0..nj {
            for i in 0..ni {
                let expect = cube[j * ni + i];
                assert!((slice.values[j * ni + i] - expect).abs() < 1e-12);
            }
        }
        // Upper layer is the log phi (0.30), distinct from the 0.25 prior.
        assert!(slice.values.iter().all(|v| (v - 0.30).abs() < 1e-9));
    }

    #[test]
    fn zone_average_is_the_hand_checked_column_mean() {
        // RESERVOIR spans all 4 layers: (0.30 + 0.30 + 0.25 + 0.25) / 4 = 0.275.
        let m = model();
        let b = m.map_bundle(&MapSpec::new().property("PORO")).unwrap();
        assert_eq!(b.zone_averages.len(), 1);
        let za = &b.zone_averages[0];
        assert_eq!(za.name, "PORO::RESERVOIR");
        for v in &za.values {
            assert!((v - 0.275).abs() < 1e-9, "zone avg {v} != 0.275");
        }
        // Non-trivial: the interval average differs from any single k-slice.
        assert!((za.range.min - 0.275).abs() < 1e-9 && (za.range.max - 0.275).abs() < 1e-9);
    }

    #[test]
    fn contact_mask_marks_crossed_columns() {
        let m = model();
        let b = m.map_bundle(&MapSpec::new()).unwrap();
        assert_eq!(b.contacts.len(), 1);
        let c = &b.contacts[0];
        assert_eq!(c.kind, "OWC");
        assert!((c.depth_m - 5025.0).abs() < 1e-9);
        // 5025 lies between every column's top (5000) and base (5040).
        assert!(
            c.crossing.iter().all(|&x| x),
            "in-column contact crosses all columns"
        );
        assert_eq!(c.crossing.len(), 16);
    }

    #[test]
    fn missing_property_and_bad_k_slice_error() {
        let m = model();
        assert!(m.map_bundle(&MapSpec::new().property("NOPE")).is_err());
        assert!(m
            .map_bundle(&MapSpec::new().property("PORO").k_slice(99))
            .is_err());
    }

    #[test]
    fn map_bundle_json_round_trips_and_keys_are_stable() {
        let m = model();
        let b = m
            .map_bundle(&MapSpec::new().property("PORO").k_slice(1))
            .unwrap();
        let json = serde_json::to_string(&b).unwrap();
        // Round-trip: a finite (box) bundle survives serialize -> deserialize.
        let back: MapBundle = serde_json::from_str(&json).unwrap();
        assert_eq!(b, back);
        // Schema snapshot: the top-level keys are the cross-codebase contract.
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = v.as_object().unwrap();
        let mut keys: Vec<&str> = obj.keys().map(String::as_str).collect();
        keys.sort_unstable();
        assert_eq!(
            keys,
            [
                "contacts",
                "frame",
                "horizons",
                "inputs_ref",
                "k_slices",
                "outline",
                "schema_version",
                "wells",
                "zone_averages",
            ]
        );
        assert_eq!(obj["schema_version"], serde_json::json!(5));
        // A layer's own key structure (legend metadata included).
        let layer = obj["horizons"][0].as_object().unwrap();
        let mut lkeys: Vec<&str> = layer.keys().map(String::as_str).collect();
        lkeys.sort_unstable();
        assert_eq!(lkeys, ["name", "range", "units", "values"]);
    }
}
