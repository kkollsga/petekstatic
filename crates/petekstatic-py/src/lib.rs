//! Thin PyO3 bindings over petekStatic — the minimal `petekstatic` Python
//! surface. Code-first: build a flat synthetic static reservoir model, then read
//! its volumes (`in_place` / `in_place_by_zone`) and the JSON view bundles
//! (`map_bundle` / `intersection_bundle` / `volume_bundle`). Logic lives in the
//! Rust core crates (`srs-model` and below); this file only marshals across the
//! Python boundary. The rich Python surface is `peteksim` (the petekSim repo);
//! this wheel is the essentials only.
//!
//! **Units (SI):** area in **m²**, depths/lengths in **m** (positive down),
//! volumes in **m³** / **Sm³**. FVF (`boi`/`bgi`) enters as a validated scalar.

use petekstatic_error::StaticError;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use srs_gridder::{Conformity, SolveOpts};
use srs_model::{
    ConstantPriors, GasFvf, InPlace, OilFvf, SectionSpec, StaticModel as CoreModel,
    StaticModelBuilder,
};
use srs_wireframe::{
    Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
};

use srs_model::BuildOpts;

/// Map a [`StaticError`] into a Python `ValueError`.
fn py_err(e: StaticError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Build a flat, single-horizon synthetic wireframe: an `n x n` top surface at a
/// constant `depth_m` over a unit square boundary, with one OWC at `owc_m`.
fn flat_wireframe(n: usize, depth_m: f64, owc_m: f64) -> Wireframe {
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
                depth_m: vec![depth_m; n * n],
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

/// A populated static reservoir model — the Python handle over the Rust
/// [`srs_model::StaticModel`]. Construct one with [`build_flat_model`].
#[pyclass(name = "StaticModel")]
pub struct StaticModel {
    inner: CoreModel,
}

impl StaticModel {
    /// Assemble a `{grv_m3, hcpv_m3, cells_in_column, ooip_sm3}` dict from an
    /// [`InPlace`] result and the caller's oil FVF; adds `ogip_sm3` when a gas FVF
    /// (`bgi`, in (0,1)) is supplied.
    fn in_place_dict(
        py: Python<'_>,
        ip: &InPlace,
        boi: f64,
        bgi: Option<f64>,
    ) -> PyResult<Py<PyDict>> {
        let d = PyDict::new(py);
        d.set_item("grv_m3", ip.grv_m3)?;
        d.set_item("hcpv_m3", ip.hcpv_m3)?;
        d.set_item("cells_in_column", ip.cells_in_column)?;
        d.set_item("ooip_sm3", ip.ooip_sm3(OilFvf::new(boi).map_err(py_err)?))?;
        if let Some(bgi) = bgi {
            d.set_item("ogip_sm3", ip.ogip_sm3(GasFvf::new(bgi).map_err(py_err)?))?;
        }
        Ok(d.unbind())
    }
}

#[pymethods]
impl StaticModel {
    /// The names of the populated property cubes.
    fn property_names(&self) -> Vec<String> {
        self.inner
            .property_names()
            .into_iter()
            .map(String::from)
            .collect()
    }

    /// Gross rock volume of the whole grid \[m³\].
    fn bulk_volume(&self) -> f64 {
        self.inner.bulk_volume()
    }

    /// Whole-column in-place volumes as a dict: `grv_m3`, `hcpv_m3`,
    /// `cells_in_column`, `ooip_sm3` (HCPV / `boi`, oil FVF `>= 1.0`, default
    /// 1.0), and — when a gas FVF `bgi` in (0,1) is supplied — `ogip_sm3`
    /// (HCPV / `bgi`).
    #[pyo3(signature = (boi=1.0, bgi=None))]
    fn in_place(&self, py: Python<'_>, boi: f64, bgi: Option<f64>) -> PyResult<Py<PyDict>> {
        let ip = self.inner.in_place().map_err(py_err)?;
        Self::in_place_dict(py, &ip, boi, bgi)
    }

    /// Summary-only in-place (no per-cell HCPV cube) — same aggregate dict shape
    /// as [`StaticModel::in_place`].
    #[pyo3(signature = (boi=1.0, bgi=None))]
    fn in_place_summary(&self, py: Python<'_>, boi: f64, bgi: Option<f64>) -> PyResult<Py<PyDict>> {
        let ip = self.inner.in_place_summary().map_err(py_err)?;
        Self::in_place_dict(py, &ip, boi, bgi)
    }

    /// Per-zone in-place with a total rollup, as a JSON string:
    /// `{"zones": [{"zone", "grv_m3", "hcpv_m3", "ooip_sm3"}...], "total": {...}}`.
    #[pyo3(signature = (boi=1.0))]
    fn in_place_by_zone(&self, boi: f64) -> PyResult<String> {
        let boi = OilFvf::new(boi).map_err(py_err)?;
        let zoned = self.inner.in_place_by_zone().map_err(py_err)?;
        let zones: Vec<_> = zoned
            .zones
            .iter()
            .map(|z| {
                serde_json::json!({
                    "zone": z.zone,
                    "grv_m3": z.in_place.grv_m3,
                    "hcpv_m3": z.in_place.hcpv_m3,
                    "ooip_sm3": z.in_place.ooip_sm3(boi),
                })
            })
            .collect();
        let out = serde_json::json!({
            "zones": zones,
            "total": {
                "grv_m3": zoned.total.grv_m3,
                "hcpv_m3": zoned.total.hcpv_m3,
                "ooip_sm3": zoned.total.ooip_sm3(boi),
            },
        });
        Ok(out.to_string())
    }

    /// The areal (plan-view) map bundle as a JSON string. Pass a `property` name
    /// to include its zone-average maps (plus a k-slice map when `k_slice` is set).
    #[pyo3(signature = (property=None, k_slice=None))]
    fn map_bundle(&self, property: Option<&str>, k_slice: Option<usize>) -> PyResult<String> {
        let mut spec = srs_model::MapSpec::new();
        if let Some(name) = property {
            spec = spec.property(name);
        }
        if let Some(k) = k_slice {
            spec = spec.k_slice(k);
        }
        let bundle = self.inner.map_bundle(&spec).map_err(py_err)?;
        serde_json::to_string(&bundle).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// A vertical cross-section along a world `[x, y]` polyline, as a JSON string.
    /// Pass a `property` name to carry its per-layer values.
    #[pyo3(signature = (line, property=None))]
    fn intersection_bundle(&self, line: Vec<[f64; 2]>, property: Option<&str>) -> PyResult<String> {
        let spec = SectionSpec::Polyline(line);
        let bundle = self
            .inner
            .intersection_bundle(&spec, property)
            .map_err(py_err)?;
        serde_json::to_string(&bundle).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// The corner-point exterior-shell volume bundle for `property`, as a JSON
    /// string (self-contained envelope with base64-wrapped binary blocks).
    fn volume_bundle(&self, property: &str) -> PyResult<String> {
        let bundle = self.inner.volume_bundle(property).map_err(py_err)?;
        serde_json::to_string(&bundle).map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

/// Build a flat synthetic single-zone static model and populate it with constant
/// priors — the minimal path from Python to a model you can read volumes off.
///
/// `n` is the top-surface node count per side; `depth_m` the flat top depth;
/// `owc_m` the oil-water contact; `area_m2` the areal footprint; `gross_height_m`
/// the column thickness; `nk` the layer count; and `porosity` / `net_to_gross` /
/// `water_saturation` the day-1 constant fraction priors.
#[pyfunction]
#[pyo3(signature = (
    n=11,
    depth_m=2000.0,
    owc_m=2100.0,
    area_m2=1_000_000.0,
    gross_height_m=50.0,
    nk=5,
    porosity=0.25,
    net_to_gross=0.8,
    water_saturation=0.3,
))]
#[allow(clippy::too_many_arguments)]
fn build_flat_model(
    n: usize,
    depth_m: f64,
    owc_m: f64,
    area_m2: f64,
    gross_height_m: f64,
    nk: usize,
    porosity: f64,
    net_to_gross: f64,
    water_saturation: f64,
) -> PyResult<StaticModel> {
    let wf = flat_wireframe(n, depth_m, owc_m);
    let opts = BuildOpts {
        area_m2,
        gross_height_m,
        nk,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity,
            net_to_gross,
            water_saturation,
        },
    };
    let inner = StaticModelBuilder::from_wireframe(&wf, opts)
        .map_err(py_err)?
        .build()
        .map_err(py_err)?;
    Ok(StaticModel { inner })
}

/// The `_petekstatic` extension module.
#[pymodule]
fn _petekstatic(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<StaticModel>()?;
    m.add_function(wrap_pyfunction!(build_flat_model, m)?)?;
    Ok(())
}
