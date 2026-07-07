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

use petekstatic::error::StaticError;
use petekstatic::gridder::{Conformity, SolveOpts};
use petekstatic::model::{
    ConstantPriors, GasFvf, Gaussian, InPlace, OilFvf, PropertyPipeline as CorePropertyPipeline,
    SectionSpec, StaticModel as CoreModel, StaticModelBuilder, UpscaleMethod,
    WellLog as CoreWellLog,
};
use petekstatic::wireframe::{
    Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
};
use petektools::{AnisotropicVariogram, VariogramModel};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use petekstatic::model::BuildOpts;

/// Map a [`StaticError`] into a Python `ValueError`.
fn py_err(e: StaticError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn parse_upscale_method(method: &str) -> PyResult<UpscaleMethod> {
    match method.trim().to_ascii_lowercase().as_str() {
        "arithmetic" | "arith" => Ok(UpscaleMethod::Arithmetic),
        "harmonic" => Ok(UpscaleMethod::Harmonic),
        "geometric" | "geo" => Ok(UpscaleMethod::Geometric),
        other => Err(PyValueError::new_err(format!(
            "unknown upscale method '{other}' (expected 'arithmetic', 'harmonic', or 'geometric')"
        ))),
    }
}

fn parse_variogram_model(model: &str) -> PyResult<VariogramModel> {
    match model.trim().to_ascii_lowercase().as_str() {
        "spherical" | "sph" => Ok(VariogramModel::Spherical),
        "exponential" | "exp" => Ok(VariogramModel::Exponential),
        "gaussian" | "gauss" => Ok(VariogramModel::Gaussian),
        other => Err(PyValueError::new_err(format!(
            "unknown variogram model '{other}' (expected 'spherical', 'exponential', or 'gaussian')"
        ))),
    }
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

/// A positioned well log for property-pipeline construction.
#[pyclass(name = "WellLog", frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct WellLog {
    inner: CoreWellLog,
}

#[pymethods]
impl WellLog {
    /// Build a positioned well log from world `x/y` and `(depth_m, value)` samples.
    #[new]
    fn new(x: f64, y: f64, samples: Vec<(f64, f64)>) -> PyResult<Self> {
        if !(x.is_finite() && y.is_finite()) {
            return Err(PyValueError::new_err("WellLog x/y must be finite"));
        }
        if samples.is_empty() {
            return Err(PyValueError::new_err("WellLog samples must not be empty"));
        }
        for (depth, value) in &samples {
            if !(depth.is_finite() && value.is_finite()) {
                return Err(PyValueError::new_err(
                    "WellLog samples must contain finite depth/value pairs",
                ));
            }
        }
        Ok(Self {
            inner: CoreWellLog::new(x, y, samples),
        })
    }

    #[getter]
    fn x(&self) -> f64 {
        self.inner.x
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.y
    }

    #[getter]
    fn samples(&self) -> Vec<(f64, f64)> {
        self.inner.samples.clone()
    }

    fn as_dict(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let d = PyDict::new(py);
        d.set_item("x", self.inner.x)?;
        d.set_item("y", self.inner.y)?;
        d.set_item("samples", self.inner.samples.clone())?;
        Ok(d.unbind())
    }

    fn __repr__(&self) -> String {
        format!(
            "WellLog(x={}, y={}, samples={})",
            self.inner.x,
            self.inner.y,
            self.inner.samples.len()
        )
    }
}

/// A Rust-backed property-population pipeline handle.
///
/// This exposes the executable construction boundary from Python: positioned well
/// logs, an upscale method, a variogram, and an SGS seed are lowered to
/// the production Rust [`PropertyPipeline`]. Applying to an arbitrary mutable grid
/// is intentionally not exposed until the grid binding exists; use
/// `apply_to_flat_model` for a small smoke execution through the existing builder.
#[pyclass(name = "PropertyPipeline", skip_from_py_object)]
#[derive(Clone)]
pub struct PropertyPipeline {
    inner: CorePropertyPipeline,
    name: String,
    well_count: usize,
    method: String,
    variogram_model: String,
    range_m: f64,
    minor_m: f64,
    vertical_m: f64,
    azimuth: f64,
    sill: f64,
    nugget: f64,
    seed: u64,
    propagate: bool,
    allow_mean_fill: bool,
    max_neighbours: Option<usize>,
    radius_m: Option<f64>,
    unbounded_search: bool,
}

struct PipelineBuildConfig<'a> {
    name: &'a str,
    wells: &'a [CoreWellLog],
    method: UpscaleMethod,
    model: VariogramModel,
    range_m: f64,
    minor_m: f64,
    vertical_m: f64,
    azimuth: f64,
    sill: f64,
    nugget: f64,
    seed: u64,
    propagate: bool,
    allow_mean_fill: bool,
    max_neighbours: Option<usize>,
    radius_m: Option<f64>,
    unbounded_search: bool,
}

impl PropertyPipeline {
    fn build_inner(config: &PipelineBuildConfig<'_>) -> PyResult<CorePropertyPipeline> {
        let variogram = AnisotropicVariogram::new(
            config.model,
            config.nugget,
            config.sill,
            config.range_m,
            config.minor_m,
            config.vertical_m,
            config.azimuth,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut pipe = CorePropertyPipeline::new(config.name.to_string())
            .upscale(config.wells.to_vec(), config.method);
        if config.propagate {
            let mut gaussian = Gaussian::new(variogram, config.seed);
            if config.allow_mean_fill {
                gaussian = gaussian.allow_mean_fill();
            }
            match (config.max_neighbours, config.radius_m) {
                (Some(max), Some(radius)) => {
                    gaussian = gaussian.with_search(max, radius);
                }
                (None, None) => {}
                _ => {
                    return Err(PyValueError::new_err(
                        "max_neighbours and radius_m must be provided together",
                    ));
                }
            }
            if config.unbounded_search {
                gaussian = gaussian.with_unbounded_search();
            }
            pipe = pipe.propagate(gaussian);
        }
        Ok(pipe)
    }
}

#[pymethods]
impl PropertyPipeline {
    /// Build a Rust property pipeline from positioned well logs and a variogram.
    /// `range_m` is the major range; omit `minor_m`, `vertical_m`, and `azimuth`
    /// for isotropic SGS. `sill` is the partial sill; `nugget` is the nugget effect.
    #[new]
    #[pyo3(signature = (
        name,
        wells,
        method,
        variogram_model,
        range_m,
        seed,
        sill=1.0,
        nugget=0.0,
        minor_m=None,
        vertical_m=None,
        azimuth=0.0,
        propagate=true,
        allow_mean_fill=false,
        max_neighbours=None,
        radius_m=None,
        unbounded_search=false,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        name: String,
        wells: Vec<PyRef<'_, WellLog>>,
        method: &str,
        variogram_model: &str,
        range_m: f64,
        seed: u64,
        sill: f64,
        nugget: f64,
        minor_m: Option<f64>,
        vertical_m: Option<f64>,
        azimuth: f64,
        propagate: bool,
        allow_mean_fill: bool,
        max_neighbours: Option<usize>,
        radius_m: Option<f64>,
        unbounded_search: bool,
    ) -> PyResult<Self> {
        let name = name.trim().to_string();
        if name.is_empty() {
            return Err(PyValueError::new_err(
                "PropertyPipeline name must be non-empty",
            ));
        }
        if wells.is_empty() {
            return Err(PyValueError::new_err(
                "PropertyPipeline requires at least one WellLog",
            ));
        }
        let minor_m = minor_m.unwrap_or(range_m);
        let vertical_m = vertical_m.unwrap_or(range_m);
        if !range_m.is_finite()
            || range_m <= 0.0
            || !minor_m.is_finite()
            || minor_m <= 0.0
            || !vertical_m.is_finite()
            || vertical_m <= 0.0
        {
            return Err(PyValueError::new_err(
                "PropertyPipeline variogram ranges must be finite and > 0",
            ));
        }
        if !azimuth.is_finite() {
            return Err(PyValueError::new_err(
                "PropertyPipeline azimuth must be finite",
            ));
        }
        if !sill.is_finite() || sill < 0.0 || !nugget.is_finite() || nugget < 0.0 {
            return Err(PyValueError::new_err(
                "PropertyPipeline sill/nugget must be finite and >= 0",
            ));
        }
        if let Some(radius) = radius_m {
            if !radius.is_finite() || radius <= 0.0 {
                return Err(PyValueError::new_err("radius_m must be finite and > 0"));
            }
        }
        if matches!(max_neighbours, Some(0)) {
            return Err(PyValueError::new_err("max_neighbours must be > 0"));
        }
        let method_value = parse_upscale_method(method)?;
        let model_value = parse_variogram_model(variogram_model)?;
        let core_wells: Vec<_> = wells.iter().map(|w| w.inner.clone()).collect();
        let inner = Self::build_inner(&PipelineBuildConfig {
            name: &name,
            wells: &core_wells,
            method: method_value,
            model: model_value,
            range_m,
            minor_m,
            vertical_m,
            azimuth,
            sill,
            nugget,
            seed,
            propagate,
            allow_mean_fill,
            max_neighbours,
            radius_m,
            unbounded_search,
        })?;
        Ok(Self {
            inner,
            name,
            well_count: core_wells.len(),
            method: method.trim().to_ascii_lowercase(),
            variogram_model: variogram_model.trim().to_ascii_lowercase(),
            range_m,
            minor_m,
            vertical_m,
            azimuth: azimuth.rem_euclid(360.0),
            sill,
            nugget,
            seed,
            propagate,
            allow_mean_fill,
            max_neighbours,
            radius_m,
            unbounded_search,
        })
    }

    #[getter]
    fn name(&self) -> &str {
        &self.name
    }

    fn config(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let d = PyDict::new(py);
        d.set_item("property", &self.name)?;
        d.set_item("well_count", self.well_count)?;
        d.set_item("method", &self.method)?;
        d.set_item("variogram_model", &self.variogram_model)?;
        d.set_item("range_m", self.range_m)?;
        d.set_item("major_m", self.range_m)?;
        d.set_item("minor_m", self.minor_m)?;
        d.set_item("vertical_m", self.vertical_m)?;
        d.set_item("azimuth", self.azimuth)?;
        d.set_item("sill", self.sill)?;
        d.set_item("nugget", self.nugget)?;
        d.set_item("seed", self.seed)?;
        d.set_item("propagate", self.propagate)?;
        d.set_item("allow_mean_fill", self.allow_mean_fill)?;
        d.set_item("max_neighbours", self.max_neighbours)?;
        d.set_item("radius_m", self.radius_m)?;
        d.set_item("unbounded_search", self.unbounded_search)?;
        Ok(d.unbind())
    }

    fn report(&self) -> String {
        let step = if self.propagate {
            "upscale+gaussian"
        } else {
            "upscale"
        };
        format!(
            "PropertyPipeline(property='{}', step={}, wells={}, method={}, variogram={} major_m={} minor_m={} vertical_m={} azimuth={} sill={} nugget={} seed={})",
            self.name,
            step,
            self.well_count,
            self.method,
            self.variogram_model,
            self.range_m,
            self.minor_m,
            self.vertical_m,
            self.azimuth,
            self.sill,
            self.nugget,
            self.seed,
        )
    }

    /// Smoke-apply this pipeline through the existing flat-model builder.
    #[pyo3(signature = (
        n=11,
        depth_m=1000.0,
        owc_m=1100.0,
        area_m2=10_000.0,
        gross_height_m=20.0,
        nk=2,
        porosity=0.25,
        net_to_gross=0.8,
        water_saturation=0.3,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn apply_to_flat_model(
        &self,
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
            .with_property(self.inner.clone())
            .build()
            .map_err(py_err)?;
        Ok(StaticModel { inner })
    }

    fn __repr__(&self) -> String {
        self.report()
    }
}

/// A populated static reservoir model — the Python handle over the Rust
/// [`petekstatic::model::StaticModel`]. Construct one with [`build_flat_model`].
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
        let mut spec = petekstatic::model::MapSpec::new();
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
    m.add_class::<WellLog>()?;
    m.add_class::<PropertyPipeline>()?;
    m.add_function(wrap_pyfunction!(build_flat_model, m)?)?;
    Ok(())
}
