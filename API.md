# petekStatic — public API

> **Status: WORKING CONTRACT.** The geomodel crates are green and the headline
> **`StaticModel`** surface (`srs-model`) has **LANDED** per the ratified seams
> (SPEC §7/§8: `decision_layer_charters`, `decision_staticmodel_regen_seam`).
> This file is the contract: changing a cross-library signature (the StaticModel
> accessors, `RealizationDraw`, the regeneration API) requires coordinator +
> consumer sign-off; library-internal signatures are petekStatic's own call. The
> file locks fully at the 0.1 release. Rust is canonical. Python exposes both
> the compact flat-model surface and the first petekStatic-owned static workflow
> facade (`Grid.from_project(...).geometry(...).horizons(...).zones(...).layers(...)`).
> The workflow facade is currently a declarative/spec layer with in-memory
> property arrays, petekTools formula evaluation, deterministic simple
> volumetrics, and lowering of `upscale(...).sgs(...)` property declarations to
> `PropertyPipelineSpec`; it is not yet the full Rust corner-point grid
> construction path.

Conventions (house style): `Result<T> = std::result::Result<T, StaticError>`;
per-cell cubes are `Vec<f64>` indexed by linear cell index; `NaN` = undefined;
depths are metres, positive-down (the optional petekIO adapter flips upstream
negative-down elevation at its ingest boundary).

---

## Python workflow facade — first vertical slice

```python
import petekstatic as pst

grid = (
    pst.Grid.from_project(project)
    .geometry(cell=(50.0, 50.0), orient=0.0, outline="ModelEdge")
    .horizons(
        [
            {
                "name": "Top reservoir",
                "surface": "Top reservoir input surface",
                "well top": "well tops/Top reservoir",
                "zone": {
                    "name": "Reservoir",
                    "sub-zones": [
                        {"zone": "Top Reservoir", "type": "constant"},
                        {"name": "Intra Shale", "well top": "Top Lower Reservoir"},
                        {"name": "Lower Reservoir", "type": "isochore"},
                    ],
                },
            },
            "Base reservoir",
            {"name": "Custom model horizon name", "surface": "input surface"},
        ],
        well_tie={"influence_radius": 800},
    )
    .layers({"Top Reservoir": pst.Layering(n=2), "Lower Reservoir": pst.Layering(n=2)})
)

p = grid.properties
p.ntg = 0.80
p.por = p.ntg * 0.85
p.sw = 0.20
p["PermXY_BC"].set(100.0)
p["PorE_BC"].set(0.25)
p.calc(
    ["RQI = $lambda * sqrt(PermXY_BC / PorE_BC)"],
    params={"lambda": 0.0314},
)

case = grid.volumes(ntg="NTG", por="POR", sw="Sw", fluid="oil", fvf=1.30)
result = case.run(progress=True)
```

Exported Python names:

```python
Grid, HorizonSpec, WellTie, Layering, Spherical,
PropertyStore, PropertyHandle, PropertyPipelineSpec,
Var, WellLogSpec, DistributionSpec, CoKriging,
UpscaleRecipeBuilder, SgsRecipe, distributions, upscale,
WellLog, PropertyPipeline, VolumeCase, VolumeResult,
StaticModel, build_flat_model, __version__
```

Current behavior:
- `Grid.from_project(project)` accepts the petekIO `Project` facade or a
  project-like fixture.
- `geometry(..., outline=...)` and `horizons([...])` validate names against
  `Project.inventory()` and mapping-like project collections. Missing assets
  raise `ValueError` with available names.
- `grid.properties` stores dense in-memory property vectors. Handles are
  available as `p["NAME"]` and convenience attributes like `p.ntg`, `p.por`,
  `p.sw`. Constants can be written as scalars (`p.ntg = 0.8` or
  `p.ntg.set(0.8)`) and are broadcast to the declared cell count; pass an
  iterable when assigning explicit per-cell values. Property handles also
  compose into assignment expressions (`p.por = p.ntg * 0.85`), routed through
  the same formula evaluator as `p.calc(...)`.
- `p.calc([...], params={...})` delegates the whole formula block to
  `petektools.evaluate_formula` and writes outputs only after the block succeeds.
- `pst.Var(model, major, minor, vertical, azimuth, sill=None, nugget=None)`
  records the canonical anisotropic variogram spec. `pst.Spherical(range_m)` is
  accepted as isotropic shorthand and converts to `Var`.
- `pst.distributions.from_logs()` records that the SGS target distribution is
  estimated from the resolved log samples. It is also the default for log-channel
  sources when no distribution is supplied.
- Assigning an `pst.upscale(source).sgs(...)` recipe, for example
  `p.por = pst.upscale(logs.PHIE(logs.NetSand > 0.50)).sgs(...)`, records the
  declaration and lowers it to `PropertyPipelineSpec`. Inspect declarations with
  `p.declarations("por")`, inspect lowered specs with `p.pipelines("por")`, and
  access the object with `p.pipeline_spec("por")`.
- `grid.volumes(...).run(progress=True|callback)` computes deterministic simple
  GRV/HCPV/in-place volumes from the declared cell size, horizon thickness, layer
  count, and named `ntg`/`por`/`sw` arrays. `VolumeResult.summary()` and
  `VolumeResult.by_zone()` return dictionaries.

Property recipe example:

```python
logs = project.logs
vgm = pst.Var(
    "spherical",
    major=1500,
    minor=700,
    vertical=20,
    azimuth=35,
    sill=1.2,
    nugget=0.05,
)

p.por = pst.upscale(logs.PHIE(logs.NetSand > 0.50)).sgs(
    variogram=vgm,
    distribution=pst.distributions.from_logs(),
    seed=12,
)

lowered = p.pipelines("por")
```

Executable recipe example, when `logs.NetSand` or the project resolver returns
positioned well logs:

```python
iso = pst.Var("spherical", major=1500, minor=1500, vertical=1500, azimuth=0)
p.ntg = pst.upscale(logs.NetSand).sgs(
    variogram=iso,
    distribution=pst.distributions.from_logs(),
    seed=11,
)

pipe = p.execute_pipeline("ntg")
config = pipe.config()
smoke_model = pipe.apply_to_flat_model()
```

Execution boundary:
- If `source.to_well_logs(project)`, `source.resolve_well_logs(project)`,
  `project.resolve_log_expression(source)`, `project.resolve_well_logs(source)`,
  or `project.resolve_log_source(source_dict)` returns positioned logs, the
  lowered spec stores `WellLogSpec` inputs.
- `p.execute_pipeline("por")` returns a Rust-backed `pst.PropertyPipeline` handle
  only when positioned wells are available, no cokriging/trend is bound,
  `distribution=pst.distributions.from_logs()`, and the variogram is isotropic
  for current Rust execution (`major == minor == vertical`, `azimuth == 0`).
- `PropertyPipeline.apply_to_flat_model(...)` is the current smoke execution
  path. Applying a pipeline to an arbitrary mutable production grid is not
  exposed through Python yet.
- Lazy unresolved log expressions, cokriging/trend binding, non-`from_logs`
  distributions, and anisotropic Rust execution raise explicit
  `NotImplementedError`s. Anisotropic `Var` specs are still serialized intact and
  `PropertyPipelineSpec.to_petektools_variogram()` lowers them to the petekTools
  anisotropic variogram object when petekTools is installed.

---

## `error`

```rust
pub enum StaticError {
    InvalidInput(String),
    Grid(String),
    OutOfRange(String),
    CrossedSurfaces { nodes: usize, worst_m: f64 },   // R1: base crosses above the top (GRV-collapsing); opt out via with_clamp_base_to_top
    #[cfg(feature = "petekio-adapter")]
    Geo(#[from] petekio::GeoError),   // optional DATA→GEOMODEL compatibility seam
}
pub type Result<T> = std::result::Result<T, StaticError>;
```

Downstream composes it: `SrsError::Static(#[from] petekstatic::error::StaticError)`.
When `petekio-adapter` is enabled, that chain also reaches `GeoError`.

## `srs-grid` — the i,j,k corner-point grid

```rust
pub struct Grid { /* geom + properties + layers, all private */ }
impl Grid {
    pub fn new(geom: CornerPointGeom) -> Grid;
    pub fn dims(&self) -> Dims;
    pub fn cell_count(&self) -> usize;
    pub fn cell(&self, c: Ijk) -> Cell;
    pub fn cells(&self) -> impl Iterator<Item = Cell> + '_;
    pub fn bulk_volume(&self) -> f64;                 // gross rock volume, m³
    pub fn properties(&self) -> &Properties;
    pub fn properties_mut(&mut self) -> &mut Properties;
    pub fn layers(&self) -> &[KLayer];
}

pub struct Cell { pub ijk: Ijk, pub corners: [Point3; 8] }
impl Cell {
    pub fn volume(&self) -> f64;
    pub fn centroid(&self) -> Point3;
    pub fn top_depth(&self) -> f64;                   // mean of 4 top corners
    pub fn bottom_depth(&self) -> f64;                // mean of 4 bottom corners
    pub fn dz(&self) -> f64;
}

pub struct Property { pub name: String, pub values: Vec<f64> }
impl Property { pub fn constant(name: impl Into<String>, value: f64, cell_count: usize) -> Property; }

pub struct Properties { /* private */ }
impl Properties {
    pub fn new(cell_count: usize) -> Properties;
    pub fn set(&mut self, prop: Property) -> Result<()>;   // errors if len != cell_count
    pub fn get(&self, name: &str) -> Option<&Property>;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    pub fn names(&self) -> impl Iterator<Item = &str>;     // unordered
}

pub struct Dims { /* ni, nj, nk */ }  pub struct Ijk { /* i, j, k */ }
pub struct CornerPointGeom { /* pillars + node z's */ }  pub struct Pillar { /* ... */ }
pub struct KLayer { /* ... */ }  pub struct Point3 { pub x: f64, pub y: f64, pub z: f64 }
pub struct Segment { /* ... */ }
pub fn build_box(spec: BoxSpec) -> Result<Grid>;   pub struct BoxSpec { /* ... */ }
pub fn hexahedron_volume(corners: &[Point3; 8]) -> f64;
```

## `srs-wireframe` — the constraining framework

```rust
pub enum Hardness { Hard, Interpolated, Assumed }
pub enum HorizonRole { Top, Base, Intermediate }
pub enum ContactKind { Owc, Goc, Gwc }
pub struct Boundary { pub ring: Vec<[f64; 2]>, pub hardness: Hardness }
pub struct GriddedDepth { pub ncol: usize, pub nrow: usize, pub depth_m: Vec<f64>, pub is_control: Vec<bool> }
pub struct Horizon { pub name: String, pub role: HorizonRole, pub surface: GriddedDepth }
pub struct Contact { pub kind: ContactKind, pub depth_m: f64, pub hardness: Hardness }
pub struct Wireframe { pub boundary: Boundary, pub horizons: Arc<Vec<Horizon>>, pub contacts: Vec<Contact> }  // horizons Arc-shared: realization-invariant, O(1) clone per MC realize (mutate via Arc::make_mut)
impl Wireframe { pub fn from_boundary(boundary: Boundary) -> Wireframe; }
```

## `srs-gridder` — the convergent gridder

```rust
pub struct Control { pub ip: usize, pub jp: usize, pub z: f64 }
pub struct SolveOpts { pub tension: f64, pub omega: f64, pub tol: f64, pub max_iter: usize }
impl Default for SolveOpts { /* tension .25, omega 1.5, tol 1e-6, max_iter 20_000 */ }
pub struct Surface { /* nx, ny, z — private */ }
impl Surface {
    pub fn constant(nx: usize, ny: usize, z: f64) -> Surface;
    pub fn offset_by(&self, dz: f64) -> Surface;      // conformable base at constant thickness
    pub fn offset_by_field(&self, dz_m: &[f64]) -> Result<Surface>;  // per-node offset (row-major jp*nx+ip) — real base relief
    pub fn nx(&self) -> usize;  pub fn ny(&self) -> usize;  pub fn z(&self, ip: usize, jp: usize) -> f64;
    // Order-repair family (SPEC §3 R1/R-c + §4a precedence). `guard_below` /
    // `repair_min_thickness` move the LOWER (base yields, top preserved); the
    // repair-precedence twins `guard_above` / `repair_min_thickness_from_below`
    // move the UPPER up (a derived surface yields to a mapped one).
    pub fn guard_below(&self, top: &Surface, clamp: bool) -> Result<Surface>;
    pub fn repair_min_thickness(&self, top: &Surface, min_thickness_m: f64) -> Result<(Surface, usize, f64)>;
    pub fn guard_above(&self, lower: &Surface, clamp: bool) -> Result<Surface>;
    pub fn repair_min_thickness_from_below(&self, lower: &Surface, min_thickness_m: f64) -> Result<(Surface, usize, f64)>;
}
pub fn solve_surface(nx: usize, ny: usize, controls: &[Control], opts: SolveOpts) -> Result<Surface>;

// Extrapolation policy beyond the data hull (SPEC §3, audit S3): explicit,
// owner-visible; default DecayToData { start_cells: 2, decay_cells: 4 }.
pub enum ExtrapolationPolicy { NaturalDip, DecayToData { start_cells: f64, decay_cells: f64 } }
impl Default for ExtrapolationPolicy { /* DecayToData { 2.0, 4.0 } */ }
impl Surface {
    pub fn taper_beyond_data(&self, controls: &[Control], policy: ExtrapolationPolicy) -> Surface;
}

// The structure-build solve entry (SPEC §3, audit S2): plane detrending (kills
// the kernel's slow affine mode) + fixed-point restarts; hard controls exact.
pub fn solve_surface_converged(nx: usize, ny: usize, controls: &[Control]) -> Result<KernelSurface>;

// The petekTools-kernel-space surface — the ONLY admissible warm-start seed.
// Sources: `flat` (constant field — a fixed point of both kernels, the safe
// bootstrap) or a prior `solve_surface_seeded` output. Deliberately NO
// `From<Surface>`: seeding the warm kernel from the cold solver's output is the
// kernel-space violation (`decision_gridder_kernel_unification`) and is now a
// compile error.
pub struct KernelSurface { /* private Surface */ }
impl KernelSurface {
    pub fn flat(nx: usize, ny: usize, z: f64) -> KernelSurface;
    pub fn surface(&self) -> &Surface;
    pub fn nx(&self) -> usize;  pub fn ny(&self) -> usize;  pub fn z(&self, ip: usize, jp: usize) -> f64;
}
// Warm-start refine (SPEC §7a): delegates the seeded SOR to petekTools'
// ConvergentGridder kernel; ~14x faster than a cold solve for a one-control
// perturbation at 50x50. No SolveOpts — petekTools owns its solver parameters.
pub fn solve_surface_seeded(seed: &KernelSurface, controls: &[Control]) -> Result<KernelSurface>;

// Layering conformity. Proportional honours `nk`; the Follow styles derive `nk`
// from geometry (ceil(max column thickness / dz_m), capped at MAX_NK) and drape
// each layer surface parallel to the top (resp. base) at a constant dz_m, with
// the deep (resp. shallow) layers TRUNCATING against the pinch-out horizon.
// A truncated cell collapses to zero thickness (zero volume — excluded from
// volumetrics, NaN-marked in the view bundles); there is no active-mask array.
pub const MAX_NK: usize = 200;
pub enum Conformity { Proportional, FollowTop { dz_m: f64 }, FollowBase { dz_m: f64 } }
// The layered grid + its report. `grid.dims().nk == nk`. `nk` is dz-derived for
// a Follow style, the passed `nk` for Proportional. `truncated_cells` counts the
// zero-thickness collapses; `nk_capped` = the derived count hit MAX_NK.
pub struct LayeredGrid { pub grid: Grid, pub nk: usize, pub truncated_cells: usize, pub nk_capped: bool }
// `nk` honoured for Proportional, ignored (dz-derived) for the Follow styles.
pub fn layer_grid(top: &Surface, base: &Surface, dx: f64, dy: f64, nk: usize, conf: Conformity) -> Result<LayeredGrid>;
```

## `srs-petro` — log upscaling

```rust
pub struct WeightedSample { pub weight: f64, pub value: f64 }
pub fn arithmetic_mean(s: &[WeightedSample]) -> f64;
pub fn geometric_mean(s: &[WeightedSample]) -> f64;
pub fn harmonic_mean(s: &[WeightedSample]) -> f64;
pub fn power_law_mean(s: &[WeightedSample], p: f64) -> f64;
pub struct NetSample { /* ... */ }  pub struct SwSample { pub length: f64, pub porosity: f64, pub water_saturation: f64 }
pub fn upscale_porosity(samples: &[WeightedSample]) -> Result<f64>;   // length-weighted
pub fn upscale_sw(samples: &[SwSample]) -> Result<f64>;               // pore-volume-weighted
pub fn net_to_gross(samples: &[NetSample]) -> f64;
pub fn perm_bounds(/* ... */) -> (f64, f64);                          // Cardwell-Parsons bracket
```

## Optional `petekio-adapter` compatibility feature

The geomodel core has no petekIO dependency. Enable `petekio-adapter` only for
the legacy model-ready-input conversion surface:

```rust
pub mod adapter;    // hardness_of, InputScalar, ModelScalars::from_summary
pub mod logs;       // petro_samples(&[WellCurveInput]) -> Vec<(f64,f64,f64)>  (tvd, φ, Sw)
pub mod wireframe;  // assemble_wireframe(&ModelInputs) -> Result<Wireframe>
pub mod petekio;    // re-export of the real petekio crate (the upstream seam)
```

The integrated product composes petekIO and petekStatic in petekSim. New core
callers should construct petekStatic-owned wireframe/model inputs directly.

## `srs-volumetrics` — GRV, in-place, FVF (relocated from petekSim 2026-07-03)

```rust
pub const PORO: &str;  pub const NTG: &str;  pub const SW: &str;   // canonical cube names

pub struct ConstantPriors { pub porosity: f64, pub net_to_gross: f64, pub water_saturation: f64 }
pub fn populate_constant(grid: &mut Grid, priors: ConstantPriors) -> Result<()>;

pub struct ZoneVolumes { pub grv_m3: f64, pub hcpv_m3: f64, pub cells: usize }  // one HC zone
pub struct InPlace {
    pub grv_m3: f64, pub hcpv_m3: f64,             // whole HC column (gas+oil)
    pub cells_in_column: usize, pub per_cell_hcpv: Vec<f64>,
    pub gas: Option<ZoneVolumes>,                    // Some only for a two-contact column
    pub oil: Option<ZoneVolumes>,
}
impl InPlace {                                        // volumes m³ internal; SI reporting scales
    pub fn grv_mcm(&self) -> f64;                     // GRV in mcm (1e6 m³)
    pub fn ooip_sm3(&self, boi: OilFvf) -> f64;       // whole column, Sm³
    pub fn oil_msm3(&self, boi: OilFvf) -> f64;       // whole column, MSm³ (oil report scale)
    pub fn ogip_sm3(&self, bgi: GasFvf) -> f64;       // whole column, Sm³
    pub fn gas_bcm(&self, bgi: GasFvf) -> f64;        // whole column, bcm (1e9 Sm³, gas report scale)
    pub fn gas_zone_ogip_sm3(&self, bgi: GasFvf) -> f64;   // gas cap only, Sm³ (0 if none)
    pub fn oil_zone_ooip_sm3(&self, boi: OilFvf) -> f64;   // oil leg only, Sm³ (0 if none)
}
pub fn compute_in_place(grid: &Grid, contact_depth_m: f64) -> Result<InPlace>;   // single contact
pub fn compute_in_place_summary(grid: &Grid, contact_depth_m: f64) -> Result<InPlace>;  // V7: same aggregates, no per-cell HCPV cube
// Two-contact (gas cap + oil rim): partitions cells gas (z<GOC) / oil
// (GOC<=z<OWC) / water. Geometry + per-zone in-place ONLY — no PVT (solution
// gas / gas-cap expansion stay in petekSim). Err if GOC deeper than OWC.
pub fn compute_in_place_two_contact(grid: &Grid, goc_m: f64, owc_m: f64, sw_gas: Option<f64>) -> Result<InPlace>;  // sw_gas: gas-cap connate-water override (R3)
pub fn compute_in_place_two_contact_summary(grid: &Grid, goc_m: f64, owc_m: f64, sw_gas: Option<f64>) -> Result<InPlace>;  // V7 summary

// FVF value types — FVF crosses the layer seam as a validated scalar INPUT
// (no PVT code); petekSim's srs-pvt keeps its own copies for the dynamic work.
// Boi/Bgi are dimensionless Rm³/Sm³ (== legacy rb/STB, rcf/scf numerically — a relabel).
pub struct OilFvf;  impl OilFvf { pub fn new(rm3_per_sm3: f64) -> Result<OilFvf>; pub fn value(self) -> f64; }   // finite, >= 1
pub struct GasFvf;  impl GasFvf { pub fn new(rm3_per_sm3: f64) -> Result<GasFvf>; pub fn value(self) -> f64; }  // finite, in (0,1)

// Physical-range predicates (H2): typed errors, never silent garbage.
pub fn validate_fraction(what: &str, x: f64) -> Result<()>;   // finite, in [0,1]
pub fn validate_positive(what: &str, x: f64) -> Result<()>;   // finite, > 0
```

## `srs-uncertainty` — the Monte Carlo toolkit (relocated from petekSim 2026-07-03)

```rust
pub enum Distribution {
    Constant(f64),
    Uniform { min: f64, max: f64 },
    Triangular { min: f64, mode: f64, max: f64 },
    Normal { mean: f64, sd: f64 },
    Lognormal { mu: f64, sigma: f64 },
}
impl Distribution {
    pub fn uniform(min: f64, max: f64) -> Result<Distribution>;          // validated ctors
    pub fn triangular(min: f64, mode: f64, max: f64) -> Result<Distribution>;
    pub fn normal(mean: f64, sd: f64) -> Result<Distribution>;
    pub fn lognormal(mu: f64, sigma: f64) -> Result<Distribution>;
    pub fn quantile(self, u: f64) -> f64;                                // inverse CDF
    pub fn sample(self, rng: &mut SplitMix64) -> f64;
}

pub struct SplitMix64;  impl SplitMix64 { pub fn new(seed: u64) -> SplitMix64; pub fn next_u64(&mut self) -> u64; pub fn next_f64(&mut self) -> f64; }
pub fn inverse_normal_cdf(p: f64) -> f64;                                // Acklam probit

pub struct Realizations { pub values: Vec<f64> }
impl Realizations { pub fn summary(&self) -> Result<PercentileSummary>; }   // Err on empty (H1)
pub fn run<F: FnMut(&mut SplitMix64) -> f64>(n: usize, seed: u64, trial: F) -> Realizations;

pub struct PercentileSummary { pub p90: f64, pub p50: f64, pub p10: f64, pub mean: f64 }
impl PercentileSummary {
    pub fn from_realizations(values: &[f64]) -> Result<PercentileSummary>;  // Err on empty (H1)
    pub fn swanson_mean(&self) -> f64;   // delegates the P90/P50/P10 digest + type-7 percentile to petektools::sampling::reservoir_summary
}
```

## `srs-model` — the `StaticModel` + the regeneration seam (LANDED 2026-07-03)

The ratified headline contract (`decision_staticmodel_regen_seam` + amendments;
`decision_layer_charters`). Re-exports the volumetrics/P-curve output surface:
`ConstantPriors`, `InPlace`, `ZoneVolumes`, `OilFvf`, `GasFvf`,
`PercentileSummary`.

### The aggregate

```rust
/// A populated static reservoir model. Owned, Clone, Send (compile-checked);
/// construction guarantees the invariants (SPEC §2) — consumers never re-validate.
pub struct StaticModel { /* framework, grid, zones, provenance — all private */ }
impl StaticModel {
    // --- read-only accessors (SPEC §7c) ---
    pub fn grid(&self) -> &Grid;
    pub fn framework(&self) -> &Wireframe;
    pub fn contacts(&self) -> &[Contact];
    pub fn property(&self, name: &str) -> Option<&Property>;
    pub fn property_names(&self) -> Vec<&str>;
    pub fn zones(&self) -> &ZoneTable;
    pub fn provenance(&self) -> &Provenance;
    pub fn georef(&self) -> Option<Georef>;              // registered WORLD frame (None = local degenerate: synthetic square/box)
    pub fn bulk_volume(&self) -> f64;                    // pure geometry (whole grid)
    // --- the volumetrics output surface (the model owns volumes) ---
    // Clipped vs the model's contact(s). A GOC + lower (OWC/FWL/GWC) contact
    // auto-returns the gas-cap + oil-rim split (InPlace::gas / ::oil); a lone
    // contact stays a generic column.
    pub fn in_place(&self) -> Result<InPlace>;
    pub fn in_place_summary(&self) -> Result<InPlace>;  // V7: same aggregates, per_cell_hcpv empty (MC hot path)
    // --- multi-zone surface (SPEC §4a; from_horizon_stack) ---
    pub fn in_place_by_zone(&self) -> Result<ZonedInPlace>;   // per-zone in-place vs EACH zone's contacts + rollup total (sum == total)
    pub fn zone_stats(&self, property: &str) -> Result<Vec<ZoneStat>>;  // per-zone count/mean/min/max over active cells
}
pub struct ZoneInPlace { pub zone: String, pub in_place: InPlace }
pub struct ZonedInPlace { pub zones: Vec<ZoneInPlace>, pub total: InPlace }  // total = summary rollup (per_cell_hcpv empty)
pub struct ZoneStat { pub zone: String, pub count: usize, pub mean: f64, pub min: f64, pub max: f64 }  // NaN aggregates when count==0
```

### Zones (SPEC §4)

```rust
pub struct Zone { pub name: String, pub color: Option<String>, pub top_horizon: String, pub base_horizon: String, pub k_range: std::ops::Range<usize> }  // color: viewer hint from StackZone.color
pub struct ZoneTable { /* ordered top→base; k-ranges partition [0, nk) */ }
impl ZoneTable {
    pub fn single(nk: usize) -> ZoneTable;                                          // whole-column degenerate zone
    pub fn from_stack(horizon_names: &[String], zone_names: &[String], zone_colors: &[Option<String>], per_zone_nk: &[usize]) -> ZoneTable;  // real multi-zone
    pub fn zones(&self) -> &[Zone];
    pub fn get(&self, name: &str) -> Option<&Zone>;
}
```

### Multi-zone horizon stack (SPEC §4a; `task_petekstatic_multizone`)

```rust
pub struct Pick { pub ip: usize, pub jp: usize, pub depth_m: f64 }                  // a well top on the node lattice
pub struct WorldPoint { pub x: f64, pub y: f64, pub depth_m: f64 }                  // a raw scatter obs in WORLD coords (positive-down depth)
pub enum HorizonSource {
    Scatter(Vec<WorldPoint>),                    // RAW world-coord scatter — the engine grids it (bilinear-conditioned, voids left NaN); only via from_scatter_stack
    Mapped(GriddedDepth),                        // a mapped surface ALREADY on the lattice — bypasses the engine's solve/conditioning fidelity (loaded grids only; raw points belong in Scatter)
    TopsOnly(Vec<Pick>),                         // no surface — draped conformally from the nearest mapped horizon above at pick-controlled thickness
}
pub struct StackHorizon { pub name: String, pub source: HorizonSource }
pub struct StackZone { pub name: String, pub color: Option<String>, pub conformity: Conformity, pub nk: usize, pub contacts: Vec<Contact> }  // name+optional colour (folds the old zone_names) + per-zone layering + its OWN contacts (empty = contactless)
impl StackZone {
    pub fn new(name: impl Into<String>, conformity: Conformity, nk: usize, contacts: Vec<Contact>) -> StackZone;  // colour = None
    pub fn with_color(self, color: impl Into<String>) -> StackZone;                  // viewer hint -> section bundle zones[{name,color}]
}
pub struct HorizonStack {                        // N horizons top→down -> N-1 zones
    pub horizons: Vec<StackHorizon>,             // first (top) must be Mapped (or Scatter, via from_scatter_stack)
    pub zone_layers: Vec<StackZone>,             // N-1 (each carries its own name/colour/layering/contacts)
}
// The areal lattice + world frame a Scatter stack conditions onto (from_scatter_stack).
pub struct StackFrame { pub ni: usize, pub nj: usize, pub georef: Georef }          // ni = √area_m2 / georef.spacing_x
// P8 per-horizon well tie: world (x,y) marker on control node (ip,jp) + measured tops per horizon.
pub struct WellTie { pub id: String, pub x: f64, pub y: f64, pub ip: usize, pub jp: usize, pub tops: Vec<(String, f64)> }
impl WellTie {
    pub fn new(id: impl Into<String>, x: f64, y: f64, ip: usize, jp: usize) -> WellTie;
    pub fn with_top(self, horizon: impl Into<String>, depth_m: f64) -> Self;   // measured formation top (positive-down)
}
// Builder entry: StaticModelBuilder::from_horizon_stack(HorizonStack, BuildOpts) — see below.

// srs-gridder: the stacked layering primitive.
pub struct ZoneLayerSpec { pub conformity: Conformity, pub requested_nk: usize }
pub struct StackedZone { pub nk: usize, pub k_start: usize, pub conformity: Conformity, pub truncated_cells: usize, pub collapsed_cells: usize }
pub struct StackedLayeredGrid { pub grid: Grid, pub nk: usize, pub zones: Vec<StackedZone>, pub truncated_cells: usize, pub collapsed_cells: usize, pub nk_capped: bool }
pub fn layer_grid_stack(surfaces: &[&Surface], dx: f64, dy: f64, zone_specs: &[ZoneLayerSpec], collapse_below_m: Option<f64>) -> Result<StackedLayeredGrid>;
// Buffer-recycling core of layer_grid_stack (the realize_into hot path): fills the caller's
// COORD + ZCORN buffers + reusable scratch in place (all cleared + fully overwritten).
pub struct LayerScratch { /* private raw/snap columns; ::new() */ }
pub struct StackLayering { pub dims: Dims, pub nk: usize, pub zones: Vec<StackedZone>, pub truncated_cells: usize, pub collapsed_cells: usize, pub nk_capped: bool }
pub fn layer_grid_stack_into(surfaces: &[&Surface], dx: f64, dy: f64, zone_specs: &[ZoneLayerSpec], collapse_below_m: Option<f64>, scratch: &mut LayerScratch, coord: &mut Vec<Pillar>, zcorn: &mut Vec<f64>) -> Result<StackLayering>;
// srs-grid recycling primitives: Grid::{take_geometry_buffers, install_geometry}, Properties::take_values, CornerPointGeom::{take_buffers, install}.
```

### Provenance (SPEC §6)

```rust
pub enum PopulationMode { Priors, Logs }
pub enum BuildWarning {                          // non-blocking build advisory
    UnusedHorizon { name: String, role: HorizonRole, reason: String },
    ThinColumnsRepaired { columns: usize, worst_m: f64 }, // R-c: with_min_thickness_m pulled thin/crossing base columns to the floor
    LayersTruncated { cells: usize },            // Follow conformity collapsed `cells` deep/shallow layers to zero thickness (informational; volume-invariant)
    CellsCollapsed { cells: usize },             // with_collapse_below_m: sub-threshold cells merged into a thicker zone-interior neighbour (volume-conserving)
    LayerCountCapped { nk: usize },              // dz-derived nk hit MAX_NK — thickest columns' deepest part not layered (coarsen dz or accept)
}
pub struct Provenance {
    pub inputs_ref: String,
    pub solve_opts: SolveOpts,
    pub conformity: Conformity,
    pub nk: usize,                              // EFFECTIVE nk (dz-derived under a Follow conformity; total across zones on the stack path)
    pub population: PopulationMode,
    pub realization: Option<RealizationDraw>,   // Some(_) for an MC realization
    pub warnings: Vec<BuildWarning>,            // e.g. supplied horizons unused by the build
    pub property_reports: Vec<PropertyReport>,  // P5: one per attached PropertyPipeline (upscale QC + propagated)
    pub stack: Option<StackProvenance>,         // the multi-zone stack record (from_horizon_stack); None on the Top+Base path
    pub well_ties: Vec<WellTieRecord>,          // P8: per-horizon tie residuals (with_well_ties); empty otherwise
}
// Multi-zone build record (SPEC §4a):
pub struct ZoneProvenance { pub name: String, pub top_horizon: String, pub base_horizon: String, pub conformity: Conformity, pub nk: usize, pub k_start: usize, pub truncated_cells: usize }
pub struct InterfaceRepair { pub interface: usize, pub columns: usize, pub worst_m: f64 }
// P8 per-horizon well ties (residual = measured_depth_m − model_depth_m at the well node, pre-tie):
pub struct HorizonTieResidual { pub horizon: String, pub measured_depth_m: f64, pub model_depth_m: f64, pub residual_m: f64 }
pub struct WellTieRecord { pub id: String, pub x: f64, pub y: f64, pub ip: usize, pub jp: usize, pub residuals: Vec<HorizonTieResidual> }
pub struct StackProvenance { pub horizons: Vec<String>, pub zones: Vec<ZoneProvenance>, pub interface_repairs: Vec<InterfaceRepair> }
```

### Per-property geostatistical pipeline (SPEC §3c; P5)

```rust
// A positioned well log for one property: world (x, y) + (tvd_m, value) samples.
pub struct WellLog { pub x: f64, pub y: f64, pub samples: Vec<(f64, f64)> }
impl WellLog { pub fn new(x: f64, y: f64, samples: Vec<(f64, f64)>) -> WellLog; }

pub enum UpscaleMethod { Arithmetic, Harmonic, Geometric }   // in-cell log average

// Monte-Carlo behaviour across realizations (decision_mc_composition).
// LevelShift adds the drawn shift to the cached pattern per cell; the fraction
// cubes PORO/NTG/SW are shift-then-clamped to [0,1] so boundary cells (NTG=0,
// SW=1) saturate instead of escaping range and tripping the per-cell H2 check
// (F9). H2 still validates the drawn INPUTS at the seam.
pub enum McMode { LevelShift /* default */, Resimulate }

// The SGS propagation spec (`ps.gaussian(...)`): data-space variogram + seed +
// optional moving-neighbourhood search + optional collocated trend cokriging.
pub struct Gaussian { /* private */ }
impl Gaussian {
    pub fn new(variogram: petektools::Variogram, seed: u64) -> Gaussian;
    // DEFAULT search is BOUNDED: 16 nodes within max(1.5*variogram.range, 4*node_spacing).
    // (Old whole-grid default -> >15 min/cube on real lattices; bounded matches it within
    // simulation tolerance since beyond-range kriging weight is ~0.)
    pub fn with_search(self, max_neighbours: usize, radius_m: f64) -> Self;
    pub fn with_unbounded_search(self) -> Self;   // opt back into the legacy whole-grid window
    pub fn allow_mean_fill(self) -> Self;         // opt into filling a data-less simulated layer with the
                                                  // conditioned mean; DEFAULT is a hard InvalidInput naming
                                                  // the property (+ zone, via the caller) — no silent constant fill
    pub fn with_trend(self, trend: TrendSurface, corr: f64) -> Self;   // collocated (Markov-1); corr=0 == plain SGS
    // A world-georeferenced trend is resampled at each column's WORLD position via the
    // MODEL georef (build with with_georef); on a local model it stays local/index-space.
    // A secondary covering < 50% of the model frame (a georef mismatch) is a hard
    // InvalidInput, NOT a silent per-node fallback to plain SGS.
}

// Model a property ONE AT A TIME: upscale logs -> per-cell conditioned values,
// then SGS-propagate (per k-layer) conditioned on those cells.
pub struct PropertyPipeline { /* private */ }
impl PropertyPipeline {
    pub fn new(name: impl Into<String>) -> PropertyPipeline;
    pub fn upscale(self, wells: Vec<WellLog>, method: UpscaleMethod) -> Self;
    pub fn propagate(self, gaussian: Gaussian) -> Self;
    pub fn name(&self) -> &str;
    // The visible upscale step: per-cell conditioned field (NaN where no log) + QC.
    // Each sample is binned against the cell interval interpolated AT THE WELL's (x,y)
    // (bilinear), not the column-centroid mean — so an off-centroid well on a dipping
    // zone boundary assigns its samples to the zone their depth truly falls in.
    pub fn upscale_cells(&self, grid: &Grid) -> Result<(Vec<f64>, UpscaleQc)>;
    // Full run against a grid (upscale + propagate), setting the cube.
    pub fn apply(&self, grid: &mut Grid) -> Result<PropertyReport>;
    // P8: run restricted to one zone's k-range, MERGING into the cube already on the
    // grid (other zones' slices untouched). apply_in_zone(0..nk) == apply().
    pub fn apply_in_zone(&self, grid: &mut Grid, k_range: Range<usize>) -> Result<PropertyReport>;
}

pub struct UpscaleQc {                              // upscaled-vs-log QC (Clone, Debug, PartialEq)
    pub property: String,
    pub conditioned_cells: usize, pub log_samples: usize,
    pub log_mean: f64, pub upscaled_mean: f64, pub upscaled_min: f64, pub upscaled_max: f64,
}
pub struct PropertyReport { pub property: String, pub upscale: UpscaleQc, pub propagated: bool }
```

> Regular / axis-aligned column lattices (the `layer_grid` box/conformable grids);
> rotated pillars are future work. Supersedes the interim `TrendSurface` multiplier
> hook below (which now only serves the deterministic `with_areal_trend` path and,
> via the shared `petektools::resample` kernel, the collocated-secondary resample).

### Areal trend (external-drift-lite; INTERIM)

> **INTERIM mechanism** (coordinator design `python-model-build-api.md`): this
> multiplier hook is the stop-gap trend conditioning until collocated cokriging
> inside per-property propagation lands (P5
> `task_petekstatic_property_modelling`), which supersedes it. Expect the
> surface to be re-ratified then.

```rust
// A gridded areal multiplier field: nearest-node resampled to the model column
// lattice, mean-normalized (field-mean preserved), applied per-column to NTG
// (and PORO if flagged) after population — lateral SHAPE only. NaN node -> 1.0.
pub struct TrendSurface { /* ncol, nrow, values, apply_porosity, georef — private */ }
impl TrendSurface {
    pub fn new(ncol: usize, nrow: usize, values: Vec<f64>) -> Result<TrendSurface>;  // finite >= 0; NaN ok
    pub fn with_porosity(self) -> Self;              // also modulate PORO (NTG always)
    pub fn with_georef(self, origin_x: f64, origin_y: f64, node_dx: f64, node_dy: f64) -> Self;  // R4: resample by world (x,y), not index fraction
    pub fn with_oriented_georef(self, origin_x: f64, origin_y: f64, node_dx: f64, node_dy: f64, rotation_deg: f64, yflip: bool) -> Self;
    pub fn is_georeferenced(&self) -> bool;
    pub fn applies_to_porosity(&self) -> bool;
}
```

### Construction — deterministic builder + the regeneration template

```rust
pub type PetroSample = (f64, f64, f64);   // (tvd_m, φ, Sw) — population input

pub struct BuildOpts {
    pub area_m2: f64, pub gross_height_m: f64,
    pub nk: usize, pub conformity: Conformity, pub solve_opts: SolveOpts, // nk honoured by Proportional only; Follow styles derive it from conformity dz_m
    pub priors: ConstantPriors,
}

// THE ONE declarative build configuration (task_petekstatic_spec_mirror), consumed
// by BOTH StaticModelBuilder and StaticModelTemplate: every with_* setter on either
// is thin sugar mutating an internal BuildSpec; with_spec installs a whole one.
// #[non_exhaustive] + #[serde(default)] (additive forward-compat — the per-zone
// isochore structural leg lands as a new optional field). Serde round-trips as a
// scenario file; compares by value. The WHOLE config layer derives serde
// (HorizonStack family incl. Scatter/WorldPoint, WellTie, Pick, BuildOpts, Georef,
// RealizationDraw, StructuralPerturbation, ZoneDraw, BuildSpec, TieSettings,
// McSettings + the srs-wireframe/-gridder/-volumetrics value types they carry).
// McInputs is EXCLUDED until petekTools' Sampler/Clamped derive serde (queued).
#[non_exhaustive]
pub struct BuildSpec {
    pub inputs_ref: Option<String>,              // None = the entry default ("flat-box"/"wireframe"/"horizon-stack")
    pub georef: Option<Georef>, pub boundary: Option<Vec<[f64; 2]>>,
    pub extrapolation: ExtrapolationPolicy,
    pub clamp_base_to_top: bool, pub min_thickness_m: Option<f64>, pub collapse_below_m: Option<f64>,
    pub sugar_cube: bool, pub sw_gas: Option<f64>,
    pub well_ties: Vec<WellTie>, pub ties: TieSettings,
}
impl BuildSpec {   // Default = every knob at the historical default; with_* sugar mirrors the consumers'
    pub fn new() -> BuildSpec;
    pub fn with_inputs_ref(self, r: impl Into<String>) -> Self;  pub fn with_georef(self, ox: f64, oy: f64, sx: f64, sy: f64) -> Self;
    pub fn with_oriented_georef(self, ox: f64, oy: f64, sx: f64, sy: f64, rotation_deg: f64, yflip: bool) -> Self;
    pub fn with_boundary(self, ring: Vec<[f64; 2]>) -> Self;     pub fn with_extrapolation(self, p: ExtrapolationPolicy) -> Self;
    pub fn with_clamp_base_to_top(self, c: bool) -> Self;        pub fn with_min_thickness_m(self, m: f64) -> Self;
    pub fn with_collapse_below_m(self, m: f64) -> Self;          pub fn with_sugar_cube(self, s: bool) -> Self;
    pub fn with_sw_gas(self, s: f64) -> Self;                    pub fn with_well_ties(self, t: Vec<WellTie>) -> Self;
    pub fn with_tie_settings(self, t: TieSettings) -> Self;
}

// How well ties fold into the mapped datum grids (the settings mirror over the
// datum-substitution tie machinery). Default = Replace (today's behaviour).
pub struct TieSettings { pub method: TieMethod }
impl TieSettings { pub fn replace() -> Self;  pub fn radius(radius_m: f64) -> Self; }
pub enum TieMethod {
    Replace,                       // measured top REPLACES the map datum at the tie node; radius 0 on a dense lattice
    Radius { radius_m: f64 },      // bounded locality: residual decays linearly 1→0 over radius_m across DEFINED datums; beyond = bit-untouched; voids stay the solve's to taper; finite>0 or typed error
}

/// Deterministic single-shot build — the relocated RefiningModel pipeline.
pub struct StaticModelBuilder { /* ... */ }
impl StaticModelBuilder {
    pub fn flat(ni: usize, nj: usize, top_depth_m: f64, contact_depth_m: f64, opts: BuildOpts) -> Result<StaticModelBuilder>;
    pub fn from_wireframe(wf: &Wireframe, opts: BuildOpts) -> Result<StaticModelBuilder>;  // wires a Base horizon's real relief; else gross_height_m offset
    pub fn from_horizon_stack(stack: HorizonStack, opts: BuildOpts) -> Result<StaticModelBuilder>;  // SPEC §4a: N horizons -> N-1 zones (per-zone conformity + contacts); opts.nk/conformity/gross_height_m unused here
    pub fn from_scatter_stack(stack: HorizonStack, opts: BuildOpts, frame: StackFrame) -> Result<StaticModelBuilder>;  // SPEC §4a: raw-scatter horizons — the engine conditions them onto `frame` (bilinear, voids NaN) then resolves as from_horizon_stack; registers frame.georef. THE single scatter-gridding authority
    pub fn condition_scatter_stack(stack: HorizonStack, frame: &StackFrame) -> Result<HorizonStack>;  // dedup seam (task_suite_scatter_perf): condition raw scatter ONCE -> the conditioned all-Mapped handle; feed to BOTH from_horizon_stack(+with_georef) and StaticModelTemplate::from_horizon_stack(+with_georef) to build a model + its MC template without re-solving the cold bilinear conditioning. Bit-identical to conditioning inside from_scatter_stack
    pub fn with_collapse_below_m(self, collapse_below_m: f64) -> Self;        // Petrel-style cell-collapse: sub-threshold cells merge volume-conservingly into a thicker zone-interior neighbour; warns CellsCollapsed
    pub fn with_logs(self, samples: Vec<PetroSample>) -> Self;
    pub fn with_areal_trend(self, trend: TrendSurface) -> Self;               // INTERIM lateral NTG (+φ) shape
    pub fn with_property(self, pipeline: PropertyPipeline) -> Self;           // P5: per-property upscale + SGS, overwriting the cube (-> Provenance.property_reports)
    pub fn with_zone_property(self, zone: impl Into<String>, pipeline: PropertyPipeline) -> Self;  // P8: per-zone upscale+SGS restricted to the zone's k-range (own variogram/trend/logs); unknown zone = build-time InvalidInput
    pub fn with_zone_priors(self, zone: impl Into<String>, priors: ConstantPriors) -> Self;        // P8: override PORO/NTG/SW across one zone's k-range (per-zone level); unknown zone = build-time InvalidInput
    pub fn with_well_ties(self, ties: Vec<WellTie>) -> Self;                  // P8: per-horizon ties per the spec's TieSettings (Replace default / Radius bounded locality) + re-run the stack resolution; pre-tie residuals -> Provenance.well_ties + map wells[].ties (from_horizon_stack)
    pub fn with_tie_settings(self, ties: TieSettings) -> Self;                 // how the ties fold in; read (with the ties) at build() time — order vs with_well_ties irrelevant here
    pub fn with_clamp_base_to_top(self, clamp: bool) -> Self;                 // R1: clamp a crossed base to zero gross instead of erroring (default: error)
    pub fn with_min_thickness_m(self, min_thickness_m: f64) -> Self;          // R-c: post-gridding repair — pull thin/crossing base to top+min (top preserved); warns ThinColumnsRepaired; precedence over clamp
    pub fn with_extrapolation(self, policy: ExtrapolationPolicy) -> Self;      // stack path: beyond-data behaviour; default DecayToData (never silent unbounded natural-dip)
    pub fn with_sw_gas(self, sw_gas: f64) -> Self;                            // R3: gas-cap connate water for a two-contact in_place split
    pub fn with_inputs_ref(self, inputs_ref: impl Into<String>) -> Self;
    pub fn with_georef(self, origin_x: f64, origin_y: f64, spacing_x: f64, spacing_y: f64) -> Self;  // register the model's WORLD frame (column (0,0) centroid + spacing) so view bundles emit ONE world frame; non-positive spacing ignored; grid geometry untouched
    pub fn with_oriented_georef(self, origin_x: f64, origin_y: f64, spacing_x: f64, spacing_y: f64, rotation_deg: f64, yflip: bool) -> Self; // CCW east→+I; grid kernels remain local
    pub fn with_boundary(self, ring: Vec<[f64; 2]>) -> Self;                  // from_horizon_stack: world outline ring for the map bundle; omitted -> a world-extent rectangle from the georef (never the unit square vs a world frame)
    pub fn with_sugar_cube(self, sugar_cube: bool) -> Self;                   // section rendering: false (default) = dip-following trapezoids; true = flat boxes -> IntersectionBundle.sugar_cube + flattened edge arrays
    pub fn with_spec(self, spec: BuildSpec) -> Self;                          // install the WHOLE declarative config in one call (the with_* above are the same sugar); bit-identical results, pinned by tests/spec_conformance.rs
    pub fn add_top_control(&mut self, ip: usize, jp: usize, depth_m: f64);   // live refine
    pub fn control_count(&self) -> usize;
    pub fn build(&self) -> Result<StaticModel>;          // warm kernel space (unified with the template, R2); base-above-top → CrossedSurfaces
}

/// The MC-regeneration template (SPEC §7a): built once, holds the reusable
/// warm-start state IN petekTools KERNEL SPACE. Send; one per worker.
pub struct StaticModelTemplate { /* ... */ }
impl StaticModelTemplate {
    // A Base horizon in `wf` is extracted as the nominal gross field g(x,y)
    // (decision_template_gross_scaling): per draw, gross = g * gross_height_m
    // / mean(g) — draw sets the LEVEL, g the SHAPE; a draw at mean(g)
    // reproduces the deterministic build on ALL column shapes (R2: builder +
    // template share one warm kernel path). No Base = constant offset (as before).
    pub fn new(wf: &Wireframe, opts: BuildOpts) -> Result<StaticModelTemplate>;
    // P8 stack-aware MC (task_petekstatic_multizone_2): resolves the multi-horizon
    // framework ONCE (surfaces are draw-invariant); realize varies per draw only the
    // areal footprint (spacing), the per-zone contacts, and the per-zone property
    // levels (RealizationDraw.zones) — bit-deterministic, realize_into-recyclable.
    pub fn from_horizon_stack(stack: HorizonStack, opts: BuildOpts) -> Result<StaticModelTemplate>;
    pub fn from_scatter_stack(stack: HorizonStack, opts: BuildOpts, frame: StackFrame) -> Result<StaticModelTemplate>;  // raw-scatter MC template — conditions + resolves scatter byte-identically to StaticModelBuilder::from_scatter_stack
    pub fn with_logs(self, samples: Vec<PetroSample>) -> Self;
    pub fn with_areal_trend(self, trend: TrendSurface) -> Self;   // INTERIM per-realization lateral shape
    // P5 property pipelines per realization (decision_mc_composition):
    pub fn with_property(self, pipeline: PropertyPipeline) -> Self;               // default McMode::LevelShift
    pub fn with_property_mode(self, pipeline: PropertyPipeline, mode: McMode) -> Self;
    // Zone-scoped pipelines (stack templates): realize honours the with_zone_property
    // upscale+SGS cube over the per-zone priors (level shift on top) — a zero-spread
    // zoned MC reproduces in_place_by_zone on every piped zone. Unknown zone = realize InvalidInput.
    pub fn with_zone_property(self, zone: impl Into<String>, pipeline: PropertyPipeline) -> Self;             // default McMode::LevelShift
    pub fn with_zone_property_mode(self, zone: impl Into<String>, pipeline: PropertyPipeline, mode: McMode) -> Self;
    // Well ties (stack templates): draw-invariant, applied ONCE at construction
    // (re-solve + repair); every draw inherits tied geometry at zero per-draw cost.
    // Err if not stack-aware / bad tie node / unknown horizon.
    pub fn with_well_ties(self, ties: Vec<WellTie>) -> Result<StaticModelTemplate>;
    pub fn with_tie_settings(self, ties: TieSettings) -> Self;    // set BEFORE with_well_ties (ties apply at that call); same tie authority as the builder — node-for-node identical
    pub fn with_sugar_cube(self, sugar_cube: bool) -> Self;       // section rendering (see StaticModelBuilder::with_sugar_cube)
    pub fn with_clamp_base_to_top(self, clamp: bool) -> Self;     // R1: clamp a crossed base per realization instead of erroring
    pub fn with_min_thickness_m(self, min_thickness_m: f64) -> Self;  // R-c: per-realization post-gridding repair to top+min (top preserved); warns ThinColumnsRepaired
    pub fn with_extrapolation(self, policy: ExtrapolationPolicy) -> Result<StaticModelTemplate>;  // stack templates: re-resolve surfaces under the policy (default DecayToData)
    pub fn with_collapse_below_m(self, collapse_below_m: f64) -> Self;  // per-realization cell-collapse (volume-conserving); warns CellsCollapsed
    pub fn with_inputs_ref(self, inputs_ref: impl Into<String>) -> Self;
    pub fn with_georef(self, origin_x: f64, origin_y: f64, spacing_x: f64, spacing_y: f64) -> Self;  // stamp the WORLD frame onto every realized model (see StaticModelBuilder::with_georef)
    pub fn with_oriented_georef(self, origin_x: f64, origin_y: f64, spacing_x: f64, spacing_y: f64, rotation_deg: f64, yflip: bool) -> Self;
    pub fn with_boundary(self, ring: Vec<[f64; 2]>) -> Self;  // stack-aware: world outline ring stamped onto each realized map bundle (see StaticModelBuilder::with_boundary)
    pub fn with_spec(self, spec: BuildSpec) -> Result<StaticModelTemplate>;   // install the WHOLE declarative config (the SAME BuildSpec the builder consumes); re-applies extrapolation + ties through the setters above — bit-identical, pinned
    /// The CHEAP per-realization call: warm-started solve, re-layer, re-populate.
    /// `&mut self` advances the warm chain; Result per draw.
    pub fn realize(&mut self, draw: &RealizationDraw) -> Result<StaticModel>;  // = new-empty-then-realize_into (allocating convenience)
    /// The buffer-recycling hot path (ratified amendment 3): realize into a REUSED
    /// model, overwriting its ZCORN + cube allocations in place (~100 MB/draw at 1M
    /// cells). Bit-identical to `realize` on the same draw + chain. Prefer on the MC loop.
    pub fn realize_into(&mut self, draw: &RealizationDraw, model: &mut StaticModel) -> Result<()>;
    /// A fresh empty model to drive a `realize_into` loop against — one per MC worker.
    pub fn reusable_model(&self) -> StaticModel;
}

/// The per-realization sampled input set — petekStatic's neutral type, filled by
/// the sampler. #[non_exhaustive]: construct via ::new + with_*. fvf EXCLUDED
/// (applied at the volumetrics output surface, never rides the draw).
#[non_exhaustive]
pub struct RealizationDraw {                              // Clone + Debug
    pub area_m2: f64, pub gross_height_m: f64, pub contact_depth_m: f64,  // contact = OWC/FWL (lower) when a GOC is set; gross = level multiplier over the template's base-relief shape when a Base horizon exists
    pub goc_depth_m: Option<f64>,                       // Some(_) => two-contact (gas cap + oil rim)
    pub porosity: f64, pub net_to_gross: f64, pub water_saturation: f64,
    pub seed_index: u64,
    pub structural: Option<StructuralPerturbation>,     // MVP explicit control-node depth shifts (2-surface path)
    pub top_structural: Option<PerturbationField>,      // correlated TOP-surface DEPTH field this draw (both paths); decision_structural_uncertainty_isochore
    pub sw_gas: Option<f64>,                             // R3: gas-cap connate-water override (with a GOC)
    pub property_shifts: Vec<(String, f64)>,            // P5: per-property additive level shift for McMode::LevelShift
    pub zones: Vec<ZoneDraw>,                           // P8: per-zone contacts + property levels (stack-aware from_horizon_stack MC); empty on the 2-surface path
}
// P8 per-zone draw for a stack-aware template. A ZoneDraw's contacts REPLACE the
// template's static contacts for that zone (neither GOC nor OWC => explicitly
// contactless: GRV, zero HC). porosity/ntg/sw override the base priors in the zone.
#[non_exhaustive]                                        // additive forward-compat (per-zone isochore leg)
pub struct ZoneDraw {                                    // Clone + Debug + PartialEq + serde
    pub zone: usize,                                     // zone index (top->down)
    pub goc_depth_m: Option<f64>, pub owc_depth_m: Option<f64>,
    pub porosity: Option<f64>, pub net_to_gross: Option<f64>, pub water_saturation: Option<f64>,
    pub isochore_structural: Option<PerturbationField>,  // correlated THICKNESS field for this zone (moves its DEEPER horizon); clamped >=0, zero-masked at exact merges; decision_structural_uncertainty_isochore
}
impl ZoneDraw {
    pub fn new(zone: usize) -> ZoneDraw;                 // no contacts / overrides
    pub fn with_owc(self, owc_depth_m: f64) -> Self;  pub fn with_goc(self, goc_depth_m: f64) -> Self;
    pub fn with_priors(self, porosity: f64, net_to_gross: f64, water_saturation: f64) -> Self;
    pub fn with_isochore_structural(self, field: PerturbationField) -> Self;   // attach this zone's isochore perturbation field
}
impl RealizationDraw {
    pub fn new(area_m2: f64, gross_height_m: f64, contact_depth_m: f64,
               porosity: f64, net_to_gross: f64, water_saturation: f64,
               seed_index: u64) -> RealizationDraw;   // goc None, structural None
    pub fn with_zone_draw(self, zone: ZoneDraw) -> Self;   // P8: attach/replace a per-zone draw (stack-aware MC)
    pub fn with_structural(self, p: StructuralPerturbation) -> Self;
    pub fn with_top_structural(self, field: PerturbationField) -> Self;   // attach the correlated TOP-surface depth perturbation field
    pub fn with_property_shift(self, property: impl Into<String>, delta: f64) -> Self;  // P5 level shift
    pub fn property_shift(&self, property: &str) -> f64;                                // 0.0 if none
    pub fn with_goc(self, goc_depth_m: f64) -> Self;    // gas cap above; must be shallower than the OWC
    pub fn with_sw_gas(self, sw_gas: f64) -> Self;       // R3: gas-cap connate water (effective only with a GOC)
    // named scalar setters (V8/R6 ergonomics) — fluent field-by-field overrides; ::new stays for compat
    pub fn with_area(self, area_m2: f64) -> Self;
    pub fn with_gross(self, gross_height_m: f64) -> Self;
    pub fn with_contact(self, contact_depth_m: f64) -> Self;
    pub fn with_porosity(self, porosity: f64) -> Self;
    pub fn with_ntg(self, net_to_gross: f64) -> Self;
    pub fn with_sw(self, water_saturation: f64) -> Self;
}
pub struct StructuralPerturbation {                       // Clone + Debug + Default
    pub control_shifts: Vec<(usize, usize, f64)>,         // (ip, jp, dz_m)
}

// A correlated structural perturbation FIELD for one horizon/isochore of a draw
// (decision_structural_uncertainty_isochore): an unconditional Gaussian random field
// (petektools sgs_unconditional) with marginal N(0, sd_m^2) and the variogram's
// spatial continuity, generated on the areal node lattice at realize() time and
// added to a TOP DEPTH surface (RealizationDraw.top_structural) or a zone THICKNESS
// (ZoneDraw.isochore_structural). The mean is variogram-INDEPENDENT (every node's
// marginal is N(0,sd^2) regardless of correlation) — only the field's shape/range
// depends on the variogram. Seed = draw.seed_index salted by the horizon index
// (bit-reproducible per seed, independent across horizons); SGS neighbourhood derived
// from the variogram range. sd_m <= 0 = a no-op (zero field). Perturbation is pinned
// to zero at well-tie nodes. Copy + Clone + Debug + PartialEq + serde.
pub struct PerturbationField { pub sd_m: f64, pub variogram: Variogram }
impl PerturbationField { pub fn new(sd_m: f64, variogram: Variogram) -> PerturbationField; }
```

```rust
// Compile-checked (ratified amendment 1):
static_assertions::assert_impl_all!(StaticModelTemplate: Send);
static_assertions::assert_impl_all!(StaticModel: Send);
```

### Structured Monte-Carlo driver (`task_peteksim_mc_structured`, LANDED)

The driver over the template: samples every uncertain input from a petekTools
`Sampler`, realizes each draw, and KEEPS the per-draw output vectors (W17). `boi`
/ `bgi` are sampled like any input but applied at the **volumetrics surface**,
never on the draw (ratified). Re-exports `Sampler` / `Correlation` /
`ReservoirSummary` from petekTools so a caller needs no direct petekTools import.

```rust
pub enum Input { Plain(Sampler), Clamped(Clamped) }   // Copy + Debug; From<Sampler>/From<Clamped>
impl Input {
    pub fn plain(sampler: Sampler) -> Input;
    pub fn clamped(sampler: Sampler, lo: f64, hi: f64) -> Result<Input>;   // petekTools .clamped (hard limiter)
}

pub struct McInputs { /* one Input per quantity; Clone + Debug */
    pub area_m2, gross_height_m, contact_depth_m, porosity, net_to_gross,
        water_saturation, boi: Input,                 // the load-bearing set
    pub goc_depth_m, sw_gas, bgi: Option<Input>,      // two-contact / gas-FVF
    pub property_shifts: Vec<(String, Input)>,        // per McMode::LevelShift property
}
impl McInputs {
    pub fn new(area_m2, gross_height_m, contact_depth_m, porosity,
               net_to_gross, water_saturation, boi: Input) -> McInputs;
    pub fn with_goc(self, Input) -> Self;  pub fn with_sw_gas(self, Input) -> Self;
    pub fn with_bgi(self, Input) -> Self;
    pub fn with_property_shift(self, property: impl Into<String>, shift: Input) -> Self;
    pub fn realize(&self, n: usize, seed: u64) -> Result<RealizedInputs>;   // Err on n == 0
}

pub struct RealizedInputs { /* per-field Vec<f64> of length n; Clone + Debug */ }  // retained for tornado reuse

pub struct McResult { /* Clone + Debug */
    pub oil_sm3: Vec<f64>,     // primary metric — two-contact oil leg, else whole column
    pub gas_sm3: Vec<f64>,     // gas cap of a two-contact run (needs bgi); else 0
    pub grv_m3: Vec<f64>,      // hydrocarbon-column GRV per draw
}
impl McResult {
    pub fn len(&self) -> usize;  pub fn is_empty(&self) -> bool;
    pub fn realized_inputs(&self) -> &RealizedInputs;
    pub fn summary(&self) -> Result<ReservoirSummary>;       // oil P90/P50/P10 (petekTools reservoir_summary)
    pub fn gas_summary(&self) -> Result<ReservoirSummary>;
    pub fn grv_summary(&self) -> Result<ReservoirSummary>;
}

// THE one structured-MC entry (task_petekstatic_spec_mirror) — the run-resources
// settings value consolidating the four historical fns (now DEPRECATED thin
// wrappers, pinned bit-identical by test). workers=1 serial (default); >1 = the
// rayon-sharded driver (sharded == serial contract). spill_dir None = in-core
// (default); Some(dir) = the spilled mode's per-shard reused f32 store under dir
// (pass std::env::temp_dir() for the old spilled(None)). Serde + PartialEq;
// #[non_exhaustive] (construct via ::new + with_*).
#[non_exhaustive]
pub struct McSettings { pub n: usize, pub seed: u64, pub workers: usize, pub spill_dir: Option<PathBuf> }
impl McSettings {
    pub fn new(n: usize, seed: u64) -> McSettings;             // serial, in-core
    pub fn with_workers(self, workers: usize) -> Self;          // clamp [1, n]; default_mc_workers() recommended
    pub fn with_spill_dir(self, dir: impl Into<PathBuf>) -> Self;  // selects the spilled mode
}

/// Seeded, reproducible: same (inputs, settings) -> bit-identical vectors.
/// FAIL-FAST error policy: the first failing draw surfaces as
/// StaticError::McDraw { index, source } (H2 typed error + the draw index).
pub fn run_mc(tmpl: &mut StaticModelTemplate, inputs: &McInputs,
              settings: &McSettings) -> Result<McResult>;

#[deprecated] pub fn run_structured_mc(..., n, seed) -> Result<McResult>;                      // = run_mc(McSettings::new(n, seed))
#[deprecated] pub fn run_structured_mc_parallel(..., n, seed, workers) -> Result<McResult>;    // = ….with_workers(workers)
#[deprecated] pub fn run_structured_mc_spilled(..., n, seed, spill_dir) -> Result<McResult>;   // = ….with_spill_dir(dir | temp_dir())
#[deprecated] pub fn run_structured_mc_parallel_spilled(...) -> Result<McResult>;              // = both

/// Field aggregation over segments' oil vectors, delegating to petekTools
/// `aggregate`: Independent (narrow) vs Comonotonic (wide downside) bracket.
pub fn aggregate_field(segments: &[&McResult], corr: Correlation) -> Vec<f64>;
```

### Tornado sensitivity (`task_peteksim_tornado`, LANDED)

One-at-a-time swings of **oil in-place Sm³**: each input pivoted at the
**realized** lo/hi percentiles (others at P50), re-realized through the template,
pre-sorted by swing descending. Pass the same `(n, seed)` as the MC run to reuse
its exact draws.

```rust
pub struct TornadoBar { pub input: String, pub lo_val, hi_val, out_lo, out_hi, swing: f64 }  // Clone + PartialEq

pub fn tornado(tmpl: &mut StaticModelTemplate, inputs: &McInputs,
               n: usize, seed: u64, lo_pct: f64, hi_pct: f64) -> Result<Vec<TornadoBar>>;
```

---

### Out-of-core — memory budget + spilled `StaticModel` (`task_petekstatic_slab_streaming`; SPEC §10)

The engine chooses its **backing-storage mode** against a declared budget (out-of-core
rulings R2/R3/R4/R5). Below budget → today's in-core path, **byte-identical**. Above
budget → geometry (ZCORN) + cubes spill to a petekTools `store` (f32 lanes, R4; COORD
resident) and the volumetric surface reads through **windowed mmap views**, streaming
one k-slab at a time (f64 accumulation). Loud switch, never OOM, temp cleanup on drop.

> **Crate origin (P10 organize wave):** these types live in the `srs-spill` crate,
> split out of `srs-model` 2026-07-05 (`task_petekstatic_organize`). `srs-model`
> re-exports the whole surface, so every path below resolves unchanged as
> `srs_model::{MemoryBudget, decide_mode, spill_grid, SpillBacking, …}`.

```rust
// The declared budget (default: DEFAULT_BUDGET_FRACTION of physical RAM).
pub const DEFAULT_BUDGET_FRACTION: f64;              // 0.5
pub struct MemoryBudget { /* private */ }            // Copy + Eq
impl MemoryBudget {
    pub fn bytes(limit_bytes: u64) -> Self;
    pub fn fraction_of_physical(fraction: f64) -> Self;   // clamped [0,1]
    pub fn unlimited() -> Self;                            // never spill (force in-core)
    pub fn limit_bytes(&self) -> u64;  pub fn is_unlimited(&self) -> bool;
}
impl Default for MemoryBudget { /* fraction_of_physical(DEFAULT_BUDGET_FRACTION) */ }
pub fn physical_ram_bytes() -> Option<u64>;          // cached; /proc/meminfo | sysctl hw.memsize
pub fn live_set_bytes(dims: Dims, n_cubes: usize) -> u64;   // ZCORN f64 + cubes f64, ×2 (warm)
pub enum BuildMode { InCore, Spilled }               // Copy + Eq
pub fn decide_mode(dims: Dims, n_cubes: usize, budget: MemoryBudget) -> (BuildMode, u64);
pub struct SpillNotice { pub cells: usize, pub budget_bytes: u64,     // the loud advisory (R5)
    pub estimate_bytes: u64, pub store_path: PathBuf }
impl SpillNotice { pub fn warn(&self); }             // → stderr on every spilled build

// Builder / template opt-in (both StaticModelBuilder and — same setters — the flow):
impl StaticModelBuilder {
    pub fn with_memory_budget(self, budget: MemoryBudget) -> Self;   // default MemoryBudget::default()
    pub fn with_spill_dir(self, dir: impl Into<PathBuf>) -> Self;    // default std::env::temp_dir()
    pub fn with_spill_persist(self, persist: bool) -> Self;          // keep the store past drop
}

// StaticModel: spilled is a backing-storage MODE, not a new type. in_place /
// in_place_summary / bulk_volume / in_place_by_zone route to the store when spilled;
// grid() / property() return a placeholder (read a spilled model through its
// volumetric methods). in_place_by_zone streams per zone (zones are contiguous
// k-bands) — v2, item 3. zone_stats stays in-core-only (typed error spilled).
impl StaticModel {
    pub fn is_spilled(&self) -> bool;
    pub fn spill_store_path(&self) -> Option<&Path>;
    pub fn dims(&self) -> Dims;                       // real dims, in-core or spilled
}

// The ONE GRV/HCPV loop (srs-volumetrics, v2 item 2): a SlabSource abstracts the
// two backings (contiguous f64 Grid vs mmap f32 slab); `compute_clipped` is the
// generic streaming core (three monomorphic loops over a Clip). GridSource is the
// in-core backing; SpillBacking::source() is the spilled one. grv.rs's public
// compute_* fns are thin wrappers over this (signatures unchanged).
pub trait SlabSource { fn dims(&self) -> Dims; fn require_cubes(&self) -> Result<()>;
    type Slab<'a>: CellSlab; fn slab(&self, k: usize) -> Result<Self::Slab<'_>>; }
pub trait CellSlab { fn centroid_z(&self, l: usize) -> f64; fn cell_volume(&self, l: usize) -> f64;
    fn poro(&self, l: usize) -> f64; fn ntg(&self, l: usize) -> f64; fn sw(&self, l: usize) -> f64; }
pub enum Clip { Bulk, Single(f64), Two { goc: f64, owc: f64, sw_gas: Option<f64> } }
pub fn compute_clipped<S: SlabSource>(src: &S, clip: Clip, k_range: Range<usize>, per_cell: bool) -> Result<InPlace>;
pub struct GridSource<'a> { /* &Grid + cube slices */ }  impl GridSource { pub fn new(grid: &Grid) -> Self; }

// Spill a built grid to a store (mmap-backed, f32 lanes); Drop removes the file
// unless detached. Used internally by the builder above budget. spill_streaming is
// the ONE store-writing path — the slab-incremental build (v2 item 1) drives it
// with an on-demand ZCORN + constant-cube producer (never a whole in-core grid);
// spill_grid_to rides it too (its producer reads a built grid).
pub fn spill_grid(grid: &Grid, dir: &Path, cleanup: bool) -> Result<SpillBacking>;
pub fn spill_grid_to(grid: &Grid, path: &Path, cleanup: bool) -> Result<SpillBacking>;
pub fn spill_streaming(path: &Path, dims: Dims, coord: &[Pillar], cube_names: &[String], cleanup: bool,
    fill_zcorn: impl FnMut(usize, &mut [f32]) -> Result<()>,
    fill_cube: impl FnMut(&str, usize, &mut [f32]) -> Result<()>) -> Result<SpillBacking>;
pub struct SpillBacking { /* private: Store + resident COORD + dims + cube names */ }
impl SpillBacking {
    pub fn dims(&self) -> Dims;  pub fn store_path(&self) -> &Path;  pub fn detach(&mut self);
    pub fn cube_names(&self) -> &[String];  pub fn bulk_volume(&self) -> Result<f64>;
    pub fn source(&self) -> SpillSource<'_>;          // the spilled SlabSource (→ compute_clipped)
}

// Spilled structured MC (R3/R4) — run_mc with McSettings::with_spill_dir(dir):
// each draw realizes in-core into the ONE reusable model, is flushed to a REUSED
// per-shard store (overwritten per draw — never a new file per draw), and its
// summary streamed from the f32 lanes → discard. Bit-deterministic within the
// spilled mode; sharded == serial at every worker count. The historical
// run_structured_mc_spilled / _parallel_spilled entries are DEPRECATED thin
// wrappers over run_mc (spill_dir None maps to std::env::temp_dir()).
```

> **Parity is tolerance-based** (R4 honesty clause): f32 storage lanes change volumes
> at a small relative bound (measured on a 2500 m fixture: GRV/bulk bit-exact, HCPV
> ≤ ~2e-9, MC oil ≤ ~3.1e-6 relative — contact-boundary sensitive; asserted ≤ 1e-5).
> **v2 (item 1):** the spilled build is now **slab-incremental** — ZCORN + constant
> cubes stream k-slab-by-k-slab into the store (never a whole in-core grid), so the
> build's heap working set is **O(slab)** and a >RAM model builds. Measured peak RSS
> (3.2M cells): in-core 90 B/cell, build-then-spill 135, **streaming 47** (≈ the f32
> store alone). Eligible: constant/per-zone-constant population, collapse off;
> otherwise build-then-spill. Deferred: streaming SGS-population + writable-store MC
> realize; true O(slab) *resident* RSS awaits a page-evicting store writer
> (`question_petektools_store_page_evict`).

---

### View bundles — the viewer export seam (`view`; SPEC §7e)

Typed inspection bundles for a **separate** viewer codebase
(`decision_viewer_home_product`): petekStatic exports, the viewer renders — it
never computes. The serialized **structure is the contract** — version-tagged and
locked by a schema-snapshot test. `f64::NAN` = undefined (JSON `null`); SI units,
positive-down depth. Bundles **stream** straight to an `io::Write` (no intermediate
`serde_json::Value` tree — the streaming path that removes the legacy ~12.7×-payload
RSS spike): map/section via `write_json`; the volume bundle via its binary-block
envelope writers (below).

```rust
pub const SCHEMA_VERSION: u32 = 6;   // v3: volume exterior-shell + binary blocks
                                     // v4: section HorizonTrace polylines + MapBundle per-well `ties`
                                     // v5: section colour-by-zone — IntersectionBundle.zones + SectionColumn.zone_ids
                                     // v6: oriented GridFrame + optional IntersectionBundle.frame
                                     //     (additive/defaulted; an older decoder ignoring unknown keys still reads v6)
                                     // map/section: plain JSON, carry the family version
impl MapBundle          { pub fn write_json<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()>; }
impl IntersectionBundle { pub fn write_json<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()>; }

// The model's registered WORLD georeference (StaticModelBuilder::with_georef):
// world (x,y) of column (0,0)'s centroid + world column spacing. Some -> the
// view frames are WORLD; None -> the local degenerate frame (synthetic square/box).
pub struct Georef { pub origin_x, origin_y, spacing_x, spacing_y: f64,
                    pub rotation_deg: f64, pub yflip: bool }
impl Georef {
    pub fn new(origin_x, origin_y, spacing_x, spacing_y: f64) -> Option<Georef>; // zero rotation, no flip
    pub fn oriented(origin_x, origin_y, spacing_x, spacing_y: f64,
                    rotation_deg: f64, yflip: bool) -> Option<Georef>;
    pub fn intrinsic_to_world(self, fi: f64, fj: f64) -> (f64, f64);
    pub fn world_to_intrinsic(self, x: f64, y: f64) -> Option<(f64, f64)>;
}

// Shared areal georeference: world (x,y) of node (i,j) = (origin + i*spacing, ...).
// origin/spacing come from the model's Georef when registered (WORLD frame — the
// raster overlays the world outline/wells and a world fence/bore section traces
// through it); otherwise from the grid's local column-centroid lattice.
pub struct GridFrame { pub origin_x, origin_y, spacing_x, spacing_y: f64,
                       pub ncol: usize, pub nrow: usize,
                       pub rotation_deg: f64, pub yflip: bool }   // ncol=ni, nrow=nj
pub struct ValueRange { pub min: f64, pub max: f64 }        // over finite entries (legend)
// A named georeferenced field on the shared frame, row-major values[j*ncol + i].
pub struct ScalarLayer { pub name: String, pub units: String,
                         pub values: Vec<f64>, pub range: ValueRange }
```

**`MapBundle` — `model.map_bundle(&MapSpec) -> Result<MapBundle>`** (areal / plan-view)

```rust
pub struct MapSpec { /* private */ }
impl MapSpec {
    pub fn new() -> Self;                        // structural surfaces + outline + contacts only
    pub fn property(self, name: impl Into<String>) -> Self;  // + zone-average map(s) for this cube
    pub fn k_slice(self, k: usize) -> Self;      // also emit a single-k-slice map per property
}

pub struct WellTieResidual { pub horizon: String, pub residual_m: f64 }                              // v4: one per horizon tied
pub struct WellMarker { pub id: String, pub x, y: f64, pub tie_residual_m: Option<f64>,              // tie_residual_m = mean of ties (None if no ties)
                        pub ties: Vec<WellTieResidual> }                                             // v4: per-horizon residuals (from Provenance.well_ties; empty if none)
pub struct ContactMask { pub kind: String, pub depth_m: f64, pub crossing: Vec<bool> }  // row-major j*ncol+i

pub struct MapBundle {
    pub schema_version: u32,
    pub inputs_ref: String,               // provenance identity
    pub frame: GridFrame,                 // shared by every areal layer below
    pub outline: Vec<Vec<[f64; 2]>>,      // boundary ring(s), world [x,y]
    pub horizons: Vec<ScalarLayer>,       // realized top + base depth surfaces (units "m")
    pub zone_averages: Vec<ScalarLayer>,  // name "{property}::{zone}"  (the useful default)
    pub k_slices: Vec<ScalarLayer>,       // name "{property}::k{n}"    (empty unless k_slice set)
    pub wells: Vec<WellMarker>,           // per-horizon tie residuals (from with_well_ties); empty if no ties
    pub contacts: Vec<ContactMask>,       // per-contact subcrop mask
}
```

Errors (`StaticError::InvalidInput`): column lattice `< 2x2` or not
axis-aligned/regular; a requested property absent; `k_slice` out of range.

**`IntersectionBundle` — `model.intersection_bundle(&SectionSpec, property: Option<&str>) -> Result<..>`** (vertical section)

The trace is marched through the areal lattice at a sub-cell step; each distinct
column crossed becomes an ordered `SectionColumn`. Raw metres — the viewer applies
vertical exaggeration. The structural top/base traces are the first/last *active*
`layer_tops` / `layer_bases`. Under a Follow conformity the active-layer count
varies per column; the arrays stay `nk`-sized and an inactive
(truncated/zero-thickness) layer is `NaN` in `layer_tops`, `layer_bases`, AND
`values` (the viewer NaN-guards).

**Interior-horizon traces (SCHEMA_VERSION 4).** `horizon_traces` carries one
polyline per *interior* framework horizon (every zone-bounding horizon strictly
between the structural top and base — `N − 2` of them for an `N`-horizon stack,
top→down); the structural top/base are NOT repeated. Each `depths` array runs
**parallel to `columns`**: `depths[c]` is that horizon's depth at `columns[c]`,
taken as the zone-top interface (the top-depth of the first cell of the zone the
horizon bounds above). A single-zone (2-surface) model emits an empty
`horizon_traces`, so the block is backward-compatible.

```rust
pub enum SectionSpec {
    Polyline(Vec<[f64; 2]>),                    // world [x,y] trace
    AlongBore { trajectory: Vec<[f64; 3]> },    // world [x,y,z] — xy traces, z overlays
}

pub struct SectionColumn {
    pub distance_m: f64, pub i: usize, pub j: usize, pub x: f64, pub y: f64,
    pub layer_tops: Vec<f64>,   // per-layer cell top depth (len nk), top->base — column CENTROID (hover/back-compat)
    pub layer_bases: Vec<f64>,  // per-layer cell base depth (len nk) — column centroid
    // v4-additive: cell top/base at the column's LEFT/RIGHT fence edges (bilinear from
    // the ZCORN corners at the section's entry/exit x,y), NaN-gapped like layer_tops.
    // With these the viewer draws a dip-following TRAPEZOID; under sugar_cube they equal
    // the centroid trace. FROZEN field names (a concurrent viewer consumes them).
    // Fence direction: for a Polyline it is the trace segment; for an AlongBore it is
    // the trace's AREAL TANGENT through the column (neighbouring column centres, not the
    // raw MD-station micro-segment — which degenerates to the centroid on a vertical /
    // densely-sampled bore). A truly vertical bore (a single areal point, one column)
    // has no tangent, so its edges stay the centroid (l == r) — the one honest flat case.
    pub layer_tops_l: Vec<f64>, pub layer_tops_r: Vec<f64>,
    pub layer_bases_l: Vec<f64>, pub layer_bases_r: Vec<f64>,
    pub values: Vec<f64>,       // per-layer property value (empty if property=None)
    // v5-additive (#[serde(default)]): per-layer zone id — index into IntersectionBundle.zones,
    // NaN-gapped in lockstep with the geometry/value arrays (an inactive layer carries the
    // sentinel SectionColumn::NO_ZONE = u16::MAX). FROZEN field name (the viewer decodes it).
    pub zone_ids: Vec<u16>,
    pub path_z: Option<f64>,    // bore path depth at this station (AlongBore only)
}
impl SectionColumn { pub const NO_ZONE: u16 = u16::MAX; }  // zone_ids gap sentinel (inactive layer)
pub struct SectionContact { pub kind: String, pub depth_m: f64 }
pub struct HorizonTrace { pub name: String, pub depths: Vec<f64> }  // v4; depths parallel to columns
pub struct SectionZone { pub name: String, pub color: Option<String> }  // v5; a zone in the id table (color from StackZone.color)

pub struct IntersectionBundle {
    pub schema_version: u32,
    pub inputs_ref: String,
    pub sugar_cube: bool,               // v4-additive: true = flat-box cells (edge arrays flattened to centroid); false (default) = dip-following trapezoids
    pub property: Option<String>,       // which property `values` carry
    pub top_name: String, pub base_name: String,   // structural-surface labels
    pub columns: Vec<SectionColumn>,    // in trace order (by distance_m)
    pub horizon_traces: Vec<HorizonTrace>,  // v4: interior-horizon polylines, top->down (empty for a single zone)
    pub zones: Vec<SectionZone>,        // v5-additive: the zone table zone_ids indexes into ([{name,color}], top->base). FROZEN
    pub contacts: Vec<SectionContact>,  // GOC / OWC / GWC depths
    pub frame: Option<GridFrame>,       // v6-additive/defaulted: world frame for columns.x/y; None when reading pre-v6
}
```

`rotation_deg` is finite degrees counter-clockwise from world +X/east to the
positive I axis and normalizes to `[0,360)`; `yflip` reverses positive J. The
zero/false fields are omitted on serialization, so legacy Georef/Frame JSON
keeps its exact member sequence. Section `frame` is appended and optional so a
pre-v6 payload deserializes and reserializes without shape drift.

**Section colour-by-zone (SCHEMA_VERSION 5, `task_suite_section_zone_color`).** The
bundle's `zones` list (`[{name, color}]`, the model's stratigraphic zones top→base)
plus each column's per-layer `zone_ids` are the payload half of the viewer's
colour-by-zone rendering. `zone_ids[k]` indexes `zones`; it is `NaN`-gapped in
lockstep with the geometry/value arrays — an inactive/truncated layer (where those
are `NaN`) carries `SectionColumn::NO_ZONE` (`u16::MAX`), which the viewer skips
exactly as it skips a `NaN` depth. `color` comes from `StackZone::color` (`None` on
the single-implicit-zone paths). Both fields are `#[serde(default)]` (additive). The
field names `zone_ids` and `zones: [{name, color}]` are **frozen** — a concurrent
viewer decodes exactly them.

Errors (`StaticError::InvalidInput`): lattice `< 2x2` or not regular;
trace with `< 2` vertices; a named property absent.

**`VolumeBundle` — `model.volume_bundle(property: &str) -> Result<VolumeBundle>`** (3-D exterior shell, SCHEMA_VERSION 4)

The corner-point cell **exterior shell** coloured by `property`: only faces
bordering an inactive/absent neighbour or the grid boundary are emitted (a cell is
active iff its property value is finite and its bulk volume is positive). For a
low-relief grid this is O(surface) not O(volume) — a 1M-cell (200×200×25) box drops
from ~12M triangles to ~200k (1.67%), and the serialized payload from ~557 B/cell
decimal-text JSON to **~6.65 B/cell** (self-contained) / **~5 B/cell** (sidecar).
Shared vertices are deduplicated by f32 position; the per-cell arrays are compacted
to the shell cells only, with a `tri_cell` index per triangle recovering cell
identity for the viewer's threshold filter / picking. The mesh builder **lives
here** — moved DOWN from petekSim's `srs-core/mesh.rs` (the DAG flows downward).

```rust
pub struct VolumeBundle {
    pub schema_version: u32,       // current family version (6; binary layout introduced in 4)
    pub inputs_ref: String,
    pub property: String,
    pub cell_count: usize,         // total grid cells (cell_values.len() = shell count)
    pub positions: Vec<f32>,       // deduped shell verts, 3 (x,y,z) each — grid-LOCAL coords
    pub indices: Vec<u32>,         // 3 per triangle, into positions
    pub tri_cell: Vec<u32>,        // compact cell index per triangle (into cell_values/zone_ids)
    pub cell_values: Vec<f32>,     // per shell cell (compact) — the filter/legend array
    pub zone_ids: Vec<u16>,        // per shell cell (compact) — index into zone_names
    pub zone_names: Vec<String>,
    pub value_range: ValueRange,   // over the shell cell values
}
impl VolumeBundle {
    // Self-contained: metadata envelope + base64-wrapped binary blocks, ONE file (save_view).
    pub fn write_self_contained<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()>;
    // Served: JSON envelope (offset,length manifest) + raw model.bin blocks, no base64 tax.
    pub fn write_sidecar<W1: std::io::Write, W2: std::io::Write>(&self, json: &mut W1, bin: &mut W2) -> std::io::Result<()>;
}

// Threshold regeneration (decision (b)): client slider is shell-only; the server
// re-cuts the shell at a cutoff (exposing revealed interior) via:
model.volume_bundle_thresholded(property: &str, cutoff: f64, keep_above: bool) -> Result<VolumeBundle>
```

Errors (`StaticError::InvalidInput`): `property` absent (both methods).

**Binary-block payload spec (SCHEMA_VERSION 4 — the viewer decode contract).** The
big arrays are emitted as raw **little-endian**, tightly packed (no padding within
or between blocks) bytes; the JSON envelope stays human-readable
(names/units/ranges/zones/provenance). dtypes: `f32` (4 B), `u32` (4 B), `u16`
(2 B); a NaN f32 is canonical `0x7FC00000`. Block order = row-major (C-order)
flatten of `shape`.

```jsonc
{ "schema_version": 4, "kind": "volume", "inputs_ref": ..., "property": ...,
  "cell_count": N, "shell_cell_count": C, "vertex_count": V, "triangle_count": T,
  "zone_names": [...], "value_range": {"min":..,"max":..},
  "encoding": "base64" | "sidecar",
  "blocks": {
    "positions":   { "dtype":"f32", "shape":[V,3], <payload> },
    "indices":     { "dtype":"u32", "shape":[T,3], <payload> },
    "tri_cell":    { "dtype":"u32", "shape":[T],   <payload> },
    "cell_values": { "dtype":"f32", "shape":[C],   <payload> },
    "zone_ids":    { "dtype":"u16", "shape":[C],   <payload> } } }
// <payload> = "data":"<base64 of the LE bytes>"   when encoding=="base64"
//           = "offset":<byte>, "length":<byte>    when encoding=="sidecar"
//             (blocks concatenated in the above order in the companion model.bin)
```

---

## Python surface

The `petekstatic` wheel exposes the Rust-backed `StaticModel` / `build_flat_model`
surface and the first notebook-facing workflow facade:

```python
petekstatic.__all__ == [
    "CoKriging", "DistributionSpec", "Grid", "HorizonSpec", "Layering",
    "PropertyHandle", "PropertyPipeline", "PropertyPipelineSpec",
    "PropertyStore", "SgsRecipe", "Spherical", "StaticModel",
    "UpscaleRecipeBuilder", "Var", "VolumeCase", "VolumeResult",
    "WellTie", "WellLog", "WellLogSpec",
    "__version__", "build_flat_model",
    "distributions", "upscale",
]
```

`build_flat_model` returns a single-zone `StaticModel` suitable for smoke tests,
volume reads, and bundle reads. `Grid.from_project(...)` owns the canonical
static declaration shape in Python: geometry, horizons, zones, layers, scalar
and formula property assignment, property recipe lowering/execution handles, and
deterministic simple volumes. Full production corner-point model building, contact scenarios,
bundles, and MC remain backed by the Rust `StaticModelBuilder` path and may
still be surfaced through petekSim's downstream product facade while the facade
lowering is completed.

---

## Change control

The cross-library seams here (`StaticModel` accessors, `RealizationDraw`, the
regeneration API, `solve_surface_seeded`/`KernelSurface`) are ratified contracts:
changing a signature requires (1) an edit to this file, (2) coordinator
(petekSuite) + consumer (petekSim) sign-off. Library-internal signatures are
petekStatic's own call. The whole file has been locked since the 0.1 release.
