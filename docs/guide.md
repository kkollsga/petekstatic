# petekStatic guide — building a populated static reservoir model

petekStatic is the **GEOMODEL layer** of the petek subsurface-modelling stack. It
turns model-ready inputs into a **populated `StaticModel`** — a corner-point grid
carrying property cubes, zones, and fluid contacts — and it **owns the volumetrics
and static-uncertainty stack** over that model (GRV / in-place, and Monte-Carlo
regeneration over model realizations). Everything downstream (dynamic simulation,
the appraisal product) reads the `StaticModel` across a one-way seam.

This is a guide, not a reference: it walks the geomodel workflow end to end and
points at the two runnable notebooks in `examples/notebooks/`. For exact
signatures see `API.md`; for the design constitution see `SPEC.md`.

> **Python surface today.** petekStatic exposes the compact Rust-backed
> `StaticModel`/`build_flat_model` surface and the first own Python workflow
> facade: `Grid.from_project(...).geometry(...).horizons(...).zones(...).layers(...)`,
> `grid.properties.calc(...)`, `upscale(...).sgs(...)` property recipes, and
> `grid.volumes(...).run(...)`. The workflow facade is a declarative/spec layer
> for notebooks and smoke tests: it validates project asset names, stores
> property arrays in memory, delegates formula blocks to
> `petektools.evaluate_formula`, lowers property recipes to
> `PropertyPipelineSpec`, and returns deterministic simple volumes. It does
> **not** yet replace the production Rust corner-point grid build, bundles,
> contact scenarios, or MC paths that are still surfaced through the `peteksim`
> facade while the lift continues.

## Where petekStatic sits

Dependencies flow one direction, downward only:

```
petekIO      DATA       → model-ready inputs (ModelInputs / .pproj)      [upstream]
   ↓
petekStatic  GEOMODEL   → populated StaticModel + volumetrics/uncertainty [here]
   ↓
petekSim     SIMULATION → dynamic/engineering + the Python product        [downstream]

petekTools   TOOLKIT    → gridding/kriging/warm-start kernels + units     [horizontal]
```

petekStatic **consumes** petekIO's model-ready inputs and **calls** petekTools'
numeric kernels; it never depends on petekSim (its consumer). No cycles, no
sideways code-sharing — conventions are shared, small types are converted at the
seam (e.g. the FVF value types that cross from PVT as validated scalars).

## The single crate, layered by module

petekStatic now publishes one Cargo crate. The historical `srs-*` crate names
survive as module/layer labels, and dependencies still point strictly downward:

| module/layer | responsibility |
|---|---|
| `error` | the one error enum (`StaticError`) + `Result<T>`; chains petekIO's `GeoError` at the ingest seam |
| `srs-grid` | the i,j,k corner-point grid: cells, corners, per-cell property cubes, gross rock volume |
| `srs-gridder` | grid construction — layering + conformity + collapse over the framework |
| `srs-wireframe` | the structural framework: boundary, horizons, contacts, the `HorizonStack` |
| `srs-data` | the petekIO → geomodel ingest boundary (depth sign flip, frame/georef registration) |
| `srs-petro` | property modelling — priors, log upscaling, geostatistical population |
| `srs-volumetrics` | GRV / in-place (OOIP / OGIP), single- and two-contact |
| `srs-uncertainty` | Monte-Carlo regeneration + tornado over static realizations |
| `srs-model` | the `StaticModel` aggregate + builder + the ratified MC-regeneration seam (top of the DAG) |

## The workflow, layer by layer

### Python declaration slice

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
p["HA_FWL"].set(1.0)
p["Jfunc"].set(1.0)
p.calc([
    "RQI = $lambda * sqrt(PermXY_BC / PorE_BC)",
    "Swirr = $SHF_c * pow(RQI, $SHF_d)",
    "Sw = if(HA_FWL == 0, 1, Swirr + (1 - Swirr) * $SHF_a * pow(Jfunc, $SHF_b))",
], params={...})

result = grid.volumes(ntg="NTG", por="POR", sw="Sw", fluid="oil", fvf=1.30).run(progress=True)
```

`progress=True` prints coarse stages; passing a callback receives event
dictionaries. `VolumeResult.summary()` returns total GRV/HCPV/in-place numbers
and `VolumeResult.by_zone()` returns the current zone dictionary.

Property recipes are declared the same way, but they produce a lowered pipeline
spec instead of an in-memory property vector:

```python
logs = project.logs
vgm = pst.Var("spherical", major=1500, minor=700, vertical=20, azimuth=35)

p.por = pst.upscale(logs.PHIE(logs.NetSand > 0.50)).sgs(
    variogram=vgm,
    distribution=pst.distributions.from_logs(),
    seed=12,
)

declared = p.declarations("por")
lowered = p.pipelines("por")
```

When the project or log source can resolve positioned wells, an isotropic recipe
can build a Rust-backed `pst.PropertyPipeline` handle:

```python
iso = pst.Var("spherical", major=1500, minor=1500, vertical=1500, azimuth=0)
p.ntg = pst.upscale(logs.NetSand).sgs(
    variogram=iso,
    distribution=pst.distributions.from_logs(),
    seed=11,
)

pipe = p.execute_pipeline("ntg")
smoke_model = pipe.apply_to_flat_model()
```

That execution is intentionally narrow today: unresolved lazy logs,
cokriging/trend binding, non-`from_logs` distributions, and anisotropic Rust
execution raise explicit errors. Anisotropic `pst.Var(...)` specs are preserved
in the serialized recipe and can still lower to petekTools' anisotropic variogram
object; arbitrary-grid pipeline application is not exposed through Python yet.
Use `PropertyPipeline.apply_to_flat_model(...)` only as the current smoke path.

### 1. Structural framework — horizons + zones

The framework is a stack of depth **horizons** bounding **zones** (one zone per
adjacent horizon pair), plus per-zone fluid **contacts**. Horizons arrive either
as gridded surfaces or as **raw scatter** (points) that the engine conditions onto
a common frame. Well **ties** pull the surfaces to picks; the family default is
the **convergent** tie method. N horizons resolve to N−1 zones, each with its own
conformity and contacts.

### 2. Grid construction — the convergent corner-point grid

`srs-gridder` builds an i,j,k **corner-point grid** on a georeferenced frame:
each zone is layered (`nk` layers, proportional/other conformity), negative
thickness is collapsed, and a minimum-thickness floor keeps degenerate cells out.
The scatter path (`from_scatter_stack`) is the *single scatter-gridding authority*
— it conditions raw horizon points onto the frame (bilinear, voids left `NaN`)
and then resolves exactly as the gridded-stack path, so a model and its MC
template are built from bit-identical geometry.

### 3. Property modelling — priors, formulas, and upscaling

`srs-petro` populates the per-cell cubes. Today: constant/day-1 **priors**,
Python-side formula assignment, and **log upscaling** lowered toward the Rust
`PropertyPipeline`. A property can be **net-only** (populated only where
net-to-gross qualifies), and each recipe carries its own variogram + seed so a
run is deterministic and reproducible.

### 4. Volumetrics — GRV and in-place, per zone

`srs-volumetrics` reads volumes off the model itself. `in_place()` gives whole-
column GRV / HCPV / OOIP / OGIP; `in_place_by_zone()` computes in-place against
**each zone's own contacts** and returns a per-zone breakdown plus a rollup total
(the zone sum equals the total). A zone with an OWC yields oil; a zone with a
GOC + FWL is a **two-contact** gas-cap-over-oil split. FVF (`Boi` / `Bgi`) enters
as a validated scalar — petekStatic holds no PVT code.

### 5. Static uncertainty — MC regeneration over realizations

`srs-uncertainty` regenerates the model under a `RealizationDraw` via
`StaticModelTemplate::realize(&draw)`. The template builds the draw-invariant
framework **once** and varies only what a draw perturbs (structure, zone contact
levels, property draws) — bit-deterministic and buffer-recyclable on the hot path.
Running many draws yields the in-place **P-curve** (P90/P50/P10); tornado
attribution is the next increment. Zero-spread draws reproduce the deterministic
model exactly (the canonical planted-truth check).

## The StaticModel seam to petekSim

The output is a `StaticModel`: framework + grid + property cubes + zones +
contacts + provenance, with `in_place*()` and the MC-regeneration seam on top. It
also emits self-contained **view bundles** as JSON — `map_bundle` (areal,
georeferenced raster + per-zone average maps), `intersection_bundle` (a vertical
cross-section along a polyline), and `volume_bundle` (the corner-point exterior
shell for 3-D view, base64-wrapped binary blocks). petekSim's facade presents
these results; it never reaches back up into geomodel internals. Changing a
seam signature (the StaticModel accessors, `RealizationDraw`, the regeneration
API) requires coordinator + consumer sign-off — the seam is a contract.

## The example notebooks

Both live in `examples/notebooks/`. They can be executed end to end on the
synthetic project tree produced by `petektools.synth_asset`, or pointed at a
user export by changing the data-source cell and replacing the visible role
literals with names from `Project.inventory()`. Each notebook keeps synthetic
generation, `Project.load`/inventory inspection, and the main modelling workflow
in separate cells so the swap point is obvious. Run the petekStatic workflow
slice with `petekio`, `petektools`, and `petekstatic`; the older product-facing
notebooks also need `peteksim` until their final build cells are fully lowered.

- **`01_stack_model_from_scatter.ipynb`** — load a project tree, print
  `Project.inventory()` before binding, declare user-replaceable literals for
  `outline="ModelEdge"`, ordered horizons, zones/subzones, contacts, and
  properties, then build a multi-zone stack model through the `peteksim` facade
  and read **per-zone volumes** (`in_place_by_zone`). Synthetic manifest checks
  are isolated in the final cell and skipped for real-data paths.
- **`02_volumes_and_bundles.ipynb`** — a per-zone GRV / in-place table, a
  **contact scenario** (replace a zone's OWC with a deeper one → more oil in that
  zone), a **bundle** peek (`volume_bundle` / `map_bundle` keys), and the current
  segment convention: run one model/MC per segment and combine segment
  realization sets with `ps.aggregate`; do not build a zone-by-segment
  single-model rollup in petekStatic.

## Conventions

Depths are metres, **positive down** (petekIO's negative-down elevation is
flipped at the `srs-data` ingest boundary). Per-cell cubes are dense vectors
indexed by linear cell index; `NaN` marks *undefined*. Every result is a
`Result<T, StaticError>`. Data in committed tests, examples, and this guide is
**synthetic only** — petekStatic never ships or references real-dataset content.
