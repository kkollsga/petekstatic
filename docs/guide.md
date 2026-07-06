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

> **Python surface today.** petekStatic has a minimal own wheel (`petekstatic`,
> exposing `StaticModel`, `build_flat_model`, and `__version__`) for the compact
> single-zone API. The **rich geomodel workflow** — project inventory, explicit
> role binding, scatter-conditioned horizon stacks, multi-zone volumes, contact
> scenarios, bundles, and MC — is driven through the **`peteksim` facade**, which
> presents petekStatic's engine to Python. Synthetic project trees come from the
> horizontal `petektools.synth_asset` unit; the notebooks treat that as a source
> of example input, not as petekStatic functionality.

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

### 3. Property modelling — priors and upscaling

`srs-petro` populates the per-cell cubes. Today: constant/day-1 **priors** and
**log upscaling**; geostatistical population (facies, variogram-driven
propagation, trends) is the growth path. A property can be **net-only**
(populated only where net-to-gross qualifies), and each property carries its own
variogram + seed so a run is deterministic and reproducible.

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
in separate cells so the swap point is obvious. Run them with a Python that has
`petektools` and `peteksim` installed.

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
