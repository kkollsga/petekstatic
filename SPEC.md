# petekStatic — design constitution (build spec)

> Repo: `Koding/Rust/petekSuite/petekStatic` · crate: `petekstatic`
> (modules preserve the historical `srs-*` layer labels). The **GEOMODEL** layer of the petek subsurface-modelling
> ecosystem. This is the design constitution (the *why*/*how*); the locked public
> contract is `API.md`; the shared conventions are the petek family house style
> (`petekSuite/dev-docs/petek-house-style.md`). Cross-library seams + lifecycle
> live in the suite planning graph (served by the `contract` MCP).

petekStatic turns **model-ready inputs** into a populated **`StaticModel`** — a
structural framework (horizons + faults +
zones), a modelling grid, per-cell property cubes — **and owns the volumetrics +
static-uncertainty stack over it**: GRV / in-place (OOIP/OGIP) off the model,
and Monte-Carlo regeneration over static-model realizations. petekSim (the
dynamic/engineering + product layer) consumes the model and its volumetric
results across the seam.

**Litmus test for what belongs here (re-scoped 2026-07-03, graph
`decision_layer_charters`):** it is a *static-modelling* concern — anything a
geomodeller or a static-uncertainty study needs: the populated grid, its
in-place volumes, and the MC loop over static realizations (tornado later).
**PVT/fluid correlations, dynamic flow, decline/material-balance, economics, and
the Python product facade are petekSim's.** FVF crosses this seam **as an
uncertain scalar input** (a validated value type), never as PVT code. Nothing
input-data-specific (parsing, QC, interpretation) lives here; that is petekIO's
job. The default core owns its input types and **calls** petekTools kernels;
petekIO conversion is an optional compatibility feature or lives in petekSim's
composition root. We **never** depend on petekSim (our downstream consumer).

---

## 1. Design constitution (this library's slice of the family house style)

1. **Strictly layered, one-way internal deps; no cycles.** The current module DAG
   inside the single `petekstatic` crate:

   ```
   error                     # the one workspace error enum (StaticError) + Result
   srs-grid                  # i,j,k corner-point grid: Grid/Cell/KLayer/Property
   srs-wireframe             # the constraining wireframe (boundary+horizons+contacts)
   srs-gridder               # convergent gridder: min-curvature surface + layering
   srs-petro                 # petrophysics: log upscaling (power-law means)
   srs-data                  # optional petekio-adapter compatibility feature
   srs-volumetrics           # GRV + in-place (OOIP/OGIP) + FVF value types + range validation
   srs-uncertainty           # Monte Carlo: distributions, SplitMix64, P90/P50/P10 (leaf)
   srs-spill                 # out-of-core backing mode: memory budget + k-slab spill
                             #   onto a petekTools store (imports srs-grid + srs-volumetrics)
   srs-model                 # TOP OF DAG: the StaticModel aggregate + builder +
                             #   the MC-regeneration template (relocated from petekSim)
   ```

   A module imports only from layers above it. `srs-volumetrics` and
   `srs-uncertainty` **relocated here from petekSim 2026-07-03** (origin SHA
   `fe6343c`, `task_relocate_refine_orchestration`) under the layer-charter
   re-scope; `srs-spill` **split out of `srs-model` 2026-07-05**
   (`task_petekstatic_organize`, P10 organize wave) as the out-of-core
   backing-storage mode below the aggregate; `srs-model` composes grid + gridder +
   petro + wireframe + volumetrics + spill into the `StaticModel` and owns the
   regeneration orchestration relocated from petekSim's `RefiningModel`. `srs-model`
   re-exports `srs-spill`'s surface, so the split is source-compatible.

2. **One error enum** (`StaticError`, `thiserror`) + `Result<T>` everywhere;
   downstream composes it with a `#[from]` variant on its own enum (house §1).
   petekSim does this: `SrsError::Static(#[from] petekstatic::error::StaticError)`.
   The optional `petekio-adapter` composes `StaticError::Geo(#[from]
   petekio::GeoError)` so legacy DATA→GEOMODEL conversions retain their source
   chain. The default geomodel core has no petekIO dependency; petekSim owns the
   integrated composition.

3. **Domain objects carry their operations** — fluent, chainable, immutable; ops
   return *new* objects, mutation is explicit `set_*`. `f64::NAN` = undefined.

4. **Compartmentalized.** One module/topic, one type/responsibility; split before a
   file owns two jobs. Boundaries are traits where backends are plural, enums where
   the set is small and closed.

5. **Compose deps, don't reinvent.** The convergent gridder / kriging / warm-start
   kernels are petekTools' job; unit conversions are petekTools'; input-data work is
   petekIO's. We map at the seam and convert small types rather than share code
   sideways. (The FVF value types in `srs-volumetrics` are the sanctioned
   duplicate-at-the-seam: petekSim's `srs-pvt` keeps its own copies for the
   dynamic work; no dependency crosses.)

6. **Physical-range validation is typed, at the seam.** Every volumetrics input
   passes `validate_fraction` (φ/NG/Sw ∈ [0,1]) / `validate_positive`
   (area/height) / the FVF constructors — a non-physical input is a
   `StaticError`, never silent garbage (`-117 MMSTB`, `inf`).

7. **PyO3-ready by construction.** Owned types on the public surface, no public
   lifetimes, plain numeric types / `ndarray`. Bindings (petekSim's facade today)
   only marshal.

---

## 2. The `StaticModel` — the contract (LANDED, `srs-model`)

A **`StaticModel` is a *populated* static reservoir model**: the single, complete,
self-describing artifact petekStatic produces. It is a value (owned, `Clone`,
`Send`), not a handle to a live builder. Its anatomy:

| part | type | what it carries |
|---|---|---|
| **structural framework** | `Wireframe` (`srs-wireframe`) | areal boundary + depth horizons (`Top`/`Intermediate`/`Base`) + fluid contacts, each tagged with `Hardness` |
| **grid** | `Grid` (`srs-grid`) | corner-point geometry + `nk` `KLayer`s. The unfaulted box is the degenerate case; the same type upgrades in place to the faulted grid |
| **property cubes** | `Properties` (`srs-grid`) | name → `f64[cell_count]` per-cell scalar fields (`PORO`, `SW`, `NTG`, …); `NaN` = undefined |
| **zones** | `ZoneTable` (`srs-model`) | named stratigraphic intervals bounding k-ranges; MVP ships the single implicit whole-column zone (§4) |
| **contacts** | `Vec<Contact>` (on the framework) | OWC / GOC / GWC depths — the model's own volumetrics clip against these |
| **provenance** | `Provenance` (`srs-model`) | inputs identity, gridder settings, population mode, and the `RealizationDraw` (with `seed_index`) for an MC realization (§6) |

**Invariant:** every property cube's `values.len() == grid.cell_count()`; every
zone's k-range lies within `[0, nk)`; the grid's areal lattice matches the
framework's. Enforced at construction (the constructor is crate-internal; only
the builder/template build models) — the consumer never re-validates.

**Units — SI everywhere, one convention** (`decision_si_units_standard`,
2026-07-04). All coordinates, depths and lengths are **metres**; areas **m²**;
volumes **m³** internally, reported in **mcm** (`1e6 m³`, GRV), **MSm³** (oil)
and **bcm** (`1e9 Sm³`, gas). Vertical is **positive-down depth** (clearly-named
`*_m` accessors), the domain-natural framing the ruling permits, and it is the
**one datum inside the `Wireframe`**: petekIO delivers surfaces (and well-curve
positions) as negative-down subsea elevation, so the optional adapter **negates
them at the ingest boundary** (`surface_depths`) onto positive-down `depth_m`, matching the
contacts petekIO already delivers positive-down. This unifies horizons and
contacts on one convention (so structural role assignment reads `Top` =
shallowest correctly) and retires the coordinate-flip deferral — its blocker,
petekIO's imperial `SummaryInputs` contacts, resolved when petekIO went SI
(2026-07-04). FVF (`OilFvf`/`GasFvf`) is
dimensionless **Rm³/Sm³** — numerically identical to the legacy rb/STB & rcf/scf,
so it is a **relabel, not a conversion**. Imperial is opt-in only (petekTools
factors), never a default.

**The model OWNS its volumetric answer** (inverted 2026-07-03 by
`decision_layer_charters`, superseding the earlier "computes no volumes"
stance): `model.in_place()` clips the column against the model's contact(s) and
returns `InPlace` (GRV + HCPV + per-cell HCPV); the caller applies FVF through
the validated `OilFvf`/`GasFvf` value types (`ooip_sm3`/`ogip_sm3`, +
`oil_msm3`/`gas_bcm` reporting scales). P-curves aggregate per-realization
in-place values via `srs-uncertainty`'s `PercentileSummary`.

**Two-contact columns (gas cap + oil rim).** A framework carrying a GOC *and* a
lower contact (OWC/FWL or GWC) makes `in_place()` partition cells by centroid
depth — gas cap (`z < GOC`) / oil leg (`GOC <= z < lower`) / water — and report
per-zone GRV + HCPV (`InPlace::gas` / `::oil`, `gas_zone_ogip_sm3` /
`oil_zone_ooip_sm3`). A lone contact stays the generic single-column answer.
**The split is geometry + in-place only**: free-gas OGIP off the cap, STOIIP off
the leg. What the model still does NOT carry: PVT correlations (solution gas Rs,
gas-cap expansion, condensate), economics, dynamic behaviour — petekSim's; FVF
still enters only as a validated scalar per zone.

---

## 3. Construction path — how a `StaticModel` is built

The build pipeline, ingest → framework → grid → population → model:

```
petekIO ModelInputs (.pproj)
   │  petekSim composition root (or optional petekio-adapter compatibility)
   ▼
Wireframe  {boundary, horizons, contacts}         + PetroSample logs [(tvd,φ,Sw)]
   │  warm_surface  (petekTools warm kernel; the SAME path the MC template uses — R2)
   ▼
solved Top/Base surfaces
   │  Surface::guard_below  (base-above-top → CrossedSurfaces, or clamp — R1)
   │  srs-gridder::layer_grid  (conformable k-layering -> corner-point Grid + report)
   │    Conformity: Proportional (equal fraction, honours nk) | FollowTop{dz_m} |
   │    FollowBase{dz_m} (layers parallel to top/base at constant dz; nk dz-derived
   │    = ceil(max thickness/dz), capped MAX_NK; deep/shallow layers TRUNCATE at the
   │    pinch-out horizon → zero-thickness cells, counted in provenance warnings)
   ▼
Grid (geometry only)
   │  srs-petro upscaling  (φ length-weighted, Sw pore-volume-weighted)
   │  population: logs-in-cell-range -> upscaled cube values; else priors
   │  + optional areal trend: lateral NTG (+φ) shape (§3b)
   ▼
StaticModel  {framework, grid, cubes, zones, contacts, provenance}
```

This pipeline is `srs_model::StaticModelBuilder` (relocated from petekSim's
`RefiningModel` 2026-07-03; petekSim's `srs-core` now keeps only a thin facade
that delegates here and applies FVF). Entry points: `from_wireframe(wf, opts)`
(the data-layer hand-off) and `flat(...)` (the model-first flat-box start);
`add_top_control` grows the control set for the live-refine loop; `build()`
converges the current state. See `API.md`.

**One declarative `BuildSpec` (the suite api-consistency contract,
`task_petekstatic_spec_mirror`).** The builder and the MC template share a single
declarative configuration value — `BuildSpec` (georef, boundary, extrapolation,
repair/clamp/collapse floors, sugar-cube, sw_gas, well ties + `TieSettings`,
inputs_ref) — held internally by both; every `with_*` setter is thin sugar
mutating it and `with_spec` installs a whole one. It serializes (a scenario is a
savable file; the whole config layer derives serde — R7 battery in
`tests/spec_conformance.rs`; `McInputs` joins once petekTools' samplers derive
serde), compares by value, and is `#[non_exhaustive]` so the structural-
uncertainty per-zone isochore leg (`decision_structural_uncertainty_isochore`)
lands additively — as do new fields on the `#[non_exhaustive]`
`RealizationDraw`/`ZoneDraw`. Well ties fold in per `TieSettings`: `Replace`
(default; datum substitution — the tie moves exactly the tied node on a dense
lattice) or `Radius { radius_m }` (bounded locality: the residual decays linearly
to zero at the radius across defined datums; beyond is bit-untouched; voids stay
the solve's to taper). One tie authority (`substitute_tie_datums`) serves both
consumers. MC run resources are the separate `McSettings { n, seed, workers,
spill_dir }` behind the single `run_mc` entry (the four historical drivers are
deprecated thin wrappers).

### 3a. Base-horizon relief (Fix S1)

`from_wireframe` wires a supplied `Base` horizon's **real relief** through the
build: its defined nodes solve a base surface the pillars follow, so gross
thickness varies spatially. `BuildOpts.gross_height_m` is the **fallback**
when no Base horizon is supplied (backward-compatible constant offset). Supplied
horizons the build cannot consume yet (Intermediate; a second or
lattice-mismatched Base) are reported non-blockingly on
`Provenance.warnings` (`BuildWarning::UnusedHorizon`), never silently dropped.
The MC template honours base relief too (`decision_template_gross_scaling`):
built from a wireframe carrying a Base horizon, it extracts the nominal per-node
gross field `g(x,y)` (base solved once in kernel space against the nominal top)
and each draw's `gross_height_m` acts as a **level** on that **shape** —
per-node gross = `g × gross_height_m / mean(g)`, so a draw at `mean(g)`
reproduces the deterministic build. Without a Base horizon the draw keeps the
constant-offset behaviour exactly (backward compatible).

**Builder/template kernel unification (R2, 2026-07-04).** Both the deterministic
build and the template's `realize` solve their surfaces through one shared
petekTools warm-kernel path (`warm_surface`: flat bootstrap at the mean control
depth, then a seeded solve). Previously `build()` used the cold `solve_surface`
(natural-dip reference kernel) while the template used the warm kernel; the two
kernels interpolate **unpinned** interior nodes differently, so on thin columns —
where gross is a small difference of two large surfaces — GRV diverged ~20× (a
fully-pinned lattice hid it, since both kernels then just return the control
values). Unifying on the warm kernel makes `realize(mean-gross draw)` reproduce
`build()` within tight tolerance on all column shapes; cold `solve_surface`
remains the accuracy reference, off the build path
(`decision_gridder_kernel_unification`).

**Base-above-top guard (R1, 2026-07-04).** After the surfaces solve, a per-node
`Surface::guard_below` asserts the base sits at/below the top (positive gross). A
crossing is a typed `StaticError::CrossedSurfaces { nodes, worst_m }` (default),
or — with `with_clamp_base_to_top(true)` — the offending columns clamp to zero
gross and the rest are untouched. Thin/crossing frameworks no longer silently
collapse GRV.

**Post-gridding order-repair (R-c, 2026-07-04).** The `CrossedSurfaces` guard is
loud but leaves the caller stuck when a real thin margin genuinely crosses:
independent gridding of Top and Base overshoots at the edges and re-introduces a
crossing a pointwise pre-repair had removed. `with_min_thickness_m(f64)` (on both
the builder and the template) opts into a remedy: `Surface::repair_min_thickness`
pulls every node where `base_z - top_z < min_thickness` **down** to exactly
`top_z + min_thickness`, **preserving the top** — the top is the better-constrained
seismic pick, so the softer base yields to it (repair direction = downward). Off
by default (the guard stays the default); when enabled, `Provenance.warnings`
gains a `BuildWarning::ThinColumnsRepaired { columns, worst_m }` recording the
repaired-node count and the worst (most-negative) original separation. It takes
precedence over `with_clamp_base_to_top` when both are set (repair, don't zero).

**Solve fidelity (2026-07-04, structure-fidelity audit S1/S2).** The stack
surface solve honours every **defined input node exactly** (hard controls,
bit-level) — the on-data contract. Between nodes, the solve is iterated to a
true **fixed point**: the kernel's per-sweep stopping rule can stall metres
short of convergence on a sparse-control lattice (the smooth modes under the
natural-dip boundary relax slowly), so `warm_surface` re-seeds the solve from
its own output until a whole re-solve moves < 1 mm (a converged re-solve exits
in ~1 sweep, so the loop is cheap once settled). Scatter→node **conditioning is
upstream of this seam**: the builder consumes `GriddedDepth` node data; if the
upstream authoring snaps scatter to nearest nodes on a coarser lattice, the
resulting on-data misfit is an authoring artifact the audit fixture quantifies
(`tests/structure_fidelity.rs` S1) — author node values by local interpolation.

**Extrapolation policy (2026-07-04, audit S3).** Every stack surface/isochore
solve is subject to an explicit [`ExtrapolationPolicy`]
(`with_extrapolation` on builder + template). Default **`DecayToData`**: within
`start_cells` (2) of the nearest datum the solve is untouched; beyond, it blends
linearly toward the **nearest datum's value** over `decay_cells` (4), so a data
void inside the model extent flattens to nearest-data instead of running the
kernel's natural-dip linear extension unbounded (or, under-converged, freezing
at the seed mean — both uncontrolled). A merged envelope's margin therefore
STAYS merged (nearest isochore datum is 0). `NaturalDip` opts back into the
legacy unbounded extension for callers who KNOW the regional dip continues.

**Repair precedence — mapped over derived (2026-07-04, `task_petekstatic_topsonly_envelope`).**
The downward repair above (base yields, top preserved) is the *mapped-vs-mapped*
rule. In a horizon stack the repair also distinguishes **mapped vs derived**
(tops-only) surfaces: a mapped horizon is measured and authoritative, so when a
*derived* surface crosses a *mapped* one the **derived surface yields** — the
`Surface::guard_above` / `repair_min_thickness_from_below` twins pull the derived
**upper** *up* to the mapped lower, leaving the mapped surface bit-unchanged. A
mapped-vs-mapped (or derived-vs-derived) crossing keeps the historical
lower-yields direction. With the build-down isochore construction (§4a) a
consistent stack never crosses, so this precedence only guards inconsistent inputs
or a `min_thickness` floor.

### 3b. Areal trend hook — external-drift-lite (Fix 2; INTERIM)

`with_areal_trend(TrendSurface)` (on both the builder and the template) applies a
gridded areal multiplier field, nearest-node resampled to the model column
lattice and **mean-normalized**, per-column to NTG (and φ if the trend flags it)
after population. Semantics (ratified logs-give-shape / draw-gives-level,
`decision_staticmodel_regen_seam`): the **trend gives lateral shape, the
prior/draw gives the level**; because multipliers are mean-normalized the
property field-mean is preserved. Undefined (`NaN`) nodes fall back to `1.0`.
Resampling is index-fraction by default; `TrendSurface::with_georef(origin,
node_dx, node_dy)` maps the field to world coordinates so population resamples by
each column's **world** `(x, y)` (from the grid cell centroids) — the trend then
lands where its data is instead of being aligned to the model lattice by luck
(R4). **This is the INTERIM trend mechanism** (coordinator design
`python-model-build-api.md`): its successor is the per-property pipeline (§3c),
where trends enter as collocated cokriging inside propagation. The hook remains
for the deterministic `with_areal_trend` path; its resampler is now the shared
`petektools::resample` kernel (the nearest-node sampler was retired 2026-07-04).

### 3c. Per-property geostatistical pipeline (P5, `task_petekstatic_property_modelling`)

Properties are modelled **one at a time** through a visible pipeline (coordinator
design `python-model-build-api.md`, owner-refined):

```
PropertyPipeline::new("PHIE")
    .upscale(wells, UpscaleMethod::Arithmetic)   // logs -> per-cell conditioned values (+ QC)
    .propagate(Gaussian::new(variogram, seed))   // SGS fills every cell, conditioned on the upscaled cells
```

- **upscale** is a first-class inspectable step: each positioned `WellLog` snaps
  to its areal column, and its samples that fall in a cell's depth range are
  upscaled (`srs-petro` power means) into that cell — a per-cell conditioned field
  (`NaN` where no log passes) plus an `UpscaleQc` report (upscaled-vs-log stats).
  The cell's depth range is interpolated **at the well's `(x, y)`** (bilinear over
  the cell corners), not the column-centroid 4-corner mean, so an off-centroid well
  on a dipping zone boundary bins each sample into the zone its depth truly falls in
  — without this, near-boundary samples were mis-assigned across the zone and the
  per-zone upscale compressed the zone proportions toward the mid-range.
- **propagate** runs petekTools' sequential Gaussian simulation
  (`geostat::sgs`) **per k-layer**, conditioned exactly on that layer's upscaled
  cells, so every cell is filled with a seeded, reproducible draw that honours the
  wells and reproduces the data histogram. `Gaussian::with_trend(surface, corr)`
  folds a trend in as a collocated (Markov-1) secondary via cokriging, the surface
  resampled to the model lattice through `petektools::resample`; `corr = 0`
  recovers plain SGS bit-for-bit. The moving-neighbourhood search defaults to a
  **bounded** window (16 nodes within `max(1.5·range, 4·spacing)`) — the old
  whole-grid default was pathologically slow on real lattices, and beyond a range a
  node's kriging weight is ~0, so the bounded field matches within simulation
  tolerance; `Gaussian::with_unbounded_search` restores the old window. A simulated
  layer with **no conditioning data** is a hard error naming the property (and, via
  the zone-scoped callers, the zone) — a silent constant mean-fill erased spatial
  structure; `Gaussian::allow_mean_fill` opts back into the fill. **Frame seam:** a
  **world**-georeferenced trend is resampled at
  each column's world position through the model `Georef` (the grid geometry is
  local; the georef maps it to world) — resampling a world trend onto the local
  lattice yields an all-`NaN` secondary that the kernel silently drops to plain SGS,
  so a secondary covering **< 50%** of the model frame is a hard error, never a
  silent no-op.
- Wired into `StaticModelBuilder::with_property` and `StaticModelTemplate::
  with_property` / `with_property_mode`; each run's `PropertyReport` (upscale QC +
  `propagated`) lands on `Provenance.property_reports`. Regular / axis-aligned
  column lattices (the `layer_grid` box/conformable grids); rotated pillars are
  future work.

**MC composition (`decision_mc_composition`).** Per property, the template chooses
a mode:
- **`LevelShift`** (default) propagates once (cached at the first realization) and
  reuses the pattern each draw, adding only the draw's per-property level shift
  (`RealizationDraw::with_property_shift`) — same spatial pattern, moved level; the
  cheap ~ms-class MC path (level uncertainty dominates this asset class). The three
  **fraction cubes** (`PORO`/`NTG`/`SW`) are **shift-then-clamped** to `[0,1]` per
  cell (F9, 2026-07-04): a conditioned real cube legitimately holds boundary cells
  (a non-net `NTG=0`, an aquifer `SW=1`), so a shift that would push them out of
  range **saturates** them instead — an `SW=1` cell stays 1 under a positive shift
  (physically right) and slides interior under a negative one; interior cells move
  by the full shift; `NaN` (undefined) cells are untouched. This is a per-cell
  *application* clamp only — H2 still validates the **drawn inputs** at the seam
  (`realize`), so a garbage sampler is still a typed error. Non-fraction cubes
  (e.g. permeability) are shifted without a clamp.
- **`Resimulate`** re-runs SGS with a fresh per-draw reseed (from `seed_index`) —
  a new pattern each draw (heterogeneity uncertainty), a simulation per draw.

Both are bit-reproducible for identical draws. Bench (50×50×20): `LevelShift`
realize ~320 µs, `Resimulate` realize ~139 ms.

---

## 4. The zones concept

Zones are the **stratigraphic index** of the model. A **zone** is a named
vertical interval of the column, bounded by two framework horizons (or the
top/contact), spanning a contiguous range of `k`-layers — the grid addressable
by geology instead of raw cell index.

Landed shape (`srs-model`): `Zone { name, top_horizon, base_horizon, k_range }`
+ `ZoneTable` (ordered top→base; k-ranges partition `[0, nk)`), with
`ZoneTable::single(nk)` the degenerate whole-column zone and
`ZoneTable::from_stack(...)` the real multi-zone table.

### 4a. The multi-zone regional framework (`task_petekstatic_multizone`, LANDED)

The canonical framework is an **ordered horizon stack** — `N` horizons top→down
define `N − 1` named zones (intra-zone splits are simply more horizons):
`StaticModelBuilder::from_horizon_stack(HorizonStack, BuildOpts)`. A
`HorizonStack` carries the ordered `horizons` and the per-zone `zone_layers`
(`StackZone { name, color, conformity, nk, contacts }` — each zone names itself and
carries an optional display colour; the old `zone_names` parallel array is folded
into `StackZone::name`, and `color` flows to the model's `Zone` and the section
bundle's colour-by-zone `zones` list).

- **Build-down via non-negative CUMULATIVE isochores, anchored at the TOP
  (2026-07-04, `task_petekstatic_topsonly_envelope` + the structure-build rider).**
  The stack is resolved by **constructing downward**, not by gridding each horizon
  independently: the **top** mapped horizon is gridded once; each deeper **mapped**
  horizon carries a **cumulative isochore vs the TOP anchor** — thickness samples
  `(z_k − z_top).max(0)` at co-located defined nodes, min-curvature gridded,
  clamped `≥ 0`, tapered per the extrapolation policy, then made **monotone** down
  the stack (`cum_k = max(cum_{k−1}, iso_k)`); the horizon sits at `top + cum_k`.
  Ordering therefore **can never invert**; where two mapped horizons merge the
  input thickness is exactly `0` and the zone collapses to **genuine zero**; and —
  the **internal-swap invariance convention** — each mapped horizon depends only on
  the TOP grid and its OWN grid, so replacing an internal mapped horizon leaves
  every other horizon (and the total envelope GRV) **bit-unchanged**, except where
  the swapped data genuinely crosses a deeper horizon (the monotone clamp then
  collapses that zone toward the data — never adds rock). Fallbacks when a horizon
  shares no defined node with the TOP: co-locate with the nearest shallower mapped
  horizon (samples lifted by its cum), else an independent solve differenced
  pointwise vs the horizon above, clamped `≥ 0`. **Well ties substitute the
  measured top into the horizon's gridded datum at the tie node and re-run this
  same resolution** — a tie on a fully-defined lattice moves exactly the tied node
  (radius of influence 0 cells; audit S4). *This changed built geometry for
  stacked models (mapped horizons now match the input scatter within gridding
  tolerance rather than being solved in isolation).*
- **The engine owns scatter gridding** (`task_petekstatic_facade_engagement`): a
  **Scatter** horizon (`HorizonSource::Scatter(Vec<WorldPoint>)`, world coords) is
  gridded **inside** the stack build — the single scatter-gridding authority — via
  `from_scatter_stack(stack, opts, StackFrame)`. The engine conditions the points
  onto the model lattice with petekTools **bilinear** off-node conditioning
  (`StackFrame.georef` maps world→node), leaving genuine data voids **`NaN`**; the
  converged solve + `ExtrapolationPolicy` + isochore build-down then act on the
  actual observations, so a data-void margin between two exactly-merged horizons
  **collapses to zero** rather than inheriting an independently-extrapolated fill.
  No caller pre-grids scatter onto the lattice. `Mapped` remains the **pre-gridded
  escape hatch** and is documented as **bypassing** the engine's solve/conditioning
  fidelity (loaded grids only).
  - **Direct-solve conditioning + the factor-once seam** (`task_suite_scatter_perf`):
    per-horizon conditioning uses petekTools' direct band-LU `MinCurvatureOperator`
    (the cap-bound iterative SOR is the degenerate-system fallback), driven through a
    crate-private `ScatterConditioner` that splits **factor** (assemble + factor the
    operator for a horizon's fixed off-node sample geometry — the dominant, once-per-
    surface cost) from **resolve** (back-substitute a fresh depth vector, ~O(n·bw)).
    Because the sample `(x, y)` — and thus the factorization — is fixed across depth
    vectors, re-seating the same geometry with new depths is a cheap resolve, not a
    re-conditioning. Conditioning is **bit-identical** to the prior one-shot path.
    The 11-horizon `condition_scatter` map stays rayon-parallel (each horizon factors
    independently). *The per-draw resolve lever has no live consumer while the MC
    template conditions each horizon once and perturbs surfaces by additive field
    (not by re-conditioning perturbed scatter); a data-perturbation MC mode that uses
    it is a coordinator decision.*
- **Three horizon source kinds** (`HorizonSource`): the raw **Scatter** above; a
  **Mapped** depth surface already on the lattice (tied or untied — the tie is a
  well-tie concern, not a source distinction); and a **TopsOnly** internal horizon
  with *no* mapped surface, defined solely by well **picks**. A tops-only internal split is
  **subordinate to its mapped envelope**: it is seated at the mapped horizon above
  it **+ `min(pick isochore, envelope isochore)`** — its own pick-thickness field
  (min-curvature gridded, constant-thickness fallback for a single pick), plainly
  **clamped** (a `min`, not a redistribution) against the zone's mapped envelope
  thickness. Consequences by construction: envelope `0` ⇒ split `0` ⇒ **both
  sub-zones collapse** to zero exactly (no phantom slab); the split lives strictly
  inside `[mapped top, mapped base]`, so a **derived surface can never displace a
  mapped one**. A *trailing* tops-only horizon with no mapped horizon beneath it
  has no envelope and falls back to the legacy absolute drape. *The geological
  statement: an untied internal split follows the mapped horizon above it at the
  thickness its well picks record, but is never allowed to breach the measured
  horizons that bound its zone.* (Redundant mapped surfaces aliasing one horizon —
  the same `GriddedDepth` given twice — is tolerated; it yields a zero-thickness
  zone by construction.)
- **Per-interface order-repair (precedence: mapped over derived).** On a
  consistent isochore-built stack no interface crosses, so the repair is a **no-op**
  safety net; it still fires under a `min_thickness` floor or genuinely inconsistent
  inputs. When it does, it honours **mapped-over-derived precedence**: a *mapped*
  horizon is authoritative, so a crossing *derived* (tops-only) surface **yields**
  (is moved), never the mapped one; a mapped-vs-mapped (or derived-vs-derived) tie
  keeps the historical behaviour (the lower yields, top preserved). Per-interface
  repair records land in provenance (`StackProvenance::interface_repairs`).
- **Per-zone layering + conformity.** Each zone is layered by its own
  `Conformity` + `dz` (the conformity engine, made per-zone); total `nk` = sum of
  the per-zone counts. The **200-layer cap applies to the TOTAL** (per-zone
  breakdown retained; proportional scale-down when exceeded).
- **Per-zone contacts (a domain statement).** The real framework assigns
  **separate** contact sets to different zones — distinct accumulations in one
  stack — and some zones have **no** contact at all. So contacts are scoped per
  zone: `in_place_by_zone()` clips each zone's hydrocarbons against *its* contacts
  (two-contact gas/oil split where a zone has both), and a **contactless zone
  contributes its gross bulk volume but ZERO hydrocarbon in-place** — no contact
  means no *known* accumulation, so it is explicitly **not** a full hydrocarbon
  column. The rollup `ZonedInPlace::total` sums the zones, so
  `sum(zone volumes) == total` to FP tolerance (conservation).
- **Per-zone stats.** `StaticModel::zone_stats(property)` reports per-zone
  active-cell count/mean/min/max over the full-grid cubes (cubes stay full-grid,
  indexed by zone). Per-zone priors / per-zone SGS variograms are the noted
  follow-up; MC draws stay global-per-property this wave (per-zone contact draws =
  the documented follow-up).
- **Provenance records the stack.** `Provenance::stack: Option<StackProvenance>`
  carries the ordered horizon names, the per-zone layering (`ZoneProvenance`:
  name, bounding horizons, conformity, nk, k_start, truncated_cells), and the
  interface repairs.

The classic Top+Base wireframe (`from_wireframe`) is the 2-horizon / single-zone
degenerate case — unchanged and backward-compatible.

### 4b. Cell-collapse threshold (Petrel-style, LANDED)

`with_collapse_below_m(f64)` (builder **and** template; default OFF): after
layering, any cell thinner than the threshold **collapses** — but
**volume-conservingly**. The layer interface is snapped so a **thicker
zone-interior neighbour absorbs the sliver's thickness** (merge into the thicker
vertical neighbour; rock is never deleted). The merge **never crosses a zone
boundary** (a zone's bounding interfaces are pinned), so each zone's volume is
individually conserved; a single-layer zone has no interior neighbour and is left
untouched. Collapse is applied per pillar, so a cell zeroes only where all four
pillars collapse it (the same max-pillar rule as truncation). Per-zone collapsed
counts land in provenance (`BuildWarning::CellsCollapsed`); a collapsed
(zero-thickness) cell is excluded from volumetrics and `NaN`-marked in the view
bundles by the existing `dz ≤ ε` rule. It bites hardest on sub-`dz` slivers under
fixed-count proportional layering in thin columns and on Follow-style truncation
partials.

---

## 5. Framework, faults, and the wireframe

The `Wireframe` is the structural skeleton the grid conforms to. Today it carries
boundary + horizons + contacts; **faults are deferred** (`question_gridder_spec` /
P5 `task_petekstatic_zones_faults`: split pillars + throw + non-neighbour
connections). The `StaticModel` framework is designed so a fault set is *additive* —
it lands as a new framework member and a grid that splits pillars, without changing
the property-cube or zones contracts. Open/closed: extend by adding, not editing.

---

## 6. Provenance — the reproducibility record (LANDED)

A `StaticModel` records **what produced it**: the input bundle identity
(`inputs_ref`), the gridder settings (`SolveOpts`, `Conformity`, `nk`), the
population mode (`Priors` | `Logs`), and — for a Monte-Carlo realization — the
full **`RealizationDraw`** including its `seed_index`. Provenance is metadata:
it never affects geometry or cubes, but it lets a consumer label a realization
and answer "which inputs gave me this P10 model?". `Hardness` per constraint
(on the wireframe) is the per-node half of this lineage.

---

## 7. The regeneration seam — RATIFIED (graph `decision_staticmodel_regen_seam`)

### 7a. Monte-Carlo regenerates the model per realization

The MC loop does **not** perturb a shared grid in place — it **regenerates a
whole `StaticModel` per realization** from that realization's sampled inputs.
petekStatic is the single owner of *how a model is built*; each realization is a
clean, independently-valid `StaticModel`. (Since the charter re-scope the MC
driver itself is also petekStatic's — `task_peteksim_mc_structured`, retargeted
— which keeps the whole loop below the product facade.)

**Warm-start makes regeneration cheap.** The reusable vs per-realization split:

| reusable across realizations (in the template) | regenerated per realization |
|---|---|
| areal lattice (dims fixed) + conformity (+ `nk` for Proportional) | sampled scalar draw (area, gross height, φ, Sw, NTG, contact, optional GOC via `with_goc`); under a Follow conformity `nk` is dz-derived per realization from the drawn thickness (deterministic in the draw → still bit-reproducible) |
| control-point **topology** (which nodes are pinned) | control-point **depths** (structural shifts, typed `Option`) |
| the previous realization's **converged surface** (the warm seed) | the re-converged surface (a few sweeps from the seed vs 20k cold) |
| attached logs / population inputs | the re-layered grid + repopulated cubes |

Landed as `StaticModelTemplate::new(wf, opts)` (holds the fixed geometry +
warm-start state) with `realize(&mut self, &RealizationDraw) ->
Result<StaticModel>` the cheap call. Measured warm-start win: **~14×** (94 ms vs
1.34 s at 50×50) — the number that makes per-realization regeneration viable.

**Stack-aware MC (P8, `task_petekstatic_multizone_2`).**
`StaticModelTemplate::from_horizon_stack(stack, opts)` runs MC over the multizone
framework: the horizon surfaces + interface repair are resolved **once** (they do
not vary per draw), so `realize`/`realize_into` vary only spacing, the **per-zone
contacts**, and the **per-zone property levels** (`RealizationDraw.zones:
Vec<ZoneDraw>`). Geometry topology is draw-invariant ⇒ bit-deterministic by
construction and `realize_into`-recyclable (the stale-buffer bit-compare holds); a
contactless zone contributes GRV, zero HC.

**The ratified amendments (all landed):**
1. **`StaticModelTemplate: Send` + `StaticModel: Send`** (Sync not required),
   compile-checked via `static_assertions` — consumers shard realizations
   one-template-per-worker; the warm chain makes `realize` serial per template.
2. **`RealizationDraw`**: `#[non_exhaustive]` + `::new(...)` (P5 structural
   fields stay additive, no re-ratification); concrete `pub seed_index: u64`;
   derives `Clone + Debug`; structural perturbation as a concretely typed
   `Option<StructuralPerturbation>` (empty for MVP). **`fvf` deliberately
   EXCLUDED** — FVF is drawn separately and applied at the volumetrics output
   surface, never rides the draw.
3. **Perf:** grid-per-realization is its own (heavier) budget — the 422µs/100k
   figure is the sampler's only. Only property-cube `Vec<f64>`s conceptually vary
   per draw; room is left for an additive
   `realize_into(&mut self, &draw, &mut StaticModel) -> Result<()>` buffer-reuse
   variant (not yet implemented).
4. **`&mut self` on `realize`** (the warm-start chain) and **`Result` per draw**
   (warm-started solves can fail to converge; `StaticError` composes into
   `SrsError` via `#[from]`).

### 7a-bis. Per-horizon correlated structural uncertainty (LANDED — `decision_structural_uncertainty_isochore`)

The MC template's structural draws are **correlated perturbation fields**, applied
**in the space the geometry is built in** so ordering + exact merges survive every
draw *by construction* (no repair rules):

- **TOP surface** takes a correlated **DEPTH** field per draw
  (`RealizationDraw::top_structural`): an unconditional Gaussian random field
  (petekTools `sgs_unconditional`, marginal `N(0, sd_m²)`, the row's variogram)
  added to the converged top. The warm-start chain still advances on the
  **unperturbed** surface (per-draw noise, not structural drift).
- **Every deeper horizon** perturbs via its zone's **ISOCHORE (thickness)** field
  (`ZoneDraw::isochore_structural`): the fixed base thickness plus a correlated
  field, **clamped `>= 0`** and **zero-masked where the base isochore is exactly 0**.
  Deeper horizons are reseated `top + Σ isochores`, so a perturbed stack can never
  invert or resurrect a collapsed zone.
- **Determinism + ties:** each field's seed = `draw.seed_index` salted by the
  horizon index (bit-reproducible per seed; independent across horizons). The stack
  structural `realize` is a **pure function of the draw** → sharded == serial and
  `realize_into` recycling stay bit-exact. Perturbation is pinned to **zero at
  well-tie nodes** (radius-0 locality).
- **Clamp bias (measured + bounded):** the `>= 0` clamp biases mean thickness up near
  pinchouts by `E[max(0, T+P)] − T = sd·φ(T/sd) − T·Φ(−T/sd) ∈ [0, sd/√(2π)]`,
  variogram-independent, → 0 as `T/sd → ∞`. Pinned by a planted-truth GRV recovery
  test (`tests/structural_uncertainty.rs`).

### 7b. Kernel constraint (from the ConvergentGridder adoption — ENFORCED BY TYPE)

The template's warm-start surface lives **in petekTools kernel space**:
`solve_surface_seeded` accepts only the `KernelSurface` newtype, whose only
sources are `KernelSurface::flat` (a constant field — a kernel fixed point, the
safe bootstrap) and a prior `solve_surface_seeded` output.

**Kernel unification landed (2026-07-04, petekTools f81b6a6).** The newtype's
original justification — warm==cold holding only *within one kernel* because
petekTools lacked the natural-dip (linear-extrapolation) boundary of the cold
`solve_surface` (12.48 ft interior sag vs 0.0 on the plane reference) — is
**retired**: petekTools now carries the natural-dip boundary, so the cold reference
and the warm kernel share one fixed point. This is gated **per-node** (not
aggregate — the old min/max/mean cross-check passed while per-node drift was
12 ft): `srs-gridder/tests/adoption_readiness.rs` now asserts per-node plane
parity, and `seeded_warm_equals_cold_of_same_kernel` re-verifies warm==cold to
1e-6 on a well-determined framework (`decision_gridder_kernel_unification`).

**The newtype is KEPT** (reassessed 2026-07-04) on its second, still-live
justification: it is a **provenance / staleness guard** on the warm chain. A seed
must be a genuine same-kernel output — never an arbitrary `Surface` (a cold solve
carries caller-chosen `SolveOpts` whose fixed point can differ from the kernel's
fixed defaults; a hand-built/foreign field is stale). The barrier is zero-cost
(open/closed) and keeps the chain sound; it intentionally has no `From<Surface>`.
A subtlety surfaced while re-verifying: a **single** flat-bootstrap solve is not
itself converged on a *sparse*-control lattice (the smooth level mode relaxes
slowly under the near-Neumann boundary), so `build()` and `realize(mean draw)`
can differ there — tracked as an R2 follow-up, not a natural-dip regression.

### 7c. What the consumer reads off a `StaticModel`

The read-only accessor contract: `model.grid()` (geometry + `bulk_volume` +
`cells()`), `model.property("PORO")` / `"SW"` / `"NTG"` + `property_names()`,
`model.contacts()`, `model.zones()`, `model.provenance()` — **plus the
volumetric answer itself**: `model.in_place()` → `InPlace`, FVF applied by the
caller (`ooip_sm3(OilFvf)` / `ogip_sm3(GasFvf)`). petekSim's facade presents
these results; it re-computes nothing.

### 7d. The structured MC driver + tornado (LANDED, `mc` / `tornado`)

The uncertainty loop over the seam (`task_peteksim_mc_structured` /
`task_peteksim_tornado`). `run_structured_mc(&mut template, &McInputs, n, seed)`
samples every uncertain input from a petekTools `Sampler` (raw or `.clamped` —
petekTools' hard limiter, so a wild fraction prior can't trip §7's H2 guard),
realizes each draw, and **keeps** the per-draw `oil_sm3`/`gas_sm3`/`grv_m3`
vectors (W17 — the realizations are the product; tornado and cross-checks reuse
them). `boi`/`bgi` are sampled but applied at the volumetrics surface, never on
the draw (amendment 2). Reproducible: one seeded stream, fixed field order,
`seed_index = seed + i` per draw.

**Error policy — fail-fast with the draw index.** The first draw whose
realization or volumetrics fails stops the run and surfaces as the typed
`StaticError::McDraw { index, source }`: the §7 H2 rejection is carried to the
caller *with the offending draw identified*, deliberately over collect-and-report
(a bad draw means the input distribution strays outside the physical range — fix
it at the sampler). `aggregate_field` brackets a multi-segment field between
petekTools' `Correlation::Independent` (narrow) and `Comonotonic` (wide downside).

`tornado(&mut template, &McInputs, n, seed, lo_pct, hi_pct)` swings each input
one-at-a-time between its **realized** percentiles (others at P50), re-realizes,
and ranks the oil-in-place swings descending. Pivots are percentiles of the
realized input vectors — not the analytic distribution — so the ranks stay
consistent with the MC P-curve; the same `(n, seed)` reuses the MC's draws.

### 7e. View bundles — the viewer export seam (LANDED, `view`; graph `decision_viewer_home_product`)

The model exports typed, JSON-stable **inspection bundles** for a *separate*
viewer codebase (the `peteksim` viewer). The division is strict — **petekStatic
exports; the viewer renders, it never computes.** Everything a map /
cross-section / 3-D view needs is pre-computed here off the populated
`StaticModel` and handed over as serde-`Serialize` value types. The viewer scales
for display (vertical exaggeration, colour ramps, threshold filtering) but runs
no reservoir computation.

Three bundles, one shared frame:

- **`MapBundle` = `model.map_bundle(&MapSpec)`** — areal (plan-view): realized
  structural depth surfaces (top + base, grid-georeferenced, named after the
  framework horizons), property maps as a single **k-slice** and the
  **zone/interval average** (the useful default), the outline ring(s) in world
  coords (a `from_horizon_stack` build with no wireframe boundary derives the world
  outline from its georef, or takes an explicit `with_boundary` ring — never the
  degenerate unit square against a world frame), well surface markers (with tie
  residuals where provenance carries them
  — none in the MVP), and per-contact **subcrop masks** (the columns the contact
  plane crosses — the honest, marked-mask form of a contact contour). Every layer
  carries its value range for a legend.
- **`IntersectionBundle` = `model.intersection_bundle(&SectionSpec)`** — a
  vertical section walked through the areal lattice (§7e below in API). Since v4
  it also emits `horizon_traces` — one depth polyline per *interior* framework
  horizon (zone boundaries between the structural top/base), parallel to `columns`.
- **`VolumeBundle` = `model.volume_bundle(...)`** — the corner-point cell mesh
  (this is where the mesh builder lives now — moved DOWN from petekSim's
  `srs-core/mesh.rs` per the layer charter; peteksim retires its copy onto this
  bundle at the viewer wave).

**Contract conventions.** All areal layers share one regular, orientable
`GridFrame` (origin + spacing + `ncol=ni` / `nrow=nj` + intrinsic rotation/J
handedness), reconstructed from the grid's column-centroid lattice (the `xy↔ij`
map — the same frame the property pipeline uses). `rotation_deg` is
counter-clockwise from world +X/east to positive I; `yflip` reverses positive J.
These are intrinsic model orientation, independent of viewer camera rotation;
row-major `values[j*ncol + i]`; `f64::NAN` = undefined/outside (JSON
`null`). SI units, positive-down depth. The serialized **key structure** is the
cross-codebase seam — documented in `API.md`, version-tagged (`SCHEMA_VERSION`),
and locked by a schema-snapshot test per bundle.

**The world/local frame discipline (a standing seam rule).** The canonical
real-data configuration is a **LOCAL-origin cell lattice + a registered world
`Georef`** (the "F5" convention): a trace/spec arrives in **world** coordinates
(the view `GridFrame`), but a cell's ZCORN `corners` live in the grid's **local**
lattice. Any geometry that consumes both a spec point *and* the raw corners —
notably `intersection_bundle`'s `fence_edge_depths` areal clip — **must convert to
one frame first** (map world through the oriented frame's exact inverse and then
into the local lattice, or map the corners into world). Local grid/population
kernels stay axis-aligned; orientation exists only at the world seam. Mixing frames does not error;
it silently produces garbage — a world line misses every local cell rectangle and
the section collapses to a flat centroid trace. This has now bitten **three** times
(the bundle-frame F5 zero-columns bug, the collocated-trend no-op, and the
section-edge centroid collapse), so it is a hard convention, not a caution.

  - **A miss is loud, not silent.** Where a frame-consistent computation cannot
    legitimately miss (e.g. a non-degenerate fence line through the cell it was
    matched to), the miss path is a `debug_assert` + a visible warning, not a quiet
    fallback. The single honest degenerate that keeps the centroid (`l == r`) is a
    zero-areal-extent trace (a single-point / vertical bore) — explicit and distinct.
  - **Test convention (mandatory): every frame-sensitive view/population test must
    carry a world-georeferenced variant** — a local-origin lattice + a world georef,
    spec in world coordinates, asserted against the local-frame image of the same
    geometry. An axis-aligned world fence is *not* sufficient (its single bounding
    slab clip is offset-invariant and hides the bug); use an **azimuthed/diagonal**
    trace. See `view::section` tests
    (`georef_world_fence_and_bore_edges_are_frame_invariant_and_dip`).

---

## 8. Contract history — ratified calls (was: open questions)

- **Q1 — Volumetrics home: RESOLVED, INVERTED (2026-07-03).**
  `decision_layer_charters` supersedes `decision_volumetrics_home_peteksim`: the
  static layer owns volumetrics AND static uncertainty (the industry split).
  `srs-volumetrics` + `srs-uncertainty` relocated here; `StaticModel` owns its
  volumes (§2); petekSim keeps `srs-pvt` for its dynamic work and the product
  facade. FVF crosses as a validated scalar input.
- **Q2 — The regeneration API: RATIFIED + LANDED** (§7a;
  `decision_staticmodel_regen_seam` with peteksim's four amendments).
- **Q3 — Draw construction: RESOLVED.** petekStatic defines the neutral
  `RealizationDraw`; whoever owns the sampler fills it (today that is also
  petekStatic — `srs-uncertainty` — but the type stays sampler-agnostic; no
  sampler dependency in `srs-model`'s seam types).
- **Q4 — Faults / zones depth (P5) do not change this contract.** Confirmed:
  the framework/zones/cube contracts absorb faults + multi-zone layering
  additively (`task_petekstatic_zones_faults`); `RealizationDraw` absorbs
  structural sampling via its `#[non_exhaustive]` + typed-Option design.

Kernel unification is closed (`decision_gridder_kernel_unification` /
`task_petektools_natural_dip_boundary`); §7b retains the type guard for
provenance/staleness rather than differing kernel fixed points.

---

## 9. What exists today vs planned

**Exists (green):** `error` (`StaticError` + `#[from] GeoError`),
`srs-grid`, `srs-gridder` (`solve_surface`, `solve_surface_seeded` +
`KernelSurface`, `layer_grid`), `srs-petro`, `srs-wireframe`, `srs-data`,
`srs-volumetrics` (in-place + FVF types + range validation), `srs-uncertainty`
(distributions, SplitMix64, fallible P90/P50/P10), `srs-model` (`StaticModel`,
`StaticModelBuilder`, `StaticModelTemplate`, `RealizationDraw`, `ZoneTable`,
`Provenance`). The full pipeline runs here; petekSim's `srs-core` keeps a thin
`RefiningModel` facade over it.

**Exists (green), continued:** the per-property geostatistical pipeline (§3c,
`task_petekstatic_property_modelling`) — `PropertyPipeline` (upscale + per-layer
SGS + collocated cokriging) on both the builder and the template, with the
`LevelShift`/`Resimulate` MC modes (`decision_mc_composition`); the structured MC
driver + tornado (§7d, `task_peteksim_mc_structured` / `task_peteksim_tornado`) —
`run_structured_mc` / `McInputs` / `McResult` / `aggregate_field` / `tornado`.

**Remaining:** the P5 faulted-gridding remainder (split pillars + throw + NNC;
`task_petekstatic_zones_faults`) and the explicit out-of-core deferrals in §10.
Boundary rings, multi-zone/per-zone population, `realize_into`, the PyO3 facade,
and the first release have landed.

---

## 10. Out-of-core — the backing-storage mode (`task_petekstatic_slab_streaming`)

Larger-than-memory models must build, populate, run MC, and answer volumetrics by
**spilling to disk**, not by dying (the P9 engine half; coordinator design
`petekSuite/dev-docs/designs/out-of-core-strategy.md`, rulings R2/R3/R4/R5). The
design is **budget-driven and loud**, and — critically — the in-core StaticModel
type keeps its API: **spilled is a backing-storage mode, not a new type.**

**The budget + the loud switch (R5).** The build/template API carries a declared
`MemoryBudget` (default: a documented fraction of physical RAM). `build()` compares
a live-set estimate (`live_set_bytes`: ZCORN f64 + cubes f64, ×2 for the warm/MC
state) against it. **Below → the in-core path, byte-identical** (zero behaviour
change; the golden/determinism suite is untouched and a generous-budget build is
bit-identical to the pre-budget one). **Above → the spilled path**, and the engine
**says so**: a `SpillNotice` on stderr naming the mode, budget, estimate, and store
path. Never an OOM kill; never a silent switch.

**The abstraction (the enum, cleanly).** A `StaticModel` carries an internal
`spill: Option<Arc<SpillBacking>>`. `None` is the in-core model — geometry + cubes
owned by the `Grid`, every accessor byte-identical to before. `Some(_)` is a spilled
model whose heavy per-cell arrays live in a memory-mapped petekTools `store` (the
R1 lane store); its `grid` is a unit placeholder. The volumetric surface
(`in_place`, `in_place_summary`, `bulk_volume`, the two-contact split) **dispatches
on the backing**: in-core reads the `Grid`; spilled **streams** the store one k-slab
at a time. Shared behind an `Arc` so the model stays `Clone`/`Send` (the mmap is
shared, not copied) and the store file is cleaned up when the last clone drops
(`Drop`), unless the caller detaches it (`with_spill_persist`). `grid()`/`property()`
on a spilled model return the placeholder — a spilled model is read through its
volumetric methods; the raw borrowed `&Grid` is the in-core representation only.

**The lane layout maps the grid's natural order (R1/R2).** ZCORN is already
k-slab-major (`zcorn[k·ni·nj·8 ..]`) and each cube is `cube[k·ni·nj ..]`, so a
k-window of any lane is **one contiguous zero-copy slice** — exactly the store's
chunked-along-k shape. ZCORN → a slab f32 lane; each cube → a slab f32 lane; COORD
(small, k-invariant) → a resident f64 flat lane. The spilled cell rebuild is
bit-parallel to the in-core vertical-lattice fast path except ZCORN is narrowed to
f32 (the sole parity delta).

**Precision + parity (R4, the honesty clause).** Storage lanes are **f32** at spill
scale — the out-of-core enabler and the MC bandwidth lever (halves the bytes on the
streaming summary read). So in-core↔spilled parity is **tolerance-based, not bit
parity**: measured on a 2500 m fixture, GRV/bulk are bit-exact, HCPV ≤ ~2e-9, and MC
oil ≤ ~3.1e-6 relative (larger because a near-contact cell can flip in/out under f32
centroid rounding) — asserted ≤ 1e-5, documented here. **Accumulations stay f64**
regardless; only storage lanes narrow. Bit-determinism holds **within** a mode
(sharded MC == serial at every worker count, in-core and spilled).

**MC never spills per-draw state (R3).** The structured-MC loop keeps its
one-reusable-model / realize → summary → discard shape. In the spilled mode each
draw realizes in-core into the one reusable model, is flushed to **one reused
per-shard store** (overwritten every draw — never a new file per draw, so no
per-draw disk accumulation), its summary streamed from the f32 lanes, then the
backing is detached (the in-core grid stays intact for the next `realize_into`).
`run_structured_mc_spilled` / `_parallel_spilled` prove the determinism +
stale-buffer contracts hold in the spilled mode.

**Slab-incremental build — true O(slab) build peak (v2, item 1).** The gridder
path now builds a spilled model **slab-incrementally**: a `StreamingLayering`
(srs-gridder) produces ZCORN one k-slab at a time from two interface planes (the
same f64 `boundary_depth` as the in-core path, narrowed to f32), and the builder
writes each ZCORN + constant-cube slab **straight into the store**
(`slab_mut_f32`), never materializing a whole in-core grid. The heap/anonymous
working set is **O(slab)** (a slab + two O(area) interface planes), so a model
whose in-core live set exceeds RAM can now be **built**. Bit-identical to
build-then-spill (same f64→f32 narrowing, same constant cubes). Measured
(`dev-docs/bench/scripts/rss_probe.sh`, 3.2M cells, `/usr/bin/time -l` max RSS):
in-core 90 B/cell (the whole f64 grid — unbuildable above RAM), v1 build-then-spill
135 B/cell (grid + store resident together), **v2 streaming 47 B/cell** (≈ the f32
store alone) — a 2.9× cut, and v2 sits *below* the in-core build. Eligibility: the
per-slab-safe population (constant priors / per-zone constant), collapse **off**
(the volume-conserving collapse pass is a whole-zone-column carry). Ineligible
spilled builds (logs/trend/SGS pipelines, or collapse) fall back to build-then-spill.

**v2 deferrals (honest).**
- **True O(slab) *resident* RSS** would need the petekTools `store` writer to shed
  written pages (`madvise(MADV_DONTNEED)`); today the mmap'd f32 store's written
  pages stay in the page cache, so total RSS is O(store_f32) = the on-disk model
  itself (file-backed, sheddable) + O(slab) heap. Raised to petekTools as
  `question_petektools_store_page_evict`.
- **Streaming SGS-pipeline population** (logs/trend/SGS builds spill via
  build-then-spill) and **writable-store MC realize** (spilled MC realizes in-core
  into the one reused model then flushes — R3: MC never spills per-draw state, so
  its peak is the bounded reused model, not a scaling wall) — both need per-slab SGS
  population + must preserve the `realize_into` budget/determinism contracts.
- Spilled `zone_stats` (arbitrary-property per-zone cube statistics) and spilled
  view-bundle windowing. (Spilled `in_place_by_zone` **shipped** in v2 — item 3.)

---

## References (suite planning graph)

- `task_relocate_refine_orchestration` (this wave) · `task_peteksim_mc_structured`
  (retargeted) · `task_peteksim_tornado` (retargeted) ·
  `task_petekstatic_zones_faults` · `task_petekstatic_property_modelling` ·
  `task_petekstatic_boundary_rings`
- `decision_layer_charters` (the 2026-07-03 re-scope) ·
  `decision_staticmodel_regen_seam` (the ratified seam) ·
  `decision_gridder_kernel_unification` · `decision_srs_core_seam_pathdeps`
- Specs: `convergent_gridder_spec`, `layer_interpolation_spec`,
  `log_upscaling_spec`, `data_to_wireframe_spec`, `grv_from_grid_spec`,
  `monte_carlo_volumetrics_spec`, `distribution_sampling_spec`
- Cross-library design: `petekSuite/dev-docs/designs/staticmodel-regen-seam.md`
