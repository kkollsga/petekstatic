# Changelog

All notable changes to petekStatic are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/); this project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.12] - 2026-07-10

### Changed
- The legacy `petekstatic::data` petekIO bridge is now behind the opt-in
  `petekio-adapter` Cargo feature. The default geomodel core and Python wheel no
  longer depend on petekIO; petekSim owns full-stack composition. Existing Rust
  callers can enable the feature during the compatibility window.
- CI now builds the ABI3 wheel once for the Python 3.10–3.14 test matrix and
  explicitly verifies that the default wheel imports and builds a model without
  petekIO installed. Rust gates exercise both the optional adapter and the
  dependency-minimal core.
- Release wheels and the source distribution now build in parallel with the
  unchanged release gates. Trusted PyPI publishing retries transient failures,
  and the workflow reports trigger-to-installable-wheel time.

## [0.1.11] - 2026-07-08

### Changed
- The petekTools floor is `0.2.7` for the topology-aware kernels and corrected
  viewer behaviour used by the geomodel workflow.

## [0.1.10] - 2026-07-08

### Changed
- Updated the Rust and Python dependency floors to `petekio` 0.3.7 and
  `petektools` 0.2.6 so petekStatic consumes the released point-edge semantics,
  structured mesh surface support, and 2-D topology-grid viewer QA.
- Hardened CI/release wheel install checks with retrying binary-only installs
  and added a PyPI visibility verification job before GitHub Release creation.

## [0.1.9] - 2026-07-08

### Changed
- Updated the Rust and Python dependency floors to `petekio` 0.3.6 and
  `petektools` 0.2.5 so petekStatic consumes the released calculated-log,
  interpolation, and standardized object-history APIs.

## [0.1.8] - 2026-07-08

### Changed
- Updated the Rust and Python dependency floor to `petekio` 0.3.5 so
  petekStatic consumes the released project organization and import/load split
  from the DATA layer.

## [0.1.7] - 2026-07-07

### Changed
- Updated the Rust and Python dependency floor to `petekio` 0.3.4 so
  petekStatic consumes the released point-set geometry/edge API from the DATA
  layer.

## [0.1.6] - 2026-07-07

### Added
- Added `HorizonSpec` and `WellTie` to the Python workflow facade. `Grid.horizons`
  now accepts per-horizon mappings such as
  `{"name": "Top reservoir", "surface": "Top reservoir input surface",
  "well top": "well tops/Top reservoir"}` plus
  `well_tie={"influence_radius": 800}`. `name` is the model horizon name;
  `surface` and `well top` bind loaded project inputs.
- Added horizon `zone` tags and inline zone declarations. `zone="Reservoir"`
  means the zone below that horizon; `zone={"name": "Reservoir", "sub-zones":
  [...]}` defines sub-zones directly from `.horizons(...)`. Sub-zone boundary
  names such as `Top Lower Reservoir` are inserted as nested grid horizons.
- Added compact mixed sub-zone sequences such as
  `["Upper Reservoir", {"name": "Intra Shale", "surface": "Top Lower Reservoir"},
  "Lower Reservoir"]`, where the mapping defines the inserted boundary horizon.
- Added sub-zone construction `type` metadata (`constant`, `conformable`,
  `isochore`, `fraction`, or `rest conformable`) and boundary mappings that bind
  to well tops, e.g. `{"name": "Intra Shale", "well top": "Top Lower Reservoir"}`.
- Added `{"name": "Intra Shale", "surface": True}` for nested model-surface
  boundaries without an external input binding.

### Changed
- Removed the old `Grid.horizons(..., tie_to_tops=..., gridding=...)`
  compatibility path and the Python `Gridding` workflow export. Well-top tying is
  now expressed only through `well_tie=...`.
- Updated petekIO/petekTools dependency floors for the release train.
- Updated CI and release workflows to current action versions and the
  Actions-owned release flow.
- Removed the stale Python perf comparison against the old `peteksim.upscale`
  facade; property recipe performance is now gated against petekStatic +
  petekIO directly.
- Tightened the LevelShift hot path by copying cached zero-shift property
  patterns directly and clarified that release perf budgets must run serially.

## [0.1.5] - 2026-07-07

### Fixed
- Declared the Python wheel runtime dependencies on `petekio` and `petektools`
  so isolated wheel installs can execute the workflow facade's log-resolution,
  formula, and variogram paths without relying on a developer venv.

## [0.1.4] - 2026-07-07

### Added
- Added the first petekStatic-owned Python static workflow facade:
  `Grid.from_project(...).geometry(...).horizons(...).zones(...).layers(...)`,
  `Gridding`, `Layering`, `Spherical`, property handles/store, scalar and
  expression assignment, `calc(...)`, and `volumes(...).run(...)` with
  `VolumeResult.summary()` / `by_zone()`.
- `grid.properties.calc(...)` delegates formula blocks to
  `petektools.evaluate_formula` and commits calculated properties atomically after
  the block validates.
- Added canonical property recipe declarations:
  `pst.upscale(logs.PHIE(logs.NetSand > 0.50)).sgs(...)`,
  `pst.Var(...)`, `pst.distributions.from_logs()`, `WellLogSpec`, and
  `PropertyPipelineSpec` lowering/inspection/execution examples.
- Added the first Python execution boundary for property recipes: lowered specs
  with resolved positioned wells, `from_logs` distribution, no trend/cokriging,
  and isotropic variograms can build a Rust-backed `PropertyPipeline`; arbitrary
  production-grid application remains unexposed except for the flat-model smoke
  path.
- Added synthetic Python tests for project asset validation, formula-backed
  calculated properties, progress callbacks, deterministic simple volumes,
  property recipe serialization/lowering, petekTools anisotropic variogram
  lowering, and Rust pipeline handle construction.

### Notes
- The workflow is currently a Python facade/spec layer for notebooks and smoke
  tests. It does not yet lower into the full Rust `StaticModelBuilder` or execute
  production corner-point geometry/property modelling against arbitrary grids.

## [0.1.2] - 2026-07-06

### Changed
- Rebuilt the petekStatic example notebooks around the ratified role-binding
  flow: synthetic project trees now come from `petektools.synth_asset`, synthetic
  generation and `Project.load` live in separate cells, the notebooks print
  `Project.inventory()` before binding, declare replaceable literals for
  outline/horizons/zones/subzones/contacts, build through the `peteksim` facade,
  and keep synthetic manifest checks in final skipped-for-real cells. Notebook 02
  now documents the current per-segment run + `ps.aggregate` pattern instead of
  implying a zone-by-segment single-model rollup.
- Fast-pathed synthetic flat-box builds through the kernel-owned flat fixed point,
  restoring the ignored release forced-spill budget without bypassing the general
  solve for user-added structural controls.

### Added
- Locked the Python wheel export surface with a test asserting
  `petekstatic.__all__ == ["StaticModel", "__version__", "build_flat_model"]`.

## [0.1.1] - 2026-07-06

### Fixed
- **Windows wheel + sdist now published.** The 0.1.0 release matrix built only
  linux (x86_64/aarch64) and macOS (x86_64/arm64) wheels — no `win_amd64` wheel
  and no sdist, so `pip install petekstatic` failed on Windows ("No matching
  distribution found"). The release workflow now builds a `windows-latest x64`
  wheel and an sdist alongside the existing four, and CI gained a permanent
  Windows build + install + smoke job so a missing release target can never go
  unnoticed again. No code changes.

## [0.1.0] - 2026-07-05

First public release: the consolidated `petekstatic` crate + the `petekstatic`
PyPI wheel (abi3-py310).

### Packaging — CONSOLIDATED the 10-crate workspace into ONE crate, `petekstatic` (owner ruling)
- **The ten published workspace crates are now a single crate, `petekstatic`
  (0.1.0).** `petekstatic-error`, `srs-grid`, `srs-gridder`, `srs-petro`,
  `srs-wireframe`, `srs-data`, `srs-volumetrics`, `srs-uncertainty`, `srs-spill`
  and `srs-model` collapse into **modules** of one crate — `petekstatic::{error,
  wireframe, grid, petro, gridder, volumetrics, uncertainty, data, spill, model}`
  — preserving today's boundaries and one-directional imports as module discipline.
- **The `srs-*` crates and `petekstatic-error` are deprecated and will be yanked
  from crates.io.** Depend on `petekstatic` instead.
- **Migration.** `srs_model::X` → `petekstatic::model::X` (or the crate-root
  re-export `petekstatic::X`); `srs_grid::X` → `petekstatic::grid::X`;
  `srs_wireframe::X` → `petekstatic::wireframe::X`; likewise for every module.
  The headline API (`StaticModelBuilder`, `StaticModelTemplate`, the `HorizonStack`
  family, `run_mc` / `McSettings`, `BuildSpec`, `StaticModel`, the view bundles,
  `StaticError`) is **re-exported at the crate root** — reach it without knowing
  the module.
- **The `petekstatic` PyPI wheel is unchanged** — the `petekstatic-py` binding
  crate simply rebinds onto the consolidated crate; the compiled extension
  (`petekstatic._petekstatic`) and its Python surface are byte-for-byte the same.
- **Behaviour-neutral packaging only** — no algorithm, tolerance, determinism seed
  or public-value change. All perf budgets, determinism tests and acceptance suites
  moved with their code and pass unchanged.

### Performance — parallelized per-layer SGS + primary-path structure solves (`task_suite_perf_round2`)
- **Per-layer zone-property SGS now runs in parallel.** `propagate_sgs_into`'s
  per-`k`-layer sequential-Gaussian sweep was serial; every layer is an independent
  2D SGS with a per-global-`k` seed (`g.seed ^ (k+1)·GOLDEN`, range-independent),
  its own conditioning, and a disjoint output slice — so the layers now simulate
  across a rayon `par_iter` and their fields scatter back in `k`-order. **Bit-for-bit
  identical** to the serial sweep (the seed is a pure function of `k`, never of
  execution order). Uses petekTools' new `sgs_seeded` to share one `&SgsParams`
  (the layer-invariant collocated secondary) across layers without per-layer clones.
  On the canonical real-model build the zone-property stage dropped **~3.5 s → ~0.9 s**.
- **Primary-path isochore solves in `resolve_stack_surfaces` now run in parallel.**
  Each mapped envelope's primary build-down solve co-locates vs the top anchor and
  depends only on the top + that envelope's own datums — never on another envelope —
  so the (dominant) warm-start solves are pre-computed across a rayon `par_iter`; the
  serial build-down consumes them and falls back to a serial solve only for the
  co-location-fallback / degenerate cases (which carry the cross-envelope cumulative
  dependency). **Bit-for-bit identical.** Surface resolution dropped **~1.0 s → ~0.25 s**
  per resolve (and the tied re-resolve likewise).
- **Net:** the canonical real-model run (build + deterministic volumes + zoned MC
  n=64 + parity MC + viewer export) went **~20.5 s → ~8.5 s** with bit-identical
  results (deterministic volumes, tie residuals, trend correlations, and per-zone +
  total P-curves all unchanged). Both levers recur across the deterministic build
  and every MC template rebuild, so they compound. No public API change.

### Changed — adopted petekTools' direct-solve `MinCurvatureOperator` in the scatter path (`task_suite_scatter_perf`)
- **The durable order-of-magnitude scatter-conditioning fix is now the kernel's
  default and adopted here.** petekTools shipped the direct band-LU
  `MinCurvatureOperator` (factor-once / solve-many), and its `grid_min_curvature*`
  entries now dispatch to it (the cap-bound ~60 s SOR is the fallback for
  degenerate systems only). So `grid_scatter`'s per-horizon bilinear conditioning —
  which calls that entry — **already runs the direct solve** (the ~60 s → ~0.44 s
  per-horizon win at the canonical ~40k-sample density lands via the dependency,
  no behaviour change beyond the direct solve *attaining* the SOR fixed point
  instead of stalling metres short of it).
- **`grid_scatter` now drives the operator explicitly** through a new crate-private
  `ScatterConditioner` (`factor` → `resolve`), pinning the fast path and exposing
  the **factor-once / solve-many** seam: the sample `(x, y)` geometry — and thus the
  band-LU factorization — is fixed across depth vectors, so re-seating the same
  geometry with new depths is a ~6 ms `resolve`, not a fresh ~0.44 s conditioning
  (the **MC lever**). The one-shot conditioning field is **bit-identical** to the
  prior `grid_min_curvature_conditioned` path (same operator); a degenerate/singular
  system still falls back to the iterative kernel. Proven bit-for-bit that a
  `resolve` on a reused factor equals conditioning the perturbed scatter from scratch
  (`resolve_on_a_reused_factor_is_bit_identical_to_conditioning_from_scratch`).
  Micro-bench (`benches/scatter_condition.rs`, 100×100 nodes / 19.6k off-node pts):
  **factor ~229 ms vs resolve ~3.0 ms** (~76× — the reuse lever).
- **Behaviour / accuracy:** results stay within petekTools' documented
  ~1e-4-of-converged-SOR change and are deterministic (fixed elimination order).
  The only pinned-value test that moved is `structure_fidelity::s2` — it asserted the
  converged wrapper must *beat* the raw flat-seeded kernel's affine-mode **stall**;
  the direct solve **eliminates** that stall (the raw kernel now reaches the fixed
  point in one solve), so S2 now asserts both the raw kernel and the converged entry
  sit at the fixed point (the wrapper agrees with the kernel rather than rescuing it).
- **The per-draw `resolve` lever has no live consumer yet.** The MC template
  conditions each scatter horizon **once** (the condition-once dedup seam) and its
  structural perturbation re-seats surfaces by **additive SGS field**, not by
  re-conditioning perturbed scatter — so nothing re-solves fixed geometry per draw
  today. The seam is in place and proven; wiring a data-perturbation MC mode onto it
  is a coordinator decision (routed via `task_suite_scatter_perf`).

### Changed — scatter-build performance: parallel conditioning + a condition-once dedup seam (`task_suite_scatter_perf`, P9-scale)
- **Parallelized per-horizon scatter conditioning.** `condition_scatter` now grids
  the 11 (N) independent scatter horizons across rayon workers instead of a serial
  loop. Each horizon is an independent, deterministic cold bilinear solve with **no
  cross-horizon reduction**, so the conditioned surfaces are **bit-identical** to the
  serial path (`structure_fidelity` scatter/determinism suite unchanged). A
  scaled-down fixture (3.6k pts/horizon) drops from ~52 s serial to ~14 s. At the
  **real canonical density (~39,668 pts/horizon)** each cold bilinear solve is
  **cap-bound at ~60 s** (burns the full ~20k SOR sweeps, never reaching TOL;
  warm-start does not help) — one 11-horizon pass is ~11 min serial, ~66 s parallel
  on a 10-core box (~10×). Every *resolve* solve is ~1 ms; the entire cost is the
  conditioning. The real `ntg_view.py` canonical run dropped ~43 min → ~29 min on a
  dev box with parallelism alone; the residual is the not-yet-deduped 3×+ redundancy
  (build + MC template + per-MC-worker) and non-conditioning stages.
- **New public dedup seam `StaticModelBuilder::condition_scatter_stack(stack,
  &frame) -> HorizonStack`.** Conditions the raw scatter **once** and returns the
  conditioned all-`Mapped` handle, so a caller building both a `StaticModel` and its
  MC `StaticModelTemplate` from the same scatter feeds the handle to
  `from_horizon_stack` (+`with_georef`) on each path **without re-running the cold
  conditioning** (conditioning is draw-invariant + identical across the two). Proven
  **bit-identical** to conditioning inside `from_scatter_stack`
  (`scatter_dedup_seam_is_bit_identical`). Removes the redundant re-solves that made
  the canonical build re-condition the full stack 3× (build + MC template + parity).
- **Release perf budget** `canonical_scatter_build_within_budget` (≤30 s; measured
  ~14 s) — a regression tripwire that trips if the parallel conditioning is lost.
- The durable order-of-magnitude fix (factor-once / direct sparse solve of the
  biharmonic+bilinear conditioning operator) belongs to the kernel owner
  (petekTools `grid_min_curvature`) and was routed to the coordinator with a precise
  spec; petekStatic never reimplements the kernel (`decision_gridder_kernel_unification`).
  **Now delivered and adopted** — see "adopted petekTools' direct-solve
  `MinCurvatureOperator`" above.

### Changed — organize wave: split `srs-spill` out of the `srs-model` god-crate (`task_petekstatic_organize`, P10)
- **New workspace crate `srs-spill`** — the out-of-core backing-storage mode
  (memory budget + build-mode decision + the k-slab spill onto a petekTools store)
  moved out of `srs-model` into a leaf crate below it (`srs-model` → `srs-spill` →
  `srs-grid`/`srs-volumetrics`/`petektools`). Behaviour is **byte-identical**;
  `srs-model` re-exports the whole surface, so the public path is unchanged
  (`srs_model::MemoryBudget`, `srs_model::spill_grid`, `srs_model::SpillBacking`,
  `srs_model::decide_mode`, … all resolve exactly as before). This removes ~670
  LOC + 2 responsibilities from the top-of-DAG aggregate. Two former `pub(crate)`
  helpers (`unique_spill_path`, `SpillBacking::to_in_core_grid`) are now `pub` on
  `srs-spill` as the cross-crate seam (documented as internal-tier helpers).
- **`tornado` swing loop is table-driven** — the seven hand-unrolled per-input
  `bar(...)` blocks collapse into a single `(name, accessor, setter)` field table;
  bar order, pivots, and every swing value are identical (no behaviour change).
- **`with_areal_trend` documented as superseded/deprecation-tracked** on both the
  builder and the MC template — the interim post-population NTG/φ multiplier hook;
  the fuller path is `PropertyPipeline` + `Gaussian::with_trend` (collocated
  cokriging). Retained (still on the build + MC-template path); no `#[deprecated]`
  attribute yet.

### Added — per-horizon correlated structural uncertainty (`task_petekstatic_structural_uncertainty`, `decision_structural_uncertainty_isochore`)
- **`PerturbationField { sd_m, variogram }`** — a correlated structural
  perturbation FIELD for one horizon/isochore of an MC draw: an unconditional
  Gaussian random field (petekTools `sgs_unconditional`) with marginal
  `N(0, sd_m²)` and the variogram's spatial continuity, generated on the areal
  node lattice at `realize` time. `Copy + Clone + Debug + PartialEq + serde`;
  construct with `PerturbationField::new(sd_m, variogram)`.
- **`RealizationDraw::top_structural: Option<PerturbationField>`** (+ `with_top_structural`)
  — a correlated **TOP-surface DEPTH** perturbation field applied per draw on
  **both** the 2-surface and stack paths. The warm-start chain still advances on
  the unperturbed converged surface (a perturbation is per-draw noise, not drift).
- **`ZoneDraw::isochore_structural: Option<PerturbationField>`** (+ `with_isochore_structural`)
  — a correlated **THICKNESS (isochore)** perturbation field for that zone,
  applied in isochore space so ordering and exact merges survive every draw **by
  construction**: the perturbed thickness is clamped `>= 0` and **zero-masked
  where the base isochore is exactly 0** (a merged/collapsed zone stays merged in
  every draw — a perturbed stack can never invert or resurrect a collapsed zone).
  Because a deeper horizon is `top + Σ isochores`, this is how a deeper horizon
  perturbs structurally.
- **Determinism + ties.** Fields are seeded from `draw.seed_index` salted by the
  horizon index → **bit-reproducible per seed**, mutually independent across
  horizons; the stack structural `realize` is a pure function of the draw, so
  **sharded == serial** and `realize_into` recycling stay bit-exact. Perturbation
  is **pinned to zero at well-tie nodes** (radius-0 locality, matching the default
  `Replace` tie). `sd_m <= 0` is a no-op. The **clamp-induced mean-GRV bias** is
  characterised analytically (`E[max(0, T+P)] = T·Φ(T/sd) + sd·φ(T/sd)`,
  variogram-independent) and pinned by a planted-truth recovery test (pinchout +
  thick-zone bound). `BuildSpec` stays `#[non_exhaustive] + #[serde(default)]`;
  the leg rides the (already `#[non_exhaustive]`) draw types additively.

### Changed — license: adopt Business Source License 1.1 (`decision_license_ratified`)
- Before its first-ever publish, petekStatic adopts the **Business Source License
  1.1** (was `license = "MIT"` in the workspace manifest, with no LICENSE file).
  The BSL text is added as [LICENSE](LICENSE) with `Licensor` / `Additional Use
  Grant` / `Change License = Apache-2.0` filled; the `{VERSION}` and Change Date
  parameters are release placeholders, filled at each release cut. The Cargo
  `[workspace.package]` `license = "MIT"` field is replaced by `license-file =
  "LICENSE"` (BSL has no SPDX expression), inherited by all nine crates via
  `license-file.workspace = true`. Each released version converts to Apache-2.0
  four years after its first publication. Never published under MIT, so no
  legacy-version note is needed.

### Added
- **`BuildSpec` — ONE declarative build configuration for BOTH `StaticModelBuilder`
  and `StaticModelTemplate`** (`task_petekstatic_spec_mirror`, suite
  api-consistency contract). The duplicated ~12-method `with_*` chains are now
  thin sugar mutating one internal `BuildSpec`; `with_spec` installs a whole
  configuration on either consumer. Fields: `inputs_ref`, `georef`, `boundary`,
  `extrapolation`, `clamp_base_to_top`, `min_thickness_m`, `collapse_below_m`,
  `sugar_cube`, `sw_gas`, `well_ties`, `ties` (a `TieSettings`). Values are
  identical either way — the determinism contracts (builder == template, sharded
  == serial, `realize_into` staleness, byte-identity) are pinned unchanged by
  `tests/spec_conformance.rs`. `#[non_exhaustive]` + `#[serde(default)]`, so the
  structural-uncertainty per-zone isochore leg
  (`decision_structural_uncertainty_isochore`) lands additively later.
- **`McSettings { n, seed, workers, spill_dir }` + `run_mc`** — the ONE structured
  MC entry consolidating `run_structured_mc` / `_parallel` / `_spilled` /
  `_parallel_spilled` (all four remain as **deprecated** thin wrappers, pinned
  bit-identical to `run_mc` by test). `workers = 1` is serial; `spill_dir:
  Some(dir)` selects the out-of-core mode (pass `std::env::temp_dir()` for the
  old `spilled(None)` behaviour).
- **`TieSettings { method: Replace | Radius { radius_m } }`** — the settings
  mirror over the datum-substitution well-tie machinery, on the spec
  (`with_tie_settings` on builder + template). `Replace` (default) = today's
  behaviour: the measured top replaces the map datum at the tie node.
  `Radius { radius_m }` = bounded locality: the tie's residual blends into every
  defined datum within `radius_m` with a linear decay (1 at the well → 0 at the
  radius); datums beyond are bit-untouched, undefined nodes stay the solve's to
  taper. One tie authority (`substitute_tie_datums`) serves builder and template
  — pinned node-for-node identical by test. On the template, set tie settings
  BEFORE `with_well_ties` (ties apply at that call).
- **Config-layer serde (R7)** — additive `Serialize`/`Deserialize` on the whole
  declarative config family: `HorizonStack`/`StackHorizon`/`HorizonSource`
  (incl. `Scatter`/`WorldPoint`)/`Pick`/`StackZone`/`StackFrame`/`WellTie`/
  `BuildOpts`/`Georef`/`RealizationDraw`/`StructuralPerturbation`/`ZoneDraw`/
  `BuildSpec`/`TieSettings`/`McSettings`, plus the value types they carry
  (`GriddedDepth`/`Contact`/`ContactKind`/`Hardness`/`Boundary`/`Horizon`/
  `HorizonRole` in srs-wireframe; `SolveOpts`/`Conformity`/`ExtrapolationPolicy`
  in srs-gridder; `ConstantPriors` in srs-volumetrics — `SolveOpts`,
  `ConstantPriors` and `BuildOpts` also gain `PartialEq`). The R7 engine-half
  battery (`tests/spec_conformance.rs`) pins round-trip == value equality on
  every type. **`McInputs` is deliberately excluded** — it wraps petekTools'
  `Sampler`/`Clamped`, which do not yet derive serde; the upstream derives are
  queued (petekTools spec wave), and the battery documents the slot.

- **Section colour-by-zone payload (rider, `task_suite_section_zone_color`;
  SCHEMA_VERSION 4 → 5, additive).** `IntersectionBundle` gains a `zones` list
  (`[{name, color}]`, the model's stratigraphic zones top→base), and each
  `SectionColumn` gains a per-layer `zone_ids: Vec<u16>` — an index into `zones`
  naming the zone each layer belongs to, `NaN`-gapped in lockstep with the
  geometry/value arrays (an inactive/truncated layer carries the sentinel
  `SectionColumn::NO_ZONE` = `u16::MAX`). This is the engine half of the viewer's
  colour-by-zone rendering; both fields are `#[serde(default)]` so older payloads
  decode unchanged. New public `SectionZone { name, color }`.
- **`StackZone` gains `name` + optional `color`** (`StackZone::new` /
  `StackZone::with_color`), and the model's `Zone` gains `color` — the display
  colour flows `StackZone::color` → `Zone::color` → the section bundle's `zones`
  list.

### Changed
- **`ZoneDraw` is now `#[non_exhaustive]`** (construct via `ZoneDraw::new` + the
  `with_*` setters — every in-tree caller already did) so the per-zone isochore
  perturbation leg (`decision_structural_uncertainty_isochore`) can land as an
  additive field.

- **⚠️ BREAKING (`HorizonStack` shape) — the `HorizonStack::zone_names` parallel
  `Vec<String>` is retired; each zone's name now lives on `StackZone::name`.**
  Construct a stack by giving every `StackZone` its `name` (and optional `color`)
  instead of a separate `zone_names` array. `ZoneTable::from_stack` gains a
  `zone_colors: &[Option<String>]` parameter.
- **⚠️ BREAKING (built geometry) — `from_horizon_stack` now BUILDS DOWN via
  non-negative isochores (`task_petekstatic_topsonly_envelope`).** Every model with
  a stacked horizon framework changes geometry. Previously each mapped horizon was
  gridded **independently** and tops-only internal splits were draped as an
  **absolute** pick-thickness offset from the horizon above; where a mapped envelope
  merged (two mapped horizons coinciding — the zone geologically absent), the
  absolute drape crossed the mapped base and the order-repair pushed the **mapped
  base down** to sit below the derived split — a **phantom** near-constant thickness
  (≈ the pick thickness) across the entire merged region, propagating down the
  stack. The build now constructs downward: the top is gridded once, each deeper
  **mapped** horizon = the horizon above **+ a clamped (`≥ 0`) gridded zone
  isochore** (so a merge samples exactly `0` and the zone collapses to genuine
  zero), and a **tops-only** split = the mapped horizon above **+ `min(pick
  isochore, envelope isochore)`** (a plain clamp inside its mapped envelope). By
  construction: ordering can never invert, mapped horizons stay **bit-authoritative**
  (a derived surface can never displace one), and a merged envelope yields **exactly
  zero** sub-zone thickness — no phantom. Mapped-horizon depths now match the input
  scatter within gridding tolerance rather than being solved in isolation. Interface
  order-repair gains **mapped-over-derived precedence** (a crossing derived surface
  yields, never a mapped one) and is a no-op on consistent inputs. New
  `Surface::guard_above` / `Surface::repair_min_thickness_from_below` (the
  repair-precedence twins). See SPEC §4a + §3 (R-c).

- **⚠️ Structure-build fidelity overhaul (owner rider, `task_petekstatic_topsonly_envelope`).**
  Driven by a committed per-stage audit fixture (`tests/structure_fidelity.rs`:
  dense off-node scatter, exact merges, a data-void margin, local ties):
  - **Solve convergence (S2)**: the flat-seeded kernel stalls metres short of its
    fixed point on sparse-control lattices (slow affine mode under the natural-dip
    boundary). New `srs_gridder::solve_surface_converged` (plane detrend — a plane
    is an exact fixed point, so superposition is exact — + fixed-point restarts)
    is now the stack build's solve entry; defined input nodes remain honoured
    exactly (hard controls).
  - **Explicit extrapolation policy (S3)**: `ExtrapolationPolicy` +
    `with_extrapolation` on builder AND template. Default `DecayToData
    { start_cells: 2, decay_cells: 4 }` — beyond the data hull the solve decays to
    the nearest datum's value; a merged envelope's margin stays merged. Legacy
    unbounded behaviour available as `NaturalDip`. New
    `Surface::taper_beyond_data`.
  - **Cumulative-from-TOP isochores (S6)**: each mapped horizon's isochore is
    anchored at the TOP horizon and made monotone down the stack, so an INTERNAL
    mapped-surface swap leaves every other horizon and the total envelope GRV
    bit-unchanged (was: chained isochores let internal swaps move the base).
  - **Ties by datum substitution (S4)**: well ties now substitute the measured top
    into the horizon's gridded datum and re-run the one stack resolution (builder
    + template) — on a fully-defined lattice a tie moves exactly the tied node
    (radius of influence 0 cells); sparse-lattice reach measured and pinned.
  - **`LayersTruncated` relabel**: the truncation warning is not Follow-specific —
    a Proportional zone over a pinched/merged envelope reports it too; docs now
    name both causes and point at per-zone provenance for attribution.

### Deprecated
- `run_structured_mc`, `run_structured_mc_parallel`, `run_structured_mc_spilled`,
  `run_structured_mc_parallel_spilled` — thin wrappers over `run_mc(tmpl, inputs,
  &McSettings { … })`; bit-identical, removal after the deprecation window.

### Added
- **First-class raw scatter: the engine owns the gridding
  (`task_petekstatic_facade_engagement`).** New `HorizonSource::Scatter(Vec<WorldPoint>)`
  + `WorldPoint { x, y, depth_m }` (world coords) + `StackFrame { ni, nj, georef }`,
  and `StaticModelBuilder::from_scatter_stack` / `StaticModelTemplate::from_scatter_stack`.
  A stack of raw world-coordinate point sets is now gridded **inside** the engine —
  the single scatter-gridding authority — instead of being pre-gridded upstream and
  handed in as fully-defined, all-control `Mapped` surfaces. The engine conditions
  the scatter onto the model lattice with petekTools **bilinear** off-node
  conditioning (`Conditioning::Bilinear`), leaving genuine data voids **`NaN`**, so
  the converged solve + `ExtrapolationPolicy::DecayToData` + isochore build-down act
  on the actual observations: a data-void margin between two exactly-merged horizons
  now **collapses to zero** (audit fixture: 4.0 m phantom → 0.05 m) instead of
  carrying independently-extrapolated fill, and on-data misfit drops to the
  lattice-representation floor (fixture: snap 0.93 m → bilinear 0.17 m rms). Scatter
  is conditioned byte-identically on the deterministic and MC-template paths.
  `HorizonSource::Mapped` is now documented as the pre-gridded escape hatch that
  **bypasses** the engine's solve/conditioning fidelity (loaded grids only). New
  `structure_fidelity` S7 tests (raw-scatter margin collapse + on-data floor,
  world-georef diagonal variant per R1, and the MC-template byte-identity cell per R2).
- **Testing doctrine R2/R4/R5 retrofit (`task_petekstatic_test_matrix`).** The six
  family testing rules are documented in `CLAUDE.md`. A **mode-matrix** acceptance
  file (`srs-model/tests/mode_matrix.rs`) declares the support matrix (in-core ×
  spilled, serial × sharded, single-zone × horizon-stack, wireframe × stack) as a
  header convention and fills its gaps: spilled map/section bundle parity, zoned ×
  sharded MC determinism, and the spilled `zone_stats` **unsupported-cell typed
  error** (the documented template). **`proptest`** is a new `srs-gridder`
  dev-dependency: degenerate-input property tests over layering / collapse /
  order-repair (zero / sub-threshold / inverted / NaN columns, single-cell +
  single-layer grids), each convergence-loop case wrapped in a **hard per-case
  timeout so a livelock FAILS** instead of hanging CI. (See *Fixed* for the bugs
  these surfaced.)
- **`StaticModelTemplate::with_zone_property` / `with_zone_property_mode` — zoned MC
  honours zone-scoped property cubes (`task_petekstatic_canonical_fixes` item 2;
  `question_zoned_mc_zone_pipe_parity`).** The stack template previously realized a
  zone that was staged via a zone-scoped pipe (`with_zone_property`) from its ZONE
  PRIORS only, ignoring the upscale+SGS cube — a zero-spread zoned MC mismatched the
  built model's `in_place_by_zone` on every piped zone (ratio ≈ prior-NTG /
  upscaled-NTG). The template now carries the per-zone pipelines and applies them in
  each realization over the per-zone constant priors, with the draw's level shift on
  top (as the whole-model path does). Acceptance: a zero-spread zoned MC now equals
  `in_place_by_zone` on a zone-piped multizone fixture (bit-level in-repo).
  *Downstream:* peteksim's `stack_template()` must thread its zone pipes into
  `with_zone_property` for the flip to reach the wheel (routed to the coordinator).
- **`StaticModelTemplate::with_well_ties` — MC draws over TIED surfaces
  (`task_petekstatic_template_ties`).** Well ties are draw-invariant, so they are
  applied **once at template construction** (control replacement + re-solve + repair,
  like `from_horizon_stack`'s once-only surface resolution); every draw inherits the
  tied geometry at zero per-draw cost, and the tie residuals ride onto each
  realization's provenance. A tied-template zero-spread realize equals a tied-builder
  build **bit-for-bit** (geometry + cubes + per-zone in-place).
- **Section "sugar-cube" toggle + dip-following trapezoid edges (frozen v4-additive
  schema).** `StaticModelBuilder`/`StaticModelTemplate` gain `with_sugar_cube(bool)`
  (default `false`) → `IntersectionBundle.sugar_cube`. Each `SectionColumn` gains four
  additive per-k arrays — `layer_tops_l`, `layer_tops_r`, `layer_bases_l`,
  `layer_bases_r` — the cell interval bilinearly interpolated from the ZCORN corners
  at the column's **left/right fence edges**, NaN-gapped like `layer_tops`. The
  centroid `layer_tops`/`layer_bases` stay (hover/back-compat). The default now draws
  each cell as a **dip-following trapezoid**; `sugar_cube=true` flattens the edge
  arrays to the centroid trace (one viewer path). SCHEMA_VERSION unchanged (additive).
- **`Gaussian::with_unbounded_search` / `Gaussian::allow_mean_fill`** — explicit
  opt-ins for the two default behaviour changes below.

### Changed
- **Bounded SGS neighbourhood is now the default (`task_petekstatic_canonical_fixes`
  item 5).** The old default searched a **whole-grid** radius (every simulated node a
  candidate → O(N²) scans, >15 min/cube on real-scale lattices). The default is now a
  **bounded** window: `max_neighbours = 16` within `radius = max(1.5·variogram.range,
  4·node_spacing)`. Beyond a covariance range a node's kriging weight is ~0, so the
  bounded field matches the unbounded one within simulation tolerance (max abs Δ
  ≈ 9e-4 on the pipeline fixture, well under 1% of the property range). Restore the
  old window with `Gaussian::with_unbounded_search`.

### Fixed
- **Map and section bundles were broken on spilled (out-of-core) models
  (`task_petekstatic_test_matrix`, R2/R4; the map/section siblings of
  `question_volume_bundle_stack_empty`).** A spilled model keeps a 1×1×1 placeholder
  grid; `map_bundle` read that placeholder for geometry AND cubes, failing with a
  misleading `"1x1 lattice"` error on every spilled model, and `intersection_bundle`
  materialized the geometry via `view_grid` but still read the property from the
  placeholder, erroring `no property 'X'` even though the cube lives in the backing.
  Both now read geometry + cubes from the materialized backing (as the volume bundle
  already did), so all three bundle kinds are supported on spilled models. Regression:
  `mode_matrix.rs` asserts spilled↔in-core parity (within the f32 lane tolerance).
- **Second `collapse_zone` livelock — a sub-threshold layer trapped between
  zero-thickness layers (`task_petekstatic_test_matrix`, R5;
  `question_collapse_zone_livelock`, second instance).** The 46a0345 fix handled a
  whole zone below the threshold; it missed a sub-threshold sliver whose *immediate*
  neighbours are both already-collapsed (zero-thickness) layers while the zone total
  exceeds the threshold — the sliver ping-ponged between adjacent zero slots forever
  (repro `[0.83, 0.0, 0.0, 1.397] @ 0.896`, **found by the new R5 collapse
  proptest**). The merge now targets the nearest **positive** neighbour, looking
  *through* any zero run, so every step strictly reduces the positive-interior-layer
  count → provable termination. Volume-conserving and existing collapse behaviour
  unchanged; held by a bounded, hard-timeout proptest.
- **Non-finite surfaces silently produced a NaN-corner grid (`task_petekstatic_test_matrix`,
  R4/R5).** `layer_grid` / `layer_grid_stack` / `StreamingLayering::prepare` did not
  validate surface finiteness; a NaN/inf depth propagated through `boundary_depth`
  into NaN ZCORN corners — a whole-grid poison no downstream volumetrics could
  distinguish from a real value. The layering seam now rejects a non-finite surface
  loudly (`StaticError::InvalidInput` naming the surface index + node); a
  **partial**-NaN surface (which slips the builder's all-NaN guard) is caught here.
- **Section fence edges collapsed to a flat centroid on world-georeferenced models
  (`task_petekstatic_section_edge_frame`; the third world/local seam bug).** On the
  canonical real-data configuration — a LOCAL-origin cell lattice + a registered
  world `Georef` — both `Polyline` and `AlongBore` sections emitted
  `layer_tops_l == layer_tops_r` **exactly** (a flat "sugar cube" everywhere). The
  trace arrives in **world** coordinates, but a cell's ZCORN corners live in the
  grid's **local** lattice; `fence_edge_depths` clipped the world point + direction
  against the local cell rectangle, so any azimuthed trace missed every cell
  (`smin > smax`) and silently fell back to the centroid. The section march now maps
  the trace point + fence direction into the local lattice frame (new
  `WorldToLattice` affine, applied on **both** the straight-march and the AlongBore
  re-tangent paths) before the clip — one frame, everywhere. A clip miss on a
  non-degenerate direction is now a `debug_assert` + a visible warning (a frame bug,
  never silent); the one honest `l == r` degenerate (a zero-areal-extent / vertical
  bore) is kept explicit and distinct. Local-frame (no-georef) and world-built-corner
  fixtures are unchanged **bit-for-bit** (the map is the identity when the view frame
  already is the local lattice). (Acceptance: a world-georef diagonal fence + bore now
  dip and match the local-frame image; a straight world fence matches direct ZCORN
  corner interpolation bit-level.)
- **AlongBore section fence edges collapsed to the centroid (flat trapezoid)
  on vertical / densely-sampled bores (`task_petekstatic_alongbore_edges`).** A
  bore column's left/right fence edges were interpolated along the raw MD-station
  window direction `(b.x−a.x, b.y−a.y)`. When that window had ~zero areal extent —
  a vertical section (including every well's vertical top before it kicks off) or
  dense stations inside one column — `fence_edge_depths` degenerated to
  `entry == exit == centroid`, so `layer_tops_l == layer_tops_r` exactly and the
  viewer drew the cell flat despite reporting trapezoid mode. The fence direction
  now comes from the **trace's areal tangent** through each column (central
  difference of neighbouring column centres, overall-azimuth fallback), so any bore
  with areal extent — even a near-vertical, slightly-drifting one — dips correctly.
  A straight bore recovers the straight-fence direction **bit-for-bit** (edges equal
  a Polyline along the same line). The one honest degenerate case is a truly vertical
  bore whose trace is a single areal point (one column, no tangent): its edges stay
  the centroid (`l == r`), documented as such. Polyline (straight-fence) edges are
  unchanged. (Acceptance: kickoff-then-deviate surface column now dips; deviated bore
  matches direct ZCORN corner interpolation bit-level; vertical convention pinned.)
- **`collapse_below_m` livelock on sub-threshold zone columns
  (`task_petekstatic_canonical_fixes` item 1; `question_collapse_zone_livelock`).** A
  zone column whose TOTAL thickness was positive but below the collapse threshold
  made `collapse_zone` loop forever — the sliver ping-ponged between interior layers
  because no single layer could ever reach the threshold (repro: interface column
  `[0, 0.3, 0.4]` @ 0.5). Such a **degenerate** zone now snaps onto its thickest
  single layer in one step (volume-conserving; other interior layers go to zero) and
  returns. Regression tests lock the exact repro (with a hard timeout) + a
  multi-layer variant; `StaticModelTemplate::realize`'s per-draw re-layer is covered
  too. Re-enables the 0.5 m collapse on the canonical framework build.
- **Stack volume bundle was empty for out-of-core (spilled) models
  (`task_petekstatic_canonical_fixes` item 3; `question_volume_bundle_stack_empty`).**
  `volume_bundle`/`intersection_bundle` read `grid()`, which is a **1×1×1
  placeholder** on a spilled model (geometry + cubes live in the mmap backing) — so
  every large model (large models are exactly what spills) exported an EMPTY shell
  ("no property" / `N cells - 0 tris` in the viewer). The shell/section exports now
  **materialize the backing** to a whole in-core grid for the (non-hot-path) export;
  a spilled stack shell now matches the in-core shell triangle-for-triangle. Added a
  stack-model shell test (in-core + spilled + realized) asserting a non-empty,
  box-consistent triangle count.
- **Silent constant mean-fill on a data-less simulated layer → loud named error
  (`task_petekstatic_canonical_fixes` items 4 & 6).** A property pipeline whose
  conditioning left a simulated layer with no data silently filled that layer with
  the conditioned mean (structureless, and any collocated trend lost) — and the
  fully-uninformed error did not name the property. Both are now **loud, named
  `InvalidInput`s** at propagate time: the empty-conditioning error names the
  property; a data-less layer errors naming the property (and, via the zone-scoped
  callers, the zone). Opt back into the fill per pipeline with
  `Gaussian::allow_mean_fill`.

### Changed (pre-existing)
- **Out-of-core v2 — the forked volumetric loop is unified behind one abstraction
  (`task_petekstatic_slab_incremental_build`, item 2).** The GRV/HCPV inner loop
  was hand-maintained as **twins** — `srs-volumetrics`' f64 in-core loops and
  srs-model's `spill.rs` f32 streaming mirrors — a silent-parity-divergence risk.
  Both now run through ONE home: a `SlabSource` / `CellSlab` trait + a generic
  `compute_clipped` core (three monomorphic streaming loops over a `Clip { Bulk,
  Single, Two }`) in `srs-volumetrics`. `GridSource` (in-core f64) and srs-model's
  `SpillSource` (spilled f32) both implement the trait; the grv.rs public functions
  and the spilled `bulk_volume` are thin wrappers (signatures unchanged). The
  `StaticModel` in-core-vs-spilled dispatch folds into one place (`clip_of` +
  `volumetrics`), removing the threefold contact-resolution duplication in
  `in_place_impl` / `zone_in_place`. **Zero behaviour change:** golden GRV tests are
  bit-exact, in-core↔spilled tolerance parity holds (≤1e-5), and the MC perf budgets
  are unchanged (monomorphization keeps the codegen monolithic).

### Fixed
- **Zoned-path correctness — three fixes from the composer's planted-truth
  validation (`task_petekstatic_zoned_fixes`).**
  - **Collocated-cokriging trend was a silent no-op on world-coordinate data.**
    The SGS resampled the collocated secondary on the grid's **local** lattice while
    a loaded trend surface is **world**-georeferenced, so the secondary came back
    all-`NaN` and the kernel silently dropped every node to plain SGS — the planted
    correlation was never applied. `Gaussian::with_trend` now resamples a
    world-georeferenced trend at each column's **world** position through the model
    `Georef` (the frame seam, not a kernel change), and a secondary covering **< 50%**
    of the model frame is a **hard `InvalidInput`** (a georef mismatch) instead of a
    silent per-node fallback. A world-georeferenced planted-ρ recovery test asserts
    the effect (recovered field-vs-trend correlation within ±0.20 of the planted 0.6;
    the bug produced ~0.11 at any ρ).
  - **Zoned map outline emitted the unit square.** `from_horizon_stack` framed its
    `map_bundle.outline` as `[0,1]×[0,1]` while the frame + wells are world
    coordinates, collapsing the viewer's content extent. The stack build now emits a
    **world** outline — a georef-derived world-extent rectangle by default, or an
    explicit ring via the new **`StaticModelBuilder::with_boundary`** /
    **`StaticModelTemplate::with_boundary`** (world `[x, y]`). Outline extent now
    matches the frame extent.
  - **Per-zone NTG upscale compressed extremes toward mid-range.** The log-upscale
    binned each sample against the cell interval interpolated at the column
    **centroid** (the 4-corner mean); an off-centroid well on a dipping zone boundary
    then mis-assigned near-boundary samples across the zone, diluting each zone's
    proportion (planted 0.45 read ~0.59 pre-fix). `upscale_cells` now bins against the
    cell interval interpolated **at the well's (x, y)** (bilinear), so a sample lands
    in the zone its depth truly falls in. A well at its column centroid is unchanged
    (bilinear at the centre == the 4-corner mean). Post-fix per-zone recovery is within
    0.05 of the planted target (tightened from the observed ~0.15 error).

### Added
- **Slab-incremental out-of-core build — true O(slab) build peak
  (`task_petekstatic_slab_incremental_build`, item 1).** A forced-spill build no
  longer materializes a whole in-core grid before spilling (the v1 O(grid)
  transient). The gridder path now builds **slab-incrementally**: a new
  `srs-gridder::StreamingLayering` produces ZCORN one k-slab at a time from two
  interface planes (the same f64 `boundary_depth` as in-core, narrowed to f32), and
  the builder writes each ZCORN + constant-cube slab **straight into the store**
  through the new `spill_streaming` writer (`slab_mut_f32` in-place views). The
  build's heap/anonymous working set is **O(slab)** (a slab + two O(area) planes), so
  a model whose in-core live set exceeds RAM can now be **built**. Bit-identical to
  build-then-spill (proved). Measured peak RSS (`dev-docs/bench/scripts/rss_probe.sh`,
  3.2M cells, `/usr/bin/time -l`): in-core **90 B/cell** (whole f64 grid, unbuildable
  above RAM), build-then-spill **135 B/cell**, **streaming 47 B/cell** (≈ the f32
  store alone) — a 2.9× cut, and streaming sits *below* the in-core build. Eligible:
  constant / per-zone-constant priors with collapse off; logs/trend/SGS pipelines or
  collapse fall back to build-then-spill. New surface: `SpillBacking::source` +
  `spill_streaming` + `srs-gridder::StreamingLayering`. **Deferrals (honest):**
  streaming SGS-population + writable-store MC realize (spilled MC realizes in-core
  into the one reused model then flushes — R3, so its peak is the bounded reused
  model); true O(slab) *resident* RSS awaits a page-evicting store writer (raised to
  petekTools). The `SpillBacking::in_place_single` / `in_place_two_contact` methods
  are removed — the model now streams the spilled surface through the unified
  `compute_clipped` over `SpillBacking::source()`.
- **Spilled `in_place_by_zone` — the v1 out-of-core gap closed
  (`task_petekstatic_slab_incremental_build`, item 3).** Per-zone volumetrics on a
  spilled (out-of-core) model now works: because zones are contiguous k-bands
  (`Zone.k_range`), each zone streams its own k-band through the unified core, so
  `in_place_by_zone` dispatches on the backing exactly like the whole-model
  `in_place*` (no more typed "not available spilled" error). Proved: in-core↔spilled
  per-zone GRV/HCPV parity ≤1e-5 (measured worst 1.0e-6) with rollup conservation, on
  the multi-zone regional fixture (single-OWC, two-contact GOC+OWC, and contactless
  zones). (Spilled `zone_stats` — arbitrary-property per-zone cube statistics —
  remains v1-scoped/deferred.)
- **Out-of-core pipeline — memory-budget mode switch + mmap-backed spilled
  `StaticModel` (`task_petekstatic_slab_streaming`, out-of-core rulings R2/R3/R4/R5).**
  The build/template API takes a declared **`MemoryBudget`** (default: a documented
  fraction of physical RAM, `MemoryBudget::default`). Below it → today's in-core
  path, **byte-identical** (proved: same-budget builds are bit-identical, and the
  golden/determinism suite is untouched). Above it → the model **spills** its
  geometry (ZCORN) + property cubes to a petekTools `store` (**f32** slab lanes,
  R4; COORD stays a resident f64 flat lane) and reads them back through **windowed
  mmap views** — the volumetric surface (`in_place`, `in_place_summary`,
  `bulk_volume`, the two-contact split) runs **streaming**, one k-slab at a time,
  with **f64 accumulation** (O(slab) read working set). The switch is **loud**
  (`SpillNotice` names mode + budget + estimate + store path; never an OOM, never
  silent). Spill stores are temp files removed on model drop unless
  `with_spill_persist`. New surface: `MemoryBudget` / `BuildMode` / `decide_mode`
  / `live_set_bytes` / `physical_ram_bytes` / `SpillNotice`, builder
  `with_memory_budget` / `with_spill_dir` / `with_spill_persist`, `StaticModel`
  `is_spilled` / `spill_store_path` / `dims`, `spill_grid` / `spill_grid_to` /
  `SpillBacking`, and spilled MC `run_structured_mc_spilled` /
  `run_structured_mc_parallel_spilled`. **Parity:** in-core↔spilled is
  tolerance-based (R4 honesty clause) — measured on a 2500 m fixture GRV/bulk
  bit-exact, HCPV ≤ ~2e-9, MC oil ≤ ~3.1e-6 relative (contact-boundary sensitive);
  asserted ≤ 1e-5. **MC (R3):** never spills per-draw state — one reusable model,
  realize → summary → discard; the spilled loop reuses **one** per-shard store
  (overwritten per draw, never a new file per draw), and sharded == serial at every
  worker count in the spilled mode. Release `perf_budgets` add a forced-spill build
  + spilled-MC rung. **v1 scope (honest):** spill is *after* the in-core build, so
  build/MC-realize peak RSS stays O(grid) (a fully slab-incremental gridder build +
  writable-store MC realize — true O(slab) *peak* — are the R2 follow-up,
  `question_petekstatic_slab_incremental_build`); spilled `zone_stats` /
  `in_place_by_zone` are deferred (whole-model `in_place*` covers spilled).
- **Stack-aware `StaticModelTemplate` — MC over the multizone framework (P8,
  `task_petekstatic_multizone_2`).** `StaticModelTemplate::from_horizon_stack(stack,
  opts)` runs Monte-Carlo over an ordered horizon stack: it resolves the
  multi-horizon surfaces + per-interface repair **once** (they are draw-invariant),
  then `realize` / `realize_into` vary per draw only the areal footprint (spacing),
  the **per-zone contacts**, and the **per-zone property levels** — new
  `RealizationDraw.zones: Vec<ZoneDraw>` (each `ZoneDraw` owns an optional GOC/OWC
  and optional per-zone φ/NTG/Sw levels; a contactless zone contributes GRV, zero
  HC). Because the geometry topology never varies, the path is **bit-deterministic
  by construction** (two fresh templates + one draw are identical) and works with
  `realize_into` buffer recycling — the **stale-buffer bit-compare** holds (two
  different stacked draws into one reused model == a fresh realize of the second).
  The 2-surface path (and its `sharded_mc_matches_serial_and_is_worker_invariant`
  determinism + release `perf_budgets`) is untouched. Out-of-range zone indices /
  non-physical per-zone contacts are typed `InvalidInput` errors, fail-fast.
- **Per-horizon well ties (P8, `task_petekstatic_multizone_2`).**
  `StaticModelBuilder::with_well_ties(Vec<WellTie>)` on the horizon-stack build:
  each `WellTie { id, x, y, ip, jp, tops: [(horizon, depth)] }` ties its **mapped**
  horizons to the measured tops (the top replaces the map control at the well node
  and the surface is re-solved — the order-repair still guarantees ordering
  afterward), and every tie's **pre-tie residual** (`measured − untied model
  surface`) is recorded in `Provenance.well_ties` (`WellTieRecord` /
  `HorizonTieResidual`). A **tops-only** horizon's tie is recorded against its
  pick-conditioned drape (already conditioned by the picks). The map bundle now
  populates `wells` from these ties: `WellMarker` gains `ties: Vec<WellTieResidual>`
  (per-horizon) and its `tie_residual_m` becomes the mean residual — the viewer's
  `wells[].ties` (SCHEMA_VERSION 4). Unknown horizon/zone names and off-lattice tie
  nodes are build-time `InvalidInput`/`Grid` errors.
- **Per-zone property population (P8, `task_petekstatic_multizone_2`).** A
  horizon-stack build can now give each zone its own distribution/variogram/logs:
  `StaticModelBuilder::with_zone_priors(zone, ConstantPriors)` overwrites
  `PORO`/`NTG`/`SW` across that zone's `k`-range (the per-zone level a sand vs shale
  zone owns), and `with_zone_property(zone, PropertyPipeline)` runs a full
  upscale→SGS pipeline **restricted to that zone's `k`-range** via the new
  `PropertyPipeline::apply_in_zone(grid, k_range)` — each zone gets its own
  variogram / trend / log-conditioning, merged into the cube so other zones' slices
  are untouched. Applied after the base population and any whole-model
  `with_property` pipelines; per-zone reports land on `Provenance::property_reports`.
  The SGS per-layer seed is per-global-`k` (range-independent), so a zone's field is
  bit-identical whether populated alone or as part of the whole grid, and
  `apply_in_zone(0..nk)` reproduces `apply()` exactly. The petekTools geostat kernels
  (SGS / collocated cokriging) remain the engines — the pipeline only scopes them.
  An unknown zone name is a build-time `InvalidInput`. Geometry (hence volume) is
  never touched, so the multizone conservation invariants hold (new coverage).
- **Section-bundle interior-horizon traces + `SCHEMA_VERSION 3 → 4`
  (`task_petekstatic_multizone_2`).** `IntersectionBundle` now carries
  `horizon_traces: Vec<HorizonTrace>` — one depth polyline per *interior* framework
  horizon (every zone-bounding horizon strictly between the structural top and base;
  `N − 2` for an `N`-horizon stack, top→down). Each `HorizonTrace { name, depths }`
  runs **parallel to `columns`**: `depths[c]` is the horizon's depth at `columns[c]`,
  taken as the zone-top interface (the top-depth of the zone's first cell). The
  structural top/base are *not* repeated (they remain the first/last active
  `layer_tops`/`layer_bases`); a single-zone (2-surface) model emits an empty
  `horizon_traces`, so the block is **additive / backward-compatible**. The viewer
  (petekTools, separate repo) is asked to render it — graph
  `question_viewer_interior_traces`.

### Changed
- **`VolumeBundle` is now the exterior shell + a binary-block payload
  (SCHEMA_VERSION 1 → 3; `task_suite_bundle_binary`, the P9 payload-killer).**
  `volume_bundle` stops shipping the full 8-vertex / 12-triangle cell soup and
  emits only the **visible faces** (a cell face is emitted iff the neighbour across
  it is inactive/absent or out of bounds); shared vertices are deduplicated and the
  per-cell arrays are compacted to the shell cells, with a per-triangle `tri_cell`
  index recovering cell identity. The big arrays serialize as raw little-endian
  binary blocks — base64-wrapped in a JSON envelope
  (`VolumeBundle::write_self_contained`, one file) or split into a `model.bin`
  sidecar with an `(offset,length)` manifest (`VolumeBundle::write_sidecar`). Both
  **stream** to an `io::Write` with no `serde_json::Value` tree. **Measured @1M
  cells (200×200×25): ~557 B/cell → ~6.65 B/cell self-contained / ~5.0 B/cell
  sidecar (≈84–111× smaller, 613 MB → ~6.7 MB — an order of magnitude under V8's
  ~512 MiB inline-script wall); triangles 12M → 200k (1.67%).** Struct changes:
  `positions`/`indices` now hold only shell geometry; new `tri_cell: Vec<u32>`;
  `cell_values` `Vec<f64>`→`Vec<f32>` and `zone_ids` `Vec<u32>`→`Vec<u16>`, both
  **compacted** to shell cells; `vertex_values` and the per-cell `active` mask
  removed (colour is per-triangle via `tri_cell`; absent cells never reach the
  shell). The authoritative binary spec is in `API.md`; the petekTools viewer
  decodes exactly it.
  - **Threshold semantics (decision (b)):** the client slider is a documented
    **shell-only** filter; true interior exposure at a cutoff is a server
    regeneration — new `StaticModel::volume_bundle_thresholded(property, cutoff,
    keep_above)` re-cuts the shell treating sub-cutoff cells as absent.
  - **Map/section bundles** gain `write_json` (the same streaming writer); their
    JSON structure is unchanged but they carry the family SCHEMA_VERSION 3.

### Added
- **`StaticModelTemplate::realize_into(&mut self, &RealizationDraw, &mut StaticModel)`
  — the ratified buffer-recycling realize variant (`task_petekstatic_realize_into`,
  `decision_staticmodel_regen_seam` amendment 3).** Realizes into a **reused**
  `StaticModel`, overwriting its geometry (ZCORN + COORD) and property-cube
  allocations in place instead of allocating a fresh ~100 MB/draw at 1M cells
  (ZCORN 64 MB + cubes 24 MB + scratch). `realize` is now
  `new-empty-then-realize_into` (one code path). New public surface:
  `StaticModelTemplate::realize_into` + `reusable_model` (a per-worker model to
  drive a realize_into loop against); in `srs-gridder`, `layer_grid_stack_into` +
  `LayerScratch` + `StackLayering` (the geometry-recycling core; `layer_grid_stack`
  is now its allocating wrapper). `run_structured_mc`/`_parallel` drive
  `realize_into` with one reused model per shard. **Determinism is unchanged** —
  the sharded MC is bit-identical to the serial run at every worker count; a
  stale-buffer test proves two consecutive draws (different gross/porosity/shift)
  into a reused model are bit-for-bit identical to a fresh `realize`.
  - **Measured @1M cells (200×200×25, LevelShift MC):** the wall is **memory
    bandwidth** (the per-draw ZCORN/cube *writes*), not allocation — the system
    allocator already recycles the large blocks — so recycling holds serial
    ≈flat (within ~5%) and lifts parallel scaling modestly (~2.05× → ~2.17× at 8
    workers) rather than cutting serial time. Release-gated perf budgets
    re-baselined accordingly (new `warm_realize_into` budget; `run_structured_mc`
    budget 25 s → 26 s).
- **Multi-zone regional framework — `StaticModelBuilder::from_horizon_stack`
  (`task_petekstatic_multizone`).** The canonical framework is now an ordered
  horizon stack: `N` horizons top→down → `N − 1` named zones, each with its own
  conformity/layering **and** its own fluid contacts. New public surface:
  `HorizonStack`, `StackHorizon`, `HorizonSource::{Mapped, TopsOnly}`, `Pick`,
  `StackZone`; `ZoneTable::from_stack`; `srs_gridder::{layer_grid_stack,
  StackedLayeredGrid, StackedZone, ZoneLayerSpec}`.
  - **Three horizon source kinds.** `Mapped` surfaces (gridded/points, tied or
    untied) and **tops-only** internal horizons (no mapped surface — defined by
    well **picks**, draped conformally from the nearest mapped horizon above at
    pick-controlled thickness; constant-thickness fallback for a single pick).
  - **Per-interface order-repair.** `with_min_thickness_m` generalizes per
    consecutive-horizon interface; repairs recorded in
    `StackProvenance::interface_repairs`.
  - **Per-zone layering + conformity.** Each zone gets its own `Conformity` + `dz`;
    total `nk` = sum of per-zone counts; the `MAX_NK` (200) cap applies to the
    **total** with the per-zone breakdown retained.
  - **Per-zone contacts + volumetrics.** `StaticModel::in_place_by_zone() ->
    ZonedInPlace` (`ZoneInPlace` per zone + rollup `total`): each zone clips its
    hydrocarbons against **its** contacts (two-contact gas/oil split per zone); a
    **contactless zone contributes gross bulk but ZERO hydrocarbon in-place** (no
    contact = no known accumulation). `sum(zone volumes) == total` (conservation).
    New `srs_volumetrics::{compute_in_place_zone, compute_in_place_two_contact_zone,
    compute_zone_bulk}`.
  - **Per-zone stats.** `StaticModel::zone_stats(property) -> Vec<ZoneStat>`
    (count/mean/min/max per zone over full-grid cubes, active cells only).
  - **Provenance records the stack** (`Provenance::stack: Option<StackProvenance>`
    with `ZoneProvenance` + `InterfaceRepair`).
  - Backward-compatible: `from_wireframe` (Top+Base) is the single-zone degenerate
    case. Follow-ups (noted): per-zone priors / SGS variograms, per-zone MC contact
    draws, section-bundle interior horizon traces.
- **Cell-collapse threshold — `with_collapse_below_m(f64)`** (builder + template,
  default OFF). Sub-threshold cells collapse **volume-conservingly**: the sliver's
  thickness merges into a thicker **zone-interior** neighbour (never across a zone
  boundary, never deleting rock). Per-zone collapsed counts →
  `BuildWarning::CellsCollapsed`; collapsed cells are `NaN`-marked in the section
  bundle by the existing `dz ≤ ε` rule.
- **Layering conformity styles — `Conformity::{Proportional, FollowTop { dz_m },
  FollowBase { dz_m }}` (`task_petekstatic_layer_conformity`).** Layering can now
  optionally follow zone edges (horizon tops/bases) instead of only proportional
  subdivision. `FollowTop` drapes each layer surface parallel to the **top** at a
  constant `dz_m` (layer `k` base = `top + (k+1)·dz`), truncating the deep layers
  against the base where the interval pinches; `FollowBase` mirrors it (onlap
  against the top). Under a Follow style `nk` is **dz-derived** — `ceil(max column
  thickness / dz_m)`, capped at `srs_gridder::MAX_NK` (200) — so the thickest
  column is fully layered and thinner columns truncate. `Proportional` (the
  default) is unchanged and still honours `nk`.
  - **Truncation = zero-thickness collapse** (no active-mask array): a truncated
    cell collapses onto the pinch-out horizon so its bulk volume is 0. Volumetrics
    already excludes it (GRV/HCPV are conformity-invariant — total column volume is
    conserved to FP), property population writes only harmless finite values into
    it (never `NaN` into an active cell), and the section bundle marks it `NaN`.
  - **`layer_grid` now returns `LayeredGrid`** `{ grid, nk, truncated_cells,
    nk_capped }` (was `Grid`) — the effective `nk`, truncated-cell count, and cap
    flag the builder/template stamp into provenance.
  - **`Provenance.nk` is now the EFFECTIVE (dz-derived) nk**; new
    `BuildWarning::{LayersTruncated { cells }, LayerCountCapped { nk }}`
    (informational). `IntersectionBundle` `SectionColumn` arrays stay `nk`-sized;
    an inactive/truncated layer is `NaN` in `layer_tops`/`layer_bases`/`values`.
  - **Replaces** the placeholder `Conformity::{Onlap, Truncation}` (fixed-`nk`,
    non-truncating) — a breaking rename to `FollowBase`/`FollowTop` with real
    dz-derived truncation.
- **View bundles — the viewer export seam (`view` module, SPEC §7e; graph
  `decision_viewer_home_product` / `task_petekstatic_view_bundles`).** The model
  now exports typed, JSON-stable inspection bundles for a *separate* viewer
  codebase: petekStatic exports, the viewer renders — it never computes. All are
  serde-`Serialize`/`Deserialize` value types on one shared, versioned schema
  (`SCHEMA_VERSION`), SI units / positive-down depth, `NaN` = undefined.
  - **`MapBundle` = `model.map_bundle(&MapSpec)`** — areal (plan-view): realized
    structural depth surfaces (top + base, grid-georeferenced), property maps as a
    single **k-slice** and the **zone/interval average** (the useful default), the
    outline ring(s) in world coords, well surface markers (with tie residuals
    where provenance carries them), and per-contact **subcrop masks**. Every layer
    carries its value range for a legend. All areal layers share one regular
    `GridFrame` (`origin`/`spacing`/`ncol=ni`/`nrow=nj`), reconstructed from the
    grid's column-centroid `xy↔ij` lattice; row-major `values[j*ncol + i]`.
  - **`IntersectionBundle` = `model.intersection_bundle(&SectionSpec, property)`**
    — a vertical cross-section walked through the areal lattice along a world
    **polyline** (`SectionSpec::Polyline`) or a **bore trajectory**
    (`SectionSpec::AlongBore`). Emits ordered `SectionColumn`s: distance-along,
    per-layer property + cell top/base depths (the horizon traces are the first
    top / last base), and — for a bore — the path's own z overlay. Section-wide
    fluid contacts carried alongside. Raw metres (the viewer scales for vertical
    exaggeration).
  - **`VolumeBundle` = `model.volume_bundle(property)`** — the corner-point cell
    mesh (8 vertices / 12 triangles per cell) + per-cell property values, zone
    ids, and an active mask (threshold filtering stays viewer-side; the value
    arrays it reads are included). **The mesh builder now lives here** — moved
    DOWN into the GEOMODEL layer from petekSim's `srs-core/mesh.rs` per the layer
    charter (the DAG flows downward; peteksim retires its copy onto this bundle at
    the viewer wave).

- **`min_thickness_m` post-gridding order-repair for crossed/thin columns (R-c).**
  Real thin margins: independent gridding of Top and Base can overshoot and
  re-introduce a crossing (base above top) a pointwise pre-repair had removed,
  leaving the caller stuck behind the (correct, loud) `CrossedSurfaces` guard.
  `StaticModelBuilder::with_min_thickness_m(f64)` and
  `StaticModelTemplate::with_min_thickness_m(f64)` opt into a repair: where the
  gridded base sits less than `min_thickness_m` below the top, the base is pulled
  **down** to exactly `top + min_thickness_m`, **preserving the top** (the
  better-constrained seismic pick). Off by default (the `CrossedSurfaces` guard
  stays the default); when enabled, `Provenance.warnings` records a new
  `BuildWarning::ThinColumnsRepaired { columns, worst_m }` (repaired-node count +
  worst, most-negative original separation). Takes precedence over
  `with_clamp_base_to_top` when both are set.

### Added
- **Rayon-sharded structured MC — `run_structured_mc_parallel(tmpl, inputs, n,
  seed, workers)` + `default_mc_workers()` (`task_petekstatic_engine_perf`).**
  Splits the `n` draws into `workers` contiguous shards, each realizing its slice
  on its own template clone (`StaticModelTemplate: Send`), and recombines the
  per-draw vectors in draw-index order. **Determinism contract:** same
  `(inputs, n, seed, workers)` → bit-identical vectors run to run (single seeded
  input stream + deterministic shard split + index-order recombination); across
  *different* worker counts the draw multiset is identical — bit-identical in the
  common no-per-draw-structural-shift case, else per-boundary draws may differ only
  within the warm-start solver tolerance. The LevelShift pattern cache is pre-warmed
  once so worker clones don't each re-propagate. `default_mc_workers()` =
  `min(6, cores)` (the MC loop is bandwidth-bound; scaling saturates below core
  count). At 1M cells / 1000 LevelShift draws: 16.2 s serial → 7.9 s at 8 workers
  (2.0×); combined with the in-place win the full MC wall is 28.9 s → 7.9 s (3.7×).

### Performance
- **In-place volumetrics ~2.9× faster on the Monte-Carlo hot path
  (`task_petekstatic_engine_perf`).** `StaticModel::in_place` / `in_place_summary`
  no longer pay a per-cell pillar-`xy` interpolation *division* on a vertical grid
  (every grid today): `CornerPointGeom` flags a vertical lattice at construction
  and builds each cell's corners straight from the pillar tops + ZCORN — **bit-for-bit
  identical** to the general path (a vertical pillar's `xy_at` returns `top.xy`). At
  1M cells `in_place_summary` drops 18.9 ms → 6.6 ms. New cheap accessors
  `Grid::cell_centroid_z_at` / `cell_volume_at` (by linear index); the volumetrics
  cores are a single fused pass that rejects below-contact cells on the cheap
  centroid-z test before touching geometry. Numerics unchanged (all GRV/HCPV/OOIP
  goldens + multi-zone conservation hold to 1e-9).
- **`#[inline]` on the `srs-grid` geometry primitives** (`Point3` algebra,
  `hexahedron_volume`, `Cell`/`CornerPointGeom` accessors) so they inline across the
  crate boundary without LTO — `realize` (populate + property traversal) ~23% faster
  at 1M (13.6 ms → 10.5 ms).

### Fixed
- **View bundles now emit ONE consistent frame for world (UTM) models — the map
  raster, outline, wells and world fence/bore sections overlay (F5-class export
  fix).** For a world-frame model the grid is an area-scaled square re-origined to
  a *local* box (e.g. centroid lattice `15..2985`), but `MapBundle.outline` (and,
  once retained, wells) pass through in *world* (UTM) coordinates — so the raster
  (`GridFrame`) and the vector layers did not overlay, and a world-coordinate
  fence/bore trace marched outside the local lattice and yielded **zero** section
  columns. The `StaticModel` now carries a registered world georeference
  (`Georef` — column `(0,0)`'s world centroid + world column spacing), set via
  `StaticModelBuilder::with_georef` / `StaticModelTemplate::with_georef` (the
  upstream seam supplies it from the source horizon lattice / CRS, the same
  `xy↔ij` the well registration uses). `GridFrame` is then emitted in that world
  frame — `map_bundle` overlays the world outline, and `intersection_bundle`
  traces world polylines / UTM bore trajectories through the same mapping. With no
  registered georeference the frame **degenerates to the grid's local
  column-centroid lattice** exactly as before (synthetic square / box —
  bit-identical, no schema change). This is a **values** fix: `SCHEMA_VERSION`
  stays `1`; the grid geometry / GRV is untouched (the georef only *labels* the
  local column lattice with its world frame).
- **`McMode::LevelShift` now works on log-conditioned real models — fraction
  cubes saturate rather than escaping `[0,1]` (F9).** A conditioned PORO/NTG/SW
  cube legitimately holds boundary cells (a non-net `NTG=0`, an aquifer `SW=1`).
  Level-shift added the drawn shift to **every** cell, pushing those boundary
  cells outside `[0,1]`, and the per-cell H2 range check then (correctly) rejected
  the draw — so property uncertainty on **any** log-conditioned real model failed
  at draw #0. The per-cell shift application now **shift-then-clamps** the three
  fraction cubes (`PORO`/`NTG`/`SW`) to `[0,1]`: boundary cells **saturate** (an
  `SW=1` aquifer cell stays 1 under a positive shift — physically right — and
  slides into the interior under a negative one), interior cells move by the full
  drawn amount, and `NaN` (undefined) cells pass through untouched. H2 still
  validates the **drawn inputs** themselves (a garbage sampler is still an error);
  only the per-cell cube application saturates.
- **Horizon z-datum unified inside `Wireframe` — latent multi-horizon role
  hazard (review Z1).** petekIO delivers surfaces as **negative-down subsea
  elevation**, but the GEOMODEL layer works in **positive-down `depth_m`** and
  `assemble_wireframe` was passing surface values through **unflipped**. Role
  assignment picks `Top` by numeric minimum (shallowest depth) — under raw
  negative-down elevation that minimum is the *most negative* = the
  **structurally deepest** surface, so on multi-horizon input the **deepest
  surface was labelled `Top`** (and the shallowest `Base`), and negative-down
  horizons were mixed with positive-down contacts in one `Wireframe`. `srs-data`
  now **negates surface z at the ingest boundary** (`surface_depths`), landing
  horizons on the same positive-down `depth_m` datum as the contacts — one datum
  inside the `Wireframe`, correct `Top`/`Base` ordering. Proven by a failing-first
  multi-horizon test. The old coordinate-flip deferral is retired (its blocker,
  petekIO's imperial `SummaryInputs` contacts, resolved when petekIO went SI).
  Well-curve positions (`logs::petro_samples`) carry the same negative-down
  elevation; the flip onto model depth there is the downstream binning seam's job
  (documented). Single-horizon models (`Top`-only) are unaffected.

### Changed
- **srs-data curve matching routes through petekio's `canonical_mnemonic`
  (review D2, dedup + correctness).** `logs::petro_samples` no longer hardcodes a
  local mnemonic set that folded total porosity (`PHIT`) into porosity and missed
  vendor variants. It now matches **effective** porosity (`PHIE`) and **effective**
  water saturation (`SW`) via petekio's one canonical alias table, so `EFFPHI`/
  `PHIEF`/`PHI` → `PHIE` and `SWE`/`SUWI`/`SW_E` → `SW` resolve, while `PHIT`/`SWT`
  (total) stay distinct and are excluded from the effective (φ, Sw) sample.
- **`srs-uncertainty` delegates its P90/P50/P10 digest to petekTools (review
  D1, dedup).** `PercentileSummary::from_realizations` now calls
  `petektools::sampling::reservoir_summary` instead of a hand-rolled type-7
  percentile kernel; the typed empty-input error (H1, FFI-visible) is kept as a
  thin wrapper. The duplicated `percentile_sorted` public fn is removed (pre-0.1;
  no external callers) — the one percentile home is now petekTools.
- **`Wireframe.horizons` is now `Arc<Vec<Horizon>>` (review P2, perf).** The
  horizon set is realization-invariant, so the MC template now **shares** it
  across every `StaticModelTemplate::realize` (an O(1) refcount bump) instead of
  deep-copying the per-node depth cubes each draw; only the two contacts are
  rebuilt per realization. The no-shift realize branch also **borrows** the
  template's control set (`Cow`) instead of cloning it. Together: ~2% faster per
  realize at 50×50×20 (1.19 → 1.17 ms), scaling with lattice size × horizon
  count. Pre-0.1 source break: construct with `Arc::new(vec![...])` and mutate
  with `Arc::make_mut(&mut wf.horizons)`; all read paths are unchanged (`Deref`).

### Added
- **Structured Monte-Carlo driver + tornado (`task_peteksim_mc_structured`,
  `task_peteksim_tornado`).** The uncertainty layer over the regeneration
  template.
  - `run_structured_mc(&mut template, &McInputs, n, seed) -> McResult`: a seeded,
    bit-reproducible run. `McInputs` carries a petekTools `Sampler` per uncertain
    input (`Input::plain` / `Input::clamped` — the latter wraps petekTools'
    `Clamped`, not a bespoke wrapper); `boi` / `bgi` are sampled but applied at
    the volumetrics surface, never on the draw (ratified). `McResult` KEEPS the
    per-draw `oil_sm3` / `gas_sm3` / `grv_m3` vectors (W17), summarises them via
    petekTools `reservoir_summary` (`summary()` = oil P90/P50/P10), and retains
    the realized inputs for tornado reuse. 10k draws (11×11×5) in ~190 ms
    single-threaded (~19 µs/draw).
  - **Error policy = fail-fast with the draw index:** the first draw whose
    realization/volumetrics fails surfaces as the new typed
    `StaticError::McDraw { index, source }` (the H2 typed error carried through
    the loop). Clamp fraction samplers and no valid configuration trips it.
  - `aggregate_field(&[&McResult], Correlation)` delegates to petekTools
    `aggregate` — the Independent (narrow) / Comonotonic (wide downside) bracket.
  - `tornado(&mut template, &McInputs, n, seed, lo_pct, hi_pct) -> Vec<TornadoBar>`:
    one-at-a-time swings of oil in-place Sm³, each input pivoted at its **realized**
    percentiles (others at P50), pre-sorted by swing descending. Same `(n, seed)`
    as the MC run reuses its exact draws.
  - `StaticError` now composes petekTools' `AlgoError` via `#[from]`
    (TOOLKIT→GEOMODEL seam) so sampling/stats `?`-chain; `srs-model` re-exports
    `Sampler` / `Correlation` / `ReservoirSummary`.
- **Per-property geostatistical population pipeline (P5,
  `task_petekstatic_property_modelling`).** Properties are modelled one at a time
  through a visible pipeline: `PropertyPipeline::new("PHIE").upscale(wells,
  UpscaleMethod::Arithmetic).propagate(Gaussian::new(variogram, seed))`.
  - **upscale** (`upscale_cells`) is a first-class inspectable step: positioned
    `WellLog`s snap to their areal column and their in-cell samples are upscaled
    (`srs-petro` power means) into per-cell conditioned values (`NaN` where no log
    passes), returning an `UpscaleQc` (upscaled-vs-log stats).
  - **propagate** runs petekTools' sequential Gaussian simulation (`geostat::sgs`)
    **per k-layer**, conditioned exactly on the upscaled cells, filling every cell
    with a seeded, reproducible draw that honours the wells and reproduces the data
    histogram. `Gaussian::with_trend(TrendSurface, corr)` steers it via collocated
    (Markov-1) cokriging, the trend resampled to the model lattice with the shared
    `petektools::resample` kernel. This supersedes the interim `TrendSurface`
    multiplier hook.
  - Wired into `StaticModelBuilder::with_property(PropertyPipeline)`; the
    `PropertyReport` (upscale QC + `propagated`) lands on
    `Provenance.property_reports`. Regular / axis-aligned column lattices (the
    `layer_grid` box/conformable grids); rotated pillars are future work.
- **Monte-Carlo property modes (`decision_mc_composition`).** `StaticModelTemplate`
  gains `with_property(pipeline)` (default `McMode::LevelShift`) and
  `with_property_mode(pipeline, mode)`:
  - **`LevelShift`** (default) propagates the field **once** (cached at the first
    realization) and reuses that pattern each draw, adding only the draw's
    per-property additive shift (`RealizationDraw::with_property_shift`) — same
    spatial pattern, moved level. Bench: a `realize` stays **~320 µs** (50×50×20).
  - **`Resimulate`** re-runs SGS with a fresh per-draw reseed (derived from
    `seed_index`) — a new pattern each draw, at its true cost (bench **~139 ms** at
    50×50×20). Both modes are bit-reproducible for identical draws.
  Additive on the ratified seam: new `RealizationDraw.property_shifts` +
  `with_property_shift` / `property_shift` (the `#[non_exhaustive]` + `with_*`
  pattern); no re-ratification.

### Changed
- **SI/metric unit standard adopted (2026-07-04, `decision_si_units_standard`).**
  petekStatic is now metric family-wide: coordinates/depths/lengths in **metres**,
  areas **m²**, volumes **m³** internally with **mcm** (GRV) / **MSm³** (oil) /
  **bcm** (gas) reporting. Renames (pre-0.1.0, breaks no downstream): `BoxSpec`
  `area_acres`/`gross_height_ft`/`top_depth_ft` → `area_m2`/`gross_height_m`/
  `top_depth_m`; `BuildOpts`/`RealizationDraw` the same; `GriddedDepth.depth_ft` →
  `depth_m`, `Contact.depth_ft` → `depth_m`, `Segment.throw_ft` → `throw_m`;
  `StaticError::CrossedSurfaces{worst_ft}` → `worst_m`; `Surface::offset_by_field`
  `dz_ft`→`dz_m`; `compute_in_place*` contact params `*_ft`→`*_m`. `InPlace`/
  `ZoneVolumes` volumes `grv_ft3`/`hcpv_ft3` → `grv_m3`/`hcpv_m3`; accessors
  `grv_acre_ft`→`grv_mcm`, `ooip_stb`→`ooip_sm3` (+ new `oil_msm3`),
  `ogip_scf`→`ogip_sm3` (+ new `gas_bcm`), `gas_zone_ogip_scf`/`oil_zone_ooip_stb`
  → `*_sm3`. **FVF `OilFvf`/`GasFvf` relabelled rb/STB & rcf/scf → Rm³/Sm³**
  (dimensionless, numerically identical — a relabel, not a conversion). Golden
  values re-expressed in SI; a `srs-volumetrics` parity test proves the SI OOIP
  converts back to the identical imperial STB to FP tolerance. Reporting factors
  come from `petektools::units`.
  - **srs-data imperial seam shim RETIRED (2026-07-04).** petekIO HEAD landed the
    SI `SummaryInputs` flip (`area_m2`, `net_pay_m`, positive-down `owc_depth_m`/
    `goc_depth_m`), so `srs-data` reads the metric fields straight through — the
    acres→m² / ft→m conversions and the `ModelScalars` `*_acres`/`*_ft` field names
    are gone (now `area_m2`/`net_pay_m`/`owc_depth_m`/`goc_depth_m`).
  - **Still deferred:** the internal negative-down elevation *coordinate* flip
    (horizon `Surface` node values are subsea elevation while our `depth_m` is
    positive-down) — contacts are now unambiguous positive-down metres, but the
    geometry-z flip stays future work.
- **Log-population presort (2026-07-04, V2).** `with_logs` (builder + template)
  now TVD-sorts the samples once (`sort_by_tvd`); `populate_from_logs` binds each
  cell's in-range window with two binary searches over the sorted array instead of
  a full linear scan per cell. Bench (40k cells × 1500 samples): **15.7 ms → 2.9 ms**
  min (~5.3×; the gap widens with denser logs). For the template the sort is paid
  once and reused across every realization.
- **Builder ↔ template kernel unification (2026-07-04, R2,
  `decision_gridder_kernel_unification`).** `StaticModelBuilder::build` now solves
  the structural surfaces through the **same petekTools warm kernel** the template
  uses (shared `warm_surface` path), not the cold `solve_surface`. The two builds
  previously diverged (~20× GRV on thin columns) because the cold reference kernel
  and the warm kernel interpolate unpinned interior nodes differently — invisible
  on fully-pinned lattices, ruinous where the thin gross is a small difference of
  two large surfaces. `template.realize(draw with gross == mean(g))` now
  reproduces the deterministic build within tight tolerance on **all** column
  shapes (thin-sparse case added alongside the existing wedge). The cold
  `srs-gridder::solve_surface` stays as the accuracy reference, off the build path.

### Added
- **`TrendSurface` world georeference (2026-07-04, R4).** `with_georef(origin_x,
  origin_y, node_dx, node_dy)` maps the trend field to world coordinates;
  population then resamples each model column by its **world** `(x, y)` (nearest
  node, from the grid's cell centroids) instead of by index fraction, so the trend
  lands where its data is rather than being aligned to the model lattice by luck.
  Non-positive spacing is ignored (stays index-space); `is_georeferenced()` reports
  the mode. Bilinear resampling still deferred to P5.
- **Summary-only in-place (2026-07-04, V7).** `StaticModel::in_place_summary` +
  `srs-volumetrics::compute_in_place_summary` /
  `compute_in_place_two_contact_summary` return the same GRV/HCPV/OGIP/OOIP
  aggregates as `in_place` but leave `InPlace::per_cell_hcpv` empty — skipping a
  `cell_count`-length allocation per realization on the MC hot path where only
  aggregates feed the P-curve. (The buffer-recycling `realize_into` half of V7 is
  still deferred — it needs `layer_grid_into` + `populate_into` in
  srs-gridder/srs-grid.)
- **`RealizationDraw` named setters (2026-07-04, V8/R6).** Fluent
  `with_area` / `with_gross` / `with_contact` / `with_porosity` / `with_ntg` /
  `with_sw` overrides (alongside the existing `with_goc` / `with_sw_gas` /
  `with_structural`) so a sampler can build a draw field-by-field from a base
  instead of positionally; `::new` stays for the load-bearing seven-scalar path.
- **Gas-cap connate-water override — `sw_gas` (2026-07-04, R3).** A single shared
  `SW` cube over-states gas-cap OGIP when the gas leg's connate water is lower than
  the oil leg's. `compute_in_place_two_contact` gains an optional `sw_gas: Option<f64>`
  applied to gas-zone cells only (the oil leg keeps the cube);
  `RealizationDraw::with_sw_gas` (additive, `#[non_exhaustive]`) and
  `StaticModelBuilder::with_sw_gas` thread it through `StaticModel::in_place`. This
  is one scalar, not a saturation-height/PVT model — transition-zone and PVT
  saturation stay in petekSim.
- **Real boundary footprints (2026-07-04, `task_petekstatic_boundary_rings`).**
  `srs-data`'s `assemble_wireframe` now reads a supplied polygon's **true exterior
  ring** via petekio `PolygonSet::rings()` (≥0.2.2) into `Boundary::ring` instead
  of rebuilding a bbox rectangle — non-rectangular field outlines are preserved
  (shoelace area matches the polygon, not its bounding box). The bbox is kept only
  as a documented fallback for a ringless (degenerate) set; the stale "no
  ring-vertex accessor" comments are corrected. Resolves `q_petekio_polygon_rings`.
- **Base-above-top guard (2026-07-04, R1).** A `Base` horizon that crosses above
  the `Top` per column (negative gross) no longer silently collapses GRV: the
  build and every `realize` now validate per-node separation and raise a typed
  `StaticError::CrossedSurfaces { nodes, worst_ft }` (offending node count + worst
  crossing). Opt into `StaticModelBuilder::with_clamp_base_to_top(true)` /
  `StaticModelTemplate::with_clamp_base_to_top(true)` to clamp the offending
  columns to zero gross instead (leaving the rest untouched). New
  `srs-gridder::Surface::guard_below(&top, clamp)`.
- **Base-horizon relief wired through the build (2026-07-03, S1 fix).**
  `StaticModelBuilder::from_wireframe` now solves a supplied `Base` horizon's
  real relief into the pillar bases (spatially varying gross thickness);
  `BuildOpts.gross_height_ft` is the fallback when no Base horizon is supplied
  (backward-compatible). `Provenance` gains `warnings: Vec<BuildWarning>`;
  supplied horizons the build cannot consume (Intermediate; a second or
  lattice-mismatched Base) raise a non-blocking `BuildWarning::UnusedHorizon`
  instead of being silently discarded. Locked by an analytic wedge test.
- **Areal trend hook — external-drift-lite (2026-07-03, INTERIM).**
  `TrendSurface` + `with_areal_trend` on both `StaticModelBuilder` and
  `StaticModelTemplate`: a gridded areal multiplier field, nearest-node
  resampled to the model column lattice, mean-normalized (field-mean
  preserved), applied per-column to NTG (and φ via `TrendSurface::
  with_porosity`) after population. Trend gives lateral *shape*, the
  prior/draw gives the *level*; `NaN` trend nodes fall back to 1.0. Interim
  until collocated cokriging lands (P5 `task_petekstatic_property_modelling`).
- **Template gross scaling over real base relief (2026-07-03,
  `decision_template_gross_scaling`).** `StaticModelTemplate` built from a
  wireframe with a `Base` horizon no longer flattens the base to
  `top.offset_by(gross_height_ft)`: it extracts the nominal per-node gross
  field `g(x,y)` (base solved once in kernel space) and each draw's
  `gross_height_ft` scales its level — per-node gross = `g × gross / mean(g)`,
  so a draw at `mean(g)` reproduces the deterministic build bit-for-bit
  (wedge-locked). No Base = the constant offset exactly (backward-compatible).
  New `srs-gridder::Surface::offset_by_field(&[f64])` (per-node offset).
- **Two-contact volumetrics — gas cap + oil rim (2026-07-03).**
  `srs-volumetrics::compute_in_place_two_contact(grid, goc_ft, owc_ft)`
  partitions cells gas cap / oil leg / water by centroid depth; `InPlace`
  gains `gas`/`oil: Option<ZoneVolumes>` + `gas_zone_ogip_scf` /
  `oil_zone_ooip_stb`. `StaticModel::in_place()` auto-splits when the
  framework carries a GOC plus a lower contact (lone contact stays generic);
  `RealizationDraw::with_goc(depth)` (additive, `#[non_exhaustive]`) realizes
  two-contact models with GOC-above-OWC validation. Geometry + in-place split
  only — solution gas / gas-cap expansion / PVT coupling stay in petekSim.
- **The static lift (2026-07-03, `task_relocate_refine_orchestration`).**
  petekStatic now owns volumetrics + static uncertainty (graph
  `decision_layer_charters`):
  - **`srs-volumetrics` + `srs-uncertainty` relocated from petekSim** (origin SHA
    `fe6343c`): GRV/in-place (OOIP/OGIP) and the Monte-Carlo toolkit
    (SplitMix64, inverse-CDF distributions, P90/P50/P10). FVF value types
    (`OilFvf`/`GasFvf`) are duplicated locally — no PVT crosses the seam.
  - **`srs-model` — the `StaticModel` contract, landed.** The aggregate
    (framework + grid + cubes + zones + contacts + provenance) with its own
    volumetrics output surface (`in_place()`); `StaticModelBuilder` (the
    relocated refine orchestration: flat-box or wireframe seed, log-driven
    population, live `add_top_control`); and the RATIFIED MC-regeneration seam:
    `StaticModelTemplate::new(wf, opts)` + `realize(&mut self, &RealizationDraw)`
    with the warm-start chain held in petekTools kernel space. `Send` on
    template + model (compile-checked). Smoke-tested at N=100 realizations.
  - **Physical-range validation (`validate_fraction`/`validate_positive`)** —
    φ/NG/Sw ∈ [0,1], positive magnitudes, typed errors instead of silent garbage;
    `PercentileSummary::from_realizations` returns a typed error on an empty set
    instead of panicking.
  - **`StaticError::Geo(#[from] petekio::GeoError)`** — `?` chains
    DATA→GEOMODEL→SIM across the seams.
  - **`KernelSurface` newtype** on `solve_surface_seeded` — seeding the warm
    kernel from a cold `solve_surface` output is now a compile error.
  - `Properties::names()` — cube-name enumeration.
- **Extraction from petekSim (2026-07-01).** petekStatic stood up as the GEOMODEL
  layer, a six-crate Rust workspace pulled out of petekSim: `petekstatic-error`
  (the single error enum), `srs-wireframe` + `srs-grid` (structural framework +
  corner-point grid construction), `srs-gridder` (convergent minimum-curvature
  gridder + stratigraphic layering), `srs-petro` (log upscaling), and `srs-data`
  (the petekIO → wireframe ingest adapter, carrying petekIO's neutral
  `Distribution` DTO). Cross-checks the gridder against petekTools' kernels.
- **`SPEC.md` + `API.md` — the StaticModel contract (2026-07-03).** Design
  constitution and provisional public API: `StaticModel` = structural framework +
  grid + property cubes + zones + contacts + provenance; the zones concept; and
  the petekSim seam (per-realization Monte-Carlo regeneration via a warm-started
  `StaticModelTemplate`).
- **`srs-gridder::solve_surface_seeded` — the warm-start refine path (SPEC §7a).**
  Re-grids a lattice from a prior converged `Surface` instead of a cold seed,
  delegating the seeded SOR to petekTools' `ConvergentGridder` kernel
  (`grid_min_curvature_seeded`). The load-bearing per-realization regeneration
  optimization: ~14x faster than a cold solve for a one-control perturbation at
  50x50 (94 ms vs 1.34 s). Not kernel-interchangeable with the cold
  `solve_surface` — see `API.md` for the boundary-treatment caveat and the open
  Q2 kernel-unification question.

### Changed
- **`solve_surface_seeded` now takes/returns `KernelSurface`** (was `Surface`) —
  the kernel-space constraint is enforced by the type; bootstrap warm chains via
  `KernelSurface::flat`.
- **`Realizations::summary()` / `PercentileSummary::from_realizations` are
  fallible** (`Result<_, StaticError>`) — empty realization sets are a typed
  error, not a panic (validation H1).
- **petekTools dependency: path → published `0.1` (crates.io).** The workspace now
  consumes petekTools `0.1.0` from crates.io rather than the sibling `../petekTools`
  path, insulating petekStatic from concurrent toolkit edits and stepping toward
  its own 0.1.0 release. petekTools is now a runtime dependency of `srs-gridder`
  (was a dev-only cross-check probe).
