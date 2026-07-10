# petekStatic — Claude Code Conventions

petekStatic is the **GEOMODEL layer** of the petek subsurface-modelling
ecosystem: a Rust library that turns model-ready inputs into a **populated
`StaticModel`** — the structural framework (horizons + faults + zones), grid
construction, property modelling (priors/log upscaling today; facies,
geostatistics, trend population planned) — **and owns volumetrics + static
uncertainty over it** (graph `decision_layer_charters`, 2026-07-03): GRV /
in-place (OOIP/OGIP) off the model itself, Monte-Carlo regeneration over
static-model realizations (`StaticModelTemplate::realize`), tornado later.

Its place in the one-way DAG (deps flow **downward only**):

```
petekIO      DATA       → model-ready inputs (ModelInputs / .pproj)   [optional adapter]
   ↓
petekStatic  GEOMODEL   → populated StaticModel + volumetrics/uncertainty    [THIS LIBRARY]
   ↓
petekSim     SIMULATION → dynamic/engineering + the Python product facade    [downstream consumer]

petekTools   TOOLKIT    → numeric kernels (gridding/kriging/warm-start) + units + pproj container  [horizontal dep]
```

The two committed sources of truth are **`SPEC.md`** (design constitution +
architecture) and **`API.md`** (the public API contract; locks fully at 0.1) at
the repo root — the petek family house style (canonically
`petekSuite/dev-docs/petek-house-style.md`), same as petekIO/petekSim. The
dev-docs + inbox + skills system below is local working state — see
`dev-docs/README.md` and `inbox/README.md` for the canonical maps.

## Where petekStatic is today

A **nine-crate Cargo workspace, built and green** (extracted from petekSim
2026-07-01; volumetrics/uncertainty relocated in 2026-07-03's static lift):
`petekstatic-error` · `srs-grid` · `srs-gridder` · `srs-petro` · `srs-wireframe`
· `srs-data` · `srs-volumetrics` · `srs-uncertainty` · `srs-model` (top of DAG:
the `StaticModel` aggregate + builder + the ratified MC-regeneration seam).
petekSim's `srs-core` consumes these across the repo seam as path deps and keeps
only a thin refine facade + the analytic box path.

The repo is a **local git repo without a GitHub remote** — `phased-plan` runs
with local commits (branch/PR/CI steps activate when a remote exists);
`release` activates at the 0.1 cut (`task_petekstatic_release_0_1`). The MC
driver over the template (`task_peteksim_mc_structured`, retargeted here) and
tornado (`task_peteksim_tornado`) are petekStatic's next uncertainty work.

## Data — test against `a local real-dataset folder`, never leak it into the repo

**You are allowed to test against real subsurface datasets under
`a local real-dataset folder`** — read them, drive them through petekStatic,
build local eval harnesses. **But never let their contents leak into the repo.**
No information derived from a dataset's *contents* (coordinates, values,
well/field names, survey rows, log/grid samples) may land in committed code,
fixtures, tests, examples, docs, commit messages, `CHANGELOG.md`, an inbox
message, the planning graph, or any published/exported output. Reference a
dataset by *path* and *format*, never by content.

- **Committed tests/examples use SYNTHETIC data** — hand-authored to format
  spec. Real-data evaluation happens in a **harness that lives in the data
  folder**, whose output also stays there (print structure/counts, never values).
- The published crate ships **no** test/example data (`Cargo.toml` `exclude`s
  `/tests` + `/examples`); keep it that way.

## Working style

- **Keep each response under 400 tokens.** For any long output, write it to a
  file (`dev-docs/temp/`, >1-day purge) and tell me the path instead of printing
  it.
- **Reproduce before fixing / claiming.** Before changing code, reproduce the
  issue and confirm the exact root cause with evidence. Don't apply a fix until
  the cause is verified. Confirm cross-library facts against the graph or the
  library itself, not assumption.
- **`API.md` is a contract, not a suggestion** (once it exists). Implement toward
  its exact signatures; changing a signature needs sign-off (see `API.md`'s
  header) and an edit to `API.md` itself — never let the code silently drift.

## Code analysis

- **Use the code-review MCP for code analysis** — `set_root_dir` to this repo,
  then Cypher over the code graph + ripgrep (`grep`). Prefer it over ad-hoc file
  reads when mapping structure, finding callers, or tracing the layered deps.
  Pair it with `SPEC.md`, the upstream contracts (petekIO's `API.md`,
  petekTools' kernel signatures), and the suite planning graph.
- **Spin up `Explore` agents** to parallelize broad sweeps and keep the
  conclusions, not the file dumps, in context.

## Architecture — the charter (from the suite)

The internal design constitution is `SPEC.md`. What's decided (by the suite, in
the graph):

- **Consumes, never reaches up.** The geomodel core calls **petekTools** numeric
  kernels (gridding/kriging/warm-start) and stays independent of petekIO.
  `petekio-adapter` is an opt-in compatibility seam for `ModelInputs`; the
  integrated conversion lives at petekSim, which already depends on both
  libraries. petekStatic **never** depends on petekSim. No cycles or sideways
  sharing — convert small types at the composition seam.
- **Produces a `StaticModel` that owns its volumes** (`decision_layer_charters`):
  framework + grid + cubes + zones + contacts + provenance, with `in_place()`
  and the MC-regeneration seam (`StaticModelTemplate::realize(&RealizationDraw)`,
  ratified `decision_staticmodel_regen_seam`) on top. petekSim's facade presents
  the results; FVF crosses as a validated scalar input, never PVT code. The
  contract lives in `SPEC.md`/`API.md` and the suite graph.
- **The kernel constraint:** the template's warm-start chain lives in petekTools
  kernel space (`KernelSurface`; `decision_gridder_kernel_unification`) — never
  seed the warm kernel from a cold `solve_surface` output.
- **House-style conventions (family-wide):** strictly layered one-way internal
  deps; a manager/substrate collection with broadcast ops (no per-item loops);
  domain objects carry their operations (fluent, chainable, immutable — ops
  return *new* objects, mutation is explicit `set_*`); open/closed (extend by
  adding, not editing); compartmentalized (one module/topic, one
  type/responsibility); compose deps, don't reinvent; `f64::NAN` = undefined;
  one error enum (`thiserror`) + `Result<T, E>` everywhere; Rust core + *thin*
  PyO3 (bindings only marshal). The specifics are locked in `SPEC.md`.

## Build & test

```bash
cargo build --all-features
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace             # Rust unit + golden/analytic tests
cargo bench                         # criterion; release build, for perf only
# Python: no petekStatic-side bindings yet — the peteksim wheel (petekSim repo)
# is the Python surface; changes here gate on petekSim's `maturin develop` too.
```

**Tooling discipline (don't relearn it the hard way):**
- **Never read a gate's status through a `tail`/`head` pipe** — a pipeline's
  exit code is the *last* command's. Run the gate bare, or `set -o pipefail`, or
  `cmd; echo "exit=$?"`.
- **After `maturin develop`, confirm it printed `Installed`.** A build error
  upstream leaves the old `.so` in place → the next `pytest` silently tests
  stale code.
- **After any Rust *behaviour* change, run `cargo test`, not just `pytest`** —
  the golden/analytic assertions live on the Rust side.

## Code health

Each pass through a file should leave it more compartmentalized than you found it.

- **No bugs left behind.** Fix a pre-existing bug you encounter in the same
  change, or surface it explicitly (a `todos.md` item) rather than stepping over
  it. First confirm it's a real defect, not deliberate behaviour (read the
  surrounding code/tests, check against `SPEC.md`/`API.md`).
- **Golden / analytic tests are the safety net.** A correctness path lands with
  the test that proves it — grid construction vs a known geometry, property
  population statistics vs an analytic expectation, log upscaling vs a hand calc,
  a `StaticModel` volumetric vs an analytic value. A new modelling step without
  such a check isn't trusted.
- **Fixing a bug — scan for the *class*.** Probe with scratch fixtures
  (`dev-docs/temp/`) before declaring scope.
- A measured perf change is only a "fix" if it measurably improves perf.

## Testing doctrine — the six rules (family-wide, `petekSuite/dev-docs/designs/testing-doctrine.md`)

Every rule names the bug class that proves its necessity; each is derived from an
actual escape. A test that violates a rule is an incomplete test.

- **R1 — Frame rule.** Every frame-sensitive test (views, population, ties,
  sections, maps) has a **world-georeferenced** variant (local-origin lattice +
  Georef + world-coordinate inputs, fictional 431000/6521000-style coords). A
  local-frame-only fixture is incomplete. *(3× world/local seam bugs.)*
- **R2 — Mode-matrix rule.** A model-level feature is tested across the mode matrix
  it claims to support — **in-core × spilled, serial × sharded, single-zone ×
  horizon-stack, wireframe × horizon-stack construction**. The matrix is **declared
  in the test-module header** (which cells the file covers); a SUPPORTED cell has ≥1
  test, an UNSUPPORTED cell is a **documented typed-error test** (the spilled
  `zone_stats` v1 gap is the template), never an untested hole. The cross-feature
  matrix lives in `srs-model/tests/mode_matrix.rs`. *(Volume bundle empty / map +
  section bundle broken on spilled models — all read the 1×1×1 placeholder.)*
- **R3 — Planted-truth rule.** Every modelling capability (population, trends,
  upscaling, MC, volumetrics) carries a planted-truth recovery test on the synthetic
  asset: plant a known value, recover it through the full pipeline, assert within a
  derived tolerance. Zero-spread MC == deterministic is the canonical instance.
- **R4 — Loudness rule.** No silent degradation, and the loudness is **tested**:
  every fallback branch either raises a typed error naming property+zone/frame or
  emits a warning the test asserts. A fallback without a loudness test is a defect.
  *(Collocated-coverage error, mean-fill error, clip-miss assert, non-finite-surface
  guard are the pattern.)*
- **R5 — Degenerate-input rule.** Kernels that iterate/converge (layering, collapse,
  order-repair, SGS, solvers) get adversarial-input **proptests**: zero /
  sub-threshold / inverted / NaN columns, empty conditioning, single-cell + single-
  layer grids — **plus a hard per-case timeout so a livelock FAILS** (times out)
  instead of hanging CI. `proptest` is a dev-dependency; convergence-loop proptests
  wrap the kernel on a worker thread with `recv_timeout`. *(collapse_below_m
  livelock — two distinct instances.)*
- **R6 — Round-trip rule.** A cross-repo feature is done only when the end-to-end
  acceptance suite passes on the canonical synthetic asset with payload invariants
  asserted (l≠r on dipping cells, non-empty watertight shells, outline == frame
  extent, ties populated) and the planted truths recovered. Homed in peteksim's test
  tree; the coordinator does not stamp a cross-repo task without it.

## Performance protocol

Before any perf-related change: baseline first (write/extend a criterion bench,
record numbers); **release build only**; trust `min` over `median` for sub-ms
benches. Heavy built grids/cubes → `dev-docs/bench/out/`; the regression rows →
`dev-docs/bench/results/results.csv`. See `dev-docs/bench/README.md`.

## Inbox hygiene

**Always use the inbox skills for cross-library communication** — never
hand-read or hand-write inbox files. Incoming → **`read-inbox`** (triage
`unread/`, lift durable info to `dev-docs/` + lean `todos.md` backlinks, route,
archive, purge). Outgoing → **`notify`** (resolve the target under `Koding/`,
compose per the schema, drop into its `inbox/unread/`). The canonical map is
`inbox/README.md`. Natural correspondents:

- **petekIO** — upstream dependency; the `ModelInputs` / `.pproj` seam we consume.
- **petekTools** — upstream horizontal dep; the numeric kernels + pproj container.
- **petekSim** — downstream consumer of the `StaticModel`.
- **petekSuite** — the coordinator; route cross-library initiatives + planning-
  graph contributions here.
- **mcp-servers** — the whole ecosystem, one inbox (never resolve a name to
  `mcp-servers/<subdir>/`).

## Planning graph — the cross-library source of truth

The petek **planning graph** (served by the `contract` MCP; homed at
`petekSuite/research/graph/research.kgl` — the coordinator) is the single source
of truth for the inter-library contracts
(the `ModelInputs` seam, the `StaticModel` seam, the layered architecture),
decisions, and open questions. Reach for it on anything cross-cutting — read the
contract before changing a shared seam; record blocking issues and choices
there, not only in local docs. Contribute **without cluttering**: runtime types
only (`Question` / `Decision` / `Artifact` / `Task` — never the managed research
nodes); **MERGE on id, never CREATE**; one node per concept; stamp `git_sha` +
`modified_by='petekstatic'`. No direct graph access → **route it through the
inbox to petekSuite** (the coordinator), who curates it in.

## Commits & releases

Commit format: `type: short description` (`feat`, `fix`, `docs`, `refactor`,
`test`, `chore`). Update `CHANGELOG.md` `[Unreleased]` for user-visible changes;
skip for internal refactors, CI, test-only, formatting.

**Pushing requires explicit, in-the-moment approval.** Default is *don't push*.
Approval is one-shot — it covers exactly that one `git push` and does not carry
to a later commit or branch.

**Exception — invoking the `release` skill IS push authorization for that
release** (the publish-triggering `main` push + its CI fix-and-push loop),
scoped to that one run. Every pre-push safeguard still applies. (First release =
`task_petekstatic_release_0_1`; a GitHub remote lands with it.)

Version source of truth: root `Cargo.toml` `[workspace.package] version` (or the
single crate's `version`) — one bump per push, all workspace members in lockstep.

## The skills (wired into these rules)

- **`phased-plan`** — run any non-trivial, multi-step change as gated phases.
  Don't use generic plan mode for large work.
- **`add-todo`** — the single authority on `todos.md` entry shape; capture work
  as a lean backlink + a `plans/` detail doc.
- **`dev-docs-cleanup`** — purge the time-boxed dirs + a todos-driven tidy. Run
  before a new phased-plan and at the end of a release.
- **`read-inbox`** / **`notify`** — the receive / send sides of the inbox.
- **`release`** — ship: goal-check, gate, reconcile `API.md`, bump, promote
  CHANGELOG, publish, tidy. Run only when asked.
