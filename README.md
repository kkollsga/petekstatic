# petekStatic

**GEOMODEL layer** of the petek subsurface-modelling ecosystem — a Rust library
that turns model-ready inputs into a populated **`StaticModel`** and owns the
**static-uncertainty** stack over it:

- **Structural framework** — horizons + zones (+ faults, planned).
- **Grid construction** — the convergent corner-point grid the properties live on.
- **Property modelling** — priors + log upscaling today; facies/geostatistics
  planned.
- **Volumetrics + static uncertainty** — GRV / in-place (OOIP/OGIP) off the
  model itself, and Monte-Carlo regeneration over static-model realizations
  (`StaticModelTemplate::realize(&RealizationDraw)`).

→ a populated `StaticModel` with its own volumes/P-curve surface, which
**petekSim** (dynamic/engineering simulation + the Python product facade)
consumes across the seam.

## Where it sits (deps flow one way, downward)

```
petekIO      DATA       → model-ready inputs (ModelInputs / .pproj)      [upstream]
   ↓
petekStatic  GEOMODEL   → populated StaticModel + volumetrics/uncertainty [here]
   ↓
petekSim     SIMULATION → dynamic/engineering + the Python product        [downstream]

petekTools   TOOLKIT    → gridding/kriging/warm-start kernels + units + pproj  [horizontal]
```

## Status — built and green (pre-0.1)

A nine-crate Cargo workspace: `petekstatic-error`, `srs-grid`, `srs-gridder`,
`srs-petro`, `srs-wireframe`, `srs-data`, `srs-volumetrics`, `srs-uncertainty`,
and `srs-model` (the `StaticModel` aggregate + the ratified MC-regeneration
seam). The volumetrics/uncertainty crates relocated here from petekSim on
2026-07-03 (the layer-charter re-scope). Design constitution: `SPEC.md`; public
contract: `API.md`. Working folders: `dev-docs/README.md`, `inbox/README.md`.

## Licensing

petekStatic is licensed under the **Business Source License 1.1** — see
[LICENSE](LICENSE). Non-production use is freely granted; production use is
permitted by the Additional Use Grant except as a competing commercial
"as-a-service" offering of the Licensed Work's functionality. Each released
version converts to the **Change License (Apache-2.0)** four years after its
first publication. For alternative licensing, contact kkollsg@gmail.com. (The
`{VERSION}` / Change Date parameters are filled in at each release cut.)
