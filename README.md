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

## Documentation

The canonical docs for the whole petek family live on the **petekSuite site**
— petekStatic's pages there:

- **[Library guide](https://peteksuite.readthedocs.io/en/latest/libraries/petekstatic/)** — the petekStatic guide.
- **Tutorial** — [Static model build (flagship)](https://peteksuite.readthedocs.io/en/latest/tutorials/static-model-build/).
- **[Notebooks](https://peteksuite.readthedocs.io/en/latest/notebooks/)** — executed examples: [stack model from scatter](https://peteksuite.readthedocs.io/en/latest/notebooks/petekstatic/01_stack_model_from_scatter/) · [volumes & bundles](https://peteksuite.readthedocs.io/en/latest/notebooks/petekstatic/02_volumes_and_bundles/).

## Where it sits (deps flow one way, downward)

```
petekIO      DATA       → model-ready inputs (ModelInputs / .pproj)      [upstream]
   ↓
petekStatic  GEOMODEL   → populated StaticModel + volumetrics/uncertainty [here]
   ↓
petekSim     SIMULATION → dynamic/engineering + the Python product        [downstream]

petekTools   TOOLKIT    → gridding/kriging/warm-start kernels + units + pproj  [horizontal]
```

## Status — built and green

A single crate, `petekstatic` (0.1.1), whose modules preserve the historical
layer boundaries and one-directional imports: `petekstatic::{error, wireframe,
grid, petro, gridder, volumetrics, uncertainty, data, spill, model}` — with the
top-of-DAG `model` surface (the `StaticModel` aggregate + the ratified
MC-regeneration seam) re-exported at the crate root. Consolidated from the
former ten-crate workspace on 2026-07-05 (owner ruling; the former `srs-*` and
`petekstatic-error` package names are retired). The volumetrics/uncertainty code
relocated here from petekSim on 2026-07-03 (the layer-charter re-scope). Design
constitution: `SPEC.md`; public contract: `API.md`. Working folders:
`dev-docs/README.md`, `inbox/README.md`.

## Licensing

petekStatic is licensed under the **Business Source License 1.1** — see
[LICENSE](LICENSE). Non-production use is freely granted; production use is
permitted by the Additional Use Grant except as a competing commercial
"as-a-service" offering of the Licensed Work's functionality. Each released
version converts to the **Change License (Apache-2.0)** four years after its
first publication. For alternative licensing, contact kkollsg@gmail.com. (The
`{VERSION}` / Change Date parameters are filled in at each release cut.)
