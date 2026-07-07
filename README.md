# petekStatic

**GEOMODEL layer** of the petek subsurface-modelling ecosystem — a Rust library
that turns model-ready inputs into a populated **`StaticModel`** and owns the
**static-uncertainty** stack over it:

- **Structural framework** — horizons + zones (+ faults, planned).
- **Grid construction** — the convergent corner-point grid the properties live on.
- **Property modelling** — priors, property formulas, log upscaling, and the
  first `PropertyPipeline` lowering path; facies/geostatistics grow from there.
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

## Python workflow facade

petekStatic now owns the first Python slice of the canonical static workflow.
The facade is deliberately small: declare the grid, assign properties, calculate
formulas through petekTools, and run deterministic smoke volumes.

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
p.ntg = 0.80                         # scalar assignment broadcasts to the grid
p.por = p.ntg * 0.85                 # assignment expression via petekTools
p.sw = 0.20
p["PermXY_BC"].set(100.0)
p["PorE_BC"].set(0.25)
p.calc(["RQI = $lambda * sqrt(PermXY_BC / PorE_BC)"], params={"lambda": 0.0314})

result = grid.volumes(ntg="NTG", por="POR", sw="Sw", fluid="oil", fvf=1.30).run(progress=True)
```

Property recipes use the same facade and lower to the Rust `PropertyPipeline`
shape:

```python
logs = project.logs
vgm = pst.Var("spherical", major=1500, minor=700, vertical=20, azimuth=35)

p.por = pst.upscale(logs.PHIE(logs.NetSand > 0.50)).sgs(
    variogram=vgm,
    distribution=pst.distributions.from_logs(),
    seed=12,
)

spec = p.pipelines("por")            # serialization-ready PropertyPipelineSpec
```

If the project or source resolves positioned wells, `p.execute_pipeline("por")`
returns a Rust-backed `pst.PropertyPipeline` handle:

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

Applying that handle to an arbitrary production grid is not exposed yet;
`PropertyPipeline.apply_to_flat_model(...)` is the current smoke execution path.
Lazy logs, cokriging/trend binding, non-`from_logs` distributions, and
anisotropic Rust execution are guarded explicitly. Anisotropic `pst.Var(...)`
specs are still preserved and can be lowered to the petekTools anisotropic
variogram object.

## Status — built and green

A single crate, `petekstatic` (0.1.7), whose modules preserve the historical
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
