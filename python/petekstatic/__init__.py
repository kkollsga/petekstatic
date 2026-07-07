"""petekStatic — the geomodel layer of the petek subsurface-modelling stack.

The compiled extension exposes the compact ``StaticModel`` surface. The Python
workflow facade adds the canonical notebook-facing declaration style:
``Grid.from_project(...).geometry(...).horizons(...).zones(...).layers(...)``,
in-memory property formula calculation, and deterministic smoke-test volumes.
See https://github.com/kkollsga/petekstatic.
"""

from ._petekstatic import (
    PropertyPipeline,
    StaticModel,
    WellLog,
    __version__,
    build_flat_model,
)
from .workflow import (
    CoKriging,
    DistributionSpec,
    Grid,
    HorizonSpec,
    Layering,
    PropertyHandle,
    PropertyPipelineSpec,
    PropertyStore,
    SgsRecipe,
    Spherical,
    UpscaleRecipeBuilder,
    Var,
    VolumeCase,
    VolumeResult,
    WellTie,
    WellLogSpec,
    distributions,
    upscale,
)

__all__ = [
    "CoKriging",
    "DistributionSpec",
    "Grid",
    "HorizonSpec",
    "Layering",
    "PropertyHandle",
    "PropertyPipeline",
    "PropertyPipelineSpec",
    "PropertyStore",
    "SgsRecipe",
    "Spherical",
    "StaticModel",
    "UpscaleRecipeBuilder",
    "Var",
    "VolumeCase",
    "VolumeResult",
    "WellTie",
    "WellLog",
    "WellLogSpec",
    "__version__",
    "build_flat_model",
    "distributions",
    "upscale",
]
