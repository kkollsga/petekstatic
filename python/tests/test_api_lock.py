import petekstatic


def test_python_wheel_all_is_locked():
    assert petekstatic.__all__ == [
        "CoKriging",
        "DistributionSpec",
        "Grid",
        "Gridding",
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
        "WellLog",
        "WellLogSpec",
        "__version__",
        "build_flat_model",
        "distributions",
        "upscale",
    ]
