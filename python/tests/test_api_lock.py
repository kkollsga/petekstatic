import petekstatic


def test_python_wheel_all_is_locked():
    assert petekstatic.__all__ == [
        "StaticModel",
        "__version__",
        "build_flat_model",
    ]
