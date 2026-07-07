from __future__ import annotations

import os
import statistics
import time
from pathlib import Path
from typing import Callable, Any

import pytest

import petekstatic as pst

petekio = pytest.importorskip("petekio")


SYNTH_FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "fixtures"
    / "wells"
)
CONFIDENTIAL_DATA = Path(
    os.environ.get(
        "PETEK_CONFIDENTIAL_DATA",
        "/Volumes/EksternalHome/Data/modellingProject/Data",
    )
)

ALIASES = {
    "PHIE": ["PHIE", "PHI", "PORO", "POR", "POROSITY"],
    "NetSand": ["NetSand", "NETSAND", "NTG", "NET", "VSH_NET"],
    "SW": ["SW", "SWT", "SUWI"],
}


def _bench_ms(fn: Callable[[], Any], *, n: int) -> dict[str, float]:
    values: list[float] = []
    for _ in range(n):
        start = time.perf_counter()
        fn()
        values.append((time.perf_counter() - start) * 1000.0)
    return {
        "min": min(values),
        "median": statistics.median(values),
        "mean": statistics.mean(values),
        "max": max(values),
    }


def _assert_under(label: str, value_ms: float, budget_ms: float) -> None:
    assert value_ms < budget_ms, (
        f"{label} regressed: {value_ms:.3f} ms exceeds budget {budget_ms:.3f} ms"
    )


def _build_recipe(project, *, module=pst):
    logs = project.wells.logs
    return module.upscale(logs.PHIE(logs.NetSand >= 0.50)).sgs(
        variogram=module.Var(
            "spherical",
            major=1500.0,
            minor=700.0,
            vertical=20.0,
            azimuth=35.0,
        ),
        distribution=module.distributions.from_logs(),
        seed=12,
    )


@pytest.mark.perf
def test_static_property_workflow_synthetic_perf_budget():
    """Cheap release-gate tripwire for the canonical workflow seam.

    This uses the tiny public petekIO fixture so it can run on every machine.
    Budgets are intentionally generous wall-clock ceilings: they catch accidental
    O(N) loops or cache bypasses without pretending to be micro-benchmarks.
    """

    project_load = _bench_ms(
        lambda: petekio.Project.load(SYNTH_FIXTURE, aliases=ALIASES),
        n=10,
    )
    project = petekio.Project.load(SYNTH_FIXTURE, aliases=ALIASES)
    recipe = _build_recipe(project, module=pst)

    first_lower = _bench_ms(lambda: recipe.lower("POR", project=project), n=1)
    spec = recipe.lower("POR", project=project)
    cached_lower = _bench_ms(lambda: recipe.lower("POR", project=project), n=100)
    construct = _bench_ms(lambda: spec.execute(), n=20)

    assert len(spec.well_logs or ()) == 1
    assert sum(len(w.samples) for w in (spec.well_logs or ())) == 5
    pipe = spec.execute()
    assert pipe.config()["minor_m"] == 700.0
    assert pipe.config()["vertical_m"] == 20.0

    _assert_under("synthetic Project.load min", project_load["min"], 25.0)
    _assert_under("synthetic first lower", first_lower["min"], 5.0)
    _assert_under("synthetic cached lower median", cached_lower["median"], 0.75)
    _assert_under("synthetic pipeline construct min", construct["min"], 2.0)


@pytest.mark.perf
def test_static_property_workflow_confidential_perf_budget():
    """Real-data release gate for the canonical property workflow.

    The default path is the local confidential data folder. CI and clean agents
    without that folder skip; agents with it assert the same budgets without
    printing or committing confidential file names or raw data values.
    """

    if not CONFIDENTIAL_DATA.exists():
        pytest.skip(
            "confidential data folder is not available; set PETEK_CONFIDENTIAL_DATA"
        )

    project_load = _bench_ms(
        lambda: petekio.Project.load(CONFIDENTIAL_DATA, aliases=ALIASES),
        n=3,
    )
    project = petekio.Project.load(CONFIDENTIAL_DATA, aliases=ALIASES)
    recipe = _build_recipe(project, module=pst)

    first_lower = _bench_ms(lambda: recipe.lower("POR", project=project), n=1)
    spec = recipe.lower("POR", project=project)
    cached_lower = _bench_ms(lambda: recipe.lower("POR", project=project), n=200)
    construct = _bench_ms(lambda: spec.execute(), n=20)

    assert len(spec.well_logs or ()) >= 1
    assert sum(len(w.samples) for w in (spec.well_logs or ())) >= 1
    pipe = spec.execute()
    assert pipe.config()["major_m"] == 1500.0
    assert pipe.config()["minor_m"] == 700.0
    assert pipe.config()["vertical_m"] == 20.0
    assert pipe.config()["azimuth"] == 35.0

    _assert_under("confidential Project.load min", project_load["min"], 750.0)
    _assert_under("confidential first lower", first_lower["min"], 60.0)
    _assert_under("confidential cached lower median", cached_lower["median"], 1.0)
    _assert_under("confidential pipeline construct min", construct["min"], 3.0)
