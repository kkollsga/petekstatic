import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest

import petekstatic as pst


@dataclass(frozen=True)
class Surface:
    values: list[float]


class SyntheticProject:
    def __init__(self):
        self.surfaces = {
            "Surfaces.Top reservoir": Surface([1000.0, 1000.0]),
            "Surfaces.Base reservoir": Surface([1020.0, 1010.0]),
            "Top reservoir input surface": Surface([900.0, 900.0]),
            "Base reservoir input surface": Surface([930.0, 920.0]),
            "input surface": Surface([950.0, 940.0]),
        }
        self.polygons = {"Polygons.ModelEdge": [[0.0, 0.0], [1.0, 0.0]]}
        self.tops = ["Top reservoir", "Base reservoir"]

    def inventory(self):
        return {
            "surfaces": list(self.surfaces),
            "polygons": list(self.polygons),
            "tops": list(self.tops),
            "points": [],
            "wells": [],
        }


@dataclass(frozen=True)
class LogExpr:
    mnemonic: str
    predicate: str | None = None

    def __call__(self, predicate):
        return LogExpr(self.mnemonic, str(predicate))

    def __gt__(self, cutoff):
        return f"{self.mnemonic} > {cutoff}"

    def as_dict(self):
        out = {"kind": "log", "mnemonic": self.mnemonic}
        if self.predicate is not None:
            out["predicate"] = self.predicate
        return out


@dataclass(frozen=True)
class Logs:
    NetSand: LogExpr = LogExpr("NetSand")
    PHIE: LogExpr = LogExpr("PHIE")


@dataclass(frozen=True)
class LazyPredicate:
    op: str
    left: str
    right: float

    def as_dict(self):
        return {
            "kind": "log_predicate",
            "op": self.op,
            "operands": [
                {"kind": "log_channel", "mnemonic": self.left, "requested": self.left},
                {"kind": "literal", "value": self.right},
            ],
        }


@dataclass(frozen=True)
class LazyLogChannel:
    mnemonic: str
    filter: LazyPredicate | None = None

    def __call__(self, predicate):
        return LazyLogChannel(self.mnemonic, predicate)

    def as_dict(self):
        out = {
            "kind": "log_channel",
            "mnemonic": self.mnemonic,
            "requested": self.mnemonic,
        }
        if self.filter is not None:
            out["filter"] = self.filter.as_dict()
        return out


class ResolvingLogChannel(LazyLogChannel):
    def to_well_logs(self, project):
        assert project is not None
        return [
            {"x": 10.0, "y": 20.0, "samples": [(1000.0, 0.18), (1005.0, 0.22)]}
        ]


def declared_grid():
    return (
        pst.Grid.from_project(SyntheticProject())
        .geometry(cell=(10.0, 10.0), orient=0.0, outline="ModelEdge")
        .horizons(
            [
                {"name": "Top reservoir", "well top": "well tops/Top reservoir"},
                "Base reservoir",
            ],
            well_tie={"influence_radius": 800},
        )
        .zones({"Reservoir": ("Top reservoir", "Base reservoir")})
        .layers({"Reservoir": pst.Layering(n=2)})
    )


def test_static_workflow_facade_calculates_properties_and_volumes():
    grid = declared_grid()
    p = grid.properties
    p["PermXY_BC"].set([100.0, 100.0, 100.0, 100.0])
    p["PorE_BC"].set([0.25, 0.25, 0.25, 0.25])
    p["HA_FWL"].set([1.0, 1.0, 1.0, 1.0])
    p["Jfunc"].set([1.0, 1.0, 1.0, 1.0])
    p.ntg.set([0.8, 0.8, 0.8, 0.8])
    p.por.set([0.25, 0.25, 0.25, 0.25])

    out = p.calc(
        [
            "RQI = $lambda * sqrt(PermXY_BC / PorE_BC)",
            "Swirr = $SHF_c * pow(RQI, $SHF_d)",
            "Sw = if(HA_FWL == 0, 1, Swirr + (1 - Swirr) * $SHF_a * pow(Jfunc, $SHF_b))",
        ],
        params={
            "lambda": 0.01,
            "SHF_a": 0.0,
            "SHF_b": 1.0,
            "SHF_c": 1.0,
            "SHF_d": 1.0,
        },
    )

    assert out["RQI"] == pytest.approx([0.2, 0.2, 0.2, 0.2])
    assert p["Sw"].values == pytest.approx([0.2, 0.2, 0.2, 0.2])

    events = []
    case = grid.volumes(ntg="NTG", por="POR", sw="Sw", fluid="oil", fvf=1.20)
    result = case.run(progress=events.append)

    assert [event["stage"] for event in events] == ["structure", "volumes", "complete"]
    assert result.summary()["grv_m3"] == pytest.approx(3000.0)
    assert result.summary()["hcpv_m3"] == pytest.approx(480.0)
    assert result.summary()["ooip_sm3"] == pytest.approx(400.0)
    assert result.by_zone()["total"]["cells"] == 4


def test_property_calc_accepts_mapping_and_volume_aliases():
    grid = declared_grid()
    grid.properties.ntg.set([0.5, 0.5])
    out = grid.properties.calc({"por": "ntg * $scale", "Sw": "0.25"}, params={"scale": 0.4})

    assert out["POR"] == pytest.approx([0.2, 0.2])
    assert out["Sw"] == pytest.approx([0.25, 0.25])

    result = grid.volumes(ntg="ntg", por="por", sw="Sw", fvf=1.0).run()
    assert result.summary()["hcpv_m3"] > 0.0


def test_property_set_broadcasts_scalars_to_declared_cells():
    grid = declared_grid()
    p = grid.properties

    p.ntg.set(0.8)
    p.por = 0.25
    p.sw.set(0.2)

    assert p.ntg.values == pytest.approx([0.8, 0.8, 0.8, 0.8])
    assert p.por.values == pytest.approx([0.25, 0.25, 0.25, 0.25])
    assert p.sw.values == pytest.approx([0.2, 0.2, 0.2, 0.2])

    result = grid.volumes(ntg="ntg", por="por", sw="sw", fvf=1.0).run()
    assert result.summary()["cells"] == 4
    assert result.summary()["hcpv_m3"] == pytest.approx(480.0)


def test_property_assignment_accepts_handle_expressions():
    grid = declared_grid()
    p = grid.properties

    p.ntg = 0.8
    p.por = p.ntg * 0.85
    p.sw = 0.2

    assert p.por.values == pytest.approx([0.68, 0.68, 0.68, 0.68])

    result = grid.volumes(ntg="ntg", por="por", sw="sw", fvf=1.0).run()
    assert result.summary()["hcpv_m3"] == pytest.approx(1305.6)


def test_variogram_distribution_and_recipe_specs_are_serializable():
    logs = Logs()
    vgm = pst.Var(
        "Spherical",
        major=1500,
        minor=700,
        vertical=20,
        azimuth=395,
        nugget=0.05,
    )
    recipe = pst.upscale(logs.PHIE(logs.NetSand > 0.5)).sgs(
        variogram=vgm,
        distribution=pst.distributions.from_logs(),
        seed=12,
        cokriging=pst.CoKriging({"map": "PoroTrend"}, rho=0.7),
    )

    assert vgm.as_dict() == {
        "kind": "variogram",
        "model": "spherical",
        "major": 1500.0,
        "minor": 700.0,
        "vertical": 20.0,
        "azimuth": 35.0,
        "nugget": 0.05,
    }
    assert recipe.as_dict() == {
        "kind": "sgs",
        "upscale": {
            "source": {
                "kind": "log",
                "mnemonic": "PHIE",
                "predicate": "NetSand > 0.5",
            },
            "method": "arithmetic",
        },
        "variogram": vgm.as_dict(),
        "distribution": {"kind": "from_logs"},
        "seed": 12,
        "cokriging": {
            "kind": "cokriging",
            "trend": {"map": "PoroTrend"},
            "rho": 0.7,
        },
        "options": {},
    }


def test_recipe_assignment_records_property_declaration_without_execution():
    grid = declared_grid()
    logs = Logs()
    p = grid.properties
    p.ntg = 0.8
    recipe = pst.upscale(logs.NetSand).sgs(
        variogram=pst.Spherical(range_m=1500),
        distribution=pst.distributions.normal(mean=0.6, std=0.1),
        seed=11,
    )

    p.ntg = recipe

    assert "NTG" in p.names()
    assert p.declarations("ntg") == [{"kind": "recipe", "args": recipe.as_dict()}]
    with pytest.raises(ValueError, match="missing grid property 'NTG'"):
        _ = p.ntg.values


def test_recipe_assignment_lowers_to_property_pipeline_spec():
    grid = declared_grid()
    p = grid.properties
    source = LazyLogChannel("PHIE")(LazyPredicate(">", "NetSand", 0.5))
    vgm = pst.Var(
        "spherical",
        major=1500,
        minor=700,
        vertical=20,
        azimuth=35,
        sill=1.2,
        nugget=0.05,
    )

    p.por = pst.upscale(source).sgs(variogram=vgm, seed=12)

    lowered = p.pipelines("por")
    assert lowered["kind"] == "property_pipeline"
    assert lowered["property"] == "POR"
    assert lowered["upscale"] == {
        "rust_type": "PropertyPipeline::upscale",
        "source": {
            "kind": "log_channel",
            "mnemonic": "PHIE",
            "requested": "PHIE",
            "filter": {
                "kind": "log_predicate",
                "op": ">",
                "operands": [
                    {
                        "kind": "log_channel",
                        "mnemonic": "NetSand",
                        "requested": "NetSand",
                    },
                    {"kind": "literal", "value": 0.5},
                ],
            },
        },
        "well_logs": None,
        "method": "arithmetic",
    }
    assert lowered["propagate"]["rust_type"] == "Gaussian"
    assert lowered["propagate"]["variogram"] == vgm.as_dict()
    assert lowered["propagate"]["distribution"] == {"kind": "from_logs"}
    assert lowered["propagate"]["seed"] == 12


def test_recipe_lowering_resolves_positioned_well_logs_when_source_supports_it():
    grid = declared_grid()
    p = grid.properties
    p.ntg = pst.upscale(ResolvingLogChannel("NetSand")).sgs(
        variogram=pst.Var("spherical", major=900, minor=500, vertical=15, azimuth=10),
        distribution=pst.distributions.from_logs(),
        seed=11,
    )

    lowered = p.pipelines("ntg")
    assert lowered["upscale"]["well_logs"] == [
        {"x": 10.0, "y": 20.0, "samples": [[1000.0, 0.18], [1005.0, 0.22]]}
    ]


def test_lowered_variogram_builds_petektools_anisotropic_object(monkeypatch):
    captured = {}

    class FakeAnisotropicVariogram:
        def __init__(self, model, **kwargs):
            captured["model"] = model
            captured.update(kwargs)

    monkeypatch.setitem(
        sys.modules,
        "petektools",
        SimpleNamespace(AnisotropicVariogram=FakeAnisotropicVariogram),
    )
    grid = declared_grid()
    grid.properties.por = pst.upscale(ResolvingLogChannel("PHIE")).sgs(
        variogram=pst.Var(
            "spherical",
            major=1500,
            minor=700,
            vertical=20,
            azimuth=35,
            sill=1.2,
            nugget=0.05,
        ),
        distribution=pst.distributions.from_logs(),
        seed=12,
    )

    grid.properties.pipeline_spec("por").to_petektools_variogram()

    assert captured == {
        "model": "spherical",
        "major": 1500.0,
        "minor": 700.0,
        "vertical": 20.0,
        "azimuth": 35.0,
        "sill": 1.2,
        "nugget": 0.05,
    }


def test_recipe_execution_boundary_builds_rust_pipeline_handle():
    grid = declared_grid()
    grid.properties.por = pst.upscale(LazyLogChannel("PHIE")).sgs(
        variogram=pst.Var("spherical", major=1500, minor=700, vertical=20, azimuth=35),
        distribution=pst.distributions.from_logs(),
        seed=12,
    )

    with pytest.raises(NotImplementedError, match="positioned WellLog inputs"):
        grid.properties.execute_pipeline("por")

    grid.properties.ntg = pst.upscale(ResolvingLogChannel("NetSand")).sgs(
        variogram=pst.Var(
            "spherical",
            major=1500,
            minor=1500,
            vertical=1500,
            azimuth=0,
            sill=1.2,
            nugget=0.05,
        ),
        distribution=pst.distributions.from_logs(),
        seed=11,
    )
    pipe = grid.properties.execute_pipeline("ntg")

    assert isinstance(pipe, pst.PropertyPipeline)
    assert pipe.name == "NTG"
    assert pipe.config() == {
        "property": "NTG",
        "well_count": 1,
        "method": "arithmetic",
        "variogram_model": "spherical",
        "range_m": 1500.0,
        "major_m": 1500.0,
        "minor_m": 1500.0,
        "vertical_m": 1500.0,
        "azimuth": 0.0,
        "sill": 1.2,
        "nugget": 0.05,
        "seed": 11,
        "propagate": True,
        "allow_mean_fill": False,
        "max_neighbours": None,
        "radius_m": None,
        "unbounded_search": False,
    }
    assert "PropertyPipeline(property='NTG'" in pipe.report()


def test_recipe_execution_uses_petekio_project_log_resolver():
    petekio = pytest.importorskip("petekio")
    root = (
        Path(__file__).resolve().parents[3]
        / "petekIO"
        / "tests"
        / "fixtures"
        / "wells_petro"
    )
    if not root.is_dir():
        pytest.skip(f"petekIO fixture not available: {root}")
    project = petekio.Project.import_data(
        root,
        aliases={
            "PHIE": ["PHI", "PHIE"],
            "NetSand": ["NTG", "NETSAND"],
        },
    )
    logs = project.wells.logs

    recipe = pst.upscale(logs.PHIE(logs.NetSand >= 0.50)).sgs(
        variogram=pst.Var(
            "spherical",
            major=1500,
            minor=1500,
            vertical=1500,
            azimuth=0,
        ),
        distribution=pst.distributions.from_logs(),
        seed=12,
    )
    spec = recipe.lower("POR", project=project)

    assert spec.well_logs is not None
    assert len(spec.well_logs) == 1
    assert spec.well_logs[0].samples == (
        (2400.0, 0.2),
        (2410.0, 0.05),
        (2420.0, 0.2),
        (2430.0, 0.2),
        (2440.0, 0.2),
    )

    pipe = spec.execute()
    assert isinstance(pipe, pst.PropertyPipeline)
    assert pipe.config()["property"] == "POR"
    assert pipe.config()["well_count"] == 1


def test_recipe_execution_builds_anisotropic_rust_pipeline_handle():
    grid = declared_grid()
    grid.properties.por = pst.upscale(ResolvingLogChannel("PHIE")).sgs(
        variogram=pst.Var(
            "spherical",
            major=1500,
            minor=700,
            vertical=20,
            azimuth=35,
            sill=1.2,
            nugget=0.05,
        ),
        distribution=pst.distributions.from_logs(),
        seed=12,
    )

    pipe = grid.properties.execute_pipeline("por")

    assert pipe.config()["property"] == "POR"
    assert pipe.config()["major_m"] == 1500.0
    assert pipe.config()["minor_m"] == 700.0
    assert pipe.config()["vertical_m"] == 20.0
    assert pipe.config()["azimuth"] == 35.0
    assert "major_m=1500" in pipe.report()


def test_recipe_validation_errors_are_loud():
    logs = Logs()

    with pytest.raises(ValueError, match="Var.major"):
        pst.Var("spherical", major=0, minor=700, vertical=20, azimuth=35)
    with pytest.raises(ValueError, match="normal distribution std"):
        pst.distributions.normal(std=0)
    with pytest.raises(ValueError, match="CoKriging.rho"):
        pst.CoKriging("trend", rho=1.5)
    with pytest.raises(TypeError, match="sgs.variogram"):
        pst.upscale(logs.NetSand).sgs(variogram="spherical", seed=11)
    with pytest.raises(ValueError, match="sgs.seed"):
        pst.upscale(logs.NetSand).sgs(variogram=pst.Var("spherical", 1, 1, 1, 0))


def test_static_workflow_facade_validates_project_asset_names_loudly():
    with pytest.raises(ValueError, match="missing outline 'MissingEdge'"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="MissingEdge")
        )

    with pytest.raises(ValueError, match="missing horizon surface 'Missing top'"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(["Missing top", "Base reservoir"])
        )


def test_static_workflow_facade_horizons_accept_well_tie_mapping():
    grid = (
        pst.Grid.from_project(SyntheticProject())
        .geometry(cell=(10.0, 10.0), outline="ModelEdge")
        .horizons(
            [
                {"name": "Top reservoir", "well top": "well tops/Top reservoir"},
                "Base reservoir",
            ],
            well_tie={"influence_radius": 800},
        )
    )

    assert grid._horizons == ["Top reservoir", "Base reservoir"]
    assert grid._horizon_specs[0] == pst.HorizonSpec(
        "Top reservoir", well_top="well tops/Top reservoir"
    )
    assert grid._horizon_specs[1] == pst.HorizonSpec("Base reservoir")
    assert grid._well_tie == pst.WellTie(influence_radius=800.0)

    with pytest.raises(TypeError, match="tie_to_tops"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(["Top reservoir", "Base reservoir"], tie_to_tops=True)
        )
    with pytest.raises(TypeError, match="gridding"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(["Top reservoir", "Base reservoir"], gridding={})
        )
    assert not hasattr(pst, "Gridding")

    with pytest.raises(ValueError, match="unknown horizon field"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons([{"name": "Top reservoir", "top": "Top reservoir"}, "Base reservoir"])
        )
    with pytest.raises(ValueError, match="WellTie.influence_radius"):
        pst.WellTie(influence_radius=0)
    with pytest.raises(ValueError, match="WellTie.influence_radius"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(["Top reservoir", "Base reservoir"], well_tie={"influence_radius": 0})
        )
    with pytest.raises(ValueError, match="unknown well_tie field"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(["Top reservoir", "Base reservoir"], well_tie={"range": 800})
        )


def test_static_workflow_facade_horizons_accept_custom_names_and_surface_bindings():
    grid = (
        pst.Grid.from_project(SyntheticProject())
        .geometry(cell=(10.0, 10.0), outline="ModelEdge")
        .horizons(
            [
                {
                    "name": "Top reservoir",
                    "surface": "Top reservoir input surface",
                    "well top": "well tops/Top reservoir",
                    "zone": "Reservoir",
                },
                "Base reservoir",
                {
                    "name": "Custom model horizon name",
                    "surface": "input surface",
                },
            ],
            well_tie={"influence_radius": 800},
        )
        .zones({"Reservoir": ("Top reservoir", "Base reservoir")})
        .layers({"Reservoir": pst.Layering(n=2)})
    )

    assert grid._horizons == [
        "Top reservoir",
        "Base reservoir",
        "Custom model horizon name",
    ]
    assert grid._horizon_specs[0] == pst.HorizonSpec(
        "Top reservoir",
        surface="Top reservoir input surface",
        well_top="well tops/Top reservoir",
        zone="Reservoir",
    )
    assert grid._horizon_specs[2] == pst.HorizonSpec(
        "Custom model horizon name",
        surface="input surface",
    )
    assert grid._declared_cell_count() == 4

    with pytest.raises(ValueError, match="missing horizon surface 'Missing surface'"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(
                [
                    {"name": "Model top", "surface": "Missing surface"},
                    "Base reservoir",
                ],
                well_tie={"influence_radius": 800},
            )
    )


def test_static_workflow_facade_zones_accept_subzone_specs_and_insert_nested_horizons():
    grid = (
        pst.Grid.from_project(SyntheticProject())
        .geometry(cell=(10.0, 10.0), outline="ModelEdge")
        .horizons(
            [
                {
                    "name": "Top reservoir",
                    "surface": "Top reservoir input surface",
                    "well top": "well tops/Top reservoir",
                    "zone": "Reservoir",
                },
                "Base reservoir",
            ],
            well_tie={"influence_radius": 800},
        )
        .zones(
            [
                {
                    "zone": "Reservoir",
                    "sub-zones": [
                        {"name": "Upper Reservoir", "base": "Top Lower Reservoir"},
                        {"name": "Lower Reservoir", "top": "Top Lower Reservoir"},
                    ],
                }
            ]
        )
    )

    assert "Top Lower Reservoir" in grid._horizons
    assert grid._horizon_specs[-1] == pst.HorizonSpec("Top Lower Reservoir")
    assert grid._zones == {
        "Upper Reservoir": ("Top reservoir", "Top Lower Reservoir"),
        "Lower Reservoir": ("Top Lower Reservoir", "Base reservoir"),
    }

    with pytest.raises(ValueError, match="zone 'Missing' has no declared bounds"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(["Top reservoir", "Base reservoir"])
            .zones([{"zone": "Missing"}])
        )


def test_static_workflow_facade_horizon_zone_tag_accepts_inline_subzones():
    grid = (
        pst.Grid.from_project(SyntheticProject())
        .geometry(cell=(10.0, 10.0), outline="ModelEdge")
        .horizons(
            [
                {
                    "name": "Top reservoir",
                    "surface": "Top reservoir input surface",
                    "well top": "well tops/Top reservoir",
                    "zone": {
                        "name": "Reservoir",
                        "sub-zones": [
                            {"name": "Upper Reservoir", "base": "Top Lower Reservoir"},
                            {"name": "Lower Reservoir", "top": "Top Lower Reservoir"},
                        ],
                    },
                },
                "Base reservoir",
            ],
            well_tie={"influence_radius": 800},
        )
    )

    assert "Top Lower Reservoir" in grid._horizons
    assert grid._horizon_specs[0].zone == {
        "name": "Reservoir",
        "sub-zones": [
            {"name": "Upper Reservoir", "base": "Top Lower Reservoir"},
            {"name": "Lower Reservoir", "top": "Top Lower Reservoir"},
        ],
    }
    assert grid._zones == {
        "Upper Reservoir": ("Top reservoir", "Top Lower Reservoir"),
        "Lower Reservoir": ("Top Lower Reservoir", "Base reservoir"),
    }


def test_static_workflow_facade_subzones_accept_names_and_boundary_surfaces():
    grid = (
        pst.Grid.from_project(SyntheticProject())
        .geometry(cell=(10.0, 10.0), outline="ModelEdge")
        .horizons(
            [
                {
                    "name": "Top reservoir",
                    "surface": "Top reservoir input surface",
                    "well top": "well tops/Top reservoir",
                    "zone": {
                        "name": "Reservoir",
                        "sub-zones": [
                            "Upper Reservoir",
                            {"name": "Intra Shale", "surface": "Top Lower Reservoir"},
                            "Lower Reservoir",
                        ],
                    },
                },
                "Base reservoir",
            ],
            well_tie={"influence_radius": 800},
        )
    )

    assert grid._horizon_specs[-1] == pst.HorizonSpec(
        "Intra Shale",
        surface="Top Lower Reservoir",
    )
    assert grid._zones == {
        "Upper Reservoir": ("Top reservoir", "Intra Shale"),
        "Lower Reservoir": ("Intra Shale", "Base reservoir"),
    }

    with pytest.raises(ValueError, match="mixed sub-zones must end with a zone name"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(
                [
                    {
                        "name": "Top reservoir",
                        "zone": {
                            "name": "Reservoir",
                            "sub-zones": [
                                "Upper Reservoir",
                                {"surface": "Top Lower Reservoir"},
                            ],
                        },
                    },
                    "Base reservoir",
                ],
                well_tie={"influence_radius": 800},
            )
        )


def test_static_workflow_facade_subzones_accept_zone_types_and_well_top_boundaries():
    grid = (
        pst.Grid.from_project(SyntheticProject())
        .geometry(cell=(10.0, 10.0), outline="ModelEdge")
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
            ],
            well_tie={"influence_radius": 800},
        )
    )

    assert grid._horizon_specs[-1] == pst.HorizonSpec(
        "Intra Shale",
        well_top="Top Lower Reservoir",
    )
    assert grid._zones == {
        "Top Reservoir": ("Top reservoir", "Intra Shale"),
        "Lower Reservoir": ("Intra Shale", "Base reservoir"),
    }
    assert grid._zone_types == {
        "Top Reservoir": "constant",
        "Lower Reservoir": "isochore",
    }

    with pytest.raises(ValueError, match="zone type"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(
                [
                    {
                        "name": "Top reservoir",
                        "zone": {
                            "name": "Reservoir",
                            "sub-zones": [
                                {"zone": "Top Reservoir", "type": "pinch"},
                                {"name": "Lower Reservoir"},
                            ],
                        },
                    },
                    "Base reservoir",
                ],
                well_tie={"influence_radius": 800},
            )
        )


def test_static_workflow_facade_subzones_accept_unbound_model_surface_boundary():
    grid = (
        pst.Grid.from_project(SyntheticProject())
        .geometry(cell=(10.0, 10.0), outline="ModelEdge")
        .horizons(
            [
                {
                    "name": "Top reservoir",
                    "surface": "Top reservoir input surface",
                    "well top": "well tops/Top reservoir",
                    "zone": {
                        "name": "Reservoir",
                        "sub-zones": [
                            "Top Reservoir",
                            {"name": "Intra Shale", "surface": True},
                            "Lower Reservoir",
                        ],
                    },
                },
                "Base reservoir",
            ],
            well_tie={"influence_radius": 800},
        )
    )

    assert grid._horizon_specs[-1] == pst.HorizonSpec("Intra Shale")
    assert grid._zones == {
        "Top Reservoir": ("Top reservoir", "Intra Shale"),
        "Lower Reservoir": ("Intra Shale", "Base reservoir"),
    }

    with pytest.raises(ValueError, match="surface=True boundary requires a 'name'"):
        (
            pst.Grid.from_project(SyntheticProject())
            .geometry(cell=(10.0, 10.0), outline="ModelEdge")
            .horizons(
                [
                    {
                        "name": "Top reservoir",
                        "zone": {
                            "name": "Reservoir",
                            "sub-zones": [
                                "Top Reservoir",
                                {"surface": True},
                                "Lower Reservoir",
                            ],
                        },
                    },
                    "Base reservoir",
                ],
                well_tie={"influence_radius": 800},
            )
        )


def test_property_calc_is_atomic_on_formula_error():
    grid = declared_grid()
    grid.properties.set("X", [1.0, 2.0])

    with pytest.raises(ValueError, match="formula parameter"):
        grid.properties.calc(["Y = X + $missing"])

    assert "Y" not in grid.properties.names()
