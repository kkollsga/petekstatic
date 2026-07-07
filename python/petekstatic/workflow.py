"""Notebook-friendly static workflow facade.

This module is the first Python vertical slice of the petekStatic workflow API.
It stores the user's static-model declarations, validates named project assets,
keeps synthetic/smoke-test property arrays in memory, delegates formula blocks
to ``petektools.evaluate_formula``, and returns deterministic simple volumes.

It is intentionally not the production Rust corner-point grid builder yet. The
public object shape is kept close to the target API so the Rust lowering can
replace the internal execution path without changing notebooks.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from math import isfinite
from typing import Any


Number = int | float
Progress = bool | Callable[[dict[str, Any]], None] | None


@dataclass(frozen=True)
class Gridding:
    """Horizon-gridding options for the workflow declaration."""

    collapse_thin: bool = False
    min_thickness_m: float = 0.0

    def __post_init__(self) -> None:
        if not isfinite(self.min_thickness_m) or self.min_thickness_m < 0.0:
            raise ValueError("Gridding.min_thickness_m must be finite and >= 0")


@dataclass(frozen=True)
class Layering:
    """Per-zone layering declaration."""

    n: int
    method: str = "proportional"

    def __post_init__(self) -> None:
        if not isinstance(self.n, int) or self.n <= 0:
            raise ValueError("Layering.n must be a positive integer")
        if not self.method:
            raise ValueError("Layering.method must be non-empty")


@dataclass(frozen=True)
class Spherical:
    """Lightweight variogram spec placeholder for property declarations."""

    range_m: float
    sill: float = 1.0
    nugget: float = 0.0

    def __post_init__(self) -> None:
        if not isfinite(self.range_m) or self.range_m <= 0.0:
            raise ValueError("Spherical.range_m must be finite and > 0")
        if not isfinite(self.sill) or self.sill < 0.0:
            raise ValueError("Spherical.sill must be finite and >= 0")
        if not isfinite(self.nugget) or self.nugget < 0.0:
            raise ValueError("Spherical.nugget must be finite and >= 0")

    def to_var(self) -> "Var":
        """Return the anisotropic canonical form with isotropic ranges."""

        return Var(
            "spherical",
            major=self.range_m,
            minor=self.range_m,
            vertical=self.range_m,
            azimuth=0.0,
            sill=self.sill,
            nugget=self.nugget,
        )

    def as_dict(self) -> dict[str, Any]:
        return self.to_var().as_dict()

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()


@dataclass(frozen=True)
class Var:
    """Serializable anisotropic variogram declaration used by SGS recipes."""

    model: str
    major: float
    minor: float
    vertical: float
    azimuth: float
    sill: float | None = None
    nugget: float | None = None

    _MODELS = frozenset({"spherical", "exponential", "gaussian"})

    def __post_init__(self) -> None:
        model = str(self.model).strip().lower()
        if model not in self._MODELS:
            raise ValueError(
                "Var.model must be one of: " + ", ".join(sorted(self._MODELS))
            )
        object.__setattr__(self, "model", model)
        for field_name in ("major", "minor", "vertical"):
            value = float(getattr(self, field_name))
            if not isfinite(value) or value <= 0.0:
                raise ValueError(f"Var.{field_name} must be finite and > 0")
            object.__setattr__(self, field_name, value)
        azimuth = float(self.azimuth)
        if not isfinite(azimuth):
            raise ValueError("Var.azimuth must be finite")
        object.__setattr__(self, "azimuth", azimuth % 360.0)
        for field_name in ("sill", "nugget"):
            value = getattr(self, field_name)
            if value is None:
                continue
            numeric = float(value)
            if not isfinite(numeric) or numeric < 0.0:
                raise ValueError(f"Var.{field_name} must be finite and >= 0")
            object.__setattr__(self, field_name, numeric)

    def as_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "kind": "variogram",
            "model": self.model,
            "major": self.major,
            "minor": self.minor,
            "vertical": self.vertical,
            "azimuth": self.azimuth,
        }
        if self.sill is not None:
            out["sill"] = self.sill
        if self.nugget is not None:
            out["nugget"] = self.nugget
        return out

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()


@dataclass(frozen=True)
class DistributionSpec:
    """Serializable distribution control for property recipes."""

    kind: str
    params: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        kind = str(self.kind).strip().lower()
        if not kind:
            raise ValueError("DistributionSpec.kind must be non-empty")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "params", dict(self.params))

    def as_dict(self) -> dict[str, Any]:
        out = {"kind": self.kind}
        out.update(
            {
                key: _serialize_recipe_value(value)
                for key, value in self.params.items()
            }
        )
        return out

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()


class _Distributions:
    """Namespace for serializable distribution specs."""

    def from_logs(self) -> DistributionSpec:
        return DistributionSpec("from_logs")

    def from_trend(self) -> DistributionSpec:
        return DistributionSpec("from_trend")

    def normal(
        self,
        *,
        mean: Number = 0.0,
        std: Number = 1.0,
        **params: Any,
    ) -> DistributionSpec:
        mean_f = float(mean)
        std_f = float(std)
        if not isfinite(mean_f):
            raise ValueError("normal distribution mean must be finite")
        if not isfinite(std_f) or std_f <= 0.0:
            raise ValueError("normal distribution std must be finite and > 0")
        clean: dict[str, Any] = {"mean": mean_f, "std": std_f}
        for key, value in params.items():
            if isinstance(value, (int, float)):
                numeric = float(value)
                if not isfinite(numeric):
                    raise ValueError(
                        f"normal distribution parameter '{key}' must be finite"
                    )
                clean[key] = numeric
            else:
                clean[key] = value
        return DistributionSpec("normal", clean)


distributions = _Distributions()


@dataclass(frozen=True)
class CoKriging:
    """Serializable collocated cokriging trend control."""

    trend: Any
    rho: float

    def __post_init__(self) -> None:
        if self.trend is None:
            raise ValueError("CoKriging.trend is required")
        rho = float(self.rho)
        if not isfinite(rho) or rho < -1.0 or rho > 1.0:
            raise ValueError("CoKriging.rho must be finite and between -1 and 1")
        object.__setattr__(self, "rho", rho)

    def as_dict(self) -> dict[str, Any]:
        return {
            "kind": "cokriging",
            "trend": _serialize_recipe_value(self.trend),
            "rho": self.rho,
        }

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()


@dataclass(frozen=True)
class UpscaleRecipeBuilder:
    """Builder returned by ``upscale(...)`` before a fill method is selected."""

    source: Any
    method: str = "arithmetic"

    def __post_init__(self) -> None:
        if self.source is None:
            raise ValueError("upscale.source is required")
        method = str(self.method).strip().lower()
        if not method:
            raise ValueError("upscale.method must be non-empty")
        object.__setattr__(self, "method", method)

    def sgs(
        self,
        variogram: Var | Spherical | None = None,
        *,
        distribution: DistributionSpec | None = None,
        seed: int | None = None,
        cokriging: CoKriging | None = None,
        **options: Any,
    ) -> "SgsRecipe":
        return SgsRecipe(
            source=self.source,
            method=self.method,
            variogram=variogram,
            distribution=distribution,
            seed=seed,
            cokriging=cokriging,
            options=dict(options),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "kind": "upscale",
            "source": _serialize_recipe_value(self.source),
            "method": self.method,
        }

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()


@dataclass(frozen=True)
class SgsRecipe:
    """Serializable upscale + SGS property declaration."""

    source: Any
    method: str
    variogram: Var | Spherical | None
    distribution: DistributionSpec | None
    seed: int | None
    cokriging: CoKriging | None = None
    options: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.source is None:
            raise ValueError("SgsRecipe.source is required")
        if isinstance(self.variogram, Spherical):
            object.__setattr__(self, "variogram", self.variogram.to_var())
        if not isinstance(self.variogram, Var):
            raise TypeError("sgs.variogram must be a Var")
        if self.distribution is not None and not isinstance(
            self.distribution, DistributionSpec
        ):
            raise TypeError("sgs.distribution must be a DistributionSpec")
        if self.seed is None:
            raise ValueError("sgs.seed is required")
        if not isinstance(self.seed, int):
            raise TypeError("sgs.seed must be an integer")
        if self.cokriging is not None and not isinstance(self.cokriging, CoKriging):
            raise TypeError("sgs.cokriging must be a CoKriging")
        method = str(self.method).strip().lower()
        if not method:
            raise ValueError("SgsRecipe.method must be non-empty")
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "options", dict(self.options))

    def as_dict(self) -> dict[str, Any]:
        return {
            "kind": "sgs",
            "upscale": {
                "source": _serialize_recipe_value(self.source),
                "method": self.method,
            },
            "variogram": self.variogram.as_dict(),
            "distribution": (
                None if self.distribution is None else self.distribution.as_dict()
            ),
            "seed": self.seed,
            "cokriging": None if self.cokriging is None else self.cokriging.as_dict(),
            "options": {
                key: _serialize_recipe_value(value)
                for key, value in sorted(self.options.items())
            },
        }

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()

    def lower(
        self,
        property_name: str,
        *,
        project: Any | None = None,
    ) -> "PropertyPipelineSpec":
        """Lower this recipe to the production ``PropertyPipeline`` shape.

        The current Python extension does not expose the Rust ``PropertyPipeline``
        class, so this returns the exact pipeline inputs petekStatic can derive
        today and leaves execution behind an explicit unsupported boundary.
        """

        return _lower_sgs_recipe(PropertyStore.key(property_name), self, project=project)


def upscale(source: Any, method: str = "arithmetic") -> UpscaleRecipeBuilder:
    """Declare a log/source upscaling recipe without executing it."""

    return UpscaleRecipeBuilder(source=source, method=method)


@dataclass(frozen=True)
class _GeometrySpec:
    cell: tuple[float, float]
    orient: float
    outline: str | None


@dataclass(frozen=True)
class _PropertyStep:
    kind: str
    args: Mapping[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "args": _serialize_recipe_value(self.args),
        }

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()


@dataclass(frozen=True)
class WellLogSpec:
    """Positioned well-log input for Rust ``WellLog::new``."""

    x: float
    y: float
    samples: Sequence[tuple[float, float]]

    def __post_init__(self) -> None:
        x = float(self.x)
        y = float(self.y)
        if not isfinite(x) or not isfinite(y):
            raise ValueError("WellLogSpec x/y must be finite")
        clean: list[tuple[float, float]] = []
        for depth, value in self.samples:
            depth_f = float(depth)
            value_f = float(value)
            if not isfinite(depth_f) or not isfinite(value_f):
                raise ValueError("WellLogSpec samples must be finite depth/value pairs")
            clean.append((depth_f, value_f))
        if not clean:
            raise ValueError("WellLogSpec requires at least one sample")
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "y", y)
        object.__setattr__(self, "samples", tuple(clean))

    def as_dict(self) -> dict[str, Any]:
        return {
            "x": self.x,
            "y": self.y,
            "samples": [[depth, value] for depth, value in self.samples],
        }

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()


@dataclass(frozen=True)
class PropertyPipelineSpec:
    """Python lowering of a recipe into the Rust ``PropertyPipeline`` path."""

    property: str
    source: Mapping[str, Any]
    method: str
    variogram: Var
    distribution: DistributionSpec
    seed: int
    cokriging: CoKriging | None = None
    options: Mapping[str, Any] = field(default_factory=dict)
    well_logs: Sequence[WellLogSpec] | None = None

    def __post_init__(self) -> None:
        key = PropertyStore.key(self.property)
        if not key:
            raise ValueError("PropertyPipelineSpec.property must be non-empty")
        method = str(self.method).strip().lower()
        if method not in {"arithmetic", "harmonic", "geometric"}:
            raise NotImplementedError(
                f"PropertyPipeline lowering does not support upscale method '{method}'"
            )
        if not isinstance(self.variogram, Var):
            raise TypeError("PropertyPipelineSpec.variogram must be a Var")
        if not isinstance(self.distribution, DistributionSpec):
            raise TypeError(
                "PropertyPipelineSpec.distribution must be a DistributionSpec"
            )
        if not isinstance(self.seed, int):
            raise TypeError("PropertyPipelineSpec.seed must be an integer")
        object.__setattr__(self, "property", key)
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "source", dict(self.source))
        object.__setattr__(self, "options", dict(self.options))
        if self.well_logs is not None:
            object.__setattr__(self, "well_logs", tuple(self.well_logs))

    def as_dict(self) -> dict[str, Any]:
        return {
            "kind": "property_pipeline",
            "property": self.property,
            "upscale": {
                "rust_type": "PropertyPipeline::upscale",
                "source": _serialize_recipe_value(self.source),
                "well_logs": (
                    None
                    if self.well_logs is None
                    else [well.as_dict() for well in self.well_logs]
                ),
                "method": self.method,
            },
            "propagate": {
                "rust_type": "Gaussian",
                "variogram": self.variogram.as_dict(),
                "distribution": self.distribution.as_dict(),
                "seed": self.seed,
                "cokriging": (
                    None if self.cokriging is None else self.cokriging.as_dict()
                ),
                "options": {
                    key: _serialize_recipe_value(value)
                    for key, value in sorted(self.options.items())
                },
            },
        }

    def to_dict(self) -> dict[str, Any]:
        return self.as_dict()

    def to_petektools_variogram(self) -> Any:
        """Build the petekTools variogram object used by the lower SGS layer."""

        try:
            import petektools as pt
        except Exception as exc:  # pragma: no cover - import environment dependent
            raise RuntimeError("petektools is required for variogram lowering") from exc
        cls = getattr(pt, "AnisotropicVariogram", None)
        if cls is None:
            raise NotImplementedError(
                "petektools.AnisotropicVariogram is required for anisotropic SGS lowering"
            )
        return cls(
            self.variogram.model,
            major=self.variogram.major,
            minor=self.variogram.minor,
            vertical=self.variogram.vertical,
            azimuth=self.variogram.azimuth,
            sill=1.0 if self.variogram.sill is None else self.variogram.sill,
            nugget=0.0 if self.variogram.nugget is None else self.variogram.nugget,
        )

    def require_executable(self) -> None:
        """Raise unless the lowered spec has enough support to execute."""

        if self.well_logs is None:
            raise NotImplementedError(
                "PropertyPipeline execution needs positioned WellLog inputs; "
                "the lazy log source was serialized but not resolved by the project"
            )
        if self.cokriging is not None:
            raise NotImplementedError(
                "PropertyPipeline execution does not yet expose Rust collocated "
                "cokriging/trend binding"
            )
        if self.distribution.kind != "from_logs":
            raise NotImplementedError(
                "PropertyPipeline execution currently supports distribution=from_logs; "
                f"got {self.distribution.kind!r}"
            )
        self._require_rust_variogram_supported()
        try:
            from . import _petekstatic
        except Exception as exc:  # pragma: no cover - import environment dependent
            raise RuntimeError("petekstatic._petekstatic is required for execution") from exc
        if not hasattr(_petekstatic, "PropertyPipeline"):
            raise NotImplementedError(
                "petekstatic._petekstatic.PropertyPipeline is not exposed to Python yet"
            )

    def execute(self) -> Any:
        """Build and return the Rust-backed ``PropertyPipeline`` handle.

        The returned object proves the executable construction boundary: Python
        recipe inputs have become the production Rust ``PropertyPipeline``. Applying
        that pipeline to an arbitrary mutable grid remains explicit until a
        production grid binding exists.
        """

        self.require_executable()
        from . import _petekstatic

        options = self._execution_options()
        wells = [
            _petekstatic.WellLog(well.x, well.y, list(well.samples))
            for well in self.well_logs or ()
        ]
        return _petekstatic.PropertyPipeline(
            self.property,
            wells,
            self.method,
            self.variogram.model,
            self.variogram.major,
            self.seed,
            sill=1.0 if self.variogram.sill is None else self.variogram.sill,
            nugget=0.0 if self.variogram.nugget is None else self.variogram.nugget,
            minor_m=self.variogram.minor,
            vertical_m=self.variogram.vertical,
            azimuth=self.variogram.azimuth,
            **options,
        )

    def _require_rust_variogram_supported(self) -> None:
        if self.variogram.model not in {"spherical", "exponential", "gaussian"}:
            raise NotImplementedError(
                "PropertyPipeline execution supports spherical, exponential, "
                f"and gaussian variograms; got {self.variogram.model!r}"
            )

    def _execution_options(self) -> dict[str, Any]:
        allowed = {"allow_mean_fill", "max_neighbours", "radius_m", "unbounded_search"}
        unknown = sorted(set(self.options) - allowed)
        if unknown:
            raise NotImplementedError(
                "PropertyPipeline execution does not support option(s): "
                + ", ".join(unknown)
            )
        out: dict[str, Any] = {
            "propagate": True,
            "allow_mean_fill": bool(self.options.get("allow_mean_fill", False)),
            "unbounded_search": bool(self.options.get("unbounded_search", False)),
        }
        if "max_neighbours" in self.options or "radius_m" in self.options:
            if "max_neighbours" not in self.options or "radius_m" not in self.options:
                raise ValueError("max_neighbours and radius_m must be set together")
            max_neighbours = self.options["max_neighbours"]
            if not isinstance(max_neighbours, int) or max_neighbours <= 0:
                raise ValueError("max_neighbours must be a positive integer")
            radius_m = float(self.options["radius_m"])
            if not isfinite(radius_m) or radius_m <= 0.0:
                raise ValueError("radius_m must be finite and > 0")
            out["max_neighbours"] = max_neighbours
            out["radius_m"] = radius_m
        return out


class Grid:
    """Chainable static workflow declaration.

    Use ``Grid.from_project(project)`` and then declare geometry, horizons,
    zones, layers, properties, and volume cases. The current implementation is a
    Python spec/execution layer for synthetic tests and early notebooks.
    """

    def __init__(self, project: Any):
        self.project = project
        self._geometry: _GeometrySpec | None = None
        self._horizons: list[str] = []
        self._tie_to_tops = False
        self._gridding = Gridding()
        self._zones: dict[str, tuple[str, str]] = {}
        self._layers: dict[str, Layering] = {}
        self._property_steps: dict[str, list[_PropertyStep]] = {}
        self._property_pipelines: dict[str, PropertyPipelineSpec] = {}
        self.properties = PropertyStore(self)

    @classmethod
    def from_project(cls, project: Any) -> "Grid":
        if project is None:
            raise ValueError("Grid.from_project requires a project")
        return cls(project)

    def geometry(
        self,
        *,
        cell: tuple[Number, Number],
        orient: Number = 0.0,
        outline: str | None = None,
    ) -> "Grid":
        dx, dy = (float(cell[0]), float(cell[1]))
        orient_f = float(orient)
        if not isfinite(dx) or not isfinite(dy) or dx <= 0.0 or dy <= 0.0:
            raise ValueError("geometry.cell must contain finite positive dx/dy")
        if not isfinite(orient_f):
            raise ValueError("geometry.orient must be finite")
        if outline is not None:
            _require_project_asset(self.project, "outline", outline, ("polygons", "outlines"))
        self._geometry = _GeometrySpec(cell=(dx, dy), orient=orient_f, outline=outline)
        return self

    def horizons(
        self,
        names: Sequence[str],
        *,
        tie_to_tops: bool = False,
        gridding: Gridding | None = None,
    ) -> "Grid":
        if len(names) < 2:
            raise ValueError("horizons requires at least a top and a base horizon")
        clean = [str(name) for name in names]
        for name in clean:
            try:
                _require_project_asset(self.project, "horizon", name, ("surfaces", "horizons"))
            except ValueError:
                if not tie_to_tops:
                    raise
                _require_project_asset(self.project, "horizon", name, ("tops", "well_tops"))
        if gridding is not None and not isinstance(gridding, Gridding):
            raise TypeError("horizons.gridding must be a Gridding instance")
        self._horizons = clean
        self._tie_to_tops = bool(tie_to_tops)
        self._gridding = gridding or Gridding()
        return self

    def zones(self, zones: Mapping[str, tuple[str, str]]) -> "Grid":
        if not zones:
            raise ValueError("zones requires at least one named zone")
        horizon_set = set(self._horizons)
        out: dict[str, tuple[str, str]] = {}
        for zone, pair in zones.items():
            if len(pair) != 2:
                raise ValueError(f"zone '{zone}' must be bounded by two horizons")
            top, base = str(pair[0]), str(pair[1])
            missing = [name for name in (top, base) if name not in horizon_set]
            if missing:
                raise ValueError(
                    f"zone '{zone}' references undeclared horizon(s): {', '.join(missing)}"
                )
            out[str(zone)] = (top, base)
        self._zones = out
        return self

    def layers(self, layers: Mapping[str, Layering | int]) -> "Grid":
        if not self._zones:
            raise ValueError("declare zones before layers")
        out: dict[str, Layering] = {}
        for zone, layering in layers.items():
            if zone not in self._zones:
                raise ValueError(f"layers references unknown zone '{zone}'")
            out[str(zone)] = layering if isinstance(layering, Layering) else Layering(int(layering))
        missing = sorted(set(self._zones) - set(out))
        if missing:
            raise ValueError(f"missing layering for zone(s): {', '.join(missing)}")
        self._layers = out
        return self

    def volumes(
        self,
        *,
        ntg: str,
        por: str,
        sw: str,
        fluid: str = "oil",
        fvf: Number | None = None,
        bo: Number | None = None,
        contacts: Mapping[str, Any] | None = None,
    ) -> "VolumeCase":
        return VolumeCase(
            self,
            ntg=PropertyStore.key(ntg),
            por=PropertyStore.key(por),
            sw=PropertyStore.key(sw),
            fluid=str(fluid),
            fvf=1.0 if fvf is None and bo is None else float(bo if bo is not None else fvf),
            contacts=dict(contacts or {}),
        )

    def _cell_bulk_volumes(self, n_cells: int) -> list[float]:
        if self._geometry is None:
            raise ValueError("declare geometry before running volumes")
        if not self._zones:
            raise ValueError("declare zones before running volumes")
        if not self._layers:
            raise ValueError("declare layers before running volumes")

        area = self._geometry.cell[0] * self._geometry.cell[1]
        vectors: list[float] = []
        for zone, (top_name, base_name) in self._zones.items():
            layers = self._layers[zone].n
            thickness = _zone_thickness(self.project, top_name, base_name)
            if thickness:
                for dz in thickness:
                    per_layer = area * max(dz, 0.0) / layers
                    vectors.extend([per_layer] * layers)
            else:
                vectors.extend([area] * layers)

        if len(vectors) == n_cells:
            return vectors
        if len(vectors) == 1:
            return vectors * n_cells
        if not vectors:
            return [area] * n_cells
        if n_cells % len(vectors) == 0:
            repeated: list[float] = []
            for value in vectors:
                repeated.extend([value] * (n_cells // len(vectors)))
            return repeated
        if len(vectors) % n_cells == 0:
            stride = len(vectors) // n_cells
            return [sum(vectors[i : i + stride]) for i in range(0, len(vectors), stride)]
        raise ValueError(
            "property array length does not match declared grid cells "
            f"({n_cells} values vs {len(vectors)} inferred cells)"
        )

    def _declared_cell_count(self) -> int | None:
        if not self._zones or not self._layers:
            return None
        total = 0
        for zone, (top_name, base_name) in self._zones.items():
            layers = self._layers[zone].n
            thickness = _zone_thickness(self.project, top_name, base_name)
            total += (len(thickness) if thickness else 1) * layers
        return total or None


class PropertyStore:
    """In-memory property arrays and declarations for a ``Grid``."""

    _ALIASES = {"ntg": "NTG", "por": "POR", "poro": "POR", "sw": "SW"}

    def __init__(self, grid: Grid):
        object.__setattr__(self, "_grid", grid)
        object.__setattr__(self, "_arrays", {})

    def __getattr__(self, name: str) -> "PropertyHandle":
        if name.startswith("_"):
            raise AttributeError(name)
        return self[self.key(name)]

    def __setattr__(
        self,
        name: str,
        value: Number | Iterable[Number] | "PropertyExpression" | SgsRecipe,
    ) -> None:
        if name.startswith("_"):
            object.__setattr__(self, name, value)
            return
        self.set(name, value)

    def __getitem__(self, name: str) -> "PropertyHandle":
        return PropertyHandle(self, self.key(name))

    def __contains__(self, name: object) -> bool:
        key = self.key(str(name))
        return key in self._arrays or key in self._grid._property_steps

    @classmethod
    def key(cls, name: str) -> str:
        return cls._ALIASES.get(str(name), str(name))

    def names(self) -> list[str]:
        return sorted(set(self._arrays) | set(self._grid._property_steps))

    def set(
        self,
        name: str,
        values: Number | Iterable[Number] | "PropertyExpression" | SgsRecipe,
    ) -> "PropertyStore":
        key = self.key(name)
        if isinstance(values, SgsRecipe):
            pipeline = values.lower(key, project=self._grid.project)
            recipe = values.as_dict()
            self._arrays.pop(key, None)
            self._grid._property_steps[key] = [_PropertyStep(kind="recipe", args=recipe)]
            self._grid._property_pipelines[key] = pipeline
            return self
        if isinstance(values, PropertyExpression):
            self.calc({key: values.formula})
            self._grid._property_steps.pop(key, None)
            self._grid._property_pipelines.pop(key, None)
            return self
        target_len = _unique_length(self._arrays.values()) or self._grid._declared_cell_count()
        array = _broadcast_if_scalar(_float_vector(values, f"property '{key}'"), target_len)
        self._arrays[key] = array
        self._grid._property_steps.pop(key, None)
        self._grid._property_pipelines.pop(key, None)
        return self

    def declarations(
        self,
        name: str | None = None,
    ) -> dict[str, list[dict[str, Any]]] | list[dict[str, Any]]:
        """Return serialization-ready property workflow declarations."""

        if name is not None:
            key = self.key(name)
            return [step.as_dict() for step in self._grid._property_steps.get(key, [])]
        return {
            key: [step.as_dict() for step in steps]
            for key, steps in sorted(self._grid._property_steps.items())
        }

    def pipelines(
        self,
        name: str | None = None,
    ) -> dict[str, dict[str, Any]] | dict[str, Any] | None:
        """Return lowered production ``PropertyPipeline`` specs."""

        if name is not None:
            key = self.key(name)
            pipeline = self._grid._property_pipelines.get(key)
            return None if pipeline is None else pipeline.as_dict()
        return {
            key: pipeline.as_dict()
            for key, pipeline in sorted(self._grid._property_pipelines.items())
        }

    def pipeline_spec(self, name: str) -> PropertyPipelineSpec:
        key = self.key(name)
        try:
            return self._grid._property_pipelines[key]
        except KeyError as exc:
            raise ValueError(f"missing lowered property pipeline '{key}'") from exc

    def execute_pipeline(self, name: str) -> Any:
        """Execute a lowered recipe once Rust pipeline bindings are available."""

        return self.pipeline_spec(name).execute()

    def values(self, name: str) -> list[float]:
        key = self.key(name)
        try:
            return list(self._arrays[key])
        except KeyError as exc:
            raise ValueError(f"missing grid property '{key}'") from exc

    def calc(
        self,
        formulas: Sequence[str | Mapping[str, str]] | Mapping[str, str],
        *,
        params: Mapping[str, Number] | None = None,
    ) -> dict[str, list[float]]:
        if not formulas:
            raise ValueError("calc requires at least one formula")
        evaluator = _petektools_formula_evaluator()
        lines = _formula_lines(formulas)
        properties = dict(self._arrays)
        for alias, canonical in self._ALIASES.items():
            if canonical in self._arrays:
                properties.setdefault(alias, self._arrays[canonical])
        # petekTools validates missing names, parameter binding, cycles and shape
        # mismatch before returning; update only after the whole block succeeds.
        out = evaluator(lines, properties, dict(params or {}))
        target_len = _unique_length(self._arrays.values())
        clean = {
            self.key(str(name)): _broadcast_if_scalar(
                _float_vector(values, f"calculated property '{name}'"),
                target_len,
            )
            for name, values in out.items()
        }
        self._arrays.update(clean)
        return {name: list(values) for name, values in clean.items()}

    def _length(self, names: Sequence[str]) -> int:
        lengths = {name: len(self.values(name)) for name in names}
        unique = set(lengths.values())
        if len(unique) != 1:
            raise ValueError(f"volume input property lengths differ: {lengths}")
        return unique.pop()


@dataclass(frozen=True)
class PropertyHandle:
    """Named property handle used for declarations and direct array access."""

    store: PropertyStore
    name: str

    @property
    def values(self) -> list[float]:
        return self.store.values(self.name)

    def set(
        self,
        values: Number | Iterable[Number] | "PropertyExpression" | SgsRecipe,
    ) -> "PropertyHandle":
        self.store.set(self.name, values)
        return self

    def _expr(self) -> "PropertyExpression":
        return PropertyExpression(self.name)

    def __add__(self, other: Any) -> "PropertyExpression":
        return self._expr() + other

    def __radd__(self, other: Any) -> "PropertyExpression":
        return _property_expr(other) + self

    def __sub__(self, other: Any) -> "PropertyExpression":
        return self._expr() - other

    def __rsub__(self, other: Any) -> "PropertyExpression":
        return _property_expr(other) - self

    def __mul__(self, other: Any) -> "PropertyExpression":
        return self._expr() * other

    def __rmul__(self, other: Any) -> "PropertyExpression":
        return _property_expr(other) * self

    def __truediv__(self, other: Any) -> "PropertyExpression":
        return self._expr() / other

    def __rtruediv__(self, other: Any) -> "PropertyExpression":
        return _property_expr(other) / self

    def upscale(self, **kwargs: Any) -> "PropertyHandle":
        self._record("upscale", kwargs)
        return self

    def propagate(self, **kwargs: Any) -> "PropertyHandle":
        self._record("propagate", kwargs)
        return self

    def _record(self, kind: str, kwargs: Mapping[str, Any]) -> None:
        steps = self.store._grid._property_steps.setdefault(self.name, [])
        steps.append(_PropertyStep(kind=kind, args=dict(kwargs)))


@dataclass(frozen=True)
class PropertyExpression:
    """Formula expression used by assignment sugar such as ``p.por = p.ntg * 0.85``."""

    formula: str

    def _binary(self, op: str, other: Any) -> "PropertyExpression":
        return PropertyExpression(f"({self.formula} {op} {_property_expr(other).formula})")

    def __add__(self, other: Any) -> "PropertyExpression":
        return self._binary("+", other)

    def __radd__(self, other: Any) -> "PropertyExpression":
        return _property_expr(other)._binary("+", self)

    def __sub__(self, other: Any) -> "PropertyExpression":
        return self._binary("-", other)

    def __rsub__(self, other: Any) -> "PropertyExpression":
        return _property_expr(other)._binary("-", self)

    def __mul__(self, other: Any) -> "PropertyExpression":
        return self._binary("*", other)

    def __rmul__(self, other: Any) -> "PropertyExpression":
        return _property_expr(other)._binary("*", self)

    def __truediv__(self, other: Any) -> "PropertyExpression":
        return self._binary("/", other)

    def __rtruediv__(self, other: Any) -> "PropertyExpression":
        return _property_expr(other)._binary("/", self)


@dataclass(frozen=True)
class VolumeCase:
    """Deferred deterministic volume calculation."""

    grid: Grid
    ntg: str
    por: str
    sw: str
    fluid: str
    fvf: float
    contacts: Mapping[str, Any] = field(default_factory=dict)

    def run(
        self,
        *,
        params: Mapping[str, Number] | None = None,
        progress: Progress = False,
        workers: str | int | None = "auto",
    ) -> "VolumeResult":
        if params:
            for name, value in params.items():
                if not isfinite(float(value)):
                    raise ValueError(f"volume parameter '{name}' must be finite")
        if not isfinite(self.fvf) or self.fvf <= 0.0:
            raise ValueError("fvf/bo must be finite and > 0")

        events: list[dict[str, Any]] = []
        emit = _progress_emitter(progress, events)
        emit("structure", "using declarative Python workflow spec")
        n_cells = self.grid.properties._length([self.ntg, self.por, self.sw])
        bulk = self.grid._cell_bulk_volumes(n_cells)
        ntg = self.grid.properties.values(self.ntg)
        por = self.grid.properties.values(self.por)
        sw = self.grid.properties.values(self.sw)

        emit("volumes", f"computing {n_cells} cells")
        grv = sum(bulk)
        hcpv = 0.0
        for i, volume in enumerate(bulk):
            hcpv += volume * ntg[i] * por[i] * (1.0 - sw[i])
        in_place = hcpv / self.fvf
        metric = "ooip_sm3" if self.fluid.lower() == "oil" else "ogip_sm3"
        total = {
            "fluid": self.fluid,
            "cells": n_cells,
            "grv_m3": grv,
            "hcpv_m3": hcpv,
            "in_place_sm3": in_place,
            metric: in_place,
            "fvf": self.fvf,
        }
        zones = {"total": dict(total)}
        emit("complete", "volume case complete")
        return VolumeResult(total=total, zones=zones, progress_events=events, workers=workers)


@dataclass(frozen=True)
class VolumeResult:
    """Volume result returned by ``VolumeCase.run``."""

    total: Mapping[str, Any]
    zones: Mapping[str, Mapping[str, Any]]
    progress_events: Sequence[Mapping[str, Any]]
    workers: str | int | None = "auto"

    def summary(self) -> dict[str, Any]:
        return dict(self.total)

    def by_zone(self) -> dict[str, dict[str, Any]]:
        return {name: dict(values) for name, values in self.zones.items()}


def _float_vector(values: Number | Iterable[Number], label: str) -> list[float]:
    if isinstance(values, (int, float)):
        values = [values]
    try:
        out = [float(value) for value in values]
    except TypeError as exc:
        raise TypeError(f"{label} must be an iterable of numbers") from exc
    if not out:
        raise ValueError(f"{label} must not be empty")
    return out


def _unique_length(vectors: Iterable[Sequence[float]]) -> int | None:
    lengths = {len(values) for values in vectors}
    if not lengths:
        return None
    if len(lengths) != 1:
        return None
    return lengths.pop()


def _broadcast_if_scalar(values: list[float], target_len: int | None) -> list[float]:
    if target_len is not None and target_len > 1 and len(values) == 1:
        return values * target_len
    return values


def _property_expr(value: Any) -> PropertyExpression:
    if isinstance(value, PropertyExpression):
        return value
    if isinstance(value, PropertyHandle):
        return value._expr()
    if isinstance(value, (int, float)):
        numeric = float(value)
        if not isfinite(numeric):
            raise ValueError("property expression constants must be finite")
        return PropertyExpression(repr(numeric))
    raise TypeError(
        "property expressions support property handles, numbers, and other expressions"
    )


def _serialize_recipe_value(value: Any) -> Any:
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, PropertyHandle):
        return {"kind": "property", "name": value.name}
    if isinstance(value, PropertyExpression):
        return {"kind": "expression", "formula": value.formula}
    if isinstance(value, Mapping):
        return {str(key): _serialize_recipe_value(item) for key, item in value.items()}
    if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        return [_serialize_recipe_value(item) for item in value]
    as_dict = getattr(value, "as_dict", None)
    if callable(as_dict):
        return _serialize_recipe_value(as_dict())
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        return _serialize_recipe_value(to_dict())
    return {"repr": repr(value)}


def _lower_sgs_recipe(
    property_name: str,
    recipe: SgsRecipe,
    *,
    project: Any | None,
) -> PropertyPipelineSpec:
    source = _log_source_dict(recipe.source)
    distribution = recipe.distribution or _default_distribution(source, recipe.cokriging)
    well_logs = _resolve_positioned_well_logs(recipe.source, project, source)
    return PropertyPipelineSpec(
        property=property_name,
        source=source,
        method=recipe.method,
        variogram=recipe.variogram,
        distribution=distribution,
        seed=recipe.seed,
        cokriging=recipe.cokriging,
        options=recipe.options,
        well_logs=well_logs,
    )


def _log_source_dict(source: Any) -> dict[str, Any]:
    data = _serialize_recipe_value(source)
    if not isinstance(data, Mapping):
        raise TypeError("upscale source must serialize to a mapping")
    out = dict(data)
    kind = out.get("kind")
    if kind in {"log", "log_channel"}:
        return out
    raise NotImplementedError(
        "PropertyPipeline lowering currently supports log-channel upscale sources"
    )


def _default_distribution(
    source: Mapping[str, Any],
    cokriging: CoKriging | None,
) -> DistributionSpec:
    if source.get("kind") in {"log", "log_channel"}:
        return distributions.from_logs()
    if cokriging is not None:
        return distributions.from_trend()
    raise NotImplementedError(
        "sgs.distribution is required when it cannot be inferred from logs or trend"
    )


def _resolve_positioned_well_logs(
    source: Any,
    project: Any | None,
    source_dict: Mapping[str, Any],
) -> tuple[WellLogSpec, ...] | None:
    for owner, method_names, args in (
        (source, ("to_well_logs", "resolve_well_logs"), (project,)),
        (project, ("resolve_log_expression", "resolve_well_logs"), (source,)),
        (project, ("resolve_log_source",), (source_dict,)),
    ):
        if owner is None:
            continue
        for method_name in method_names:
            method = getattr(owner, method_name, None)
            if not callable(method):
                continue
            result = _call_resolver(method, args)
            if result is not None:
                return tuple(_coerce_well_log_spec(item) for item in result)
    return None


def _call_resolver(method: Callable[..., Any], args: tuple[Any, ...]) -> Any:
    clean_args = tuple(arg for arg in args if arg is not None)
    try:
        return method(*clean_args)
    except TypeError:
        return method()


def _coerce_well_log_spec(value: Any) -> WellLogSpec:
    if isinstance(value, WellLogSpec):
        return value
    if isinstance(value, Mapping):
        return WellLogSpec(
            x=value["x"],
            y=value["y"],
            samples=_coerce_samples(value["samples"]),
        )
    x = getattr(value, "x", None)
    y = getattr(value, "y", None)
    samples = getattr(value, "samples", None)
    if x is None or y is None or samples is None:
        raise TypeError("positioned well logs need x, y and samples")
    return WellLogSpec(x=x, y=y, samples=_coerce_samples(samples))


def _coerce_samples(samples: Any) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    for sample in samples:
        if isinstance(sample, Mapping):
            depth = sample.get("depth_m", sample.get("tvd", sample.get("depth")))
            value = sample.get("value")
            out.append((depth, value))
        else:
            depth, value = sample
            out.append((depth, value))
    return out


def _petektools_formula_evaluator():
    try:
        import petektools as pt
    except Exception as exc:  # pragma: no cover - import environment dependent
        raise RuntimeError("petektools is required for property formula evaluation") from exc
    evaluator = getattr(pt, "evaluate_formula", None)
    if evaluator is None:
        raise RuntimeError(
            "petektools.evaluate_formula is required for property formula evaluation"
        )
    return evaluator


def _formula_lines(
    formulas: Sequence[str | Mapping[str, str]] | Mapping[str, str],
) -> list[str]:
    if isinstance(formulas, Mapping):
        return [f"{lhs} = {rhs}" for lhs, rhs in formulas.items()]
    lines: list[str] = []
    for item in formulas:
        if isinstance(item, Mapping):
            lines.extend(f"{lhs} = {rhs}" for lhs, rhs in item.items())
        else:
            lines.append(str(item))
    return lines


def _progress_emitter(progress: Progress, events: list[dict[str, Any]]):
    def emit(stage: str, message: str) -> None:
        event = {"stage": stage, "message": message}
        events.append(event)
        if progress is True:
            print(f"petekstatic: {stage}: {message}")
        elif callable(progress):
            progress(dict(event))

    return emit


def _project_inventory(project: Any) -> Mapping[str, Any]:
    inv = getattr(project, "inventory", None)
    if callable(inv):
        try:
            data = inv()
            if isinstance(data, Mapping):
                return data
        except Exception:
            return {}
    return {}


def _asset_names(project: Any, kinds: Sequence[str]) -> set[str]:
    names: set[str] = set()
    inv = _project_inventory(project)
    for kind in kinds:
        value = inv.get(kind)
        if isinstance(value, str):
            names.add(value)
        elif isinstance(value, Iterable):
            names.update(str(item) for item in value)

        attr = getattr(project, kind, None)
        if attr is None:
            continue
        names_method = getattr(attr, "names", None)
        if callable(names_method):
            names.update(str(item) for item in names_method())
        elif isinstance(attr, Mapping):
            names.update(str(item) for item in attr.keys())
        elif isinstance(attr, Iterable) and not isinstance(attr, (str, bytes)):
            names.update(str(item) for item in attr)
    return names


def _require_project_asset(
    project: Any,
    role: str,
    name: str,
    kinds: Sequence[str],
) -> None:
    names = _asset_names(project, kinds)
    if _matches_asset_name(name, names):
        return
    else:
        available = ", ".join(sorted(names)) or "none"
        raise ValueError(
            f"project is missing {role} '{name}' in {', '.join(kinds)} "
            f"(available: {available})"
        )


def _matches_asset_name(name: str, names: set[str]) -> bool:
    if name in names:
        return True
    suffix_matches = [candidate for candidate in names if candidate.rsplit(".", 1)[-1] == name]
    return len(suffix_matches) == 1


def _zone_thickness(project: Any, top_name: str, base_name: str) -> list[float]:
    top = _surface_values(project, top_name)
    base = _surface_values(project, base_name)
    if not top or not base:
        return []
    if len(top) != len(base):
        raise ValueError(
            f"horizon surface size mismatch for '{top_name}' ({len(top)}) and "
            f"'{base_name}' ({len(base)})"
        )
    return [base[i] - top[i] for i in range(len(top))]


def _surface_values(project: Any, name: str) -> list[float]:
    surface = _surface_object(project, name)
    if surface is None:
        return []
    for attr in ("values", "depth_m", "z", "data"):
        value = getattr(surface, attr, None)
        if value is not None and not callable(value):
            return _flatten_numbers(value)
    if isinstance(surface, Mapping):
        for key in ("values", "depth_m", "z", "data"):
            if key in surface:
                return _flatten_numbers(surface[key])
    if isinstance(surface, Iterable) and not isinstance(surface, (str, bytes)):
        return _flatten_numbers(surface)
    return []


def _surface_object(project: Any, name: str) -> Any:
    for attr_name in ("surfaces", "horizons"):
        collection = getattr(project, attr_name, None)
        if collection is None:
            continue
        if isinstance(collection, Mapping) and name in collection:
            return collection[name]
        if isinstance(collection, Mapping):
            matches = [key for key in collection if str(key).rsplit(".", 1)[-1] == name]
            if len(matches) == 1:
                return collection[matches[0]]
        get = getattr(collection, "get", None)
        if callable(get):
            value = get(name)
            if value is not None:
                return value
            names = _asset_names(project, (attr_name,))
            matches = [candidate for candidate in names if candidate.rsplit(".", 1)[-1] == name]
            if len(matches) == 1:
                value = get(matches[0])
                if value is not None:
                    return value
        try:
            return collection[name]
        except Exception:
            pass
    surface_method = getattr(project, "surface", None)
    if callable(surface_method):
        try:
            return surface_method(name)
        except Exception:
            return None
    return None


def _flatten_numbers(value: Any) -> list[float]:
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        out: list[float] = []
        for item in value:
            if isinstance(item, Iterable) and not isinstance(item, (str, bytes)):
                out.extend(_flatten_numbers(item))
            else:
                out.append(float(item))
        return out
    return []
