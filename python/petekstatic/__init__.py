"""petekStatic — the geomodel layer of the petek subsurface-modelling stack.

A thin Python surface over the Rust `petekStatic` library: build a populated
static reservoir model and read its volumes (GRV / HCPV / OOIP / OGIP) and the
JSON view bundles (map / cross-section / volume). The rich Python product is
`peteksim`; this wheel is the essentials only.
See https://github.com/kkollsga/petekstatic.
"""

from ._petekstatic import (
    StaticModel,
    __version__,
    build_flat_model,
)

__all__ = [
    "StaticModel",
    "__version__",
    "build_flat_model",
]
