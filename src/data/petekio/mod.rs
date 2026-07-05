//! Re-export facade for petekio's model-ready input contract.
//!
//! Was a hand-written stub of the locked shape; now that petekio 0.2.1 is published
//! it re-exports the **real** crate's types, so the adapter/wireframe consume the
//! genuine contract (`task_integrate_petekio`). The `crate::data::petekio::X` import path is
//! kept so the rest of srs-data is unchanged by the swap.
//!
//! Construction lives in petekio (`GeoData::new(...).load_*(...).model_inputs()`);
//! srs-data only *consumes* what `model_inputs()` returns.

pub use petekio::{
    Distribution, GeoData, GridGeometry, HorizonInput, ModelInputs, PolygonSet, Provenance,
    SpatialInputs, SummaryInputs, Surface, Uncertain, Unit, WellCurveInput,
};

// The one canonical-mnemonic alias table (petekio owns curve-name normalization);
// srs-data routes its PHIE/SW curve matching through it rather than hardcoding a
// second, divergent alias set (D2).
pub use petekio::analysis::normalize::canonical_mnemonic;
