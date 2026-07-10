//! `petekstatic` — the GEOMODEL layer of the petek subsurface-modelling suite,
//! consolidated into **one crate**.
//!
//! Packaging note (2026-07-05): the ten historical workspace crates were merged
//! into this single `petekstatic` crate. Today's boundaries are preserved as
//! **modules** with the same one-directional import discipline they had as crates:
//!
//! ```text
//! error → wireframe → grid → petro → gridder → volumetrics → uncertainty → data → spill → model
//! ```
//!
//! - [`error`] — the one workspace error type ([`StaticError`]).
//! - [`wireframe`] — the constraining wireframe (boundary + horizons + contacts).
//! - [`grid`] — the i,j,k corner-point grid data model.
//! - [`petro`] — petrophysics (log upscaling, facies).
//! - [`gridder`] — the convergent gridder + conformable layering.
//! - [`volumetrics`] — GRV / in-place volumetrics + FVF.
//! - [`uncertainty`] — the Monte Carlo toolkit (distributions, P90/P50/P10).
//! - [`data`] — the thin petekio adapter (model-ready inputs → srs input types).
//! - [`spill`] — the out-of-core backing-storage mode.
//! - [`model`] — the top of the DAG: the [`StaticModel`] aggregate + MC template.
//!
//! The headline API is re-exported at the crate root, so callers reach the common
//! types (`StaticModelBuilder`, `StaticModelTemplate`, the `HorizonStack` family,
//! `run_mc` / `McSettings`, `BuildSpec`, `StaticModel`, the view bundles, …)
//! without needing to know which module they live in.

#[cfg(feature = "petekio-adapter")]
pub mod data;
pub mod error;
pub mod grid;
pub mod gridder;
pub mod model;
pub mod petro;
pub mod spill;
pub mod uncertainty;
pub mod volumetrics;
pub mod wireframe;

// The single workspace error type, at the crate root.
pub use error::{Result, StaticError};

// The full model surface (the top of the DAG) at the crate root — this carries the
// headline API: StaticModelBuilder, StaticModelTemplate, the HorizonStack family,
// run_mc / McSettings, BuildSpec, StaticModel, the view bundles, and the
// re-exported volumetrics / uncertainty / spill / sampling seam types.
pub use model::*;
