//! srs-data — the **thin petekio adapter**.
//!
//! petekio owns *all* input-data work (ingest · normalize · validate · interpret ·
//! characterise · surface ops) and hands us a model-ready [`petekio::ModelInputs`]
//! bundle. srs-data does **zero** data processing — it only *maps* that bundle into
//! srs input types and assembles the constraining [`crate::wireframe::Wireframe`].
//!
//! ## Integrated against petekio 0.2.1
//! Built stub-first, now integrated: [`petekio`] re-exports the **real**
//! `petekio = "0.2.1"` crate (verified end-to-end — see `tests/integration_petekio.rs`,
//! petekio's GATE-0 verification). One contract gap surfaced: `PolygonSet` exposes no
//! ring-vertex accessor, so [`wireframe`] reconstructs the boundary from its bounding
//! box (interim) — raised as `q_petekio_polygon_rings`.
//!
//! Discipline (the boundary rule — "whose type is the output"): if anything here is
//! tempted to parse/alias/harmonise/cutoff/QC/default/reproject, that is a petikio
//! gap → raise a Question / notify, never a local workaround.

#![forbid(unsafe_code)]

pub mod adapter;
pub mod logs;
pub mod petekio;
pub mod wireframe;
