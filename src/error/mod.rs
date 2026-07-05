//! The petekStatic workspace error type.
//!
//! petekStatic owns one error enum, [`StaticError`], surfaced by every geomodel
//! crate — the analog of petekSim's `SrsError` (petek house style §1). Downstream
//! libraries compose it across the seam with a `#[from]` variant on their own enum.

#[allow(clippy::module_inception)]
mod error;

pub use error::{Result, StaticError};
