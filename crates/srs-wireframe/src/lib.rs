//! srs-wireframe — the **constraining wireframe** contract.
//!
//! The wireframe is the structural skeleton the convergent gridder and the 3D
//! grid build conform to: an areal **boundary**, a set of depth **horizons**
//! (gridded on the model lattice), and fluid **contacts** — every constraint
//! tagged with how *hard* it is (measured data vs interpolated vs assumed).
//!
//! ## Where it comes from (GATE-0 seam)
//! The wireframe is assembled (in `srs-data`/`srs-wireframe` builders) from
//! petekio's model-ready inputs — a thin mapping, no data processing on our side:
//!
//! | petekio `ModelInputs`                    | wireframe |
//! |------------------------------------------|-----------|
//! | `SpatialInputs.boundary`                 | [`Boundary`] |
//! | `SpatialInputs.horizons` (resampled to our lattice) | [`Horizon`] |
//! | `SummaryInputs.owc_depth_m` / `goc_depth_m` | [`Contact`] |
//! | `Provenance` (hard/interpolated/…)       | [`Hardness`] |
//!
//! This crate is self-contained (no petekio dependency); the adapter bridges.

#![forbid(unsafe_code)]

use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// How hard a constraint is — the data-vs-interpolated flag the gridder needs to
/// know which nodes to honour exactly. Mirrors petekio `Provenance`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Hardness {
    /// Honoured exactly — measured data (a well pick, a mapped point).
    Hard,
    /// Interpolated/gridded between hard data.
    Interpolated,
    /// Assumed where no data constrains it.
    Assumed,
}

/// The structural role of a horizon in the stratigraphic column.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum HorizonRole {
    Top,
    Base,
    Intermediate,
}

/// A fluid contact type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ContactKind {
    /// Oil–water contact.
    Owc,
    /// Gas–oil contact.
    Goc,
    /// Gas–water contact.
    Gwc,
}

/// Areal outline of the modelled region (closed ring of `[x, y]`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Boundary {
    pub ring: Vec<[f64; 2]>,
    pub hardness: Hardness,
}

/// A depth surface gridded on the model lattice (row-major `ncol × nrow`).
/// `depth_m[k]` is `NaN` where undefined; `is_control[k]` marks nodes honoured
/// exactly to hard data.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GriddedDepth {
    pub ncol: usize,
    pub nrow: usize,
    /// Depth (m, positive down), row-major; `NaN` = undefined.
    pub depth_m: Vec<f64>,
    /// Per-node: was this node pinned to hard data during gridding?
    pub is_control: Vec<bool>,
}

/// A structural horizon constraint.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Horizon {
    pub name: String,
    pub role: HorizonRole,
    pub surface: GriddedDepth,
}

/// A fluid contact constraint (a depth scalar, metres positive-down).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Contact {
    pub kind: ContactKind,
    pub depth_m: f64,
    pub hardness: Hardness,
}

/// The constraining wireframe consumed by the convergent gridder + 3D grid build.
///
/// `horizons` is reference-counted (`Arc`): it is **realization-invariant**
/// structural interpretation, so the MC template shares one horizon set across
/// every `realize` (an `O(1)` refcount bump per realization) instead of
/// deep-copying the per-node depth cubes each draw. Only `contacts` (two small
/// scalars) change per realization. Mutate the horizon set in place with
/// `Arc::make_mut` (copy-on-write if it is shared).
#[derive(Debug, Clone, PartialEq)]
pub struct Wireframe {
    pub boundary: Boundary,
    pub horizons: Arc<Vec<Horizon>>,
    pub contacts: Vec<Contact>,
}

impl Wireframe {
    /// An empty wireframe with just an areal boundary — the minimal seed that
    /// later `add_*` passes (wells, surfaces) tighten.
    pub fn from_boundary(boundary: Boundary) -> Self {
        Self {
            boundary,
            horizons: Arc::new(Vec::new()),
            contacts: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seed_from_boundary_is_empty() {
        let b = Boundary {
            ring: vec![[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]],
            hardness: Hardness::Hard,
        };
        let wf = Wireframe::from_boundary(b);
        assert!(wf.horizons.is_empty());
        assert!(wf.contacts.is_empty());
        assert_eq!(wf.boundary.ring.len(), 4);
    }
}
