//! `srs-grid` — the i,j,k corner-point grid data model.
//!
//! The spine type is [`Grid`]: corner-point [`geometry`](CornerPointGeom) plus
//! per-cell [`Properties`] plus [`KLayer`] layering. The box grid
//! ([`build_box`]) is the degenerate, unfaulted case (Ponting 1989); the same
//! structure carries the full convergent grid later, so it upgrades in place.

mod box_construction;
mod cell;
mod geometry;
#[allow(clippy::module_inception)]
mod grid;
mod index;
mod klayer;
mod point;
mod property;
mod segment;
mod volume;

pub use box_construction::{build_box, BoxSpec};
pub use cell::Cell;
pub use geometry::{CornerPointGeom, Pillar};
pub use grid::Grid;
pub use index::{Dims, Ijk};
pub use klayer::KLayer;
pub use point::Point3;
pub use property::{Properties, Property};
pub use segment::Segment;
pub use volume::hexahedron_volume;
