//! Segments: fault-bounded compartments of the grid. The box grid is a single
//! unfaulted segment; the convergent gridder (Phase 7) introduces fault offset
//! and throw, splitting pillars into multiple segments. Minimal here by design.

/// A fault-bounded compartment. `throw_m` is the vertical offset across its
/// bounding fault, in metres (0 for the unfaulted box).
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub id: usize,
    pub name: String,
    pub throw_m: f64,
}

impl Segment {
    /// The single unfaulted segment covering the whole box grid.
    #[must_use]
    pub fn whole_grid() -> Self {
        Self {
            id: 0,
            name: "segment_0".to_string(),
            throw_m: 0.0,
        }
    }
}
