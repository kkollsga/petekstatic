//! K-layers: the stratigraphic layering of the grid. In the box grid each `k`
//! is one uniform layer; `layer_interpolation` (Phase 7) splits/curves these.

/// One stratigraphic layer, identified by its `k` index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KLayer {
    pub k: usize,
    pub name: String,
}

impl KLayer {
    /// A layer with a default name `layer_{k}`.
    #[must_use]
    pub fn new(k: usize) -> Self {
        Self {
            k,
            name: format!("layer_{k}"),
        }
    }
}
