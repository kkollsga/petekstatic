//! The single petekStatic error type. Every geomodel crate surfaces failures as
//! [`StaticError`].

use thiserror::Error;

/// petekStatic-wide result alias.
pub type Result<T> = std::result::Result<T, StaticError>;

/// Errors raised across the petekStatic geomodel crates.
#[derive(Debug, Error)]
pub enum StaticError {
    /// A caller-supplied value was outside its valid range (with a reason).
    #[error("invalid input: {0}")]
    InvalidInput(String),

    /// A grid dimension or index was out of bounds.
    #[error("grid error: {0}")]
    Grid(String),

    /// A value fell outside a correlation/spec's documented validity range.
    #[error("outside validity range: {0}")]
    OutOfRange(String),

    /// The base surface crosses above the top at one or more lattice nodes, so
    /// the gross column has negative (crossed) thickness there — a thin/crossing
    /// framework that would silently collapse GRV. Reports the offending node
    /// count and the worst (most negative) separation in metres. Opt into
    /// `clamp_base_to_top` to zero those columns instead of erroring.
    #[error("base surface crosses above the top at {nodes} node(s); worst separation {worst_m} m (thin/crossing surfaces collapse GRV)")]
    CrossedSurfaces {
        /// Number of lattice nodes where the base sits above the top.
        nodes: usize,
        /// The worst (most negative) base-minus-top separation encountered, metres.
        worst_m: f64,
    },

    /// A Monte-Carlo realization draw failed. **Fail-fast policy** (the MC driver,
    /// `task_peteksim_mc_structured`): the loop stops at the first bad draw and
    /// reports its `index` alongside the underlying cause, so the H2 typed error
    /// (an out-of-range priors draw, a crossed base, an off-lattice shift) reaches
    /// the caller *with the offending draw identified* rather than silently
    /// aborting the whole run. `source()` reaches the original failure.
    #[error("Monte-Carlo draw #{index} failed: {source}")]
    McDraw {
        /// The index of the failing draw in the realization loop.
        index: usize,
        /// The underlying failure (validation, grid, crossed surfaces, …).
        #[source]
        source: Box<StaticError>,
    },

    /// A failure from the petekIO DATA layer (ingest/normalize/interpret),
    /// composed across the seam so `?` chains DATA→GEOMODEL and `source()`
    /// reaches the origin (house-style §1). This transitively lets petekSim's
    /// `SrsError` reach `GeoError` too (it already has `#[from] StaticError`).
    #[error(transparent)]
    Geo(#[from] petekio::GeoError),

    /// A failure from the petekTools TOOLKIT layer (sampling / stats / numeric
    /// kernels), composed across the horizontal seam so `?` chains
    /// TOOLKIT→GEOMODEL. The MC driver's percentile / `reservoir_summary` /
    /// `aggregate` calls surface here.
    #[error(transparent)]
    Algo(#[from] petektools::AlgoError),
}
