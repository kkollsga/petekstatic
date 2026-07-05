//! `srs-uncertainty` — the Monte Carlo toolkit (an independent leaf).
//!
//! [`Distribution`]s sampled by inverse-CDF transform of a seeded [`SplitMix64`]
//! stream, driven by [`run`], summarised as a petroleum [`PercentileSummary`]
//! (P90 low / P50 / P10 high). Inputs-as-distributions per
//! `monte_carlo_volumetrics_spec`; the volumetric coupling lives in `srs-core`.
//!
//! Deferred (post-skeleton): Latin-Hypercube sampling and Iman-Conover rank
//! correlation (the spec's optional refinements).

mod distribution;
mod gaussian;
mod monte_carlo;
mod rng;
mod summary;

pub use distribution::Distribution;
pub use gaussian::inverse_normal_cdf;
pub use monte_carlo::{run, Realizations};
pub use rng::SplitMix64;
pub use summary::PercentileSummary;
