//! The Monte Carlo driver: run a seeded sampling function `N` times and collect
//! the realizations. Generic over what each trial computes, so the volumetric
//! coupling (draw A,h,phi,Sw,FVF -> V) lives in the orchestrator (`srs-core`),
//! keeping this crate an independent leaf.

use crate::error::StaticError;
use crate::uncertainty::rng::SplitMix64;
use crate::uncertainty::summary::PercentileSummary;

/// The realizations produced by a Monte Carlo run.
#[derive(Debug, Clone)]
pub struct Realizations {
    /// One value per trial, in draw order.
    pub values: Vec<f64>,
}

impl Realizations {
    /// The P90/P50/P10 + mean summary.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if there are no realizations (H1) — an empty
    /// set has no percentiles; the FFI surfaces this as a clean `ValueError`.
    pub fn summary(&self) -> Result<PercentileSummary, StaticError> {
        PercentileSummary::from_realizations(&self.values)
    }
}

/// Run `n` trials with a seeded generator, calling `trial` each time.
///
/// Same `seed` + same `trial` => identical realizations (reproducibility).
pub fn run<F>(n: usize, seed: u64, mut trial: F) -> Realizations
where
    F: FnMut(&mut SplitMix64) -> f64,
{
    let mut rng = SplitMix64::new(seed);
    let values = (0..n).map(|_| trial(&mut rng)).collect();
    Realizations { values }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::uncertainty::distribution::Distribution;

    #[test]
    fn reproducible_for_same_seed() {
        let d = Distribution::lognormal(3.0, 0.4).unwrap();
        let a = run(5_000, 99, |r| d.sample(r));
        let b = run(5_000, 99, |r| d.sample(r));
        assert_eq!(a.values, b.values);
    }

    #[test]
    fn recovers_known_normal_percentiles() {
        // Normal(mean=100, sd=10): P50~100, P10(high)~100+1.2816*10, P90(low)~100-1.2816*10.
        let d = Distribution::normal(100.0, 10.0).unwrap();
        let s = run(200_000, 7, |r| d.sample(r)).summary().unwrap();
        assert!((s.p50 - 100.0).abs() < 0.3, "p50={}", s.p50);
        assert!((s.p10 - 112.816).abs() < 0.3, "p10={}", s.p10);
        assert!((s.p90 - 87.184).abs() < 0.3, "p90={}", s.p90);
        assert!((s.mean - 100.0).abs() < 0.2, "mean={}", s.mean);
    }

    #[test]
    fn swanson_mean_near_true_mean_for_symmetric() {
        let d = Distribution::normal(50.0, 5.0).unwrap();
        let s = run(200_000, 3, |r| d.sample(r)).summary().unwrap();
        assert!((s.swanson_mean() - 50.0).abs() < 0.3);
    }
}
