//! Input distributions, each sampled by inverse-CDF transform of a uniform.
//!
//! `distribution_sampling_spec`: lognormal is the usual default for area/perm,
//! triangular for sparsely-known parameters. A `Constant` covers deterministic
//! inputs so every input can be expressed as a distribution.

use crate::gaussian::inverse_normal_cdf;
use crate::rng::SplitMix64;
use petekstatic_error::StaticError;

/// A 1-D probability distribution over a volumetric input.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Distribution {
    /// A fixed value (deterministic input).
    Constant(f64),
    /// Uniform on `[min, max]`.
    Uniform { min: f64, max: f64 },
    /// Triangular with `min <= mode <= max`.
    Triangular { min: f64, mode: f64, max: f64 },
    /// Normal with given mean and standard deviation.
    Normal { mean: f64, sd: f64 },
    /// Lognormal: `exp(Normal(mu, sigma))`, parameters in log-space.
    Lognormal { mu: f64, sigma: f64 },
}

impl Distribution {
    /// Validate a uniform distribution.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] unless `min < max` and both finite.
    pub fn uniform(min: f64, max: f64) -> Result<Self, StaticError> {
        if min.is_finite() && max.is_finite() && min < max {
            Ok(Self::Uniform { min, max })
        } else {
            Err(StaticError::InvalidInput(format!(
                "uniform requires finite min < max, got [{min}, {max}]"
            )))
        }
    }

    /// Validate a triangular distribution.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] unless `min <= mode <= max`, `min < max`, finite.
    pub fn triangular(min: f64, mode: f64, max: f64) -> Result<Self, StaticError> {
        if [min, mode, max].iter().all(|v| v.is_finite()) && min <= mode && mode <= max && min < max
        {
            Ok(Self::Triangular { min, mode, max })
        } else {
            Err(StaticError::InvalidInput(format!(
                "triangular requires finite min <= mode <= max and min < max, got [{min}, {mode}, {max}]"
            )))
        }
    }

    /// Validate a normal distribution.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] unless `sd > 0` and parameters finite.
    pub fn normal(mean: f64, sd: f64) -> Result<Self, StaticError> {
        if mean.is_finite() && sd.is_finite() && sd > 0.0 {
            Ok(Self::Normal { mean, sd })
        } else {
            Err(StaticError::InvalidInput(format!(
                "normal requires finite mean and sd > 0, got mean={mean}, sd={sd}"
            )))
        }
    }

    /// Validate a lognormal distribution (log-space mu, sigma).
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] unless `sigma > 0` and parameters finite.
    pub fn lognormal(mu: f64, sigma: f64) -> Result<Self, StaticError> {
        if mu.is_finite() && sigma.is_finite() && sigma > 0.0 {
            Ok(Self::Lognormal { mu, sigma })
        } else {
            Err(StaticError::InvalidInput(format!(
                "lognormal requires finite mu and sigma > 0, got mu={mu}, sigma={sigma}"
            )))
        }
    }

    /// The quantile (inverse CDF) at cumulative probability `u` in `[0, 1)`.
    #[must_use]
    pub fn quantile(self, u: f64) -> f64 {
        match self {
            Self::Constant(c) => c,
            Self::Uniform { min, max } => min + u * (max - min),
            Self::Triangular { min, mode, max } => {
                let fc = (mode - min) / (max - min);
                if u < fc {
                    min + (u * (max - min) * (mode - min)).sqrt()
                } else {
                    max - ((1.0 - u) * (max - min) * (max - mode)).sqrt()
                }
            }
            Self::Normal { mean, sd } => mean + sd * inverse_normal_cdf(u),
            Self::Lognormal { mu, sigma } => (mu + sigma * inverse_normal_cdf(u)).exp(),
        }
    }

    /// Draw one sample using the generator.
    pub fn sample(self, rng: &mut SplitMix64) -> f64 {
        self.quantile(rng.next_f64())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructors_validate() {
        assert!(Distribution::uniform(2.0, 1.0).is_err());
        assert!(Distribution::triangular(0.0, 5.0, 4.0).is_err());
        assert!(Distribution::normal(0.0, -1.0).is_err());
        assert!(Distribution::lognormal(0.0, 0.0).is_err());
        assert!(Distribution::uniform(1.0, 2.0).is_ok());
    }

    #[test]
    fn uniform_quantiles_are_linear() {
        let d = Distribution::uniform(10.0, 20.0).unwrap();
        assert!((d.quantile(0.0) - 10.0).abs() < 1e-12);
        assert!((d.quantile(0.5) - 15.0).abs() < 1e-12);
    }

    #[test]
    fn triangular_endpoints_and_mode_band() {
        let d = Distribution::triangular(0.0, 3.0, 6.0).unwrap();
        assert!((d.quantile(0.0) - 0.0).abs() < 1e-9);
        assert!((d.quantile(1.0 - 1e-12) - 6.0).abs() < 1e-3);
        // Symmetric triangle: median == mode.
        assert!((d.quantile(0.5) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn normal_median_is_mean() {
        let d = Distribution::normal(100.0, 15.0).unwrap();
        assert!((d.quantile(0.5) - 100.0).abs() < 1e-7);
    }

    #[test]
    fn lognormal_median_is_exp_mu() {
        // Median of lognormal(mu, sigma) is exp(mu).
        let d = Distribution::lognormal(2.0, 0.5).unwrap();
        assert!((d.quantile(0.5) - (2.0_f64).exp()).abs() < 1e-6);
    }
}
