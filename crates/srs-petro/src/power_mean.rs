//! Generalized (power-law) mean over length-weighted samples
//! (`log_upscaling_spec`).
//!
//! `M_p = ( sum(h_i * v_i^p) / sum(h_i) )^(1/p)`, with the geometric limit at
//! `p -> 0`: `exp( sum(h_i ln v_i) / sum(h_i) )`. `p=+1` arithmetic (porosity,
//! horizontal k), `p=-1` harmonic (vertical k), `p=0` geometric (isotropic k).

use petekstatic_error::StaticError;

/// A sample with a length weight along the well trajectory.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WeightedSample {
    /// Sample length \[ft\] (the weight).
    pub length: f64,
    /// Property value.
    pub value: f64,
}

impl WeightedSample {
    /// Construct a weighted sample.
    #[must_use]
    pub fn new(length: f64, value: f64) -> Self {
        Self { length, value }
    }
}

/// Length-weighted power-law mean with exponent `p` (`p == 0.0` => geometric).
///
/// # Errors
/// [`StaticError::InvalidInput`] if there are no samples, total length is not
/// positive, any length is negative/non-finite, or a value is non-positive when
/// `p <= 0` (geometric/harmonic need strictly positive values).
pub fn power_law_mean(samples: &[WeightedSample], p: f64) -> Result<f64, StaticError> {
    if samples.is_empty() {
        return Err(StaticError::InvalidInput("no samples to upscale".into()));
    }
    let total: f64 = samples.iter().map(|s| s.length).sum();
    if !(total.is_finite() && total > 0.0) {
        return Err(StaticError::InvalidInput(format!(
            "total sample length must be finite and > 0, got {total}"
        )));
    }
    for s in samples {
        if !s.length.is_finite() || s.length < 0.0 {
            return Err(StaticError::InvalidInput(format!(
                "sample length must be finite and >= 0, got {}",
                s.length
            )));
        }
        if p <= 0.0 && !(s.value.is_finite() && s.value > 0.0) {
            return Err(StaticError::InvalidInput(format!(
                "geometric/harmonic mean needs values > 0, got {}",
                s.value
            )));
        }
    }

    if p == 0.0 {
        let acc: f64 = samples.iter().map(|s| s.length * s.value.ln()).sum();
        Ok((acc / total).exp())
    } else {
        let acc: f64 = samples.iter().map(|s| s.length * s.value.powf(p)).sum();
        Ok((acc / total).powf(1.0 / p))
    }
}

/// Length-weighted arithmetic mean (`p = +1`): porosity, horizontal k.
///
/// # Errors
/// As [`power_law_mean`].
pub fn arithmetic_mean(samples: &[WeightedSample]) -> Result<f64, StaticError> {
    power_law_mean(samples, 1.0)
}

/// Length-weighted harmonic mean (`p = -1`): vertical k.
///
/// # Errors
/// As [`power_law_mean`].
pub fn harmonic_mean(samples: &[WeightedSample]) -> Result<f64, StaticError> {
    power_law_mean(samples, -1.0)
}

/// Length-weighted geometric mean (`p -> 0`): isotropic k.
///
/// # Errors
/// As [`power_law_mean`].
pub fn geometric_mean(samples: &[WeightedSample]) -> Result<f64, StaticError> {
    power_law_mean(samples, 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(len: f64, v: f64) -> WeightedSample {
        WeightedSample::new(len, v)
    }

    #[test]
    fn equal_weight_arithmetic() {
        let x = arithmetic_mean(&[s(1.0, 0.2), s(1.0, 0.3)]).unwrap();
        assert!((x - 0.25).abs() < 1e-12);
    }

    #[test]
    fn length_weighting() {
        // 3 ft at 0.1, 1 ft at 0.3 -> (0.3+0.3)/4 = 0.15.
        let x = arithmetic_mean(&[s(3.0, 0.1), s(1.0, 0.3)]).unwrap();
        assert!((x - 0.15).abs() < 1e-12);
    }

    #[test]
    fn harmonic_le_geometric_le_arithmetic() {
        // Power-mean inequality for positive, non-constant values.
        let d = [s(1.0, 10.0), s(1.0, 1000.0)];
        let h = harmonic_mean(&d).unwrap();
        let g = geometric_mean(&d).unwrap();
        let a = arithmetic_mean(&d).unwrap();
        assert!(h < g && g < a, "h={h} g={g} a={a}");
        // Geometric of 10 and 1000 is 100.
        assert!((g - 100.0).abs() < 1e-9);
        // Harmonic of 10 and 1000 ~ 19.8.
        assert!((h - 2.0 / (0.1 + 0.001)).abs() < 1e-9);
    }

    #[test]
    fn nonpositive_value_rejected_for_geometric() {
        assert!(geometric_mean(&[s(1.0, 0.0)]).is_err());
        assert!(harmonic_mean(&[s(1.0, -2.0)]).is_err());
        // arithmetic tolerates zero/negative.
        assert!(arithmetic_mean(&[s(1.0, 0.0)]).is_ok());
    }

    #[test]
    fn empty_rejected() {
        assert!(arithmetic_mean(&[]).is_err());
    }
}
