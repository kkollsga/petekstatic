//! The petroleum P90/P50/P10 summary — a thin, typed wrapper over petekTools'
//! `reservoir_summary`.
//!
//! Petroleum convention: **P90 is the low (conservative) estimate** — the value
//! exceeded with 90% probability, i.e. the 10th percentile of the sorted values;
//! **P10 is the high** (90th percentile). P50 is the median. The percentile
//! kernel itself (type-7 / Hyndman-Fan) lives once in petekTools; this crate only
//! adds the typed empty-input error the FFI surface needs (H1).

use petekstatic_error::StaticError;

/// The P90/P50/P10 + mean summary of a realization set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PercentileSummary {
    /// P90 — low estimate (10th percentile).
    pub p90: f64,
    /// P50 — median.
    pub p50: f64,
    /// P10 — high estimate (90th percentile).
    pub p10: f64,
    /// Arithmetic mean.
    pub mean: f64,
}

impl PercentileSummary {
    /// Summarise realizations (order-independent; a sorted copy is made).
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if `values` is empty — a percentile of zero
    /// realizations is undefined. This is the typed error the FFI surfaces as a
    /// clean `ValueError` when a caller asks for `realizations = 0` (H1), rather
    /// than the historical panic that aborted the interpreter.
    pub fn from_realizations(values: &[f64]) -> Result<Self, StaticError> {
        if values.is_empty() {
            return Err(StaticError::InvalidInput(
                "cannot summarise zero realizations (need at least 1)".into(),
            ));
        }
        // Delegate the P90/P50/P10/mean digest (and its type-7 percentile kernel)
        // to petekTools — the one home for the reservoir summary. The empty guard
        // above keeps the FFI-visible typed error (H1); for a non-empty input this
        // never errors (the `?` maps petekTools' `AlgoError` via `#[from]`).
        let s = petektools::sampling::reservoir_summary(values)?;
        Ok(Self {
            p90: s.p90,
            p50: s.p50,
            p10: s.p10,
            mean: s.mean,
        })
    }

    /// Swanson's mean check `0.3*P10 + 0.4*P50 + 0.3*P90` (spec sanity test).
    #[must_use]
    pub fn swanson_mean(&self) -> f64 {
        0.3 * self.p10 + 0.4 * self.p50 + 0.3 * self.p90
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn petroleum_convention_ordering() {
        let v: Vec<f64> = (1..=1000).map(f64::from).collect();
        let s = PercentileSummary::from_realizations(&v).unwrap();
        // P90 (low) < P50 < P10 (high).
        assert!(s.p90 < s.p50 && s.p50 < s.p10);
        assert!((s.p90 - 100.9).abs() < 1.0);
        assert!((s.p10 - 900.1).abs() < 1.0);
        assert!((s.mean - 500.5).abs() < 1e-9);
    }

    #[test]
    fn zero_realizations_is_a_typed_error_not_a_panic() {
        // H1: an empty realization set must surface a typed error (which the FFI
        // maps to a clean ValueError), never an interpreter-aborting panic.
        let err = PercentileSummary::from_realizations(&[]).unwrap_err();
        assert!(matches!(err, StaticError::InvalidInput(_)));
    }
}
