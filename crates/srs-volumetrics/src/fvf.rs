//! Formation volume factors (FVF) — the reservoir-to-surface volume ratios that
//! convert hydrocarbon pore volume into in-place surface volumes.
//!
//! FVF enters the static-uncertainty layer as an **uncertain scalar input** (the
//! layer charter, graph `decision_layer_charters`): no PVT code crosses the seam
//! down from petekSim. These are small, self-validating value types — the petek
//! house-style "duplicate a small type at the seam" rather than depend sideways on
//! petekSim's `srs-pvt` (which keeps its own copy for its dynamic/PVT-correlation
//! work). The validation here is also the volumetrics half of the physical-range
//! hardening: a non-physical FVF is a typed error, not silent `inf`.

use petekstatic_error::StaticError;

/// Oil formation volume factor `Boi` [reservoir m³ / standard m³ = Rm³/Sm³].
/// Reservoir oil shrinks to the tank, so physically `Boi >= 1`.
///
/// Note: `Boi` in the legacy rb/STB is **numerically identical** to Rm³/Sm³ —
/// both are reservoir-volume / surface-volume ratios (dimensionless), so this
/// is a **relabel, not a conversion** (family SI standard,
/// `decision_si_units_standard`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OilFvf(f64);

/// Gas formation volume factor `Bgi` [reservoir m³ / standard m³ = Rm³/Sm³].
/// Reservoir gas expands at surface, so physically `0 < Bgi < 1`. Like `Boi`
/// this is dimensionless: the rcf/scf value relabels to Rm³/Sm³ unchanged.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasFvf(f64);

impl OilFvf {
    /// Construct an oil FVF.
    ///
    /// # Errors
    /// Returns [`StaticError::InvalidInput`] unless finite and `>= 1.0`.
    pub fn new(rm3_per_sm3: f64) -> Result<Self, StaticError> {
        if rm3_per_sm3.is_finite() && rm3_per_sm3 >= 1.0 {
            Ok(Self(rm3_per_sm3))
        } else {
            Err(StaticError::InvalidInput(format!(
                "oil FVF (Boi) must be finite and >= 1.0 Rm³/Sm³, got {rm3_per_sm3}"
            )))
        }
    }

    /// The value in Rm³/Sm³.
    #[must_use]
    pub fn value(self) -> f64 {
        self.0
    }
}

impl GasFvf {
    /// Construct a gas FVF.
    ///
    /// # Errors
    /// Returns [`StaticError::InvalidInput`] unless finite and in `(0, 1)`.
    pub fn new(rm3_per_sm3: f64) -> Result<Self, StaticError> {
        if rm3_per_sm3.is_finite() && rm3_per_sm3 > 0.0 && rm3_per_sm3 < 1.0 {
            Ok(Self(rm3_per_sm3))
        } else {
            Err(StaticError::InvalidInput(format!(
                "gas FVF (Bgi) must be finite and in (0,1) Rm³/Sm³, got {rm3_per_sm3}"
            )))
        }
    }

    /// The value in Rm³/Sm³.
    #[must_use]
    pub fn value(self) -> f64 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_fvfs_round_trip() {
        assert!((OilFvf::new(1.25).unwrap().value() - 1.25).abs() < 1e-12);
        assert!((GasFvf::new(0.005).unwrap().value() - 0.005).abs() < 1e-12);
    }

    #[test]
    fn oil_fvf_rejects_below_one_and_nonfinite() {
        assert!(OilFvf::new(0.9).is_err());
        assert!(OilFvf::new(0.0).is_err());
        assert!(OilFvf::new(f64::NAN).is_err());
    }

    #[test]
    fn gas_fvf_rejects_out_of_unit_interval() {
        assert!(GasFvf::new(0.0).is_err());
        assert!(GasFvf::new(1.0).is_err());
        assert!(GasFvf::new(f64::INFINITY).is_err());
    }
}
