//! Physical-range validation for volumetric inputs — the canonical predicates
//! every volumetrics path (grid population, in-place, and the consumer's analytic
//! box tail) checks so a non-physical input is a **typed error**, never silent
//! garbage (validation finding H2: φ=-0.1 → -117.7 MMSTB; Sw=1.2 → negative;
//! FVF=0 → `inf`).

use petekstatic_error::StaticError;

/// A saturation / porosity / net-to-gross value must be a finite fraction in
/// `[0, 1]`.
///
/// # Errors
/// [`StaticError::InvalidInput`] if `x` is not finite or is outside `[0, 1]`.
pub fn validate_fraction(what: &str, x: f64) -> Result<(), StaticError> {
    if x.is_finite() && (0.0..=1.0).contains(&x) {
        Ok(())
    } else {
        Err(StaticError::InvalidInput(format!(
            "{what} must be a fraction in [0,1], got {x}"
        )))
    }
}

/// A magnitude (area, gross height, …) must be finite and strictly positive.
///
/// # Errors
/// [`StaticError::InvalidInput`] if `x` is not finite or is `<= 0`.
pub fn validate_positive(what: &str, x: f64) -> Result<(), StaticError> {
    if x.is_finite() && x > 0.0 {
        Ok(())
    } else {
        Err(StaticError::InvalidInput(format!(
            "{what} must be finite and > 0, got {x}"
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fraction_accepts_unit_interval_endpoints() {
        assert!(validate_fraction("phi", 0.0).is_ok());
        assert!(validate_fraction("phi", 1.0).is_ok());
        assert!(validate_fraction("phi", 0.25).is_ok());
    }

    #[test]
    fn fraction_rejects_the_h2_garbage_cases() {
        // The exact garbage the validator drove through the wheel.
        assert!(validate_fraction("porosity", -0.1).is_err()); // -> was -117.7 MMSTB
        assert!(validate_fraction("water_saturation", 1.2).is_err()); // -> was negative
        assert!(validate_fraction("net_to_gross", f64::NAN).is_err());
        assert!(validate_fraction("porosity", f64::INFINITY).is_err());
    }

    #[test]
    fn positive_rejects_zero_negative_and_nonfinite() {
        assert!(validate_positive("area", 100.0).is_ok());
        assert!(validate_positive("area", 0.0).is_err());
        assert!(validate_positive("gross_height", -5.0).is_err());
        assert!(validate_positive("area", f64::INFINITY).is_err());
    }
}
