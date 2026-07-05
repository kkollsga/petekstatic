//! Inverse standard-normal CDF (probit) via Acklam's rational approximation.
//! Accurate to ~1.15e-9 over `(0, 1)`; used to inverse-CDF sample the normal and
//! lognormal distributions.

// Acklam (2003) coefficients.
const A: [f64; 6] = [
    -3.969_683_028_665_376e1,
    2.209_460_984_245_205e2,
    -2.759_285_104_469_687e2,
    1.383_577_518_672_69e2,
    -3.066_479_806_614_716e1,
    2.506_628_277_459_239e0,
];
const B: [f64; 5] = [
    -5.447_609_879_822_406e1,
    1.615_858_368_580_409e2,
    -1.556_989_798_598_866e2,
    6.680_131_188_771_972e1,
    -1.328_068_155_288_572e1,
];
const C: [f64; 6] = [
    -7.784_894_002_430_293e-3,
    -3.223_964_580_411_365e-1,
    -2.400_758_277_161_838e0,
    -2.549_732_539_343_734e0,
    4.374_664_141_464_968e0,
    2.938_163_982_698_783e0,
];
const D: [f64; 4] = [
    7.784_695_709_041_462e-3,
    3.224_671_290_700_398e-1,
    2.445_134_137_142_996e0,
    3.754_408_661_907_416e0,
];

/// Inverse CDF of the standard normal distribution for `p` in `(0, 1)`.
///
/// Returns `-inf`/`+inf` at the boundaries `p <= 0` / `p >= 1`.
#[must_use]
pub fn inverse_normal_cdf(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    // Break-points where the central rational region gives way to the tails.
    const P_LOW: f64 = 0.024_25;
    const P_HIGH: f64 = 1.0 - P_LOW;

    if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= P_HIGH {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_is_zero() {
        assert!(inverse_normal_cdf(0.5).abs() < 1e-9);
    }

    #[test]
    fn symmetric_quantiles() {
        // Phi^-1(0.975) ~ 1.959964, the canonical 95% z-value.
        assert!((inverse_normal_cdf(0.975) - 1.959_963_98).abs() < 1e-6);
        assert!((inverse_normal_cdf(0.025) + 1.959_963_98).abs() < 1e-6);
    }

    #[test]
    fn one_sigma() {
        // Phi^-1(0.8413) ~ 1.0.
        assert!((inverse_normal_cdf(0.841_344_746) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn boundaries_are_infinite() {
        assert_eq!(inverse_normal_cdf(0.0), f64::NEG_INFINITY);
        assert_eq!(inverse_normal_cdf(1.0), f64::INFINITY);
    }
}
