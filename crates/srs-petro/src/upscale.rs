//! Property-appropriate upscaling, each weighted by what it conserves
//! (`log_upscaling_spec`): porosity by length, Sw by pore volume, NTG by net
//! length. This avoids double-counting net pore volume downstream in GRV.

use crate::power_mean::{arithmetic_mean, harmonic_mean, WeightedSample};
use petekstatic_error::StaticError;

/// A net-flagged sample (e.g. above a porosity/Vshale cutoff).
#[derive(Debug, Clone, Copy)]
pub struct NetSample {
    pub length: f64,
    pub is_net: bool,
}

/// A sample carrying both porosity and water saturation.
#[derive(Debug, Clone, Copy)]
pub struct SwSample {
    pub length: f64,
    pub porosity: f64,
    pub water_saturation: f64,
}

/// Porosity upscales as the length-weighted arithmetic mean.
///
/// # Errors
/// As [`arithmetic_mean`].
pub fn upscale_porosity(samples: &[WeightedSample]) -> Result<f64, StaticError> {
    arithmetic_mean(samples)
}

/// Net-to-gross = net length / gross length, in `[0, 1]`.
///
/// # Errors
/// [`StaticError::InvalidInput`] if empty, any length invalid, or gross <= 0.
pub fn net_to_gross(samples: &[NetSample]) -> Result<f64, StaticError> {
    if samples.is_empty() {
        return Err(StaticError::InvalidInput("no samples for NTG".into()));
    }
    let mut gross = 0.0;
    let mut net = 0.0;
    for s in samples {
        if !s.length.is_finite() || s.length < 0.0 {
            return Err(StaticError::InvalidInput(format!(
                "sample length must be finite and >= 0, got {}",
                s.length
            )));
        }
        gross += s.length;
        if s.is_net {
            net += s.length;
        }
    }
    if gross <= 0.0 {
        return Err(StaticError::InvalidInput("gross length must be > 0".into()));
    }
    Ok(net / gross)
}

/// Pore-volume-weighted water saturation: `Sw = sum(h phi Sw) / sum(h phi)`.
///
/// # Errors
/// [`StaticError::InvalidInput`] if empty, a length/porosity is invalid, or total
/// pore volume is zero.
pub fn upscale_sw(samples: &[SwSample]) -> Result<f64, StaticError> {
    if samples.is_empty() {
        return Err(StaticError::InvalidInput("no samples for Sw".into()));
    }
    let mut pv = 0.0;
    let mut num = 0.0;
    for s in samples {
        if !s.length.is_finite() || s.length < 0.0 || !(0.0..=1.0).contains(&s.porosity) {
            return Err(StaticError::InvalidInput(format!(
                "invalid Sw sample: length={}, porosity={}",
                s.length, s.porosity
            )));
        }
        let w = s.length * s.porosity;
        pv += w;
        num += w * s.water_saturation;
    }
    if pv <= 0.0 {
        return Err(StaticError::InvalidInput(
            "total pore volume is zero; Sw undefined".into(),
        ));
    }
    Ok(num / pv)
}

/// Cardwell-Parsons bounds bracketing effective permeability:
/// `(harmonic, arithmetic)` = (lower, upper).
///
/// # Errors
/// As the underlying means (values must be > 0).
pub fn perm_bounds(samples: &[WeightedSample]) -> Result<(f64, f64), StaticError> {
    Ok((harmonic_mean(samples)?, arithmetic_mean(samples)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ntg_counts_net_length() {
        let s = [
            NetSample {
                length: 10.0,
                is_net: true,
            },
            NetSample {
                length: 30.0,
                is_net: false,
            },
            NetSample {
                length: 10.0,
                is_net: true,
            },
        ];
        // 20 net / 50 gross = 0.4.
        assert!((net_to_gross(&s).unwrap() - 0.4).abs() < 1e-12);
    }

    #[test]
    fn sw_is_pore_volume_weighted() {
        // High-phi wet sand dominates over low-phi tight oil.
        let s = [
            SwSample {
                length: 1.0,
                porosity: 0.30,
                water_saturation: 0.8,
            },
            SwSample {
                length: 1.0,
                porosity: 0.05,
                water_saturation: 0.2,
            },
        ];
        // (0.30*0.8 + 0.05*0.2) / (0.30 + 0.05) = 0.25/0.35.
        assert!((upscale_sw(&s).unwrap() - 0.25 / 0.35).abs() < 1e-12);
    }

    #[test]
    fn perm_bounds_bracket() {
        let d = [
            WeightedSample::new(1.0, 50.0),
            WeightedSample::new(1.0, 500.0),
        ];
        let (lo, hi) = perm_bounds(&d).unwrap();
        assert!(lo < hi);
        assert!((hi - 275.0).abs() < 1e-9); // arithmetic
    }

    #[test]
    fn sw_rejects_zero_pore_volume() {
        let s = [SwSample {
            length: 1.0,
            porosity: 0.0,
            water_saturation: 0.5,
        }];
        assert!(upscale_sw(&s).is_err());
    }
}
