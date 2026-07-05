//! The thin adapter: map petekio's summary [`SummaryInputs`] into srs input types.
//!
//! **Mapping only — no data processing.** Each petikio [`Uncertain`] becomes an
//! [`InputScalar`] carrying the deterministic `value` (for the box/grid path), the
//! neutral petekio [`Distribution`] DTO (for the downstream Monte Carlo P-curve), and
//! the [`Hardness`] lineage (for the trust surface). The DTO→sampler mapping is
//! **not** srs-data's job — it lives in the consuming simulation layer (petekSim's
//! srs-core), so this crate never depends on the sampler (`srs-uncertainty`). Spatial
//! → wireframe assembly is Phase 3 ([`crate::data::wireframe`]).

use crate::data::petekio::{Distribution, Provenance, SummaryInputs, Uncertain};
use crate::wireframe::Hardness;

/// petikio `Provenance` → wireframe `Hardness`. The one classification mapping the
/// gridder/trust surface need; petikio owns the underlying determination.
#[must_use]
pub fn hardness_of(p: Provenance) -> Hardness {
    match p {
        Provenance::HardData => Hardness::Hard,
        Provenance::Interpolated => Hardness::Interpolated,
        // Defaulted/Assumed are both "not constrained by data".
        Provenance::Defaulted | Provenance::Assumed => Hardness::Assumed,
    }
}

/// An adapted scalar input: deterministic value + its (neutral) sampling
/// distribution DTO + hardness lineage. The srs-side form of a petikio
/// [`Uncertain`]. The [`Distribution`] is petekio's neutral DTO carried straight
/// through; the DTO→sampler conversion happens downstream in petekSim (srs-core).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InputScalar {
    pub value: f64,
    pub distribution: Distribution,
    pub hardness: Hardness,
}

impl InputScalar {
    /// Adapt a single petikio [`Uncertain`] — pure field mapping, the distribution
    /// DTO passes straight through.
    #[must_use]
    pub fn from_uncertain(u: Uncertain) -> Self {
        Self {
            value: u.value,
            distribution: u.distribution,
            hardness: hardness_of(u.provenance),
        }
    }
}

/// The adapted summary inputs — the srs-side mirror of [`SummaryInputs`].
///
/// SI throughout, mirroring petekio ≥0.3.0's metric `SummaryInputs`: `area_m2`,
/// `net_pay_m`, and positive-down `*_depth_m` contacts pass straight through
/// (the acres→m² / ft→m seam shim retired 2026-07-04 when petekio went SI).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModelScalars {
    pub area_m2: InputScalar,
    pub net_pay_m: InputScalar,
    pub porosity_frac: InputScalar,
    pub water_saturation_frac: InputScalar,
    pub net_to_gross_frac: InputScalar,
    pub owc_depth_m: Option<InputScalar>,
    pub goc_depth_m: Option<InputScalar>,
}

impl ModelScalars {
    /// Adapt petikio's [`SummaryInputs`]. Pure field-by-field mapping (petekio is
    /// SI-native, so no unit conversion happens here).
    #[must_use]
    pub fn from_summary(s: &SummaryInputs) -> Self {
        let opt = |o: Option<Uncertain>| o.map(InputScalar::from_uncertain);
        Self {
            area_m2: InputScalar::from_uncertain(s.area_m2),
            net_pay_m: InputScalar::from_uncertain(s.net_pay_m),
            porosity_frac: InputScalar::from_uncertain(s.porosity_frac),
            water_saturation_frac: InputScalar::from_uncertain(s.water_saturation_frac),
            net_to_gross_frac: InputScalar::from_uncertain(s.net_to_gross_frac),
            owc_depth_m: opt(s.owc_depth_m),
            goc_depth_m: opt(s.goc_depth_m),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::petekio::SummaryInputs;

    /// A small SummaryInputs built directly (all fields pub) — controlled inputs for
    /// the mapping logic, no file load needed.
    fn summary() -> SummaryInputs {
        let u = |value, distribution, provenance| Uncertain {
            value,
            distribution,
            provenance,
        };
        SummaryInputs {
            area_m2: u(
                2_509_000.0,
                Distribution::Triangular {
                    lo: 2_023_000.0,
                    mode: 2_509_000.0,
                    hi: 3_157_000.0,
                },
                Provenance::Interpolated,
            ),
            net_pay_m: u(
                25.0,
                Distribution::Normal {
                    mean: 25.0,
                    std: 3.0,
                },
                Provenance::HardData,
            ),
            porosity_frac: u(
                0.22,
                Distribution::Normal {
                    mean: 0.22,
                    std: 0.02,
                },
                Provenance::HardData,
            ),
            water_saturation_frac: u(
                0.30,
                Distribution::Normal {
                    mean: 0.30,
                    std: 0.04,
                },
                Provenance::HardData,
            ),
            net_to_gross_frac: u(
                0.80,
                Distribution::Triangular {
                    lo: 0.6,
                    mode: 0.8,
                    hi: 0.92,
                },
                Provenance::Interpolated,
            ),
            owc_depth_m: Some(u(2511.6, Distribution::Deterministic, Provenance::HardData)),
            goc_depth_m: None,
        }
    }

    #[test]
    fn provenance_maps_to_hardness() {
        assert_eq!(hardness_of(Provenance::HardData), Hardness::Hard);
        assert_eq!(
            hardness_of(Provenance::Interpolated),
            Hardness::Interpolated
        );
        assert_eq!(hardness_of(Provenance::Defaulted), Hardness::Assumed);
        assert_eq!(hardness_of(Provenance::Assumed), Hardness::Assumed);
    }

    #[test]
    fn summary_adapts_field_by_field() {
        let s = ModelScalars::from_summary(&summary());
        assert_eq!(s.area_m2.value, 2_509_000.0);
        assert_eq!(s.area_m2.hardness, Hardness::Interpolated);
        assert_eq!(s.net_pay_m.hardness, Hardness::Hard);
        // owc present, goc absent — the Option seam survives the mapping.
        assert!(s.owc_depth_m.is_some());
        assert!(s.goc_depth_m.is_none());
        // The neutral DTO is carried straight through, unconverted.
        assert_eq!(
            s.owc_depth_m.unwrap().distribution,
            Distribution::Deterministic
        );
    }
}
