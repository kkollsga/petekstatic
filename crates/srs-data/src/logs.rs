//! Extract positioned petro samples from petekio's well curves, for grid population.
//!
//! Pairs each well's **effective** porosity (`PHIE`) and **effective** water-
//! saturation (`SW`) curves into `(z, φ, Sw)` triples using `WellCurveInput.xyz`,
//! skipping NaN values and unpositioned `[NaN;3]` samples. The position `z` is
//! petekio's **negative-down subsea elevation** (its `xyz()` convention, ≥0.3.0),
//! passed through unchanged — pure mapping. The downstream consumer (srs-core
//! `RefiningModel::with_logs`) does the cell binning + upscaling, aligning `z` to
//! the model's positive-down `depth_m` datum at that seam.
//!
//! ## Curve selection (D2)
//! Mnemonics are matched through petekio's [`canonical_mnemonic`] — the one alias
//! table — rather than a hardcoded local set. So vendor variants resolve
//! (`EFFPHI`/`PHIEF`/`PHI` → `PHIE`; `SWE`/`SUWI`/`SW_E` → `SW`), and the
//! **effective/total distinction is preserved**: total porosity (`PHIT`) and
//! total water saturation (`SWT`) canonicalize to distinct names and are *not*
//! folded into the effective sample.

use crate::petekio::{canonical_mnemonic, WellCurveInput};

/// True if `mnemonic` resolves to **effective** porosity (`PHIE`). Total porosity
/// (`PHIT`) is a physically distinct curve and is deliberately excluded.
fn is_effective_porosity(mnemonic: &str) -> bool {
    canonical_mnemonic(mnemonic) == "PHIE"
}

/// True if `mnemonic` resolves to **effective** water saturation (`SW`). Total
/// water saturation (`SWT`) is deliberately excluded.
fn is_effective_sw(mnemonic: &str) -> bool {
    canonical_mnemonic(mnemonic) == "SW"
}

/// Paired `(z_m, porosity, water_saturation)` samples across all wells. The well
/// position `z` is petekio's negative-down subsea elevation, passed through
/// unchanged (see the module note); the flip onto the model's positive-down
/// `depth_m` datum is the downstream binning seam's job (srs-core), mirroring the
/// surface-ingest flip in [`crate::wireframe`].
///
/// Assumes a well's `PHIE`/`SW` curves share the same sample grid (same LAS) — they
/// are paired by index over the shorter of the two. Samples with a non-finite value
/// or an unpositioned `[NaN;3]` location are dropped.
#[must_use]
pub fn petro_samples(curves: &[WellCurveInput]) -> Vec<(f64, f64, f64)> {
    let mut well_ids: Vec<&str> = Vec::new();
    for c in curves {
        if !well_ids.contains(&c.well_id.as_str()) {
            well_ids.push(&c.well_id);
        }
    }

    let mut out = Vec::new();
    for wid in well_ids {
        let phi = curves
            .iter()
            .find(|c| c.well_id == wid && is_effective_porosity(&c.mnemonic));
        let sw = curves
            .iter()
            .find(|c| c.well_id == wid && is_effective_sw(&c.mnemonic));
        let (Some(phi), Some(sw)) = (phi, sw) else {
            continue; // need both to form a (φ, Sw) sample
        };
        let n = phi.values.len().min(sw.values.len()).min(phi.xyz.len());
        for i in 0..n {
            let [x, y, z] = phi.xyz[i];
            let (p, s) = (phi.values[i], sw.values[i]);
            if x.is_finite() && y.is_finite() && z.is_finite() && p.is_finite() && s.is_finite() {
                out.push((z, p, s));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::petekio::Provenance;

    fn curve(mnemonic: &str, vals: Vec<f64>, zs: Vec<f64>) -> WellCurveInput {
        WellCurveInput {
            well_id: "A1".to_string(),
            mnemonic: mnemonic.to_string(),
            md: zs.clone(),
            values: vals,
            xyz: zs.iter().map(|&z| [100.0, 200.0, z]).collect(),
            provenance: Provenance::HardData,
        }
    }

    #[test]
    fn pairs_phie_and_sw_skipping_nan() {
        let phie = curve(
            "PHIE",
            vec![0.24, f64::NAN, 0.21],
            vec![8000.0, 8010.0, 8020.0],
        );
        let sw = curve("SW", vec![0.28, 0.30, 0.31], vec![8000.0, 8010.0, 8020.0]);
        let s = petro_samples(&[phie, sw]);
        // middle sample dropped (φ NaN); two survive
        assert_eq!(s, vec![(8000.0, 0.24, 0.28), (8020.0, 0.21, 0.31)]);
    }

    #[test]
    fn unpositioned_xyz_is_dropped() {
        let mut phie = curve("PHIE", vec![0.24, 0.22], vec![8000.0, 8010.0]);
        phie.xyz[1] = [f64::NAN; 3]; // unpositioned
        let sw = curve("SW", vec![0.28, 0.30], vec![8000.0, 8010.0]);
        let s = petro_samples(&[phie, sw]);
        assert_eq!(s, vec![(8000.0, 0.24, 0.28)]);
    }

    #[test]
    fn no_sw_curve_yields_nothing() {
        let phie = curve("PHIE", vec![0.24], vec![8000.0]);
        assert!(petro_samples(&[phie]).is_empty());
    }

    #[test]
    fn resolves_vendor_mnemonic_variants_through_canonical() {
        // D2: EFFPHI → PHIE and SUWI → SW via petekio's canonical table, so a
        // well logged with vendor variants still pairs (no hardcoded local set).
        let phi = curve("EFFPHI", vec![0.24, 0.21], vec![8000.0, 8020.0]);
        let sw = curve("SUWI", vec![0.28, 0.31], vec![8000.0, 8020.0]);
        let s = petro_samples(&[phi, sw]);
        assert_eq!(s, vec![(8000.0, 0.24, 0.28), (8020.0, 0.21, 0.31)]);
    }

    #[test]
    fn total_porosity_is_not_folded_into_effective() {
        // D2: PHIT (total) canonicalizes to a distinct name, so it must NOT be
        // taken as the effective-porosity curve — a well with only PHIT + SW
        // yields no effective (φ, Sw) sample.
        let phit = curve("PHIT", vec![0.30, 0.28], vec![8000.0, 8020.0]);
        let sw = curve("SW", vec![0.28, 0.31], vec![8000.0, 8020.0]);
        assert!(petro_samples(&[phit, sw]).is_empty());
    }
}
