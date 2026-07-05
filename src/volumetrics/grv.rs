//! GRV and in-place volumetrics from a populated grid
//! (`grv_from_grid_spec`, Cosentino / Ringrose-Bentley / Leverett).
//!
//! Per-cell hydrocarbon pore volume: `HCPV_i = V_i * NTG_i * phi_i * (1 - Sw_i)`.
//! GRV is the bulk volume of cells in the hydrocarbon column. All volumes are
//! m³ internally (family SI standard). Surface in-place: `STOIIP = sum HCPV /
//! Boi`, `GIIP = sum HCPV / Bgi` — HCPV is reservoir m³, Boi/Bgi are Rm³/Sm³, so
//! the quotient is Sm³ directly (no imperial ft³→rb step). Reporting accessors
//! present mcm / MSm³ (oil) / bcm (gas).
//!
//! MVP uses a **hard contact** (a cell is in or out by its centroid depth).
//! Saturation-height partial-fill through the transition zone is deferred
//! (see `consider-for-future.md`).
//!
//! ## Two-contact columns (gas cap + oil rim)
//! [`compute_in_place_two_contact`] partitions cells by centroid depth against a
//! GOC (upper) and an OWC/FWL (lower): gas cap (`z < GOC`), oil leg
//! (`GOC <= z < OWC`), water (`z >= OWC`). Each hydrocarbon zone gets its own
//! GRV + HCPV so the caller can apply the right FVF per zone
//! ([`InPlace::gas_zone_ogip_sm3`] / [`InPlace::oil_zone_ooip_sm3`]). This is a
//! **geometry + in-place split only** — free-gas OGIP off the gas cap and STOIIP
//! off the oil leg. Solution gas (Rs), gas-cap expansion and any other PVT
//! coupling stay in petekSim's PVT domain; no fluid physics crosses this seam.

use crate::error::StaticError;
use crate::grid::Grid;
use crate::volumetrics::backing::{compute_clipped, Clip, GridSource};
use crate::volumetrics::fvf::{GasFvf, OilFvf};
use core::ops::Range;
use petektools::units::{m3_to_bcm, m3_to_mcm, m3_to_msm3};

/// GRV + HCPV of one hydrocarbon zone (a gas cap or an oil leg).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZoneVolumes {
    /// Gross rock volume of this zone \[m³\].
    pub grv_m3: f64,
    /// Hydrocarbon pore volume of this zone \[m³\].
    pub hcpv_m3: f64,
    /// Number of cells in this zone.
    pub cells: usize,
}

/// Deterministic in-place result. Volumes in m³; HCPV is reservoir volume.
///
/// `grv_m3` / `hcpv_m3` are the whole hydrocarbon column (gas + oil). For a
/// two-contact column [`InPlace::gas`] / [`InPlace::oil`] carry the per-zone
/// split; for a single-contact column both are `None` (the column is a generic
/// hydrocarbon leg the caller resolves with the FVF it applies).
#[derive(Debug, Clone)]
pub struct InPlace {
    /// Gross rock volume of the hydrocarbon column \[m³\].
    pub grv_m3: f64,
    /// Total hydrocarbon pore volume \[m³\].
    pub hcpv_m3: f64,
    /// Number of cells in the hydrocarbon column.
    pub cells_in_column: usize,
    /// Per-cell HCPV \[m³\], linear cell order (0 outside the column).
    pub per_cell_hcpv: Vec<f64>,
    /// Gas-cap volumes (two-contact column only; `None` single-contact).
    pub gas: Option<ZoneVolumes>,
    /// Oil-leg volumes (two-contact column only; `None` single-contact).
    pub oil: Option<ZoneVolumes>,
}

impl InPlace {
    /// GRV in **mcm** (million m³).
    #[must_use]
    pub fn grv_mcm(&self) -> f64 {
        m3_to_mcm(self.grv_m3)
    }

    /// Stock-tank oil initially in place off the **whole** column \[Sm³\]. HCPV
    /// is reservoir m³ and `Boi` is Rm³/Sm³, so the quotient is Sm³. For a
    /// two-contact column prefer [`InPlace::oil_zone_ooip_sm3`] (oil leg only).
    #[must_use]
    pub fn ooip_sm3(&self, boi: OilFvf) -> f64 {
        self.hcpv_m3 / boi.value()
    }

    /// Stock-tank oil in place off the **whole** column in **MSm³** (million Sm³)
    /// — the oil reporting scale.
    #[must_use]
    pub fn oil_msm3(&self, boi: OilFvf) -> f64 {
        m3_to_msm3(self.ooip_sm3(boi))
    }

    /// Gas initially in place off the **whole** column \[Sm³\]. For a two-contact
    /// column prefer [`InPlace::gas_zone_ogip_sm3`] (gas cap only).
    #[must_use]
    pub fn ogip_sm3(&self, bgi: GasFvf) -> f64 {
        self.hcpv_m3 / bgi.value()
    }

    /// Gas in place off the **whole** column in **bcm** (billion Sm³) — the gas
    /// reporting scale.
    #[must_use]
    pub fn gas_bcm(&self, bgi: GasFvf) -> f64 {
        m3_to_bcm(self.ogip_sm3(bgi))
    }

    /// Free gas in place off the **gas cap** only \[Sm³\]; `0` if there is no
    /// gas cap (single-contact column).
    #[must_use]
    pub fn gas_zone_ogip_sm3(&self, bgi: GasFvf) -> f64 {
        self.gas.map_or(0.0, |z| z.hcpv_m3) / bgi.value()
    }

    /// Stock-tank oil in place off the **oil leg** only \[Sm³\]; `0` if there is
    /// no distinct oil leg (single-contact column).
    #[must_use]
    pub fn oil_zone_ooip_sm3(&self, boi: OilFvf) -> f64 {
        self.oil.map_or(0.0, |z| z.hcpv_m3) / boi.value()
    }
}

/// Compute GRV + HCPV from a populated grid with a hard hydrocarbon contact.
///
/// Cells whose centroid is shallower than `contact_depth_m` (smaller z) are in
/// the hydrocarbon column. Requires `PORO`, `NTG`, `SW` properties.
///
/// # Errors
/// Returns [`StaticError::InvalidInput`] if a required property is missing or a
/// value is not a fraction in `[0, 1]`, or [`StaticError::OutOfRange`] if the
/// contact depth is not finite.
pub fn compute_in_place(grid: &Grid, contact_depth_m: f64) -> Result<InPlace, StaticError> {
    let nk = grid.dims().nk;
    compute_clipped(
        &GridSource::new(grid),
        Clip::Single(contact_depth_m),
        0..nk,
        true,
    )
}

/// Per-zone single-contact in-place: [`compute_in_place`] restricted to the cells
/// whose k-layer falls in `k_range` (a zone's layer band). Summary-only (no
/// per-cell HCPV cube). The rollup over a zone partition of `[0, nk)` reproduces
/// the whole-grid total (each cell is counted by exactly one zone).
///
/// # Errors
/// Same as [`compute_in_place`].
pub fn compute_in_place_zone(
    grid: &Grid,
    contact_depth_m: f64,
    k_range: Range<usize>,
) -> Result<InPlace, StaticError> {
    compute_clipped(
        &GridSource::new(grid),
        Clip::Single(contact_depth_m),
        k_range,
        false,
    )
}

/// Per-zone two-contact (gas-cap + oil-leg) in-place: [`compute_in_place_two_contact`]
/// restricted to `k_range`. Summary-only.
///
/// # Errors
/// Same as [`compute_in_place_two_contact`].
pub fn compute_in_place_two_contact_zone(
    grid: &Grid,
    goc_m: f64,
    owc_m: f64,
    sw_gas: Option<f64>,
    k_range: Range<usize>,
) -> Result<InPlace, StaticError> {
    compute_clipped(
        &GridSource::new(grid),
        Clip::Two {
            goc: goc_m,
            owc: owc_m,
            sw_gas,
        },
        k_range,
        false,
    )
}

/// The gross bulk volume of a zone's cells (its `k_range`) as an [`InPlace`] with
/// `grv_m3` = summed active-cell bulk volume and `hcpv_m3` = 0 — the volumetric
/// contribution of a **contactless** zone. A zone with no known fluid contact has
/// no known accumulation, so it contributes gross rock but **zero hydrocarbon
/// in-place** (explicitly not treated as a full hydrocarbon column). Truncated
/// (zero-volume) cells are excluded.
#[must_use]
pub fn compute_zone_bulk(grid: &Grid, k_range: Range<usize>) -> InPlace {
    // Bulk reads only geometry (no cubes), so the unified core is infallible here.
    compute_clipped(&GridSource::new(grid), Clip::Bulk, k_range, false)
        .expect("bulk volumetrics reads only geometry — infallible in-core")
}

/// Summary-only variant of [`compute_in_place`]: the same GRV/HCPV totals without
/// materializing the per-cell HCPV cube (V7). `per_cell_hcpv` is left empty — use
/// this on the MC hot path where only aggregates feed the P-curve.
///
/// # Errors
/// Same as [`compute_in_place`].
pub fn compute_in_place_summary(grid: &Grid, contact_depth_m: f64) -> Result<InPlace, StaticError> {
    let nk = grid.dims().nk;
    compute_clipped(
        &GridSource::new(grid),
        Clip::Single(contact_depth_m),
        0..nk,
        false,
    )
}

/// Compute per-zone GRV + HCPV for a **two-contact** column: a gas cap above
/// `goc_m` and an oil leg between `goc_m` and `owc_m` (an OWC or FWL). A cell
/// is gas if its centroid `z < goc_m`, oil if `goc_m <= z < owc_m`, water
/// otherwise. The returned [`InPlace`] carries the split in `gas` / `oil`; its
/// `grv_m3` / `hcpv_m3` totals are gas + oil.
///
/// This is a **geometry + in-place partition only** — no fluid physics. Apply
/// free-gas FVF to the gas cap ([`InPlace::gas_zone_ogip_sm3`]) and oil FVF to
/// the oil leg ([`InPlace::oil_zone_ooip_sm3`]); solution gas and any PVT
/// coupling belong to petekSim.
///
/// `sw_gas` (R3) is an optional gas-cap connate-water override: with `Some(s)`
/// the gas-cap cells use the single scalar `s` for `(1 - Sw)` instead of the
/// shared `SW` cube (a shared cube over-states gas-cap OGIP when the gas leg's
/// connate water is lower than the oil leg's). The oil leg always uses the cube.
/// This is still one scalar, not a saturation-height model — full transition-zone
/// / PVT saturation stays in petekSim.
///
/// # Errors
/// [`StaticError::OutOfRange`] if a contact depth is not finite;
/// [`StaticError::InvalidInput`] if `goc_m > owc_m` (the gas cap must sit above
/// the oil-water contact), `sw_gas` is not a fraction, a required property is
/// missing, or a value is not a fraction in `[0, 1]`.
pub fn compute_in_place_two_contact(
    grid: &Grid,
    goc_m: f64,
    owc_m: f64,
    sw_gas: Option<f64>,
) -> Result<InPlace, StaticError> {
    let nk = grid.dims().nk;
    compute_clipped(
        &GridSource::new(grid),
        Clip::Two {
            goc: goc_m,
            owc: owc_m,
            sw_gas,
        },
        0..nk,
        true,
    )
}

/// Summary-only variant of [`compute_in_place_two_contact`]: the same per-zone
/// totals without the per-cell HCPV cube (V7). `per_cell_hcpv` is left empty.
///
/// # Errors
/// Same as [`compute_in_place_two_contact`].
pub fn compute_in_place_two_contact_summary(
    grid: &Grid,
    goc_m: f64,
    owc_m: f64,
    sw_gas: Option<f64>,
) -> Result<InPlace, StaticError> {
    let nk = grid.dims().nk;
    compute_clipped(
        &GridSource::new(grid),
        Clip::Two {
            goc: goc_m,
            owc: owc_m,
            sw_gas,
        },
        0..nk,
        false,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::{build_box, BoxSpec, Dims};
    use crate::volumetrics::populate::{populate_constant, ConstantPriors};
    use petektools::units::{ft3_to_rb, sm3_to_stb, FT_TO_M};

    fn populated(area: f64, h: f64, top: f64, dims: Dims, p: ConstantPriors) -> Grid {
        let mut g = build_box(BoxSpec {
            area_m2: area,
            gross_height_m: h,
            dims,
            top_depth_m: top,
            aspect_ratio: 1.0,
        })
        .unwrap();
        populate_constant(&mut g, p).unwrap();
        g
    }

    #[test]
    fn golden_ooip_full_column() {
        // 400_000 m², 50 m, phi=0.25, NTG=0.8, Sw=0.3, Boi=1.25, contact below base.
        let p = ConstantPriors {
            porosity: 0.25,
            net_to_gross: 0.8,
            water_saturation: 0.3,
        };
        let grid = populated(400_000.0, 50.0, 5000.0, Dims::new(10, 10, 5).unwrap(), p);
        let res = compute_in_place(&grid, 6000.0).unwrap();

        // Whole grid is above contact -> GRV == bulk volume.
        let bulk = 400_000.0 * 50.0;
        assert!((res.grv_m3 - bulk).abs() / bulk < 1e-9);

        // HCPV = bulk * NTG * phi * (1-Sw).
        let hcpv = bulk * 0.8 * 0.25 * 0.7;
        assert!((res.hcpv_m3 - hcpv).abs() / hcpv < 1e-9);

        // OOIP [Sm³] = HCPV(reservoir m³) / Boi(Rm³/Sm³).
        let ooip = hcpv / 1.25;
        assert!((res.ooip_sm3(OilFvf::new(1.25).unwrap()) - ooip).abs() / ooip < 1e-9);
    }

    #[test]
    fn si_ooip_converts_back_to_the_identical_imperial_stb() {
        // Parity proof (decision_si_units_standard): the SAME physical reservoir
        // computed in SI, then converted Sm³→STB, equals the OLD imperial path
        // (dimension the box in ft, ft³→rb, /Boi). Proves numerical equivalence
        // old-imperial == new-SI × factor to FP tolerance.
        let boi = OilFvf::new(1.25).unwrap();
        let (phi, ntg, sw) = (0.25, 0.8, 0.3);
        // SI box: 400_000 m² × 50 m, whole column above contact.
        let grid = populated(
            400_000.0,
            50.0,
            5000.0,
            Dims::new(10, 10, 5).unwrap(),
            ConstantPriors {
                porosity: phi,
                net_to_gross: ntg,
                water_saturation: sw,
            },
        );
        let ip = compute_in_place(&grid, 6000.0).unwrap();
        // New SI answer, reported in STB via the pure geometric factor.
        let si_stb = sm3_to_stb(ip.ooip_sm3(boi));

        // Old imperial path: the identical reservoir dimensioned in feet.
        let hcpv_m3 = 400_000.0 * 50.0 * ntg * phi * (1.0 - sw);
        let hcpv_ft3 = hcpv_m3 / (FT_TO_M * FT_TO_M * FT_TO_M); // m³ -> ft³
        let imperial_stb = ft3_to_rb(hcpv_ft3) / boi.value(); // rb -> STB / Boi

        // Numbers (this case): HCPV 2.8e6 m³, OOIP = 2.24e6 Sm³ ≈ 14_088_647 STB.
        assert!(
            (si_stb - imperial_stb).abs() / imperial_stb < 1e-9,
            "SI STB {si_stb} != imperial STB {imperial_stb}"
        );
        assert!(
            (ip.ooip_sm3(boi) - 2.24e6).abs() / 2.24e6 < 1e-9,
            "OOIP == 2.24 MSm³, got {}",
            ip.ooip_sm3(boi)
        );
    }

    #[test]
    fn contact_excludes_deeper_cells() {
        let p = ConstantPriors {
            porosity: 0.2,
            net_to_gross: 1.0,
            water_saturation: 0.2,
        };
        // 4 layers over 5000-5100 m; contact at 5050 -> top 2 layers in column.
        let grid = populated(100.0, 100.0, 5000.0, Dims::new(2, 2, 4).unwrap(), p);
        let res = compute_in_place(&grid, 5050.0).unwrap();
        assert_eq!(res.cells_in_column, 2 * 2 * 2);
        let full = compute_in_place(&grid, 9999.0).unwrap();
        assert!((res.grv_m3 - full.grv_m3 / 2.0).abs() / full.grv_m3 < 1e-9);
    }

    #[test]
    fn two_contact_splits_gas_cap_and_oil_leg() {
        // 400_000 m², 100 m over 5000-5100, 10 layers of 10 m (centroids
        // 5005..5095). GOC 5030 -> gas cells 5005/5015/5025 (3 layers); OWC 5070
        // -> oil cells 5035/5045/5055/5065 (4 layers); the rest is water.
        let p = ConstantPriors {
            porosity: 0.2,
            net_to_gross: 1.0,
            water_saturation: 0.2,
        };
        let grid = populated(400_000.0, 100.0, 5000.0, Dims::new(2, 2, 10).unwrap(), p);
        let res = compute_in_place_two_contact(&grid, 5030.0, 5070.0, None).unwrap();

        let gas = res.gas.expect("gas zone");
        let oil = res.oil.expect("oil zone");
        assert_eq!(gas.cells, 2 * 2 * 3, "gas cap = 3 layers");
        assert_eq!(oil.cells, 2 * 2 * 4, "oil leg = 4 layers");

        let per_layer = 400_000.0 * 10.0;
        assert!((gas.grv_m3 - per_layer * 3.0).abs() / (per_layer * 3.0) < 1e-9);
        assert!((oil.grv_m3 - per_layer * 4.0).abs() / (per_layer * 4.0) < 1e-9);

        // Totals are gas + oil; HCPV ratio matches the 3:4 cell split (uniform).
        assert!((res.hcpv_m3 - (gas.hcpv_m3 + oil.hcpv_m3)).abs() < 1e-6);
        assert!((gas.hcpv_m3 / oil.hcpv_m3 - 3.0 / 4.0).abs() < 1e-9);

        // Per-zone surface volumes come off the right zone with the right FVF.
        let bgi = GasFvf::new(0.004).unwrap();
        let boi = OilFvf::new(1.25).unwrap();
        assert!(
            (res.gas_zone_ogip_sm3(bgi) - gas.hcpv_m3 / 0.004).abs() / (gas.hcpv_m3 / 0.004) < 1e-9
        );
        let expect_ooip = oil.hcpv_m3 / 1.25;
        assert!((res.oil_zone_ooip_sm3(boi) - expect_ooip).abs() / expect_ooip < 1e-9);

        // Cross-check: the total column equals a single contact clipped at OWC.
        let single = compute_in_place(&grid, 5070.0).unwrap();
        assert!((res.hcpv_m3 - single.hcpv_m3).abs() / single.hcpv_m3 < 1e-9);
        assert!(single.gas.is_none() && single.oil.is_none());
    }

    #[test]
    fn two_contact_rejects_goc_below_owc() {
        let p = ConstantPriors {
            porosity: 0.2,
            net_to_gross: 1.0,
            water_saturation: 0.2,
        };
        let grid = populated(100.0, 100.0, 5000.0, Dims::new(2, 2, 4).unwrap(), p);
        assert!(compute_in_place_two_contact(&grid, 5080.0, 5030.0, None).is_err());
        assert!(compute_in_place_two_contact(&grid, f64::NAN, 5030.0, None).is_err());
        // A non-fraction sw_gas is rejected (H2).
        assert!(compute_in_place_two_contact(&grid, 5030.0, 5080.0, Some(1.5)).is_err());
    }

    #[test]
    fn sw_gas_override_lowers_gas_cap_hcpv() {
        // Cube Sw = 0.4 everywhere; a lower gas-cap connate water 0.1 must raise
        // gas HCPV in the (1-Sw) ratio (0.9/0.6) and leave the oil leg untouched.
        let p = ConstantPriors {
            porosity: 0.2,
            net_to_gross: 1.0,
            water_saturation: 0.4,
        };
        let grid = populated(100.0, 100.0, 5000.0, Dims::new(2, 2, 10).unwrap(), p);
        let base = compute_in_place_two_contact(&grid, 5030.0, 5070.0, None).unwrap();
        let over = compute_in_place_two_contact(&grid, 5030.0, 5070.0, Some(0.1)).unwrap();
        let (gb, go) = (base.gas.unwrap(), over.gas.unwrap());
        // Oil leg identical; gas HCPV scales by (1-0.1)/(1-0.4) = 1.5.
        assert!((base.oil.unwrap().hcpv_m3 - over.oil.unwrap().hcpv_m3).abs() < 1e-6);
        assert!(
            (go.hcpv_m3 / gb.hcpv_m3 - (0.9 / 0.6)).abs() < 1e-9,
            "gas HCPV ratio {} != 1.5",
            go.hcpv_m3 / gb.hcpv_m3
        );
    }

    #[test]
    fn missing_property_errors() {
        let grid = build_box(BoxSpec::square(40.0, 20.0, Dims::new(2, 2, 2).unwrap())).unwrap();
        assert!(compute_in_place(&grid, 100.0).is_err());
    }

    #[test]
    fn ogip_from_gas_fvf() {
        let p = ConstantPriors {
            porosity: 0.2,
            net_to_gross: 0.9,
            water_saturation: 0.25,
        };
        let grid = populated(50.0, 40.0, 8000.0, Dims::new(5, 5, 4).unwrap(), p);
        let res = compute_in_place(&grid, 9000.0).unwrap();
        let expected = res.hcpv_m3 / 0.004;
        assert!((res.ogip_sm3(GasFvf::new(0.004).unwrap()) - expected).abs() / expected < 1e-9);
    }
}
