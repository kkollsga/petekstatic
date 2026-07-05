//! Zones — the stratigraphic index of the model (SPEC §4). A [`Zone`] is a named
//! vertical interval bounded by two framework horizons, spanning a contiguous
//! `k`-range, so the grid is addressable by geology instead of raw cell index.
//!
//! MVP: a **single implicit zone** (the whole column, Top→contact). The seam is
//! defined now; the multi-zone build (P5 `task_petekstatic_zones_faults`) slots in
//! without changing the consumer contract.

use std::ops::Range;

/// A named stratigraphic interval spanning a contiguous range of `k`-layers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Zone {
    pub name: String,
    /// Optional display colour (a viewer hint, e.g. `"#ffcc00"`), carried from the
    /// source [`crate::StackZone::color`] on a horizon-stack build. `None` on the
    /// single-implicit-zone (from-wireframe / flat) paths. Surfaced in the section
    /// bundle's per-zone `zones` list for colour-by-zone rendering.
    pub color: Option<String>,
    pub top_horizon: String,
    pub base_horizon: String,
    pub k_range: Range<usize>,
}

/// The ordered top→base zone table; `k`-ranges partition `[0, nk)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ZoneTable {
    zones: Vec<Zone>,
}

impl ZoneTable {
    /// The MVP degenerate whole-column zone (Top→contact over all `nk` layers).
    #[must_use]
    pub fn single(nk: usize) -> Self {
        Self {
            zones: vec![Zone {
                name: "RESERVOIR".to_string(),
                color: None,
                top_horizon: "Top".to_string(),
                base_horizon: "Contact".to_string(),
                k_range: 0..nk,
            }],
        }
    }

    /// Build the real multi-zone table from an ordered horizon stack: `N`
    /// horizon names (top→down) and the per-zone layer counts (`N - 1` of them),
    /// producing `N - 1` zones whose `k`-ranges partition `[0, sum(nk))`
    /// top→base. Zone `z` is bounded by horizons `z` (top) and `z + 1` (base) and
    /// carries the caller-supplied zone name.
    ///
    /// The builder constructs this from the resolved framework stack; a consumer
    /// then addresses the grid by geology (`get(name)` → its `k`-range).
    #[must_use]
    pub fn from_stack(
        horizon_names: &[String],
        zone_names: &[String],
        zone_colors: &[Option<String>],
        per_zone_nk: &[usize],
    ) -> Self {
        debug_assert_eq!(zone_names.len(), per_zone_nk.len());
        debug_assert_eq!(zone_names.len(), zone_colors.len());
        debug_assert_eq!(horizon_names.len(), zone_names.len() + 1);
        let mut zones = Vec::with_capacity(zone_names.len());
        let mut k_start = 0usize;
        for (z, nk) in per_zone_nk.iter().copied().enumerate() {
            zones.push(Zone {
                name: zone_names[z].clone(),
                color: zone_colors[z].clone(),
                top_horizon: horizon_names[z].clone(),
                base_horizon: horizon_names[z + 1].clone(),
                k_range: k_start..k_start + nk,
            });
            k_start += nk;
        }
        Self { zones }
    }

    /// The zones, ordered top→base.
    #[must_use]
    pub fn zones(&self) -> &[Zone] {
        &self.zones
    }

    /// A zone by name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Zone> {
        self.zones.iter().find(|z| z.name == name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_zone_partitions_the_column() {
        let zt = ZoneTable::single(8);
        assert_eq!(zt.zones().len(), 1);
        assert_eq!(zt.zones()[0].k_range, 0..8);
        assert!(zt.get("RESERVOIR").is_some());
        assert!(zt.get("MISSING").is_none());
    }

    #[test]
    fn stack_zones_partition_the_column_top_to_base() {
        let horizons: Vec<String> = ["H0", "H1", "H2", "H3"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let names: Vec<String> = ["ZA", "ZB", "ZC"].iter().map(|s| s.to_string()).collect();
        let colors = vec![None, Some("#0af".to_string()), None];
        let zt = ZoneTable::from_stack(&horizons, &names, &colors, &[3, 2, 4]);
        assert_eq!(zt.zones().len(), 3);
        assert_eq!(zt.zones()[0].k_range, 0..3);
        assert_eq!(zt.zones()[1].k_range, 3..5);
        assert_eq!(zt.zones()[2].k_range, 5..9);
        // Bounding horizons + names line up.
        assert_eq!(zt.get("ZB").unwrap().top_horizon, "H1");
        assert_eq!(zt.get("ZB").unwrap().base_horizon, "H2");
        // The ranges tile [0, 9) with no gap/overlap.
        let covered: usize = zt.zones().iter().map(|z| z.k_range.len()).sum();
        assert_eq!(covered, 9);
    }
}
