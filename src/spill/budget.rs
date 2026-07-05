//! Memory budget + build-mode decision for the out-of-core pipeline (rulings
//! **R2/R4/R5**, `petekSuite/dev-docs/designs/out-of-core-strategy.md`).
//!
//! The engine chooses its backing-storage mode **loudly** (R5): a declared
//! [`MemoryBudget`] (default: a documented fraction of physical RAM) is compared
//! against a live-set estimate; below it the pipeline stays on today's pure
//! in-core path (byte-identical), above it it switches to the k-slab-streaming
//! spill (geometry + cubes onto a petekTools store, f32 at spill scale) and
//! **says so** ([`SpillNotice::warn`]). Never an OOM kill; never a silent switch.

use crate::grid::Dims;
use std::path::PathBuf;
use std::sync::OnceLock;

/// Default budget fraction of physical RAM when the caller does not override it
/// (R5). Conservative: half of RAM leaves headroom for the OS, the warm-start
/// state, and the reusable MC model (the live set roughly doubles under MC).
pub const DEFAULT_BUDGET_FRACTION: f64 = 0.5;

/// Bytes per in-core geometry/cube element (`f64`). ZCORN is `8` of these per
/// cell; each property cube is `1` per cell.
const ELEM_BYTES: u64 = 8;

/// The warm-start template + reusable MC model roughly **double** the bare
/// grid's live set (`out-of-core-strategy.md` sizing table). The estimate folds
/// this in so a model that only *just* fits a bare in-core grid but blows the
/// budget under MC still spills.
const WARM_FACTOR: f64 = 2.0;

/// Fallback physical-RAM assumption when the platform query fails — deliberately
/// modest (8 GiB) so an unknown machine errs toward spilling a genuinely large
/// model rather than risking an OOM.
const FALLBACK_PHYSICAL_BYTES: u64 = 8 * 1024 * 1024 * 1024;

/// Physical RAM in bytes, queried once per process and cached. `None` when the
/// platform is unsupported or the query fails (callers fall back to
/// [`FALLBACK_PHYSICAL_BYTES`]).
#[must_use]
pub fn physical_ram_bytes() -> Option<u64> {
    static CACHE: OnceLock<Option<u64>> = OnceLock::new();
    *CACHE.get_or_init(query_physical_ram)
}

#[cfg(target_os = "linux")]
fn query_physical_ram() -> Option<u64> {
    let txt = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in txt.lines() {
        if let Some(rest) = line.strip_prefix("MemTotal:") {
            let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

#[cfg(target_os = "macos")]
fn query_physical_ram() -> Option<u64> {
    // `sysctl -n hw.memsize` → total physical bytes. One process spawn per
    // process (cached), off any hot path.
    let out = std::process::Command::new("sysctl")
        .arg("-n")
        .arg("hw.memsize")
        .output()
        .ok()?;
    String::from_utf8(out.stdout).ok()?.trim().parse().ok()
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
fn query_physical_ram() -> Option<u64> {
    None
}

/// A declared memory budget for the build/template API (R5). Below it the
/// pipeline stays in-core; above it it spills.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MemoryBudget {
    limit_bytes: u64,
}

impl MemoryBudget {
    /// An explicit byte budget.
    #[must_use]
    pub fn bytes(limit_bytes: u64) -> Self {
        Self { limit_bytes }
    }

    /// A fraction of physical RAM (clamped to `[0, 1]`); uses
    /// [`FALLBACK_PHYSICAL_BYTES`] when the RAM query fails.
    #[must_use]
    pub fn fraction_of_physical(fraction: f64) -> Self {
        let phys = physical_ram_bytes().unwrap_or(FALLBACK_PHYSICAL_BYTES);
        let frac = fraction.clamp(0.0, 1.0);
        Self {
            limit_bytes: (phys as f64 * frac) as u64,
        }
    }

    /// An effectively unlimited budget — the pipeline never spills (the explicit
    /// opt-out; also what tests use to force the in-core path).
    #[must_use]
    pub fn unlimited() -> Self {
        Self {
            limit_bytes: u64::MAX,
        }
    }

    /// The budget in bytes.
    #[must_use]
    pub fn limit_bytes(&self) -> u64 {
        self.limit_bytes
    }

    /// Whether the budget is effectively unlimited.
    #[must_use]
    pub fn is_unlimited(&self) -> bool {
        self.limit_bytes == u64::MAX
    }
}

impl Default for MemoryBudget {
    /// The documented default (R5): [`DEFAULT_BUDGET_FRACTION`] of physical RAM.
    fn default() -> Self {
        Self::fraction_of_physical(DEFAULT_BUDGET_FRACTION)
    }
}

/// The backing-storage mode the engine picks for a build (R5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BuildMode {
    /// Today's pure in-core path — geometry + cubes are owned `Vec`s. Byte-identical
    /// to the pre-out-of-core behaviour.
    InCore,
    /// The out-of-core path — geometry + cubes are streamed to a petekTools store
    /// (f32 lanes at spill scale, ruling R4) and the model reads through windowed
    /// mmap views.
    Spilled,
}

/// The in-core live-set estimate for a grid of `dims` with `n_cubes` property
/// cubes, in bytes (ZCORN `f64` + cubes `f64`, scaled by [`WARM_FACTOR`]). COORD
/// is `O(area)` and negligible against the per-cell arrays, so it is omitted.
#[must_use]
pub fn live_set_bytes(dims: Dims, n_cubes: usize) -> u64 {
    let cells = dims.cell_count() as u64;
    let zcorn = cells.saturating_mul(8).saturating_mul(ELEM_BYTES);
    let cubes = (n_cubes as u64)
        .saturating_mul(cells)
        .saturating_mul(ELEM_BYTES);
    (((zcorn.saturating_add(cubes)) as f64) * WARM_FACTOR) as u64
}

/// Decide the build mode for `dims` / `n_cubes` against `budget`, returning the
/// mode and the live-set estimate that drove it (the estimate is surfaced in the
/// loud [`SpillNotice`], R5).
#[must_use]
pub fn decide_mode(dims: Dims, n_cubes: usize, budget: MemoryBudget) -> (BuildMode, u64) {
    let estimate = live_set_bytes(dims, n_cubes);
    let mode = if estimate > budget.limit_bytes() {
        BuildMode::Spilled
    } else {
        BuildMode::InCore
    };
    (mode, estimate)
}

/// The loud mode-switch advisory (R5): naming the mode, the budget, the computed
/// live-set estimate, and the store path. Emitted to stderr by [`SpillNotice::warn`]
/// on every spilled build — never a silent switch.
#[derive(Debug, Clone)]
pub struct SpillNotice {
    /// Total cell count of the spilled model.
    pub cells: usize,
    /// The declared budget that was exceeded, bytes.
    pub budget_bytes: u64,
    /// The live-set estimate that drove the switch, bytes.
    pub estimate_bytes: u64,
    /// Where the spill store was written.
    pub store_path: PathBuf,
}

impl SpillNotice {
    /// Emit the advisory to stderr.
    pub fn warn(&self) {
        eprintln!(
            "petekstatic: OUT-OF-CORE mode — live-set estimate {} MiB exceeds budget {} MiB \
             ({} cells); spilling geometry + cubes (f32) to {}",
            self.estimate_bytes >> 20,
            self.budget_bytes >> 20,
            self.cells,
            self.store_path.display(),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unlimited_never_spills() {
        let dims = Dims::new(1000, 1000, 100).unwrap(); // 100M cells
        let (mode, est) = decide_mode(dims, 3, MemoryBudget::unlimited());
        assert_eq!(mode, BuildMode::InCore);
        assert!(est > 0);
    }

    #[test]
    fn tiny_budget_forces_spill() {
        let dims = Dims::new(100, 100, 100).unwrap(); // 1M cells
        let (mode, est) = decide_mode(dims, 3, MemoryBudget::bytes(1024));
        assert_eq!(mode, BuildMode::Spilled);
        // 1M cells: ZCORN 8*8 = 64 MB + 3 cubes 24 MB, doubled ≈ 176 MB.
        assert!(est > 100 * 1024 * 1024, "estimate {est} looks too small");
    }

    #[test]
    fn generous_budget_stays_in_core() {
        let dims = Dims::new(50, 50, 10).unwrap(); // 25k cells — kilobytes
        let (mode, _) = decide_mode(dims, 3, MemoryBudget::bytes(1024 * 1024 * 1024));
        assert_eq!(mode, BuildMode::InCore);
    }

    #[test]
    fn fraction_is_clamped_and_positive() {
        assert!(MemoryBudget::fraction_of_physical(0.5).limit_bytes() > 0);
        // Over-unit fraction clamps to the whole (never panics / overflows).
        let whole = MemoryBudget::fraction_of_physical(2.0).limit_bytes();
        let phys = physical_ram_bytes().unwrap_or(FALLBACK_PHYSICAL_BYTES);
        assert_eq!(whole, phys);
    }
}
