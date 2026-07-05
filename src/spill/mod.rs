//! `srs-spill` — the out-of-core **backing-storage mode** for a `StaticModel`
//! (rulings **R1/R2/R4/R5**, `petekSuite/dev-docs/designs/out-of-core-strategy.md`).
//!
//! Split out of `srs-model` (organize wave P10): a leaf crate below the model
//! aggregate that owns two responsibilities, both self-contained over `srs-grid`
//! geometry + `srs-volumetrics` streaming:
//!
//! - [`budget`] — the loud memory-budget / build-mode decision ([`MemoryBudget`],
//!   [`decide_mode`], [`SpillNotice`]): below the budget the pipeline stays on the
//!   byte-identical in-core path; above it it spills and **says so** (R5).
//! - [`spill`] — the k-slab spill itself: geometry (ZCORN) + property cubes stream
//!   onto a petekTools `store` (f32 lanes at spill scale, R4) and read back through
//!   windowed mmap views ([`SpillBacking`] / [`SpillSource`]). Peak working set is
//!   `O(slab)`, not `O(grid)` (R2).
//!
//! `srs-model` re-exports this crate's surface, so the public path
//! (`crate::model::MemoryBudget`, `crate::model::spill_grid`, …) is unchanged by the split.

mod budget;
#[allow(clippy::module_inception)]
mod spill;

pub use budget::{
    decide_mode, live_set_bytes, physical_ram_bytes, BuildMode, MemoryBudget, SpillNotice,
    DEFAULT_BUDGET_FRACTION,
};
pub use spill::{
    spill_grid, spill_grid_to, spill_streaming, unique_spill_path, SpillBacking, SpillSlab,
    SpillSource,
};
