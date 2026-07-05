//! Stratigraphic layering (`layer_interpolation_spec`, high): subdivide the
//! volume between a top and base surface into conformable k-layers, writing the
//! corner-point ZCORN so the grid honours structure.
//!
//! ## Conformity styles
//! - [`Conformity::Proportional`] — each layer takes an equal fraction of the
//!   local top→base thickness; `nk` is honoured as given (the historical default).
//! - [`Conformity::FollowTop`] — each layer surface is draped **parallel to the
//!   top** at a constant `dz_m` (layer `k` base = `top + (k+1)·dz`). Where a
//!   column is thinner than the full layer stack, the deep layers **truncate
//!   against the base** (erosional/truncation geometry).
//! - [`Conformity::FollowBase`] — the mirror: layer surfaces parallel to the
//!   **base**, pinching (onlapping) against the top.
//!
//! ## Truncation = zero-thickness collapse (the chosen representation)
//! A truncated (inactive) cell is represented **geometrically**, by collapsing
//! its ZCORN onto the pinch-out horizon so the cell has **zero bulk volume** —
//! exactly how a Proportional pinch-out already behaves. There is no parallel
//! `active` mask array: a zero-volume cell is already respected everywhere
//! downstream — volumetrics excludes it (its bulk/pore volume is 0, so GRV/HCPV
//! totals are conformity-invariant), property population writes harmless finite
//! values into it (never NaN into an *active* cell), and the view bundles mark a
//! collapsed cell with `NaN` on a pure-geometric `dz ≤ ε` test so their schemas
//! stay `nk`-sized and stable. [`LayeredGrid::truncated_cells`] reports the count
//! for the caller's provenance (informational, never an error).
//!
//! ## dz-derived nk (Follow styles)
//! Under a Follow style `nk` is **derived** from geometry — `ceil(max column
//! thickness / dz)` — so the thickest column is fully layered and thinner columns
//! truncate. The passed `nk` argument is ignored for Follow styles (honoured only
//! by `Proportional`). The derived count is capped at [`MAX_NK`] (a finer dz than
//! that allows over the thickest column sets [`LayeredGrid::nk_capped`] and lets
//! the deepest part of the thickest columns truncate against the cap rather than
//! exploding the layer count).

use crate::error::StaticError;
use crate::grid::{CornerPointGeom, Dims, Grid, Pillar, Point3};
use crate::gridder::surface::Surface;
use serde::{Deserialize, Serialize};

/// Reusable per-pillar interface-depth scratch for [`layer_grid_stack_into`]: the
/// `raw` (conformity result) and `snap` (post-collapse) columns, each
/// `pillar_count × (nk + 1)` f64. Owned by the caller (one per MC worker) and
/// refilled each draw, so the layering scratch is allocated once, not per
/// realization (`StaticModelTemplate::realize_into`).
#[derive(Debug, Default)]
pub struct LayerScratch {
    raw: Vec<f64>,
    snap: Vec<f64>,
}

impl LayerScratch {
    /// An empty scratch; it grows to fit on the first [`layer_grid_stack_into`].
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

// Scratch is a reusable buffer, not logical state: a clone starts **empty** (each
// owner grows its own on first use), so cloning a template one-per-worker never
// copies megabytes of scratch.
impl Clone for LayerScratch {
    fn clone(&self) -> Self {
        Self::default()
    }
}

/// The layering result of [`layer_grid_stack_into`]: the grid dimensions plus the
/// same per-zone report as [`StackedLayeredGrid`], but **without** a built [`Grid`]
/// — the caller installed the refilled geometry into a recycled grid.
#[derive(Debug, Clone)]
pub struct StackLayering {
    /// Grid dimensions (`nk` = total layer count).
    pub dims: Dims,
    /// Total k-layer count = sum of the per-zone counts (after the total cap).
    pub nk: usize,
    /// Per-zone layering, top→down; k-ranges partition `[0, nk)`.
    pub zones: Vec<StackedZone>,
    /// Total cells collapsed to zero thickness by truncation across all zones.
    pub truncated_cells: usize,
    /// Total cells collapsed to zero thickness by the cell-collapse pass.
    pub collapsed_cells: usize,
    /// `true` if the summed per-zone count exceeded [`MAX_NK`] and the zones were
    /// scaled down to fit.
    pub nk_capped: bool,
}

/// Upper bound on the dz-derived k-layer count under a Follow style. Matches the
/// fine-layering (~1 m default dz) work: a dz that would need more than this many
/// layers to span the thickest column is capped here.
pub const MAX_NK: usize = 200;

/// A cell whose collapsed thickness is at or below this (metres) counts as
/// truncated/inactive (zero volume).
const TRUNC_EPS: f64 = 1e-9;

/// How layer boundaries are distributed between top and base.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Conformity {
    /// Proportional: each layer takes an equal fraction of local thickness; `nk`
    /// is honoured as passed.
    Proportional,
    /// Follow-top (truncation/erosional): layer surfaces drape parallel to the
    /// top at constant `dz_m` (layer `k` base = `top + (k+1)·dz`); deep layers
    /// truncate against the base where the column pinches. `nk` is dz-derived.
    FollowTop { dz_m: f64 },
    /// Follow-base (onlap): the mirror — layer surfaces parallel to the base,
    /// pinching against the top. `nk` is dz-derived.
    FollowBase { dz_m: f64 },
}

impl Conformity {
    /// The constant layer thickness of a Follow style, if any.
    fn dz_m(self) -> Option<f64> {
        match self {
            Conformity::Proportional => None,
            Conformity::FollowTop { dz_m } | Conformity::FollowBase { dz_m } => Some(dz_m),
        }
    }
}

/// A layered corner-point grid plus the layering report the caller stamps into
/// provenance: the effective (possibly dz-derived) `nk`, the count of cells
/// collapsed to zero thickness by truncation, and whether the derived count hit
/// [`MAX_NK`].
#[derive(Debug, Clone)]
pub struct LayeredGrid {
    /// The built grid; its `dims().nk` equals [`LayeredGrid::nk`].
    pub grid: Grid,
    /// Effective k-layer count — dz-derived for a Follow style, the passed `nk`
    /// for `Proportional`.
    pub nk: usize,
    /// Cells collapsed to zero thickness by truncation against the pinch-out
    /// horizon (informational; conformity conserves total volume regardless).
    pub truncated_cells: usize,
    /// `true` if the dz-derived `nk` was clamped to [`MAX_NK`] (dz finer than the
    /// thickest column can carry within the cap; its deepest part truncates).
    pub nk_capped: bool,
}

/// Depth of layer-boundary level `level` (0..=nk) for one pillar.
fn boundary_depth(conf: Conformity, zt: f64, zb: f64, level: usize, nk: usize) -> f64 {
    let (lo, hi) = (zt.min(zb), zt.max(zb));
    let z = match conf {
        Conformity::Proportional => zt + (level as f64 / nk as f64) * (zb - zt),
        // Parallel to the base: constant dz measured up from the base; the top
        // (level 0) side onlaps.
        Conformity::FollowBase { dz_m } => zb - (nk - level) as f64 * dz_m,
        // Parallel to the top: constant dz measured down from the top; the base
        // (level nk) side truncates.
        Conformity::FollowTop { dz_m } => zt + level as f64 * dz_m,
    };
    z.clamp(lo, hi)
}

/// Derive the effective k-layer count for `conformity`: the passed `nk` for
/// `Proportional`, or `ceil(max column thickness / dz)` (≥ 1, capped at
/// [`MAX_NK`]) for a Follow style. Returns `(nk, capped)`.
fn derive_nk(
    top: &Surface,
    base: &Surface,
    conformity: Conformity,
    requested_nk: usize,
) -> (usize, bool) {
    let Some(dz) = conformity.dz_m() else {
        return (requested_nk, false);
    };
    let mut max_t = 0.0_f64;
    for jp in 0..top.ny() {
        for ip in 0..top.nx() {
            max_t = max_t.max((base.z(ip, jp) - top.z(ip, jp)).abs());
        }
    }
    let n = (max_t / dz).ceil() as usize;
    let n = n.max(1);
    if n > MAX_NK {
        (MAX_NK, true)
    } else {
        (n, false)
    }
}

/// Build a k-layered corner-point grid between `top` and `base` surfaces sharing
/// a `(ni+1) x (nj+1)` node lattice, with uniform areal spacing `dx, dy`.
///
/// `nk` is honoured for [`Conformity::Proportional`] and **ignored** for the
/// Follow styles (there it is dz-derived — see [`LayeredGrid::nk`]). The return
/// carries the layering report (effective `nk`, truncated-cell count, cap flag).
///
/// This is the **2-surface degenerate case** of [`layer_grid_stack`] (one zone
/// between one top and one base); the multi-zone stack path reuses the same
/// per-cell ZCORN machinery.
///
/// # Errors
/// [`StaticError::InvalidInput`] if the surfaces' lattices differ, `nk == 0` (for
/// `Proportional`), a Follow `dz_m` is not finite and positive, or `dx`/`dy` are
/// not positive; [`StaticError::Grid`] on a degenerate lattice.
pub fn layer_grid(
    top: &Surface,
    base: &Surface,
    dx: f64,
    dy: f64,
    nk: usize,
    conformity: Conformity,
) -> Result<LayeredGrid, StaticError> {
    let stacked = layer_grid_stack(
        &[top, base],
        dx,
        dy,
        &[ZoneLayerSpec {
            conformity,
            requested_nk: nk,
        }],
        None,
    )?;
    Ok(LayeredGrid {
        grid: stacked.grid,
        nk: stacked.nk,
        truncated_cells: stacked.truncated_cells,
        nk_capped: stacked.nk_capped,
    })
}

/// The layering scheme for one zone of a horizon stack: its [`Conformity`] plus
/// the requested layer count (honoured only by [`Conformity::Proportional`]; a
/// Follow style derives `nk` from the zone's own `dz` and thickness).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZoneLayerSpec {
    /// This zone's conformity/layering style.
    pub conformity: Conformity,
    /// Requested layer count — used only for `Proportional` (Follow derives it).
    pub requested_nk: usize,
}

/// The per-zone slice of a stacked layering: this zone's effective layer count,
/// its first global `k`, its own conformity, and the cells it truncated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StackedZone {
    /// Effective k-layer count for this zone (dz-derived for a Follow style,
    /// the requested count for `Proportional`, after the total cap is applied).
    pub nk: usize,
    /// First global k-layer this zone owns (its zones partition `[0, total_nk)`).
    pub k_start: usize,
    /// This zone's conformity/layering style.
    pub conformity: Conformity,
    /// Cells in this zone collapsed to zero thickness by truncation.
    pub truncated_cells: usize,
    /// Cells in this zone collapsed to zero thickness by the sub-threshold
    /// **cell-collapse** pass (their sliver thickness merged into a thicker
    /// zone-interior neighbour; volume-conserving). `0` when collapse is off.
    pub collapsed_cells: usize,
}

impl StackedZone {
    /// The global k-range `[k_start, k_start + nk)` this zone owns.
    #[must_use]
    pub fn k_range(&self) -> core::ops::Range<usize> {
        self.k_start..self.k_start + self.nk
    }
}

/// A stacked multi-zone corner-point grid plus its per-zone layering report:
/// N ordered surfaces (top→down) define N−1 zones, each layered by its own
/// [`ZoneLayerSpec`]; the grid is their vertical concatenation in k.
#[derive(Debug, Clone)]
pub struct StackedLayeredGrid {
    /// The built grid; `dims().nk` equals [`StackedLayeredGrid::nk`] (the total).
    pub grid: Grid,
    /// Total k-layer count = sum of the per-zone counts (after the total cap).
    pub nk: usize,
    /// Per-zone layering, top→down; k-ranges partition `[0, nk)`.
    pub zones: Vec<StackedZone>,
    /// Total cells collapsed to zero thickness by truncation across all zones.
    pub truncated_cells: usize,
    /// Total cells collapsed to zero thickness by the cell-collapse pass across all
    /// zones (volume merged into a neighbour); `0` when collapse is off.
    pub collapsed_cells: usize,
    /// `true` if the summed per-zone count exceeded [`MAX_NK`] and the zones were
    /// scaled down to fit (the per-zone breakdown is in [`StackedLayeredGrid::zones`]).
    pub nk_capped: bool,
}

/// Build a stacked multi-zone corner-point grid from `surfaces` — N ordered depth
/// surfaces (top→down, sharing one `(ni+1) x (nj+1)` node lattice) — with one
/// [`ZoneLayerSpec`] per consecutive-surface interval (`zone_specs.len() ==
/// surfaces.len() - 1`). Each zone is layered by its own conformity between its
/// bounding surfaces; the resulting sub-grids are concatenated in k, so the total
/// layer count is the sum of the per-zone counts.
///
/// The surfaces are the resolved framework horizons (a tops-only horizon has
/// already been draped by the caller); ordering between consecutive surfaces is
/// the caller's contract (the builder's per-interface order-repair enforces it).
/// Where a surface pair crosses, the corresponding cells simply truncate to zero
/// volume (as a single-zone Follow pinch-out already does).
///
/// ## Total layer cap
/// The [`MAX_NK`] cap applies to the **total** (`sum` of per-zone counts). If the
/// sum exceeds it the per-zone counts are scaled down proportionally (each kept
/// `>= 1`) and [`StackedLayeredGrid::nk_capped`] is set; the per-zone breakdown
/// stays in [`StackedLayeredGrid::zones`] for the caller's warning.
///
/// ## Cell-collapse threshold (`collapse_below_m`)
/// With `Some(threshold)`, after layering any cell thinner than `threshold`
/// **collapses**: the layer interface is snapped so a thicker zone-interior
/// neighbour absorbs the sliver's thickness — **volume-conserving** (rock is
/// merged, never deleted). The merge is into the thicker vertical neighbour and
/// **never crosses a zone boundary** (a zone's bounding interfaces are pinned), so
/// zone volumes are individually conserved. Collapse is applied per pillar, so a
/// cell zeroes only where all four pillars collapse it (the same max-pillar rule
/// as truncation); per-zone counts land in [`StackedZone::collapsed_cells`]. `None`
/// = off (the default).
///
/// # Errors
/// [`StaticError::InvalidInput`] if fewer than 2 surfaces are given, the spec
/// count is not `surfaces.len() - 1`, the lattices differ, a `Proportional` zone
/// requests `nk == 0`, a Follow `dz_m` is not finite-positive, `dx`/`dy` are not
/// positive, `collapse_below_m` is not finite-positive, or the zone count exceeds
/// [`MAX_NK`]; [`StaticError::Grid`] on a degenerate lattice.
pub fn layer_grid_stack(
    surfaces: &[&Surface],
    dx: f64,
    dy: f64,
    zone_specs: &[ZoneLayerSpec],
    collapse_below_m: Option<f64>,
) -> Result<StackedLayeredGrid, StaticError> {
    // The allocating convenience path: fresh scratch + geometry buffers, filled by
    // the recycling core, then wrapped in a Grid. One code path with `realize_into`.
    let mut scratch = LayerScratch::new();
    let mut coord = Vec::new();
    let mut zcorn = Vec::new();
    let lay = layer_grid_stack_into(
        surfaces,
        dx,
        dy,
        zone_specs,
        collapse_below_m,
        &mut scratch,
        &mut coord,
        &mut zcorn,
    )?;
    Ok(StackedLayeredGrid {
        grid: Grid::new(CornerPointGeom::new(lay.dims, coord, zcorn)),
        nk: lay.nk,
        zones: lay.zones,
        truncated_cells: lay.truncated_cells,
        collapsed_cells: lay.collapsed_cells,
        nk_capped: lay.nk_capped,
    })
}

/// [`layer_grid_stack`] writing its geometry into **caller-owned** buffers instead
/// of freshly allocating — the allocation-recycling core of
/// `StaticModelTemplate::realize_into`. `scratch` (the per-pillar interface
/// columns), `coord` (the pillar lattice), and `zcorn` (the per-cell corner
/// depths, the 64 MB/1M-cell dominant) are **cleared and fully overwritten**, so
/// carrying stale buffers from a prior draw is safe and the result is
/// bit-identical to the allocating path. Returns the [`StackLayering`] report; the
/// caller installs `coord`/`zcorn` into a recycled grid
/// ([`crate::grid::Grid::install_geometry`]).
///
/// # Errors
/// Identical to [`layer_grid_stack`].
#[allow(clippy::too_many_arguments)]
pub fn layer_grid_stack_into(
    surfaces: &[&Surface],
    dx: f64,
    dy: f64,
    zone_specs: &[ZoneLayerSpec],
    collapse_below_m: Option<f64>,
    scratch: &mut LayerScratch,
    coord: &mut Vec<Pillar>,
    zcorn: &mut Vec<f64>,
) -> Result<StackLayering, StaticError> {
    if surfaces.len() < 2 {
        return Err(StaticError::InvalidInput(format!(
            "a horizon stack needs at least 2 surfaces, got {}",
            surfaces.len()
        )));
    }
    if zone_specs.len() != surfaces.len() - 1 {
        return Err(StaticError::InvalidInput(format!(
            "expected {} zone specs for {} surfaces, got {}",
            surfaces.len() - 1,
            surfaces.len(),
            zone_specs.len()
        )));
    }
    let nx = surfaces[0].nx();
    let ny = surfaces[0].ny();
    for s in surfaces {
        if s.nx() != nx || s.ny() != ny {
            return Err(StaticError::InvalidInput(
                "all stack surfaces must share one lattice".into(),
            ));
        }
    }
    require_finite_surfaces(surfaces)?;
    for spec in zone_specs {
        if let Some(dz) = spec.conformity.dz_m() {
            if !(dz.is_finite() && dz > 0.0) {
                return Err(StaticError::InvalidInput(format!(
                    "conformity dz_m must be finite and > 0, got {dz}"
                )));
            }
        } else if spec.requested_nk == 0 {
            return Err(StaticError::InvalidInput("nk must be >= 1".into()));
        }
    }
    if !(dx.is_finite() && dx > 0.0 && dy.is_finite() && dy > 0.0) {
        return Err(StaticError::InvalidInput(format!(
            "dx, dy must be finite and > 0, got {dx}, {dy}"
        )));
    }
    if zone_specs.len() > MAX_NK {
        return Err(StaticError::InvalidInput(format!(
            "{} zones exceed the {MAX_NK}-layer total cap (a zone needs >= 1 layer)",
            zone_specs.len()
        )));
    }
    if let Some(t) = collapse_below_m {
        if !(t.is_finite() && t > 0.0) {
            return Err(StaticError::InvalidInput(format!(
                "collapse_below_m must be finite and > 0, got {t}"
            )));
        }
    }

    // Per-zone effective layer counts (dz-derived for Follow, requested for
    // Proportional), then the TOTAL cap: scale the whole stack down if the sum
    // exceeds MAX_NK, keeping every zone >= 1 layer.
    let mut zone_capped = false;
    let per_zone: Vec<usize> = zone_specs
        .iter()
        .zip(surfaces.windows(2))
        .map(|(spec, pair)| {
            let (n, capped) = derive_nk(pair[0], pair[1], spec.conformity, spec.requested_nk);
            zone_capped |= capped;
            n
        })
        .collect();
    let (per_zone, total_capped) = cap_total_nk(&per_zone);
    let nk_capped = zone_capped || total_capped;

    // k -> owning zone map + per-zone k_start.
    let total_nk: usize = per_zone.iter().sum();
    let mut k_start = Vec::with_capacity(per_zone.len());
    let mut zone_of_k = Vec::with_capacity(total_nk);
    let mut acc = 0usize;
    for (z, &n) in per_zone.iter().enumerate() {
        k_start.push(acc);
        for _ in 0..n {
            zone_of_k.push(z);
        }
        acc += n;
    }

    let dims = Dims::new(nx - 1, ny - 1, total_nk)?;

    // Per-pillar interface depths: for each areal node, the `total_nk + 1` layer
    // boundary levels 0..=total_nk (level 0 = surfaces[0], level total_nk =
    // surfaces[last]). `raw` is the conformity result; `snap` is the same after the
    // volume-conserving cell-collapse pass (== raw when collapse is off). Corner
    // depths read straight off these arrays, so per-pillar collapse stays
    // watertight (adjacent cells share a pillar's interface column).
    let ncorner = total_nk + 1;
    let np = nx * ny;
    // Recycle the caller's scratch, sized to exactly `np * ncorner`. The loop below
    // writes EVERY element (each pillar's levels `0..=total_nk` in `raw`, and the
    // whole column via `copy_from_slice` in `snap`), so stale content carried from a
    // prior draw is fully overwritten before any read — no zero-init needed.
    // Deliberately NOT `resize(_, 0.0)`: that memsets ~np·ncorner f64 per draw
    // (~17 MB at 1M cells) which the fill then overwrites — pure waste that measured
    // as a ~10% serial regression against the allocating path (whose `alloc_zeroed`
    // got its zero pages from the OS). `size_scratch` is a no-op `truncate` on the
    // steady-state (constant-nk) path.
    // Move the reused buffers into owned locals for the hot fill loops. Writing
    // through the `&mut Vec` parameters forces the compiler to reload the backing
    // pointer/len each iteration (it cannot prove they don't alias the surface
    // slices), which measured as a real serial regression vs the allocating path's
    // local `Vec`s. Take → fill locals → write back: all O(1), capacity preserved,
    // codegen restored.
    let mut raw = std::mem::take(&mut scratch.raw);
    let mut snap = std::mem::take(&mut scratch.snap);
    let mut coord_buf = std::mem::take(coord);
    let mut zcorn_buf = std::mem::take(zcorn);
    size_scratch(&mut raw, np * ncorner);
    size_scratch(&mut snap, np * ncorner);
    for jp in 0..ny {
        for ip in 0..nx {
            let base = (jp * nx + ip) * ncorner;
            for (z, &znk) in per_zone.iter().enumerate() {
                let (zt, zb) = (surfaces[z].z(ip, jp), surfaces[z + 1].z(ip, jp));
                for local in 0..=znk {
                    raw[base + k_start[z] + local] =
                        boundary_depth(zone_specs[z].conformity, zt, zb, local, znk);
                }
            }
            let col = &mut snap[base..base + ncorner];
            col.copy_from_slice(&raw[base..base + ncorner]);
            if let Some(threshold) = collapse_below_m {
                for (z, &znk) in per_zone.iter().enumerate() {
                    collapse_zone(col, k_start[z], k_start[z] + znk, threshold);
                }
            }
        }
    }

    // COORD: pillars span the whole stack (surfaces[0] top -> surfaces[last] base).
    // Recycle the caller's buffer (clear + reserve reuses capacity; every pillar is
    // pushed, so no stale entry survives).
    let base = surfaces[surfaces.len() - 1];
    coord_buf.clear();
    coord_buf.reserve(dims.pillar_count());
    for jp in 0..ny {
        for ip in 0..nx {
            let x = ip as f64 * dx;
            let y = jp as f64 * dy;
            coord_buf.push(Pillar {
                top: Point3::new(x, y, surfaces[0].z(ip, jp)),
                bottom: Point3::new(x, y, base.z(ip, jp)),
            });
        }
    }

    // ZCORN: the 8 corners of cell `c` read the four pillars' interface columns at
    // levels `c.k` (top) and `c.k + 1` (bottom). A cell truncates (was already zero
    // pre-collapse) or collapses (zeroed by the collapse pass) only where ALL four
    // pillars agree — the max-pillar rule.
    let corner_z = |arr: &[f64], ip: usize, jp: usize, level: usize| -> f64 {
        arr[(jp * nx + ip) * ncorner + level]
    };
    // Recycle the caller's ZCORN buffer (the 64 MB/1M-cell dominant): clear +
    // reserve reuses its capacity; every cell's 8 corners are pushed, a full
    // overwrite.
    zcorn_buf.clear();
    zcorn_buf.reserve(dims.cell_count() * 8);
    let mut trunc_by_zone = vec![0usize; per_zone.len()];
    let mut collapse_by_zone = vec![0usize; per_zone.len()];
    for c in dims.iter() {
        let z = zone_of_k[c.k];
        let mut corners = [0.0_f64; 8];
        for (corner, slot) in corners.iter_mut().enumerate() {
            let di = corner & 1;
            let dj = (corner >> 1) & 1;
            let dk = (corner >> 2) & 1;
            *slot = corner_z(snap.as_slice(), c.i + di, c.j + dj, c.k + dk);
        }
        let max_dz = (0..4)
            .map(|p| corners[4 + p] - corners[p])
            .fold(0.0_f64, f64::max);
        let max_dz_raw = (0..4)
            .map(|p| {
                let di = p & 1;
                let dj = (p >> 1) & 1;
                corner_z(raw.as_slice(), c.i + di, c.j + dj, c.k + 1)
                    - corner_z(raw.as_slice(), c.i + di, c.j + dj, c.k)
            })
            .fold(0.0_f64, f64::max);
        if max_dz_raw <= TRUNC_EPS {
            trunc_by_zone[z] += 1; // already zero before collapse: a truncation
        } else if max_dz <= TRUNC_EPS {
            collapse_by_zone[z] += 1; // zeroed by the collapse pass
        }
        zcorn_buf.extend_from_slice(&corners);
    }

    // Write the reused buffers back into the caller's scratch + out-params (O(1)).
    scratch.raw = raw;
    scratch.snap = snap;
    *coord = coord_buf;
    *zcorn = zcorn_buf;

    let truncated_cells: usize = trunc_by_zone.iter().sum();
    let zones = per_zone
        .iter()
        .enumerate()
        .map(|(z, &nk)| StackedZone {
            nk,
            k_start: k_start[z],
            conformity: zone_specs[z].conformity,
            truncated_cells: trunc_by_zone[z],
            collapsed_cells: collapse_by_zone[z],
        })
        .collect();

    Ok(StackLayering {
        dims,
        nk: total_nk,
        zones,
        truncated_cells,
        collapsed_cells: collapse_by_zone.iter().sum(),
        nk_capped,
    })
}

/// A **streaming** producer of a stacked layering's ZCORN — one k-slab at a time,
/// **never materializing the whole grid** (the out-of-core slab-incremental build,
/// R2 follow-up). Where [`layer_grid_stack`] allocates the full `cell_count·8`
/// ZCORN Vec (the 64 MB/1M-cell dominant), this computes each k-slab's corner
/// depths on demand from two areal interface planes (levels `k` and `k+1`), so the
/// build's peak interface state is `O(area)`, not `O(grid)`.
///
/// The corner depths are computed by the **same** [`boundary_depth`] as the
/// in-core path, so a slab produced here is **bit-identical** to the corresponding
/// slab of a fully-built grid (before any f32 narrowing the caller applies).
///
/// **Scope:** collapse (`collapse_below_m`) is **not** supported here — the
/// volume-conserving cell-collapse pass is inherently a whole-zone-column operation
/// (it needs every interface level of a pillar's zone at once), so it cannot run
/// from a two-plane window. A spilled build that requests collapse falls back to
/// the build-then-spill path. Truncation (pinch-out to zero volume) **is** honoured
/// and counted, since it is a per-cell test on the two bounding planes.
#[derive(Debug, Clone)]
pub struct StreamingLayering {
    dims: Dims,
    nk: usize,
    nk_capped: bool,
    surfaces: Vec<Surface>,
    conformities: Vec<Conformity>,
    per_zone: Vec<usize>,
    k_start: Vec<usize>,
    zone_of_k: Vec<usize>,
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
}

impl StreamingLayering {
    /// Prepare a streaming layering over `surfaces` (N top→down surfaces sharing one
    /// `(ni+1)×(nj+1)` lattice) with one [`ZoneLayerSpec`] per interval. Validates
    /// the inputs and resolves the per-zone layer counts + total cap up front (all
    /// `O(area)`); no ZCORN is built until [`StreamingLayering::fill_zcorn_slab`].
    ///
    /// # Errors
    /// As [`layer_grid_stack`] (minus collapse, which this path does not accept).
    pub fn prepare(
        surfaces: &[&Surface],
        dx: f64,
        dy: f64,
        zone_specs: &[ZoneLayerSpec],
    ) -> Result<Self, StaticError> {
        if surfaces.len() < 2 {
            return Err(StaticError::InvalidInput(format!(
                "a horizon stack needs at least 2 surfaces, got {}",
                surfaces.len()
            )));
        }
        if zone_specs.len() != surfaces.len() - 1 {
            return Err(StaticError::InvalidInput(format!(
                "expected {} zone specs for {} surfaces, got {}",
                surfaces.len() - 1,
                surfaces.len(),
                zone_specs.len()
            )));
        }
        let nx = surfaces[0].nx();
        let ny = surfaces[0].ny();
        for s in surfaces {
            if s.nx() != nx || s.ny() != ny {
                return Err(StaticError::InvalidInput(
                    "all stack surfaces must share one lattice".into(),
                ));
            }
        }
        require_finite_surfaces(surfaces)?;
        for spec in zone_specs {
            if let Some(dz) = spec.conformity.dz_m() {
                if !(dz.is_finite() && dz > 0.0) {
                    return Err(StaticError::InvalidInput(format!(
                        "conformity dz_m must be finite and > 0, got {dz}"
                    )));
                }
            } else if spec.requested_nk == 0 {
                return Err(StaticError::InvalidInput("nk must be >= 1".into()));
            }
        }
        if !(dx.is_finite() && dx > 0.0 && dy.is_finite() && dy > 0.0) {
            return Err(StaticError::InvalidInput(format!(
                "dx, dy must be finite and > 0, got {dx}, {dy}"
            )));
        }
        if zone_specs.len() > MAX_NK {
            return Err(StaticError::InvalidInput(format!(
                "{} zones exceed the {MAX_NK}-layer total cap (a zone needs >= 1 layer)",
                zone_specs.len()
            )));
        }

        let mut zone_capped = false;
        let per_zone: Vec<usize> = zone_specs
            .iter()
            .zip(surfaces.windows(2))
            .map(|(spec, pair)| {
                let (n, capped) = derive_nk(pair[0], pair[1], spec.conformity, spec.requested_nk);
                zone_capped |= capped;
                n
            })
            .collect();
        let (per_zone, total_capped) = cap_total_nk(&per_zone);
        let nk_capped = zone_capped || total_capped;

        let total_nk: usize = per_zone.iter().sum();
        let mut k_start = Vec::with_capacity(per_zone.len());
        let mut zone_of_k = Vec::with_capacity(total_nk);
        let mut acc = 0usize;
        for (z, &n) in per_zone.iter().enumerate() {
            k_start.push(acc);
            for _ in 0..n {
                zone_of_k.push(z);
            }
            acc += n;
        }
        let dims = Dims::new(nx - 1, ny - 1, total_nk)?;

        Ok(Self {
            dims,
            nk: total_nk,
            nk_capped,
            surfaces: surfaces.iter().map(|s| (*s).clone()).collect(),
            conformities: zone_specs.iter().map(|s| s.conformity).collect(),
            per_zone,
            k_start,
            zone_of_k,
            nx,
            ny,
            dx,
            dy,
        })
    }

    /// The grid dimensions of the streamed build (`nk` = total layer count).
    #[must_use]
    pub fn dims(&self) -> Dims {
        self.dims
    }

    /// Whether the summed layer count hit [`MAX_NK`] and was scaled down.
    #[must_use]
    pub fn nk_capped(&self) -> bool {
        self.nk_capped
    }

    /// The per-zone layering report (k-ranges, conformity). `truncated_cells` /
    /// `collapsed_cells` are `0` here — truncation is counted per slab by
    /// [`StreamingLayering::fill_zcorn_slab`] (collapse is unsupported, always `0`).
    #[must_use]
    pub fn zones(&self) -> Vec<StackedZone> {
        self.per_zone
            .iter()
            .enumerate()
            .map(|(z, &nk)| StackedZone {
                nk,
                k_start: self.k_start[z],
                conformity: self.conformities[z],
                truncated_cells: 0,
                collapsed_cells: 0,
            })
            .collect()
    }

    /// Fill the pillar-lattice COORD into `out` (recycled): each pillar spans the
    /// top surface (`surfaces[0]`) to the base (`surfaces[last]`), areal spacing
    /// `dx, dy`. `O(area)`, k-invariant — computed once for the whole build.
    pub fn fill_coord(&self, out: &mut Vec<Pillar>) {
        out.clear();
        out.reserve(self.dims.pillar_count());
        let base = &self.surfaces[self.surfaces.len() - 1];
        for jp in 0..self.ny {
            for ip in 0..self.nx {
                out.push(Pillar {
                    top: Point3::new(
                        ip as f64 * self.dx,
                        jp as f64 * self.dy,
                        self.surfaces[0].z(ip, jp),
                    ),
                    bottom: Point3::new(ip as f64 * self.dx, jp as f64 * self.dy, base.z(ip, jp)),
                });
            }
        }
    }

    /// Interface-boundary depth of layer-boundary `level` (0..=nk) at pillar
    /// `(ip, jp)` — the [`boundary_depth`] of the zone that owns that level.
    #[inline]
    fn level_depth(&self, ip: usize, jp: usize, level: usize) -> f64 {
        let (z, local) = if level >= self.nk {
            let z = self.per_zone.len() - 1;
            (z, self.per_zone[z])
        } else {
            let z = self.zone_of_k[level];
            (z, level - self.k_start[z])
        };
        boundary_depth(
            self.conformities[z],
            self.surfaces[z].z(ip, jp),
            self.surfaces[z + 1].z(ip, jp),
            local,
            self.per_zone[z],
        )
    }

    /// Fill the areal interface plane for boundary `level` into `plane` (length
    /// `nx·ny`). `O(area)`.
    fn fill_plane(&self, level: usize, plane: &mut [f64]) {
        for jp in 0..self.ny {
            for ip in 0..self.nx {
                plane[jp * self.nx + ip] = self.level_depth(ip, jp, level);
            }
        }
    }

    /// Fill k-slab `k`'s ZCORN corner depths (as **f32**) into `out` (length
    /// `ni·nj·8`, cell-major `j·ni+i`, 8 local corners), computing the two bounding
    /// interface planes into the caller's `plane_top`/`plane_bot` scratch (each
    /// `nx·ny`). Returns this slab's **truncated** (zero-thickness) cell count.
    ///
    /// The f64 corner depths are bit-identical to a fully-built grid's; narrowing to
    /// f32 here matches the spill store's f32 ZCORN lane (R4), so a streamed spilled
    /// build is bit-identical to build-then-spill.
    pub fn fill_zcorn_slab(
        &self,
        k: usize,
        plane_top: &mut [f64],
        plane_bot: &mut [f64],
        out: &mut [f32],
    ) -> usize {
        self.fill_plane(k, plane_top);
        self.fill_plane(k + 1, plane_bot);
        let (ni, nj) = (self.dims.ni, self.dims.nj);
        let mut truncated = 0usize;
        for j in 0..nj {
            for i in 0..ni {
                let local = j * ni + i;
                let base = local * 8;
                let mut max_dz = 0.0f64;
                for corner in 0..8 {
                    let di = corner & 1;
                    let dj = (corner >> 1) & 1;
                    let dk = (corner >> 2) & 1;
                    let plane: &[f64] = if dk == 0 { plane_top } else { plane_bot };
                    out[base + corner] = plane[(j + dj) * self.nx + (i + di)] as f32;
                }
                for p in 0..4 {
                    let di = p & 1;
                    let dj = (p >> 1) & 1;
                    let idx = (j + dj) * self.nx + (i + di);
                    max_dz = max_dz.max(plane_bot[idx] - plane_top[idx]);
                }
                if max_dz <= TRUNC_EPS {
                    truncated += 1;
                }
            }
        }
        truncated
    }
}

/// Reject a non-finite surface node up front. A NaN/inf depth silently propagates
/// through [`boundary_depth`] into NaN ZCORN corners — a whole-grid poison no
/// downstream volumetrics can distinguish from a real value — so the layering seam
/// refuses it **loudly** (R4 loudness / R5 degenerate-input), naming the offending
/// surface index and node. The builder guards *all-NaN* surfaces upstream, but a
/// **partial**-NaN surface (only some nodes undefined) slips that check and would
/// poison exactly those cells; this catches it at the kernel seam that every layering
/// path funnels through.
fn require_finite_surfaces(surfaces: &[&Surface]) -> Result<(), StaticError> {
    for (s_idx, s) in surfaces.iter().enumerate() {
        for jp in 0..s.ny() {
            for ip in 0..s.nx() {
                let z = s.z(ip, jp);
                if !z.is_finite() {
                    return Err(StaticError::InvalidInput(format!(
                        "surface {s_idx} has a non-finite depth {z} at node ({ip}, {jp}); \
                         layering needs finite top/base surfaces"
                    )));
                }
            }
        }
    }
    Ok(())
}

/// Size a reused scratch buffer to exactly `n` elements **without** re-zeroing the
/// reused region — the caller overwrites every element before reading it. Grows with
/// a one-off `resize` (only the new tail is written) and shrinks/keeps with a
/// write-free `truncate`, so the steady-state (constant-`n`) path does no memset.
fn size_scratch(v: &mut Vec<f64>, n: usize) {
    if v.len() < n {
        v.resize(n, 0.0);
    } else {
        v.truncate(n);
    }
}

/// Volume-conserving cell-collapse on one pillar's interface column, restricted to
/// one zone's interface band `[lo, hi]` (`hi - lo` layers). Any zone-interior layer
/// thinner than `threshold` (but still positive — an already-zero truncated layer
/// is skipped) has its sliver thickness merged into the nearest **positive** vertical
/// neighbour (looking through any run of already-zeroed layers) by snapping the
/// intervening interfaces, so the sliver goes to zero thickness while total zone
/// thickness is preserved. The zone's bounding interfaces (`lo`, `hi`) are never
/// moved, so collapse never crosses a zone boundary. A single-layer zone
/// (`hi - lo == 1`) has no interior neighbour, so it is left untouched.
fn collapse_zone(col: &mut [f64], lo: usize, hi: usize, threshold: f64) {
    if hi <= lo + 1 {
        return; // single-layer zone: nothing to merge into
    }
    // Degenerate zone: the whole column is thinner than the threshold, so NO
    // arrangement can lift any single layer to/above the threshold. Whichever
    // interior layer holds the thickness is itself a sub-threshold sliver with an
    // interior neighbour, so the merge loop below would ping-pong the thickness
    // between interior layers forever (question_collapse_zone_livelock, repro
    // [0.0, 0.3, 0.4] @ 0.5). Snap the whole zone onto its thickest single layer
    // in one step — volume-conserving (total zone thickness preserved on that
    // layer, every other interior layer goes to zero) — and return.
    let (top, base) = (col[lo], col[hi]);
    let total = base - top;
    if total > TRUNC_EPS && total < threshold {
        let mut survivor = lo;
        let mut thickest = col[lo + 1] - col[lo];
        for k in (lo + 1)..hi {
            let t = col[k + 1] - col[k];
            if t > thickest {
                thickest = t;
                survivor = k;
            }
        }
        // Interfaces above the survivor snap to the zone top, those below to the
        // zone base; the survivor layer then spans the full [lo, hi] thickness.
        for c in col.iter_mut().take(survivor + 1).skip(lo + 1) {
            *c = top;
        }
        for c in col.iter_mut().take(hi).skip(survivor + 1) {
            *c = base;
        }
        return;
    }
    // Repeatedly collapse the thinnest sub-threshold, still-positive layer into the
    // nearest **positive** neighbour, until none remain.
    //
    // Progress guarantee (why this cannot livelock): the sliver is merged into a
    // neighbour that is itself positive, looking *through* any run of already-zeroed
    // layers to the nearest positive layer/rock on each side. The chosen neighbour
    // stays positive (it only grows) and the sliver goes to zero permanently, so
    // every iteration strictly reduces the count of positive interior layers → the
    // loop terminates in at most `hi - lo` steps.
    //
    // The earlier "merge into the immediate thicker neighbour" rule livelocked when a
    // sliver's immediate neighbours were both already zero-thickness: the sliver just
    // ping-ponged between adjacent zero slots without the positive-layer count ever
    // falling (repro [0.83, 0.0, 0.0, 1.397] @ 0.896, found by the R5 collapse
    // proptest; question_collapse_zone_livelock, second instance). Since total >=
    // threshold here (the whole-zone-degenerate case returned above), a sliver <
    // threshold always leaves positive rock elsewhere, so a positive neighbour exists.
    loop {
        let mut target: Option<usize> = None;
        let mut thinnest = f64::INFINITY;
        for k in lo..hi {
            let t = col[k + 1] - col[k];
            if t > TRUNC_EPS && t < threshold && t < thinnest {
                thinnest = t;
                target = Some(k);
            }
        }
        let Some(k) = target else { break };
        // Nearest positive layer to the left / right (scanning through zeroed layers).
        let left = (lo..k).rev().find(|&l| col[l + 1] - col[l] > TRUNC_EPS);
        let right = ((k + 1)..hi).find(|&r| col[r + 1] - col[r] > TRUNC_EPS);
        let merge_left = match (left, right) {
            (Some(l), Some(r)) => (col[l + 1] - col[l]) >= (col[r + 1] - col[r]),
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => break, // no positive neighbour (degenerate handled above)
        };
        if merge_left {
            let l = left.unwrap();
            // Push the sliver (and the zero run between) down into layer l: snap every
            // interface l+1..=k to the sliver's base, so l grows and l+1..k go to zero.
            let anchor = col[k + 1];
            for c in col.iter_mut().take(k + 1).skip(l + 1) {
                *c = anchor;
            }
        } else {
            let r = right.unwrap();
            // Push the sliver (and the zero run between) up into layer r: snap every
            // interface k+1..=r to the sliver's top, so r grows and k..r-1 go to zero.
            let anchor = col[k];
            for c in col.iter_mut().take(r + 1).skip(k + 1) {
                *c = anchor;
            }
        }
    }
}

/// Apply the [`MAX_NK`] cap to the **total** of per-zone layer counts. Returns the
/// (possibly reduced) per-zone counts and whether the cap bit. Scaling is
/// proportional with every zone kept `>= 1`; any rounding residue over the cap is
/// trimmed from the thickest zones. Callers guarantee `counts.len() <= MAX_NK`.
fn cap_total_nk(counts: &[usize]) -> (Vec<usize>, bool) {
    let total: usize = counts.iter().sum();
    if total <= MAX_NK {
        return (counts.to_vec(), false);
    }
    let scale = MAX_NK as f64 / total as f64;
    let mut scaled: Vec<usize> = counts
        .iter()
        .map(|&n| ((n as f64 * scale).floor() as usize).max(1))
        .collect();
    // Rounding can leave the sum a little over MAX_NK; trim the largest zones
    // (never below 1) until it fits.
    let mut over = scaled.iter().sum::<usize>().saturating_sub(MAX_NK);
    while over > 0 {
        let Some(idx) = scaled
            .iter()
            .enumerate()
            .filter(|(_, &n)| n > 1)
            .max_by_key(|(_, &n)| n)
            .map(|(i, _)| i)
        else {
            break; // every zone already at 1 — cannot trim further
        };
        scaled[idx] -= 1;
        over -= 1;
    }
    (scaled, true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gridder::surface::{solve_surface, Control, SolveOpts};

    /// A flat surface at constant depth over an `nx x ny` lattice.
    fn flat(nx: usize, ny: usize, z: f64) -> Surface {
        let c: Vec<Control> = (0..nx)
            .flat_map(|ip| (0..ny).map(move |jp| Control { ip, jp, z }))
            .collect();
        solve_surface(nx, ny, &c, SolveOpts::default()).unwrap()
    }

    #[test]
    fn flat_surfaces_proportional_recover_a_box() {
        // 10x10 cells, top 5000, base 5050, 5 layers, 100x80 ft cells.
        let top = flat(11, 11, 5000.0);
        let base = flat(11, 11, 5050.0);
        let lg = layer_grid(&top, &base, 100.0, 80.0, 5, Conformity::Proportional).unwrap();
        let grid = &lg.grid;
        assert_eq!(lg.nk, 5);
        assert_eq!(lg.truncated_cells, 0);
        assert!(!lg.nk_capped);
        assert_eq!(grid.cell_count(), 10 * 10 * 5);
        // Bulk volume == footprint area * gross thickness.
        let expected = (10.0 * 100.0) * (10.0 * 80.0) * 50.0;
        assert!((grid.bulk_volume() - expected).abs() / expected < 1e-9);
        // Each layer is 10 ft thick.
        assert!((grid.cell(crate::grid::Ijk::new(0, 0, 0)).dz() - 10.0).abs() < 1e-6);
    }

    #[test]
    fn proportional_layers_split_thickness_evenly() {
        let top = flat(6, 6, 4000.0);
        let base = flat(6, 6, 4100.0);
        let lg = layer_grid(&top, &base, 50.0, 50.0, 4, Conformity::Proportional).unwrap();
        for k in 0..4 {
            let dz = lg.grid.cell(crate::grid::Ijk::new(0, 0, k)).dz();
            assert!((dz - 25.0).abs() < 1e-6, "layer {k} dz={dz}");
        }
    }

    /// A wedge: flat top at 5000, base dipping only in i so node `ip` sits at
    /// `5000 + (ip+1)·10` (thickness 10 at ip=0 → 110 at ip=10), constant in j.
    fn wedge(nx: usize) -> (Surface, Surface) {
        let top = flat(nx, nx, 5000.0);
        let ctrl: Vec<Control> = (0..nx)
            .flat_map(|ip| {
                (0..nx).map(move |jp| Control {
                    ip,
                    jp,
                    z: 5000.0 + (ip as f64 + 1.0) * 10.0,
                })
            })
            .collect();
        let base = solve_surface(
            nx,
            nx,
            &ctrl,
            SolveOpts {
                tol: 1e-9,
                max_iter: 60_000,
                ..SolveOpts::default()
            },
        )
        .unwrap();
        (top, base)
    }

    #[test]
    fn follow_top_derives_nk_and_layers_parallel_to_the_top() {
        // dz=10, max thickness 110 -> nk = ceil(110/10) = 11.
        let (top, base) = wedge(11);
        let lg = layer_grid(
            &top,
            &base,
            10.0,
            10.0,
            999,
            Conformity::FollowTop { dz_m: 10.0 },
        )
        .unwrap();
        assert_eq!(lg.nk, 11, "nk derived from max thickness / dz");
        assert!(!lg.nk_capped);
        // Fully-active shallow layers are EXACTLY dz thick everywhere (parallel to
        // the flat top) — the truncation only bites the deep layers.
        for k in 0..1 {
            for i in [0usize, 5, 9] {
                let dz = lg.grid.cell(crate::grid::Ijk::new(i, 0, k)).dz();
                assert!((dz - 10.0).abs() < 1e-6, "cell ({i},0,{k}) dz={dz} != 10");
            }
        }
    }

    #[test]
    fn follow_top_truncation_count_is_exact() {
        // Closed form (max-pillar rule): cell col i (nodes i,i+1; deeper base at
        // i+1 = 5000+(i+2)·10) is active for k with k·dz < (i+2)·dz i.e. k<=i+1,
        // so active = min(i+2, nk); truncated per col = nk - active. nk=11.
        let (top, base) = wedge(11);
        let lg = layer_grid(
            &top,
            &base,
            10.0,
            10.0,
            0,
            Conformity::FollowTop { dz_m: 10.0 },
        )
        .unwrap();
        let nk = lg.nk as i64;
        let nj = 10i64;
        let expected: i64 = (0..10).map(|i| (nk - (i + 2).min(nk)).max(0)).sum::<i64>() * nj;
        assert_eq!(expected, 450, "analytic truncated count");
        assert_eq!(lg.truncated_cells as i64, expected);
    }

    #[test]
    fn volume_conserved_across_conformity_styles() {
        // Total bulk volume is conformity-invariant (all styles fill exactly the
        // top->base wedge). FollowTop / FollowBase == Proportional to FP.
        let (top, base) = wedge(11);
        let prop = layer_grid(&top, &base, 10.0, 10.0, 7, Conformity::Proportional)
            .unwrap()
            .grid
            .bulk_volume();
        let ftop = layer_grid(
            &top,
            &base,
            10.0,
            10.0,
            0,
            Conformity::FollowTop { dz_m: 10.0 },
        )
        .unwrap()
        .grid
        .bulk_volume();
        let fbase = layer_grid(
            &top,
            &base,
            10.0,
            10.0,
            0,
            Conformity::FollowBase { dz_m: 10.0 },
        )
        .unwrap()
        .grid
        .bulk_volume();
        assert!(
            (ftop - prop).abs() / prop < 1e-9,
            "FollowTop {ftop} != {prop}"
        );
        assert!(
            (fbase - prop).abs() / prop < 1e-9,
            "FollowBase {fbase} != {prop}"
        );
    }

    #[test]
    fn follow_base_mirrors_follow_top() {
        // FollowBase drapes parallel to the (dipping) base; the deepest layer is
        // dz thick everywhere, the top onlaps.
        let (top, base) = wedge(11);
        let lg = layer_grid(
            &top,
            &base,
            10.0,
            10.0,
            0,
            Conformity::FollowBase { dz_m: 10.0 },
        )
        .unwrap();
        assert_eq!(lg.nk, 11);
        // Deepest layer sits against the base at dz thickness in the thick columns.
        let dz = lg.grid.cell(crate::grid::Ijk::new(9, 0, lg.nk - 1)).dz();
        assert!((dz - 10.0).abs() < 1e-6, "deepest layer dz={dz} != 10");
        // Truncation happens near the TOP (mirror of FollowTop), same total count.
        assert_eq!(lg.truncated_cells, 450);
    }

    #[test]
    fn follow_top_caps_nk_and_flags_it() {
        // dz tiny -> derived nk would be huge; capped at MAX_NK with the flag set.
        let top = flat(4, 4, 5000.0);
        let base = flat(4, 4, 5000.0 + (MAX_NK as f64 + 50.0)); // needs > MAX_NK layers
        let lg = layer_grid(
            &top,
            &base,
            10.0,
            10.0,
            0,
            Conformity::FollowTop { dz_m: 1.0 },
        )
        .unwrap();
        assert_eq!(lg.nk, MAX_NK);
        assert!(lg.nk_capped, "cap flag set");
    }

    #[test]
    fn rejects_nonpositive_dz() {
        let top = flat(4, 4, 5000.0);
        let base = flat(4, 4, 5050.0);
        assert!(layer_grid(
            &top,
            &base,
            10.0,
            10.0,
            0,
            Conformity::FollowTop { dz_m: 0.0 }
        )
        .is_err());
        assert!(layer_grid(
            &top,
            &base,
            10.0,
            10.0,
            0,
            Conformity::FollowBase { dz_m: f64::NAN }
        )
        .is_err());
    }

    #[test]
    fn tilted_top_follows_structure() {
        // Top dips 2 ft per node in i; base flat. Cells should thicken updip.
        let plane = |ip: usize, jp: usize| 5000.0 + 2.0 * ip as f64 + 0.0 * jp as f64;
        let ctrl: Vec<Control> = [(0usize, 0usize), (10, 0), (0, 10), (10, 10), (5, 5)]
            .iter()
            .map(|&(ip, jp)| Control {
                ip,
                jp,
                z: plane(ip, jp),
            })
            .collect();
        let top = solve_surface(
            11,
            11,
            &ctrl,
            SolveOpts {
                tol: 1e-9,
                max_iter: 60_000,
                ..SolveOpts::default()
            },
        )
        .unwrap();
        let base = flat(11, 11, 5100.0);
        let grid = layer_grid(&top, &base, 100.0, 100.0, 1, Conformity::Proportional)
            .unwrap()
            .grid;
        let updip = grid.cell(crate::grid::Ijk::new(0, 0, 0)).dz();
        let downdip = grid.cell(crate::grid::Ijk::new(9, 0, 0)).dz();
        assert!(
            updip > downdip,
            "updip {updip} should be thicker than downdip {downdip}"
        );
    }

    #[test]
    fn follow_top_pinches_out_flat_column() {
        // Top 5000, base 5030 (flat, thickness 30), FollowTop dz=20 -> nk =
        // ceil(30/20) = 2; layer0 = 20, layer1 truncates to 10 (clipped at base).
        let top = flat(4, 4, 5000.0);
        let base = flat(4, 4, 5030.0);
        let lg = layer_grid(
            &top,
            &base,
            50.0,
            50.0,
            0,
            Conformity::FollowTop { dz_m: 20.0 },
        )
        .unwrap();
        assert_eq!(lg.nk, 2);
        let l0 = lg.grid.cell(crate::grid::Ijk::new(0, 0, 0)).dz();
        let l1 = lg.grid.cell(crate::grid::Ijk::new(0, 0, 1)).dz();
        assert!((l0 - 20.0).abs() < 1e-6, "layer0 {l0}");
        assert!((l1 - 10.0).abs() < 1e-6, "layer1 {l1} (clipped at base)");
        // A partial-thickness cell is NOT zero-volume, so nothing is truncated here.
        assert_eq!(lg.truncated_cells, 0);
    }

    #[test]
    fn rejects_mismatched_lattice() {
        let top = flat(5, 5, 1.0);
        let base = flat(6, 6, 2.0);
        assert!(layer_grid(&top, &base, 10.0, 10.0, 2, Conformity::Proportional).is_err());
    }

    // --- multi-zone stack layering ---

    fn prop(nk: usize) -> ZoneLayerSpec {
        ZoneLayerSpec {
            conformity: Conformity::Proportional,
            requested_nk: nk,
        }
    }

    #[test]
    fn stack_concatenates_zones_and_partitions_k() {
        // Three flat surfaces (5000, 5050, 5090) -> two zones of 4 + 3 layers.
        let s0 = flat(6, 6, 5000.0);
        let s1 = flat(6, 6, 5050.0);
        let s2 = flat(6, 6, 5090.0);
        let sg = layer_grid_stack(&[&s0, &s1, &s2], 50.0, 50.0, &[prop(4), prop(3)], None).unwrap();
        assert_eq!(sg.nk, 7, "total nk = 4 + 3");
        assert_eq!(sg.zones.len(), 2);
        assert_eq!(sg.zones[0].k_range(), 0..4);
        assert_eq!(sg.zones[1].k_range(), 4..7);
        assert_eq!(sg.truncated_cells, 0);
        // Zone 0 layers are 50/4 = 12.5 m; zone 1 layers are 40/3 m.
        assert!((sg.grid.cell(crate::grid::Ijk::new(0, 0, 0)).dz() - 12.5).abs() < 1e-6);
        assert!((sg.grid.cell(crate::grid::Ijk::new(0, 0, 4)).dz() - 40.0 / 3.0).abs() < 1e-6);
        // Bulk volume == footprint * total gross (90 m over a 250 m square).
        let expected = 250.0 * 250.0 * 90.0;
        assert!((sg.grid.bulk_volume() - expected).abs() / expected < 1e-9);
    }

    #[test]
    fn stack_two_surface_matches_single_layer_grid() {
        // The 2-surface stack must equal the old single-zone layer_grid bit-for-bit.
        let (top, base) = wedge(11);
        let a = layer_grid(&top, &base, 10.0, 10.0, 7, Conformity::Proportional).unwrap();
        let b = layer_grid_stack(&[&top, &base], 10.0, 10.0, &[prop(7)], None).unwrap();
        assert_eq!(a.nk, b.nk);
        assert_eq!(a.truncated_cells, b.truncated_cells);
        assert_eq!(a.grid.bulk_volume(), b.grid.bulk_volume());
    }

    #[test]
    fn stack_mixes_conformity_per_zone() {
        // Zone 0 proportional, zone 1 FollowTop with its own dz -> per-zone nk.
        let s0 = flat(6, 6, 5000.0);
        let s1 = flat(6, 6, 5030.0);
        let s2 = flat(6, 6, 5090.0); // zone1 thickness 60, dz 20 -> 3 layers
        let sg = layer_grid_stack(
            &[&s0, &s1, &s2],
            50.0,
            50.0,
            &[
                prop(2),
                ZoneLayerSpec {
                    conformity: Conformity::FollowTop { dz_m: 20.0 },
                    requested_nk: 0,
                },
            ],
            None,
        )
        .unwrap();
        assert_eq!(sg.zones[0].nk, 2);
        assert_eq!(sg.zones[1].nk, 3, "60 m / 20 m dz");
        assert_eq!(sg.nk, 5);
    }

    #[test]
    fn stack_total_cap_scales_zones_down() {
        // Two zones each wanting 150 proportional layers -> 300 > MAX_NK -> capped.
        let s0 = flat(4, 4, 0.0);
        let s1 = flat(4, 4, 100.0);
        let s2 = flat(4, 4, 200.0);
        let sg =
            layer_grid_stack(&[&s0, &s1, &s2], 10.0, 10.0, &[prop(150), prop(150)], None).unwrap();
        assert!(sg.nk_capped);
        assert!(sg.nk <= MAX_NK, "total capped to <= MAX_NK, got {}", sg.nk);
        assert!(sg.zones.iter().all(|z| z.nk >= 1));
    }

    #[test]
    fn stack_rejects_bad_arity() {
        let s0 = flat(4, 4, 0.0);
        let s1 = flat(4, 4, 10.0);
        assert!(layer_grid_stack(&[&s0], 10.0, 10.0, &[], None).is_err()); // < 2 surfaces
        assert!(layer_grid_stack(&[&s0, &s1], 10.0, 10.0, &[prop(2), prop(2)], None).is_err());
        // wrong spec count
    }

    #[test]
    fn collapse_conserves_volume_and_stays_within_zones() {
        // Two flat zones: zone0 thickness 30 in 6 proportional layers -> 5 m each;
        // zone1 thickness 40 in 4 layers -> 10 m each. Collapse at 8 m zeroes zone0's
        // 5 m slivers (merged into neighbours) but never zone1's 10 m layers, and
        // total bulk volume is unchanged.
        let s0 = flat(4, 4, 0.0);
        let s1 = flat(4, 4, 30.0);
        let s2 = flat(4, 4, 70.0);
        let no = layer_grid_stack(&[&s0, &s1, &s2], 10.0, 10.0, &[prop(6), prop(4)], None).unwrap();
        let yes =
            layer_grid_stack(&[&s0, &s1, &s2], 10.0, 10.0, &[prop(6), prop(4)], Some(8.0)).unwrap();
        // (a) volume conserved to FP tolerance.
        let (vn, vy) = (no.grid.bulk_volume(), yes.grid.bulk_volume());
        assert!(
            (vn - vy).abs() / vn < 1e-12,
            "collapse conserved volume: {vn} vs {vy}"
        );
        // (b) zone0's 5 m layers all collapse (6 layers, interior-mergeable), zone1's
        //     10 m layers never do; per-zone counts exact.
        let ncell_zone = 3 * 3; // 3x3 cells per layer
        assert!(yes.zones[0].collapsed_cells > 0);
        assert_eq!(
            yes.zones[1].collapsed_cells, 0,
            "zone1 layers exceed threshold"
        );
        // (c) no collapse crosses a zone boundary: zone1 fully intact.
        for k in yes.zones[1].k_range() {
            for j in 0..3 {
                for i in 0..3 {
                    let dz = yes.grid.cell(crate::grid::Ijk::new(i, j, k)).dz();
                    assert!((dz - 10.0).abs() < 1e-9, "zone1 cell dz {dz} != 10");
                }
            }
        }
        // A single collapsed layer in zone0 zeroes exactly one layer's worth of cells.
        assert_eq!(yes.zones[0].collapsed_cells % ncell_zone, 0);
    }

    #[test]
    fn collapse_subthreshold_zone_terminates_and_conserves() {
        // Regression: a zone whose TOTAL thickness is positive but below the
        // threshold used to livelock — whichever interior layer held the zone's
        // thickness was itself a sub-threshold sliver, so the merge shuffled the
        // thickness between interior layers forever (question_collapse_zone_livelock).
        // The exact synthetic repro: interface column [0.0, 0.3, 0.4], threshold 0.5.
        // Runs on a worker thread with a hard timeout so a regression FAILS (hangs
        // manifest as a timeout) instead of wedging the whole test binary.
        use std::sync::mpsc;
        use std::time::Duration;
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let mut col = [0.0f64, 0.3, 0.4];
            collapse_zone(&mut col, 0, 2, 0.5);
            tx.send(col).ok();
        });
        let col = rx
            .recv_timeout(Duration::from_secs(5))
            .expect("collapse_zone must terminate on a sub-threshold zone column");
        // Volume conserved: total zone thickness (0.4) preserved.
        let total = col[2] - col[0];
        assert!(
            (total - 0.4).abs() < 1e-12,
            "total thickness {total} != 0.4"
        );
        // Degenerate zone snaps onto a single layer: exactly one layer holds the
        // whole 0.4, the other is zero thickness.
        let l0 = col[1] - col[0];
        let l1 = col[2] - col[1];
        assert!(l0 >= 0.0 && l1 >= 0.0, "no negative layer: {l0}, {l1}");
        assert!(
            ((l0 - 0.4).abs() < 1e-12 && l1 < TRUNC_EPS)
                || ((l1 - 0.4).abs() < 1e-12 && l0 < TRUNC_EPS),
            "one layer holds the full thickness: {l0}, {l1}"
        );
    }

    #[test]
    fn collapse_subthreshold_zone_multilayer_terminates() {
        // A multi-layer variant of the same degenerate column: total 0.4 spread
        // over many interior layers must still terminate and conserve in one snap.
        use std::sync::mpsc;
        use std::time::Duration;
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            // 30-layer zone, uniform 0.4/30 slivers, total 0.4 < threshold 0.5.
            let nk = 30usize;
            let mut col: Vec<f64> = (0..=nk).map(|k| 0.4 * k as f64 / nk as f64).collect();
            collapse_zone(&mut col, 0, nk, 0.5);
            tx.send(col).ok();
        });
        let col = rx
            .recv_timeout(Duration::from_secs(5))
            .expect("multi-layer sub-threshold collapse must terminate");
        let total = col[col.len() - 1] - col[0];
        assert!((total - 0.4).abs() < 1e-12, "total {total} != 0.4");
        // Exactly one layer non-zero.
        let nonzero = col.windows(2).filter(|w| w[1] - w[0] > TRUNC_EPS).count();
        assert_eq!(nonzero, 1, "degenerate zone snaps to a single layer");
    }

    #[test]
    fn collapse_single_layer_zone_is_untouched() {
        // A zone with one layer has no interior neighbour to merge into.
        let s0 = flat(4, 4, 0.0);
        let s1 = flat(4, 4, 2.0); // 2 m single layer, below a 5 m threshold
        let sg = layer_grid_stack(&[&s0, &s1], 10.0, 10.0, &[prop(1)], Some(5.0)).unwrap();
        assert_eq!(sg.collapsed_cells, 0, "single-layer zone cannot collapse");
        assert_eq!(sg.nk, 1);
    }

    // --- R4 loudness: a non-finite surface node is refused, never a silent
    // NaN-corner grid (`require_finite_surfaces`). ---

    #[test]
    fn nan_surface_is_a_typed_error_not_a_silent_nan_grid() {
        // A NaN base depth silently propagated into NaN ZCORN (poisoning every
        // downstream volume as NaN) — now a loud typed error naming the surface.
        let top = flat(6, 6, 5000.0);
        let base = top.offset_by(f64::NAN);
        let err = layer_grid(&top, &base, 50.0, 50.0, 4, Conformity::Proportional)
            .expect_err("a non-finite surface must be a typed error");
        assert!(
            matches!(err, StaticError::InvalidInput(ref m) if m.contains("non-finite")),
            "message must name the non-finite surface, got: {err}"
        );
    }

    #[test]
    fn nan_surface_in_stack_names_the_offending_surface() {
        // A 3-surface stack whose MIDDLE surface (index 1) is undefined: the error
        // names surface 1, not a generic failure.
        let s0 = flat(4, 4, 5000.0);
        let s1 = flat(4, 4, 5050.0).offset_by(f64::NAN);
        let s2 = flat(4, 4, 5090.0);
        let err = layer_grid_stack(&[&s0, &s1, &s2], 50.0, 50.0, &[prop(2), prop(2)], None)
            .expect_err("a non-finite stack surface must be a typed error");
        assert!(
            matches!(err, StaticError::InvalidInput(ref m) if m.contains("surface 1") && m.contains("non-finite")),
            "message must name surface 1, got: {err}"
        );
    }

    // --- R5 degenerate-input property tests: the collapse convergence loop over
    // random column profiles, each on a worker thread with a HARD TIMEOUT so a
    // livelock FAILS (times out) instead of hanging CI (the collapse_below_m
    // livelock class, question_collapse_zone_livelock). ---

    /// Run `f` on a worker thread; fail (not hang) if it does not finish in `secs`.
    fn bounded<T: Send + 'static>(secs: u64, f: impl FnOnce() -> T + Send + 'static) -> T {
        use std::sync::mpsc;
        use std::time::Duration;
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(f());
        });
        rx.recv_timeout(Duration::from_secs(secs))
            .expect("kernel must terminate within the hard timeout (livelock guard)")
    }

    proptest::proptest! {
        #![proptest_config(proptest::prelude::ProptestConfig::with_cases(200))]

        /// `collapse_zone` over an arbitrary non-decreasing interface column (any mix
        /// of zero / sub-threshold / above-threshold layers) TERMINATES and conserves
        /// volume: the zone endpoints are pinned (total thickness preserved) and the
        /// column stays non-decreasing (no inversion, no negative layer).
        #[test]
        fn prop_collapse_zone_terminates_and_conserves(
            gaps in proptest::collection::vec(0.0f64..1.5, 1..=10),
            threshold in 0.05f64..1.2,
        ) {
            let n = gaps.len();
            let mut col = Vec::with_capacity(n + 1);
            let mut acc = 4000.0f64;
            col.push(acc);
            for g in &gaps {
                acc += *g;
                col.push(acc);
            }
            let total_before = col[n] - col[0];
            let (top, base) = (col[0], col[n]);
            let out = bounded(5, move || {
                let mut c = col;
                collapse_zone(&mut c, 0, n, threshold);
                c
            });
            // (a) endpoints pinned → volume conserved.
            proptest::prop_assert!((out[0] - top).abs() < 1e-9);
            proptest::prop_assert!((out[n] - base).abs() < 1e-9);
            proptest::prop_assert!((( out[n] - out[0]) - total_before).abs() < 1e-9);
            // (b) non-decreasing: no interface crosses another (no negative layer).
            for k in 0..n {
                proptest::prop_assert!(
                    out[k + 1] - out[k] >= -1e-9,
                    "inverted layer {k}: {} -> {}", out[k], out[k + 1]
                );
            }
        }

        /// `layer_grid` over a random flat top + arbitrary uniform offset (inverted,
        /// zero, sub-threshold, or thick) and a random layer count / conformity: the
        /// call either succeeds with a FINITE, non-negative bulk volume or returns a
        /// typed error — never a panic, never a NaN grid, always within the timeout.
        #[test]
        fn prop_layer_grid_degenerate_offsets_are_finite_or_typed(
            offset in -80.0f64..250.0,
            nk in 1usize..8,
            follow in proptest::bool::ANY,
            dz in 0.5f64..40.0,
        ) {
            let conf = if follow { Conformity::FollowTop { dz_m: dz } } else { Conformity::Proportional };
            let res = bounded(5, move || {
                let top = flat(5, 5, 5000.0);
                let base = top.offset_by(offset);
                layer_grid(&top, &base, 50.0, 50.0, nk, conf).map(|lg| lg.grid.bulk_volume())
            });
            // Ok(finite volume) or a typed error are both acceptable — never a
            // panic, a NaN grid, or a hang.
            if let Ok(v) = res {
                proptest::prop_assert!(v.is_finite() && v >= -1e-6, "bulk volume {v}");
            }
        }
    }

    #[test]
    fn single_cell_single_layer_grid_is_built() {
        // The smallest possible grid: one areal cell (2x2 nodes), one layer.
        let top = flat(2, 2, 5000.0);
        let base = flat(2, 2, 5020.0);
        let lg = layer_grid(&top, &base, 25.0, 25.0, 1, Conformity::Proportional).unwrap();
        assert_eq!(lg.grid.cell_count(), 1);
        assert_eq!(lg.nk, 1);
        assert!(lg.grid.bulk_volume() > 0.0);
    }
}
