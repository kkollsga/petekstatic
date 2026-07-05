//! The **one home** for the bandwidth-bound GRV/HCPV inner loop.
//!
//! GRV + in-place volumetrics are structurally identical whether the per-cell
//! geometry + cubes live in a contiguous in-core [`Grid`] (f64) or in a
//! memory-mapped, k-slab-major spill store (f32, `srs-model`'s out-of-core
//! backing). Before this module the two were **hand-maintained twins** (grv.rs
//! and srs-model's spill.rs), a silent-parity-divergence risk. Here the loop has
//! a single home and one place to tune: a [`SlabSource`] abstracts the two
//! backings, and three monomorphic streaming cores ([`stream_single`],
//! [`stream_two`], [`stream_bulk`]) run over any source **one k-slab at a time**.
//!
//! Monomorphization keeps the codegen monolithic: each concrete backing
//! (`GridSource` here, `SpillSource` in srs-model) instantiates its own tight,
//! branch-free-on-clip loop — no vtable, no per-cell dispatch. The in-core
//! iteration order (`for k { for local }`, i.e. ascending linear cell index) is
//! byte-identical to the previous `dims.iter()` loop, so in-core results stay
//! **bit-exact**. f32 backings widen to f64 at the cell boundary; **accumulations
//! are always f64** (out-of-core ruling R4 honesty clause).

use crate::error::StaticError;
use crate::grid::{Dims, Grid};
use crate::volumetrics::grv::{InPlace, ZoneVolumes};
use crate::volumetrics::names::{NTG, PORO, SW};
use crate::volumetrics::valid::validate_fraction;
use core::ops::Range;

/// A source of per-cell volumetric inputs, addressed **one k-slab at a time**.
/// Implemented by [`GridSource`] (in-core f64) and by srs-model's spilled mmap
/// backing (f32). The generic streaming cores below iterate any `SlabSource`.
pub trait SlabSource {
    /// Grid dimensions of the backing.
    fn dims(&self) -> Dims;

    /// Validate that the volumetrics cubes (`PORO`/`NTG`/`SW`) are present — the
    /// hydrocarbon clips call this once up front; the bulk clip does not.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if a required cube is missing.
    fn require_cubes(&self) -> Result<(), StaticError>;

    /// The per-cell view of one k-slab (zero-copy where possible).
    type Slab<'a>: CellSlab
    where
        Self: 'a;

    /// Borrow k-slab `k` as a [`CellSlab`].
    ///
    /// # Errors
    /// [`StaticError::Algo`] on a store read failure (spilled backings only).
    fn slab(&self, k: usize) -> Result<Self::Slab<'_>, StaticError>;
}

/// The per-cell accessors over one k-slab. `local` is the slab-local cell index
/// (`0..ni·nj`); property values are widened to f64 by f32 backings.
pub trait CellSlab {
    /// Centroid depth of slab-local cell `local` (the cheap contact test).
    fn centroid_z(&self, local: usize) -> f64;
    /// Gross (hexahedron) volume of slab-local cell `local` \[m³\].
    fn cell_volume(&self, local: usize) -> f64;
    /// `PORO` of slab-local cell `local`.
    fn poro(&self, local: usize) -> f64;
    /// `NTG` of slab-local cell `local`.
    fn ntg(&self, local: usize) -> f64;
    /// `SW` of slab-local cell `local`.
    fn sw(&self, local: usize) -> f64;
}

/// The fluid-contact clip that classifies a cell into the accumulation buckets —
/// the shared resolution the in-core and spilled paths previously each re-derived.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Clip {
    /// Contactless: gross rock volume only, **zero** hydrocarbon (a zone with no
    /// known contact). No cube reads.
    Bulk,
    /// A single hard contact at `depth_m`: a cell is in-column iff its centroid is
    /// shallower (smaller z).
    Single(f64),
    /// Two contacts (gas cap above `goc`, oil leg `goc..owc`, water below `owc`);
    /// `sw_gas` optionally overrides the gas-cap connate water (R3).
    Two {
        /// Gas–oil contact depth \[m\].
        goc: f64,
        /// Oil–water / free-water contact depth \[m\].
        owc: f64,
        /// Optional gas-cap connate-water override.
        sw_gas: Option<f64>,
    },
}

impl Clip {
    /// Whether this clip reads the `PORO`/`NTG`/`SW` cubes (all but [`Clip::Bulk`]).
    #[must_use]
    pub fn needs_cubes(&self) -> bool {
        !matches!(self, Clip::Bulk)
    }
}

/// Per-cell hydrocarbon pore volume `V·NTG·φ·(1 - Sw)`; `sw_override` (gas-cap
/// connate water, R3) replaces the cube `Sw`. Inputs are already f64 (f32 backings
/// widen at the [`CellSlab`] boundary); validation is per-cell (H2).
#[inline]
fn cell_hcpv(
    volume: f64,
    poro: f64,
    ntg: f64,
    sw: f64,
    sw_override: Option<f64>,
) -> Result<f64, StaticError> {
    let phi = frac(PORO, poro)?;
    let n = frac(NTG, ntg)?;
    let water = match sw_override {
        Some(s) => s,
        None => frac(SW, sw)?,
    };
    Ok(volume * n * phi * (1.0 - water))
}

#[inline]
fn frac(what: &str, x: f64) -> Result<f64, StaticError> {
    validate_fraction(what, x)?;
    Ok(x)
}

/// Dispatch a [`Clip`] to its monomorphic streaming core (validating the clip and
/// the cube presence first). The single entry the grv.rs wrappers and the model's
/// backing dispatch both call — the loop's one home.
///
/// `k_range` restricts the slab sweep to a zone's layer band (`0..nk` = whole
/// grid); `per_cell` gates the `cell_count`-length HCPV cube (the MC summary path
/// leaves it empty).
///
/// # Errors
/// [`StaticError::OutOfRange`] on a non-finite contact; [`StaticError::InvalidInput`]
/// on `goc > owc`, a non-fraction `sw_gas`, or a missing/non-physical cube;
/// [`StaticError::Algo`] on a store read failure.
pub fn compute_clipped<S: SlabSource>(
    src: &S,
    clip: Clip,
    k_range: Range<usize>,
    per_cell: bool,
) -> Result<InPlace, StaticError> {
    match clip {
        Clip::Bulk => stream_bulk(src, k_range),
        Clip::Single(depth) => {
            if !depth.is_finite() {
                return Err(StaticError::OutOfRange(format!(
                    "contact depth must be finite, got {depth}"
                )));
            }
            src.require_cubes()?;
            stream_single(src, depth, k_range, per_cell)
        }
        Clip::Two { goc, owc, sw_gas } => {
            if !goc.is_finite() || !owc.is_finite() {
                return Err(StaticError::OutOfRange(format!(
                    "contact depths must be finite, got GOC {goc}, OWC {owc}"
                )));
            }
            if goc > owc {
                return Err(StaticError::InvalidInput(format!(
                    "GOC ({goc}) must be shallower than OWC ({owc})"
                )));
            }
            if let Some(s) = sw_gas {
                validate_fraction("sw_gas", s)?;
            }
            src.require_cubes()?;
            stream_two(src, goc, owc, sw_gas, k_range, per_cell)
        }
    }
}

/// Streaming single-contact GRV/HCPV over `k_range`. Monomorphic per backing.
fn stream_single<S: SlabSource>(
    src: &S,
    contact_depth_m: f64,
    k_range: Range<usize>,
    per_cell: bool,
) -> Result<InPlace, StaticError> {
    let dims = src.dims();
    let per_slab = dims.ni * dims.nj;
    let mut per_cell_hcpv = if per_cell {
        vec![0.0; dims.cell_count()]
    } else {
        Vec::new()
    };
    let (mut grv, mut hcpv, mut cells) = (0.0f64, 0.0f64, 0usize);
    for k in k_range {
        let slab = src.slab(k)?;
        let kbase = k * per_slab;
        for local in 0..per_slab {
            if slab.centroid_z(local) >= contact_depth_m {
                continue; // below contact — rejected before any corner build
            }
            let v = slab.cell_volume(local);
            if v <= 0.0 {
                continue; // truncated/pinched cell
            }
            let h = cell_hcpv(v, slab.poro(local), slab.ntg(local), slab.sw(local), None)?;
            if per_cell {
                per_cell_hcpv[kbase + local] = h;
            }
            grv += v;
            hcpv += h;
            cells += 1;
        }
    }
    Ok(InPlace {
        grv_m3: grv,
        hcpv_m3: hcpv,
        cells_in_column: cells,
        per_cell_hcpv,
        gas: None,
        oil: None,
    })
}

/// Streaming two-contact (gas cap + oil leg) GRV/HCPV over `k_range`.
fn stream_two<S: SlabSource>(
    src: &S,
    goc_m: f64,
    owc_m: f64,
    sw_gas: Option<f64>,
    k_range: Range<usize>,
    per_cell: bool,
) -> Result<InPlace, StaticError> {
    let dims = src.dims();
    let per_slab = dims.ni * dims.nj;
    let mut per_cell_hcpv = if per_cell {
        vec![0.0; dims.cell_count()]
    } else {
        Vec::new()
    };
    let (mut gg, mut gh, mut gc) = (0.0f64, 0.0f64, 0usize);
    let (mut og, mut oh, mut oc) = (0.0f64, 0.0f64, 0usize);
    for k in k_range {
        let slab = src.slab(k)?;
        let kbase = k * per_slab;
        for local in 0..per_slab {
            let z = slab.centroid_z(local);
            if z >= owc_m {
                continue; // water leg
            }
            let v = slab.cell_volume(local);
            if v <= 0.0 {
                continue;
            }
            let is_gas = z < goc_m;
            let h = cell_hcpv(
                v,
                slab.poro(local),
                slab.ntg(local),
                slab.sw(local),
                if is_gas { sw_gas } else { None },
            )?;
            if per_cell {
                per_cell_hcpv[kbase + local] = h;
            }
            if is_gas {
                gg += v;
                gh += h;
                gc += 1;
            } else {
                og += v;
                oh += h;
                oc += 1;
            }
        }
    }
    Ok(InPlace {
        grv_m3: gg + og,
        hcpv_m3: gh + oh,
        cells_in_column: gc + oc,
        per_cell_hcpv,
        gas: Some(ZoneVolumes {
            grv_m3: gg,
            hcpv_m3: gh,
            cells: gc,
        }),
        oil: Some(ZoneVolumes {
            grv_m3: og,
            hcpv_m3: oh,
            cells: oc,
        }),
    })
}

/// Streaming gross bulk volume over `k_range` (the contactless zone): summed
/// active-cell volume, zero hydrocarbon. Reads no cubes.
fn stream_bulk<S: SlabSource>(src: &S, k_range: Range<usize>) -> Result<InPlace, StaticError> {
    let dims = src.dims();
    let per_slab = dims.ni * dims.nj;
    let (mut grv, mut cells) = (0.0f64, 0usize);
    for k in k_range {
        let slab = src.slab(k)?;
        for local in 0..per_slab {
            let v = slab.cell_volume(local);
            if v <= 0.0 {
                continue;
            }
            grv += v;
            cells += 1;
        }
    }
    Ok(InPlace {
        grv_m3: grv,
        hcpv_m3: 0.0,
        cells_in_column: cells,
        per_cell_hcpv: Vec::new(),
        gas: None,
        oil: None,
    })
}

// --- the in-core backing: a contiguous f64 Grid ---------------------------

/// The in-core [`SlabSource`]: a contiguous f64 [`Grid`] + its volumetrics cube
/// slices. Constructing it is infallible (cube presence is checked lazily by
/// [`SlabSource::require_cubes`], so the bulk path needs no cubes).
pub struct GridSource<'a> {
    grid: &'a Grid,
    poro: Option<&'a [f64]>,
    ntg: Option<&'a [f64]>,
    sw: Option<&'a [f64]>,
}

impl<'a> GridSource<'a> {
    /// Wrap a grid, resolving the volumetrics cube slices (absent → `None`).
    #[must_use]
    pub fn new(grid: &'a Grid) -> Self {
        let get = |n: &str| grid.properties().get(n).map(|p| p.values.as_slice());
        Self {
            grid,
            poro: get(PORO),
            ntg: get(NTG),
            sw: get(SW),
        }
    }
}

/// One k-slab view of a [`GridSource`] — the grid + cube slices offset to the
/// slab's base cell index (`k·ni·nj`).
pub struct GridSlab<'a> {
    grid: &'a Grid,
    poro: Option<&'a [f64]>,
    ntg: Option<&'a [f64]>,
    sw: Option<&'a [f64]>,
    kbase: usize,
}

impl SlabSource for GridSource<'_> {
    fn dims(&self) -> Dims {
        self.grid.dims()
    }

    fn require_cubes(&self) -> Result<(), StaticError> {
        for (name, present) in [
            (PORO, self.poro.is_some()),
            (NTG, self.ntg.is_some()),
            (SW, self.sw.is_some()),
        ] {
            if !present {
                return Err(StaticError::InvalidInput(format!(
                    "grid is missing required property '{name}'"
                )));
            }
        }
        Ok(())
    }

    type Slab<'s>
        = GridSlab<'s>
    where
        Self: 's;

    fn slab(&self, k: usize) -> Result<GridSlab<'_>, StaticError> {
        let d = self.grid.dims();
        Ok(GridSlab {
            grid: self.grid,
            poro: self.poro,
            ntg: self.ntg,
            sw: self.sw,
            kbase: k * d.ni * d.nj,
        })
    }
}

impl CellSlab for GridSlab<'_> {
    #[inline]
    fn centroid_z(&self, local: usize) -> f64 {
        self.grid.cell_centroid_z_at(self.kbase + local)
    }
    #[inline]
    fn cell_volume(&self, local: usize) -> f64 {
        self.grid.cell_volume_at(self.kbase + local)
    }
    #[inline]
    fn poro(&self, local: usize) -> f64 {
        self.poro.expect("require_cubes verified PORO")[self.kbase + local]
    }
    #[inline]
    fn ntg(&self, local: usize) -> f64 {
        self.ntg.expect("require_cubes verified NTG")[self.kbase + local]
    }
    #[inline]
    fn sw(&self, local: usize) -> f64 {
        self.sw.expect("require_cubes verified SW")[self.kbase + local]
    }
}
