//! The out-of-core backing-storage mode (rulings **R1/R2/R4**,
//! `petekSuite/dev-docs/designs/out-of-core-strategy.md`).
//!
//! A `StaticModel` (up in `srs-model`) above the memory budget spills its
//! **geometry (ZCORN)** and **property cubes** onto a petekTools `store` — a
//! chunked, memory-mapped, k-slab-major lane file (petekTools ruling R1) — and
//! reads them back through **windowed mmap views**. COORD (the pillar lattice) is
//! `O(area)` and stays resident; only the per-cell arrays spill. Storage lanes are
//! **f32** at spill scale (R4: halves the bytes — the out-of-core enabler and the
//! MC bandwidth lever), so in-core↔spilled parity is TOLERANCE-based (~1e-7 on
//! volumes); **accumulations stay f64** regardless (R4 honesty clause).
//!
//! The spilled surface reproduces the in-core volumetrics loops
//! (`crate::volumetrics::grv`) **streaming**, one k-slab at a time: peak working set
//! is `O(ncol·nrow)`, not `O(grid)` (R2). The lane layout maps the grid's natural
//! k-slab-major order exactly (`zcorn[k·ni·nj·8 ..]`, `cube[k·ni·nj ..]`), so a
//! k-window is one contiguous zero-copy slice.

use crate::error::StaticError;
use crate::grid::{hexahedron_volume, Dims, Grid, Pillar, Point3};
use crate::volumetrics::{compute_clipped, CellSlab, Clip, SlabSource, NTG, PORO, SW};
use petektools::store::{Dtype, LaneSpec, Store, StoreSchema, StoreWriter};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// The ZCORN lane name (8 corner depths per cell, k-slab-major).
const ZCORN: &str = "ZCORN";
/// The COORD flat lane name (pillar lattice: top+bottom xyz per pillar).
const COORD: &str = "COORD";
/// Elements per pillar in the flat COORD lane: `top{x,y,z}` + `bottom{x,y,z}`.
const COORD_ELEMS_PER_PILLAR: usize = 6;

/// A unique-per-process spill file name (pid + a monotonic counter + nanos), so
/// concurrent MC workers and repeated builds never collide on one temp path.
///
/// `pub` for the `srs-model` builder/MC loops that pick a per-shard spill path;
/// an internal helper of the out-of-core mode, not a headline API.
pub fn unique_spill_path(dir: &Path) -> PathBuf {
    static SEQ: AtomicU64 = AtomicU64::new(0);
    let seq = SEQ.fetch_add(1, Ordering::Relaxed);
    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    dir.join(format!("petekstatic-spill-{pid}-{seq}-{nanos}.pts"))
}

/// The property cubes a spilled model always carries (the volumetrics inputs).
/// Extra cubes are spilled too but volumetrics reads exactly these three.
fn cube_lane_names(grid: &Grid) -> Vec<String> {
    let mut names: Vec<String> = grid.properties().names().map(str::to_owned).collect();
    names.sort(); // deterministic lane order → byte-deterministic store layout
    names
}

/// Spill a fully-built in-core `grid` to a petekTools store at a chosen location
/// under `dir` (unique file name), returning the mmap-backed [`SpillBacking`].
///
/// Writes ZCORN (f32 slab lane) + every property cube (f32 slab lane) + COORD
/// (f64 flat lane) k-slab-by-k-slab, then `finalize`s the seal and re-`open`s the
/// file read-only. The caller then drops the in-core grid buffers; the model's
/// steady-state resident set is `O(slab)` (the mmap pages in on demand).
///
/// # Errors
/// [`StaticError::InvalidInput`] on a non-vertical lattice (unsupported by spill
/// v1) or a missing volumetrics cube; [`StaticError::Algo`] on any store I/O.
pub fn spill_grid(grid: &Grid, dir: &Path, cleanup: bool) -> Result<SpillBacking, StaticError> {
    std::fs::create_dir_all(dir).map_err(|e| StaticError::Algo(e.into()))?;
    let path = unique_spill_path(dir);
    spill_grid_to(grid, &path, cleanup)
}

/// Spill `grid` to an **exact** store path (overwriting it), returning the
/// mmap-backed [`SpillBacking`]. The reused-path variant behind [`spill_grid`]:
/// the structured-MC spilled loop overwrites one per-shard store per draw (never
/// a new file per draw — ruling R3's "MC never spills per-draw state" = no
/// per-draw accumulation on disk).
///
/// # Errors
/// As [`spill_grid`].
pub fn spill_grid_to(grid: &Grid, path: &Path, cleanup: bool) -> Result<SpillBacking, StaticError> {
    let geom = grid.geometry();
    if !geom.is_vertical() {
        return Err(StaticError::InvalidInput(
            "out-of-core spill supports vertical lattices only (every grid today)".into(),
        ));
    }
    let dims = grid.dims();
    let per_slab_zcorn = dims.ni * dims.nj * 8;
    let per_slab_cube = dims.ni * dims.nj;
    let zcorn = geom.zcorn();
    let cube_names = cube_lane_names(grid);
    // Cube value slices, borrowed once (the write below is a slab-major stream).
    let cubes: Vec<(&str, &[f64])> = cube_names
        .iter()
        .map(|n| {
            (
                n.as_str(),
                grid.properties()
                    .get(n)
                    .map(|p| p.values.as_slice())
                    .expect("named cube present"),
            )
        })
        .collect();

    // Reuse the single streaming writer: the built grid's ZCORN/cube slabs are the
    // slab producers (narrowed f64 -> f32). Identical bytes to the old inline path.
    spill_streaming(
        path,
        dims,
        geom.coord(),
        &cube_names,
        cleanup,
        |k, out| {
            let z0 = k * per_slab_zcorn;
            for (dst, src) in out.iter_mut().zip(&zcorn[z0..z0 + per_slab_zcorn]) {
                *dst = *src as f32;
            }
            Ok(())
        },
        |name, k, out| {
            let vals = cubes
                .iter()
                .find(|(n, _)| *n == name)
                .map(|(_, v)| *v)
                .expect("named cube present");
            let c0 = k * per_slab_cube;
            for (dst, src) in out.iter_mut().zip(&vals[c0..c0 + per_slab_cube]) {
                *dst = *src as f32;
            }
            Ok(())
        },
    )
}

/// The store schema for a spilled `StaticModel`: an f32 ZCORN slab lane, an f64
/// COORD flat lane (the resident pillar lattice), and one f32 slab lane per cube.
fn store_schema(dims: Dims, pillar_count: usize, cube_names: &[String]) -> StoreSchema {
    let per_slab_zcorn = dims.ni * dims.nj * 8;
    let per_slab_cube = dims.ni * dims.nj;
    let mut lanes = vec![
        LaneSpec::slab(ZCORN, Dtype::F32, per_slab_zcorn as u64),
        LaneSpec::flat(
            COORD,
            Dtype::F64,
            (pillar_count * COORD_ELEMS_PER_PILLAR) as u64,
        ),
    ];
    for name in cube_names {
        lanes.push(LaneSpec::slab(name, Dtype::F32, per_slab_cube as u64));
    }
    let app = serde_json::json!({
        "kind": "petekstatic-static-model",
        "ni": dims.ni, "nj": dims.nj, "nk": dims.nk,
        "cubes": cube_names,
    });
    StoreSchema::new(dims.nk as u64, lanes).with_app(app)
}

/// The ONE store-writing path: create the store, write COORD, then fill **each
/// k-slab in place** through the writer's zero-copy `slab_mut_f32` views —
/// `fill_zcorn(k, &mut[f32])` for the ZCORN slab and `fill_cube(name, k, &mut[f32])`
/// per cube slab — then finalize and re-open read-only, returning the mmap-backed
/// [`SpillBacking`]. Peak writer working set is **one slab** (the mmap pages are
/// written directly), so a **slab-incremental build** that produces geometry +
/// cubes on demand achieves `O(slab)` build peak (never a whole in-core grid). The
/// post-hoc [`spill_grid_to`] rides the same path (its producer reads a built grid).
///
/// # Errors
/// [`StaticError::Algo`] on any store I/O; propagates a producer error.
#[allow(clippy::too_many_arguments)]
pub fn spill_streaming(
    path: &Path,
    dims: Dims,
    coord: &[Pillar],
    cube_names: &[String],
    cleanup: bool,
    mut fill_zcorn: impl FnMut(usize, &mut [f32]) -> Result<(), StaticError>,
    mut fill_cube: impl FnMut(&str, usize, &mut [f32]) -> Result<(), StaticError>,
) -> Result<SpillBacking, StaticError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| StaticError::Algo(e.into()))?;
    }
    let schema = store_schema(dims, coord.len(), cube_names);
    let mut w = StoreWriter::create(path, schema)?;

    // COORD flat lane: pack each pillar as top{x,y,z}, bottom{x,y,z} in place.
    {
        let flat = w.flat_mut_f64(COORD)?;
        for (p, slot) in coord
            .iter()
            .zip(flat.chunks_exact_mut(COORD_ELEMS_PER_PILLAR))
        {
            slot.copy_from_slice(&[
                p.top.x, p.top.y, p.top.z, p.bottom.x, p.bottom.y, p.bottom.z,
            ]);
        }
    }

    for k in 0..dims.nk {
        fill_zcorn(k, w.slab_mut_f32(ZCORN, k as u64)?)?;
        for name in cube_names {
            fill_cube(name, k, w.slab_mut_f32(name, k as u64)?)?;
        }
    }
    w.finalize()?;

    let store = petektools::store::open(path)?;
    Ok(SpillBacking {
        store,
        path: path.to_path_buf(),
        cleanup,
        dims,
        coord: coord.to_vec(),
        cube_names: cube_names.to_vec(),
    })
}

/// The mmap-backed store behind a spilled `StaticModel` (up in `srs-model`) —
/// the out-of-core "backing-storage mode". Holds the read-only store (ZCORN +
/// cubes, f32) plus the resident COORD lattice. All heavy reads window the store.
///
/// **Drop semantics (R5):** unless [`SpillBacking::detach`]ed, the store file is
/// removed on drop (temp cleanup). Detach to keep it (a caller-owned location).
#[derive(Debug)]
pub struct SpillBacking {
    store: Store,
    path: PathBuf,
    cleanup: bool,
    dims: Dims,
    coord: Vec<Pillar>,
    cube_names: Vec<String>,
}

impl Drop for SpillBacking {
    fn drop(&mut self) {
        if self.cleanup {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

impl SpillBacking {
    /// Grid dimensions of the spilled model.
    #[must_use]
    pub fn dims(&self) -> Dims {
        self.dims
    }

    /// The on-disk store path.
    #[must_use]
    pub fn store_path(&self) -> &Path {
        &self.path
    }

    /// Retain the store file past drop (a caller-owned location) — the "detach"
    /// side of the cleanup contract.
    pub fn detach(&mut self) {
        self.cleanup = false;
    }

    /// The names of the spilled property cubes.
    #[must_use]
    pub fn cube_names(&self) -> &[String] {
        &self.cube_names
    }

    /// One ZCORN k-slab (`ni·nj·8` f32, zero-copy mmap view).
    fn zcorn_slab(&self, k: usize) -> Result<&[f32], StaticError> {
        Ok(self.store.slab_f32(ZCORN, k as u64)?)
    }

    /// One property-cube k-slab (`ni·nj` f32, zero-copy), or `None` if the cube is
    /// not spilled (the bulk path reads no cubes).
    fn cube_slab_opt(&self, name: &str, k: usize) -> Result<Option<&[f32]>, StaticError> {
        if self.cube_names.iter().any(|n| n == name) {
            Ok(Some(self.store.slab_f32(name, k as u64)?))
        } else {
            Ok(None)
        }
    }

    /// A [`SlabSource`] view over this backing — the streaming f32 mirror of
    /// `crate::volumetrics::GridSource`. The whole spilled volumetric surface
    /// (`bulk_volume`, the model's `in_place*` / `in_place_by_zone`) runs through
    /// the ONE unified core (`crate::volumetrics::compute_clipped`) over this source,
    /// so there is no hand-maintained twin of the GRV/HCPV loop.
    #[must_use]
    pub fn source(&self) -> SpillSource<'_> {
        SpillSource { backing: self }
    }

    /// Total gross rock volume \[m³\] streamed slab-by-slab — the spilled
    /// `StaticModel::bulk_volume` (up in `srs-model`).
    ///
    /// # Errors
    /// [`StaticError::Algo`] on a store read failure.
    pub fn bulk_volume(&self) -> Result<f64, StaticError> {
        Ok(compute_clipped(&self.source(), Clip::Bulk, 0..self.dims.nk, false)?.grv_m3)
    }

    /// Reconstruct a **whole in-core [`Grid`]** from the mmap store — geometry
    /// (COORD + ZCORN, widened f32→f64) plus every spilled property cube — for the
    /// **non-hot-path view/shell exports**, which need random cell access the
    /// streaming volumetric surface does not provide. This allocates an `O(grid)`
    /// transient (the very thing out-of-core avoids on the MC hot path), so it is
    /// only for one-shot exports (`StaticModel::volume_bundle` on a spilled model),
    /// never the realization loop. Geometry carries the store's f32 precision (the
    /// documented ~1e-7 in-core↔spilled parity delta).
    ///
    /// `pub` for the `srs-model` view/shell exporters that need a random-access grid;
    /// an internal helper of the out-of-core mode, not a headline API.
    ///
    /// # Errors
    /// [`StaticError`] on a store read or grid-construction failure.
    pub fn to_in_core_grid(&self) -> Result<Grid, StaticError> {
        let dims = self.dims;
        let mut grid = crate::grid::build_box(crate::grid::BoxSpec::square(
            1.0,
            1.0,
            Dims::new(1, 1, 1).expect("1x1x1 dims are valid"),
        ))?;
        // ZCORN, k-slab-major (matches the in-core layout the spill was written from).
        let mut zcorn = Vec::with_capacity(dims.cell_count() * 8);
        for k in 0..dims.nk {
            zcorn.extend(self.zcorn_slab(k)?.iter().map(|&z| z as f64));
        }
        grid.install_geometry(dims, self.coord.clone(), zcorn);
        // Every spilled cube, k-slab-major.
        for name in self.cube_names().to_vec() {
            let mut values = Vec::with_capacity(dims.cell_count());
            for k in 0..dims.nk {
                let slab = self
                    .cube_slab_opt(&name, k)?
                    .expect("cube listed in cube_names is spilled");
                values.extend(slab.iter().map(|&v| v as f64));
            }
            grid.properties_mut()
                .set(crate::grid::Property { name, values })?;
        }
        Ok(grid)
    }
}

/// The spilled [`SlabSource`]: streams one f32 k-slab of the mmap store at a time,
/// widening geometry + cubes to f64 at the cell boundary (accumulations stay f64,
/// R4). The f32 narrowing of ZCORN is the sole parity delta vs the in-core path.
pub struct SpillSource<'a> {
    backing: &'a SpillBacking,
}

impl SlabSource for SpillSource<'_> {
    fn dims(&self) -> Dims {
        self.backing.dims
    }

    fn require_cubes(&self) -> Result<(), StaticError> {
        for name in [PORO, NTG, SW] {
            if !self.backing.cube_names.iter().any(|n| n == name) {
                return Err(StaticError::InvalidInput(format!(
                    "spilled grid is missing required property '{name}'"
                )));
            }
        }
        Ok(())
    }

    type Slab<'s>
        = SpillSlab<'s>
    where
        Self: 's;

    fn slab(&self, k: usize) -> Result<SpillSlab<'_>, StaticError> {
        Ok(SpillSlab {
            zslab: self.backing.zcorn_slab(k)?,
            poro: self.backing.cube_slab_opt(PORO, k)?,
            ntg: self.backing.cube_slab_opt(NTG, k)?,
            sw: self.backing.cube_slab_opt(SW, k)?,
            coord: &self.backing.coord,
            dims: self.backing.dims,
        })
    }
}

/// One f32 k-slab of a [`SpillSource`]: the zero-copy ZCORN + cube mmap windows for
/// slab `k` plus the resident pillar tops. Cell accessors widen to f64.
pub struct SpillSlab<'a> {
    zslab: &'a [f32],
    poro: Option<&'a [f32]>,
    ntg: Option<&'a [f32]>,
    sw: Option<&'a [f32]>,
    coord: &'a [Pillar],
    dims: Dims,
}

impl SpillSlab<'_> {
    /// The 8 corners of slab-local cell `local` from its f32 ZCORN + the resident
    /// pillar tops — bit-parallel to `CornerPointGeom::cell_corners_at` on a
    /// vertical lattice, except ZCORN is narrowed to f32 (the sole parity delta).
    #[inline]
    fn corners(&self, local: usize) -> [Point3; 8] {
        let ni = self.dims.ni;
        let j = local / ni;
        let i = local % ni;
        let z = &self.zslab[local * 8..local * 8 + 8];
        let mut corners = [Point3::new(0.0, 0.0, 0.0); 8];
        for (idx, slot) in corners.iter_mut().enumerate() {
            let di = idx & 1;
            let dj = (idx >> 1) & 1;
            let top = self.coord[self.dims.pillar_linear(i + di, j + dj)].top;
            *slot = Point3::new(top.x, top.y, z[idx] as f64);
        }
        corners
    }
}

impl CellSlab for SpillSlab<'_> {
    /// Mean of the 8 f32 corner depths widened to f64 — matching
    /// `CornerPointGeom::centroid_z_at`.
    #[inline]
    fn centroid_z(&self, local: usize) -> f64 {
        let z = &self.zslab[local * 8..local * 8 + 8];
        let s: f64 = z.iter().map(|&v| v as f64).sum();
        s / 8.0
    }
    #[inline]
    fn cell_volume(&self, local: usize) -> f64 {
        hexahedron_volume(&self.corners(local))
    }
    #[inline]
    fn poro(&self, local: usize) -> f64 {
        self.poro.expect("require_cubes verified PORO")[local] as f64
    }
    #[inline]
    fn ntg(&self, local: usize) -> f64 {
        self.ntg.expect("require_cubes verified NTG")[local] as f64
    }
    #[inline]
    fn sw(&self, local: usize) -> f64 {
        self.sw.expect("require_cubes verified SW")[local] as f64
    }
}
