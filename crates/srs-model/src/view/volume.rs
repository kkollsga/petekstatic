//! [`VolumeBundle`] — the 3-D corner-point **exterior-shell** inspection bundle
//! ([`StaticModel::volume_bundle`]).
//!
//! Instead of the full 8-vertex / 12-triangle cell soup (which made the viewer
//! payload O(volume) and blew past V8's inline-script wall at ~0.9M cells), the
//! bundle ships only the **visible faces**: for each active cell, a cell face is
//! emitted iff the neighbour across it is inactive/absent or out of bounds. For a
//! low-relief reservoir grid this is O(surface) — a full box collapses to exactly
//! its outer shell (~`2·ni·nj` caps + skin instead of `ni·nj·nk·12` triangles).
//!
//! Shared vertices are deduplicated by f32 position, and the per-cell arrays are
//! **compacted to the shell cells only** (a `tri_cell` index per triangle recovers
//! cell identity for the viewer's threshold filter / picking). The arrays serialize
//! as raw little-endian binary blocks — see [`super::wire`] for the decode spec.
//!
//! ## Threshold semantics (the honest MVP — `task_suite_bundle_binary` decision (b))
//! The default bundle is the shell of the **full active set**. A client-side
//! threshold slider can only *hide* shell cells it already holds — it cannot expose
//! a cell's interior faces, because those triangles are not in the payload. So the
//! viewer's slider is documented as a **shell-only** filter. To truly re-cut the
//! shell at a property cutoff (exposing the interior that the cutoff reveals), the
//! server regenerates via [`StaticModel::volume_bundle_thresholded`] — live viewers
//! re-request it; a `save_view` can pre-compute a few cutoff steps.
//!
//! ## Mesh home (layer charter, `decision_layer_charters`)
//! The corner-point mesh builder **lives here** — moved DOWN into the GEOMODEL
//! layer from petekSim's `srs-core/mesh.rs`; the DAG flows downward, so the mesh
//! belongs beside the grid it meshes.

use super::frame::ValueRange;
use super::wire::{self, Block, BlockData, Head};
use super::SCHEMA_VERSION;
use crate::model::StaticModel;
use petekstatic_error::StaticError;
use serde::{Deserialize, Serialize};
use srs_grid::{Dims, Ijk, Point3};
use std::collections::HashMap;
use std::io::{self, Write};

/// The 6 quad faces of a cell as corner indices (corner = `di + 2*dj + 4*dk`),
/// each split into two triangles at emit time. Aligned index-for-index with
/// [`NEIGHBORS`] (the `(di,dj,dk)` step to the cell sharing that face).
const FACES: [[usize; 4]; 6] = [
    [0, 1, 3, 2], // k- (top)
    [4, 5, 7, 6], // k+ (bottom)
    [0, 1, 5, 4], // j-
    [2, 3, 7, 6], // j+
    [0, 2, 6, 4], // i-
    [1, 3, 7, 5], // i+
];

/// The neighbouring cell across each face in [`FACES`], as a signed `(di,dj,dk)`
/// step. A face is emitted iff this neighbour is out of bounds or not present.
const NEIGHBORS: [(isize, isize, isize); 6] = [
    (0, 0, -1), // k-
    (0, 0, 1),  // k+
    (0, -1, 0), // j-
    (0, 1, 0),  // j+
    (-1, 0, 0), // i-
    (1, 0, 0),  // i+
];

/// The 3-D volume **exterior-shell** bundle: the visible cell faces coloured by a
/// property, with compact per-shell-cell values + zone ids. Positions are grid
/// **local** corner coordinates (metres; z positive-down depth) — the viewer
/// renders the volume in its own space and needs no world georeference.
///
/// Layout (all indices/counts derive from the arrays):
/// - `positions`: deduplicated shell vertices, 3 `f32` `(x,y,z)` each.
/// - `indices`: 3 `u32` per triangle, into `positions`.
/// - `tri_cell`: one `u32` per triangle — the **compact** cell index (into
///   `cell_values` / `zone_ids`) the triangle belongs to (viewer picking / filter).
/// - `cell_values`: one `f32` per shell cell (compact order), the property value.
/// - `zone_ids`: one `u16` per shell cell (compact order), index into `zone_names`.
///
/// The bundle serializes via [`Self::write_self_contained`] (base64, one file) or
/// [`Self::write_sidecar`] (JSON envelope + `model.bin`); the `serde` derive is a
/// convenience for Rust-side round-trips, not the viewer wire format.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VolumeBundle {
    pub schema_version: u32,
    pub inputs_ref: String,
    /// The property carried in `cell_values`.
    pub property: String,
    /// Total grid cell count (reference; `cell_values.len()` is the shell count).
    pub cell_count: usize,
    /// Deduplicated shell-vertex positions, 3 floats `(x, y, z)` per vertex.
    /// Length `vertex_count * 3`.
    pub positions: Vec<f32>,
    /// Triangle indices into the vertex list, 3 per triangle. Length
    /// `triangle_count * 3`.
    pub indices: Vec<u32>,
    /// Per-triangle **compact** cell index (into `cell_values` / `zone_ids`).
    /// Length `triangle_count` (both triangles of a quad face share it).
    pub tri_cell: Vec<u32>,
    /// Per-shell-cell property value in compact (first-emitted) order — the array
    /// the viewer's threshold filter reads. Length = shell cell count. `f32`
    /// (display precision); `NaN`-valued cells never reach the shell.
    pub cell_values: Vec<f32>,
    /// Per-shell-cell zone id (index into `zone_names`), compact order. Length =
    /// shell cell count.
    pub zone_ids: Vec<u16>,
    /// Zone names indexed by `zone_ids`.
    pub zone_names: Vec<String>,
    /// Value span over the shell cell values (colour-ramp legend).
    pub value_range: ValueRange,
}

impl StaticModel {
    /// Export the 3-D corner-point **exterior shell** ([`VolumeBundle`]) coloured
    /// by `property`. Only faces bordering an inactive/absent neighbour or the grid
    /// boundary are emitted; a cell is *active* iff its property value is finite and
    /// its bulk volume is positive (a collapsed / pinched cell is absent).
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if `property` is absent.
    pub fn volume_bundle(&self, property: &str) -> Result<VolumeBundle, StaticError> {
        self.extract_shell(property, |v| v.is_finite())
    }

    /// Re-cut the exterior shell treating only cells whose property value passes the
    /// `cutoff` as active — the server-side threshold regeneration path (decision
    /// (b)): hiding sub-cutoff cells exposes the interior faces the cutoff reveals.
    /// `keep_above == true` keeps cells with `value >= cutoff`, else `value <= cutoff`.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if `property` is absent.
    pub fn volume_bundle_thresholded(
        &self,
        property: &str,
        cutoff: f64,
        keep_above: bool,
    ) -> Result<VolumeBundle, StaticError> {
        self.extract_shell(property, |v| {
            v.is_finite() && if keep_above { v >= cutoff } else { v <= cutoff }
        })
    }

    /// Shared shell extraction. `keep` decides, from a cell's property value, whether
    /// the cell is a candidate active cell (it must additionally have positive bulk
    /// volume — a collapsed cell is never present, so its neighbours expose faces).
    fn extract_shell(
        &self,
        property: &str,
        keep: impl Fn(f64) -> bool,
    ) -> Result<VolumeBundle, StaticError> {
        // The real grid whether in-core or spilled: a spilled (out-of-core) model's
        // `grid()` is a 1×1×1 placeholder with no cubes, so the shell must materialize
        // the backing — otherwise every large (spilled) model exported an EMPTY shell
        // ("no property" / 0 tris in the viewer; `question_volume_bundle_stack_empty`).
        let grid = self.view_grid()?;
        let dims = grid.dims();
        let cell_count = dims.cell_count();
        let cube = &grid
            .properties()
            .get(property)
            .ok_or_else(|| {
                StaticError::InvalidInput(format!("volume_bundle: no property '{property}'"))
            })?
            .values;

        // k -> zone id lookup (zones partition [0, nk)).
        let zones = self.zones().zones();
        let zone_names: Vec<String> = zones.iter().map(|z| z.name.clone()).collect();
        let mut zone_of_k = vec![0u16; dims.nk];
        for (zid, z) in zones.iter().enumerate() {
            for k in z.k_range.clone() {
                zone_of_k[k] = zid as u16;
            }
        }

        // Presence mask (one pass): active iff the property value passes `keep`
        // and the cell has positive bulk volume (collapsed cells are absent).
        let present: Vec<bool> = (0..cell_count)
            .map(|c| keep(cube[c]) && grid.cell_volume_at(c) > 0.0)
            .collect();

        let shell = extract_faces(
            dims,
            &present,
            |c| grid.cell(ijk_of(c, dims)).corners,
            |c| cube[c] as f32,
            |c| zone_of_k[ijk_of(c, dims).k],
        );

        let value_range = ValueRange::of(shell.cell_values.iter().map(|&v| v as f64));

        Ok(VolumeBundle {
            schema_version: SCHEMA_VERSION,
            inputs_ref: self.provenance().inputs_ref.clone(),
            property: property.to_string(),
            cell_count,
            positions: shell.positions,
            indices: shell.indices,
            tri_cell: shell.tri_cell,
            cell_values: shell.cell_values,
            zone_ids: shell.zone_ids,
            zone_names,
            value_range,
        })
    }
}

/// `(i,j,k)` of a linear cell index (`i` fastest), for the extraction closures.
#[inline]
fn ijk_of(c: usize, dims: Dims) -> Ijk {
    let per_layer = dims.ni * dims.nj;
    let k = c / per_layer;
    let rem = c % per_layer;
    Ijk::new(rem % dims.ni, rem / dims.ni, k)
}

/// The raw shell arrays produced by [`extract_faces`].
struct Shell {
    positions: Vec<f32>,
    indices: Vec<u32>,
    tri_cell: Vec<u32>,
    cell_values: Vec<f32>,
    zone_ids: Vec<u16>,
}

/// The pure exterior-shell extraction: over the presence mask, emit each active
/// cell's faces whose neighbour is absent or out of bounds, deduplicating vertices
/// by f32 position and compacting the per-cell arrays to the shell cells only.
/// Driven by closures so it is testable off any geometry — `corners(c)` gives a
/// cell's 8 corners (order `di + 2*dj + 4*dk`), `value(c)` / `zone(c)` its scalars.
fn extract_faces(
    dims: Dims,
    present: &[bool],
    corners: impl Fn(usize) -> [Point3; 8],
    value: impl Fn(usize) -> f32,
    zone: impl Fn(usize) -> u16,
) -> Shell {
    let (ni, nj, nk) = (dims.ni, dims.nj, dims.nk);
    let mut positions: Vec<f32> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let mut tri_cell: Vec<u32> = Vec::new();
    let mut cell_values: Vec<f32> = Vec::new();
    let mut zone_ids: Vec<u16> = Vec::new();
    // Dedup shell vertices by exact f32 position bits (correct for conformal AND
    // split-node geometry: identical position => one vertex; colour is per-triangle
    // via `tri_cell`, so geometric merging is safe).
    let mut vert_map: HashMap<[u32; 3], u32> = HashMap::new();

    for k in 0..nk {
        for j in 0..nj {
            for i in 0..ni {
                let c = (k * nj + j) * ni + i;
                if !present[c] {
                    continue;
                }
                let mut cid = u32::MAX; // compact id, assigned on the first face
                let cell = corners(c);
                for (f, face) in FACES.iter().enumerate() {
                    if neighbor_present(present, dims, (i, j, k), NEIGHBORS[f]) {
                        continue; // interior shared face — vanishes
                    }
                    if cid == u32::MAX {
                        cid = cell_values.len() as u32;
                        cell_values.push(value(c));
                        zone_ids.push(zone(c));
                    }
                    let mut v = [0u32; 4];
                    for (n, &corner) in face.iter().enumerate() {
                        let p = cell[corner];
                        let (x, y, z) = (p.x as f32, p.y as f32, p.z as f32);
                        let key = [x.to_bits(), y.to_bits(), z.to_bits()];
                        v[n] = *vert_map.entry(key).or_insert_with(|| {
                            let id = (positions.len() / 3) as u32;
                            positions.extend_from_slice(&[x, y, z]);
                            id
                        });
                    }
                    // Two triangles per quad, keeping the face's winding.
                    indices.extend_from_slice(&[v[0], v[1], v[2], v[0], v[2], v[3]]);
                    tri_cell.push(cid);
                    tri_cell.push(cid);
                }
            }
        }
    }
    Shell {
        positions,
        indices,
        tri_cell,
        cell_values,
        zone_ids,
    }
}

/// Is the neighbour of cell `(i,j,k)` across step `(di,dj,dk)` in bounds AND
/// present? Out of bounds (grid boundary) counts as absent → the face is emitted.
#[inline]
fn neighbor_present(
    present: &[bool],
    dims: Dims,
    (i, j, k): (usize, usize, usize),
    (di, dj, dk): (isize, isize, isize),
) -> bool {
    let (ni, nj, nk) = (dims.ni as isize, dims.nj as isize, dims.nk as isize);
    let (ii, jj, kk) = (i as isize + di, j as isize + dj, k as isize + dk);
    if ii < 0 || ii >= ni || jj < 0 || jj >= nj || kk < 0 || kk >= nk {
        return false;
    }
    let c = (kk * nj + jj) * ni + ii;
    present[c as usize]
}

impl VolumeBundle {
    /// The five binary blocks in wire order (positions, indices, tri_cell,
    /// cell_values, zone_ids). See [`super::wire`] for the decode spec.
    fn blocks(&self) -> [Block<'_>; 5] {
        let vcount = self.positions.len() / 3;
        let tcount = self.indices.len() / 3;
        let scount = self.cell_values.len();
        [
            Block {
                name: "positions",
                shape: vec![vcount, 3],
                data: BlockData::F32(&self.positions),
            },
            Block {
                name: "indices",
                shape: vec![tcount, 3],
                data: BlockData::U32(&self.indices),
            },
            Block {
                name: "tri_cell",
                shape: vec![tcount],
                data: BlockData::U32(&self.tri_cell),
            },
            Block {
                name: "cell_values",
                shape: vec![scount],
                data: BlockData::F32(&self.cell_values),
            },
            Block {
                name: "zone_ids",
                shape: vec![scount],
                data: BlockData::U16(&self.zone_ids),
            },
        ]
    }

    /// Envelope head fields borrowed for the streaming writers.
    fn head(&self) -> Head<'_> {
        Head {
            schema_version: self.schema_version,
            inputs_ref: &self.inputs_ref,
            property: &self.property,
            cell_count: self.cell_count,
            shell_cell_count: self.cell_values.len(),
            vertex_count: self.positions.len() / 3,
            triangle_count: self.indices.len() / 3,
            zone_names: &self.zone_names,
            value_range: &self.value_range,
        }
    }

    /// Stream the **self-contained** payload (metadata envelope + base64-wrapped
    /// binary blocks) to `w` — one file, the `save_view` export. Peak memory stays
    /// ~1x the payload (streamed base64, no whole-string intermediate).
    pub fn write_self_contained<W: Write>(&self, w: &mut W) -> io::Result<()> {
        wire::write_self_contained(&self.head(), &self.blocks(), w)
    }

    /// Stream the **served** payload: the JSON envelope (with a `(offset,length)`
    /// manifest) to `json` and the raw `model.bin` block bytes to `bin` — no base64
    /// overhead when HTTP serves the two files.
    pub fn write_sidecar<W1: Write, W2: Write>(
        &self,
        json: &mut W1,
        bin: &mut W2,
    ) -> io::Result<()> {
        wire::write_sidecar(&self.head(), &self.blocks(), json, bin)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BuildOpts, ConstantPriors, StaticModelBuilder};
    use base64::engine::general_purpose::STANDARD;
    use base64::Engine;
    use srs_gridder::{Conformity, SolveOpts};

    // ni×nj×nk flat box populated with constant priors (PORO = 0.25 everywhere).
    fn model_nk(ni: usize, nj: usize, nk: usize) -> StaticModel {
        let opts = BuildOpts {
            area_m2: 10_000.0,
            gross_height_m: 40.0,
            nk,
            conformity: Conformity::Proportional,
            solve_opts: SolveOpts::default(),
            priors: ConstantPriors {
                porosity: 0.25,
                net_to_gross: 0.8,
                water_saturation: 0.3,
            },
        };
        StaticModelBuilder::flat(ni, nj, 5000.0, 5100.0, opts)
            .unwrap()
            .build()
            .unwrap()
    }

    fn model() -> StaticModel {
        model_nk(3, 2, 4)
    }

    // Expected outer-shell triangle count of a solid ni×nj×nk box: 2 triangles per
    // face, 2 faces on each of the 3 axis pairs.
    fn box_shell_tris(ni: usize, nj: usize, nk: usize) -> usize {
        2 * (ni * nj + ni * nk + nj * nk) * 2
    }

    #[test]
    fn full_box_is_exactly_the_outer_shell() {
        let (ni, nj, nk) = (3, 2, 4);
        let b = model_nk(ni, nj, nk).volume_bundle("PORO").unwrap();
        assert_eq!(b.cell_count, ni * nj * nk);
        // Only the outer skin: interior shared faces vanish.
        assert_eq!(b.indices.len() / 3, box_shell_tris(ni, nj, nk));
        assert_eq!(b.tri_cell.len(), b.indices.len() / 3);
        // A watertight box hull: unique vertices == (ni+1)(nj+1)(nk+1) corner
        // lattice minus the interior lattice nodes no shell face touches. Just
        // assert the dedup did real work (far fewer than 4 verts/triangle).
        let vcount = b.positions.len() / 3;
        assert!(
            vcount < b.indices.len(),
            "vertices deduplicated across faces"
        );
        // Every index references an existing vertex.
        assert!(b.indices.iter().all(|&i| (i as usize) < vcount));
    }

    #[test]
    fn shell_cells_carry_their_value_and_zone() {
        let b = model().volume_bundle("PORO").unwrap();
        // Constant priors: every shell cell value is 0.25.
        assert!(b.cell_values.iter().all(|&v| (v - 0.25).abs() < 1e-6));
        assert_eq!(b.zone_ids.len(), b.cell_values.len());
        assert!(b.zone_ids.iter().all(|&z| z == 0));
        assert_eq!(b.zone_names, vec!["RESERVOIR".to_string()]);
        assert!((b.value_range.min - 0.25).abs() < 1e-6);
        assert!((b.value_range.max - 0.25).abs() < 1e-6);
        // tri_cell indexes the compact arrays.
        assert!(b
            .tri_cell
            .iter()
            .all(|&c| (c as usize) < b.cell_values.len()));
    }

    // Unit-cube corners for cell (i,j,k), order `di + 2*dj + 4*dk` — a trivial
    // analytic geometry to exercise the pure `extract_faces` face counting.
    fn unit_cube(dims: Dims) -> impl Fn(usize) -> [Point3; 8] {
        move |c| {
            let ijk = ijk_of(c, dims);
            let mut cn = [Point3::new(0.0, 0.0, 0.0); 8];
            for (idx, slot) in cn.iter_mut().enumerate() {
                let (di, dj, dk) = (idx & 1, (idx >> 1) & 1, (idx >> 2) & 1);
                *slot = Point3::new(
                    (ijk.i + di) as f64,
                    (ijk.j + dj) as f64,
                    (ijk.k + dk) as f64,
                );
            }
            cn
        }
    }

    #[test]
    fn extract_faces_full_box_is_the_outer_shell_exactly() {
        let (ni, nj, nk) = (3, 3, 3);
        let dims = Dims::new(ni, nj, nk).unwrap();
        let present = vec![true; ni * nj * nk];
        let s = extract_faces(dims, &present, unit_cube(dims), |_| 1.0, |_| 0);
        assert_eq!(s.indices.len() / 3, box_shell_tris(ni, nj, nk)); // 108
                                                                     // The hull dedups to the 4×4×4 corner lattice minus the 2×2×2 fully
                                                                     // interior nodes no shell face touches: 64 - 8 = 56.
        assert_eq!(s.positions.len() / 3, 56);
        assert_eq!(s.cell_values.len(), 26); // every cell but the centre is on the skin
    }

    #[test]
    fn collapsed_interior_cell_exposes_its_neighbours_faces() {
        // 3×3×3 box with the dead-centre cell absent: its 6 neighbours each gain
        // the one face that was shared with it -> +6 faces = +12 triangles, and the
        // centre contributes no shell cell (26 of 27 cells on the shell).
        let (ni, nj, nk) = (3, 3, 3);
        let dims = Dims::new(ni, nj, nk).unwrap();
        let mut present = vec![true; ni * nj * nk];
        let centre = (nj + 1) * ni + 1; // (i,j,k) = (1,1,1)
        present[centre] = false;
        let s = extract_faces(dims, &present, unit_cube(dims), |_| 1.0, |_| 0);
        assert_eq!(s.indices.len() / 3, box_shell_tris(ni, nj, nk) + 6 * 2); // 120
        assert_eq!(s.cell_values.len(), 26);
        assert_eq!(s.tri_cell.len(), s.indices.len() / 3);
        assert!(s
            .tri_cell
            .iter()
            .all(|&c| (c as usize) < s.cell_values.len()));
    }

    #[test]
    fn thresholded_recut_gates_cells_by_cutoff() {
        // Constant 0.25 model: a cutoff above the value drops every cell (empty
        // shell); a cutoff below keeps them all (== the full shell). This exercises
        // the server-side regeneration predicate end to end.
        let m = model_nk(4, 4, 4);
        let full = m.volume_bundle("PORO").unwrap();
        let all = m.volume_bundle_thresholded("PORO", 0.20, true).unwrap();
        assert_eq!(all.indices.len(), full.indices.len());
        assert_eq!(all.cell_values.len(), full.cell_values.len());
        let none = m.volume_bundle_thresholded("PORO", 0.30, true).unwrap();
        assert!(none.indices.is_empty());
        assert!(none.cell_values.is_empty());
        // keep_below keeps the sub-cutoff set (here: everything at/under 0.30).
        let below = m.volume_bundle_thresholded("PORO", 0.30, false).unwrap();
        assert_eq!(below.cell_values.len(), full.cell_values.len());
    }

    #[test]
    fn missing_property_errors() {
        assert!(model().volume_bundle("NOPE").is_err());
        assert!(model()
            .volume_bundle_thresholded("NOPE", 0.0, true)
            .is_err());
    }

    // --- binary round-trip: encode -> decode == source arrays, bit-exact ---

    fn read_f32_le(bytes: &[u8]) -> Vec<f32> {
        bytes
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
            .collect()
    }
    fn read_u32_le(bytes: &[u8]) -> Vec<u32> {
        bytes
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
            .collect()
    }
    fn read_u16_le(bytes: &[u8]) -> Vec<u16> {
        bytes
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes(c.try_into().unwrap()))
            .collect()
    }

    #[test]
    fn self_contained_base64_round_trips_bit_exact() {
        let b = model_nk(4, 3, 5).volume_bundle("PORO").unwrap();
        let mut out = Vec::new();
        b.write_self_contained(&mut out).unwrap();
        let env: serde_json::Value = serde_json::from_slice(&out).unwrap();
        assert_eq!(env["schema_version"], serde_json::json!(5));
        assert_eq!(env["encoding"], serde_json::json!("base64"));
        let blk = &env["blocks"];
        let dec = |name: &str| {
            STANDARD
                .decode(blk[name]["data"].as_str().unwrap())
                .unwrap()
        };
        assert_eq!(read_f32_le(&dec("positions")), b.positions);
        assert_eq!(read_u32_le(&dec("indices")), b.indices);
        assert_eq!(read_u32_le(&dec("tri_cell")), b.tri_cell);
        assert_eq!(read_f32_le(&dec("cell_values")), b.cell_values);
        assert_eq!(read_u16_le(&dec("zone_ids")), b.zone_ids);
    }

    #[test]
    fn sidecar_manifest_round_trips_bit_exact() {
        let b = model_nk(4, 3, 5).volume_bundle("PORO").unwrap();
        let (mut json, mut bin) = (Vec::new(), Vec::new());
        b.write_sidecar(&mut json, &mut bin).unwrap();
        let env: serde_json::Value = serde_json::from_slice(&json).unwrap();
        assert_eq!(env["encoding"], serde_json::json!("sidecar"));
        let blk = &env["blocks"];
        let slice = |name: &str| {
            let off = blk[name]["offset"].as_u64().unwrap() as usize;
            let len = blk[name]["length"].as_u64().unwrap() as usize;
            &bin[off..off + len]
        };
        assert_eq!(read_f32_le(slice("positions")), b.positions);
        assert_eq!(read_u32_le(slice("indices")), b.indices);
        assert_eq!(read_u32_le(slice("tri_cell")), b.tri_cell);
        assert_eq!(read_f32_le(slice("cell_values")), b.cell_values);
        assert_eq!(read_u16_le(slice("zone_ids")), b.zone_ids);
        // The manifest exactly tiles the .bin (tight packing, no gaps).
        let total: usize = [
            "positions",
            "indices",
            "tri_cell",
            "cell_values",
            "zone_ids",
        ]
        .iter()
        .map(|n| blk[n]["length"].as_u64().unwrap() as usize)
        .sum();
        assert_eq!(total, bin.len());
    }

    #[test]
    fn schema_snapshot_v4_envelope() {
        let b = model().volume_bundle("PORO").unwrap();
        let mut out = Vec::new();
        b.write_self_contained(&mut out).unwrap();
        let env: serde_json::Value = serde_json::from_slice(&out).unwrap();
        let mut keys: Vec<&str> = env
            .as_object()
            .unwrap()
            .keys()
            .map(String::as_str)
            .collect();
        keys.sort_unstable();
        assert_eq!(
            keys,
            [
                "blocks",
                "cell_count",
                "encoding",
                "inputs_ref",
                "kind",
                "property",
                "schema_version",
                "shell_cell_count",
                "triangle_count",
                "value_range",
                "vertex_count",
                "zone_names",
            ]
        );
        assert_eq!(env["kind"], serde_json::json!("volume"));
        assert_eq!(env["schema_version"], serde_json::json!(5));
        // Block names + dtypes are the decode contract.
        let blocks = env["blocks"].as_object().unwrap();
        let mut bnames: Vec<&str> = blocks.keys().map(String::as_str).collect();
        bnames.sort_unstable();
        assert_eq!(
            bnames,
            [
                "cell_values",
                "indices",
                "positions",
                "tri_cell",
                "zone_ids"
            ]
        );
        assert_eq!(
            env["blocks"]["positions"]["dtype"],
            serde_json::json!("f32")
        );
        assert_eq!(env["blocks"]["indices"]["dtype"], serde_json::json!("u32"));
        assert_eq!(env["blocks"]["zone_ids"]["dtype"], serde_json::json!("u16"));
    }
}
