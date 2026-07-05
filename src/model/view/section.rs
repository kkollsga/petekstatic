//! [`IntersectionBundle`] — the vertical cross-section inspection bundle
//! ([`StaticModel::intersection_bundle`]).
//!
//! A section trace (an arbitrary world **polyline**, or a **bore trajectory**) is
//! marched through the grid's areal lattice; each column the trace crosses emits
//! an ordered [`SectionColumn`] carrying its distance-along-section, per-layer
//! property value, per-layer cell top/base depths (so the horizon traces are the
//! first top and last base), and — for a bore — the path's own z overlay. Raw
//! metres throughout; the viewer applies vertical exaggeration.

use super::frame::GridFrame;
use super::SCHEMA_VERSION;
use crate::error::StaticError;
use crate::grid::{Cell, Ijk};
use crate::model::model::StaticModel;
use crate::model::pipeline::areal_lattice;
use petektools::Lattice;
use serde::{Deserialize, Serialize};

/// A cell thinner than this (metres) is treated as an inactive/truncated layer in
/// a section column and rendered as `NaN` (see [`StaticModel::intersection_bundle`]).
const LAYER_ACTIVE_EPS_M: f64 = 1e-9;

/// Where a section is taken: a world-`(x, y)` **polyline**, or **along a bore**
/// trajectory (world `(x, y, z)` stations — the areal `(x, y)` drives the trace,
/// the `z` becomes each column's [`SectionColumn::path_z`] overlay).
#[derive(Debug, Clone, PartialEq)]
pub enum SectionSpec {
    /// A world-`[x, y]` polyline trace.
    Polyline(Vec<[f64; 2]>),
    /// A bore trajectory: world `[x, y, z]` stations (positive-down z).
    AlongBore { trajectory: Vec<[f64; 3]> },
}

/// One column the section trace crosses. `distance_m` is cumulative along the
/// trace at the column entry; `(x, y)` is the trace sample position (world). The
/// per-layer arrays are length `nk`, top→base. The **horizon traces** are the
/// first *active* `layer_tops` (structural top) and the last active `layer_bases`
/// (structural base).
///
/// ## Inactive layers = `NaN` (stable `nk`-sized schema)
/// Under a Follow conformity (`FollowTop`/`FollowBase`) the number of *active*
/// layers **varies per column** — thin columns truncate deep (resp. shallow)
/// layers against the pinch-out horizon. To keep the schema stable, every
/// `SectionColumn` still carries `nk`-sized arrays; a truncated (zero-thickness)
/// layer is written as `f64::NAN` in `layer_tops`, `layer_bases`, **and** `values`.
/// The viewer NaN-guards these (draws no cell there). Consumers must skip NaN
/// entries rather than assume every index is a live cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SectionColumn {
    pub distance_m: f64,
    pub i: usize,
    pub j: usize,
    pub x: f64,
    pub y: f64,
    /// Per-layer cell top depth \[m, positive-down\], top→base (length `nk`);
    /// `NaN` for an inactive/truncated layer. This is the column **centroid** trace
    /// (4-corner mean) — kept for hover + back-compat.
    pub layer_tops: Vec<f64>,
    /// Per-layer cell base depth \[m\], top→base (length `nk`); `NaN` for an
    /// inactive/truncated layer. Column centroid trace (see [`Self::layer_tops`]).
    pub layer_bases: Vec<f64>,
    /// Per-layer cell top depth at the column's **left** fence edge (bilinear from
    /// the cell's top ZCORN corners at the section's entry point through the cell),
    /// top→base (length `nk`); `NaN`-gapped like [`Self::layer_tops`]. With the four
    /// `*_l`/`*_r` arrays a cell renders as a **trapezoid following dip within the
    /// column** (not a flat "sugar cube"). Under sugar-cube mode these equal the
    /// centroid `layer_tops` (one viewer code path). (SCHEMA_VERSION 4, additive.)
    pub layer_tops_l: Vec<f64>,
    /// Per-layer cell top depth at the column's **right** fence edge. See
    /// [`Self::layer_tops_l`].
    pub layer_tops_r: Vec<f64>,
    /// Per-layer cell base depth at the column's **left** fence edge. See
    /// [`Self::layer_tops_l`].
    pub layer_bases_l: Vec<f64>,
    /// Per-layer cell base depth at the column's **right** fence edge. See
    /// [`Self::layer_tops_l`].
    pub layer_bases_r: Vec<f64>,
    /// Per-layer property value (length `nk`); empty when no property requested,
    /// `NaN` for an inactive/truncated layer.
    pub values: Vec<f64>,
    /// Per-layer zone id (length `nk`, top→base) — an index into
    /// [`IntersectionBundle::zones`] naming the stratigraphic zone each layer belongs
    /// to (the payload half of section colour-by-zone). Runs parallel to
    /// [`Self::values`] / [`Self::layer_tops`]; an inactive/truncated layer (where
    /// those are `NaN`) carries the gap sentinel [`SectionColumn::NO_ZONE`]
    /// (`u16::MAX`), which the viewer skips exactly as it skips a `NaN` depth.
    /// (SCHEMA_VERSION 5, additive — absent on older payloads.)
    #[serde(default)]
    pub zone_ids: Vec<u16>,
    /// The bore path's own depth at this station (`AlongBore` only; else `None`).
    pub path_z: Option<f64>,
}

impl SectionColumn {
    /// The [`Self::zone_ids`] gap sentinel: an inactive/truncated layer belongs to no
    /// zone (its geometry + value are `NaN`). Equals `u16::MAX`.
    pub const NO_ZONE: u16 = u16::MAX;
}

/// One stratigraphic zone referenced by a section's [`SectionColumn::zone_ids`]
/// (SCHEMA_VERSION 5): its name + an optional display colour (viewer hint). The
/// zone's id is its index in [`IntersectionBundle::zones`]. Sourced from the model's
/// zone table (`StackZone` name/colour on a horizon-stack build).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SectionZone {
    /// The zone's name (matches [`crate::model::Zone::name`]).
    pub name: String,
    /// Optional display colour (e.g. `"#ffcc00"`); `None` = viewer default.
    pub color: Option<String>,
}

/// A fluid contact along the section — a flat depth scalar (metres, positive-down)
/// the viewer draws as a horizontal line across the section.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SectionContact {
    pub kind: String,
    pub depth_m: f64,
}

/// One **interior** framework-horizon trace along the section (SCHEMA_VERSION 4,
/// `task_petekstatic_multizone_2`) — the depth polyline of a zone-bounding horizon
/// *between* the structural top and base, sampled at every column the section
/// crosses. The viewer draws it as a labelled polyline; `depths[c]` is the horizon
/// depth at [`IntersectionBundle::columns`]`[c]`, so the two arrays run parallel.
///
/// The top and base structural horizons are **not** repeated here — they are the
/// first active `layer_tops` / last active `layer_bases` of each column. Only the
/// `N − 2` interior horizons of an `N`-horizon stack are emitted, top→down; a
/// single-zone (2-surface) model has none, so `horizon_traces` is empty (the
/// additive block stays backward-compatible).
///
/// A horizon depth is the zone-top interface at that column (the top-depth of the
/// first cell of the zone the horizon bounds above). It is `NaN` only where the
/// column itself does not reach that interface (fully truncated), matching the
/// section's `NaN`-inactive convention.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HorizonTrace {
    /// The framework horizon's name.
    pub name: String,
    /// Depth \[m, positive-down\] at each section column, parallel to
    /// [`IntersectionBundle::columns`]; `NaN` where the column does not reach it.
    pub depths: Vec<f64>,
}

/// The vertical cross-section bundle: ordered columns plus the section-wide
/// contacts and the structural-surface labels.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IntersectionBundle {
    pub schema_version: u32,
    pub inputs_ref: String,
    /// **Sugar-cube mode** (SCHEMA_VERSION 4, additive): when `true`, cells are drawn
    /// as flat boxes and the per-column `*_l`/`*_r` edge arrays are flattened to the
    /// centroid trace (`layer_tops_l == layer_tops_r == layer_tops`, likewise bases).
    /// Default `false` — the engine geometry is corner-point (zone-following), so the
    /// edge arrays carry the true left/right cell depths and the viewer draws each
    /// cell as a dip-following trapezoid. Set via `with_sugar_cube` on the
    /// builder/template. Absent on older (pre-v4-edge) payloads → the viewer falls
    /// back to the flat path.
    pub sugar_cube: bool,
    /// Which property the per-column `values` carry (`None` = geometry only).
    pub property: Option<String>,
    /// The structural top / base horizon labels (the `layer_tops[0]` /
    /// `layer_bases[nk-1]` traces).
    pub top_name: String,
    pub base_name: String,
    /// Columns in trace order (by `distance_m`).
    pub columns: Vec<SectionColumn>,
    /// Per-interior-horizon depth traces (SCHEMA_VERSION 4), one per zone-bounding
    /// framework horizon between the structural top and base, top→down; each
    /// `depths` array runs parallel to `columns`. Empty for a single-zone model.
    pub horizon_traces: Vec<HorizonTrace>,
    /// The stratigraphic zones of the model, top→base (SCHEMA_VERSION 5) — the id
    /// table [`SectionColumn::zone_ids`] indexes into. Each carries `{name, color}`;
    /// the payload half of section colour-by-zone. A single-implicit-zone model still
    /// lists its one zone. (Additive — absent on older payloads.)
    #[serde(default)]
    pub zones: Vec<SectionZone>,
    /// Section-wide fluid contacts (GOC / OWC / GWC depths).
    pub contacts: Vec<SectionContact>,
}

impl IntersectionBundle {
    /// Stream this bundle to `w` as JSON (no intermediate `Value` tree) — the same
    /// streaming path the volume bundle's envelope uses (`view::wire`). Section
    /// bundles are small, so they stay plain JSON.
    pub fn write_json<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()> {
        super::wire::write_json(self, w)
    }
}

/// A section-trace vertex: world `(x, y)` and, for a bore, the path depth `z`.
struct Vert {
    x: f64,
    y: f64,
    z: Option<f64>,
}

impl SectionSpec {
    /// The trace vertices (xy, + optional path z).
    fn verts(&self) -> Vec<Vert> {
        match self {
            SectionSpec::Polyline(pts) => pts
                .iter()
                .map(|p| Vert {
                    x: p[0],
                    y: p[1],
                    z: None,
                })
                .collect(),
            SectionSpec::AlongBore { trajectory } => trajectory
                .iter()
                .map(|p| Vert {
                    x: p[0],
                    y: p[1],
                    z: Some(p[2]),
                })
                .collect(),
        }
    }
}

impl StaticModel {
    /// Export the vertical cross-section ([`IntersectionBundle`]) along `spec`,
    /// carrying `property`'s per-layer values (pass `None` for a geometry-only
    /// section). The trace is marched through the areal lattice at a sub-cell step;
    /// each distinct column crossed becomes an ordered [`SectionColumn`].
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the column lattice is smaller than `2x2`
    /// or not axis-aligned/regular, the trace has fewer than two vertices, or the
    /// named property is absent.
    pub fn intersection_bundle(
        &self,
        spec: &SectionSpec,
        property: Option<&str>,
    ) -> Result<IntersectionBundle, StaticError> {
        // Spilled models keep a 1×1×1 placeholder grid, so materialize the backing
        // for the (non-hot-path) section export — the same fix the shell export uses.
        let grid = self.view_grid()?;
        let frame = GridFrame::of_grid(&grid, self.georef())?;
        // The trace is authored in `frame` (world when georeferenced); the cell ZCORN
        // corners live in the grid's LOCAL lattice. `w2l` maps the trace point +
        // direction into that local frame before the fence clip, so
        // `fence_edge_depths` runs in ONE frame (identity when the frame already IS
        // the local lattice). See `WorldToLattice` / `task_petekstatic_section_edge_frame`.
        let w2l = WorldToLattice::new(&frame, &areal_lattice(&grid)?);
        let dims = grid.dims();
        let (ni, nj, nk) = (dims.ni, dims.nj, dims.nk);
        let sugar_cube = self.provenance().sugar_cube;

        let verts = spec.verts();
        if verts.len() < 2 {
            return Err(StaticError::InvalidInput(
                "intersection_bundle: a section trace needs at least two vertices".into(),
            ));
        }
        // Read the cube from the MATERIALIZED grid (not `self.property`, which is the
        // 1×1×1 placeholder on a spilled model — the same escape the geometry read
        // above avoids via `view_grid`). Otherwise every spilled section errored with
        // a spurious "no property" even though the cube lives in the backing.
        let prop = property
            .map(|name| {
                grid.properties().get(name).ok_or_else(|| {
                    StaticError::InvalidInput(format!("intersection_bundle: no property '{name}'"))
                })
            })
            .transpose()?;

        // Interior framework horizons (SCHEMA_VERSION 4): every zone except the
        // first contributes its top interface — a horizon strictly between the
        // structural top and base. Its depth at a column is the top-depth of the
        // first cell of the zone it bounds above (`k_start`). A single-zone model
        // has no interior horizon, so `interior` is empty and no trace is emitted.
        let interior: Vec<(String, usize)> = self
            .zones()
            .zones()
            .iter()
            .skip(1)
            .map(|z| (z.top_horizon.clone(), z.k_range.start))
            .collect();
        let mut trace_depths: Vec<Vec<f64>> = vec![Vec::new(); interior.len()];

        // Colour-by-zone payload (SCHEMA_VERSION 5): the section's zone table + the
        // per-layer zone id. `zone_of_k[k]` is the index (into `zones` below) of the
        // stratigraphic zone whose `k`-range contains layer `k`; a `k` outside every
        // zone range (should not happen for a partitioned column) stays `NO_ZONE`.
        let section_zones: Vec<SectionZone> = self
            .zones()
            .zones()
            .iter()
            .map(|z| SectionZone {
                name: z.name.clone(),
                color: z.color.clone(),
            })
            .collect();
        let mut zone_of_k = vec![SectionColumn::NO_ZONE; nk];
        for (zi, z) in self.zones().zones().iter().enumerate() {
            let id = u16::try_from(zi).unwrap_or(SectionColumn::NO_ZONE);
            let lo = z.k_range.start.min(nk);
            let hi = z.k_range.end.min(nk);
            for slot in &mut zone_of_k[lo..hi] {
                *slot = id;
            }
        }

        // Sub-cell march step: half the finer node spacing (honest cell sequencing
        // on the axis-aligned regular lattice this module is designed for).
        let step = 0.5 * frame.spacing_x.min(frame.spacing_y);
        let mut columns: Vec<SectionColumn> = Vec::new();
        let mut last_cell: Option<(usize, usize)> = None;
        let mut cum = 0.0f64;

        for w in verts.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            let seg = (b.x - a.x).hypot(b.y - a.y);
            let n_steps = (seg / step).ceil().max(1.0) as usize;
            for s in 0..=n_steps {
                let t = s as f64 / n_steps as f64;
                let (x, y) = (a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
                let d = cum + t * seg;
                let Some((fi, fj)) = frame_xy_to_ij(&frame, x, y) else {
                    continue;
                };
                let (ri, rj) = (fi.round(), fj.round());
                if ri < 0.0 || rj < 0.0 {
                    continue;
                }
                let (i, j) = (ri as usize, rj as usize);
                if i >= ni || j >= nj || last_cell == Some((i, j)) {
                    continue;
                }
                last_cell = Some((i, j));
                let path_z = match (a.z, b.z) {
                    (Some(za), Some(zb)) => Some(za + t * (zb - za)),
                    _ => None,
                };
                let mut layer_tops = Vec::with_capacity(nk);
                let mut layer_bases = Vec::with_capacity(nk);
                let mut layer_tops_l = Vec::with_capacity(nk);
                let mut layer_tops_r = Vec::with_capacity(nk);
                let mut layer_bases_l = Vec::with_capacity(nk);
                let mut layer_bases_r = Vec::with_capacity(nk);
                let mut values = Vec::new();
                let mut zone_ids = Vec::with_capacity(nk);
                // The section's local direction through this column (for the fence
                // entry/exit edge interpolation). For a `Polyline` the window vector
                // is the fence direction; for an `AlongBore` it can degenerate to
                // ~zero areal extent on a vertical / densely-sampled section, so the
                // edges below are RECOMPUTED from the trace's areal tangent in the
                // post-march pass (see `retangent_alongbore_edges`).
                let (segdx, segdy) = (b.x - a.x, b.y - a.y);
                // Map the (world/view-frame) trace point + fence direction into the
                // grid's LOCAL lattice frame — the frame the cell corners live in — so
                // the fence clip below never mixes frames (the section-edge defect).
                let (lpx, lpy) = w2l.point(x, y);
                let (ldx, ldy) = w2l.dir(segdx, segdy);
                for k in 0..nk {
                    let cell = grid.cell(Ijk::new(i, j, k));
                    // A layer truncated (or pinched) to zero thickness in this
                    // column is inactive: emit NaN for its geometry AND its
                    // property so the bundle stays nk-sized and stable while the
                    // viewer (which NaN-guards) skips it. `dz()` is the honest
                    // per-column cell thickness — zero exactly for a collapsed
                    // Follow-conformity truncation cell.
                    let (top_d, base_d, active) = if cell.dz() <= LAYER_ACTIVE_EPS_M {
                        (f64::NAN, f64::NAN, false)
                    } else {
                        (cell.top_depth(), cell.bottom_depth(), true)
                    };
                    layer_tops.push(top_d);
                    layer_bases.push(base_d);
                    // Left/right fence-edge depths (dip-following trapezoid). Under
                    // sugar-cube mode these collapse to the centroid trace (one
                    // viewer path). Inactive layers stay NaN-gapped.
                    let (tl, bl, tr, br) = if !active {
                        (f64::NAN, f64::NAN, f64::NAN, f64::NAN)
                    } else if sugar_cube {
                        (top_d, base_d, top_d, base_d)
                    } else {
                        fence_edge_depths(&cell, lpx, lpy, ldx, ldy)
                    };
                    layer_tops_l.push(tl);
                    layer_bases_l.push(bl);
                    layer_tops_r.push(tr);
                    layer_bases_r.push(br);
                    if let Some(p) = &prop {
                        values.push(if active {
                            p.values[(k * nj + j) * ni + i]
                        } else {
                            f64::NAN
                        });
                    }
                    // Zone id per layer, NaN-gapped like the geometry/value arrays: an
                    // active layer carries its zone index, an inactive one the sentinel.
                    zone_ids.push(if active {
                        zone_of_k.get(k).copied().unwrap_or(SectionColumn::NO_ZONE)
                    } else {
                        SectionColumn::NO_ZONE
                    });
                }
                // Interior-horizon depths at this column: the zone-top interface,
                // i.e. the top-depth of the zone's first cell. Even a fully
                // truncated (zero-thickness) zone-top cell carries the collapsed
                // interface location as its top-depth, so the trace stays a
                // continuous horizon polyline the viewer can draw.
                for (t, (_, k_start)) in trace_depths.iter_mut().zip(&interior) {
                    t.push(grid.cell(Ijk::new(i, j, *k_start)).top_depth());
                }
                columns.push(SectionColumn {
                    distance_m: d,
                    i,
                    j,
                    x,
                    y,
                    layer_tops,
                    layer_bases,
                    layer_tops_l,
                    layer_tops_r,
                    layer_bases_l,
                    layer_bases_r,
                    values,
                    zone_ids,
                    path_z,
                });
            }
            cum += seg;
        }

        // AlongBore fence direction (defect fix, `task_petekstatic_alongbore_edges`):
        // a bore column is emitted the first sample it enters a new cell, carrying
        // that MD-station window's raw `(segdx, segdy)`. On a vertical or
        // densely-sampled section that window has ~zero areal extent, so
        // `fence_edge_depths` degenerates to the centroid and the trapezoid draws
        // flat. Recompute the edges from the TRACE's areal tangent through each
        // column instead (Polyline edges are already the fence direction and stay
        // untouched — bit-for-bit).
        if matches!(spec, SectionSpec::AlongBore { .. }) && !sugar_cube {
            retangent_alongbore_edges(&grid, &w2l, nk, &mut columns);
        }

        let horizon_traces: Vec<HorizonTrace> = interior
            .into_iter()
            .zip(trace_depths)
            .map(|((name, _), depths)| HorizonTrace { name, depths })
            .collect();

        let contacts = self
            .contacts()
            .iter()
            .map(|c| SectionContact {
                kind: format!("{:?}", c.kind).to_uppercase(),
                depth_m: c.depth_m,
            })
            .collect();

        let name_for = |role: crate::wireframe::HorizonRole, fallback: &str| -> String {
            self.framework()
                .horizons
                .iter()
                .find(|h| h.role == role)
                .map_or_else(|| fallback.to_string(), |h| h.name.clone())
        };

        Ok(IntersectionBundle {
            schema_version: SCHEMA_VERSION,
            inputs_ref: self.provenance().inputs_ref.clone(),
            sugar_cube,
            property: property.map(String::from),
            top_name: name_for(crate::wireframe::HorizonRole::Top, "TOP"),
            base_name: name_for(crate::wireframe::HorizonRole::Base, "BASE"),
            columns,
            horizon_traces,
            zones: section_zones,
            contacts,
        })
    }
}

/// Areal-tangent floor (m): a trace tangent shorter than this is treated as having
/// no areal extent — a single-point (truly vertical) bore, the one honest degenerate
/// case that keeps the centroid trace.
const TANGENT_EPS_M: f64 = 1e-9;

/// Recompute an `AlongBore` section's left/right fence edges from the **trace's
/// areal tangent** through each column, replacing the raw MD-station micro-segment
/// direction (`task_petekstatic_alongbore_edges`).
///
/// The fence direction for column `c` is the central difference of its neighbouring
/// column centres (`columns[c+1] − columns[c-1]`, forward/backward at the ends) — the
/// areal tangent of the trace polyline through the column, robust to a vertical or
/// densely-sampled bore segment whose adjacent MD stations coincide areally. Where
/// the local tangent is itself degenerate (a trace that doubles back onto the column)
/// it falls back to the section's **overall azimuth** (first→last column centre).
///
/// **The one honest degenerate case:** a truly vertical bore whose trace is a single
/// areal point produces a single column with no tangent and no azimuth; there is no
/// meaningful fence direction, so the centroid trace (`l == r`) is kept — the cell is
/// drawn flat because the section has no areal extent to dip along. Any bore with
/// areal extent (even a near-vertical, slightly-drifting one) gets a real tangent and
/// dips correctly.
fn retangent_alongbore_edges(
    grid: &crate::grid::Grid,
    w2l: &WorldToLattice,
    nk: usize,
    columns: &mut [SectionColumn],
) {
    let n = columns.len();
    if n < 2 {
        return; // single (or no) column: no areal extent → keep centroid (honest).
    }
    // Overall azimuth (first→last column centre) — the fallback tangent. Column `x`/`y`
    // are in the view frame (world when georeferenced); the tangent + point are mapped
    // into the LOCAL lattice frame before the clip, exactly as the march path does.
    let (ax, ay) = (
        columns[n - 1].x - columns[0].x,
        columns[n - 1].y - columns[0].y,
    );
    for c in 0..n {
        let (px, py) = (columns[c].x, columns[c].y);
        let (i, j) = (columns[c].i, columns[c].j);
        let prev = &columns[c.saturating_sub(1)];
        let next = &columns[(c + 1).min(n - 1)];
        let (mut dx, mut dy) = (next.x - prev.x, next.y - prev.y);
        if dx.hypot(dy) < TANGENT_EPS_M {
            (dx, dy) = (ax, ay); // local tangent degenerate → overall azimuth
        }
        if dx.hypot(dy) < TANGENT_EPS_M {
            continue; // no areal extent anywhere → keep centroid (honest degenerate)
        }
        let (lpx, lpy) = w2l.point(px, py);
        let (ldx, ldy) = w2l.dir(dx, dy);
        for k in 0..nk {
            if !columns[c].layer_tops[k].is_finite() {
                continue; // inactive/truncated layer stays NaN-gapped
            }
            let cell = grid.cell(Ijk::new(i, j, k));
            let (tl, bl, tr, br) = fence_edge_depths(&cell, lpx, lpy, ldx, ldy);
            columns[c].layer_tops_l[k] = tl;
            columns[c].layer_bases_l[k] = bl;
            columns[c].layer_tops_r[k] = tr;
            columns[c].layer_bases_r[k] = br;
        }
    }
}

/// The cell's `(top_l, base_l, top_r, base_r)` depths at the section's **left** and
/// **right** fence edges — where the section line through `(px, py)` in direction
/// `(dx, dy)` enters and leaves the cell's areal rectangle. Each edge depth is the
/// **bilinear** interpolation of the cell's top (and bottom) ZCORN corners at that
/// entry/exit `(x, y)`. On a flat cell all four corners share a depth, so
/// `left == right == centroid`; on a cell dipping along the section they differ,
/// giving the trapezoid the viewer draws. The lattice is axis-aligned/regular (the
/// module limitation), so the areal quad is a rectangle.
///
/// ## Frame contract
/// `(px, py)` and `(dx, dy)` **must be in the cell's own (local grid lattice)
/// frame** — the same frame as `cell.corners`. Callers map a world/view-frame trace
/// through [`WorldToLattice`] first; feeding a world point against a local rectangle
/// is exactly the `task_petekstatic_section_edge_frame` defect (the clip misses and
/// the section collapses to a flat centroid trace). A `smin > smax` clip miss on a
/// non-degenerate direction therefore means a FRAME BUG, not a legitimate case, and
/// is flagged loudly (`debug_assert` + a `warn`) rather than silently swallowed.
///
/// Corners are ordered `di + 2·dj + 4·dk`: top face `0..4` = nodes (i,j), (i+1,j),
/// (i,j+1), (i+1,j+1); bottom face `4..8` in the same areal order.
fn fence_edge_depths(cell: &Cell, px: f64, py: f64, dx: f64, dy: f64) -> (f64, f64, f64, f64) {
    let centroid = || {
        let (t, b) = (cell.top_depth(), cell.bottom_depth());
        (t, b, t, b)
    };
    // The ONE legitimate `l == r` case: a trace with no areal extent through this
    // column (a single-point / truly vertical bore). No fence direction ⇒ no dip to
    // follow ⇒ keep the centroid. This is explicit and distinct from a clip miss.
    if dx.hypot(dy) < TANGENT_EPS_M {
        return centroid();
    }
    let c = &cell.corners;
    let (x0, x1) = (c[0].x, c[1].x);
    let (y0, y1) = (c[0].y, c[2].y);
    // Clip the parametric line P(s) = (px, py) + s·(dx, dy) to the cell rectangle.
    let (mut smin, mut smax) = (f64::NEG_INFINITY, f64::INFINITY);
    for &(p, d, lo, hi) in &[
        (px, dx, x0.min(x1), x0.max(x1)),
        (py, dy, y0.min(y1), y0.max(y1)),
    ] {
        if d.abs() < 1e-12 {
            continue; // parallel to this slab; P is inside, so it does not bound s
        }
        let (t0, t1) = ((lo - p) / d, (hi - p) / d);
        let (t0, t1) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        smin = smin.max(t0);
        smax = smax.min(t1);
    }
    // Bilinear interpolation of a 4-corner face at local (xx, yy).
    let interp = |xx: f64, yy: f64, face: &[crate::grid::Point3]| -> f64 {
        let frac = |q: f64, a: f64, b: f64| {
            let d = b - a;
            if d.abs() > 1e-12 {
                ((q - a) / d).clamp(0.0, 1.0)
            } else {
                0.5
            }
        };
        let (u, v) = (frac(xx, x0, x1), frac(yy, y0, y1));
        (1.0 - u) * (1.0 - v) * face[0].z
            + u * (1.0 - v) * face[1].z
            + (1.0 - u) * v * face[2].z
            + u * v * face[3].z
    };
    if !(smin.is_finite() && smax.is_finite()) || smin > smax {
        // A non-degenerate line that misses the cell it was matched to is a FRAME
        // MISMATCH (a point/direction not in the cell's local lattice frame), never a
        // legitimate section. Fail loudly in debug; warn + fall back to the centroid in
        // release so the degradation is visible, not silent.
        debug_assert!(
            false,
            "fence_edge_depths: line misses cell rectangle — frame bug (not the local \
             lattice frame?): p=({px}, {py}) dir=({dx}, {dy}) rect x[{x0}, {x1}] y[{y0}, {y1}]"
        );
        eprintln!(
            "warning: fence_edge_depths clip miss (frame mismatch?) — section edge fell \
             back to the centroid: p=({px}, {py}) dir=({dx}, {dy}) cell x[{x0}, {x1}] y[{y0}, {y1}]"
        );
        return centroid();
    }
    let (lx, ly) = (px + smin * dx, py + smin * dy);
    let (rx, ry) = (px + smax * dx, py + smax * dy);
    (
        interp(lx, ly, &c[0..4]),
        interp(lx, ly, &c[4..8]),
        interp(rx, ry, &c[0..4]),
        interp(rx, ry, &c[4..8]),
    )
}

/// The affine **world→lattice** areal map, so the fence clip runs in ONE frame.
///
/// The section trace is authored in the model's *view* frame — the world
/// [`GridFrame`] when a georef is registered (the standard real-data F5 config: a
/// LOCAL-origin cell lattice + a world georef). But a cell's ZCORN `corners` live
/// in the grid's own **local** lattice. [`fence_edge_depths`] clips the trace line
/// against a cell rectangle built from those corners, so the trace point +
/// direction MUST be mapped into the local lattice first — otherwise an azimuthed
/// world line misses every local rectangle (`smin > smax`) and silently degrades to
/// the centroid (the `task_petekstatic_section_edge_frame` defect).
///
/// Both frames are regular/axis-aligned, so the map is per-axis affine:
/// `local = lat_ori + (world − frame_ori) · (lat_inc / frame_inc)`. When the view
/// frame *is* the local lattice (no georef, or a grid whose corners are themselves
/// world-valued) origins + spacings coincide and this reduces to the identity — the
/// local fixtures stay bit-for-bit.
struct WorldToLattice {
    frame_ox: f64,
    frame_oy: f64,
    lat_ox: f64,
    lat_oy: f64,
    scale_x: f64,
    scale_y: f64,
}

impl WorldToLattice {
    /// Build the map from the view [`GridFrame`] to the grid's local [`Lattice`].
    fn new(frame: &GridFrame, lat: &Lattice) -> Self {
        let scale = |lat_inc: f64, frame_inc: f64| {
            if frame_inc != 0.0 {
                lat_inc / frame_inc
            } else {
                1.0
            }
        };
        Self {
            frame_ox: frame.origin_x,
            frame_oy: frame.origin_y,
            lat_ox: lat.xori,
            lat_oy: lat.yori,
            scale_x: scale(lat.xinc, frame.spacing_x),
            scale_y: scale(lat.yinc, frame.spacing_y),
        }
    }

    /// A view-frame `(x, y)` in the local lattice frame.
    fn point(&self, x: f64, y: f64) -> (f64, f64) {
        (
            self.lat_ox + (x - self.frame_ox) * self.scale_x,
            self.lat_oy + (y - self.frame_oy) * self.scale_y,
        )
    }

    /// A view-frame direction `(dx, dy)` in the local lattice frame (origin-free).
    fn dir(&self, dx: f64, dy: f64) -> (f64, f64) {
        (dx * self.scale_x, dy * self.scale_y)
    }
}

/// World `(x, y)` → fractional `(i, j)` on the areal frame (regular axis-aligned).
fn frame_xy_to_ij(frame: &GridFrame, x: f64, y: f64) -> Option<(f64, f64)> {
    if frame.spacing_x == 0.0 || frame.spacing_y == 0.0 {
        return None;
    }
    Some((
        (x - frame.origin_x) / frame.spacing_x,
        (y - frame.origin_y) / frame.spacing_y,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gridder::{Conformity, SolveOpts};
    use crate::model::{BuildOpts, ConstantPriors, StaticModelBuilder};
    use crate::wireframe::{
        Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
    };

    fn opts() -> BuildOpts {
        BuildOpts {
            area_m2: 100.0,
            gross_height_m: 50.0,
            nk: 5,
            conformity: Conformity::Proportional,
            solve_opts: SolveOpts::default(),
            priors: ConstantPriors {
                porosity: 0.25,
                net_to_gross: 0.8,
                water_saturation: 0.3,
            },
        }
    }

    fn wedge_model(n: usize, thin: f64, thick: f64) -> StaticModel {
        StaticModelBuilder::from_wireframe(&wedge_wf(n, thin, thick), opts())
            .unwrap()
            .build()
            .unwrap()
    }

    /// A wedge: flat Top at 5000, a Base dipping in i from 5000+thin (col 0) to
    /// 5000+thick (col n-1); constant in j. Contact deep (whole column counts).
    fn wedge_wf(n: usize, thin: f64, thick: f64) -> Wireframe {
        let top = vec![5000.0; n * n];
        let mut base = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                base[r * n + c] = 5000.0 + thin + (thick - thin) * (c as f64 / (n - 1) as f64);
            }
        }
        Wireframe {
            boundary: Boundary {
                ring: vec![
                    [0.0, 0.0],
                    [10.0, 0.0],
                    [10.0, 10.0],
                    [0.0, 10.0],
                    [0.0, 0.0],
                ],
                hardness: Hardness::Interpolated,
            },
            horizons: std::sync::Arc::new(vec![
                Horizon {
                    name: "TopRes".into(),
                    role: HorizonRole::Top,
                    surface: GriddedDepth {
                        ncol: n,
                        nrow: n,
                        depth_m: top,
                        is_control: vec![true; n * n],
                    },
                },
                Horizon {
                    name: "BaseRes".into(),
                    role: HorizonRole::Base,
                    surface: GriddedDepth {
                        ncol: n,
                        nrow: n,
                        depth_m: base,
                        is_control: vec![true; n * n],
                    },
                },
            ]),
            contacts: vec![Contact {
                kind: ContactKind::Owc,
                depth_m: 6000.0,
                hardness: Hardness::Hard,
            }],
        }
    }

    #[test]
    fn section_edge_arrays_follow_dip_bitwise_and_sugar_cube_flattens() {
        // The section runs along +x at y = 0.5 (each cell's j-centroid), so the fence
        // left/right edges are the cell's node-x faces (x0, x1) at v = 0.5 — i.e. the
        // means of the two left / two right ZCORN corners. On the dipping wedge these
        // differ (trapezoid); on a flat model and under sugar-cube mode they collapse
        // to the centroid trace.
        let spec = SectionSpec::Polyline(vec![[0.5, 0.5], [9.5, 0.5]]);

        // (a) Dipping wedge: edges follow dip and match direct ZCORN interpolation
        //     bit-for-bit.
        let m = wedge_model(11, 20.0, 120.0);
        let b = m.intersection_bundle(&spec, None).unwrap();
        assert!(!b.sugar_cube, "trapezoid mode is the default");
        assert_eq!(b.schema_version, SCHEMA_VERSION);
        let col = &b.columns[5];
        let (i, j, k) = (col.i, col.j, 2usize);
        let c = &m.grid().cell(Ijk::new(i, j, k)).corners;
        // At u = 0 (left face) / u = 1 (right face), v = 0.5: the bilinear reduces to
        // the mean of the two corners on that face — computed here EXACTLY as the
        // encoder's interp does, so the compare is bit-level.
        let (top_l, top_r) = (0.5 * c[0].z + 0.5 * c[2].z, 0.5 * c[1].z + 0.5 * c[3].z);
        let (base_l, base_r) = (0.5 * c[4].z + 0.5 * c[6].z, 0.5 * c[5].z + 0.5 * c[7].z);
        assert_eq!(col.layer_tops_l[k], top_l, "left top edge bit-level");
        assert_eq!(col.layer_tops_r[k], top_r, "right top edge bit-level");
        assert_eq!(col.layer_bases_l[k], base_l, "left base edge bit-level");
        assert_eq!(col.layer_bases_r[k], base_r, "right base edge bit-level");
        // The wedge dips along the section, so left != right (a real trapezoid).
        assert!(
            (col.layer_bases_l[k] - col.layer_bases_r[k]).abs() > 1e-3,
            "base edge dips across the column"
        );
        assert!(
            (col.layer_tops_l[k] - col.layer_tops_r[k]).abs() > 1e-3,
            "interior top edge dips across the column"
        );
        // Centroid trace stays the 4-corner mean, between the two edges.
        assert!((col.layer_tops[k] - 0.5 * (top_l + top_r)).abs() < 1e-9);

        // (b) Flat model: left == right == centroid.
        let flat = wedge_model(11, 50.0, 50.0);
        let fb = flat.intersection_bundle(&spec, None).unwrap();
        for col in &fb.columns {
            for k in 0..col.layer_tops.len() {
                if !col.layer_tops[k].is_finite() {
                    continue;
                }
                assert_eq!(col.layer_tops_l[k], col.layer_tops[k]);
                assert_eq!(col.layer_tops_r[k], col.layer_tops[k]);
                assert_eq!(col.layer_bases_l[k], col.layer_bases[k]);
                assert_eq!(col.layer_bases_r[k], col.layer_bases[k]);
            }
        }

        // (c) Sugar-cube mode on the SAME dipping wedge: edges flatten to the centroid.
        let sc = StaticModelBuilder::from_wireframe(&wedge_wf(11, 20.0, 120.0), opts())
            .unwrap()
            .with_sugar_cube(true)
            .build()
            .unwrap();
        let sb = sc.intersection_bundle(&spec, None).unwrap();
        assert!(sb.sugar_cube, "flag flows into the bundle");
        for col in &sb.columns {
            for k in 0..col.layer_tops.len() {
                if !col.layer_tops[k].is_finite() {
                    continue;
                }
                assert_eq!(col.layer_tops_l[k], col.layer_tops[k], "sugar: l==centroid");
                assert_eq!(col.layer_tops_r[k], col.layer_tops[k], "sugar: r==centroid");
                assert_eq!(col.layer_bases_l[k], col.layer_bases[k]);
                assert_eq!(col.layer_bases_r[k], col.layer_bases[k]);
            }
        }
    }

    #[test]
    fn section_columns_are_ordered_and_track_the_wedge() {
        // 11-node horizon -> ni=nj=10, area 100 -> side 10, dx=dy=1, col-0 centroid
        // x=0.5. A trace along i at y=0.5 crosses columns j=0, i=0..9 in order.
        let m = wedge_model(11, 20.0, 120.0);
        let spec = SectionSpec::Polyline(vec![[0.5, 0.5], [9.5, 0.5]]);
        let b = m.intersection_bundle(&spec, Some("PORO")).unwrap();
        assert_eq!(b.property.as_deref(), Some("PORO"));
        assert_eq!(b.top_name, "TopRes");
        assert_eq!(b.base_name, "BaseRes");
        // One column per i, ordered by distance and by i.
        assert_eq!(b.columns.len(), 10);
        for (n, col) in b.columns.iter().enumerate() {
            assert_eq!(col.i, n);
            assert_eq!(col.j, 0);
            if n > 0 {
                assert!(
                    col.distance_m > b.columns[n - 1].distance_m,
                    "distance monotone"
                );
            }
            // Per-layer arrays are nk long; property carried.
            assert_eq!(col.layer_tops.len(), 5);
            assert_eq!(col.values.len(), 5);
            assert!(col.path_z.is_none());
        }
        // Horizon trace: structural top is flat 5000; base thickens downdip.
        let first = &b.columns[0];
        let last = &b.columns[9];
        assert!(
            (first.layer_tops[0] - 5000.0).abs() < 1e-6,
            "flat top updip"
        );
        assert!(
            (last.layer_tops[0] - 5000.0).abs() < 1e-6,
            "flat top downdip"
        );
        let thick_updip = first.layer_bases[4] - first.layer_tops[0];
        let thick_downdip = last.layer_bases[4] - last.layer_tops[0];
        // Analytic wedge column thickness ~ 20 (updip) -> ~120 (downdip).
        assert!(
            (thick_updip - 24.0).abs() < 6.0,
            "updip thickness {thick_updip} ~ 20"
        );
        assert!(
            (thick_downdip - 116.0).abs() < 8.0,
            "downdip thickness {thick_downdip} ~ 120"
        );
        assert!(thick_downdip > thick_updip * 3.0, "wedge thickens downdip");
        // The contact is carried section-wide.
        assert_eq!(b.contacts.len(), 1);
        assert_eq!(b.contacts[0].kind, "OWC");
    }

    #[test]
    fn along_bore_overlays_the_path_z() {
        let m = wedge_model(11, 20.0, 120.0);
        // A bore descending from 5000 to 5100 as it steps in i.
        let traj = vec![[0.5, 0.5, 5000.0], [9.5, 0.5, 5100.0]];
        let b = m
            .intersection_bundle(&SectionSpec::AlongBore { trajectory: traj }, Some("PORO"))
            .unwrap();
        assert!(
            b.columns.iter().all(|c| c.path_z.is_some()),
            "bore z overlaid"
        );
        let z0 = b.columns.first().unwrap().path_z.unwrap();
        let z1 = b.columns.last().unwrap().path_z.unwrap();
        assert!(z0 < z1, "path descends: {z0} -> {z1}");
        assert!((z0 - 5000.0).abs() < 15.0 && (z1 - 5100.0).abs() < 15.0);
    }

    #[test]
    fn along_bore_edges_follow_trace_tangent_and_vertical_convention() {
        // `task_petekstatic_alongbore_edges`: an AlongBore section's fence direction
        // must come from the TRACE's areal tangent, not the raw MD-station
        // micro-segment (which degenerates to the centroid on a vertical / densely
        // sampled window and draws the trapezoid flat).
        let m = wedge_model(11, 20.0, 120.0);
        let k = 2usize; // an interior (dipping) layer of the wedge

        // (a) VERTICAL KICKOFF then deviate: the surface column is entered via a
        //     zero-areal-extent first window. Pre-fix it collapsed to the centroid
        //     (l == r, flat); the trace tangent now gives it real dip.
        let mut kick = vec![[0.5, 0.5, 5000.0], [0.5, 0.5, 5010.0]];
        for s in 1..=9 {
            kick.push([0.5 + s as f64, 0.5, 5010.0 + s as f64 * 10.0]);
        }
        let kb = m
            .intersection_bundle(&SectionSpec::AlongBore { trajectory: kick }, None)
            .unwrap();
        let c0 = &kb.columns[0];
        assert_eq!(c0.i, 0, "first column is the surface (kickoff) cell");
        assert!(
            (c0.layer_bases_l[k] - c0.layer_bases_r[k]).abs() > 1e-3,
            "kickoff column dips across the trace tangent (was flat pre-fix): \
             l={} r={}",
            c0.layer_bases_l[k],
            c0.layer_bases_r[k]
        );

        // (b) A straight deviated bore along +x recovers the straight-fence
        //     direction: its edge arrays match a Polyline along the same line
        //     BIT-FOR-BIT, and match direct ZCORN corner interpolation at v = 0.5.
        let ab = m
            .intersection_bundle(
                &SectionSpec::AlongBore {
                    trajectory: vec![[0.5, 0.5, 5000.0], [9.5, 0.5, 5100.0]],
                },
                None,
            )
            .unwrap();
        let pl = m
            .intersection_bundle(&SectionSpec::Polyline(vec![[0.5, 0.5], [9.5, 0.5]]), None)
            .unwrap();
        assert_eq!(ab.columns.len(), pl.columns.len());
        let mut dipping = 0;
        for (a, p) in ab.columns.iter().zip(&pl.columns) {
            assert_eq!((a.i, a.j), (p.i, p.j));
            for kk in 0..a.layer_tops_l.len() {
                if !a.layer_tops_l[kk].is_finite() {
                    continue;
                }
                assert_eq!(
                    a.layer_tops_l[kk], p.layer_tops_l[kk],
                    "tangent == fence l top"
                );
                assert_eq!(
                    a.layer_tops_r[kk], p.layer_tops_r[kk],
                    "tangent == fence r top"
                );
                assert_eq!(a.layer_bases_l[kk], p.layer_bases_l[kk]);
                assert_eq!(a.layer_bases_r[kk], p.layer_bases_r[kk]);
                if (a.layer_bases_l[kk] - a.layer_bases_r[kk]).abs() > 1e-3 {
                    dipping += 1;
                }
            }
        }
        assert!(dipping > 0, "deviated bore dips on dipping cells (l != r)");
        // Direct ZCORN corner interpolation along +x at v = 0.5 (bit-level): the base
        // edges are the means of the two left / two right bottom-face corners.
        let col = &ab.columns[5];
        let cc = &m.grid().cell(Ijk::new(col.i, col.j, k)).corners;
        let (bl, br) = (0.5 * cc[4].z + 0.5 * cc[6].z, 0.5 * cc[5].z + 0.5 * cc[7].z);
        assert_eq!(
            col.layer_bases_l[k], bl,
            "bore base_l == direct ZCORN interp"
        );
        assert_eq!(
            col.layer_bases_r[k], br,
            "bore base_r == direct ZCORN interp"
        );

        // (c) The one honest degenerate case: a perfectly VERTICAL bore (single areal
        //     point) has no trace tangent, so the edges stay the centroid (l == r).
        let vb = m
            .intersection_bundle(
                &SectionSpec::AlongBore {
                    trajectory: vec![[3.5, 3.5, 5000.0], [3.5, 3.5, 5100.0]],
                },
                None,
            )
            .unwrap();
        assert_eq!(vb.columns.len(), 1, "single areal point -> one column");
        let v = &vb.columns[0];
        for kk in 0..v.layer_tops.len() {
            if !v.layer_tops[kk].is_finite() {
                continue;
            }
            assert_eq!(
                v.layer_tops_l[kk], v.layer_tops[kk],
                "vertical: l == centroid"
            );
            assert_eq!(
                v.layer_tops_r[kk], v.layer_tops[kk],
                "vertical: r == centroid"
            );
            assert_eq!(v.layer_bases_l[kk], v.layer_bases[kk]);
            assert_eq!(v.layer_bases_r[kk], v.layer_bases[kk]);
        }
    }

    #[test]
    fn geometry_only_section_and_error_paths() {
        let m = wedge_model(11, 20.0, 120.0);
        let spec = SectionSpec::Polyline(vec![[0.5, 0.5], [9.5, 0.5]]);
        // No property -> empty per-column values.
        let b = m.intersection_bundle(&spec, None).unwrap();
        assert!(b.property.is_none());
        assert!(b.columns.iter().all(|c| c.values.is_empty()));
        // Degenerate trace + missing property error.
        assert!(m
            .intersection_bundle(&SectionSpec::Polyline(vec![[0.5, 0.5]]), None)
            .is_err());
        assert!(m.intersection_bundle(&spec, Some("NOPE")).is_err());
    }

    // A UTM-origin (world-frame) wedge, built through `from_wireframe` and given
    // its registered world georeference. Grid stays a local area-scaled square
    // (side 300, dx=30); the georef labels column (0,0)'s world centroid.
    const S_UTM_X0: f64 = 552_000.0;
    const S_UTM_Y0: f64 = 6_805_000.0;
    const S_UTM_INC: f64 = 30.0;

    fn utm_wedge() -> StaticModel {
        let n = 11usize;
        let top = vec![5000.0; n * n];
        let mut base = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                base[r * n + c] = 5000.0 + 20.0 + 100.0 * (c as f64 / (n - 1) as f64);
            }
        }
        let wf = Wireframe {
            boundary: Boundary {
                ring: vec![
                    [S_UTM_X0, S_UTM_Y0],
                    [S_UTM_X0 + 300.0, S_UTM_Y0],
                    [S_UTM_X0 + 300.0, S_UTM_Y0 + 300.0],
                    [S_UTM_X0, S_UTM_Y0 + 300.0],
                    [S_UTM_X0, S_UTM_Y0],
                ],
                hardness: Hardness::Hard,
            },
            horizons: std::sync::Arc::new(vec![
                Horizon {
                    name: "TopRes".into(),
                    role: HorizonRole::Top,
                    surface: GriddedDepth {
                        ncol: n,
                        nrow: n,
                        depth_m: top,
                        is_control: vec![true; n * n],
                    },
                },
                Horizon {
                    name: "BaseRes".into(),
                    role: HorizonRole::Base,
                    surface: GriddedDepth {
                        ncol: n,
                        nrow: n,
                        depth_m: base,
                        is_control: vec![true; n * n],
                    },
                },
            ]),
            contacts: vec![Contact {
                kind: ContactKind::Owc,
                depth_m: 6000.0,
                hardness: Hardness::Hard,
            }],
        };
        let mut o = opts();
        o.area_m2 = 90_000.0; // side 300 -> dx = 30 over 10 columns
        StaticModelBuilder::from_wireframe(&wf, o)
            .unwrap()
            .with_georef(
                S_UTM_X0 + S_UTM_INC / 2.0,
                S_UTM_Y0 + S_UTM_INC / 2.0,
                S_UTM_INC,
                S_UTM_INC,
            )
            .build()
            .unwrap()
    }

    #[test]
    fn utm_world_fence_yields_ordered_columns() {
        // F5-class fix: a fence polyline in WORLD (UTM) coordinates now traces the
        // world-georeferenced frame and yields non-empty, ordered columns. Under
        // the old local frame the UTM trace marched outside a 0..300 lattice and
        // produced ZERO columns.
        let m = utm_wedge();
        let f = m.map_bundle(&crate::model::MapSpec::new()).unwrap().frame;
        let ymid = f.origin_y + (f.nrow as f64 - 1.0) * f.spacing_y / 2.0;
        let x0 = f.origin_x;
        let x1 = f.origin_x + (f.ncol as f64 - 1.0) * f.spacing_x;
        let line = SectionSpec::Polyline(vec![[x0, ymid], [x1, ymid]]);
        let b = m.intersection_bundle(&line, Some("PORO")).unwrap();
        assert!(!b.columns.is_empty(), "world fence produced no columns");
        assert_eq!(b.columns.len(), 10, "one column per crossed i");
        for (n, col) in b.columns.iter().enumerate() {
            assert_eq!(col.i, n);
            // Sample position is carried in WORLD coordinates.
            assert!(col.x > 500_000.0, "column x is world: {}", col.x);
            if n > 0 {
                assert!(col.distance_m > b.columns[n - 1].distance_m);
            }
            assert_eq!(col.layer_tops.len(), 5);
            assert_eq!(col.values.len(), 5);
        }
        // Wedge still tracks: base thickens downdip in i.
        let thin = b.columns[0].layer_bases[4] - b.columns[0].layer_tops[0];
        let thick = b.columns[9].layer_bases[4] - b.columns[9].layer_tops[0];
        assert!(
            thick > thin * 3.0,
            "wedge thickens downdip: {thin} -> {thick}"
        );
    }

    // A LOCAL-origin lattice (boundary at origin 0, so the cell CORNERS are local
    // ~0..10) carrying a registered WORLD georeference — the standard real-data F5
    // configuration (petekIO hands a local lattice + a world georef). This is
    // *distinct* from `utm_wedge`, whose boundary ring is itself world-valued (a
    // single consistent world frame, which is why its corners == its georef). Here
    // the corners are LOCAL while the georef labels the WORLD lattice, so a section
    // trace given in WORLD coordinates (exactly as peteksim drives it) exercises the
    // world→lattice frame seam that `fence_edge_depths` sits on.
    const W_X0: f64 = 500_000.0;
    const W_Y0: f64 = 6_800_000.0;

    fn georef_local_wedge(thin: f64, thick: f64) -> StaticModel {
        // wedge_wf: boundary [0,0]..[10,10], area 100 -> ni=nj=10, local dx=dy=1;
        // column (0,0)'s local centroid is (0.5, 0.5). The georef labels that node's
        // WORLD position (W_X0+0.5, W_Y0+0.5) with a matching 1 m increment, so a
        // world trace maps onto the same cells the local fixtures use.
        StaticModelBuilder::from_wireframe(&wedge_wf(11, thin, thick), opts())
            .unwrap()
            .with_georef(W_X0 + 0.5, W_Y0 + 0.5, 1.0, 1.0)
            .build()
            .unwrap()
    }

    #[test]
    fn georef_world_fence_and_bore_edges_are_frame_invariant_and_dip() {
        // REGRESSION (`task_petekstatic_section_edge_frame`, the third world/local
        // silent seam bug after the bundle-frame F5 and the collocated no-op): on a
        // world-georeferenced model the section trace arrives in WORLD coordinates
        // while the cell ZCORN corners live in the LOCAL lattice. `fence_edge_depths`
        // clipped the WORLD point + direction against the LOCAL cell rectangle; for
        // any azimuthed (non-axis-aligned) trace the x- and y-slab crossings fall at
        // wildly different s, so `smin > smax` — the clip MISSES every cell and the
        // fallback silently degraded to entry==exit==centroid -> layer_tops_l ==
        // layer_tops_r EXACTLY (a flat "sugar cube" everywhere).
        //
        // Two checks, both on a DIAGONAL fence (a straight axis-aligned fence hides the
        // bug — its single bounding slab's clip is offset-invariant, which is why the
        // earlier axis-aligned `utm_*` fixtures passed while the real azimuthed model
        // did not):
        //   1. The world-georef section DIPS (l != r) instead of collapsing to the
        //      centroid — the direct reproduction: pre-fix every pair was l == r EXACTLY.
        //   2. FRAME CONSISTENCY: the SAME wedge geometry expressed in the local frame
        //      (no georef) vs the local lattice + a world georef produces the SAME edge
        //      geometry. Not bit-identical — the world path does its trace arithmetic at
        //      ~5e5 magnitude and loses ~1e-10 of mantissa vs the local path — so the
        //      match is to a tight absolute tolerance, well below any real dip.
        let local = wedge_model(11, 20.0, 120.0); // boundary [0,0]..[10,10], no georef
        let world = georef_local_wedge(20.0, 120.0); // SAME corners + world georef
        let k = 2usize; // an interior (dipping) layer of the wedge
        const FRAME_TOL: f64 = 1e-6;

        // Diagonal trace: local (0.5,0.5)->(9.5,9.5), and its WORLD image under the
        // georef (origin W_X0+0.5 / W_Y0+0.5, 1 m increment) is a pure +offset.
        let l0 = [0.5, 0.5];
        let l1 = [9.5, 9.5];
        let w0 = [W_X0 + l0[0], W_Y0 + l0[1]];
        let w1 = [W_X0 + l1[0], W_Y0 + l1[1]];

        let assert_dips_and_consistent = |lb: &IntersectionBundle, wb: &IntersectionBundle| {
            assert!(!wb.columns.is_empty(), "world trace produced no columns");
            assert_eq!(
                lb.columns.len(),
                wb.columns.len(),
                "same geometry -> same crossed columns"
            );
            let mut dipping = 0;
            for (lc, wc) in lb.columns.iter().zip(&wb.columns) {
                assert_eq!((lc.i, lc.j), (wc.i, wc.j), "same (i,j) sequence");
                for kk in 0..lc.layer_tops_l.len() {
                    if !lc.layer_tops_l[kk].is_finite() {
                        continue;
                    }
                    // Frame consistency: world-georef edges == local edges (tight tol).
                    assert!((wc.layer_tops_l[kk] - lc.layer_tops_l[kk]).abs() < FRAME_TOL);
                    assert!((wc.layer_tops_r[kk] - lc.layer_tops_r[kk]).abs() < FRAME_TOL);
                    assert!((wc.layer_bases_l[kk] - lc.layer_bases_l[kk]).abs() < FRAME_TOL);
                    assert!((wc.layer_bases_r[kk] - lc.layer_bases_r[kk]).abs() < FRAME_TOL);
                    // The world-georef section must DIP (would be exactly 0 pre-fix).
                    if (wc.layer_bases_l[kk] - wc.layer_bases_r[kk]).abs() > 1e-3 {
                        dipping += 1;
                    }
                }
            }
            assert!(
                dipping > 0,
                "world-georef section must dip (l != r); pre-fix it collapsed to the centroid"
            );
        };

        // (a) straight (diagonal) Polyline.
        let lb = local
            .intersection_bundle(&SectionSpec::Polyline(vec![l0, l1]), None)
            .unwrap();
        let wb = world
            .intersection_bundle(&SectionSpec::Polyline(vec![w0, w1]), None)
            .unwrap();
        assert!(
            wb.columns[0].x > 400_000.0,
            "world sample x: {}",
            wb.columns[0].x
        );
        assert_dips_and_consistent(&lb, &wb);

        // (b) AlongBore with the same diagonal areal path (world trajectory).
        let lab = local
            .intersection_bundle(
                &SectionSpec::AlongBore {
                    trajectory: vec![[l0[0], l0[1], 5000.0], [l1[0], l1[1], 5100.0]],
                },
                None,
            )
            .unwrap();
        let wab = world
            .intersection_bundle(
                &SectionSpec::AlongBore {
                    trajectory: vec![[w0[0], w0[1], 5000.0], [w1[0], w1[1], 5100.0]],
                },
                None,
            )
            .unwrap();
        assert!(
            wab.columns.iter().all(|c| c.path_z.is_some()),
            "bore z overlaid"
        );
        assert_dips_and_consistent(&lab, &wab);

        // BIT-LEVEL vs direct ZCORN corner interpolation: a straight world +x fence at a
        // j-centroid maps (through the world→lattice frame) to the local cell's node-x
        // faces at v = 0.5, so each edge is EXACTLY the mean of the two left / two right
        // ZCORN corners — computed here as the encoder's `interp` does. This anchors the
        // frame mapping to the true corner geometry (a wrong scale/offset would shift u/v
        // off these exact values).
        let (bx0, bx1, bymid) = (W_X0 + 0.5, W_X0 + 9.5, W_Y0 + 0.5);
        let sb = world
            .intersection_bundle(
                &SectionSpec::Polyline(vec![[bx0, bymid], [bx1, bymid]]),
                None,
            )
            .unwrap();
        let col = &sb.columns[5];
        let cc = &world.grid().cell(Ijk::new(col.i, col.j, k)).corners;
        let (top_l, top_r) = (0.5 * cc[0].z + 0.5 * cc[2].z, 0.5 * cc[1].z + 0.5 * cc[3].z);
        let (base_l, base_r) = (0.5 * cc[4].z + 0.5 * cc[6].z, 0.5 * cc[5].z + 0.5 * cc[7].z);
        assert_eq!(
            col.layer_tops_l[k], top_l,
            "world fence top_l == direct ZCORN interp"
        );
        assert_eq!(
            col.layer_tops_r[k], top_r,
            "world fence top_r == direct ZCORN interp"
        );
        assert_eq!(
            col.layer_bases_l[k], base_l,
            "world fence base_l == direct ZCORN interp"
        );
        assert_eq!(
            col.layer_bases_r[k], base_r,
            "world fence base_r == direct ZCORN interp"
        );
    }

    #[test]
    fn utm_along_bore_trajectory_sections() {
        // A bore with a UTM trajectory sections the world frame and overlays its z.
        let m = utm_wedge();
        let f = m.map_bundle(&crate::model::MapSpec::new()).unwrap().frame;
        let ymid = f.origin_y + (f.nrow as f64 - 1.0) * f.spacing_y / 2.0;
        let x0 = f.origin_x;
        let x1 = f.origin_x + (f.ncol as f64 - 1.0) * f.spacing_x;
        let traj = vec![[x0, ymid, 5000.0], [x1, ymid, 5100.0]];
        let b = m
            .intersection_bundle(&SectionSpec::AlongBore { trajectory: traj }, Some("PORO"))
            .unwrap();
        assert!(!b.columns.is_empty(), "UTM bore produced no columns");
        assert!(
            b.columns.iter().all(|c| c.path_z.is_some()),
            "bore z overlaid"
        );
        let z0 = b.columns.first().unwrap().path_z.unwrap();
        let z1 = b.columns.last().unwrap().path_z.unwrap();
        assert!(z0 < z1, "path descends: {z0} -> {z1}");
    }

    #[test]
    fn section_json_round_trips_and_keys_are_stable() {
        let m = wedge_model(11, 20.0, 120.0);
        let spec = SectionSpec::Polyline(vec![[0.5, 0.5], [9.5, 0.5]]);
        let b = m.intersection_bundle(&spec, Some("PORO")).unwrap();
        let json = serde_json::to_string(&b).unwrap();
        let back: IntersectionBundle = serde_json::from_str(&json).unwrap();
        assert_eq!(b, back);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let mut keys: Vec<&str> = v.as_object().unwrap().keys().map(String::as_str).collect();
        keys.sort_unstable();
        assert_eq!(
            keys,
            [
                "base_name",
                "columns",
                "contacts",
                "horizon_traces",
                "inputs_ref",
                "property",
                "schema_version",
                "sugar_cube",
                "top_name",
                "zones",
            ]
        );
        let col = v["columns"][0].as_object().unwrap();
        let mut ck: Vec<&str> = col.keys().map(String::as_str).collect();
        ck.sort_unstable();
        assert_eq!(
            ck,
            [
                "distance_m",
                "i",
                "j",
                "layer_bases",
                "layer_bases_l",
                "layer_bases_r",
                "layer_tops",
                "layer_tops_l",
                "layer_tops_r",
                "path_z",
                "values",
                "x",
                "y",
                "zone_ids",
            ]
        );
    }

    /// A 3-horizon / 2-zone flat stack model, UPPER coloured, LOWER uncoloured — the
    /// fixture for the section colour-by-zone rider.
    fn two_zone_stack_model() -> StaticModel {
        use crate::model::{HorizonSource, HorizonStack, StackHorizon, StackZone};
        const N: usize = 11;
        let flat = |d: f64| GriddedDepth {
            ncol: N,
            nrow: N,
            depth_m: vec![d; N * N],
            is_control: vec![true; N * N],
        };
        let stack = HorizonStack {
            horizons: vec![
                StackHorizon {
                    name: "TOP".into(),
                    source: HorizonSource::Mapped(flat(5000.0)),
                },
                StackHorizon {
                    name: "MID".into(),
                    source: HorizonSource::Mapped(flat(5030.0)),
                },
                StackHorizon {
                    name: "BASE".into(),
                    source: HorizonSource::Mapped(flat(5060.0)),
                },
            ],
            zone_layers: vec![
                StackZone::new("UPPER", Conformity::Proportional, 4, Vec::new())
                    .with_color("#ffcc00"),
                StackZone::new("LOWER", Conformity::Proportional, 4, Vec::new()),
            ],
        };
        StaticModelBuilder::from_horizon_stack(stack, opts())
            .unwrap()
            .with_boundary(vec![
                [0.0, 0.0],
                [10.0, 0.0],
                [10.0, 10.0],
                [0.0, 10.0],
                [0.0, 0.0],
            ])
            .build()
            .unwrap()
    }

    #[test]
    fn section_carries_zone_ids_and_colour_by_zone_table() {
        // The colour-by-zone rider (SCHEMA_VERSION 5, task_suite_section_zone_color):
        // the bundle carries the zone table {name, color} and each column a per-layer
        // `zone_ids` array, NaN-gapped in lockstep with the geometry/value arrays.
        let m = two_zone_stack_model();
        let spec = SectionSpec::Polyline(vec![[0.5, 0.5], [9.5, 0.5]]);
        let b = m.intersection_bundle(&spec, Some("PORO")).unwrap();

        // Zone table mirrors the model's zones, in order, with colour carried through.
        assert_eq!(
            b.zones,
            vec![
                SectionZone {
                    name: "UPPER".into(),
                    color: Some("#ffcc00".into()),
                },
                SectionZone {
                    name: "LOWER".into(),
                    color: None,
                },
            ]
        );

        assert!(!b.columns.is_empty());
        let mut saw_upper = false;
        let mut saw_lower = false;
        for col in &b.columns {
            assert_eq!(
                col.zone_ids.len(),
                col.layer_tops.len(),
                "zone_ids is nk-sized"
            );
            for k in 0..col.zone_ids.len() {
                if col.layer_tops[k].is_nan() {
                    // NaN-gapped in lockstep: an inactive layer carries the sentinel.
                    assert_eq!(
                        col.zone_ids[k],
                        SectionColumn::NO_ZONE,
                        "inactive layer -> NO_ZONE"
                    );
                } else {
                    let id = col.zone_ids[k];
                    assert!(
                        (id as usize) < b.zones.len(),
                        "active layer -> valid zone id"
                    );
                    // Layers 0..4 are UPPER (id 0), 4..8 are LOWER (id 1).
                    let expect = if k < 4 { 0 } else { 1 };
                    assert_eq!(id, expect, "layer {k} zone id");
                    saw_upper |= id == 0;
                    saw_lower |= id == 1;
                }
            }
        }
        assert!(
            saw_upper && saw_lower,
            "both zones appear along the section"
        );

        // Round-trips through serde intact (the wire contract the viewer decodes).
        let json = serde_json::to_string(&b).unwrap();
        let back: IntersectionBundle = serde_json::from_str(&json).unwrap();
        assert_eq!(b, back);
        assert_eq!(b.schema_version, 5);
    }
}
