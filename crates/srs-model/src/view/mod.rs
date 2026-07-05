//! `view` — typed, JSON-stable **inspection bundles** the model exports for a
//! separate viewer codebase (graph `decision_viewer_home_product`,
//! `task_petekstatic_view_bundles`). The division of labour is strict:
//!
//! > **petekStatic exports; the viewer renders — it never computes.**
//!
//! Everything a map / cross-section / 3-D volume view needs is pre-computed here
//! off the already-populated [`StaticModel`] (grid geometry + property cubes +
//! framework horizons + contacts + zones + provenance) and handed over as plain
//! serde-`Serialize` value types. The viewer scales for display (vertical
//! exaggeration, colour ramps, threshold filtering) but performs no reservoir
//! computation.
//!
//! ## The three bundles
//! - [`MapBundle`] ([`StaticModel::map_bundle`]) — areal (plan-view) layers:
//!   per-horizon depth grids, property maps (k-slice + zone-average), the areal
//!   outline, well markers, and per-contact subcrop masks.
//! - [`IntersectionBundle`] ([`StaticModel::intersection_bundle`]) — a vertical
//!   cross-section walked along a polyline or a bore trajectory: ordered columns
//!   of per-layer property + cell top/base depths, horizon traces, contacts, and
//!   (along a bore) the path's own z overlay.
//! - [`VolumeBundle`] ([`StaticModel::volume_bundle`]) — the corner-point cell
//!   **exterior shell** (only the faces bordering an inactive/absent neighbour or
//!   the grid boundary) + compact per-shell-cell property values + zone ids. The
//!   big arrays serialize as raw little-endian **binary blocks**, base64-wrapped in
//!   a JSON envelope (self-contained) or split into a `model.bin` sidecar (served).
//!   See [`wire`] for the block spec the viewer decodes.
//!
//! ## Contract stability (the serialized *shape* is the seam)
//! The viewer is a separate codebase, so the serialized **structure** is a
//! contract — documented per bundle in `API.md` and locked by a schema-snapshot
//! test. Every bundle carries [`SCHEMA_VERSION`]; a breaking change bumps it. The
//! map/section bundles are plain JSON; the volume bundle's envelope wraps binary
//! blocks (SCHEMA_VERSION 4, [`wire`]).
//!
//! ## Streaming serialization ([`wire`])
//! All bundles serialize straight to an [`std::io::Write`] with no intermediate
//! `serde_json::Value` tree (the legacy Value-tree path was the ~12.7x-payload RSS
//! spike + the ~3 MB/s serializer). Map/section stream via [`wire::write_json`];
//! the volume bundle streams its envelope + binary blocks via
//! [`VolumeBundle::write_self_contained`] / [`VolumeBundle::write_sidecar`].
//!
//! ## Conventions (family SI + house style)
//! - **Units:** SI throughout — depths/distances in metres (positive-down),
//!   coordinates in the model's world frame; each layer names its unit label.
//! - **Georeference:** areal layers share one regular [`GridFrame`] reconstructed
//!   from the grid's column-centroid lattice (the `xy↔ij` map). Row-major indexing
//!   is `values[j * ncol + i]`, `i` along x (`ncol == ni`), `j` along y
//!   (`nrow == nj`); `NaN` marks an undefined / outside-boundary node.
//! - **`f64::NAN` = undefined** (serializes to JSON `null`).

mod frame;
mod map;
mod section;
mod volume;
pub(crate) mod wire;

pub use frame::{GridFrame, ScalarLayer, ValueRange};
pub use map::{ContactMask, MapBundle, MapSpec, WellMarker, WellTieResidual};
pub use section::{
    HorizonTrace, IntersectionBundle, SectionColumn, SectionContact, SectionSpec, SectionZone,
};
pub use volume::VolumeBundle;

/// The view-bundle schema version. Bumped on a breaking change to any bundle's
/// serialized **structure** (the cross-codebase contract). **v3** introduced the
/// volume bundle's exterior-shell + binary-block payload (see [`wire`]). **v4**
/// (`task_petekstatic_multizone_2`) adds **additive** blocks — all new fields on
/// existing bundles, so a v3 decoder that ignores unknown keys still reads a v4
/// payload (the viewer is asked to render the new blocks):
/// - the section bundle's per-interior-horizon [`section::HorizonTrace`] polylines;
/// - the map bundle's per-well per-horizon tie residuals;
/// - the section bundle's **sugar-cube** block: [`section::IntersectionBundle::sugar_cube`]
///   (root `bool`) + four per-column edge arrays [`section::SectionColumn::layer_tops_l`]
///   / `layer_tops_r` / `layer_bases_l` / `layer_bases_r` (cell interval bilinearly
///   interpolated at the column's left/right fence edges, NaN-gapped like `layer_tops`).
///   The default renders dip-following trapezoids; `sugar_cube=true` flattens the edge
///   arrays to the centroid trace. **Frozen field names** (a concurrent viewer consumes
///   them) — no rename without a coordinated bump.
pub const SCHEMA_VERSION: u32 = 5;
