//! `wire` — the streaming serializer + the volume bundle's **binary-block payload
//! spec** (SCHEMA_VERSION 3). This is the authoritative decode contract for the
//! separate viewer codebase (petekTools): every array below is emitted as raw
//! little-endian bytes; the viewer reconstructs typed arrays (`ArrayBuffer`s)
//! from exactly this layout.
//!
//! ## Why binary (the P9 payload-killer, `task_suite_bundle_binary`)
//! The legacy viewer payload was decimal-text JSON at ~519–620 B/cell — a 1M-cell
//! model serialized to ~613 MB, over V8's `kMaxStringLength` (~512 MiB) wall, so
//! it never loaded. Binary f32/u32/u16 blocks + the volume bundle's exterior-shell
//! reduction (see [`super::volume`]) collapse that by two orders of magnitude.
//!
//! ## The block spec (decode contract — DO NOT change without a SCHEMA_VERSION bump)
//! * **Endianness:** little-endian, always (independent of host).
//! * **Packing:** tight — no padding or alignment between elements or between
//!   blocks. A block of `n` elements of a `k`-byte dtype is exactly `n*k` bytes.
//! * **dtypes:** `"f32"` (4 B IEEE-754), `"u32"` (4 B), `"u16"` (2 B). A `NaN`
//!   f32 is the canonical quiet NaN `0x7FC00000`.
//! * **Order:** the row-major (C-order) flattening of `shape`. E.g. `positions`
//!   with `shape [V, 3]` is `x0,y0,z0, x1,y1,z1, …`.
//!
//! ## The two envelopes
//! Both are a JSON object carrying the human-readable metadata
//! (names/units/ranges/zones/provenance) + a `blocks` map. They differ only in
//! where the block bytes live:
//! * **`encoding: "base64"`** (self-contained) — each block object carries
//!   `"data": "<base64 of the raw LE bytes>"`. One file; a `save_view` export.
//!   [`super::VolumeBundle::write_self_contained`].
//! * **`encoding: "sidecar"`** (served) — each block object carries `"offset"`
//!   and `"length"` (bytes) into a companion `model.bin`; the blocks are written
//!   to `model.bin` tightly concatenated in the envelope-declared order. No base64
//!   overhead when HTTP serves the two files.
//!   [`super::VolumeBundle::write_sidecar`].
//!
//! ## Streaming (the RSS win)
//! Every writer emits straight to an [`Write`] — metadata via
//! [`serde_json::to_writer`], block bytes via a small fixed scratch buffer, base64
//! via [`base64::write::EncoderWriter`]. There is never a whole-payload `String`
//! or `serde_json::Value` in memory, so peak RSS stays ~1x the (already tiny)
//! payload instead of the legacy ~12.7x spike.

use super::frame::ValueRange;
use base64::engine::general_purpose::STANDARD;
use base64::write::EncoderWriter;
use serde::Serialize;
use std::io::{self, Write};

/// A typed view over one of the bundle's big arrays, streamable as little-endian
/// bytes into any [`Write`] (the raw sidecar sink, or a base64 [`EncoderWriter`]).
#[derive(Clone, Copy)]
pub(crate) enum BlockData<'a> {
    F32(&'a [f32]),
    U32(&'a [u32]),
    U16(&'a [u16]),
}

impl BlockData<'_> {
    /// The wire dtype tag.
    pub(crate) fn dtype(&self) -> &'static str {
        match self {
            BlockData::F32(_) => "f32",
            BlockData::U32(_) => "u32",
            BlockData::U16(_) => "u16",
        }
    }

    /// Total byte length of the block on the wire (tight little-endian packing).
    pub(crate) fn byte_len(&self) -> usize {
        match self {
            BlockData::F32(v) => v.len() * 4,
            BlockData::U32(v) => v.len() * 4,
            BlockData::U16(v) => v.len() * 2,
        }
    }

    /// Stream the block as tightly packed little-endian bytes. A fixed 8 KiB
    /// scratch bounds the extra memory regardless of block size (the streaming
    /// RSS guarantee); the caller passes either the raw sidecar writer or a
    /// base64 encoder wrapping the envelope writer.
    fn write_le<W: Write>(&self, w: &mut W) -> io::Result<()> {
        // 8192 is a multiple of both 4 (f32/u32) and 2 (u16), so no element ever
        // straddles a flush boundary.
        let mut buf = [0u8; 8192];
        macro_rules! stream {
            ($slice:expr, $bytes:expr) => {{
                let mut n = 0usize;
                for &x in $slice {
                    let b = x.to_le_bytes();
                    buf[n..n + $bytes].copy_from_slice(&b);
                    n += $bytes;
                    if n == buf.len() {
                        w.write_all(&buf)?;
                        n = 0;
                    }
                }
                if n > 0 {
                    w.write_all(&buf[..n])?;
                }
            }};
        }
        match self {
            BlockData::F32(v) => stream!(*v, 4),
            BlockData::U32(v) => stream!(*v, 4),
            BlockData::U16(v) => stream!(*v, 2),
        }
        Ok(())
    }
}

/// One named binary block: its wire name, logical `shape` (row-major), and the
/// typed data view.
pub(crate) struct Block<'a> {
    pub(crate) name: &'static str,
    pub(crate) shape: Vec<usize>,
    pub(crate) data: BlockData<'a>,
}

/// Write `"key":<json(val)>` (no surrounding comma) — value escaping via serde.
fn kv<W: Write, T: Serialize + ?Sized>(w: &mut W, key: &str, val: &T) -> io::Result<()> {
    serde_json::to_writer(&mut *w, key).map_err(io::Error::other)?;
    w.write_all(b":")?;
    serde_json::to_writer(&mut *w, val).map_err(io::Error::other)?;
    Ok(())
}

/// Stream a bundle straight to `w` as plain JSON (no intermediate `Value` tree) —
/// the map/section path (small bundles, no binary blocks) and the streaming
/// primitive the whole `view` seam routes through.
pub fn write_json<T: Serialize, W: Write>(bundle: &T, w: &mut W) -> io::Result<()> {
    serde_json::to_writer(w, bundle).map_err(io::Error::other)
}

/// Open the envelope and stream every metadata field (leaving the object open for
/// the caller to append `,"blocks":…`).
#[allow(clippy::too_many_arguments)]
fn write_head<W: Write>(
    w: &mut W,
    schema_version: u32,
    inputs_ref: &str,
    property: &str,
    cell_count: usize,
    shell_cell_count: usize,
    vertex_count: usize,
    triangle_count: usize,
    zone_names: &[String],
    value_range: &impl Serialize,
    encoding: &str,
) -> io::Result<()> {
    w.write_all(b"{")?;
    kv(w, "schema_version", &schema_version)?;
    w.write_all(b",")?;
    kv(w, "kind", "volume")?;
    w.write_all(b",")?;
    kv(w, "inputs_ref", inputs_ref)?;
    w.write_all(b",")?;
    kv(w, "property", property)?;
    w.write_all(b",")?;
    kv(w, "cell_count", &cell_count)?;
    w.write_all(b",")?;
    kv(w, "shell_cell_count", &shell_cell_count)?;
    w.write_all(b",")?;
    kv(w, "vertex_count", &vertex_count)?;
    w.write_all(b",")?;
    kv(w, "triangle_count", &triangle_count)?;
    w.write_all(b",")?;
    kv(w, "zone_names", zone_names)?;
    w.write_all(b",")?;
    kv(w, "value_range", value_range)?;
    w.write_all(b",")?;
    kv(w, "encoding", encoding)?;
    Ok(())
}

/// Envelope head fields, borrowed from a [`super::VolumeBundle`].
pub(crate) struct Head<'a> {
    pub(crate) schema_version: u32,
    pub(crate) inputs_ref: &'a str,
    pub(crate) property: &'a str,
    pub(crate) cell_count: usize,
    pub(crate) shell_cell_count: usize,
    pub(crate) vertex_count: usize,
    pub(crate) triangle_count: usize,
    pub(crate) zone_names: &'a [String],
    pub(crate) value_range: &'a ValueRange,
}

/// Self-contained envelope: metadata + base64-wrapped binary blocks, one file.
pub(crate) fn write_self_contained<W: Write>(
    head: &Head,
    blocks: &[Block],
    w: &mut W,
) -> io::Result<()> {
    write_head(
        w,
        head.schema_version,
        head.inputs_ref,
        head.property,
        head.cell_count,
        head.shell_cell_count,
        head.vertex_count,
        head.triangle_count,
        head.zone_names,
        head.value_range,
        "base64",
    )?;
    w.write_all(b",\"blocks\":{")?;
    for (n, blk) in blocks.iter().enumerate() {
        if n > 0 {
            w.write_all(b",")?;
        }
        serde_json::to_writer(&mut *w, blk.name).map_err(io::Error::other)?;
        w.write_all(b":{")?;
        kv(w, "dtype", blk.data.dtype())?;
        w.write_all(b",")?;
        kv(w, "shape", &blk.shape)?;
        w.write_all(b",\"data\":\"")?;
        {
            let mut enc = EncoderWriter::new(&mut *w, &STANDARD);
            blk.data.write_le(&mut enc)?;
            enc.finish()?;
        }
        w.write_all(b"\"}")?;
    }
    w.write_all(b"}}")?;
    Ok(())
}

/// Sidecar envelope: metadata + a `(offset,length)` manifest into `bin`, with the
/// raw block bytes streamed to `bin` in the declared order.
pub(crate) fn write_sidecar<W1: Write, W2: Write>(
    head: &Head,
    blocks: &[Block],
    json: &mut W1,
    bin: &mut W2,
) -> io::Result<()> {
    write_head(
        json,
        head.schema_version,
        head.inputs_ref,
        head.property,
        head.cell_count,
        head.shell_cell_count,
        head.vertex_count,
        head.triangle_count,
        head.zone_names,
        head.value_range,
        "sidecar",
    )?;
    json.write_all(b",\"blocks\":{")?;
    let mut offset = 0usize;
    for (n, blk) in blocks.iter().enumerate() {
        if n > 0 {
            json.write_all(b",")?;
        }
        serde_json::to_writer(&mut *json, blk.name).map_err(io::Error::other)?;
        json.write_all(b":{")?;
        kv(json, "dtype", blk.data.dtype())?;
        json.write_all(b",")?;
        kv(json, "shape", &blk.shape)?;
        json.write_all(b",")?;
        kv(json, "offset", &offset)?;
        json.write_all(b",")?;
        kv(json, "length", &blk.data.byte_len())?;
        json.write_all(b"}")?;
        offset += blk.data.byte_len();
    }
    json.write_all(b"}}")?;
    // Block bytes, tightly concatenated in the same order the manifest declares.
    for blk in blocks {
        blk.data.write_le(bin)?;
    }
    Ok(())
}
