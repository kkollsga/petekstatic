//! `volume_bench` — measure the `VolumeBundle` payload (exterior shell + binary
//! blocks) at a given grid size: bytes, effective **B/cell**, extraction time and
//! serialize time. Run under `/usr/bin/time -l` to capture peak RSS.
//!
//! ```text
//! cargo run --release --example volume_bench -- <ni> <nj> <nk> [self|sidecar]
//! /usr/bin/time -l cargo run --release --example volume_bench -- 200 200 25
//! ```
//! Synthetic data only (constant priors) — no dataset content is touched.

use srs_gridder::{Conformity, SolveOpts};
use srs_model::{BuildOpts, ConstantPriors, StaticModelBuilder};
use std::io::{self, Write};
use std::time::Instant;

/// A sink that counts bytes and discards them (payload-size probe without holding
/// the serialized output in memory).
struct CountingSink(u64);
impl Write for CountingSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0 += buf.len() as u64;
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let ni: usize = a.get(1).and_then(|s| s.parse().ok()).unwrap_or(200);
    let nj: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(200);
    let nk: usize = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(25);
    let mode = a.get(4).map(String::as_str).unwrap_or("self");
    let cells = ni * nj * nk;

    let opts = BuildOpts {
        area_m2: 1.0e7,
        gross_height_m: 100.0,
        nk,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.25,
            net_to_gross: 0.8,
            water_saturation: 0.3,
        },
    };
    let t = Instant::now();
    let model = StaticModelBuilder::flat(ni, nj, 2000.0, 2100.0, opts)
        .unwrap()
        .build()
        .unwrap();
    let build_ms = t.elapsed().as_secs_f64() * 1e3;

    let t = Instant::now();
    let vb = model.volume_bundle("PORO").unwrap();
    let extract_ms = t.elapsed().as_secs_f64() * 1e3;

    let mut sink = CountingSink(0);
    let t = Instant::now();
    match mode {
        "sidecar" => {
            let mut binsink = CountingSink(0);
            vb.write_sidecar(&mut sink, &mut binsink).unwrap();
            sink.0 += binsink.0;
        }
        _ => vb.write_self_contained(&mut sink).unwrap(),
    }
    let serialize_ms = t.elapsed().as_secs_f64() * 1e3;

    let bytes = sink.0;
    let bpc = bytes as f64 / cells as f64;
    let tris = vb.indices.len() / 3;
    let verts = vb.positions.len() / 3;
    let shell = vb.cell_values.len();
    println!(
        "grid {ni}x{nj}x{nk} = {cells} cells | mode {mode}\n\
         shell: {shell} cells, {verts} verts, {tris} tris ({:.2}% of full soup tris)\n\
         payload: {bytes} B = {:.2} MB | {bpc:.3} B/cell\n\
         time: build {build_ms:.1} ms | extract {extract_ms:.1} ms | serialize {serialize_ms:.1} ms",
        100.0 * tris as f64 / (cells * 12) as f64,
        bytes as f64 / 1.0e6,
    );
}
