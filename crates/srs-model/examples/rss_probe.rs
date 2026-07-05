//! Peak-RSS probe for the out-of-core build (`task_petekstatic_slab_incremental_build`).
//!
//! Builds one flat constant-prior model in one of three modes and exits, so an
//! external `/usr/bin/time -l` reads the process's maximum resident set size.
//! `incore` = in-core build (whole f64 grid); `v1` = build-then-spill (whole f64
//! grid + f32 store resident together); `v2` = slab-incremental streaming spilled
//! build (ZCORN + cubes stream k-slab-by-k-slab into the store, no whole grid).
//!
//! Usage: `rss_probe <incore|v1|v2> <ni> <nj> <nk>`. Build release:
//! `cargo build -p srs-model --release --example rss_probe`.

use srs_gridder::{Conformity, SolveOpts};
use srs_model::{spill_grid_to, BuildOpts, ConstantPriors, MemoryBudget, StaticModelBuilder};

fn opts(nk: usize) -> BuildOpts {
    BuildOpts {
        area_m2: 1.0e6,
        gross_height_m: 40.0,
        nk,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.25,
            net_to_gross: 0.8,
            water_saturation: 0.3,
        },
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 5 {
        eprintln!("usage: rss_probe <incore|v1|v2> <ni> <nj> <nk>");
        std::process::exit(2);
    }
    let (mode, ni, nj, nk) = (
        a[1].as_str(),
        a[2].parse().unwrap(),
        a[3].parse().unwrap(),
        a[4].parse().unwrap(),
    );
    let cells = ni * nj * nk;
    match mode {
        "incore" => {
            let m = StaticModelBuilder::flat(ni, nj, 2500.0, 9000.0, opts(nk))
                .unwrap()
                .with_memory_budget(MemoryBudget::unlimited())
                .build()
                .unwrap();
            std::hint::black_box(m.bulk_volume());
        }
        "v1" => {
            // The v1 forced-spill peak: build the whole in-core grid, THEN spill it.
            let m = StaticModelBuilder::flat(ni, nj, 2500.0, 9000.0, opts(nk))
                .unwrap()
                .with_memory_budget(MemoryBudget::unlimited())
                .build()
                .unwrap();
            let path = std::env::temp_dir().join(format!("rss-v1-{}.pts", std::process::id()));
            let backing = spill_grid_to(m.grid(), &path, true).unwrap();
            std::hint::black_box(backing.bulk_volume().unwrap());
        }
        "v2" => {
            // Slab-incremental streaming spilled build — a tiny budget forces spill.
            let m = StaticModelBuilder::flat(ni, nj, 2500.0, 9000.0, opts(nk))
                .unwrap()
                .with_memory_budget(MemoryBudget::bytes(1024))
                .build()
                .unwrap();
            assert!(m.is_spilled());
            std::hint::black_box(m.bulk_volume());
        }
        other => {
            eprintln!("unknown mode: {other}");
            std::process::exit(2);
        }
    }
    eprintln!("mode={mode} cells={cells}");
}
