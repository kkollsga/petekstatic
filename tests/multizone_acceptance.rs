//! Acceptance fixture for the multi-zone regional-framework build
//! (`task_petekstatic_multizone`). Per graph policy `decision_public_fixtures_synthetic`
//! the framework is **procedurally generated** — seeded, deterministic
//! (bit-reproducible), from our own geostatistics, at an arbitrary fictional study
//! area. The owner recipe (a fixture factory `build_framework(n, intervals, seed)`):
//!
//!   1. a synthetic TOP surface = a regional dip + a closure, plus **spatially
//!      correlated** Gaussian noise (petektools SGS, not white noise);
//!   2. one **Gaussian-simulated ISOCHORE** per interval (SGS with a per-interval
//!      mean/variogram), **clamped ≥ 0** so zeros become natural pinch-outs;
//!   3. **build down**: horizon `k+1 = horizon k + isochore k`, stacking from the
//!      top — ordered by construction (no crossing possible), realistic thickness
//!      variation, and pinch-outs that exercise truncation + the collapse threshold.
//!
//! SGS is unconditional here, so the normal-score transform is anchored on a handful
//! of **seeded pseudo-points** drawn (deterministically) from each field's target
//! distribution — the documented nscore-anchoring trick.
//!
//! It mirrors the owner-spec SHAPE — 11 horizons → 10 zones — with one untied
//! internal horizon that is **tops-only** (4 fictional well picks, no mapped
//! surface), a per-zone conformity mix, and per-zone contacts (a single-OWC zone, a
//! two-contact GOC+OWC zone, and contactless zones). All names/ids are neutral
//! fictional placeholders — no dataset content. Since the surfaces are
//! non-conformable, per-zone volumes are hand-checked against an **independent
//! in-test geometry reference** + conservation invariants.
//!
//! (This is also an indirect exercise of the SGS kernel — see the module tail for a
//! note on the one bit of kernel friction hit while wiring it up.)
//!
//! **Mode matrix (doctrine R2).** This file owns the **horizon-stack (multizone)**
//! column of the support matrix: per-zone `in_place`/`by_zone` (in-core `:924`,
//! spilled `:1000`), zone-aware map/section/volume bundles (`:1067`), and stack
//! construction/determinism. The full cross-feature matrix + the single-zone,
//! spilled-MC, and unsupported-cell cells live in `mode_matrix.rs`.

use petekstatic::grid::Ijk;
use petekstatic::gridder::{Conformity, SolveOpts};
use petekstatic::model::{
    BuildOpts, ConstantPriors, HorizonSource, HorizonStack, MapSpec, MemoryBudget, Pick,
    PropertyPipeline, SectionSpec, StackHorizon, StackZone, StaticModel, StaticModelBuilder,
    UpscaleMethod, WellLog, WellTie,
};
use petekstatic::volumetrics::{NTG, PORO, SW};
use petekstatic::wireframe::{Contact, ContactKind, GriddedDepth, Hardness};

use petektools::geostat::{sgs, SgsParams};
use petektools::{Lattice, Variogram, VariogramModel};

const SEED: u64 = 20_260_704; // fixed → the procedural framework reproduces bit-for-bit
const N_NODES: usize = 9; // 8×8 cells
const NI: usize = N_NODES - 1;
const AREA_M2: f64 = 1_000_000.0; // side 1000 m
const BASE_DEPTH: f64 = 2000.0;
const SIDE: f64 = 1000.0;
const DX: f64 = SIDE / NI as f64;
// A fictional study-area coordinate window (NOT the 43xxxx/652xxxx vicinity used
// elsewhere).
const ORIGIN_X: f64 = 700_000.0;
const ORIGIN_Y: f64 = 7_100_000.0;

// The 10 interval (isochore) targets: (mean thickness, sd). Zone 8 is deliberately
// thin + high-variance so its clamped isochore pinches out in places (truncation).
const INTERVALS: [(f64, f64); 10] = [
    (12.0, 3.0),
    (13.0, 3.0),
    (9.0, 2.5),
    (14.0, 3.5),
    (12.0, 3.0),
    (11.0, 3.0),
    (14.0, 3.5),
    (7.0, 2.0),
    (5.0, 4.0), // pinch-prone
    (12.0, 3.0),
];

const PHI: f64 = 0.2;
const NTG_V: f64 = 1.0;
const SW_V: f64 = 0.25;
const HC: f64 = NTG_V * PHI * (1.0 - SW_V); // 0.15 hydrocarbon per m³ of GRV

fn lattice() -> Lattice {
    Lattice::regular(ORIGIN_X, ORIGIN_Y, DX, DX, N_NODES, N_NODES)
}

/// A spatially-correlated Gaussian field over the lattice via petektools SGS, in
/// **data space** ~ N(mean, sd) with range `range`. The normal-score transform is
/// anchored on 9 **seeded pseudo-points** whose values are a deterministic spread
/// around the target (the documented nscore-anchoring trick for unconditional
/// simulation). Returned row-major `jp * N_NODES + ip`.
fn sgs_field(mean: f64, sd: f64, range: f64, seed: u64) -> Vec<f64> {
    let lat = lattice();
    // Deterministic spread → a non-degenerate data distribution for nscore.
    let spread = [-1.5, -0.9, -0.4, -0.1, 0.1, 0.4, 0.9, 1.5, 0.0];
    let locs = [
        (0.1, 0.1),
        (0.9, 0.1),
        (0.5, 0.5),
        (0.1, 0.9),
        (0.9, 0.9),
        (0.3, 0.7),
        (0.7, 0.3),
        (0.5, 0.1),
        (0.5, 0.9),
    ];
    let coords: Vec<[f64; 3]> = locs
        .iter()
        .zip(spread.iter())
        .map(|(&(fx, fy), &s)| [ORIGIN_X + fx * SIDE, ORIGIN_Y + fy * SIDE, mean + sd * s])
        .collect();
    let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, range).unwrap();
    let params = SgsParams::new(vgm, 16, SIDE * 2.0, seed).unwrap();
    let field = sgs(&coords, &lat, &params).expect("sgs field");
    let mut out = vec![0.0; N_NODES * N_NODES];
    for jp in 0..N_NODES {
        for ip in 0..N_NODES {
            out[jp * N_NODES + ip] = field[[ip, jp]];
        }
    }
    out
}

/// Crest-to-spill relief of the top-surface closure [m].
const CREST_RELIEF: f64 = 80.0;

/// The top-surface trend (no noise): a **4-way dip closure** — an elongated (~5:3)
/// elliptical dome closed in all directions, with a gentle superimposed regional
/// tilt. `BASE_DEPTH` is the flat regional level; the crest sits `CREST_RELIEF`
/// shallower, the Gaussian decaying to the regional plane well inside the extent so
/// the closing contours (and the spill point) lie within the study area.
fn top_trend(ip: usize, jp: usize) -> f64 {
    let (x, y) = (ip as f64 / NI as f64, jp as f64 / NI as f64);
    let tilt = 8.0 * x + 4.0 * y; // gentle regional tilt (~12 m across; << relief)
    let (dxr, dyr) = (x - 0.5, y - 0.5);
    // Elliptical dome, ~5:3 aspect (sx : sy ≈ 1.68).
    let (ex, ey) = (dxr / 0.32, dyr / 0.19);
    let dome = -CREST_RELIEF * (-(ex * ex + ey * ey)).exp();
    BASE_DEPTH + tilt + dome
}

/// Build the framework node-fields down from the top (the reusable factory). Returns
/// `n_horizons` surfaces (row-major node depth fields), ordered by construction:
/// `surf[k+1] = surf[k] + max(isochore_k, 0)`.
fn build_framework(intervals: &[(f64, f64)], seed: u64) -> Vec<Vec<f64>> {
    // (1) TOP = the 4-way dip closure + SUBTLE spatially-correlated SGS noise
    //     (sigma ~4% of crest-to-spill relief so texture reads without lumpiness).
    let mut top = vec![0.0; N_NODES * N_NODES];
    for jp in 0..N_NODES {
        for ip in 0..N_NODES {
            top[jp * N_NODES + ip] = top_trend(ip, jp);
        }
    }
    let noise = sgs_field(0.0, 0.04 * CREST_RELIEF, 300.0, seed);
    for (t, n) in top.iter_mut().zip(&noise) {
        *t += n;
    }

    // (2)+(3) Gaussian-simulated isochores, clamped ≥ 0, built down.
    let mut surfaces = vec![top];
    for (k, &(mean, sd)) in intervals.iter().enumerate() {
        let iso = sgs_field(mean, sd, 400.0, seed + 101 + k as u64);
        let prev = surfaces.last().unwrap();
        // build down: clamp isochore ≥ 0 → pinch-outs.
        let next: Vec<f64> = prev.iter().zip(&iso).map(|(p, i)| p + i.max(0.0)).collect();
        surfaces.push(next);
    }
    surfaces
}

fn gridded(field: &[f64]) -> GriddedDepth {
    GriddedDepth {
        ncol: N_NODES,
        nrow: N_NODES,
        depth_m: field.to_vec(),
        is_control: vec![true; N_NODES * N_NODES],
    }
}

fn owc(depth: f64) -> Contact {
    Contact {
        kind: ContactKind::Owc,
        depth_m: depth,
        hardness: Hardness::Interpolated,
    }
}
fn goc(depth: f64) -> Contact {
    Contact {
        kind: ContactKind::Goc,
        depth_m: depth,
        hardness: Hardness::Interpolated,
    }
}

// ---------------------------------------------------------------------------
// The full synthetic asset: framework + WELLS + TOPS + zone TARGETS + LOGS +
// contacts. The INTERFACE is the deliverable here (the dedicated synth wave
// upgrades the generators in-place behind these signatures).
// ---------------------------------------------------------------------------

/// Fixture instantiation parameters (factory params, per the owner extension).
const N_WELLS: usize = 6;
const MAX_RESIDUAL_M: f64 = 8.0; // uniform ± this, injected into every well top

/// A tiny seeded deterministic RNG (splitmix64 + Box-Muller) — enough for the
/// first-cut generators; the synth wave swaps in the petektools samplers behind
/// the same `synth_*` signatures.
struct Rng64(u64);
impl Rng64 {
    fn new(seed: u64) -> Self {
        Self(
            seed.wrapping_mul(0x9E37_79B9_7F4A_7C15)
                .wrapping_add(0x1234_5678),
        )
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn next_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
    fn next_gauss(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// A fictional synthetic well at an interior lattice node.
struct SynthWell {
    id: String,
    ip: usize,
    jp: usize,
}

impl SynthWell {
    /// Local grid coordinates (the frame the grid/property pipeline live in): the
    /// **centroid** of areal column `(ip, jp)` — the well sits in the middle of its
    /// column, not on a node corner. The log-upscale now bins each sample against the
    /// cell interval interpolated *at the well's position*, so a well at the column
    /// centroid bins exactly as the 4-corner-mean interval `zone_depth_grid` authors
    /// its samples against (`task_petekstatic_zoned_fixes` finding 3).
    fn local_xy(&self) -> (f64, f64) {
        ((self.ip as f64 + 0.5) * DX, (self.jp as f64 + 0.5) * DX)
    }
}

/// Per-zone petrophysical targets — the fixture's spec object (sand vs shale
/// contrast). The generators + the acceptance assertions both read this.
#[derive(Debug, Clone, Copy)]
struct ZoneTargets {
    por_mean: f64,
    por_std: f64,
    ntg_mean: f64,
    ntg_std: f64,
}

const SAND: ZoneTargets = ZoneTargets {
    por_mean: 0.28,
    por_std: 0.05,
    ntg_mean: 0.85,
    ntg_std: 0.10,
};
const SHALE: ZoneTargets = ZoneTargets {
    por_mean: 0.08,
    por_std: 0.03,
    ntg_mean: 0.15,
    ntg_std: 0.08,
};

/// Per-zone targets: the two hydrocarbon zones (Z2, Z7) are sands; shaly zones
/// interleave for contrast.
fn zone_targets() -> [ZoneTargets; 10] {
    [
        SAND, SHALE, SAND, SHALE, SAND, SHALE, SAND, SAND, SHALE, SAND,
    ]
}

/// Seeded random well placement on INTERIOR nodes (factory params: `n_wells`,
/// `seed`). Wells 0 and 1 are pinned by design: one **near the crest**, one **on
/// the flank** (so ties + contacts are exercised across the structure); the rest
/// are seeded-random, deduplicated.
fn place_wells(n_wells: usize, seed: u64) -> Vec<SynthWell> {
    let interior = 1..N_NODES - 1;
    // Crest = shallowest interior node of the trend; flank = deepest.
    let (mut crest, mut crest_d) = ((1, 1), f64::INFINITY);
    let (mut flank, mut flank_d) = ((1, 1), f64::NEG_INFINITY);
    for jp in interior.clone() {
        for ip in interior.clone() {
            let d = top_trend(ip, jp);
            if d < crest_d {
                crest_d = d;
                crest = (ip, jp);
            }
            if d > flank_d {
                flank_d = d;
                flank = (ip, jp);
            }
        }
    }
    let mut nodes = vec![crest, flank];
    let mut rng = Rng64::new(seed ^ 0x5745_4C4C); // "WELL"
    while nodes.len() < n_wells {
        let ip = 1 + (rng.next_u64() as usize) % (N_NODES - 2);
        let jp = 1 + (rng.next_u64() as usize) % (N_NODES - 2);
        if !nodes.contains(&(ip, jp)) {
            nodes.push((ip, jp));
        }
    }
    nodes
        .into_iter()
        .enumerate()
        .map(|(w, (ip, jp))| SynthWell {
            id: format!("99/{}-{}", w + 1, 1),
            ip,
            jp,
        })
        .collect()
}

/// TOPS generated FROM the synthetic surfaces: sample each horizon at the well
/// node and inject a seeded uniform ±`max_residual_m` depth residual. Returns
/// `(tops, residuals)`, both `[well][horizon]`.
fn well_tops(
    wells: &[SynthWell],
    surfaces: &[Vec<f64>],
    seed: u64,
    max_residual_m: f64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut rng = Rng64::new(seed ^ 0x7075);
    let mut tops = Vec::with_capacity(wells.len());
    let mut residuals = Vec::with_capacity(wells.len());
    for w in wells {
        let (mut t, mut r) = (Vec::new(), Vec::new());
        for surf in surfaces {
            let res = rng.next_range(-max_residual_m, max_residual_m);
            t.push(surf[w.jp * N_NODES + w.ip] + res);
            r.push(res);
        }
        tops.push(t);
        residuals.push(r);
    }
    (tops, residuals)
}

/// First-cut synthetic porosity curve for one zone at one well: an AR(1)-smoothed
/// (spatially correlated, NOT white-noise) Gaussian series hitting the zone
/// targets, clamped to `[0, 1]`. `depth_grid` is the sample depths (metres); one
/// value per depth. **Clean seam** — the dedicated synth wave upgrades this
/// in-place (same signature).
fn synth_por(targets: &ZoneTargets, depth_grid: &[f64], seed: u64) -> Vec<f64> {
    synth_curve(targets.por_mean, targets.por_std, depth_grid.len(), seed)
}

/// First-cut synthetic NTG curve — same contract as [`synth_por`].
fn synth_ntg(targets: &ZoneTargets, depth_grid: &[f64], seed: u64) -> Vec<f64> {
    synth_curve(
        targets.ntg_mean,
        targets.ntg_std,
        depth_grid.len(),
        seed ^ 0x004E_5447, // "NTG"
    )
}

/// The shared AR(1) generator: `g_t = a·g_{t-1} + sqrt(1-a²)·ε_t` (stationary unit
/// variance), `x = mean + std·g`, clamped to `[0, 1]`.
fn synth_curve(mean: f64, std: f64, n: usize, seed: u64) -> Vec<f64> {
    const A: f64 = 0.7;
    let mut rng = Rng64::new(seed);
    let mut g = rng.next_gauss();
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        out.push((mean + std * g).clamp(0.0, 1.0));
        g = A * g + (1.0 - A * A).sqrt() * rng.next_gauss();
    }
    out
}

/// The depth sample grid for one zone at one well: 0.25 m steps across the zone's
/// depth band at the well's areal COLUMN — the 4-node average of the bounding
/// horizons over the column's corner nodes, matching how a corner-point cell
/// defines its depth range (`Cell::top_depth` = mean of the 4 top corners), so the
/// engine's depth binning assigns the samples to the intended zone.
fn zone_depth_grid(surfaces: &[Vec<f64>], z: usize, well: &SynthWell) -> Vec<f64> {
    // An interior node (ip, jp) snaps to areal column (ip, jp): corners at nodes
    // {ip, ip+1} × {jp, jp+1}.
    let col_avg = |surf: &[f64]| {
        let (i, j) = (well.ip, well.jp);
        (surf[j * N_NODES + i]
            + surf[j * N_NODES + i + 1]
            + surf[(j + 1) * N_NODES + i]
            + surf[(j + 1) * N_NODES + i + 1])
            / 4.0
    };
    let (top, base) = (col_avg(&surfaces[z]), col_avg(&surfaces[z + 1]));
    let mut d = top;
    let mut out = Vec::new();
    while d <= base {
        out.push(d);
        d += 0.25;
    }
    out
}

/// The complete synthetic asset: framework + wells + tops (+injected residuals) +
/// per-zone targets. H5's tops-only picks come from the first 4 (tie) wells' tops
/// — REAL residual-bearing conditioning, so the drape pins H5 to the measured tops.
struct SyntheticAsset {
    stack: HorizonStack,
    surfaces: Vec<Vec<f64>>,
    wells: Vec<SynthWell>,
    tops: Vec<Vec<f64>>,
    residuals: Vec<Vec<f64>>,
    targets: [ZoneTargets; 10],
}

impl SyntheticAsset {
    /// Indices of the wells whose H5 tops condition the tops-only drape.
    fn tie_wells(&self) -> std::ops::Range<usize> {
        0..4.min(self.wells.len())
    }

    /// One [`WellLog`] per well for `property` (`PORO`/`NTG`), concatenating the
    /// per-zone synthetic curves down-hole (targets per zone).
    fn well_logs(&self, property: &str, seed: u64) -> Vec<WellLog> {
        self.wells
            .iter()
            .enumerate()
            .map(|(w, well)| {
                let mut samples = Vec::new();
                for z in 0..10 {
                    let grid = zone_depth_grid(&self.surfaces, z, well);
                    let zseed = seed ^ ((w as u64) << 32) ^ ((z as u64) << 16);
                    let vals = match property {
                        PORO => synth_por(&self.targets[z], &grid, zseed),
                        _ => synth_ntg(&self.targets[z], &grid, zseed),
                    };
                    samples.extend(grid.iter().copied().zip(vals));
                }
                let (x, y) = well.local_xy();
                WellLog::new(x, y, samples)
            })
            .collect()
    }
}

/// The asset factory (`n_wells`, `seed`, `max_residual_m` are the owner's factory
/// params).
fn synthetic_asset(n_wells: usize, seed: u64, max_residual_m: f64) -> SyntheticAsset {
    let surfaces = build_framework(&INTERVALS, seed);
    let wells = place_wells(n_wells, seed);
    let (tops, residuals) = well_tops(&wells, &surfaces, seed, max_residual_m);
    let mean_h7_top = mean_of(&surfaces[7]);

    let mut horizons = Vec::with_capacity(11);
    for (i, surf) in surfaces.iter().enumerate() {
        if i == 5 {
            // The untied internal split: tops-only, conditioned on the tie wells'
            // MEASURED (residual-bearing) H5 tops.
            let picks = (0..4.min(wells.len()))
                .map(|w| Pick {
                    ip: wells[w].ip,
                    jp: wells[w].jp,
                    depth_m: tops[w][5],
                })
                .collect();
            horizons.push(StackHorizon {
                name: "H5".into(),
                source: HorizonSource::TopsOnly(picks),
            });
        } else {
            horizons.push(StackHorizon {
                name: format!("H{i}"),
                source: HorizonSource::Mapped(gridded(surf)),
            });
        }
    }

    let mut zone_layers = Vec::with_capacity(10);
    for (z, &(mean, _)) in INTERVALS.iter().enumerate() {
        let conf = if z == 3 {
            Conformity::FollowTop { dz_m: 1.0 }
        } else {
            Conformity::Proportional
        };
        let contacts = match z {
            2 => vec![owc(1.0e7)], // below all → whole zone oil
            7 => vec![goc(mean_h7_top + 3.0), owc(mean_h7_top + 5.0)], // splits inside zone 7
            _ => vec![],           // contactless
        };
        // ~1 m proportional layering (owner default); zone 8 may pinch out.
        let nk = (mean.round() as usize).max(1);
        zone_layers.push(StackZone {
            name: format!("Z{z}"),
            color: None,
            conformity: conf,
            nk,
            contacts,
        });
    }

    SyntheticAsset {
        stack: HorizonStack {
            horizons,
            zone_layers,
        },
        surfaces,
        wells,
        tops,
        residuals,
        targets: zone_targets(),
    }
}

/// The canonical fixture stack (the asset's framework at the default params).
fn regional_stack() -> HorizonStack {
    synthetic_asset(N_WELLS, SEED, MAX_RESIDUAL_M).stack
}

fn mean_of(field: &[f64]) -> f64 {
    field.iter().sum::<f64>() / field.len() as f64
}

fn opts() -> BuildOpts {
    BuildOpts {
        area_m2: AREA_M2,
        gross_height_m: 100.0,
        nk: 1,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: PHI,
            net_to_gross: NTG_V,
            water_saturation: SW_V,
        },
    }
}

fn build() -> StaticModel {
    build_with_budget(MemoryBudget::unlimited())
}

fn build_with_budget(budget: MemoryBudget) -> StaticModel {
    StaticModelBuilder::from_horizon_stack(regional_stack(), opts())
        .unwrap()
        // A small min-thickness keeps the tops-only drape ordered inside its zone.
        .with_min_thickness_m(0.0)
        .with_georef(ORIGIN_X, ORIGIN_Y, DX, DX)
        .with_memory_budget(budget)
        .build()
        .unwrap()
}

fn approx(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs().max(1.0)
}

/// Independent in-test geometry reference: per-zone bulk + gas/oil GRV, recomputed
/// straight off the built grid by centroid classification — a check on
/// `in_place_by_zone` independent of it (works for any, non-conformable, geometry).
fn reference(
    m: &StaticModel,
    k_range: std::ops::Range<usize>,
    goc: Option<f64>,
    owc: Option<f64>,
) -> (f64, f64, f64) {
    let dims = m.grid().dims();
    let (mut bulk, mut gas, mut oil) = (0.0, 0.0, 0.0);
    for k in k_range {
        for j in 0..dims.nj {
            for i in 0..dims.ni {
                let c = m.grid().cell(Ijk::new(i, j, k));
                let (z, v) = (c.centroid().z, c.volume());
                if v <= 0.0 {
                    continue;
                }
                bulk += v;
                match (goc, owc) {
                    (Some(g), Some(w)) => {
                        if z < g {
                            gas += v
                        } else if z < w {
                            oil += v
                        }
                    }
                    (None, Some(w)) => {
                        if z < w {
                            oil += v
                        }
                    }
                    _ => oil += v,
                }
            }
        }
    }
    (bulk, gas, oil)
}

#[test]
fn top_surface_is_a_four_way_dip_closure() {
    // The crest (shallowest node) is INTERIOR, and every boundary node is deeper than
    // a mid-relief contour — so that contour is a closed ring inside the extent,
    // closing the structure in all four directions (a proper 4-way trap). Tested on
    // the noiseless trend so the geometry is exact; the subtle SGS noise preserves it.
    let mut crest = f64::INFINITY;
    let mut crest_ij = (0, 0);
    for jp in 0..N_NODES {
        for ip in 0..N_NODES {
            let d = top_trend(ip, jp);
            if d < crest {
                crest = d;
                crest_ij = (ip, jp);
            }
        }
    }
    let (ci, cj) = crest_ij;
    assert!(
        ci > 0 && ci < N_NODES - 1 && cj > 0 && cj < N_NODES - 1,
        "crest must be interior, got {crest_ij:?}"
    );
    // Shallowest boundary node ≈ the spill level.
    let mut spill = f64::INFINITY;
    for jp in 0..N_NODES {
        for ip in 0..N_NODES {
            if ip == 0 || jp == 0 || ip == N_NODES - 1 || jp == N_NODES - 1 {
                spill = spill.min(top_trend(ip, jp));
            }
        }
    }
    assert!(
        spill - crest > 40.0,
        "generous closure relief: spill {spill} crest {crest}"
    );
    // A contour at mid relief is deeper than EVERY boundary node → a closed interior
    // ring encircling the crest (4-way closure).
    let contour = crest + 0.5 * (spill - crest);
    for jp in 0..N_NODES {
        for ip in 0..N_NODES {
            if ip == 0 || jp == 0 || ip == N_NODES - 1 || jp == N_NODES - 1 {
                assert!(
                    top_trend(ip, jp) > contour,
                    "boundary node ({ip},{jp}) breaches the closing contour"
                );
            }
        }
    }
}

#[test]
fn procedural_framework_is_ordered_by_construction() {
    let surfaces = build_framework(&INTERVALS, SEED);
    assert_eq!(surfaces.len(), 11);
    // Build-down (isochore ≥ 0) → every horizon sits at/below the one above at every
    // node: ordered by construction, no crossing possible.
    for k in 0..10 {
        for (n, (&lo, &hi)) in surfaces[k].iter().zip(&surfaces[k + 1]).enumerate() {
            assert!(
                hi >= lo - 1e-9,
                "horizon {k}->{} crosses at node {n}",
                k + 1
            );
        }
    }
    // The pinch-prone interval (zone 8) actually reaches ~0 somewhere.
    let z8_min_iso = (0..N_NODES * N_NODES)
        .map(|n| surfaces[9][n] - surfaces[8][n])
        .fold(f64::INFINITY, f64::min);
    assert!(
        z8_min_iso < 1.0,
        "zone 8 isochore pinches (min {z8_min_iso})"
    );
}

#[test]
fn procedural_stack_structure_and_zone_table() {
    let m = build();
    assert_eq!(m.zones().zones().len(), 10);
    let mut start = 0;
    for (z, zone) in m.zones().zones().iter().enumerate() {
        assert_eq!(zone.name, format!("Z{z}"));
        assert_eq!(zone.top_horizon, format!("H{z}"));
        assert_eq!(zone.base_horizon, format!("H{}", z + 1));
        assert_eq!(zone.k_range.start, start);
        start = zone.k_range.end;
    }
    assert_eq!(start, m.provenance().nk, "zones tile [0, nk)");
    assert_eq!(m.framework().horizons.len(), 11);

    // Provenance records the stack; zone 3 carries its FollowTop conformity.
    let stack = m.provenance().stack.as_ref().expect("stack provenance");
    assert_eq!(stack.zones.len(), 10);
    assert_eq!(
        stack.zones[3].conformity,
        Conformity::FollowTop { dz_m: 1.0 }
    );

    // Tops-only H5 honours its (residual-bearing) tie-well picks exactly and sits
    // within its zone (H4 ≤ H5 ≤ H6; the order-repair guarantees it).
    let asset = synthetic_asset(N_WELLS, SEED, MAX_RESIDUAL_M);
    let h5 = m
        .framework()
        .horizons
        .iter()
        .find(|h| h.name == "H5")
        .unwrap();
    for w in asset.tie_wells() {
        let well = &asset.wells[w];
        let got = h5.surface.depth_m[well.jp * N_NODES + well.ip];
        assert!(
            approx(got, asset.tops[w][5], 1e-4),
            "H5 pick at {} = {got} vs top {}",
            well.id,
            asset.tops[w][5]
        );
    }
    let (h4, h6) = (
        &m.framework().horizons[4].surface,
        &m.framework().horizons[6].surface,
    );
    for n in 0..N_NODES * N_NODES {
        assert!(h5.surface.depth_m[n] >= h4.depth_m[n] - 1e-6);
        assert!(h5.surface.depth_m[n] <= h6.depth_m[n] + 1e-6);
    }
}

// ---------------------------------------------------------------------------
// Wells + tops + ties + logs acceptance (the full synthetic asset)
// ---------------------------------------------------------------------------

#[test]
fn wells_cover_crest_and_flank_and_tie_residuals_verify() {
    let asset = synthetic_asset(N_WELLS, SEED, MAX_RESIDUAL_M);
    assert_eq!(asset.wells.len(), N_WELLS);

    // Placement: well 0 near the crest, well 1 on the flank — a real structural
    // spread (crest markedly shallower than flank on the trend).
    let d = |w: &SynthWell| top_trend(w.ip, w.jp);
    assert!(
        d(&asset.wells[1]) - d(&asset.wells[0]) > 30.0,
        "crest well {} vs flank well {}: insufficient structural spread",
        d(&asset.wells[0]),
        d(&asset.wells[1])
    );
    // All wells interior.
    for w in &asset.wells {
        assert!(w.ip >= 1 && w.ip <= N_NODES - 2 && w.jp >= 1 && w.jp <= N_NODES - 2);
    }

    let m = StaticModelBuilder::from_horizon_stack(asset.stack.clone(), opts())
        .unwrap()
        .with_min_thickness_m(0.0)
        .build()
        .unwrap();

    // The tie report, computed off the built framework: residual = measured top −
    // model surface at the well node.
    let model_surf =
        |h: usize, w: &SynthWell| m.framework().horizons[h].surface.depth_m[w.jp * N_NODES + w.ip];

    // (a) MAPPED horizons are pinned to their own (clean) surfaces, so the tie
    //     residual must equal the INJECTED residual exactly — real, non-zero.
    let mut saw_nonzero = false;
    for (w, well) in asset.wells.iter().enumerate() {
        for h in 0..11 {
            if h == 5 {
                continue;
            }
            let residual = asset.tops[w][h] - model_surf(h, well);
            assert!(
                approx(residual, asset.residuals[w][h], 1e-4),
                "{} H{h}: residual {residual} != injected {}",
                well.id,
                asset.residuals[w][h]
            );
            if residual.abs() > 1.0 {
                saw_nonzero = true;
            }
        }
    }
    assert!(
        saw_nonzero,
        "the injected residuals must be genuinely non-zero"
    );

    // (b) H5 (tops-only) is CONDITIONED on the tie wells' measured tops — the
    //     warm-start conditioning pins them, so the residual there is ~0 even
    //     though the injected residual was not.
    for w in asset.tie_wells() {
        let well = &asset.wells[w];
        let residual = asset.tops[w][5] - model_surf(5, well);
        assert!(
            residual.abs() < 1e-4,
            "{} H5: tie residual {residual} not pinned (injected {})",
            well.id,
            asset.residuals[w][5]
        );
        assert!(
            asset.residuals[w][5].abs() > 1e-3,
            "the pinning is only meaningful if the injected residual was non-zero"
        );
    }
}

#[test]
fn synthetic_logs_hit_zone_targets_and_are_correlated() {
    // Generator acceptance: per zone per well the AR-smoothed series hits the
    // target statistics and is spatially correlated (NOT white noise).
    let asset = synthetic_asset(N_WELLS, SEED, MAX_RESIDUAL_M);
    for (w, well) in asset.wells.iter().enumerate() {
        for z in 0..10 {
            let grid = zone_depth_grid(&asset.surfaces, z, well);
            if grid.len() < 30 {
                continue; // pinched zone at this well — too short to test stats
            }
            let zseed = SEED ^ ((w as u64) << 32) ^ ((z as u64) << 16);
            let t = &asset.targets[z];
            for (vals, mean, std) in [
                (synth_por(t, &grid, zseed), t.por_mean, t.por_std),
                (synth_ntg(t, &grid, zseed), t.ntg_mean, t.ntg_std),
            ] {
                let n = vals.len() as f64;
                let m = vals.iter().sum::<f64>() / n;
                // AR(1) a=0.7 shrinks the effective sample size ~5.7×, so the mean
                // tolerance scales with the target sigma.
                assert!(
                    (m - mean).abs() < 0.03 + 0.6 * std,
                    "{} Z{z}: series mean {m} vs target {mean} (sd {std})",
                    well.id
                );
                let var = vals.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / n;
                assert!(
                    var.sqrt() > 0.3 * std,
                    "{} Z{z}: series not variable",
                    well.id
                );
                assert!(
                    vals.iter().all(|v| (0.0..=1.0).contains(v)),
                    "clamped to [0,1]"
                );
                // Lag-1 autocorrelation well above white noise (AR(1), a = 0.7).
                let mut cov = 0.0;
                for k in 1..vals.len() {
                    cov += (vals[k] - m) * (vals[k - 1] - m);
                }
                let ac1 = cov / ((n - 1.0) * var);
                assert!(
                    ac1 > 0.15,
                    "{} Z{z}: lag-1 autocorr {ac1} — white noise is not acceptable",
                    well.id
                );
            }
        }
    }
}

#[test]
fn per_zone_upscaled_log_means_hit_zone_targets() {
    // Engine acceptance: feed the synthetic logs through the per-property upscale
    // (no propagate — only well-column cells conditioned), then check the per-zone
    // upscaled means (zone_stats over the conditioned cells) against the targets.
    let asset = synthetic_asset(N_WELLS, SEED, MAX_RESIDUAL_M);
    let m = StaticModelBuilder::from_horizon_stack(asset.stack.clone(), opts())
        .unwrap()
        .with_min_thickness_m(0.0)
        .with_property(
            PropertyPipeline::new(PORO)
                .upscale(asset.well_logs(PORO, SEED), UpscaleMethod::Arithmetic),
        )
        .with_property(
            PropertyPipeline::new(NTG)
                .upscale(asset.well_logs(NTG, SEED), UpscaleMethod::Arithmetic),
        )
        .build()
        .unwrap();

    // zone_stats sees only the conditioned (finite) cells — the upscaled well
    // columns — so its per-zone mean is the upscaled zone mean at the wells.
    let por = m.zone_stats(PORO).unwrap();
    let ntg = m.zone_stats(NTG).unwrap();
    let mut zones_checked = 0;
    for z in 0..10 {
        let t = &asset.targets[z];
        if por[z].count < 8 {
            continue; // pinched zone: too few conditioned cells to average
        }
        // Zones 4 and 5 are bounded by the DRAPED (untied) H5, whose model
        // position legitimately differs from the synthetic truth between the tie
        // wells — the boundary smear is real behavior, so those two get a looser
        // band; every mapped-bounded zone is held tight.
        let (por_tol, ntg_tol) = if z == 4 || z == 5 {
            (0.12, 0.16)
        } else {
            (0.06, 0.08)
        };
        assert!(
            (por[z].mean - t.por_mean).abs() < por_tol,
            "Z{z}: upscaled PORO mean {} vs target {}",
            por[z].mean,
            t.por_mean
        );
        assert!(
            (ntg[z].mean - t.ntg_mean).abs() < ntg_tol,
            "Z{z}: upscaled NTG mean {} vs target {}",
            ntg[z].mean,
            t.ntg_mean
        );
        zones_checked += 1;
    }
    assert!(
        zones_checked >= 8,
        "most zones must be checkable, got {zones_checked}"
    );

    // The sand/shale CONTRAST survives upscaling: the sand zones' upscaled PORO
    // means sit far above the shale zones'.
    let sand_mean = por[2].mean; // Z2 sand
    let shale_mean = por[3].mean; // Z3 shale
    assert!(
        sand_mean - shale_mean > 0.12,
        "sand/shale contrast lost: {sand_mean} vs {shale_mean}"
    );
}

#[test]
fn per_zone_in_place_hand_checked_vs_reference_and_conserves() {
    let m = build();
    let zoned = m.in_place_by_zone().unwrap();
    assert_eq!(zoned.zones.len(), 10);
    let by = |name: &str| {
        zoned
            .zones
            .iter()
            .find(|z| z.zone == name)
            .map(|z| &z.in_place)
            .unwrap()
    };
    let krange = |name: &str| m.zones().get(name).unwrap().k_range.clone();

    // Contactless zones: in_place GRV == the independent bulk reference, HCPV == 0.
    for z in [0usize, 1, 4, 5, 6, 8, 9] {
        let name = format!("Z{z}");
        let (bulk, _, _) = reference(&m, krange(&name), None, None);
        let ip = by(&name);
        assert!(
            approx(ip.grv_m3, bulk, 1e-9),
            "{name} bulk {} vs reference {bulk}",
            ip.grv_m3
        );
        assert_eq!(ip.hcpv_m3, 0.0, "contactless {name} has zero hydrocarbon");
    }

    // Zone 2 all-oil (OWC below the zone): GRV == full zone bulk, HCPV == GRV·HC.
    let (z2_bulk, _, _) = reference(&m, krange("Z2"), None, None);
    let z2 = by("Z2");
    assert!(approx(z2.grv_m3, z2_bulk, 1e-9));
    assert!(approx(z2.hcpv_m3, z2_bulk * HC, 1e-9));

    // Zone 7 genuine two-contact split — cross-checked against the geometry reference.
    let mean_h7 = mean_of(&build_framework(&INTERVALS, SEED)[7]);
    let (_, o_gas, o_oil) = reference(&m, krange("Z7"), Some(mean_h7 + 3.0), Some(mean_h7 + 5.0));
    let z7 = by("Z7");
    let (gas, oil) = (z7.gas.expect("gas"), z7.oil.expect("oil"));
    assert!(
        gas.grv_m3 > 0.0 && oil.grv_m3 > 0.0,
        "genuine gas+oil split"
    );
    assert!(
        approx(gas.grv_m3, o_gas, 1e-9),
        "gas {} vs reference {o_gas}",
        gas.grv_m3
    );
    assert!(
        approx(oil.grv_m3, o_oil, 1e-9),
        "oil {} vs reference {o_oil}",
        oil.grv_m3
    );
    assert!(approx(z7.grv_m3, gas.grv_m3 + oil.grv_m3, 1e-9));

    // Rollup conservation: total == sum of the per-zone GRV / HCPV.
    let sum_grv: f64 = zoned.zones.iter().map(|z| z.in_place.grv_m3).sum();
    let sum_hcpv: f64 = zoned.zones.iter().map(|z| z.in_place.hcpv_m3).sum();
    assert!(approx(zoned.total.grv_m3, sum_grv, 1e-12));
    assert!(approx(zoned.total.hcpv_m3, sum_hcpv, 1e-12));

    // Partition conservation: the per-zone bulk sums to the whole-grid bulk.
    let sum_bulk: f64 = (0..10)
        .map(|z| reference(&m, krange(&format!("Z{z}")), None, None).0)
        .sum();
    assert!(approx(sum_bulk, m.bulk_volume(), 1e-9));

    // Per-zone stats: constant priors → mean/min/max == the prior over active cells.
    for s in m.zone_stats(PORO).unwrap() {
        assert!(s.count > 0);
        assert!(
            approx(s.mean, PHI, 1e-12) && approx(s.min, PHI, 1e-12) && approx(s.max, PHI, 1e-12)
        );
    }
}

#[test]
fn spilled_in_place_by_zone_matches_in_core_within_tolerance() {
    // The v2 spilled by-zone gap closed (item 3): per-zone volumetrics on a spilled
    // (out-of-core) model streams each zone's contiguous k-band through the SAME
    // unified core as in-core. Assert per-zone GRV/HCPV parity within the f32
    // tolerance (R4), cell counts integer-identical (topology unchanged), and the
    // rollup still conserves.
    let core = build_with_budget(MemoryBudget::unlimited());
    let spill = build_with_budget(MemoryBudget::bytes(1024));
    assert!(
        !core.is_spilled() && spill.is_spilled(),
        "budgets drive the mode"
    );

    let zc = core.in_place_by_zone().unwrap();
    let zs = spill.in_place_by_zone().unwrap();
    assert_eq!(zc.zones.len(), zs.zones.len());

    let mut worst = 0.0f64;
    for (c, s) in zc.zones.iter().zip(&zs.zones) {
        assert_eq!(c.zone, s.zone, "zone order identical");
        // NB cell COUNT is not a stable f32 parity metric for a truncation-heavy
        // (FollowTop / pinch-out) zone: cells clamped to exactly zero thickness
        // in-core pick up tiny f32 noise and flip to marginally-positive volume. Their
        // volume is ~0, so the volume-weighted GRV/HCPV parity (the R4 contract) stays
        // tight — that is the meaningful assertion, not the integer count.
        let g = rel(s.in_place.grv_m3, c.in_place.grv_m3);
        let h = rel(s.in_place.hcpv_m3, c.in_place.hcpv_m3);
        worst = worst.max(g).max(h);
        assert!(g <= 1e-5, "{}: GRV parity {g:.2e}", c.zone);
        assert!(h <= 1e-5, "{}: HCPV parity {h:.2e}", c.zone);
        // Two-contact zone (Z7): the gas/oil split streams too.
        if let (Some(gc), Some(gs)) = (c.in_place.gas, s.in_place.gas) {
            assert!(rel(gs.hcpv_m3, gc.hcpv_m3) <= 1e-5, "{}: gas HCPV", c.zone);
        }
    }
    eprintln!("spilled↔in-core by-zone worst GRV/HCPV relative error: {worst:.2e}");

    // Rollup conservation holds in the spilled mode.
    let sum_grv: f64 = zs.zones.iter().map(|z| z.in_place.grv_m3).sum();
    let sum_hcpv: f64 = zs.zones.iter().map(|z| z.in_place.hcpv_m3).sum();
    assert!(approx(zs.total.grv_m3, sum_grv, 1e-12));
    assert!(approx(zs.total.hcpv_m3, sum_hcpv, 1e-12));
}

/// Relative error (|a-b|/|b|, guarding b==0) for the tolerance-parity asserts.
fn rel(a: f64, b: f64) -> f64 {
    if b == 0.0 {
        a.abs()
    } else {
        (a - b).abs() / b.abs()
    }
}

#[test]
fn procedural_stack_is_deterministic() {
    let a = build();
    let b = build();
    assert_eq!(a.bulk_volume(), b.bulk_volume());
    let (ia, ib) = (a.in_place_by_zone().unwrap(), b.in_place_by_zone().unwrap());
    assert_eq!(ia.total.grv_m3, ib.total.grv_m3);
    assert_eq!(ia.total.hcpv_m3, ib.total.hcpv_m3);
    for (za, zb) in ia.zones.iter().zip(&ib.zones) {
        assert_eq!(za.in_place.grv_m3, zb.in_place.grv_m3);
    }
}

#[test]
fn zone_aware_bundles_render_real_multi_zone() {
    let m = build();
    let nk = m.provenance().nk;

    let vb = m.volume_bundle(PORO).unwrap();
    assert_eq!(vb.zone_names.len(), 10);
    assert_eq!(vb.zone_names[9], "Z9");
    assert!(vb.zone_ids.contains(&9));

    let mb = m.map_bundle(&MapSpec::new().property(PORO)).unwrap();
    assert_eq!(mb.zone_averages.len(), 10);
    assert!(mb.zone_averages.iter().any(|l| l.name.contains("Z7")));

    let sec = m
        .intersection_bundle(
            &SectionSpec::Polyline(vec![
                [ORIGIN_X + 60.0, ORIGIN_Y + 60.0],
                [ORIGIN_X + 940.0, ORIGIN_Y + 60.0],
            ]),
            Some(PORO),
        )
        .unwrap();
    assert!(!sec.columns.is_empty());
    for col in &sec.columns {
        assert_eq!(col.layer_tops.len(), nk);
        for k in 0..nk {
            if !col.layer_tops[k].is_nan() {
                assert!(col.layer_bases[k] >= col.layer_tops[k]);
                assert!(col.values[k].is_finite());
            }
        }
    }

    // SCHEMA_VERSION 4: interior-horizon traces. An 11-horizon stack has 9 interior
    // horizons (H1..H9), each a depth polyline parallel to columns and ordered
    // top→down; H0/H10 are the structural top/base and are NOT repeated here.
    assert_eq!(sec.horizon_traces.len(), 9, "H1..H9 interior traces");
    for (t, expect) in sec.horizon_traces.iter().zip(1..=9) {
        assert_eq!(t.name, format!("H{expect}"));
        assert_eq!(
            t.depths.len(),
            sec.columns.len(),
            "trace runs parallel to columns"
        );
        assert!(t.depths.iter().all(|d| d.is_finite()));
    }
    // Interior horizons are ordered top→down at every column, bracketed by the
    // structural top/base of that column.
    for (c, col) in sec.columns.iter().enumerate() {
        let top = col
            .layer_tops
            .iter()
            .copied()
            .find(|d| !d.is_nan())
            .unwrap();
        let base = col
            .layer_bases
            .iter()
            .rev()
            .copied()
            .find(|d| !d.is_nan())
            .unwrap();
        let mut prev = top;
        for t in &sec.horizon_traces {
            assert!(t.depths[c] >= prev - 1e-6, "interior horizon ordered");
            prev = t.depths[c];
        }
        assert!(
            base >= prev - 1e-6,
            "base below the deepest interior horizon"
        );
    }
}

// --- cell-collapse on a controlled procedural thinning wedge ---

const WN: usize = 11; // 10×10 cells

/// A procedural 2-zone wedge: a thin PLATEAU (whole cells uniformly thin → the
/// max-pillar rule collapses them) ramping to a thick edge, over a flat slab whose
/// layers always exceed the threshold.
fn wedge_stack() -> HorizonStack {
    let flat = |d: f64| GriddedDepth {
        ncol: WN,
        nrow: WN,
        depth_m: vec![d; WN * WN],
        is_control: vec![true; WN * WN],
    };
    let mut mid = vec![0.0; WN * WN];
    for r in 0..WN {
        for c in 0..WN {
            let thickness = if c <= 3 {
                2.0
            } else {
                2.0 + 40.0 * ((c - 3) as f64 / (WN as f64 - 1.0 - 3.0))
            };
            mid[r * WN + c] = 5000.0 + thickness;
        }
    }
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "WT".into(),
                source: HorizonSource::Mapped(flat(5000.0)),
            },
            StackHorizon {
                name: "WM".into(),
                source: HorizonSource::Mapped(GriddedDepth {
                    ncol: WN,
                    nrow: WN,
                    depth_m: mid,
                    is_control: vec![true; WN * WN],
                }),
            },
            StackHorizon {
                name: "WB".into(),
                source: HorizonSource::Mapped(flat(5100.0)),
            },
        ],
        zone_layers: vec![
            StackZone {
                name: "WEDGE".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 10,
                contacts: vec![],
            },
            StackZone {
                name: "SLAB".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 8,
                contacts: vec![],
            },
        ],
    }
}

#[test]
fn collapse_conserves_volume_marks_cells_and_spares_thick_zone() {
    let no = StaticModelBuilder::from_horizon_stack(wedge_stack(), opts())
        .unwrap()
        .build()
        .unwrap();
    let yes = StaticModelBuilder::from_horizon_stack(wedge_stack(), opts())
        .unwrap()
        .with_collapse_below_m(0.5)
        .build()
        .unwrap();

    // (a) volume conserved to FP tolerance before/after collapse.
    assert!(
        approx(no.bulk_volume(), yes.bulk_volume(), 1e-12),
        "collapse conserves volume"
    );

    // Collapse fired, only in the thin WEDGE zone (no truncation involved).
    let stack = yes.provenance().stack.as_ref().unwrap();
    assert_eq!(stack.zones[0].name, "WEDGE");
    assert_eq!(stack.zones[0].truncated_cells, 0);
    assert!(yes.provenance().warnings.iter().any(
        |w| matches!(w, petekstatic::model::BuildWarning::CellsCollapsed { cells } if *cells > 0)
    ));

    // (c) the thick SLAB zone is untouched.
    let slab = yes.zones().get("SLAB").unwrap();
    let dims = yes.grid().dims();
    for k in slab.k_range.clone() {
        for j in 0..dims.nj {
            for i in 0..dims.ni {
                assert!(yes.grid().cell(Ijk::new(i, j, k)).dz() > 0.5);
            }
        }
    }

    // (d) property population unaffected on surviving cells — no NaN.
    let poro = yes.property(PORO).unwrap().values.as_slice();
    let ntg = yes.property(NTG).unwrap().values.as_slice();
    let sw = yes.property(SW).unwrap().values.as_slice();
    for c in dims.iter() {
        if yes.grid().cell(c).dz() <= 1e-9 {
            continue;
        }
        let idx = (c.k * dims.nj + c.j) * dims.ni + c.i;
        assert!(poro[idx].is_finite() && ntg[idx].is_finite() && sw[idx].is_finite());
    }

    // (e) collapsed cells are NaN-marked in the section bundle.
    let sec = yes
        .intersection_bundle(
            &SectionSpec::Polyline(vec![[15.0, 15.0], [285.0, 15.0]]),
            Some(PORO),
        )
        .unwrap();
    let mut saw_nan = false;
    for col in &sec.columns {
        for k in 0..col.layer_tops.len() {
            if col.layer_tops[k].is_nan() {
                saw_nan = true;
                assert!(col.layer_bases[k].is_nan() && col.values[k].is_nan());
            }
        }
    }
    assert!(saw_nan, "collapsed cells must be NaN-marked in the section");
}

// --- per-horizon well ties (P8 `task_petekstatic_multizone_2`) ---

#[test]
fn well_ties_condition_mapped_horizons_and_record_residuals() {
    let asset = synthetic_asset(N_WELLS, SEED, MAX_RESIDUAL_M);
    // Every well ties every horizon at its measured (residual-bearing) top.
    let ties: Vec<WellTie> = asset
        .wells
        .iter()
        .enumerate()
        .map(|(w, well)| {
            let (x, y) = well.local_xy();
            let mut t = WellTie::new(well.id.clone(), x, y, well.ip, well.jp);
            for h in 0..11 {
                t = t.with_top(format!("H{h}"), asset.tops[w][h]);
            }
            t
        })
        .collect();
    let m = StaticModelBuilder::from_horizon_stack(asset.stack.clone(), opts())
        .unwrap()
        .with_min_thickness_m(0.0)
        .with_well_ties(ties)
        .build()
        .unwrap();

    // (a) provenance carries one record per well, per-horizon residuals; the pre-tie
    //     residual on a MAPPED horizon equals the INJECTED residual (untied mapped
    //     surface == the clean synthetic surface).
    let recs = &m.provenance().well_ties;
    assert_eq!(recs.len(), asset.wells.len());
    for (w, rec) in recs.iter().enumerate() {
        assert_eq!(rec.id, asset.wells[w].id);
        assert_eq!(rec.residuals.len(), 11);
        let mut saw_nonzero = false;
        for h in 0..11 {
            if h == 5 {
                continue; // tops-only: recorded against the drape, ~0
            }
            let r = rec
                .residuals
                .iter()
                .find(|r| r.horizon == format!("H{h}"))
                .unwrap();
            assert!(
                approx(r.residual_m, asset.residuals[w][h], 1e-3),
                "well {w} H{h}: recorded {} vs injected {}",
                r.residual_m,
                asset.residuals[w][h]
            );
            if r.residual_m.abs() > 1.0 {
                saw_nonzero = true;
            }
        }
        assert!(
            saw_nonzero,
            "the recorded residuals must be genuinely non-zero"
        );
    }

    // (b) after tying, the surfaces move to honour the wells. The TOP horizon (H0)
    //     is never pushed by the order-repair, so it lands exactly on the measured
    //     top; an interior mapped horizon is honoured too, except where the
    //     order-repair must pull it DOWN to stay below the horizon above (independent
    //     noisy tops can cross) — so it is only ever at/deeper than the measured top,
    //     never shallower. That is a genuine tie: the untied build sits at the clean
    //     surface, off by the injected residual.
    for (w, well) in asset.wells.iter().enumerate() {
        let node = well.jp * N_NODES + well.ip;
        let h0 = m.framework().horizons[0].surface.depth_m[node];
        assert!(
            (h0 - asset.tops[w][0]).abs() < 1e-3,
            "well {w} H0: tied top {h0} not at measured {}",
            asset.tops[w][0]
        );
        let mut honoured = 0;
        for h in 0..11 {
            if h == 5 {
                continue;
            }
            let got = m.framework().horizons[h].surface.depth_m[node];
            assert!(
                got >= asset.tops[w][h] - 1e-3,
                "well {w} H{h}: tied surface {got} shallower than measured {}",
                asset.tops[w][h]
            );
            if (got - asset.tops[w][h]).abs() < 1e-3 {
                honoured += 1;
            }
        }
        assert!(
            honoured >= 6,
            "well {w}: too few horizons honoured ({honoured})"
        );
    }

    // (c) the map bundle surfaces the ties per well (wells[].ties + summary residual).
    let mb = m.map_bundle(&MapSpec::new()).unwrap();
    assert_eq!(mb.wells.len(), asset.wells.len());
    assert!(mb
        .wells
        .iter()
        .all(|w| w.ties.len() == 11 && w.tie_residual_m.is_some()));
    assert!(mb
        .wells
        .iter()
        .any(|w| w.ties.iter().any(|t| t.horizon == "H2")));

    // (d) an unknown horizon name in a tie is a build-time error.
    let bad = vec![WellTie::new("99/x-1", 0.0, 0.0, 1, 1).with_top("NOPE", 2000.0)];
    assert!(
        StaticModelBuilder::from_horizon_stack(asset.stack.clone(), opts())
            .unwrap()
            .with_well_ties(bad)
            .build()
            .is_err()
    );
}

// --- per-zone property population (P8 `task_petekstatic_multizone_2`) ---

#[test]
fn per_zone_priors_set_distinct_zone_levels_and_conserve_volume() {
    // Two-zone wedge; give WEDGE a sand prior and SLAB a shale prior. Each zone's
    // PORO/NTG level is its own, and the geometry (hence volume) is untouched — a
    // per-zone-population conservation check.
    let sand = ConstantPriors {
        porosity: 0.28,
        net_to_gross: 0.85,
        water_saturation: 0.20,
    };
    let shale = ConstantPriors {
        porosity: 0.08,
        net_to_gross: 0.15,
        water_saturation: 0.60,
    };
    let plain = StaticModelBuilder::from_horizon_stack(wedge_stack(), opts())
        .unwrap()
        .build()
        .unwrap();
    let m = StaticModelBuilder::from_horizon_stack(wedge_stack(), opts())
        .unwrap()
        .with_zone_priors("WEDGE", sand)
        .with_zone_priors("SLAB", shale)
        .build()
        .unwrap();

    // (a) geometry / volume conserved (population never moves a corner).
    assert!(approx(m.bulk_volume(), plain.bulk_volume(), 1e-12));

    // (b) each zone carries its OWN constant level over its active cells.
    let por = m.zone_stats(PORO).unwrap();
    let ntg = m.zone_stats(NTG).unwrap();
    let by = |v: &[petekstatic::model::ZoneStat], name: &str| {
        v.iter().find(|s| s.zone == name).unwrap().clone()
    };
    let (wp, sp) = (by(&por, "WEDGE"), by(&por, "SLAB"));
    assert!(approx(wp.mean, sand.porosity, 1e-12) && approx(wp.min, sand.porosity, 1e-12));
    assert!(approx(sp.mean, shale.porosity, 1e-12) && approx(sp.max, shale.porosity, 1e-12));
    assert!(approx(by(&ntg, "WEDGE").mean, sand.net_to_gross, 1e-12));
    assert!(approx(by(&ntg, "SLAB").mean, shale.net_to_gross, 1e-12));

    // (c) an unknown zone name is a build-time error.
    assert!(
        StaticModelBuilder::from_horizon_stack(wedge_stack(), opts())
            .unwrap()
            .with_zone_priors("NOPE", sand)
            .build()
            .is_err()
    );
}

#[test]
fn per_zone_property_populates_only_its_zone_over_a_zone_prior_baseline() {
    // A per-zone SGS pipeline on WEDGE simulates that zone from its own logs while
    // SLAB keeps its (distinct) constant prior — per-zone variogram/log-conditioning.
    let sand = ConstantPriors {
        porosity: 0.28,
        net_to_gross: 0.85,
        water_saturation: 0.20,
    };
    let shale = ConstantPriors {
        porosity: 0.08,
        net_to_gross: 0.15,
        water_saturation: 0.60,
    };
    // Two wells spanning the WEDGE depth band (WT=5000 down); local column lattice
    // dx = side/10 = 100, so column (0,0)/(9,9) centroids are at 50/950.
    let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 400.0).unwrap();
    let lo = WellLog::new(50.0, 50.0, vec![(5000.5, 0.24), (5001.0, 0.26)]);
    let hi = WellLog::new(950.0, 950.0, vec![(5000.5, 0.32), (5001.0, 0.30)]);
    let m = StaticModelBuilder::from_horizon_stack(wedge_stack(), opts())
        .unwrap()
        .with_zone_priors("WEDGE", sand)
        .with_zone_priors("SLAB", shale)
        .with_zone_property(
            "WEDGE",
            PropertyPipeline::new(PORO)
                .upscale(vec![lo, hi], UpscaleMethod::Arithmetic)
                // The WEDGE zone truncates, so deeper simulated layers legitimately
                // carry no conditioning data — opt into the mean-fill rather than the
                // default hard error (item 4).
                .propagate(petekstatic::model::Gaussian::new(vgm, 7).allow_mean_fill()),
        )
        .build()
        .unwrap();

    let dims = m.grid().dims();
    let poro = m.property(PORO).unwrap().values.as_slice();
    let wedge = m.zones().get("WEDGE").unwrap().k_range.clone();
    let slab = m.zones().get("SLAB").unwrap().k_range.clone();

    // SLAB untouched by the WEDGE pipeline: still exactly the shale prior.
    for k in slab {
        for j in 0..dims.nj {
            for i in 0..dims.ni {
                if m.grid().cell(Ijk::new(i, j, k)).dz() <= 1e-9 {
                    continue;
                }
                let v = poro[(k * dims.nj + j) * dims.ni + i];
                assert!(approx(v, shale.porosity, 1e-12), "SLAB cell moved: {v}");
            }
        }
    }
    // WEDGE simulated: finite everywhere, honours the wells, and genuinely varies
    // (not the flat sand baseline) — the SGS field is present.
    let mut distinct = false;
    for k in wedge {
        for j in 0..dims.nj {
            for i in 0..dims.ni {
                if m.grid().cell(Ijk::new(i, j, k)).dz() <= 1e-9 {
                    continue;
                }
                let v = poro[(k * dims.nj + j) * dims.ni + i];
                assert!(v.is_finite());
                if (v - sand.porosity).abs() > 1e-6 {
                    distinct = true;
                }
            }
        }
    }
    assert!(
        distinct,
        "WEDGE was not simulated (stayed the flat baseline)"
    );
    // Conservation: per-zone in-place rollup still balances.
    let zoned = m.in_place_by_zone().unwrap();
    let sum_grv: f64 = zoned.zones.iter().map(|z| z.in_place.grv_m3).sum();
    assert!(approx(zoned.total.grv_m3, sum_grv, 1e-12));
}

// SGS kernel friction note (finding, for the coordinator): `petektools::geostat::sgs`
// is CONDITIONAL-only — it requires ≥ 1 conditioning datum and derives the
// normal-score transform from the supplied data values (`NormalScore::fit`). There
// is no "unconditional simulation with target mean/variance" entry point, so an
// unconditional field must be faked by anchoring on seeded pseudo-points drawn from
// the target distribution (done above). A first-class `sgs_unconditional(lattice,
// mean, variance, variogram, seed)` would remove this dance for synthetic-fixture
// use; noted for petekTools.
