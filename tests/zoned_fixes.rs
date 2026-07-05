//! Rust-level acceptance for `task_petekstatic_zoned_fixes` — the three zoned-path
//! correctness fixes from the composer's planted-truth validation (2026-07-04):
//! (1) collocated-cokriging trend on WORLD-coordinate data (was a silent no-op),
//! (2) the zoned map outline (was the unit square vs a world frame), and
//! (3) per-zone NTG upscale compression (boundary cells averaging across zones).
//! All fixtures are synthetic (no dataset content) at a fictional study area.

use petekstatic::gridder::{Conformity, SolveOpts};
use petekstatic::model::{
    BuildOpts, ConstantPriors, Gaussian, HorizonSource, HorizonStack, MapSpec, PropertyPipeline,
    StackHorizon, StackZone, StaticModelBuilder, TrendSurface, UpscaleMethod, WellLog,
};
use petekstatic::volumetrics::{NTG, PORO};
use petekstatic::wireframe::{
    Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
};
use petektools::{Variogram, VariogramModel};

// A fictional UTM-magnitude study window (NOT any real field's coordinates).
const WX0: f64 = 431_000.0;
const WY0: f64 = 6_521_000.0;

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma).powi(2);
        sbb += (y - mb).powi(2);
    }
    sab / (saa.sqrt() * sbb.sqrt())
}

// ---- Finding 1: collocated cokriging on WORLD-coordinate data recovers the planted ρ ----

/// A flat world-georeferenced top surface on an `(n×n)`-node lattice at world origin.
fn world_wireframe(n: usize) -> Wireframe {
    let side = 100.0 * (n as f64 - 1.0); // dx = 100 m
    Wireframe {
        boundary: Boundary {
            ring: vec![
                [WX0, WY0],
                [WX0 + side, WY0],
                [WX0 + side, WY0 + side],
                [WX0, WY0 + side],
                [WX0, WY0],
            ],
            hardness: Hardness::Hard,
        },
        horizons: std::sync::Arc::new(vec![Horizon {
            name: "TopRes".into(),
            role: HorizonRole::Top,
            surface: GriddedDepth {
                ncol: n,
                nrow: n,
                depth_m: vec![2000.0; n * n],
                is_control: vec![true; n * n],
            },
        }]),
        contacts: vec![Contact {
            kind: ContactKind::Owc,
            depth_m: 2025.0,
            hardness: Hardness::Hard,
        }],
    }
}

fn world_opts(nk: usize) -> BuildOpts {
    BuildOpts {
        area_m2: (100.0 * 20.0_f64).powi(2), // side 2000 -> dx = 100 over 20 columns
        gross_height_m: 40.0,
        nk,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.22,
            net_to_gross: 0.7,
            water_saturation: 0.3,
        },
    }
}

/// A world-georeferenced trend ramp increasing in i+j over the model's 20×20 columns,
/// its georef == the grid's world column lattice (origin (0,0)-centroid + dx spacing).
fn world_trend(ni: usize, nj: usize, dx: f64) -> TrendSurface {
    let values: Vec<f64> = (0..ni * nj)
        .map(|k| (k % ni + k / ni) as f64) // i + j ramp
        .collect();
    TrendSurface::new(ni, nj, values)
        .unwrap()
        .with_georef(WX0 + dx / 2.0, WY0 + dx / 2.0, dx, dx)
}

fn recovered_rho(corr: f64) -> f64 {
    let (ni, nj, dx) = (20usize, 20usize, 100.0);
    // Two sparse conditioning wells at opposite corners so the trend, not the wells,
    // governs the areal pattern. Well logs are positioned in the model's LOCAL grid
    // frame (column centroids 50..1950) — it is the *trend* that is world-georeferenced
    // (finding 1: the georef maps the local grid to world so the world trend overlaps).
    let wells = vec![
        WellLog::new(50.0, 50.0, vec![(2005.0, 0.20)]),
        WellLog::new(1950.0, 1950.0, vec![(2005.0, 0.24)]),
    ];
    let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 1200.0).unwrap();
    let pipe = PropertyPipeline::new(PORO)
        .upscale(wells, UpscaleMethod::Arithmetic)
        // Some simulated layers carry no well data here; opt into the mean-fill
        // rather than the default hard error (item 4) — the collocated recovery is
        // measured on the conditioned layers.
        .propagate(
            Gaussian::new(vgm, 42)
                .with_trend(world_trend(ni, nj, dx), corr)
                .allow_mean_fill(),
        );
    let m = StaticModelBuilder::from_wireframe(&world_wireframe(21), world_opts(3))
        .unwrap()
        .with_georef(WX0 + dx / 2.0, WY0 + dx / 2.0, dx, dx)
        .with_property(pipe)
        .build()
        .unwrap();
    // Correlate the k=0 slice with the trend ramp over the columns.
    let cube = &m.property(PORO).unwrap().values;
    let dims = m.grid().dims();
    let mut field = Vec::new();
    let mut trend = Vec::new();
    for j in 0..dims.nj {
        for i in 0..dims.ni {
            field.push(cube[j * dims.ni + i]);
            trend.push((i + j) as f64);
        }
    }
    pearson(&field, &trend)
}

#[test]
fn collocated_trend_recovers_planted_correlation_on_world_data() {
    // The BUG (finding 1): on world-coordinate data the trend resampled onto the grid's
    // LOCAL lattice was all-NaN, so every node fell back to plain SGS and the recovered
    // field-vs-trend correlation was INDEPENDENT of the requested ρ (~0.11 at any ρ).
    // The fix samples the world trend at each column's world position (via the model
    // georef), so the planted ρ is recovered.
    let r0 = recovered_rho(0.0);
    let r6 = recovered_rho(0.6);
    // ρ=0 must NOT track the trend (plain SGS with two corner wells still leaves a weak
    // diagonal from the conditioning, so bound it loosely well below the planted 0.6).
    assert!(
        r0.abs() < 0.45,
        "ρ=0 should not track the trend, got {r0:.3}"
    );
    // ρ=0.6: recovered field-vs-trend Pearson within a justified tolerance of 0.6.
    // Tolerance derivation: Markov-1 collocated cokriging reproduces the secondary
    // correlation only approximately, and the Pearson estimate over a 20×20 = 400-node
    // field carries SGS sampling noise ~1/sqrt(400) ≈ 0.05; a ±0.20 band absorbs both
    // the Markov-1 bias and the finite-size noise while still proving ρ is recovered
    // (the bug produced ~0.11 regardless — far outside this band).
    assert!(
        (r6 - 0.6).abs() < 0.20,
        "ρ=0.6 planted, recovered field-vs-trend correlation {r6:.3} (ρ=0 gave {r0:.3})"
    );
    // And the effect is monotone: more collocation => stronger trend correlation.
    assert!(r6 > r0 + 0.15, "ρ=0.6 ({r6:.3}) must beat ρ=0 ({r0:.3})");
}

// ---- Finding 2: the zoned (from_horizon_stack) map outline is world, not the unit square ----

fn flat_gd(n: usize, depth: f64) -> GriddedDepth {
    GriddedDepth {
        ncol: n,
        nrow: n,
        depth_m: vec![depth; n * n],
        is_control: vec![true; n * n],
    }
}

fn two_zone_stack(n: usize) -> HorizonStack {
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "H0".into(),
                source: HorizonSource::Mapped(flat_gd(n, 2000.0)),
            },
            StackHorizon {
                name: "H1".into(),
                source: HorizonSource::Mapped(flat_gd(n, 2030.0)),
            },
            StackHorizon {
                name: "H2".into(),
                source: HorizonSource::Mapped(flat_gd(n, 2060.0)),
            },
        ],
        zone_layers: vec![
            StackZone {
                name: "Z0".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 4,
                contacts: vec![],
            },
            StackZone {
                name: "Z1".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 4,
                contacts: vec![],
            },
        ],
    }
}

#[test]
fn zoned_map_outline_is_world_not_unit_square() {
    // The BUG (finding 2): from_horizon_stack emitted map_bundle.outline as the unit
    // square [0,1]×[0,1] while the frame + wells are world coords, collapsing the
    // viewer's content extent. The fix carries a world outline (georef-derived extent,
    // or an explicit with_boundary ring).
    let (n, dx) = (11usize, 100.0_f64); // 10×10 cells, side 1000
    let opts = BuildOpts {
        area_m2: (dx * 10.0).powi(2),
        gross_height_m: 0.0,
        nk: 0,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.2,
            net_to_gross: 0.7,
            water_saturation: 0.3,
        },
    };
    let m = StaticModelBuilder::from_horizon_stack(two_zone_stack(n), opts)
        .unwrap()
        .with_georef(WX0 + dx / 2.0, WY0 + dx / 2.0, dx, dx)
        .build()
        .unwrap();
    let b = m.map_bundle(&MapSpec::new()).unwrap();
    let f = &b.frame;
    // The frame is world (via the georef).
    assert!(f.origin_x > 400_000.0, "frame is world: {f:?}");

    // The outline is NOT the degenerate unit square, and every outline point lies
    // inside the frame's world cell-edge extent — so the viewer's content extent no
    // longer collapses (the outline extent ≈ the frame extent).
    let xmin = f.origin_x - f.spacing_x / 2.0;
    let xmax = f.origin_x + (f.ncol as f64 - 0.5) * f.spacing_x;
    let ymin = f.origin_y - f.spacing_y / 2.0;
    let ymax = f.origin_y + (f.nrow as f64 - 0.5) * f.spacing_y;
    assert!(!b.outline.is_empty());
    let (mut ox0, mut ox1) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut oy0, mut oy1) = (f64::INFINITY, f64::NEG_INFINITY);
    for ring in &b.outline {
        for p in ring {
            assert!(
                p[0] >= xmin - 1e-6
                    && p[0] <= xmax + 1e-6
                    && p[1] >= ymin - 1e-6
                    && p[1] <= ymax + 1e-6,
                "outline point {p:?} outside frame extent x[{xmin},{xmax}] y[{ymin},{ymax}]"
            );
            ox0 = ox0.min(p[0]);
            ox1 = ox1.max(p[0]);
            oy0 = oy0.min(p[1]);
            oy1 = oy1.max(p[1]);
        }
    }
    // Outline extent ≈ frame extent (world), not the collapsed unit square.
    assert!(
        (ox0 - xmin).abs() < 1e-6 && (ox1 - xmax).abs() < 1e-6,
        "x extent"
    );
    assert!(
        (oy0 - ymin).abs() < 1e-6 && (oy1 - ymax).abs() < 1e-6,
        "y extent"
    );
    assert!(ox1 - ox0 > 900.0, "outline spans ~1000 m, not ~1");

    // An explicit with_boundary ring is carried through verbatim (world shape fidelity).
    let ring = vec![
        [WX0 + 100.0, WY0 + 100.0],
        [WX0 + 900.0, WY0 + 200.0],
        [WX0 + 500.0, WY0 + 900.0],
        [WX0 + 100.0, WY0 + 100.0],
    ];
    let opts2 = BuildOpts {
        area_m2: (dx * 10.0).powi(2),
        gross_height_m: 0.0,
        nk: 0,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.2,
            net_to_gross: 0.7,
            water_saturation: 0.3,
        },
    };
    let m2 = StaticModelBuilder::from_horizon_stack(two_zone_stack(n), opts2)
        .unwrap()
        .with_georef(WX0 + dx / 2.0, WY0 + dx / 2.0, dx, dx)
        .with_boundary(ring.clone())
        .build()
        .unwrap();
    let b2 = m2.map_bundle(&MapSpec::new()).unwrap();
    assert_eq!(
        b2.outline,
        vec![ring],
        "explicit world ring carried through"
    );
}

// ---- Finding 3: per-zone NTG upscale no longer compresses off-centroid wells ----

/// A 2-zone stack with the internal boundary tilted in i (steep) so a well off the
/// column centroid sees a zone boundary depth that differs from the centroid.
fn tilted_stack(n: usize, slope: f64) -> HorizonStack {
    let top = vec![2000.0; n * n];
    let base = vec![2060.0; n * n];
    let mut mid = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            mid[j * n + i] = 2030.0 + slope * (i as f64 - (n as f64 - 1.0) / 2.0);
        }
    }
    let gd = |d: Vec<f64>| GriddedDepth {
        ncol: n,
        nrow: n,
        depth_m: d,
        is_control: vec![true; n * n],
    };
    HorizonStack {
        horizons: vec![
            StackHorizon {
                name: "H0".into(),
                source: HorizonSource::Mapped(gd(top)),
            },
            StackHorizon {
                name: "H1".into(),
                source: HorizonSource::Mapped(gd(mid)),
            },
            StackHorizon {
                name: "H2".into(),
                source: HorizonSource::Mapped(gd(base)),
            },
        ],
        zone_layers: vec![
            StackZone {
                name: "Z0".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 4,
                contacts: vec![],
            },
            StackZone {
                name: "Z1".into(),
                color: None,
                conformity: Conformity::Proportional,
                nk: 4,
                contacts: vec![],
            },
        ],
    }
}

/// A net-flag NTG well log at LOCAL grid `(wx, wy)`: samples across [2000, 2060] with
/// the per-zone planted proportions (Z0 = 0.45, Z1 = 0.80), the zone split taken at the
/// TRUE boundary depth at the well (bilinear on the tilted mid horizon, node index
/// `wx/dx`) — the physical truth a real LAS records at the well's own position.
fn net_log(wx: f64, wy: f64, dx: f64, slope: f64, n: usize) -> Vec<WellLog> {
    let fi = wx / dx; // fractional node index in i (local node 0 at x=0)
    let mid = 2030.0 + slope * (fi - (n as f64 - 1.0) / 2.0);
    let mut samples = Vec::new();
    let (mut n0, mut n1) = (0usize, 0usize);
    let mut d = 2000.5;
    while d < 2059.5 {
        let (cnt, frac) = if d >= mid {
            (&mut n1, 0.80)
        } else {
            (&mut n0, 0.45)
        };
        let flag =
            ((*cnt as f64 * frac).floor() != ((*cnt as f64 + 1.0) * frac).floor()) as i32 as f64;
        *cnt += 1;
        samples.push((d, flag));
        d += 0.25;
    }
    vec![WellLog::new(wx, wy, samples)]
}

fn zoned_ntg(wx: f64, wy: f64) -> (f64, f64) {
    // slope 5 m/node keeps the tilted mid boundary strictly inside (2000, 2060) while
    // still offsetting the boundary depth between a column centroid and its corner.
    let (n, dx, slope) = (9usize, 100.0_f64, 5.0_f64);
    let opts = BuildOpts {
        area_m2: (dx * (n as f64 - 1.0)).powi(2),
        gross_height_m: 0.0,
        nk: 0,
        conformity: Conformity::Proportional,
        solve_opts: SolveOpts::default(),
        priors: ConstantPriors {
            porosity: 0.2,
            net_to_gross: 0.5,
            water_saturation: 0.3,
        },
    };
    let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 400.0).unwrap();
    let m = StaticModelBuilder::from_horizon_stack(tilted_stack(n, slope), opts)
        .unwrap()
        .with_min_thickness_m(0.0)
        .with_georef(WX0 + dx / 2.0, WY0 + dx / 2.0, dx, dx)
        .with_zone_property(
            "Z0",
            PropertyPipeline::new(NTG)
                .upscale(net_log(wx, wy, dx, slope, n), UpscaleMethod::Arithmetic)
                .propagate(Gaussian::new(vgm, 1)),
        )
        .with_zone_property(
            "Z1",
            PropertyPipeline::new(NTG)
                .upscale(net_log(wx, wy, dx, slope, n), UpscaleMethod::Arithmetic)
                .propagate(Gaussian::new(vgm, 2)),
        )
        .build()
        .unwrap();
    let s = m.zone_stats(NTG).unwrap();
    (s[0].mean, s[1].mean)
}

#[test]
fn per_zone_ntg_recovers_planted_targets_for_offcentroid_wells() {
    // The BUG (finding 3): a well off its column centroid on a dipping zone boundary had
    // its near-boundary net-flag samples mis-binned across the zone boundary (the cell
    // interval was interpolated at the column CENTROID, not the well), so the per-zone
    // upscaled NTG compressed toward the mid-range. Reproduced pre-fix: planted Z0=0.45
    // read ~0.59 (a ~0.14 error, matching the composer's ~0.15). The fix bins each
    // sample against the cell interval interpolated AT THE WELL, recovering the targets.
    //
    // A well AT the column centroid was always fine (bilinear at centre == 4-corner
    // mean); the off-centroid well is the regression witness.
    // Local column 1 centroid is at x = 150; a corner-offset well at local x = 100.
    let (z0_off, z1_off) = zoned_ntg(100.0, 150.0);
    let (z0_ctr, z1_ctr) = zoned_ntg(150.0, 150.0);

    // Post-fix per-zone tolerance: 0.05 (tightened from the observed ~0.15 pre-fix
    // error). The residual is SGS reproduction noise + the net-flag quantization over a
    // ~60-sample zone column, NOT boundary cross-contamination.
    for (z0, z1, who) in [
        (z0_off, z1_off, "off-centroid"),
        (z0_ctr, z1_ctr, "centroid"),
    ] {
        assert!(
            (z0 - 0.45).abs() < 0.05,
            "{who}: Z0 NTG {z0:.3} vs planted 0.45"
        );
        assert!(
            (z1 - 0.80).abs() < 0.05,
            "{who}: Z1 NTG {z1:.3} vs planted 0.80"
        );
    }
}
