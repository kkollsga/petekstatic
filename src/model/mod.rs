//! `srs-model` — the top of the geomodel DAG: the [`StaticModel`] aggregate, the
//! deterministic [`StaticModelBuilder`], and the Monte-Carlo regeneration seam
//! ([`StaticModelTemplate`] + [`RealizationDraw`], graph
//! `decision_staticmodel_regen_seam`, ratified 2026-07-03).
//!
//! A `StaticModel` is a *populated* static reservoir model — framework + grid +
//! property cubes + zones + contacts + provenance — and, per the layer charter
//! (graph `decision_layer_charters`), it **owns its volumetrics output surface**:
//! [`StaticModel::in_place`] → [`InPlace`] (+ FVF via [`OilFvf`]/[`GasFvf`]),
//! with P-curves aggregated over realizations via [`PercentileSummary`].
//!
//! The relocated `RefiningModel` orchestration (petekSim `srs-core`,
//! `task_relocate_refine_orchestration`) lives in [`builder`]; petekSim's
//! `srs-core` now consumes this crate across the repo seam.

mod builder;
mod draw;
mod mc;
#[allow(clippy::module_inception)]
mod model;
mod pipeline;
mod population;
mod provenance;
mod spec;
mod template;
mod tornado;
mod trend;
mod view;
mod zones;

pub use builder::{
    BuildOpts, HorizonSource, HorizonStack, Pick, StackFrame, StackHorizon, StackZone,
    StaticModelBuilder, WellTie, WorldPoint,
};
pub use draw::{PerturbationField, RealizationDraw, StructuralPerturbation, ZoneDraw};
pub use mc::{
    aggregate_field, default_mc_workers, run_mc, Input, McInputs, McResult, McSettings,
    RealizedInputs,
};
// The out-of-core backing-storage mode lives in `srs-spill` (split out P10); its
// surface is re-exported here so the public path (`crate::model::MemoryBudget`,
// `crate::model::spill_grid`, …) is unchanged by the crate split.
pub use crate::spill::{
    decide_mode, live_set_bytes, physical_ram_bytes, spill_grid, spill_grid_to, BuildMode,
    MemoryBudget, SpillBacking, SpillNotice, DEFAULT_BUDGET_FRACTION,
};
// The four historical MC entries — deprecated thin wrappers over `run_mc`.
#[allow(deprecated)]
pub use mc::{
    run_structured_mc, run_structured_mc_parallel, run_structured_mc_parallel_spilled,
    run_structured_mc_spilled,
};
pub use model::{Georef, StaticModel, ZoneInPlace, ZoneStat, ZonedInPlace};
pub use pipeline::{
    Gaussian, McMode, PropertyPipeline, PropertyReport, UpscaleMethod, UpscaleQc, WellLog,
};
pub use population::PetroSample;
pub use provenance::{
    BuildWarning, HorizonTieResidual, InterfaceRepair, PopulationMode, Provenance, StackProvenance,
    WellTieRecord, ZoneProvenance,
};
pub use spec::{BuildSpec, TieMethod, TieSettings};
pub use template::StaticModelTemplate;
pub use tornado::{tornado, TornadoBar};
pub use trend::TrendSurface;
pub use view::{
    ContactMask, GridFrame, HorizonTrace, IntersectionBundle, MapBundle, MapSpec, ScalarLayer,
    SectionColumn, SectionContact, SectionSpec, SectionZone, ValueRange, VolumeBundle, WellMarker,
    WellTieResidual, SCHEMA_VERSION,
};
pub use zones::{Zone, ZoneTable};

// The volumetrics / P-curve output surface (the model owns volumes now).
// `OilFvf`/`GasFvf` are the FVF **seam scalars** (`crate::volumetrics::fvf`): the
// house-style "duplicate a small type at the seam" — petekSim's `srs-pvt` keeps its
// own copy for the dynamic/PVT work, and no PVT code crosses down into this layer
// (graph `decision_layer_charters`). FVF enters MC as a validated uncertain scalar,
// applied at the volumetrics conversion, never inside a `RealizationDraw`.
pub use crate::uncertainty::PercentileSummary;
pub use crate::volumetrics::{ConstantPriors, GasFvf, InPlace, OilFvf, ZoneVolumes};

// The structured-MC seam re-exports petekTools' sampler + correlation + digest
// types so a caller builds an [`McInputs`] and reads [`McResult::summary`]
// without a direct petekTools import (the driver owns the seam).
pub use petektools::sampling::{Correlation, ReservoirSummary, Sampler};

// Ratified amendment 1: template + model are Send (Sync not required) — the
// consumer shards realizations one-template-per-worker. Compile-time checked.
static_assertions::assert_impl_all!(StaticModelTemplate: Send);
static_assertions::assert_impl_all!(StaticModel: Send);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gridder::{Conformity, SolveOpts};
    use crate::wireframe::{
        Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
    };

    /// A wireframe with an `n×n` top surface at a constant `depth_m` and one OWC.
    fn flat_wireframe(n: usize, depth_m: f64, owc_m: f64) -> Wireframe {
        Wireframe {
            boundary: Boundary {
                ring: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
                hardness: Hardness::Interpolated,
            },
            horizons: std::sync::Arc::new(vec![Horizon {
                name: "top".into(),
                role: HorizonRole::Top,
                surface: GriddedDepth {
                    ncol: n,
                    nrow: n,
                    depth_m: vec![depth_m; n * n],
                    is_control: vec![true; n * n],
                },
            }]),
            contacts: vec![Contact {
                kind: ContactKind::Owc,
                depth_m: owc_m,
                hardness: Hardness::Hard,
            }],
        }
    }

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

    // --- per-property geostatistical pipeline (P5) ---

    #[test]
    fn builder_property_pipeline_populates_a_conditioned_cube() {
        use petektools::{Variogram, VariogramModel};
        // 11×11 flat top at 5000, 100 m² (side 10, dx=dy=1), gross 50 over nk=5
        // (dz=10). Columns centroid x = (i+0.5); layer k spans [5000+10k, 5010+10k].
        let wf = flat_wireframe(11, 5000.0, 5100.0);
        let low = WellLog::new(
            0.5,
            0.5,
            vec![
                (5005.0, 0.10),
                (5015.0, 0.12),
                (5025.0, 0.14),
                (5035.0, 0.16),
                (5045.0, 0.18),
            ],
        );
        let high = WellLog::new(
            9.5,
            9.5,
            vec![
                (5005.0, 0.26),
                (5015.0, 0.25),
                (5025.0, 0.24),
                (5035.0, 0.23),
                (5045.0, 0.22),
            ],
        );
        let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 5.0).unwrap();
        let pipe = PropertyPipeline::new("PHIE")
            .upscale(vec![low, high], UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(vgm, 42));

        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .with_property(pipe)
            .build()
            .unwrap();

        // The property landed, fully populated (no NaN), conditioned on the wells.
        let prop = m.property("PHIE").expect("PHIE cube present");
        assert_eq!(prop.values.len(), 10 * 10 * 5);
        assert!(prop.values.iter().all(|v| v.is_finite()));
        // A pipeline flips population to Logs and records a report.
        assert_eq!(m.provenance().population, PopulationMode::Logs);
        assert_eq!(m.provenance().property_reports.len(), 1);
        let report = &m.provenance().property_reports[0];
        assert_eq!(report.property, "PHIE");
        assert!(report.propagated);
        assert_eq!(report.upscale.conditioned_cells, 10); // 2 wells × 5 layers
    }

    // --- MC modes for the property pipeline (decision_mc_composition) ---

    fn phie_pipeline(seed: u64) -> PropertyPipeline {
        use petektools::{Variogram, VariogramModel};
        let low = WellLog::new(
            0.5,
            0.5,
            vec![
                (5005.0, 0.10),
                (5015.0, 0.12),
                (5025.0, 0.14),
                (5035.0, 0.16),
                (5045.0, 0.18),
            ],
        );
        let high = WellLog::new(
            9.5,
            9.5,
            vec![
                (5005.0, 0.26),
                (5015.0, 0.25),
                (5025.0, 0.24),
                (5035.0, 0.23),
                (5045.0, 0.22),
            ],
        );
        let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 5.0).unwrap();
        PropertyPipeline::new("PHIE")
            .upscale(vec![low, high], UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(vgm, seed))
    }

    fn base_draw(seed_index: u64) -> RealizationDraw {
        RealizationDraw::new(100.0, 50.0, 5100.0, 0.25, 0.8, 0.3, seed_index)
    }

    /// A `PORO` cube conditioned to the physical boundary: a non-net well at 0.0
    /// and a fully-porous well at 1.0, so SGS honours 0.0 and 1.0 cells exactly —
    /// the situation a real log-conditioned model produces (a non-net NTG=0 / an
    /// aquifer SW=1 boundary cell). Targets `PORO` so `in_place` per-cell
    /// H2-validates the cube (F9).
    fn poro_boundary_pipeline(seed: u64) -> PropertyPipeline {
        use petektools::{Variogram, VariogramModel};
        let zero = WellLog::new(
            0.5,
            0.5,
            vec![
                (5005.0, 0.0),
                (5015.0, 0.0),
                (5025.0, 0.0),
                (5035.0, 0.0),
                (5045.0, 0.0),
            ],
        );
        let one = WellLog::new(
            9.5,
            9.5,
            vec![
                (5005.0, 1.0),
                (5015.0, 1.0),
                (5025.0, 1.0),
                (5035.0, 1.0),
                (5045.0, 1.0),
            ],
        );
        let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 5.0).unwrap();
        PropertyPipeline::new("PORO")
            .upscale(vec![zero, one], UpscaleMethod::Arithmetic)
            .propagate(Gaussian::new(vgm, seed))
    }

    #[test]
    fn level_shift_saturates_a_fraction_cube_holding_zeros_and_ones() {
        // F9: a conditioned PORO cube legitimately holds 0.0 and 1.0 cells. A level
        // shift adds the drawn amount to every cell; boundary cells must SATURATE at
        // [0,1] (shift-then-clamp) rather than escape the range and trip the per-cell
        // H2 check — else property uncertainty on any log-conditioned model fails.
        let wf = flat_wireframe(11, 5000.0, 5100.0);
        let mut t = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property(poro_boundary_pipeline(42));
        // Draw 0, no shift: the raw conditioned pattern; it holds a ~0 and a ~1 cell.
        let m0 = t.realize(&base_draw(0)).unwrap();
        let pat = m0.property("PORO").unwrap().values.clone();
        assert!(
            pat.iter().any(|&v| v <= 1e-9),
            "cube has a non-net (0.0) cell"
        );
        assert!(
            pat.iter().any(|&v| v >= 1.0 - 1e-9),
            "cube has a fully-saturated (1.0) cell"
        );

        // A large positive shift: the 1.0 cells would go > 1 without the clamp.
        let shift = 0.1;
        let m1 = t
            .realize(&base_draw(1).with_property_shift("PORO", shift))
            .unwrap();
        let cube = m1.property("PORO").unwrap().values.clone();
        let mut saturated = 0usize;
        let (mut d_interior, mut n_interior) = (0.0f64, 0usize);
        for (v0, v1) in pat.iter().zip(&cube) {
            let expect = (v0 + shift).clamp(0.0, 1.0);
            assert!(
                (v1 - expect).abs() < 1e-12,
                "shift-then-clamp per cell: {v0} -> {v1}"
            );
            assert!((0.0..=1.0).contains(v1), "shifted cell out of [0,1]: {v1}");
            if v0 + shift > 1.0 {
                saturated += 1;
                assert!((v1 - 1.0).abs() < 1e-12, "boundary cell saturates at 1.0");
            } else {
                d_interior += v1 - v0;
                n_interior += 1;
            }
        }
        assert!(
            saturated > 0,
            "a 1.0 cell saturated under the positive shift"
        );
        // Over the interior (non-saturated) cells the field moves by ~the shift.
        assert!(
            (d_interior / n_interior as f64 - shift).abs() < 1e-12,
            "interior mean shifts by the drawn amount"
        );
        // in_place succeeds — the conditioned model is no longer rejected at draw #0.
        assert!(
            m1.in_place().is_ok(),
            "level-shifted conditioned cube is a valid model"
        );
    }

    #[test]
    fn level_shift_conditioned_cube_survives_1000_draws_but_inputs_stay_h2() {
        use crate::error::StaticError;
        let wf = flat_wireframe(11, 5000.0, 5100.0);
        let mut t = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property(poro_boundary_pipeline(42));
        // 1000 draws over a wide shift band (±0.5) that would push boundary cells
        // far outside [0,1] without the clamp — every realization stays valid.
        for i in 0..1000u64 {
            let shift = -0.5 + (i as f64) * (1.0 / 999.0);
            let m = t
                .realize(&base_draw(i).with_property_shift("PORO", shift))
                .expect("level-shift realization is valid");
            let cube = &m.property("PORO").unwrap().values;
            assert!(
                cube.iter().all(|v| (0.0..=1.0).contains(v)),
                "draw {i}: cube in [0,1]"
            );
            assert!(m.in_place().is_ok(), "draw {i}: valid in-place");
        }
        // H2 still rejects a garbage DRAWN INPUT — the clamp is a per-cell cube
        // application, NOT a licence for a garbage sampler. A porosity prior outside
        // [0,1] still errors.
        let bad = RealizationDraw::new(100.0, 50.0, 5100.0, 1.5, 0.8, 0.3, 0);
        assert!(
            matches!(t.realize(&bad), Err(StaticError::InvalidInput(_))),
            "garbage input prior still rejected (H2)"
        );
    }

    #[test]
    fn level_shift_keeps_the_pattern_and_moves_the_mean() {
        let wf = flat_wireframe(11, 5000.0, 5100.0);
        // Default mode is LevelShift.
        let mut t = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property(phie_pipeline(42));
        // First draw: no shift -> the raw once-propagated pattern.
        let m0 = t.realize(&base_draw(0)).unwrap();
        let c0 = m0.property("PHIE").unwrap().values.clone();
        // Second draw: +0.05 level shift, different seed_index (must NOT re-propagate).
        let m1 = t
            .realize(&base_draw(1).with_property_shift("PHIE", 0.05))
            .unwrap();
        let c1 = m1.property("PHIE").unwrap().values.clone();
        // Same pattern everywhere, shifted by exactly 0.05.
        for (a, b) in c0.iter().zip(&c1) {
            assert!(
                (b - a - 0.05).abs() < 1e-12,
                "not a pure level shift: {a} {b}"
            );
        }
        let mean0 = c0.iter().sum::<f64>() / c0.len() as f64;
        let mean1 = c1.iter().sum::<f64>() / c1.len() as f64;
        assert!(
            (mean1 - mean0 - 0.05).abs() < 1e-9,
            "mean not shifted by 0.05"
        );
    }

    #[test]
    fn resimulate_redraws_a_new_pattern_per_seed() {
        let wf = flat_wireframe(11, 5000.0, 5100.0);
        let mut t = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property_mode(phie_pipeline(42), McMode::Resimulate);
        let a = t
            .realize(&base_draw(1))
            .unwrap()
            .property("PHIE")
            .unwrap()
            .values
            .clone();
        let b = t
            .realize(&base_draw(2))
            .unwrap()
            .property("PHIE")
            .unwrap()
            .values
            .clone();
        assert_ne!(
            a, b,
            "resimulate must give a different pattern per seed_index"
        );
        // Both still honour the well at column (0,0,0) exactly (SGS conditioning);
        // cell (0,0,0) is row-major index 0.
        assert!((a[0] - 0.10).abs() < 1e-6 && (b[0] - 0.10).abs() < 1e-6);
    }

    #[test]
    fn mc_property_pipeline_is_bit_reproducible_both_modes() {
        let wf = flat_wireframe(11, 5000.0, 5100.0);
        // LevelShift: identical draws -> identical cubes.
        let mut tl = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property(phie_pipeline(42));
        let l0 = tl
            .realize(&base_draw(7))
            .unwrap()
            .property("PHIE")
            .unwrap()
            .values
            .clone();
        let l1 = tl
            .realize(&base_draw(7))
            .unwrap()
            .property("PHIE")
            .unwrap()
            .values
            .clone();
        assert_eq!(
            l0, l1,
            "LevelShift must be bit-reproducible for identical draws"
        );
        // Resimulate: same seed_index -> identical pattern.
        let mut tr = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property_mode(phie_pipeline(42), McMode::Resimulate);
        let r0 = tr
            .realize(&base_draw(5))
            .unwrap()
            .property("PHIE")
            .unwrap()
            .values
            .clone();
        let r1 = tr
            .realize(&base_draw(5))
            .unwrap()
            .property("PHIE")
            .unwrap()
            .values
            .clone();
        assert_eq!(
            r0, r1,
            "Resimulate must be bit-reproducible for the same seed_index"
        );
    }

    // --- deterministic builder (the relocated pipeline) ---

    #[test]
    fn flat_build_matches_box_volume() {
        // 11×11 flat top, contact below the base -> full column == area*height.
        let wf = flat_wireframe(11, 5000.0, 5100.0);
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();
        let expected_grv_m3 = 100.0 * 50.0;
        let ip = m.in_place().unwrap();
        assert!((ip.grv_m3 - expected_grv_m3).abs() / expected_grv_m3 < 1e-6);
        // The model owns volumes: OOIP off the model's own in-place + FVF input.
        let ooip = ip.ooip_sm3(OilFvf::new(1.25).unwrap());
        assert!(ooip > 0.0);
        // Cubes + zones + provenance are populated.
        assert_eq!(m.property_names().len(), 3);
        assert_eq!(m.zones().zones().len(), 1);
        assert!(m.provenance().realization.is_none());
    }

    #[test]
    fn structural_high_via_control_raises_in_place() {
        // Contact mid-column; a shallow crest pulls more rock above it (H3's
        // closure behaviour, regression-locked in the relocated pipeline).
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut b = StaticModelBuilder::from_wireframe(&wf, opts()).unwrap();
        let before = b.build().unwrap().in_place().unwrap().hcpv_m3;
        b.add_top_control(5, 5, 4980.0);
        let after = b.build().unwrap().in_place().unwrap().hcpv_m3;
        assert!(
            after > before,
            "structural high should raise in-place: {before} -> {after}"
        );
    }

    #[test]
    fn logs_populate_cells_in_their_depth_range() {
        // Column spans 5000..5050; logs only in the upper half with φ=0.30
        // distinct from the 0.25 prior.
        let wf = flat_wireframe(11, 5000.0, 5100.0);
        let samples: Vec<PetroSample> = (0..=25)
            .map(|i| (5000.0 + f64::from(i), 0.30, 0.20))
            .collect();
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .with_logs(samples)
            .build()
            .unwrap();
        let poro = &m.property("PORO").unwrap().values;
        assert!(poro.iter().any(|&v| (v - 0.30).abs() < 1e-9), "log cells");
        assert!(poro.iter().any(|&v| (v - 0.25).abs() < 1e-9), "prior cells");
        assert!(matches!(m.provenance().population, PopulationMode::Logs));
    }

    #[test]
    fn builder_needs_a_top_and_a_contact() {
        let mut wf = flat_wireframe(4, 5000.0, 5100.0);
        std::sync::Arc::make_mut(&mut wf.horizons).clear();
        assert!(StaticModelBuilder::from_wireframe(&wf, opts()).is_err());
        let mut wf2 = flat_wireframe(4, 5000.0, 5100.0);
        wf2.contacts.clear();
        assert!(StaticModelBuilder::from_wireframe(&wf2, opts()).is_err());
    }

    // --- Fix 1: real base-horizon relief wired through the build ---

    /// A wedge: flat Top at `top_m`, a `Base` horizon dipping in i from
    /// `top_m + thin_m` (column 0) to `top_m + thick_m` (column n-1),
    /// constant in j; contact deep so the whole column counts.
    fn wedge_wireframe(n: usize, top_m: f64, thin_m: f64, thick_m: f64) -> Wireframe {
        let mut base_depth = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                let frac = c as f64 / (n - 1) as f64;
                base_depth[r * n + c] = top_m + thin_m + (thick_m - thin_m) * frac;
            }
        }
        let mut wf = flat_wireframe(n, top_m, top_m + 1000.0);
        std::sync::Arc::make_mut(&mut wf.horizons).push(Horizon {
            name: "base".into(),
            role: HorizonRole::Base,
            surface: GriddedDepth {
                ncol: n,
                nrow: n,
                depth_m: base_depth,
                is_control: vec![true; n * n],
            },
        });
        wf
    }

    #[test]
    fn base_horizon_relief_drives_spatially_varying_gross() {
        use crate::grid::Ijk;
        // Thickness 20 ft (updip, column 0) -> 120 ft (downdip, column ni). Mean
        // thickness of a linear wedge = (20+120)/2 = 70 ft, so the analytic bulk
        // is footprint_area * 70 — NOT the constant gross_height_m (50) offset.
        let wf = wedge_wireframe(11, 5000.0, 20.0, 120.0);
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();

        let footprint = 100.0;
        let expected = footprint * 70.0;
        assert!(
            (m.bulk_volume() - expected).abs() / expected < 1e-6,
            "wedge bulk {} != analytic {expected}",
            m.bulk_volume()
        );
        // The constant-offset (fallback) answer would be footprint*50 — the base
        // horizon must move the number well away from it.
        let constant_offset = footprint * 50.0;
        assert!(
            (m.bulk_volume() - constant_offset).abs() / constant_offset > 0.3,
            "base relief was ignored: bulk {} ~= constant-offset {constant_offset}",
            m.bulk_volume()
        );
        // Gross varies across i: the top layer thickens from column 0 to column ni-1.
        let dz_updip = m.grid().cell(Ijk::new(0, 0, 0)).dz();
        let dz_downdip = m.grid().cell(Ijk::new(9, 0, 0)).dz();
        assert!(
            dz_downdip > dz_updip * 3.0,
            "gross should thicken downdip: updip {dz_updip} vs downdip {dz_downdip}"
        );
        // A cleanly-consumed Base raises no advisory.
        assert!(
            m.provenance().warnings.is_empty(),
            "clean wedge build has no warnings"
        );
    }

    #[test]
    fn no_base_horizon_falls_back_to_constant_offset() {
        // Backward compat: a Top-only wireframe still gives the gross_height_m box.
        let wf = flat_wireframe(11, 5000.0, 6000.0);
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();
        let expected = 100.0 * 50.0; // gross_height_m
        assert!((m.bulk_volume() - expected).abs() / expected < 1e-9);
        assert!(m.provenance().warnings.is_empty());
    }

    #[test]
    fn unused_intermediate_horizon_is_warned() {
        use crate::model::BuildWarning;
        let mut wf = flat_wireframe(11, 5000.0, 6000.0);
        std::sync::Arc::make_mut(&mut wf.horizons).push(Horizon {
            name: "mid".into(),
            role: HorizonRole::Intermediate,
            surface: GriddedDepth {
                ncol: 11,
                nrow: 11,
                depth_m: vec![5030.0; 11 * 11],
                is_control: vec![true; 11 * 11],
            },
        });
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();
        assert!(
            m.provenance().warnings.iter().any(|w| matches!(
                w,
                BuildWarning::UnusedHorizon {
                    role: HorizonRole::Intermediate,
                    ..
                }
            )),
            "intermediate horizon should raise an UnusedHorizon warning: {:?}",
            m.provenance().warnings
        );
        // Base still falls back to the constant offset.
        let expected = 100.0 * 50.0;
        assert!((m.bulk_volume() - expected).abs() / expected < 1e-9);
    }

    // --- template gross scaling over real base relief
    //     (`decision_template_gross_scaling`) ---

    #[test]
    fn wedge_template_at_mean_gross_matches_deterministic_build() {
        // Wedge g(x): 20 -> 120 ft, mean(g) = 70. A draw with gross == mean(g)
        // (scale 1) must reproduce the deterministic Fix-1 build exactly: the
        // fully-pinned lattice makes both solvers return the control values, so
        // the surfaces — and the volumes — coincide bit-for-bit.
        let wf = wedge_wireframe(11, 5000.0, 20.0, 120.0);
        let det = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let draw = RealizationDraw::new(100.0, 70.0, 6000.0, 0.25, 0.8, 0.3, 0);
        let real = t.realize(&draw).unwrap();
        let (a, b) = (det.bulk_volume(), real.bulk_volume());
        assert!(
            (a - b).abs() <= 1e-12 * a,
            "template at mean gross diverged from deterministic build: {a} vs {b}"
        );
        // Per-cell geometry matches too (spot-check both wedge ends, top layer).
        use crate::grid::Ijk;
        for ijk in [Ijk::new(0, 0, 0), Ijk::new(9, 0, 0)] {
            let (dd, rd) = (det.grid().cell(ijk).dz(), real.grid().cell(ijk).dz());
            assert!(
                (dd - rd).abs() <= 1e-12 * dd.max(1.0),
                "cell dz {dd} vs {rd}"
            );
        }
    }

    #[test]
    fn wedge_template_gross_scales_level_preserving_shape() {
        use crate::grid::Ijk;
        let wf = wedge_wireframe(11, 5000.0, 20.0, 120.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let at_mean = t
            .realize(&RealizationDraw::new(
                100.0, 70.0, 6000.0, 0.25, 0.8, 0.3, 0,
            ))
            .unwrap();
        let doubled = t
            .realize(&RealizationDraw::new(
                100.0, 140.0, 6000.0, 0.25, 0.8, 0.3, 1,
            ))
            .unwrap();
        // Level: bulk doubles; every column's gross doubles.
        assert!(
            (doubled.bulk_volume() - 2.0 * at_mean.bulk_volume()).abs()
                <= 1e-9 * doubled.bulk_volume()
        );
        let dz = |m: &StaticModel, i| m.grid().cell(Ijk::new(i, 0, 0)).dz();
        for i in [0usize, 5, 9] {
            let (one, two) = (dz(&at_mean, i), dz(&doubled, i));
            assert!(
                (two - 2.0 * one).abs() <= 1e-9 * two,
                "col {i}: {one} -> {two}"
            );
        }
        // Shape: the updip/downdip thickness ratio is preserved.
        let r1 = dz(&at_mean, 9) / dz(&at_mean, 0);
        let r2 = dz(&doubled, 9) / dz(&doubled, 0);
        assert!((r1 - r2).abs() < 1e-9, "shape ratio drifted: {r1} vs {r2}");
        // Two-contact in_place works identically over the scaled base.
        let split = t
            .realize(&RealizationDraw::new(100.0, 70.0, 5100.0, 0.25, 0.8, 0.3, 2).with_goc(5030.0))
            .unwrap()
            .in_place()
            .unwrap();
        let (g, o) = (split.gas.unwrap(), split.oil.unwrap());
        assert!(g.cells > 0 && o.cells > 0);
        assert!((split.hcpv_m3 - (g.hcpv_m3 + o.hcpv_m3)).abs() <= 1e-9 * split.hcpv_m3);
    }

    /// A **sparse** top (only 5 defined nodes; the interior interpolates) with a
    /// thin base a constant `sep_m` below it at those same nodes. Sparse +
    /// interpolated is exactly the shape the old cold-builder / warm-template split
    /// diverged on (R2): a fully-pinned lattice hid the kernel mismatch because
    /// both solvers returned control values.
    fn sparse_thin_wireframe(n: usize, sep_m: f64) -> Wireframe {
        let mut top = vec![f64::NAN; n * n];
        let mut base = vec![f64::NAN; n * n];
        let defined = [
            (0, 0, 5000.0),
            (n - 1, 0, 5020.0),
            (0, n - 1, 5010.0),
            (n - 1, n - 1, 5030.0),
            (n / 2, n / 2, 4995.0),
        ];
        for &(c, r, z) in &defined {
            top[r * n + c] = z;
            base[r * n + c] = z + sep_m;
        }
        let mut wf = flat_wireframe(n, 5000.0, 9000.0);
        std::sync::Arc::make_mut(&mut wf.horizons)[0].surface = GriddedDepth {
            ncol: n,
            nrow: n,
            depth_m: top,
            is_control: vec![true; n * n],
        };
        std::sync::Arc::make_mut(&mut wf.horizons).push(Horizon {
            name: "base".into(),
            role: HorizonRole::Base,
            surface: GriddedDepth {
                ncol: n,
                nrow: n,
                depth_m: base,
                is_control: vec![true; n * n],
            },
        });
        wf
    }

    #[test]
    fn thin_sparse_template_at_mean_gross_matches_deterministic_build() {
        // R2 regression lock: on a NON-fully-pinned, thin (2 ft) column the unified
        // kernel path makes builder and template coincide. A constant 2 ft base
        // offset over sparse controls means the kernel reproduces it exactly, so
        // mean(g) == 2.0 and a scale-1 draw reproduces the build.
        let wf = sparse_thin_wireframe(11, 2.0);
        let det = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let draw = RealizationDraw::new(100.0, 2.0, 9000.0, 0.25, 0.8, 0.3, 0);
        let real = t.realize(&draw).unwrap();
        let (a, b) = (det.bulk_volume(), real.bulk_volume());
        assert!(a > 0.0, "thin column has positive bulk");
        assert!(
            (a - b).abs() <= 1e-9 * a,
            "thin sparse template diverged from build: {a} vs {b}"
        );
        // Per-cell geometry matches across the interpolated interior too.
        use crate::grid::Ijk;
        for ijk in [Ijk::new(0, 0, 0), Ijk::new(5, 5, 0), Ijk::new(9, 9, 0)] {
            let (dd, rd) = (det.grid().cell(ijk).dz(), real.grid().cell(ijk).dz());
            assert!(
                (dd - rd).abs() <= 1e-9 * dd.max(1.0),
                "cell dz {dd} vs {rd}"
            );
        }
    }

    // --- R1: base-above-top guard (crossing surfaces collapse GRV) ---

    /// A crossing base: flat Top at `top_m`; the Base sits 10 ft ABOVE the top on
    /// the left half (a crossing that collapses gross) and 10 ft below on the right.
    fn crossing_wireframe(n: usize, top_m: f64) -> Wireframe {
        let mut base_depth = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                base_depth[r * n + c] = if c < n / 2 {
                    top_m - 10.0
                } else {
                    top_m + 10.0
                };
            }
        }
        let mut wf = flat_wireframe(n, top_m, top_m + 1000.0);
        std::sync::Arc::make_mut(&mut wf.horizons).push(Horizon {
            name: "base".into(),
            role: HorizonRole::Base,
            surface: GriddedDepth {
                ncol: n,
                nrow: n,
                depth_m: base_depth,
                is_control: vec![true; n * n],
            },
        });
        wf
    }

    #[test]
    fn crossing_base_errors_by_default() {
        use crate::error::StaticError;
        let wf = crossing_wireframe(11, 5000.0);
        let err = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap_err();
        match err {
            StaticError::CrossedSurfaces { nodes, worst_m } => {
                assert!(nodes > 0, "reports offending nodes");
                assert!(worst_m < 0.0, "worst crossing is negative: {worst_m}");
            }
            other => panic!("expected CrossedSurfaces, got {other}"),
        }
        // The template path guards identically per realization.
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let draw = RealizationDraw::new(100.0, 20.0, 6000.0, 0.25, 0.8, 0.3, 0);
        assert!(matches!(
            t.realize(&draw),
            Err(StaticError::CrossedSurfaces { .. })
        ));
    }

    #[test]
    fn crossing_base_clamps_offending_columns_only() {
        use crate::grid::Ijk;
        // Opt in: the crossed (left) columns clamp to zero gross; the good (right)
        // columns keep their thickness untouched.
        let wf = crossing_wireframe(11, 5000.0);
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .with_clamp_base_to_top(true)
            .build()
            .unwrap();
        let dz_left = m.grid().cell(Ijk::new(0, 0, 0)).dz();
        let dz_right = m.grid().cell(Ijk::new(9, 0, 0)).dz();
        assert!(dz_left.abs() < 1e-6, "crossed column zeroed: dz={dz_left}");
        assert!(dz_right > 1.0, "good column preserved: dz={dz_right}");
    }

    #[test]
    fn min_thickness_repairs_crossing_and_records_a_warning() {
        use crate::error::StaticError;
        use crate::grid::Ijk;
        use crate::model::provenance::BuildWarning;
        let wf = crossing_wireframe(11, 5000.0);
        // Default: the crossing wedge still errors (the guard stays the default).
        assert!(matches!(
            StaticModelBuilder::from_wireframe(&wf, opts())
                .unwrap()
                .build(),
            Err(StaticError::CrossedSurfaces { .. })
        ));

        // Opt in: post-gridding repair pulls the crossed (left) columns down to a
        // minimum thickness (top preserved) and builds cleanly.
        let min_t = 2.0;
        let nk = 5; // opts().nk
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .with_min_thickness_m(min_t)
            .build()
            .unwrap();
        let col_thick = |model: &StaticModel, i: usize| {
            (0..nk)
                .map(|k| model.grid().cell(Ijk::new(i, 0, k)).dz())
                .sum::<f64>()
        };
        // A fully-crossed left column is repaired to exactly the minimum thickness.
        assert!(
            (col_thick(&m, 0) - min_t).abs() < 1e-6,
            "crossed column repaired to min thickness: {}",
            col_thick(&m, 0)
        );
        // A non-crossing right column is bit-unchanged vs the (untouched) guarded
        // geometry — the clamp build leaves non-crossing columns alone identically.
        let clean = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .with_clamp_base_to_top(true)
            .build()
            .unwrap();
        assert_eq!(
            col_thick(&m, 9),
            col_thick(&clean, 9),
            "non-crossing column bit-unchanged"
        );
        // The warning carries the repaired-column count + the worst (negative,
        // a true crossing) violation.
        let (columns, worst_m) = m
            .provenance()
            .warnings
            .iter()
            .find_map(|w| match w {
                BuildWarning::ThinColumnsRepaired { columns, worst_m } => {
                    Some((*columns, *worst_m))
                }
                _ => None,
            })
            .expect("a ThinColumnsRepaired warning was recorded");
        assert!(columns > 0, "repaired-column count populated");
        assert!(
            worst_m < 0.0,
            "worst violation is a true crossing (negative): {worst_m}"
        );
    }

    #[test]
    fn template_min_thickness_repairs_per_realization() {
        use crate::model::provenance::BuildWarning;
        let wf = crossing_wireframe(11, 5000.0);
        // The template path honours the same opt-in remedy per realization.
        let mut t = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_min_thickness_m(2.0);
        let draw = RealizationDraw::new(100.0, 20.0, 6000.0, 0.25, 0.8, 0.3, 0);
        let m = t
            .realize(&draw)
            .expect("min-thickness realization builds cleanly");
        let (columns, worst_m) = m
            .provenance()
            .warnings
            .iter()
            .find_map(|w| match w {
                BuildWarning::ThinColumnsRepaired { columns, worst_m } => {
                    Some((*columns, *worst_m))
                }
                _ => None,
            })
            .expect("a ThinColumnsRepaired warning was recorded on the realization");
        assert!(columns > 0, "repaired-column count populated");
        assert!(worst_m < 0.0, "worst violation negative: {worst_m}");
    }

    #[test]
    fn no_base_template_keeps_constant_offset() {
        // Backward compat: without a Base horizon the draw's gross_height_m is
        // the exact constant offset it always was.
        let wf = flat_wireframe(11, 5000.0, 6000.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let m = t
            .realize(&RealizationDraw::new(
                100.0, 50.0, 6000.0, 0.25, 0.8, 0.3, 0,
            ))
            .unwrap();
        let expected = 100.0 * 50.0;
        assert!((m.bulk_volume() - expected).abs() / expected < 1e-9);
    }

    // --- Fix 2: minimal areal trend hook (external-drift-lite) ---

    #[test]
    fn uniform_trend_is_a_noop() {
        // A trend of all-1.0 multipliers must leave NTG at the prior everywhere.
        let wf = flat_wireframe(11, 5000.0, 6000.0);
        let trend = TrendSurface::new(10, 10, vec![1.0; 100]).unwrap();
        let plain = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();
        let trended = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .with_areal_trend(trend)
            .build()
            .unwrap();
        let a = &plain.property("NTG").unwrap().values;
        let b = &trended.property("NTG").unwrap().values;
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b) {
            assert!(
                (x - y).abs() < 1e-12,
                "uniform trend changed NTG: {x} vs {y}"
            );
        }
    }

    #[test]
    fn step_trend_shifts_ntg_preserving_mean() {
        use crate::grid::Ijk;
        // 10x10 columns; a step trend (left half 1.0, right half 2.0) has mean
        // 1.5, so multipliers are 2/3 (left) and 4/3 (right). Base NTG 0.5 ->
        // 0.333 / 0.667, and the field mean is preserved at 0.5.
        let mut o = opts();
        o.priors.net_to_gross = 0.5;
        let wf = flat_wireframe(11, 5000.0, 6000.0);
        let mut vals = vec![0.0; 100];
        for r in 0..10 {
            for c in 0..10 {
                vals[r * 10 + c] = if c < 5 { 1.0 } else { 2.0 };
            }
        }
        let trend = TrendSurface::new(10, 10, vals).unwrap();
        let m = StaticModelBuilder::from_wireframe(&wf, o)
            .unwrap()
            .with_areal_trend(trend)
            .build()
            .unwrap();
        let ntg = &m.property("NTG").unwrap().values;
        let dims = m.grid().dims();
        let at = |i, j, k| ntg[dims.linear(Ijk::new(i, j, k)).unwrap()];
        assert!(
            (at(0, 0, 0) - 0.5 * (2.0 / 3.0)).abs() < 1e-9,
            "left NTG {}",
            at(0, 0, 0)
        );
        assert!(
            (at(9, 0, 0) - 0.5 * (4.0 / 3.0)).abs() < 1e-9,
            "right NTG {}",
            at(9, 0, 0)
        );
        // Field mean preserved at the prior level.
        let mean = ntg.iter().sum::<f64>() / ntg.len() as f64;
        assert!(
            (mean - 0.5).abs() < 1e-9,
            "trend shifted the field mean to {mean}"
        );
    }

    #[test]
    fn trend_with_porosity_flag_modulates_both_cubes() {
        let mut o = opts();
        o.priors.net_to_gross = 0.4;
        o.priors.porosity = 0.2;
        let wf = flat_wireframe(11, 5000.0, 6000.0);
        let mut vals = vec![1.0; 100];
        for r in 0..10 {
            vals[r * 10 + 9] = 1.5; // one column boosted
        }
        let trend = TrendSurface::new(10, 10, vals).unwrap().with_porosity();
        let m = StaticModelBuilder::from_wireframe(&wf, o)
            .unwrap()
            .with_areal_trend(trend)
            .build()
            .unwrap();
        // Porosity must vary now (flag on) — not a single constant.
        let poro = &m.property("PORO").unwrap().values;
        let pmin = poro.iter().cloned().fold(f64::INFINITY, f64::min);
        let pmax = poro.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(pmax > pmin + 1e-6, "porosity should vary with the flag on");
    }

    // --- Fix 3: two-contact volumetrics (gas cap + oil rim) ---

    #[test]
    fn two_contact_framework_splits_in_place() {
        // Column 5000-5050 (opts gross 50, 5 layers of 10 ft; centroids
        // 5005..5045). GOC 5020 -> gas 5005/5015; OWC 5040 -> oil 5025/5035.
        let mut wf = flat_wireframe(11, 5000.0, 5040.0); // the OWC
        wf.contacts.insert(
            0,
            Contact {
                kind: ContactKind::Goc,
                depth_m: 5020.0,
                hardness: Hardness::Hard,
            },
        );
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();
        let ip = m.in_place().unwrap();
        let gas = ip.gas.expect("gas zone present");
        let oil = ip.oil.expect("oil zone present");
        assert_eq!(gas.cells, 10 * 10 * 2, "2 gas layers");
        assert_eq!(oil.cells, 10 * 10 * 2, "2 oil layers");
        assert!((ip.hcpv_m3 - (gas.hcpv_m3 + oil.hcpv_m3)).abs() < 1e-6);
        // Per-zone surface volumes both positive with the right FVF.
        assert!(ip.gas_zone_ogip_sm3(GasFvf::new(0.004).unwrap()) > 0.0);
        assert!(ip.oil_zone_ooip_sm3(OilFvf::new(1.25).unwrap()) > 0.0);
    }

    #[test]
    fn sw_gas_override_lowers_gas_zone_ogip() {
        // R3: with_sw_gas gives the gas cap a lower connate water than the shared
        // SW cube, raising gas-cap HCPV/OGIP and leaving the oil leg untouched.
        let mut wf = flat_wireframe(11, 5000.0, 5040.0);
        wf.contacts.insert(
            0,
            Contact {
                kind: ContactKind::Goc,
                depth_m: 5020.0,
                hardness: Hardness::Hard,
            },
        );
        let base = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap()
            .in_place()
            .unwrap();
        let over = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .with_sw_gas(0.1) // cube Sw is 0.3
            .build()
            .unwrap()
            .in_place()
            .unwrap();
        assert!(
            over.gas.unwrap().hcpv_m3 > base.gas.unwrap().hcpv_m3,
            "lower gas Sw raises gas HCPV"
        );
        // Oil leg identical; only the gas cap moved.
        assert!(
            (over.oil.unwrap().hcpv_m3 - base.oil.unwrap().hcpv_m3).abs() < 1e-6,
            "oil leg untouched by sw_gas"
        );
        // The template path honours draw.sw_gas too (draw carries the GOC).
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let draw = RealizationDraw::new(100.0, 50.0, 5040.0, 0.25, 0.8, 0.3, 0)
            .with_goc(5020.0)
            .with_sw_gas(0.1);
        let tip = t.realize(&draw).unwrap().in_place().unwrap();
        assert!(tip.gas.unwrap().hcpv_m3 > base.gas.unwrap().hcpv_m3);
    }

    #[test]
    fn in_place_summary_matches_full_aggregates() {
        // V7: the summary path returns identical aggregates but leaves the
        // per-cell HCPV cube empty (single- and two-contact).
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let m = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap();
        let full = m.in_place().unwrap();
        let sum = m.in_place_summary().unwrap();
        assert!((full.grv_m3 - sum.grv_m3).abs() < 1e-9);
        assert!((full.hcpv_m3 - sum.hcpv_m3).abs() < 1e-9);
        assert_eq!(full.cells_in_column, sum.cells_in_column);
        assert!(!full.per_cell_hcpv.is_empty());
        assert!(
            sum.per_cell_hcpv.is_empty(),
            "summary skips the per-cell cube"
        );

        // Two-contact path too.
        let mut wf2 = flat_wireframe(11, 5000.0, 5040.0);
        wf2.contacts.insert(
            0,
            Contact {
                kind: ContactKind::Goc,
                depth_m: 5020.0,
                hardness: Hardness::Hard,
            },
        );
        let m2 = StaticModelBuilder::from_wireframe(&wf2, opts())
            .unwrap()
            .build()
            .unwrap();
        let f2 = m2.in_place().unwrap();
        let s2 = m2.in_place_summary().unwrap();
        assert!((f2.gas.unwrap().hcpv_m3 - s2.gas.unwrap().hcpv_m3).abs() < 1e-9);
        assert!((f2.oil.unwrap().hcpv_m3 - s2.oil.unwrap().hcpv_m3).abs() < 1e-9);
        assert!(s2.per_cell_hcpv.is_empty());
    }

    #[test]
    fn single_contact_framework_stays_generic() {
        // Backward compat: a lone OWC yields no gas/oil split (generic column).
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let ip = StaticModelBuilder::from_wireframe(&wf, opts())
            .unwrap()
            .build()
            .unwrap()
            .in_place()
            .unwrap();
        assert!(ip.gas.is_none() && ip.oil.is_none());
        assert!(ip.ooip_sm3(OilFvf::new(1.25).unwrap()) > 0.0);
    }

    #[test]
    fn draw_with_goc_realizes_a_two_contact_model() {
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        // Gas cap above 5020, oil rim down to the drawn OWC 5040.
        let draw = RealizationDraw::new(100.0, 50.0, 5040.0, 0.25, 0.8, 0.3, 0).with_goc(5020.0);
        let ip = t.realize(&draw).unwrap().in_place().unwrap();
        assert!(
            ip.gas.is_some() && ip.oil.is_some(),
            "GOC draw splits the column"
        );
        assert!(ip.gas.unwrap().cells > 0 && ip.oil.unwrap().cells > 0);
        // A GOC below the OWC is a typed error.
        let bad = RealizationDraw::new(100.0, 50.0, 5030.0, 0.25, 0.8, 0.3, 1).with_goc(5060.0);
        assert!(t.realize(&bad).is_err());
    }

    // --- the regeneration seam (template + draw) ---

    #[test]
    fn template_realizes_n100_and_aggregates_a_p_curve() {
        // THE seam smoke test: one template, 100 draws through `realize`,
        // in-place aggregated to a P-curve — proving template reuse, the warm
        // chain, per-draw provenance, and the volumetrics output surface.
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let boi = OilFvf::new(1.25).unwrap();
        let mut ooip = Vec::with_capacity(100);
        for i in 0..100u64 {
            // A deterministic pseudo-sampler: small scalar perturbations around
            // the nominal inputs (the real sampler fills draws upstream).
            let f = (i % 10) as f64 / 10.0; // 0.0 .. 0.9
            let draw = RealizationDraw::new(
                90.0 + 20.0 * f,   // area
                45.0 + 10.0 * f,   // gross height
                5020.0 + 10.0 * f, // contact
                0.20 + 0.10 * f,   // porosity
                0.75 + 0.10 * f,   // NTG
                0.25 + 0.10 * f,   // Sw
                i,                 // seed index
            );
            let m = t.realize(&draw).unwrap();
            assert_eq!(m.provenance().realization.as_ref().unwrap().seed_index, i);
            ooip.push(m.in_place().unwrap().ooip_sm3(boi));
        }
        let s = PercentileSummary::from_realizations(&ooip).unwrap();
        assert!(s.p90 < s.p50 && s.p50 < s.p10, "P-curve ordered: {s:?}");
        assert!(s.p90 > 0.0);
    }

    #[test]
    fn template_areal_trend_shifts_ntg_preserving_mean() {
        use crate::grid::Ijk;
        // The template analog of the builder step-trend test: a step trend on the
        // MC path shifts NTG laterally while preserving the field mean.
        let wf = flat_wireframe(11, 5000.0, 6000.0);
        let mut vals = vec![0.0; 100];
        for r in 0..10 {
            for c in 0..10 {
                vals[r * 10 + c] = if c < 5 { 1.0 } else { 2.0 };
            }
        }
        let trend = TrendSurface::new(10, 10, vals).unwrap();
        let mut t = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_areal_trend(trend);
        let draw = RealizationDraw::new(100.0, 50.0, 6000.0, 0.25, 0.5, 0.3, 0);
        let m = t.realize(&draw).unwrap();
        let ntg = &m.property("NTG").unwrap().values;
        let dims = m.grid().dims();
        let at = |i, j, k| ntg[dims.linear(Ijk::new(i, j, k)).unwrap()];
        assert!(
            (at(0, 0, 0) - 0.5 * (2.0 / 3.0)).abs() < 1e-9,
            "left {}",
            at(0, 0, 0)
        );
        assert!(
            (at(9, 0, 0) - 0.5 * (4.0 / 3.0)).abs() < 1e-9,
            "right {}",
            at(9, 0, 0)
        );
        let mean = ntg.iter().sum::<f64>() / ntg.len() as f64;
        assert!((mean - 0.5).abs() < 1e-9, "field mean drifted to {mean}");
    }

    #[test]
    fn identical_draws_realize_identical_models() {
        // Warm==cold within the kernel: chaining realize() with the same draw
        // must reproduce the same converged surface, hence the same volumes.
        let wf = flat_wireframe(9, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let draw = RealizationDraw::new(100.0, 50.0, 5025.0, 0.25, 0.8, 0.3, 7);
        let a = t.realize(&draw).unwrap().in_place().unwrap().hcpv_m3;
        let b = t.realize(&draw).unwrap().in_place().unwrap().hcpv_m3;
        assert!(
            (a - b).abs() <= 1e-6 * a.abs(),
            "chained identical draws diverged: {a} vs {b}"
        );
    }

    #[test]
    fn structural_perturbation_shifts_a_control() {
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let nominal = RealizationDraw::new(100.0, 50.0, 5025.0, 0.25, 0.8, 0.3, 0);
        let before = t.realize(&nominal).unwrap().in_place().unwrap().hcpv_m3;
        // Raise the centre node 20 ft (shallower) -> more rock above the contact.
        let crest = nominal.clone().with_structural(StructuralPerturbation {
            control_shifts: vec![(5, 5, -20.0)],
        });
        let after = t.realize(&crest).unwrap().in_place().unwrap().hcpv_m3;
        assert!(
            after > before,
            "crest shift should raise in-place: {before} -> {after}"
        );
        // Off-lattice shift is a typed error.
        let bad = nominal.with_structural(StructuralPerturbation {
            control_shifts: vec![(99, 99, -20.0)],
        });
        assert!(t.realize(&bad).is_err());
    }

    /// Bit-compare a recycled model against a freshly-realized one: geometry (ZCORN
    /// via per-cell centroid-z + volume) and every cube, exact `f64` bits.
    fn assert_models_bit_identical(a: &StaticModel, b: &StaticModel, ctx: &str) {
        assert_eq!(
            a.grid().cell_count(),
            b.grid().cell_count(),
            "{ctx}: cell count"
        );
        for lin in 0..a.grid().cell_count() {
            assert_eq!(
                a.grid().cell_centroid_z_at(lin).to_bits(),
                b.grid().cell_centroid_z_at(lin).to_bits(),
                "{ctx}: ZCORN centroid-z diverged at cell {lin}"
            );
            assert_eq!(
                a.grid().cell_volume_at(lin).to_bits(),
                b.grid().cell_volume_at(lin).to_bits(),
                "{ctx}: cell volume diverged at cell {lin}"
            );
        }
        for name in ["PORO", "SW", "NTG"] {
            assert_eq!(
                a.property(name).map(|p| &p.values),
                b.property(name).map(|p| &p.values),
                "{ctx}: cube {name} diverged"
            );
        }
    }

    #[test]
    fn realize_into_recycles_buffers_without_staleness() {
        // Trap #5: `realize_into` reuses the passed model's ZCORN + cube buffers, so a
        // second draw with a DIFFERENT gross/porosity must FULLY overwrite the first
        // draw's geometry + cubes — never accumulate or leave a stale tail. Proof: run
        // two consecutive draws (A then B) into one reused model and bit-compare the
        // result against a fresh `realize` of the SAME chain history (A discarded, then
        // B) — identical chains, so any difference is a recycling bug.
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        // A PORO LevelShift pipeline + per-draw shifts also exercises the recycled
        // shift-cube buffer.
        let poro = || {
            use petektools::{Variogram, VariogramModel};
            let low = WellLog::new(0.5, 0.5, vec![(5005.0, 0.10), (5045.0, 0.18)]);
            let high = WellLog::new(9.5, 9.5, vec![(5005.0, 0.26), (5045.0, 0.22)]);
            let vgm = Variogram::new(VariogramModel::Spherical, 0.0, 1.0, 5.0).unwrap();
            PropertyPipeline::new("PORO")
                .upscale(vec![low, high], UpscaleMethod::Arithmetic)
                // Per-draw geometry can leave a simulated layer data-less; opt into
                // the mean-fill rather than the default hard error (item 4).
                .propagate(Gaussian::new(vgm, 42).allow_mean_fill())
        };
        // Different gross (50 vs 80) -> different base/layering geometry; different
        // porosity + PORO shift -> different cubes.
        let draw_a = RealizationDraw::new(100.0, 50.0, 5025.0, 0.20, 0.80, 0.30, 1)
            .with_property_shift("PORO", -0.03);
        let draw_b = RealizationDraw::new(120.0, 80.0, 5030.0, 0.26, 0.85, 0.25, 2)
            .with_property_shift("PORO", 0.04);

        // Reference: fresh allocations, chain A→B.
        let mut t_ref = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property(poro());
        let model_a = t_ref.realize(&draw_a).unwrap();
        let model_ref_b = t_ref.realize(&draw_b).unwrap();

        // Recycled: one reused model, chain A→B.
        let mut t_rec = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property(poro());
        let mut m = t_rec.reusable_model();
        t_rec.realize_into(&draw_a, &mut m).unwrap();
        // After draw A the buffers hold A's geometry + cubes (the stale content that B
        // must overwrite) — guard that A and B genuinely differ, else the test is inert.
        assert_ne!(
            m.property("PORO").unwrap().values,
            model_ref_b.property("PORO").unwrap().values,
            "draws A and B must differ for the staleness test to bite"
        );
        assert_models_bit_identical(&model_a, &m, "draw A: recycled vs fresh");

        t_rec.realize_into(&draw_b, &mut m).unwrap();
        // The reused buffers now carry B EXACTLY — no stale A tail, bit-for-bit.
        assert_models_bit_identical(&model_ref_b, &m, "draw B: recycled (post-A) vs fresh");
    }

    #[test]
    fn realize_rejects_a_garbage_draw() {
        // H2 at the seam.
        let wf = flat_wireframe(9, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut bad = RealizationDraw::new(100.0, 50.0, 5025.0, -0.1, 0.8, 0.3, 0);
        assert!(t.realize(&bad).is_err()); // φ = -0.1
        bad.porosity = 0.25;
        bad.water_saturation = 1.2;
        assert!(t.realize(&bad).is_err()); // Sw = 1.2
        bad.water_saturation = 0.3;
        bad.area_m2 = 0.0;
        assert!(t.realize(&bad).is_err()); // area = 0
        bad.area_m2 = 100.0;
        bad.contact_depth_m = f64::INFINITY;
        assert!(t.realize(&bad).is_err()); // non-finite contact
    }
}
