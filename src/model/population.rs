//! Property population — the population half of the relocated `RefiningModel`
//! pipeline (petekSim `srs-core/src/refine.rs`, moved here 2026-07-03 per
//! `task_relocate_refine_orchestration`).
//!
//! For each cell, positioned petro samples whose TVD falls in the cell's depth
//! range are upscaled (φ length-weighted, Sw pore-volume-weighted via srs-petro);
//! cells with no samples fall back to the constant priors. NTG is the prior
//! everywhere (logs carry φ/Sw; net flags are petekIO's interpret step).

use crate::error::StaticError;
use crate::grid::{Grid, Property};
use crate::model::model::Georef;
use crate::model::trend::TrendSurface;
use crate::petro::{upscale_porosity, upscale_sw, SwSample, WeightedSample};
use crate::volumetrics::{populate_constant, validate_fraction, ConstantPriors, NTG, PORO, SW};

/// A positioned petrophysical sample for grid population: model-internal
/// positive-down depth \[m\] + porosity + water saturation (fractions). Built by
/// the data layer (srs-data) from petekIO's positioned well curves
/// (`WellCurveInput.xyz`, negative-down elevation) and aligned to this
/// positive-down `depth_m` datum at the srs-core binning seam; binned into cells
/// here. Lateral placement is deferred — every column sees the same samples
/// (a single-well layer-cake); areal heterogeneity (georeferenced placement +
/// kriging) is the P5 property-modelling work.
pub type PetroSample = (f64, f64, f64);

/// Sort petro samples by TVD ascending — the invariant [`populate_from_logs`]
/// relies on to binary-search each cell's depth range instead of scanning all
/// samples (V2). The `with_logs` setters call this once so the cost is paid a
/// single time, not per realization.
pub(crate) fn sort_by_tvd(mut samples: Vec<PetroSample>) -> Vec<PetroSample> {
    samples.sort_by(|a, b| a.0.total_cmp(&b.0));
    samples
}

/// Populate `PORO`/`SW`/`NTG`: constants from `priors`, or upscaled logs where
/// samples cover a cell's depth range (priors elsewhere).
///
/// # Errors
/// [`StaticError::InvalidInput`] if a prior is not a fraction in `[0, 1]` (H2),
/// or if upscaling rejects a sample set.
pub(crate) fn populate(
    grid: &mut Grid,
    priors: ConstantPriors,
    logs: Option<&[PetroSample]>,
    trend: Option<&TrendSurface>,
    georef: Option<Georef>,
) -> Result<(), StaticError> {
    match logs {
        None => populate_constant(grid, priors),
        Some(samples) => populate_from_logs(grid, priors, samples),
    }?;
    if let Some(trend) = trend {
        apply_areal_trend(grid, trend, georef)?;
    }
    Ok(())
}

/// Apply an areal trend as a per-column multiplier on the NTG cube (and, if the
/// trend so flags, porosity): the trend gives lateral *shape*, the already-set
/// prior/log value gives the *level*. Multipliers are mean-normalized (see
/// [`TrendSurface`]) so the property field-mean is preserved; results are
/// clamped to `[0, 1]`.
fn apply_areal_trend(
    grid: &mut Grid,
    trend: &TrendSurface,
    georef: Option<Georef>,
) -> Result<(), StaticError> {
    let dims = grid.dims();
    let (ni, nj, nk) = (dims.ni, dims.nj, dims.nk);
    // Resample the trend to the model areal lattice through the SHARED kernel
    // (`petektools::resample`, via `column_multipliers_on`): a georeferenced trend
    // lands by world coordinate, a bare one is stretched to the column extent. The
    // lattice is reconstructed from the top-layer cell centroids.
    let local = crate::model::pipeline::areal_lattice(grid)?;
    let lattice = match georef {
        Some(g) if trend.is_georeferenced() => petektools::Lattice {
            xori: g.origin_x,
            yori: g.origin_y,
            xinc: g.spacing_x,
            yinc: g.spacing_y,
            ncol: ni,
            nrow: nj,
            rotation_deg: g.rotation_deg,
            yflip: g.yflip,
        },
        _ => local,
    };
    let mult = trend.column_multipliers_on(&lattice)?;

    let mut targets = vec![NTG];
    if trend.applies_to_porosity() {
        targets.push(PORO);
    }
    for name in targets {
        let mut prop = grid.properties().get(name).cloned().ok_or_else(|| {
            StaticError::InvalidInput(format!("areal trend target cube '{name}' is not populated"))
        })?;
        for k in 0..nk {
            for j in 0..nj {
                for i in 0..ni {
                    let idx = (k * nj + j) * ni + i;
                    prop.values[idx] = (prop.values[idx] * mult[j * ni + i]).clamp(0.0, 1.0);
                }
            }
        }
        grid.properties_mut().set(prop)?;
    }
    Ok(())
}

/// Overwrite `PORO`/`SW`/`NTG` with constant `priors` across one zone's `k`-range
/// (P8 per-zone population, `task_petekstatic_multizone_2`) — the per-zone
/// distribution level a stack zone owns, applied after the base population. Every
/// cell in the range (including truncated ones, whose values are ignored by
/// volumetrics) is set, so the override is uniform across the zone.
///
/// # Errors
/// [`StaticError::InvalidInput`] if a prior is not a fraction in `[0, 1]` (H2), or a
/// target cube is not yet populated.
pub(crate) fn override_zone_priors(
    grid: &mut Grid,
    priors: ConstantPriors,
    k_range: core::ops::Range<usize>,
) -> Result<(), StaticError> {
    validate_fraction("porosity", priors.porosity)?;
    validate_fraction("net_to_gross", priors.net_to_gross)?;
    validate_fraction("water_saturation", priors.water_saturation)?;
    let dims = grid.dims();
    let (ni, nj) = (dims.ni, dims.nj);
    for (name, val) in [
        (PORO, priors.porosity),
        (SW, priors.water_saturation),
        (NTG, priors.net_to_gross),
    ] {
        let mut prop = grid.properties().get(name).cloned().ok_or_else(|| {
            StaticError::InvalidInput(format!("zone-priors target cube '{name}' is not populated"))
        })?;
        for k in k_range.clone() {
            for j in 0..nj {
                for i in 0..ni {
                    prop.values[(k * nj + j) * ni + i] = val;
                }
            }
        }
        grid.properties_mut().set(prop)?;
    }
    Ok(())
}

fn populate_from_logs(
    grid: &mut Grid,
    priors: ConstantPriors,
    samples: &[PetroSample],
) -> Result<(), StaticError> {
    // The prior is the fallback for uncovered cells — it must be physical too.
    validate_fraction("porosity", priors.porosity)?;
    validate_fraction("net_to_gross", priors.net_to_gross)?;
    validate_fraction("water_saturation", priors.water_saturation)?;

    let n = grid.cell_count();
    // Recycle the existing PORO/SW/NTG cube buffers (take → refill → reinstall):
    // taking them out releases the mutable borrow so the per-cell loop can read the
    // grid geometry, and reusing their capacity means no fresh cube alloc on the
    // steady-state `realize_into` path. A cold build takes empty Vecs (allocates
    // once on `reserve`) — identical outcome.
    let (mut poro, mut sw, mut ntg) = {
        let props = grid.properties_mut();
        (
            props.take_values(PORO),
            props.take_values(SW),
            props.take_values(NTG),
        )
    };
    poro.clear();
    poro.reserve(n);
    sw.clear();
    sw.reserve(n);
    ntg.clear();
    ntg.reserve(n);
    for cell in grid.cells() {
        let (lo, hi) = {
            let (t, b) = (cell.top_depth(), cell.bottom_depth());
            (t.min(b), t.max(b))
        };
        // V2: `samples` is TVD-sorted (guaranteed by `with_logs` → `sort_by_tvd`),
        // so the in-range window is a contiguous slice found by two binary
        // searches instead of a full scan per cell.
        let start = samples.partition_point(|(tvd, _, _)| *tvd < lo);
        let end = samples.partition_point(|(tvd, _, _)| *tvd <= hi);
        let in_range = &samples[start..end];
        if in_range.is_empty() {
            poro.push(priors.porosity);
            sw.push(priors.water_saturation);
        } else {
            let phi: Vec<WeightedSample> = in_range
                .iter()
                .map(|(_, p, _)| WeightedSample::new(1.0, *p))
                .collect();
            poro.push(upscale_porosity(&phi)?);
            let sws: Vec<SwSample> = in_range
                .iter()
                .map(|(_, p, s)| SwSample {
                    length: 1.0,
                    porosity: *p,
                    water_saturation: *s,
                })
                .collect();
            sw.push(upscale_sw(&sws)?);
        }
        ntg.push(priors.net_to_gross);
    }
    let props = grid.properties_mut();
    props.set(Property {
        name: PORO.to_string(),
        values: poro,
    })?;
    props.set(Property {
        name: SW.to_string(),
        values: sw,
    })?;
    props.set(Property {
        name: NTG.to_string(),
        values: ntg,
    })?;
    Ok(())
}
