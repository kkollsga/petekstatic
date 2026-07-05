//! `tornado` — one-at-a-time sensitivity of oil in-place to each uncertain input
//! (`task_peteksim_tornado`, owned by petekStatic).
//!
//! A tornado holds every input at its P50 and swings **one** input between two
//! percentiles, re-realizing the model each time, so the output swing measures
//! that input's leverage. Ranked by swing, the bars form the familiar tornado
//! chart — the top bar is the input worth resolving first.
//!
//! ## Pivots at the *realized* percentiles (design note (2))
//! The pivot values are percentiles of the **realized input vectors** (the very
//! draws a matching [`run_structured_mc`](crate::model::run_structured_mc) would make
//! under the same `(n, seed)`), not of the analytic distribution. That is what
//! keeps the swing ranks consistent with the MC P-curve (a clamped/truncated
//! sampler's realized P10 differs from its nominal one). Call `tornado` with the
//! same `(n, seed)` as the MC run and the pivots reuse those exact draws.
//!
//! ## Output metric
//! Oil in-place \[Sm³\] (the two-contact oil leg, or the whole column for a
//! single-contact run) — the same primary metric as [`McResult`](crate::model::McResult).
//! A gas or GRV variant would swap the extracted output in [`realize_oil`]; oil
//! is shipped as the decision metric.

use crate::error::StaticError;
use crate::model::draw::RealizationDraw;
use crate::model::mc::{McInputs, RealizedInputs};
use crate::model::template::StaticModelTemplate;
use petektools::stats::percentile;
use std::cmp::Ordering;

/// One tornado bar: an input swung from its low to its high pivot, with the
/// resulting oil in-place at each and the absolute swing between them.
#[derive(Debug, Clone, PartialEq)]
pub struct TornadoBar {
    /// The input's name (e.g. `"net_to_gross"`).
    pub input: String,
    /// The low pivot value (input at `lo_pct`).
    pub lo_val: f64,
    /// The high pivot value (input at `hi_pct`).
    pub hi_val: f64,
    /// Oil in-place \[Sm³\] with the input at `lo_val`, others at P50.
    pub out_lo: f64,
    /// Oil in-place \[Sm³\] with the input at `hi_val`, others at P50.
    pub out_hi: f64,
    /// `|out_hi - out_lo|` — the sensitivity magnitude the chart ranks by.
    pub swing: f64,
}

/// A full scalar input set for one realization — the tornado's working point.
/// Every input at a fixed value; a pivot overrides exactly one field.
#[derive(Debug, Clone)]
struct Scalars {
    area: f64,
    gross: f64,
    contact: f64,
    goc: Option<f64>,
    porosity: f64,
    ntg: f64,
    sw: f64,
    sw_gas: Option<f64>,
    boi: f64,
    bgi: Option<f64>,
    shifts: Vec<(String, f64)>,
}

impl Scalars {
    fn to_draw(&self, seed_index: u64) -> RealizationDraw {
        let mut d = RealizationDraw::new(
            self.area,
            self.gross,
            self.contact,
            self.porosity,
            self.ntg,
            self.sw,
            seed_index,
        );
        if let Some(g) = self.goc {
            d = d.with_goc(g);
        }
        if let Some(s) = self.sw_gas {
            d = d.with_sw_gas(s);
        }
        for (name, delta) in &self.shifts {
            d = d.with_property_shift(name.clone(), *delta);
        }
        d
    }
}

/// `percentile` mapped into a [`StaticError`] (petekTools stats error → seam).
fn pct(data: &[f64], p: f64) -> Result<f64, StaticError> {
    Ok(percentile(data, p)?)
}

/// Build the P50 working point from the realized input vectors.
fn p50_scalars(r: &RealizedInputs) -> Result<Scalars, StaticError> {
    Ok(Scalars {
        area: pct(&r.area_m2, 50.0)?,
        gross: pct(&r.gross_height_m, 50.0)?,
        contact: pct(&r.contact_depth_m, 50.0)?,
        goc: r.goc_depth_m.as_ref().map(|v| pct(v, 50.0)).transpose()?,
        porosity: pct(&r.porosity, 50.0)?,
        ntg: pct(&r.net_to_gross, 50.0)?,
        sw: pct(&r.water_saturation, 50.0)?,
        sw_gas: r.sw_gas.as_ref().map(|v| pct(v, 50.0)).transpose()?,
        boi: pct(&r.boi, 50.0)?,
        bgi: r.bgi.as_ref().map(|v| pct(v, 50.0)).transpose()?,
        shifts: r
            .property_shifts
            .iter()
            .map(|(n, v)| Ok::<_, StaticError>((n.clone(), pct(v, 50.0)?)))
            .collect::<Result<_, _>>()?,
    })
}

/// The oil in-place \[Sm³\] of the model realized at `s` (the tornado metric).
/// Uses a **fixed** `seed_index` so a `Resimulate` property's pattern is held
/// constant across pivots — only the swung scalar moves the output.
fn realize_oil(
    tmpl: &mut StaticModelTemplate,
    s: &Scalars,
    seed_index: u64,
) -> Result<f64, StaticError> {
    let model = tmpl.realize(&s.to_draw(seed_index))?;
    let (oil, _gas, _grv) = crate::model::mc::outputs(&model, s.boi, s.bgi)?;
    Ok(oil)
}

/// Compute one bar: pivot `field` (via `set`) between `lo_val`/`hi_val`, others
/// held at the base P50.
fn bar(
    tmpl: &mut StaticModelTemplate,
    base: &Scalars,
    name: &str,
    lo_val: f64,
    hi_val: f64,
    seed_index: u64,
    set: impl Fn(&mut Scalars, f64),
) -> Result<TornadoBar, StaticError> {
    let mut lo_s = base.clone();
    set(&mut lo_s, lo_val);
    let mut hi_s = base.clone();
    set(&mut hi_s, hi_val);
    let out_lo = realize_oil(tmpl, &lo_s, seed_index)?;
    let out_hi = realize_oil(tmpl, &hi_s, seed_index)?;
    Ok(TornadoBar {
        input: name.to_string(),
        lo_val,
        hi_val,
        out_lo,
        out_hi,
        swing: (out_hi - out_lo).abs(),
    })
}

/// One-at-a-time tornado of oil in-place over the inputs of `inputs`, pivoting
/// each at the realized `lo_pct`/`hi_pct` percentiles (statistical percentiles in
/// `[0, 100]`, e.g. 10 and 90) with the others held at P50. Bars are returned
/// **pre-sorted by swing, descending** (the top driver first).
///
/// `n` / `seed` draw the realized input vectors the pivots come from — pass the
/// same values as the [`run_structured_mc`](crate::model::run_structured_mc) run to
/// reuse its exact draws. Deterministic: identical `(inputs, n, seed, lo_pct,
/// hi_pct)` give identical bars.
///
/// # Errors
/// [`StaticError`] if a pivot realization fails, `n == 0`, or a percentile is
/// out of range (propagated from petekTools stats).
pub fn tornado(
    tmpl: &mut StaticModelTemplate,
    inputs: &McInputs,
    n: usize,
    seed: u64,
    lo_pct: f64,
    hi_pct: f64,
) -> Result<Vec<TornadoBar>, StaticError> {
    let r = inputs.realize(n, seed)?;
    let base = p50_scalars(&r)?;
    // Hold Resimulate patterns fixed across every pivot realize.
    let sidx = seed;

    // One bar per uncertain input, in a fixed, table-driven order (the field table
    // is the single declaration of the scalar input set for the swing loop). Each
    // row is `(name, accessor, setter)`; the always-present rows always emit a bar,
    // the optional rows only when their realized vector is present. Order and
    // per-bar computation are identical to the hand-unrolled form.
    #[allow(clippy::type_complexity)]
    let required: &[(
        &str,
        fn(&RealizedInputs) -> &Vec<f64>,
        fn(&mut Scalars, f64),
    )] = &[
        ("area_m2", |r| &r.area_m2, |s, x| s.area = x),
        ("gross_height_m", |r| &r.gross_height_m, |s, x| s.gross = x),
        (
            "contact_depth_m",
            |r| &r.contact_depth_m,
            |s, x| s.contact = x,
        ),
        ("porosity", |r| &r.porosity, |s, x| s.porosity = x),
        ("net_to_gross", |r| &r.net_to_gross, |s, x| s.ntg = x),
        ("water_saturation", |r| &r.water_saturation, |s, x| s.sw = x),
        ("boi", |r| &r.boi, |s, x| s.boi = x),
    ];
    #[allow(clippy::type_complexity)]
    let optional: &[(
        &str,
        fn(&RealizedInputs) -> Option<&Vec<f64>>,
        fn(&mut Scalars, f64),
    )] = &[
        (
            "goc_depth_m",
            |r| r.goc_depth_m.as_ref(),
            |s, x| s.goc = Some(x),
        ),
        ("sw_gas", |r| r.sw_gas.as_ref(), |s, x| s.sw_gas = Some(x)),
        ("bgi", |r| r.bgi.as_ref(), |s, x| s.bgi = Some(x)),
    ];

    let mut bars = Vec::with_capacity(required.len() + optional.len() + r.property_shifts.len());
    for (name, get, set) in required {
        let v = get(&r);
        bars.push(bar(
            tmpl,
            &base,
            name,
            pct(v, lo_pct)?,
            pct(v, hi_pct)?,
            sidx,
            *set,
        )?);
    }
    for (name, get, set) in optional {
        if let Some(v) = get(&r) {
            bars.push(bar(
                tmpl,
                &base,
                name,
                pct(v, lo_pct)?,
                pct(v, hi_pct)?,
                sidx,
                *set,
            )?);
        }
    }
    // Per-property level shifts.
    for (name, v) in &r.property_shifts {
        let lo = pct(v, lo_pct)?;
        let hi = pct(v, hi_pct)?;
        let target = name.clone();
        bars.push(bar(tmpl, &base, name, lo, hi, sidx, move |s, x| {
            if let Some((_, d)) = s.shifts.iter_mut().find(|(nm, _)| nm == &target) {
                *d = x;
            }
        })?);
    }

    // Pre-sort by swing, descending (the tornado ordering).
    bars.sort_by(|a, b| b.swing.partial_cmp(&a.swing).unwrap_or(Ordering::Equal));
    Ok(bars)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gridder::{Conformity, SolveOpts};
    use crate::model::mc::Input;
    use crate::model::{BuildOpts, ConstantPriors};
    use crate::wireframe::{
        Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
    };
    use petektools::sampling::Sampler;

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

    fn tri(min: f64, mode: f64, max: f64) -> Input {
        Input::plain(Sampler::new_triangular(min, mode, max).unwrap())
    }

    /// A near-degenerate uniform (a pinned input with negligible spread).
    fn pinned(v: f64) -> Input {
        Input::plain(Sampler::new_uniform(v - 1e-4, v + 1e-4).unwrap())
    }

    #[test]
    fn dominant_input_ranks_first_and_bars_are_ordered() {
        // Porosity swings wide; everything else is pinned. Porosity must top the
        // chart, and the bars come pre-sorted by descending swing.
        let inputs = McInputs::new(
            pinned(100.0),
            pinned(50.0),
            pinned(5025.0),
            tri(0.05, 0.20, 0.35), // porosity — the wide one
            pinned(0.80),
            pinned(0.30),
            pinned(1.30),
        );
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let bars = tornado(&mut t, &inputs, 500, 42, 10.0, 90.0).unwrap();
        assert_eq!(
            bars[0].input, "porosity",
            "porosity must dominate: {bars:#?}"
        );
        for w in bars.windows(2) {
            assert!(w[0].swing >= w[1].swing, "bars not descending: {bars:#?}");
        }
        assert!(bars[0].swing > 0.0);
    }

    #[test]
    fn tornado_is_deterministic_under_seed() {
        let inputs = McInputs::new(
            tri(90.0, 100.0, 110.0),
            tri(45.0, 50.0, 55.0),
            tri(5020.0, 5025.0, 5030.0),
            tri(0.18, 0.22, 0.26),
            tri(0.70, 0.80, 0.90),
            tri(0.25, 0.30, 0.35),
            tri(1.20, 1.30, 1.45),
        );
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t1 = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut t2 = StaticModelTemplate::new(&wf, opts()).unwrap();
        let a = tornado(&mut t1, &inputs, 300, 7, 10.0, 90.0).unwrap();
        let b = tornado(&mut t2, &inputs, 300, 7, 10.0, 90.0).unwrap();
        assert_eq!(a, b, "tornado not deterministic under a fixed seed");
    }

    #[test]
    fn pivot_values_match_realized_percentiles() {
        let inputs = McInputs::new(
            tri(90.0, 100.0, 110.0),
            tri(45.0, 50.0, 55.0),
            tri(5020.0, 5025.0, 5030.0),
            tri(0.18, 0.22, 0.26),
            tri(0.70, 0.80, 0.90),
            tri(0.25, 0.30, 0.35),
            tri(1.20, 1.30, 1.45),
        );
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let (n, seed, lo, hi) = (400usize, 99u64, 15.0, 85.0);
        let bars = tornado(&mut t, &inputs, n, seed, lo, hi).unwrap();
        // Re-derive the realized vectors identically and check one bar's pivots.
        let realized = inputs.realize(n, seed).unwrap();
        let area_bar = bars.iter().find(|b| b.input == "area_m2").unwrap();
        assert!((area_bar.lo_val - percentile(&realized.area_m2, lo).unwrap()).abs() < 1e-12);
        assert!((area_bar.hi_val - percentile(&realized.area_m2, hi).unwrap()).abs() < 1e-12);
        let ntg_bar = bars.iter().find(|b| b.input == "net_to_gross").unwrap();
        assert!((ntg_bar.lo_val - percentile(&realized.net_to_gross, lo).unwrap()).abs() < 1e-12);
    }
}
