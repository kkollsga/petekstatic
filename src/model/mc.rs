//! `mc` — the structured Monte-Carlo driver over the [`StaticModelTemplate`]
//! regeneration seam (`task_peteksim_mc_structured`, owned by petekStatic per
//! graph `decision_layer_charters`).
//!
//! The template makes one *realization* cheap; this module makes a *run* cheap:
//! it samples every uncertain input from a [petekTools `Sampler`](petektools::sampling::Sampler),
//! feeds each draw through [`StaticModelTemplate::realize`], and keeps the full
//! per-draw output vectors (the **W17 lesson**: never throw the realizations
//! away — P-curves, tornado reuse and cross-checks all need them).
//!
//! ## Boi at the output surface, not in the draw (ratified)
//! The oil/gas formation-volume factors (`boi`/`bgi`) are sampled like any other
//! uncertain input but applied **at the volumetrics conversion** (`HCPV / Boi`),
//! never inside the geometry/property draw — PVT never rides the
//! [`RealizationDraw`] (graph `decision_staticmodel_regen_seam`). Keeping them in
//! [`McInputs`] spares the caller a parallel sampler.
//!
//! ## Error policy — fail-fast with the failing draw index
//! A draw that fails validation (H2: out-of-range priors, a crossed base, an
//! off-lattice structural shift) is a **typed** error to the caller:
//! [`run_structured_mc`] stops at the first bad draw and returns
//! [`StaticError::McDraw`] carrying that draw's index and the underlying cause
//! (`source()` reaches the origin). This is deliberate over collect-and-report:
//! a bad draw means the *input distribution* strays outside the physical range,
//! which the caller should fix at the sampler (clamp it — see [`Input::clamped`])
//! rather than silently average over. Clamp your fraction samplers and no valid
//! configuration ever trips this.

use crate::error::StaticError;
use crate::model::draw::RealizationDraw;
use crate::model::model::StaticModel;
use crate::model::template::StaticModelTemplate;
use crate::spill::{spill_grid_to, unique_spill_path};
use crate::volumetrics::{GasFvf, OilFvf};
use petektools::sampling::{
    aggregate, reservoir_summary, seeded_rng, Clamped, Correlation, ReservoirSummary, Sampler,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// The settings of a structured MC run — the ONE run-resources value consolidating
/// the four historical entry points (`run_structured_mc` / `_parallel` / `_spilled`
/// / `_parallel_spilled`, now deprecated thin wrappers) behind [`run_mc`]
/// (`task_petekstatic_spec_mirror`; suite ruling: run resources are a settings
/// value). Serializable (a run recipe is a savable file), value-comparable.
///
/// - `workers == 1` (the [`McSettings::new`] default) = the serial driver;
///   `> 1` = the rayon-sharded driver (one template clone per shard; see the
///   sharding determinism contract on [`run_mc`]). Clamped to `[1, n]`.
/// - `spill_dir: None` (default) = the in-core mode; `Some(dir)` = the spilled
///   (out-of-core) mode with each shard's reused f32 store written under `dir`
///   (rulings R3/R4). *The directory selects the mode* — pass
///   [`std::env::temp_dir()`] for the historical "spilled at the platform temp
///   dir" behaviour.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
#[serde(default)]
pub struct McSettings {
    /// Number of draws (`> 0`).
    pub n: usize,
    /// The run seed — same `(inputs, n, seed, workers)` → bit-identical vectors.
    pub seed: u64,
    /// Worker count; `1` = serial. See [`default_mc_workers`] for the recommended
    /// parallel sizing (the loop is memory-bandwidth-bound).
    pub workers: usize,
    /// `Some(dir)` = spilled (out-of-core) mode, stores under `dir`; `None` = in-core.
    pub spill_dir: Option<PathBuf>,
}

impl Default for McSettings {
    fn default() -> Self {
        Self {
            n: 1,
            seed: 0,
            workers: 1,
            spill_dir: None,
        }
    }
}

impl McSettings {
    /// A serial, in-core run of `n` draws at `seed` — add the optionals with the
    /// `with_*` sugar.
    #[must_use]
    pub fn new(n: usize, seed: u64) -> Self {
        Self {
            n,
            seed,
            ..Self::default()
        }
    }

    /// Shard across `workers` template clones ([`default_mc_workers`] for the
    /// recommended sizing).
    #[must_use]
    pub fn with_workers(mut self, workers: usize) -> Self {
        self.workers = workers;
        self
    }

    /// Run in the spilled (out-of-core) mode, each shard's reused f32 store under
    /// `dir` (rulings R3/R4).
    #[must_use]
    pub fn with_spill_dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.spill_dir = Some(dir.into());
        self
    }
}

/// Run a seeded, reproducible structured MC per [`McSettings`] — THE single MC
/// entry (`task_petekstatic_spec_mirror`): sample every input in [`McInputs`],
/// realize each draw through `tmpl` (serial, or sharded across
/// `settings.workers` template clones), in-core or spilled per
/// `settings.spill_dir`, and keep the oil/gas/GRV vectors (W17).
///
/// ## Determinism (the consolidated contract)
/// - Same `(inputs, settings)` → bit-identical vectors, run to run — one seeded
///   stream, fixed field order, `seed_index = seed + i` per draw, deterministic
///   shard boundaries, recombination by draw index.
/// - Across worker counts: the same draw multiset; identical vectors in the
///   common case (no per-draw structural shifts) — the sharded == serial
///   contract pinned by the driver tests.
/// - The spilled mode is bit-deterministic within itself (f32 quantization is a
///   deterministic per-cell function of the draw, shard-split-independent).
///
/// The template's [`crate::model::McMode::LevelShift`] pattern cache is pre-warmed once
/// (draw 0 on `tmpl`) before any shards clone it.
///
/// # Errors
/// [`StaticError::McDraw`] (carrying the failing draw index) on the first draw
/// whose realization or volumetrics fails (fail-fast, see the [module
/// docs](self)); [`StaticError::InvalidInput`] if `settings.n == 0`.
pub fn run_mc(
    tmpl: &mut StaticModelTemplate,
    inputs: &McInputs,
    settings: &McSettings,
) -> Result<McResult, StaticError> {
    let realized = inputs.realize(settings.n, settings.seed)?;
    let n = settings.n;
    let dir = settings.spill_dir.as_deref();
    let workers = settings.workers.clamp(1, n);
    if workers == 1 {
        let triples = run_draws_impl(tmpl, &realized, 0..n, dir)?;
        return Ok(assemble(triples, realized));
    }
    // Pre-warm the LevelShift pattern cache on the shared template (draw 0) so the
    // per-worker clones inherit it — otherwise each of `workers` clones would
    // re-propagate the field (up to `workers`× the one-time propagate cost, which
    // dominates a run). Its output is recomputed cheaply (cached) inside its shard.
    let _ = run_draws_impl(tmpl, &realized, 0..1, dir)?;

    let ranges = shard_ranges(n, workers);
    let realized_ref = &realized;
    let template_ref = &*tmpl;
    // Each shard is one heavy sequential task; with `workers` shards, rayon runs at
    // most `workers` in parallel (the intended 4–6-way width) on the global pool.
    let shards: Result<Vec<Vec<OutTriple>>, StaticError> = ranges
        .into_par_iter()
        .map(|r| {
            let mut t = template_ref.clone();
            run_draws_impl(&mut t, realized_ref, r, dir)
        })
        .collect();

    // Contiguous shards recombine to draw-index order by concatenation.
    let triples: Vec<OutTriple> = shards?.into_iter().flatten().collect();
    Ok(assemble(triples, realized))
}

/// A per-draw output triple `(oil Sm³, gas Sm³, GRV m³)` — the sharded driver's
/// unit of work, recombined into the [`McResult`] vectors.
type OutTriple = (f64, f64, f64);

/// One uncertain input's distribution: a raw petekTools [`Sampler`] or a
/// [`Clamped`] one (its draws snapped into `[lo, hi]`).
///
/// **Clamp fraction inputs.** φ / NTG / Sw must land in `[0, 1]` or
/// [`StaticModelTemplate::realize`] rejects the draw (H2). A `Normal` prior can
/// stray outside, so wrap it with [`Input::clamped`] — this uses petekTools'
/// own [`Clamped`] (a hard limiter), *not* a bespoke wrapper. When a quantity is
/// genuinely bounded (a saturation), prefer a
/// [`Sampler::new_truncated_normal`] which reshapes the density instead of
/// piling mass on the bounds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Input {
    /// A raw sampler — draws are used as-is.
    Plain(Sampler),
    /// A clamped sampler — every draw is snapped into `[lo, hi]`.
    Clamped(Clamped),
}

impl Input {
    /// A raw (unclamped) input from `sampler`.
    #[must_use]
    pub fn plain(sampler: Sampler) -> Self {
        Input::Plain(sampler)
    }

    /// A clamped input: `sampler`'s draws snapped into `[lo, hi]` via petekTools'
    /// [`Clamped`].
    ///
    /// # Errors
    /// [`StaticError::Algo`] (from [`Sampler::clamped`]) unless `lo < hi` and both
    /// are finite.
    pub fn clamped(sampler: Sampler, lo: f64, hi: f64) -> Result<Self, StaticError> {
        Ok(Input::Clamped(sampler.clamped(lo, hi)?))
    }

    /// Draw `n` samples from `rng` (dispatches to the inner sampler).
    fn sample_n<R: rand::Rng>(&self, n: usize, rng: &mut R) -> Vec<f64> {
        match self {
            Input::Plain(s) => s.sample_n(n, rng),
            Input::Clamped(c) => c.sample_n(n, rng),
        }
    }
}

impl From<Sampler> for Input {
    fn from(s: Sampler) -> Self {
        Input::Plain(s)
    }
}

impl From<Clamped> for Input {
    fn from(c: Clamped) -> Self {
        Input::Clamped(c)
    }
}

/// The uncertain inputs of a structured MC run — one [`Input`] per quantity.
///
/// The seven load-bearing scalars are always present; `goc_depth_m` /`sw_gas`
/// (two-contact gas cap) and `bgi` (gas FVF) are optional, and
/// `property_shifts` supplies the per-draw additive level shift for each
/// [`crate::model::McMode::LevelShift`] property attached to the template.
/// [`crate::model::McMode::Resimulate`] properties need **no** entry here — the
/// template redraws their pattern from the draw's `seed_index`.
#[derive(Debug, Clone)]
pub struct McInputs {
    /// Areal footprint \[m²\].
    pub area_m2: Input,
    /// Gross column thickness / base-relief level \[m\].
    pub gross_height_m: Input,
    /// The (lower, OWC/FWL) hydrocarbon contact depth \[m\].
    pub contact_depth_m: Input,
    /// Optional gas–oil contact depth \[m\] — `Some` makes every draw a
    /// two-contact (gas-cap + oil-rim) column.
    pub goc_depth_m: Option<Input>,
    /// Porosity prior (fraction; clamp to `[0, 1]`).
    pub porosity: Input,
    /// Net-to-gross prior (fraction; clamp to `[0, 1]`).
    pub net_to_gross: Input,
    /// Water-saturation prior (fraction; clamp to `[0, 1]`).
    pub water_saturation: Input,
    /// Optional gas-cap connate-water override (fraction) for two-contact draws.
    pub sw_gas: Option<Input>,
    /// Oil formation-volume factor \[Rm³/Sm³\], applied at the output surface.
    pub boi: Input,
    /// Optional gas formation-volume factor \[Rm³/Sm³\] for the gas-cap OGIP of a
    /// two-contact draw; `None` leaves `gas_sm3` at zero.
    pub bgi: Option<Input>,
    /// Per-property additive level shifts `(property, Input)` for
    /// [`crate::model::McMode::LevelShift`] properties on the template.
    pub property_shifts: Vec<(String, Input)>,
}

impl McInputs {
    /// A bare inputs set — the seven load-bearing scalars, no gas cap / FVF-gas /
    /// property shifts. Add the optionals with the `with_*` builders.
    #[must_use]
    pub fn new(
        area_m2: Input,
        gross_height_m: Input,
        contact_depth_m: Input,
        porosity: Input,
        net_to_gross: Input,
        water_saturation: Input,
        boi: Input,
    ) -> Self {
        Self {
            area_m2,
            gross_height_m,
            contact_depth_m,
            goc_depth_m: None,
            porosity,
            net_to_gross,
            water_saturation,
            sw_gas: None,
            boi,
            bgi: None,
            property_shifts: Vec::new(),
        }
    }

    /// Add a gas–oil contact input (every draw becomes a two-contact column).
    #[must_use]
    pub fn with_goc(mut self, goc_depth_m: Input) -> Self {
        self.goc_depth_m = Some(goc_depth_m);
        self
    }

    /// Add a gas-cap connate-water override input (two-contact draws).
    #[must_use]
    pub fn with_sw_gas(mut self, sw_gas: Input) -> Self {
        self.sw_gas = Some(sw_gas);
        self
    }

    /// Add a gas formation-volume-factor input so two-contact draws report
    /// gas-cap OGIP in `gas_sm3`.
    #[must_use]
    pub fn with_bgi(mut self, bgi: Input) -> Self {
        self.bgi = Some(bgi);
        self
    }

    /// Add a per-draw level-shift input for a [`crate::model::McMode::LevelShift`]
    /// property. Replaces any prior shift for the same property name.
    #[must_use]
    pub fn with_property_shift(mut self, property: impl Into<String>, shift: Input) -> Self {
        let property = property.into();
        self.property_shifts.retain(|(p, _)| p != &property);
        self.property_shifts.push((property, shift));
        self
    }

    /// Draw the `n` realized input vectors from a single seeded stream (fixed
    /// field order → bit-reproducible for a given `seed`). Each field is an
    /// independent draw batch; the vectors are retained on the [`McResult`] and
    /// re-derived identically by [`crate::model::tornado`] under the same `(n, seed)`.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if `n == 0`.
    pub fn realize(&self, n: usize, seed: u64) -> Result<RealizedInputs, StaticError> {
        if n == 0 {
            return Err(StaticError::InvalidInput(
                "structured MC needs at least one draw (n = 0)".into(),
            ));
        }
        let mut rng = seeded_rng(seed);
        // Fixed order — the sole guarantor of reproducibility across the stream.
        let area_m2 = self.area_m2.sample_n(n, &mut rng);
        let gross_height_m = self.gross_height_m.sample_n(n, &mut rng);
        let contact_depth_m = self.contact_depth_m.sample_n(n, &mut rng);
        let goc_depth_m = self.goc_depth_m.as_ref().map(|s| s.sample_n(n, &mut rng));
        let porosity = self.porosity.sample_n(n, &mut rng);
        let net_to_gross = self.net_to_gross.sample_n(n, &mut rng);
        let water_saturation = self.water_saturation.sample_n(n, &mut rng);
        let sw_gas = self.sw_gas.as_ref().map(|s| s.sample_n(n, &mut rng));
        let boi = self.boi.sample_n(n, &mut rng);
        let bgi = self.bgi.as_ref().map(|s| s.sample_n(n, &mut rng));
        let property_shifts = self
            .property_shifts
            .iter()
            .map(|(name, s)| (name.clone(), s.sample_n(n, &mut rng)))
            .collect();
        Ok(RealizedInputs {
            area_m2,
            gross_height_m,
            contact_depth_m,
            goc_depth_m,
            porosity,
            net_to_gross,
            water_saturation,
            sw_gas,
            boi,
            bgi,
            property_shifts,
            seed,
        })
    }
}

/// The realized (sampled) input vectors of a run — one `Vec<f64>` of length `n`
/// per field. Retained on [`McResult`] so [`crate::model::tornado`] can pivot at the
/// **realized** percentiles (design note (2): realized, not analytic, percentiles
/// keep the swing ranks consistent with the MC).
#[derive(Debug, Clone)]
pub struct RealizedInputs {
    /// Areal footprint draws \[m²\].
    pub area_m2: Vec<f64>,
    /// Gross-thickness draws \[m\].
    pub gross_height_m: Vec<f64>,
    /// Contact-depth (OWC/FWL) draws \[m\].
    pub contact_depth_m: Vec<f64>,
    /// Gas–oil-contact draws \[m\] (two-contact runs only).
    pub goc_depth_m: Option<Vec<f64>>,
    /// Porosity draws (fraction).
    pub porosity: Vec<f64>,
    /// Net-to-gross draws (fraction).
    pub net_to_gross: Vec<f64>,
    /// Water-saturation draws (fraction).
    pub water_saturation: Vec<f64>,
    /// Gas-cap connate-water draws (fraction; two-contact runs only).
    pub sw_gas: Option<Vec<f64>>,
    /// Oil-FVF draws \[Rm³/Sm³\].
    pub boi: Vec<f64>,
    /// Gas-FVF draws \[Rm³/Sm³\] (two-contact runs only).
    pub bgi: Option<Vec<f64>>,
    /// Per-property level-shift draws `(property, deltas)`.
    pub property_shifts: Vec<(String, Vec<f64>)>,
    seed: u64,
}

impl RealizedInputs {
    /// Number of draws.
    #[must_use]
    pub fn len(&self) -> usize {
        self.area_m2.len()
    }

    /// Whether the input set is empty (never true after a successful `realize`).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.area_m2.is_empty()
    }

    /// The [`RealizationDraw`] for draw `i` (geometry + priors; the draw's
    /// `seed_index = seed + i` gives every draw a reproducible, distinct pattern
    /// seed for [`crate::model::McMode::Resimulate`] properties).
    fn draw_at(&self, i: usize) -> RealizationDraw {
        let mut d = RealizationDraw::new(
            self.area_m2[i],
            self.gross_height_m[i],
            self.contact_depth_m[i],
            self.porosity[i],
            self.net_to_gross[i],
            self.water_saturation[i],
            self.seed.wrapping_add(i as u64),
        );
        if let Some(g) = &self.goc_depth_m {
            d = d.with_goc(g[i]);
        }
        if let Some(s) = &self.sw_gas {
            d = d.with_sw_gas(s[i]);
        }
        for (name, deltas) in &self.property_shifts {
            d = d.with_property_shift(name.clone(), deltas[i]);
        }
        d
    }
}

/// The kept result of a structured MC run — the **per-draw realization vectors**
/// (oil Sm³, gas Sm³, GRV m³), plus the realized inputs for tornado reuse.
///
/// The vectors are the product (W17): summarise them via [`McResult::summary`]
/// (P90/P50/P10 in the oil-industry convention), aggregate several segments via
/// [`aggregate_field`], or feed them to [`crate::model::tornado`].
#[derive(Debug, Clone)]
pub struct McResult {
    /// Oil in-place per draw \[Sm³\] — the two-contact oil leg, or the whole
    /// column for a single-contact run. The **primary output metric**.
    pub oil_sm3: Vec<f64>,
    /// Free gas in-place per draw \[Sm³\] — the gas cap of a two-contact run
    /// (needs `bgi`); zero for a single-contact run.
    pub gas_sm3: Vec<f64>,
    /// Gross rock volume of the hydrocarbon column per draw \[m³\].
    pub grv_m3: Vec<f64>,
    realized: RealizedInputs,
}

impl McResult {
    /// Number of draws (`oil_sm3.len()`).
    #[must_use]
    pub fn len(&self) -> usize {
        self.oil_sm3.len()
    }

    /// Whether the result is empty (never after a successful run).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.oil_sm3.is_empty()
    }

    /// The realized input vectors (for tornado reuse / auditing).
    #[must_use]
    pub fn realized_inputs(&self) -> &RealizedInputs {
        &self.realized
    }

    /// P90/P50/P10 + mean of the **oil** realizations (the primary metric), in the
    /// oil-industry exceedance convention (`p90 ≤ p50 ≤ p10`) via petekTools'
    /// [`reservoir_summary`].
    ///
    /// # Errors
    /// [`StaticError::Algo`] on an empty realization set.
    pub fn summary(&self) -> Result<ReservoirSummary, StaticError> {
        Ok(reservoir_summary(&self.oil_sm3)?)
    }

    /// P-curve summary of the **gas** realizations \[Sm³\].
    ///
    /// # Errors
    /// [`StaticError::Algo`] on an empty realization set.
    pub fn gas_summary(&self) -> Result<ReservoirSummary, StaticError> {
        Ok(reservoir_summary(&self.gas_sm3)?)
    }

    /// P-curve summary of the **GRV** realizations \[m³\].
    ///
    /// # Errors
    /// [`StaticError::Algo`] on an empty realization set.
    pub fn grv_summary(&self) -> Result<ReservoirSummary, StaticError> {
        Ok(reservoir_summary(&self.grv_m3)?)
    }
}

/// Run a seeded, reproducible structured MC: sample every input in [`McInputs`],
/// realize each draw through `tmpl`, and keep the oil/gas/GRV vectors.
///
/// Reproducible: the same `(inputs, n, seed)` gives bit-identical vectors (one
/// seeded stream, fixed field order, `seed_index = seed + i` per draw).
///
/// # Errors
/// [`StaticError::McDraw`] (carrying the failing draw index) on the first draw
/// whose realization or volumetrics fails — the fail-fast policy (see the
/// [module docs](self)). [`StaticError::InvalidInput`] if `n == 0`.
#[deprecated(since = "0.1.0", note = "use `run_mc` with `McSettings::new(n, seed)`")]
pub fn run_structured_mc(
    tmpl: &mut StaticModelTemplate,
    inputs: &McInputs,
    n: usize,
    seed: u64,
) -> Result<McResult, StaticError> {
    run_mc(tmpl, inputs, &McSettings::new(n, seed))
}

/// The recommended default worker count for [`run_structured_mc_parallel`]:
/// `min(6, available_parallelism)`. The MC loop is **memory-bandwidth-bound**
/// (~90 MB churned per draw at 1M cells), so parallel scaling saturates well
/// before core count — the P9 baseline measured 4.78× at 8 threads and ~2.7× at
/// 4, so 4–6 workers capture essentially all the win while leaving headroom.
#[must_use]
pub fn default_mc_workers() -> usize {
    let cores = std::thread::available_parallelism()
        .map(std::num::NonZeroUsize::get)
        .unwrap_or(1);
    cores.clamp(1, 6)
}

/// Rayon-sharded [`run_structured_mc`]: split the `n` draws across `workers`
/// shards, each realizing its contiguous slice on its **own template clone**
/// (`StaticModelTemplate: Send`; the warm-start chain is serial *within* a shard),
/// then recombine the per-draw vectors **in draw-index order**. `workers` is
/// clamped to `[1, n]`; pass [`default_mc_workers`] for the recommended sizing.
///
/// ## Determinism (the sharding contract)
/// - **Same `(inputs, n, seed, workers)` → bit-identical vectors, run to run.**
///   Every draw's inputs come from the *single* seeded stream
///   [`McInputs::realize`] (identical to the serial path), shard boundaries are a
///   deterministic function of `(n, workers)`, and recombination is by draw index
///   — no scheduling nondeterminism reaches the result.
/// - **Across *different* `workers` → the same draw multiset.** A draw's inputs
///   and its per-draw pattern seed (`seed_index = seed + i`, worker-independent)
///   do not depend on the shard split, and the structural solve converges to the
///   control-determined surface regardless of its warm-start guess — so in the
///   common case (no *per-draw* structural shifts) the vectors are identical to
///   the serial run. Where per-draw `structural` shifts make the warm-start chain
///   path-dependent, a draw that lands at a shard boundary may differ from the
///   serial value **within the solver's convergence tolerance** (a different warm
///   start, same fixed point) — the draw *set* is unchanged. Summaries (P90/P50/P10)
///   are invariant either way.
///
/// The template's [`McMode::LevelShift`] pattern cache is **pre-warmed once**
/// (draw 0 on `tmpl`) before the shards clone it, so the one-time field
/// propagation is paid once rather than once per worker.
///
/// # Errors
/// [`StaticError::McDraw`] (carrying the failing draw index) on the first draw
/// whose realization or volumetrics fails; [`StaticError::InvalidInput`] if
/// `n == 0`.
#[deprecated(
    since = "0.1.0",
    note = "use `run_mc` with `McSettings::new(n, seed).with_workers(workers)`"
)]
pub fn run_structured_mc_parallel(
    tmpl: &mut StaticModelTemplate,
    inputs: &McInputs,
    n: usize,
    seed: u64,
    workers: usize,
) -> Result<McResult, StaticError> {
    run_mc(
        tmpl,
        inputs,
        &McSettings::new(n, seed).with_workers(workers),
    )
}

/// Spilled (out-of-core) structured MC (ruling R3/R4): like [`run_structured_mc`],
/// but each draw's summary streams from an f32 spill store instead of the in-core
/// cubes — the mode for a model whose per-draw live set exceeds the memory budget.
/// The reusable model is realized in-core then flushed to **one reused per-shard
/// store** (overwritten per draw, never a new file per draw); the P90/P50/P10 are
/// tolerance-equal to the in-core run (f32 storage, R4 honesty clause).
///
/// `spill_dir` is where the reused store is written (`None` = the platform temp
/// dir). Bit-deterministic **within** the spilled mode: same `(inputs, n, seed)` →
/// identical vectors.
///
/// # Errors
/// As [`run_structured_mc`]; plus [`StaticError::Algo`] on a store I/O failure.
#[deprecated(
    since = "0.1.0",
    note = "use `run_mc` with `McSettings::new(n, seed).with_spill_dir(dir)` \
            (pass `std::env::temp_dir()` for the old `None` behaviour)"
)]
pub fn run_structured_mc_spilled(
    tmpl: &mut StaticModelTemplate,
    inputs: &McInputs,
    n: usize,
    seed: u64,
    spill_dir: Option<PathBuf>,
) -> Result<McResult, StaticError> {
    let dir = spill_dir.unwrap_or_else(std::env::temp_dir);
    run_mc(tmpl, inputs, &McSettings::new(n, seed).with_spill_dir(dir))
}

/// Rayon-sharded spilled structured MC — the parallel sibling of
/// [`run_structured_mc_spilled`] (each shard its own template clone + its own reused
/// store). Proves the sharding determinism contract holds in the spilled mode too:
/// **sharded == serial at every worker count** (the f32 quantization is a
/// deterministic per-cell function of the draw, independent of the shard split).
///
/// # Errors
/// As [`run_structured_mc_parallel`]; plus [`StaticError::Algo`] on store I/O.
#[deprecated(
    since = "0.1.0",
    note = "use `run_mc` with `McSettings::new(n, seed).with_workers(workers).with_spill_dir(dir)`"
)]
pub fn run_structured_mc_parallel_spilled(
    tmpl: &mut StaticModelTemplate,
    inputs: &McInputs,
    n: usize,
    seed: u64,
    workers: usize,
    spill_dir: Option<PathBuf>,
) -> Result<McResult, StaticError> {
    let dir = spill_dir.unwrap_or_else(std::env::temp_dir);
    run_mc(
        tmpl,
        inputs,
        &McSettings::new(n, seed)
            .with_workers(workers)
            .with_spill_dir(dir),
    )
}

/// Even contiguous split of `0..n` into `workers` ranges, remainder on the first
/// shards — a deterministic function of `(n, workers)` (the determinism contract).
fn shard_ranges(n: usize, workers: usize) -> Vec<std::ops::Range<usize>> {
    let base = n / workers;
    let rem = n % workers;
    let mut ranges = Vec::with_capacity(workers);
    let mut start = 0;
    for w in 0..workers {
        let len = base + usize::from(w < rem);
        ranges.push(start..start + len);
        start += len;
    }
    ranges
}

/// The per-draw body of the serial and sharded drivers, in-core (`spill_dir` None)
/// or **spilled** (`Some(dir)`, ruling R3). `range` runs on `tmpl`'s warm-start
/// chain (serial within the range).
///
/// **In-core (`None`):** one reusable model per shard — `realize_into` recycles
/// its ZCORN + cube buffers in place, so the ~100 MB/draw churn that caps parallel
/// scaling is paid once per shard, and `in_place_summary` consumes it without
/// materializing more. **MC never spills per-draw state** (R3): one model,
/// realize → summary → discard, nothing written to disk.
///
/// **Spilled (`Some(dir)`):** the reusable model is realized in-core, then each
/// draw is flushed to **one per-shard store** (a unique path, *overwritten* every
/// draw — not a new file per draw, R3) and its summary streamed from the f32 lanes
/// (the R4 bandwidth lever). The store is attached to the reused model only for the
/// summary (the in-core grid stays intact for the next `realize_into`), then
/// detached. Bit-deterministic within the mode (f32 quantization is deterministic),
/// so sharded == serial at every worker count.
fn run_draws_impl(
    tmpl: &mut StaticModelTemplate,
    realized: &RealizedInputs,
    range: std::ops::Range<usize>,
    spill_dir: Option<&Path>,
) -> Result<Vec<OutTriple>, StaticError> {
    let mut out = Vec::with_capacity(range.len());
    let mut model = tmpl.reusable_model();
    // One reused per-shard store path for the spilled mode (overwritten per draw).
    let shard_store: Option<PathBuf> = spill_dir.map(unique_spill_path);
    for i in range {
        let draw = realized.draw_at(i);
        tmpl.realize_into(&draw, &mut model)
            .map_err(|e| at_draw(i, e))?;
        let boi = realized.boi[i];
        let bgi = realized.bgi.as_ref().map(|b| b[i]);
        let triple = match &shard_store {
            None => outputs(&model, boi, bgi).map_err(|e| at_draw(i, e))?,
            Some(path) => {
                // Flush the just-realized in-core grid to the reused store, attach
                // it so `in_place_summary` streams the f32 lanes, then detach.
                let backing =
                    spill_grid_to(model.grid(), path, false).map_err(|e| at_draw(i, e))?;
                model.set_spill(Some(Arc::new(backing)));
                let triple = outputs(&model, boi, bgi).map_err(|e| at_draw(i, e));
                model.set_spill(None);
                triple?
            }
        };
        out.push(triple);
    }
    if let Some(path) = &shard_store {
        let _ = std::fs::remove_file(path); // reused store is transient; clean at shard end
    }
    Ok(out)
}

/// Assemble the per-draw `(oil, gas, grv)` triples (in draw-index order) into an
/// [`McResult`] with the realized inputs retained.
fn assemble(triples: Vec<OutTriple>, realized: RealizedInputs) -> McResult {
    let n = triples.len();
    let mut oil_sm3 = Vec::with_capacity(n);
    let mut gas_sm3 = Vec::with_capacity(n);
    let mut grv_m3 = Vec::with_capacity(n);
    for (oil, gas, grv) in triples {
        oil_sm3.push(oil);
        gas_sm3.push(gas);
        grv_m3.push(grv);
    }
    McResult {
        oil_sm3,
        gas_sm3,
        grv_m3,
        realized,
    }
}

/// Aggregate several segments' **oil** realizations into a field total under an
/// explicit dependence assumption, delegating to petekTools' [`aggregate`].
///
/// ## The Independent / Comonotonic bracketing pattern (documented once)
/// The library offers no inter-segment correlation model, so a field total is
/// **bracketed** between two bounding cases (both share the same mean):
/// - [`Correlation::Independent`] — segments sampled from their own streams sum
///   index-wise; the *narrow* bound (partial cancellation of independent draws).
/// - [`Correlation::Comonotonic`] — perfect positive rank dependence (sorted,
///   summed rank-for-rank); the *wide* downside bound (everything low together /
///   high together), so its P90 is lower and its P10 higher than the independent
///   case. Carry Independent as the base and quote Comonotonic as the downside-P90
///   stress. The truth sits between (shared seismic/depth-conversion/analogue
///   couplings pull toward Comonotonic; local geometry toward Independent).
///
/// The result length is the shortest segment's; an empty input gives an empty
/// `Vec` (petekTools semantics).
#[must_use]
pub fn aggregate_field(segments: &[&McResult], corr: Correlation) -> Vec<f64> {
    let vecs: Vec<&[f64]> = segments.iter().map(|s| s.oil_sm3.as_slice()).collect();
    aggregate(&vecs, corr)
}

/// The per-draw output triple `(oil Sm³, gas Sm³, GRV m³)` from a realized model.
/// Two-contact: oil = the oil leg, gas = the gas cap (needs `bgi`). Single
/// contact: oil = the whole column, gas = 0.
///
/// Uses [`StaticModel::in_place_summary`] (V7): only the GRV/HCPV/OOIP/OGIP
/// aggregates feed the P-curve, so the per-cell HCPV cube is never materialized —
/// avoiding a `cell_count`-length alloc-and-discard on every MC draw / tornado
/// pivot.
pub(crate) fn outputs(
    model: &StaticModel,
    boi: f64,
    bgi: Option<f64>,
) -> Result<(f64, f64, f64), StaticError> {
    let ip = model.in_place_summary()?;
    let two_contact = ip.gas.is_some() && ip.oil.is_some();
    let oil = if two_contact {
        ip.oil_zone_ooip_sm3(OilFvf::new(boi)?)
    } else {
        ip.ooip_sm3(OilFvf::new(boi)?)
    };
    let gas = match (two_contact, bgi) {
        (true, Some(bgi)) => ip.gas_zone_ogip_sm3(GasFvf::new(bgi)?),
        _ => 0.0,
    };
    Ok((oil, gas, ip.grv_m3))
}

/// Wrap a per-draw failure with its draw index (fail-fast typed error).
fn at_draw(index: usize, source: StaticError) -> StaticError {
    StaticError::McDraw {
        index,
        source: Box::new(source),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gridder::{Conformity, SolveOpts};
    use crate::model::{
        BuildOpts, ConstantPriors, Gaussian, McMode, PropertyPipeline, UpscaleMethod, WellLog,
    };
    use crate::wireframe::{
        Boundary, Contact, ContactKind, GriddedDepth, Hardness, Horizon, HorizonRole, Wireframe,
    };

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

    /// A nominal single-contact inputs set with modest spreads (all fractions
    /// safely inside `[0,1]` so no clamp is needed).
    fn nominal_inputs() -> McInputs {
        McInputs::new(
            tri(90.0, 100.0, 110.0),     // area
            tri(45.0, 50.0, 55.0),       // gross
            tri(5020.0, 5025.0, 5030.0), // contact
            tri(0.18, 0.22, 0.26),       // porosity
            tri(0.70, 0.80, 0.90),       // ntg
            tri(0.25, 0.30, 0.35),       // sw
            tri(1.20, 1.30, 1.45),       // boi
        )
    }

    #[test]
    fn run_is_bit_reproducible() {
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let inputs = nominal_inputs();
        let mut t1 = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut t2 = StaticModelTemplate::new(&wf, opts()).unwrap();
        let a = run_mc(&mut t1, &inputs, &McSettings::new(200, 42)).unwrap();
        let b = run_mc(&mut t2, &inputs, &McSettings::new(200, 42)).unwrap();
        assert_eq!(a.oil_sm3, b.oil_sm3, "oil vector not bit-reproducible");
        assert_eq!(a.gas_sm3, b.gas_sm3);
        assert_eq!(a.grv_m3, b.grv_m3);
        // A different seed gives a different stream.
        let c = run_mc(&mut t1, &inputs, &McSettings::new(200, 43)).unwrap();
        assert_ne!(a.oil_sm3, c.oil_sm3);
    }

    #[test]
    fn sharded_mc_matches_serial_and_is_worker_invariant() {
        // The sharding contract (`run_structured_mc_parallel`): for a run with no
        // per-draw structural shifts, the parallel driver reproduces the serial
        // vectors EXACTLY at every worker count, and is bit-reproducible run to run.
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let inputs = nominal_inputs();
        let (n, seed) = (250usize, 42u64);

        let mut ts = StaticModelTemplate::new(&wf, opts()).unwrap();
        let serial = run_mc(&mut ts, &inputs, &McSettings::new(n, seed)).unwrap();

        for workers in [1usize, 2, 3, 5, 8] {
            let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
            let par = run_mc(
                &mut t,
                &inputs,
                &McSettings::new(n, seed).with_workers(workers),
            )
            .unwrap();
            assert_eq!(par.len(), n, "workers={workers} length");
            assert_eq!(
                par.oil_sm3, serial.oil_sm3,
                "workers={workers}: oil vector diverged from serial"
            );
            assert_eq!(par.gas_sm3, serial.gas_sm3, "workers={workers}: gas");
            assert_eq!(par.grv_m3, serial.grv_m3, "workers={workers}: grv");
        }

        // Same (seed, n, workers) twice -> bit-identical (no scheduling nondeterminism).
        let mut ta = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut tb = StaticModelTemplate::new(&wf, opts()).unwrap();
        let a = run_mc(&mut ta, &inputs, &McSettings::new(n, seed).with_workers(4)).unwrap();
        let b = run_mc(&mut tb, &inputs, &McSettings::new(n, seed).with_workers(4)).unwrap();
        assert_eq!(
            a.oil_sm3, b.oil_sm3,
            "same (seed,n,workers) not reproducible"
        );

        // Workers is clamped to [1, n]: an over-large count still runs correctly.
        let mut tc = StaticModelTemplate::new(&wf, opts()).unwrap();
        let over = run_mc(
            &mut tc,
            &inputs,
            &McSettings::new(n, seed).with_workers(10_000),
        )
        .unwrap();
        assert_eq!(over.oil_sm3, serial.oil_sm3, "over-large worker count");
    }

    #[test]
    fn result_vectors_have_length_n() {
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let r = run_mc(&mut t, &nominal_inputs(), &McSettings::new(137, 1)).unwrap();
        assert_eq!(r.len(), 137);
        assert_eq!(r.oil_sm3.len(), 137);
        assert_eq!(r.gas_sm3.len(), 137);
        assert_eq!(r.grv_m3.len(), 137);
        assert!(r.oil_sm3.iter().all(|v| v.is_finite() && *v > 0.0));
        let s = r.summary().unwrap();
        assert!(s.p90 <= s.p50 && s.p50 <= s.p10, "P-curve ordered: {s:?}");
    }

    #[test]
    fn clamped_garbage_sampler_never_aborts() {
        // Wild Normal priors that WITHOUT a clamp would routinely leave [0,1] and
        // trip realize's H2 guard. Clamped into a valid sub-range, every one of
        // 1000 draws stays valid — the run completes, no McDraw error.
        let wild = |mean: f64, sd: f64, lo: f64, hi: f64| {
            Input::clamped(Sampler::new_normal(mean, sd).unwrap(), lo, hi).unwrap()
        };
        let inputs = McInputs::new(
            wild(100.0, 50.0, 10.0, 200.0), // area (positive)
            wild(50.0, 40.0, 5.0, 120.0),   // gross (positive)
            tri(5020.0, 5025.0, 5030.0),    // contact
            wild(0.20, 0.30, 0.02, 0.35),   // porosity — sd 0.30 wildly overspills
            wild(0.60, 0.40, 0.05, 0.95),   // ntg
            wild(0.30, 0.30, 0.02, 0.60),   // sw
            wild(1.30, 0.50, 1.05, 1.80),   // boi (>= 1)
        );
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let r = run_mc(&mut t, &inputs, &McSettings::new(1000, 7)).unwrap();
        assert_eq!(r.len(), 1000);
        assert!(r.oil_sm3.iter().all(|v| v.is_finite() && *v > 0.0));
    }

    #[test]
    fn two_contact_run_produces_gas_and_oil_vectors() {
        // A gas cap above the drawn GOC + a gas FVF -> both legs populated.
        let wf = flat_wireframe(11, 5000.0, 5040.0);
        let inputs = McInputs::new(
            tri(90.0, 100.0, 110.0),
            tri(45.0, 50.0, 55.0),
            tri(5038.0, 5040.0, 5042.0), // OWC
            tri(0.18, 0.22, 0.26),
            tri(0.70, 0.80, 0.90),
            tri(0.25, 0.30, 0.35),
            tri(1.20, 1.30, 1.45),
        )
        .with_goc(tri(5018.0, 5020.0, 5022.0))
        .with_bgi(tri(0.0035, 0.004, 0.0045));
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let r = run_mc(&mut t, &inputs, &McSettings::new(100, 11)).unwrap();
        assert!(r.oil_sm3.iter().all(|v| *v > 0.0), "oil leg positive");
        assert!(r.gas_sm3.iter().all(|v| *v > 0.0), "gas cap positive");
    }

    #[test]
    fn level_shift_and_resimulate_both_drive_through_the_loop() {
        use petektools::{Variogram, VariogramModel};
        // Target the PORO cube itself so the geostatistical pipeline *drives
        // volumetrics* — otherwise a side-property (e.g. PHIE) leaves the oil
        // metric untouched (in-place reads PORO/NTG/SW). With PORO the two MC
        // composition modes produce genuinely different oil realization sets.
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let poro = |seed: u64| {
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
            PropertyPipeline::new("PORO")
                .upscale(vec![low, high], UpscaleMethod::Arithmetic)
                .propagate(Gaussian::new(vgm, seed))
        };

        // LevelShift template + a per-draw shift input on PORO.
        let mut tl = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property(poro(42));
        let ls_inputs = nominal_inputs().with_property_shift("PORO", tri(-0.02, 0.0, 0.02));
        // The shift is plumbed onto each draw the loop builds.
        let realized = ls_inputs.realize(100, 3).unwrap();
        assert!(
            (realized.draw_at(0).property_shift("PORO") - realized.property_shifts[0].1[0]).abs()
                < 1e-15,
            "the sampled PORO shift must reach the draw"
        );
        let ls = run_mc(&mut tl, &ls_inputs, &McSettings::new(100, 3)).unwrap();
        assert_eq!(ls.len(), 100);
        assert!(ls.oil_sm3.iter().all(|v| v.is_finite() && *v > 0.0));

        // Resimulate template — no property-shift input needed (seed_index drives it).
        let mut tr = StaticModelTemplate::new(&wf, opts())
            .unwrap()
            .with_property_mode(poro(42), McMode::Resimulate);
        let rs = run_mc(&mut tr, &nominal_inputs(), &McSettings::new(100, 3)).unwrap();
        assert_eq!(rs.len(), 100);
        assert!(rs.oil_sm3.iter().all(|v| v.is_finite() && *v > 0.0));
        // The two composition modes give genuinely different realization sets.
        assert_ne!(ls.oil_sm3, rs.oil_sm3);
    }

    #[test]
    fn aggregation_brackets_are_ordered_sensibly() {
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t1 = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut t2 = StaticModelTemplate::new(&wf, opts()).unwrap();
        let seg1 = run_mc(&mut t1, &nominal_inputs(), &McSettings::new(500, 101)).unwrap();
        let seg2 = run_mc(&mut t2, &nominal_inputs(), &McSettings::new(500, 202)).unwrap();

        let ind = aggregate_field(&[&seg1, &seg2], Correlation::Independent);
        let com = aggregate_field(&[&seg1, &seg2], Correlation::Comonotonic);
        let s_ind = reservoir_summary(&ind).unwrap();
        let s_com = reservoir_summary(&com).unwrap();
        // Same mean (a hard identity of the two summing modes).
        assert!((s_ind.mean - s_com.mean).abs() < 1e-6, "means must match");
        // Comonotonic is the wider (downside) bracket: lower P90, higher P10.
        assert!(
            s_com.p90 <= s_ind.p90 + 1e-9,
            "comonotonic P90 not <= independent"
        );
        assert!(
            s_com.p10 >= s_ind.p10 - 1e-9,
            "comonotonic P10 not >= independent"
        );
    }

    #[test]
    fn n_zero_is_a_typed_error() {
        let wf = flat_wireframe(9, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let err = run_mc(&mut t, &nominal_inputs(), &McSettings::new(0, 1)).unwrap_err();
        assert!(matches!(err, StaticError::InvalidInput(_)));
    }

    #[test]
    fn bad_draw_surfaces_as_typed_mcdraw_with_index() {
        // A contact sampler pinned far ABOVE the top so every column is fully
        // below contact -> zero HCPV is fine; instead force a real H2 failure by
        // an unclamped porosity that overspills [0,1]. The first offending draw
        // must surface as McDraw carrying an index (fail-fast).
        let inputs = McInputs::new(
            tri(90.0, 100.0, 110.0),
            tri(45.0, 50.0, 55.0),
            tri(5020.0, 5025.0, 5030.0),
            Input::plain(Sampler::new_uniform(0.9, 1.5).unwrap()), // porosity > 1 sometimes
            tri(0.70, 0.80, 0.90),
            tri(0.25, 0.30, 0.35),
            tri(1.20, 1.30, 1.45),
        );
        let wf = flat_wireframe(9, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let err = run_mc(&mut t, &inputs, &McSettings::new(500, 5)).unwrap_err();
        match err {
            StaticError::McDraw { index, source } => {
                assert!(index < 500, "index in range");
                assert!(matches!(*source, StaticError::InvalidInput(_)));
            }
            other => panic!("expected McDraw, got {other}"),
        }
    }

    #[test]
    fn smoke_10k_draws_with_timing() {
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let mut t = StaticModelTemplate::new(&wf, opts()).unwrap();
        let start = std::time::Instant::now();
        let r = run_mc(&mut t, &nominal_inputs(), &McSettings::new(10_000, 2024)).unwrap();
        let dt = start.elapsed();
        assert_eq!(r.len(), 10_000);
        let s = r.summary().unwrap();
        assert!(s.p90 <= s.p50 && s.p50 <= s.p10);
        eprintln!(
            "[mc smoke] 10k draws (11x11x5) in {:.3?} ({:.1} us/draw); oil P90/P50/P10 = {:.3}/{:.3}/{:.3} Sm3",
            dt,
            dt.as_secs_f64() * 1e6 / 10_000.0,
            s.p90,
            s.p50,
            s.p10
        );
    }
    #[test]
    #[allow(deprecated)]
    fn deprecated_wrappers_match_run_mc() {
        // The four historical entries are thin wrappers over `run_mc` — pinned
        // bit-identical so the deprecation window cannot drift.
        let wf = flat_wireframe(11, 5000.0, 5025.0);
        let inputs = nominal_inputs();
        let (n, seed) = (60usize, 9u64);
        let dir = std::env::temp_dir();

        let mut a = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut b = StaticModelTemplate::new(&wf, opts()).unwrap();
        let new = run_mc(&mut a, &inputs, &McSettings::new(n, seed)).unwrap();
        let old = run_structured_mc(&mut b, &inputs, n, seed).unwrap();
        assert_eq!(new.oil_sm3, old.oil_sm3, "serial wrapper");

        let mut a = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut b = StaticModelTemplate::new(&wf, opts()).unwrap();
        let new = run_mc(&mut a, &inputs, &McSettings::new(n, seed).with_workers(3)).unwrap();
        let old = run_structured_mc_parallel(&mut b, &inputs, n, seed, 3).unwrap();
        assert_eq!(new.oil_sm3, old.oil_sm3, "parallel wrapper");

        let mut a = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut b = StaticModelTemplate::new(&wf, opts()).unwrap();
        let new = run_mc(
            &mut a,
            &inputs,
            &McSettings::new(n, seed).with_spill_dir(dir.clone()),
        )
        .unwrap();
        let old = run_structured_mc_spilled(&mut b, &inputs, n, seed, None).unwrap();
        assert_eq!(new.oil_sm3, old.oil_sm3, "spilled wrapper");

        let mut a = StaticModelTemplate::new(&wf, opts()).unwrap();
        let mut b = StaticModelTemplate::new(&wf, opts()).unwrap();
        let new = run_mc(
            &mut a,
            &inputs,
            &McSettings::new(n, seed).with_workers(2).with_spill_dir(dir),
        )
        .unwrap();
        let old = run_structured_mc_parallel_spilled(&mut b, &inputs, n, seed, 2, None).unwrap();
        assert_eq!(new.oil_sm3, old.oil_sm3, "parallel-spilled wrapper");
    }
}
