//! The [`StaticModel`] aggregate — a *populated* static reservoir model: the
//! single, complete, self-describing artifact this library produces (SPEC §2).
//! It is a value (owned, `Clone`), not a handle to a live builder; construction
//! (via the builder/template, crate-internal) guarantees the invariants, so
//! consumers never re-validate.
//!
//! Per the layer charter (graph `decision_layer_charters`) the model **owns its
//! volumetrics output surface**: GRV / in-place come off the model itself
//! ([`StaticModel::in_place`]), with FVF applied by the caller through
//! srs-volumetrics' validated value types ([`crate::volumetrics::OilFvf`] /
//! [`crate::volumetrics::GasFvf`]). P-curves aggregate per-realization results via
//! srs-uncertainty (`PercentileSummary`).

use crate::error::StaticError;
use crate::grid::{build_box, BoxSpec, Dims, Grid, Property};
use crate::gridder::{Conformity, SolveOpts};
use crate::model::provenance::{PopulationMode, Provenance};
use crate::model::zones::ZoneTable;
use crate::spill::SpillBacking;
use crate::volumetrics::{compute_clipped, Clip, GridSource, InPlace, ZoneVolumes};
use crate::wireframe::{Contact, ContactKind, Wireframe};
use std::path::Path;
use std::sync::Arc;

/// The model's registered **world** georeference: the world `(x, y)` of column
/// `(0, 0)`'s centroid and the column spacing. It is the mapping from the grid's
/// local column lattice `(i, j)` to world coordinates — the same `xy↔ij` the
/// upstream well registration uses — so the view bundles can emit their shared
/// areal frame in ONE consistent world frame (outline + wells + raster overlay,
/// world fence/bore sections). `None` on the model is the **degenerate local
/// case** (a synthetic square / box with no world georeference), for which the
/// frame falls back to the grid's own local column-centroid lattice.
///
/// Column `(i, j)`'s world centroid is
/// `(origin_x + i * spacing_x, origin_y + j * spacing_y)`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Georef {
    /// World `x` of column `(0, 0)`'s centroid.
    pub origin_x: f64,
    /// World `y` of column `(0, 0)`'s centroid.
    pub origin_y: f64,
    /// World column spacing along `x` (metres).
    pub spacing_x: f64,
    /// World column spacing along `y` (metres).
    pub spacing_y: f64,
}

impl Georef {
    /// Build a world georeference, returning `None` (→ the local degenerate case)
    /// when the spacing is not finite-positive — so a caller may pass through an
    /// unvalidated georeference and get the safe local fallback, mirroring
    /// [`crate::model::TrendSurface::with_georef`].
    #[must_use]
    pub fn new(origin_x: f64, origin_y: f64, spacing_x: f64, spacing_y: f64) -> Option<Self> {
        if origin_x.is_finite()
            && origin_y.is_finite()
            && spacing_x.is_finite()
            && spacing_x > 0.0
            && spacing_y.is_finite()
            && spacing_y > 0.0
        {
            Some(Self {
                origin_x,
                origin_y,
                spacing_x,
                spacing_y,
            })
        } else {
            None
        }
    }
}

/// A populated static reservoir model: structural framework + corner-point grid
/// + property cubes + zones + contacts + provenance.
#[derive(Debug, Clone)]
pub struct StaticModel {
    framework: Wireframe,
    grid: Grid,
    zones: ZoneTable,
    provenance: Provenance,
    /// Optional gas-cap connate-water override (R3): applied to gas-zone cells in
    /// the two-contact `in_place` split so a single shared `SW` cube does not
    /// over-state gas-cap OGIP. `None` = the gas cap uses the `SW` cube.
    sw_gas: Option<f64>,
    /// The registered world georeference (see [`Georef`]); `None` = the grid is
    /// its own local frame (a synthetic square / box) and the view frame
    /// degenerates to the local column-centroid lattice.
    georef: Option<Georef>,
    /// Per-zone fluid contacts, parallel to `zones().zones()` — the multi-zone
    /// (`from_horizon_stack`) framework where each zone has its **own** accumulation
    /// (or none). Drives [`StaticModel::in_place_by_zone`]. `None` for a legacy
    /// single-contact model (the whole-model contacts on the framework govern).
    zone_contacts: Option<Vec<Vec<Contact>>>,
    /// The out-of-core **backing-storage mode** (rulings R1/R2/R4; `srs_spill`).
    /// `None` (the default) is the in-core model, where geometry and cubes are owned
    /// by `grid` and every accessor is byte-identical to the pre-out-of-core
    /// behaviour. `Some(_)` is a **spilled** model whose heavy per-cell arrays (ZCORN
    /// and cubes, f32) live in a memory-mapped petekTools store with `grid` reduced
    /// to a unit placeholder; the volumetric surface (`in_place`, `bulk_volume`)
    /// streams through the store's windowed views while every other field (framework,
    /// zones, contacts, provenance, georef) stays resident. Shared behind an `Arc` so
    /// the model stays `Clone` (the mmap is shared, not copied) and the store file is
    /// cleaned up when the last clone drops (unless the backing was detached).
    /// `grid()` and `property()` on a spilled model return the placeholder, so read a
    /// spilled model through its volumetric methods; raw cube borrows are the in-core
    /// representation only.
    spill: Option<Arc<SpillBacking>>,
}

impl StaticModel {
    /// Crate-internal constructor — only the builder/template build models, so
    /// the SPEC §2 invariants (cube lengths == cell count, zone k-ranges within
    /// nk, lattice consistency) hold by construction.
    pub(crate) fn new(
        framework: Wireframe,
        grid: Grid,
        zones: ZoneTable,
        provenance: Provenance,
        sw_gas: Option<f64>,
    ) -> Self {
        Self {
            framework,
            grid,
            zones,
            provenance,
            sw_gas,
            georef: None,
            zone_contacts: None,
            spill: None,
        }
    }

    /// Convert a freshly-built in-core model into a **spilled** one (crate-internal;
    /// the builder/template call it above the memory budget, ruling R5). The heavy
    /// in-core `grid` is dropped — replaced by a unit placeholder — and the real
    /// geometry + cubes now live in the memory-mapped `backing` store. Every other
    /// field (framework, zones, contacts, provenance, sw_gas, georef, zone_contacts)
    /// is preserved, so the model's non-volumetric surface is unchanged; the
    /// volumetric surface routes through the store.
    pub(crate) fn into_spilled(mut self, backing: SpillBacking) -> Self {
        self.grid = build_box(BoxSpec::square(
            1.0,
            1.0,
            Dims::new(1, 1, 1).expect("1x1x1 dims are valid"),
        ))
        .expect("unit box always builds");
        self.spill = Some(Arc::new(backing));
        self
    }

    /// Construct a **spilled** model directly from a store backing, **without ever
    /// building a whole in-core grid** — the slab-incremental build path
    /// ([`crate::model::StaticModelBuilder`] above budget, R2). The geometry + cubes already
    /// live in `backing`'s mmap store (written slab-by-slab); the `grid` is a unit
    /// placeholder. Equivalent to `new(..).with_georef_opt(..).into_spilled(..)` but
    /// never allocates the `O(grid)` transient the post-hoc spill path does.
    pub(crate) fn spilled(
        framework: Wireframe,
        zones: ZoneTable,
        provenance: Provenance,
        sw_gas: Option<f64>,
        georef: Option<Georef>,
        backing: SpillBacking,
    ) -> Self {
        let placeholder = build_box(BoxSpec::square(
            1.0,
            1.0,
            Dims::new(1, 1, 1).expect("1x1x1 dims are valid"),
        ))
        .expect("unit box always builds");
        Self::new(framework, placeholder, zones, provenance, sw_gas)
            .with_georef_opt(georef)
            .into_spilled(backing)
    }

    /// Attach (or clear) a spill backing **without** disturbing the in-core grid —
    /// crate-internal, for the structured-MC spilled loop (ruling R3). The reusable
    /// model's in-core grid buffers stay intact (the next `realize_into` recycles
    /// them); attaching a backing only re-routes the volumetric surface
    /// ([`StaticModel::in_place_impl`]) to stream the just-spilled f32 store, so the
    /// summary reads through mmap windows. Cleared again after the summary so the
    /// next draw realizes in-core. Distinct from [`StaticModel::into_spilled`],
    /// which *replaces* the grid with a placeholder (the persistent build path).
    pub(crate) fn set_spill(&mut self, spill: Option<Arc<SpillBacking>>) {
        self.spill = spill;
    }

    /// Whether this model is spilled (out-of-core backing-storage mode).
    #[must_use]
    pub fn is_spilled(&self) -> bool {
        self.spill.is_some()
    }

    /// The spill store path, if this model is spilled (ruling R5: the store
    /// location is observable).
    #[must_use]
    pub fn spill_store_path(&self) -> Option<&Path> {
        self.spill.as_ref().map(|b| b.store_path())
    }

    /// A whole in-core [`Grid`] for the **non-hot-path view exports** — the model's
    /// own grid in-core, or one reconstructed from the mmap backing when spilled (a
    /// spilled model's `grid()` is a 1×1×1 placeholder, so the shell/section exports
    /// must materialize the real geometry + cubes to read them). The `Owned` branch
    /// allocates an `O(grid)` transient, so this is for one-shot exports only, never
    /// the realization loop.
    pub(crate) fn view_grid(&self) -> Result<std::borrow::Cow<'_, Grid>, StaticError> {
        match &self.spill {
            Some(backing) => Ok(std::borrow::Cow::Owned(backing.to_in_core_grid()?)),
            None => Ok(std::borrow::Cow::Borrowed(&self.grid)),
        }
    }

    /// Grid dimensions — the real dims whether in-core or spilled (a spilled
    /// model's `grid()` is a placeholder, so read dims through this).
    #[must_use]
    pub fn dims(&self) -> Dims {
        match &self.spill {
            Some(b) => b.dims(),
            None => self.grid.dims(),
        }
    }

    /// Attach the registered world [`Georef`] (crate-internal; the builder /
    /// template thread it through from their own `with_georef`). `None` leaves the
    /// model in the local degenerate frame.
    pub(crate) fn with_georef_opt(mut self, georef: Option<Georef>) -> Self {
        self.georef = georef;
        self
    }

    /// Attach the per-zone contact sets (crate-internal; the stack build threads
    /// them through parallel to `zones()`). `None` leaves the model on the legacy
    /// whole-model contact path.
    pub(crate) fn with_zone_contacts(mut self, zone_contacts: Option<Vec<Vec<Contact>>>) -> Self {
        self.zone_contacts = zone_contacts;
        self
    }

    /// A structurally-empty single-cell model to be filled by
    /// [`crate::model::StaticModelTemplate::realize_into`] — every field is overwritten on
    /// the first realize (the grid geometry + cube buffers then grow in place and are
    /// recycled thereafter). The allocating [`crate::model::StaticModelTemplate::realize`]
    /// starts here, so both share one build path. `framework` is a throwaway seed
    /// (the realize overwrites it); the placeholder provenance never escapes.
    pub(crate) fn empty(framework: Wireframe) -> Self {
        let grid = build_box(BoxSpec::square(
            1.0,
            1.0,
            Dims::new(1, 1, 1).expect("1x1x1 dims are valid"),
        ))
        .expect("unit box always builds");
        Self {
            framework,
            grid,
            zones: ZoneTable::single(1),
            provenance: Provenance {
                inputs_ref: String::new(),
                solve_opts: SolveOpts::default(),
                conformity: Conformity::Proportional,
                nk: 1,
                population: PopulationMode::Priors,
                realization: None,
                warnings: Vec::new(),
                property_reports: Vec::new(),
                stack: None,
                well_ties: Vec::new(),
                sugar_cube: false,
            },
            sw_gas: None,
            georef: None,
            zone_contacts: None,
            spill: None,
        }
    }

    /// Mutable access to the corner-point grid — crate-internal, for the
    /// buffer-recycling realize path ([`crate::model::StaticModelTemplate::realize_into`]).
    pub(crate) fn grid_mut(&mut self) -> &mut Grid {
        &mut self.grid
    }

    /// Overwrite the reused model's non-geometry state after its grid was recycled in
    /// place (crate-internal; [`crate::model::StaticModelTemplate::realize_into`]). Sets
    /// exactly the fields `StaticModel::new(..).with_georef_opt(..)` would for a
    /// single-contact realization, so a recycled model is indistinguishable from a
    /// freshly built one.
    pub(crate) fn reset_state(
        &mut self,
        framework: Wireframe,
        zones: ZoneTable,
        provenance: Provenance,
        sw_gas: Option<f64>,
        georef: Option<Georef>,
        zone_contacts: Option<Vec<Vec<Contact>>>,
    ) {
        self.framework = framework;
        self.zones = zones;
        self.provenance = provenance;
        self.sw_gas = sw_gas;
        self.georef = georef;
        self.zone_contacts = zone_contacts;
        // A recycled realize_into target is always in-core (its grid buffers were
        // just refilled in place); clear any spill backing from a prior life.
        self.spill = None;
    }

    /// The registered world georeference, if any (see [`Georef`]). `None` = the
    /// grid is its own local frame (a synthetic square / box) — the view frame
    /// then degenerates to the grid's local column-centroid lattice.
    #[must_use]
    pub fn georef(&self) -> Option<Georef> {
        self.georef
    }

    // --- read-only accessors (the consumer surface, SPEC §7b) ---

    /// The corner-point grid (geometry + cubes + `cells()`).
    #[must_use]
    pub fn grid(&self) -> &Grid {
        &self.grid
    }

    /// The structural framework (boundary + horizons + contacts).
    #[must_use]
    pub fn framework(&self) -> &Wireframe {
        &self.framework
    }

    /// The fluid contacts (convenience for `framework().contacts`).
    #[must_use]
    pub fn contacts(&self) -> &[Contact] {
        &self.framework.contacts
    }

    /// A property cube by name (`PORO` / `SW` / `NTG` / …).
    #[must_use]
    pub fn property(&self, name: &str) -> Option<&Property> {
        self.grid.properties().get(name)
    }

    /// The names of the populated cubes (unordered).
    #[must_use]
    pub fn property_names(&self) -> Vec<&str> {
        self.grid.properties().names().collect()
    }

    /// The stratigraphic zone table.
    #[must_use]
    pub fn zones(&self) -> &ZoneTable {
        &self.zones
    }

    /// Per-zone summary statistics of a property cube (count / mean / min / max),
    /// over the **active** cells (non-truncated, `dz > 0`) of each zone's `k`-range
    /// — the geology-addressed view of a full-grid cube (the cubes stay full-grid;
    /// this indexes them by zone). Returns one [`ZoneStat`] per zone, top→base; a
    /// zone with no active cells reports `count == 0` and `NaN` aggregates.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the named cube is not populated.
    pub fn zone_stats(&self, property: &str) -> Result<Vec<ZoneStat>, StaticError> {
        if self.spill.is_some() {
            // Per-zone cube statistics window the cubes randomly across all
            // k-slabs; the spilled read surface is v1-scoped to the volumetric
            // re-cuts (in_place*). Route zone stats through an in-core build.
            return Err(StaticError::InvalidInput(
                "zone_stats is not available on a spilled (out-of-core) model in v1".into(),
            ));
        }
        let cube = self.grid.properties().get(property).ok_or_else(|| {
            StaticError::InvalidInput(format!("grid is missing property '{property}'"))
        })?;
        let dims = self.grid.dims();
        let mut out = Vec::with_capacity(self.zones.zones().len());
        for zone in self.zones.zones() {
            let (mut n, mut sum, mut min, mut max) =
                (0usize, 0.0f64, f64::INFINITY, f64::NEG_INFINITY);
            for k in zone.k_range.clone() {
                for j in 0..dims.nj {
                    for i in 0..dims.ni {
                        let c = crate::grid::Ijk::new(i, j, k);
                        if self.grid.cell(c).dz() <= 1e-9 {
                            continue; // inactive/truncated/collapsed
                        }
                        let v = cube.values[(k * dims.nj + j) * dims.ni + i];
                        if !v.is_finite() {
                            continue;
                        }
                        n += 1;
                        sum += v;
                        min = min.min(v);
                        max = max.max(v);
                    }
                }
            }
            out.push(if n == 0 {
                ZoneStat {
                    zone: zone.name.clone(),
                    count: 0,
                    mean: f64::NAN,
                    min: f64::NAN,
                    max: f64::NAN,
                }
            } else {
                ZoneStat {
                    zone: zone.name.clone(),
                    count: n,
                    mean: sum / n as f64,
                    min,
                    max,
                }
            });
        }
        Ok(out)
    }

    /// What produced this model (inputs, gridder settings, population mode, and
    /// the realization draw if this is an MC realization).
    #[must_use]
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }

    /// Gross rock volume of the whole grid \[m³\] (pure geometry). In-core: sums
    /// the grid. Spilled: streams the ZCORN slabs through the store (`NaN` on the
    /// rare store-read failure — never a panic).
    #[must_use]
    pub fn bulk_volume(&self) -> f64 {
        match &self.spill {
            Some(b) => b.bulk_volume().unwrap_or(f64::NAN),
            None => self.grid.bulk_volume(),
        }
    }

    // --- the volumetrics output surface (`decision_layer_charters`) ---

    /// GRV + HCPV of the hydrocarbon column: the model's own volumes, clipped
    /// against its fluid contact(s). Apply FVF on the result
    /// ([`InPlace::ooip_sm3`] / [`InPlace::ogip_sm3`]) — FVF is an uncertain
    /// scalar *input* at this seam, never PVT code.
    ///
    /// **Two-contact columns.** When the framework carries a GOC *and* a lower
    /// contact (OWC/FWL or GWC), this returns the gas-cap + oil-rim split
    /// (`InPlace::gas` / `InPlace::oil`); read the per-zone in-place with
    /// [`InPlace::gas_zone_ogip_sm3`] / [`InPlace::oil_zone_ooip_sm3`].
    /// Otherwise it clips the single (first) contact as a generic hydrocarbon
    /// column. A gas-cap connate-water override (`with_sw_gas`, R3) applies to the
    /// gas-zone cells of the two-contact split.
    ///
    /// # Errors
    /// [`StaticError::Grid`] if the model carries no fluid contact;
    /// [`StaticError::InvalidInput`] if a required cube is missing or
    /// non-physical (H2), or the drawn GOC sits below the OWC.
    pub fn in_place(&self) -> Result<InPlace, StaticError> {
        self.in_place_impl(true)
    }

    /// Summary-only in-place (V7): the same GRV/HCPV/OGIP/OOIP aggregates as
    /// [`StaticModel::in_place`] but **without materializing the per-cell HCPV
    /// cube** ([`InPlace::per_cell_hcpv`] is left empty). Use this on the MC hot
    /// path, where only the aggregates feed the P-curve and the per-cell map is
    /// never read — it avoids a `cell_count`-length allocation per realization.
    ///
    /// # Errors
    /// Same as [`StaticModel::in_place`].
    pub fn in_place_summary(&self) -> Result<InPlace, StaticError> {
        self.in_place_impl(false)
    }

    /// Resolve the whole-model contact geometry and dispatch to the unified
    /// volumetrics core (`per_cell` gates the per-cell HCPV cube). A model with no
    /// fluid contact is an error (unlike a contactless *zone*, which is bulk).
    fn in_place_impl(&self, per_cell: bool) -> Result<InPlace, StaticError> {
        let clip = clip_of(self.contacts(), self.sw_gas)
            .ok_or_else(|| StaticError::Grid("static model has no fluid contact".into()))?;
        self.volumetrics(clip, 0..self.dims().nk, per_cell)
    }

    /// Dispatch a [`Clip`] over `k_range` to the ONE unified GRV/HCPV core
    /// (`crate::volumetrics::compute_clipped`), picking the backing: in-core reads the
    /// contiguous f64 [`Grid`]; spilled STREAMS the mmap store one k-slab at a time
    /// (R2/R4 — same aggregates, f64 accumulation, tolerance-parity). The single
    /// place the in-core-vs-spilled fork lives.
    fn volumetrics(
        &self,
        clip: Clip,
        k_range: core::ops::Range<usize>,
        per_cell: bool,
    ) -> Result<InPlace, StaticError> {
        match &self.spill {
            Some(backing) => compute_clipped(&backing.source(), clip, k_range, per_cell),
            None => compute_clipped(&GridSource::new(&self.grid), clip, k_range, per_cell),
        }
    }

    /// Per-zone in-place with a total rollup (SPEC §4/§8) — the authoritative
    /// multi-zone answer for a stack-built model where each zone carries its **own**
    /// fluid contacts. Every zone's hydrocarbons are clipped against *its* contacts
    /// over *its* `k`-range: a two-contact zone splits gas-cap / oil-leg, a
    /// single-contact zone clips one leg, and a **contactless** zone contributes its
    /// gross bulk volume but **zero** hydrocarbon in-place (no contact = no known
    /// accumulation). [`ZonedInPlace::total`] is the sum over zones (GRV / HCPV /
    /// gas / oil), so `sum(zone volumes) == total` to FP tolerance (conservation).
    ///
    /// A legacy single-contact model (no per-zone contacts) resolves every zone
    /// against the whole-model framework contacts — for the single implicit zone
    /// this reproduces [`StaticModel::in_place_summary`].
    ///
    /// # Errors
    /// [`StaticError`] if a zone's contact geometry is invalid (e.g. a drawn GOC
    /// below its OWC) or a required cube is missing/non-physical.
    pub fn in_place_by_zone(&self) -> Result<ZonedInPlace, StaticError> {
        // Zones are contiguous k-bands (`Zone.k_range`), so per-zone volumetrics is a
        // k-range-restricted stream through the same unified core — spilled works
        // identically to in-core (each zone's slabs are windowed from the store).
        let mut zones = Vec::with_capacity(self.zones.zones().len());
        let (mut t_grv, mut t_hcpv, mut t_cells) = (0.0, 0.0, 0usize);
        let (mut t_gas, mut t_oil): (Option<ZoneVolumes>, Option<ZoneVolumes>) = (None, None);
        for (z, zone) in self.zones.zones().iter().enumerate() {
            let contacts: &[Contact] = match &self.zone_contacts {
                Some(zc) => &zc[z],
                None => self.contacts(),
            };
            let ip = self.zone_in_place(contacts, zone.k_range.clone())?;
            t_grv += ip.grv_m3;
            t_hcpv += ip.hcpv_m3;
            t_cells += ip.cells_in_column;
            t_gas = add_zone_volumes(t_gas, ip.gas);
            t_oil = add_zone_volumes(t_oil, ip.oil);
            zones.push(ZoneInPlace {
                zone: zone.name.clone(),
                in_place: ip,
            });
        }
        Ok(ZonedInPlace {
            zones,
            total: InPlace {
                grv_m3: t_grv,
                hcpv_m3: t_hcpv,
                cells_in_column: t_cells,
                per_cell_hcpv: Vec::new(),
                gas: t_gas,
                oil: t_oil,
            },
        })
    }

    /// Resolve one zone's contact geometry and compute its in-place over `k_range`
    /// through the unified dispatch: two-contact (GOC + a lower OWC/GWC) split,
    /// single-contact clip, or — with no contact — gross bulk with zero hydrocarbon
    /// ([`Clip::Bulk`]). Backing-agnostic: works in-core and spilled.
    fn zone_in_place(
        &self,
        contacts: &[Contact],
        k_range: core::ops::Range<usize>,
    ) -> Result<InPlace, StaticError> {
        let clip = clip_of(contacts, self.sw_gas).unwrap_or(Clip::Bulk);
        self.volumetrics(clip, k_range, false)
    }
}

/// Resolve a contact list to its volumetric [`Clip`] — the ONE home of the
/// two-contact (GOC + deepest OWC/GWC, `goc <= owc`) vs single-contact resolution
/// that `in_place_impl` and `zone_in_place` previously each re-derived. `None` =
/// no fluid contact (the caller decides: an error for the whole model, `Bulk` for
/// a zone).
fn clip_of(contacts: &[Contact], sw_gas: Option<f64>) -> Option<Clip> {
    let goc = contacts
        .iter()
        .find(|c| c.kind == ContactKind::Goc)
        .map(|c| c.depth_m);
    // The lower (base) contact: the deepest OWC/FWL or GWC present.
    let lower = contacts
        .iter()
        .filter(|c| matches!(c.kind, ContactKind::Owc | ContactKind::Gwc))
        .map(|c| c.depth_m)
        .fold(None, |acc: Option<f64>, d| {
            Some(acc.map_or(d, |a| a.max(d)))
        });
    if let (Some(g), Some(w)) = (goc, lower) {
        if g <= w {
            return Some(Clip::Two {
                goc: g,
                owc: w,
                sw_gas,
            });
        }
    }
    contacts.first().map(|c| Clip::Single(c.depth_m))
}

/// Per-zone summary statistics of a property cube (from [`StaticModel::zone_stats`]),
/// over the active cells of one zone. `count == 0` (with `NaN` aggregates) when the
/// zone has no active cells.
#[derive(Debug, Clone, PartialEq)]
pub struct ZoneStat {
    /// The zone's name.
    pub zone: String,
    /// Active-cell count in the zone.
    pub count: usize,
    /// Arithmetic mean of the cube over the zone's active cells.
    pub mean: f64,
    /// Minimum value over the zone's active cells.
    pub min: f64,
    /// Maximum value over the zone's active cells.
    pub max: f64,
}

/// One zone's in-place result (from [`StaticModel::in_place_by_zone`]).
#[derive(Debug, Clone)]
pub struct ZoneInPlace {
    /// The zone's name.
    pub zone: String,
    /// The zone's GRV / HCPV (+ gas/oil split if it is a two-contact zone; zero
    /// hydrocarbon if it is contactless).
    pub in_place: InPlace,
}

/// Per-zone in-place plus the total rollup (from [`StaticModel::in_place_by_zone`]).
/// `total` is the sum over `zones` (GRV / HCPV / gas / oil), so
/// `sum(zone volumes) == total` to FP tolerance.
#[derive(Debug, Clone)]
pub struct ZonedInPlace {
    /// Per-zone results, ordered top→base.
    pub zones: Vec<ZoneInPlace>,
    /// The rollup total (summary-only: `per_cell_hcpv` is empty).
    pub total: InPlace,
}

/// Accumulate an optional per-zone gas/oil [`ZoneVolumes`] into a running total.
fn add_zone_volumes(acc: Option<ZoneVolumes>, add: Option<ZoneVolumes>) -> Option<ZoneVolumes> {
    match (acc, add) {
        (a, None) => a,
        (None, Some(v)) => Some(v),
        (Some(a), Some(v)) => Some(ZoneVolumes {
            grv_m3: a.grv_m3 + v.grv_m3,
            hcpv_m3: a.hcpv_m3 + v.hcpv_m3,
            cells: a.cells + v.cells,
        }),
    }
}
