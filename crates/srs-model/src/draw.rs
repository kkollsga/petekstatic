//! [`RealizationDraw`] — the per-realization sampled input set the MC loop feeds
//! [`crate::StaticModelTemplate::realize`].
//!
//! petekStatic **owns** this neutral type; the sampler (petekSim, or petekStatic's
//! own MC driver in `task_peteksim_mc_structured`) *fills* it per draw — no sampler
//! dependency crosses down into the geomodel layer (the one-way DAG; Q3 of the
//! ratified regeneration seam).
//!
//! Ratified shape (graph `decision_staticmodel_regen_seam`): `#[non_exhaustive]` +
//! a [`RealizationDraw::new`] constructor (so P5 structural fields stay additive
//! without re-ratifying the seam), a concrete `pub seed_index: u64`, derives
//! `Clone + Debug`, and a concretely-typed structural [`Option`] (empty for the
//! MVP scalar draw). **`fvf` is deliberately EXCLUDED** — PVT enters volumetrics as
//! a separate uncertain scalar the facade supplies; it never rides the draw.

use petektools::Variogram;

/// A per-realization structural perturbation of the template's control lattice.
/// Empty for the MVP scalar-only draw; the concrete shape lets P5 structural
/// sampling land additively (the `#[non_exhaustive]` draw need not change).
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct StructuralPerturbation {
    /// Per-node top-surface depth shifts this realization: `(ip, jp, dz_m)`.
    pub control_shifts: Vec<(usize, usize, f64)>,
}

/// A **correlated structural perturbation field** for one horizon/isochore of a
/// Monte-Carlo draw (`decision_structural_uncertainty_isochore`): an unconditional
/// Gaussian random field (petekTools `sgs_unconditional`) with marginal
/// `N(0, sd_m²)` and the spatial continuity of `variogram`, generated on the areal
/// node lattice at [`crate::StaticModelTemplate::realize`] time and added to the
/// surface (a TOP **depth** field) or the zone's **thickness** (a deeper-horizon
/// **isochore** field). Because every node's marginal is `N(0, sd_m²)` regardless
/// of the correlation, the **mean** perturbation is variogram-independent — only
/// the field's *shape* (range) depends on the variogram.
///
/// The RNG seed derives from [`RealizationDraw::seed_index`] salted by the horizon
/// index, so a field is **bit-reproducible per seed** and mutually independent
/// across horizons; the SGS search neighbourhood is derived from the variogram
/// range (the same bounded default the property pipelines use). `sd_m <= 0` is a
/// no-op (zero field). Rides the `#[non_exhaustive]` [`RealizationDraw`] /
/// [`ZoneDraw`] additively.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PerturbationField {
    /// Standard deviation \[m\] of the perturbation (a TOP-surface **depth** sd, or a
    /// zone **thickness/isochore** sd). `<= 0` disables the field.
    pub sd_m: f64,
    /// Spatial-continuity model of the field (the `range` is metres; the sill is
    /// irrelevant — the field is rescaled to `sd_m`, only the shape/range matter).
    pub variogram: Variogram,
}

impl PerturbationField {
    /// A perturbation field of magnitude `sd_m` \[m\] with the spatial continuity of
    /// `variogram`.
    #[must_use]
    pub fn new(sd_m: f64, variogram: Variogram) -> Self {
        Self { sd_m, variogram }
    }
}

/// A per-realization draw for **one zone** of a stack-aware
/// [`crate::StaticModelTemplate::from_horizon_stack`] MC (P8,
/// `task_petekstatic_multizone_2`): this zone's optional fluid contacts and optional
/// per-zone property-level overrides. A zone absent from
/// [`RealizationDraw::zones`] uses the template's static per-zone contacts and the
/// draw's base priors; a zone with **no contacts** (neither here nor on the
/// template) contributes GRV but **zero hydrocarbon** in-place.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct ZoneDraw {
    /// Zone index (into the stack's `zone_layers`, top→down).
    pub zone: usize,
    /// This zone's gas–oil contact \[m, positive-down\]; `Some` + a lower OWC makes
    /// it a two-contact (gas-cap + oil-leg) zone. Must be shallower than `owc_depth_m`.
    pub goc_depth_m: Option<f64>,
    /// This zone's oil/gas–water contact \[m, positive-down\]; `None` (with no GOC)
    /// = a contactless zone (GRV, zero hydrocarbon).
    pub owc_depth_m: Option<f64>,
    /// Per-zone porosity level override (fraction); `None` = the draw's base prior.
    pub porosity: Option<f64>,
    /// Per-zone net-to-gross level override (fraction); `None` = the base prior.
    pub net_to_gross: Option<f64>,
    /// Per-zone water-saturation level override (fraction); `None` = the base prior.
    pub water_saturation: Option<f64>,
    /// Optional **isochore (thickness) perturbation field** for this zone
    /// (`decision_structural_uncertainty_isochore`): a correlated `N(0, sd²)` field
    /// added to the zone's base thickness this draw, clamped `>= 0` and **zero-masked
    /// where the base isochore is exactly 0** (a merged zone stays merged in every
    /// draw). Because the deeper horizon of the stack is `top + Σ isochores`, this is
    /// how a **deeper horizon** perturbs structurally. `None` = the zone's thickness
    /// is draw-invariant. Stack-aware templates only.
    pub isochore_structural: Option<PerturbationField>,
}

impl ZoneDraw {
    /// A zone draw for zone index `zone` with no contacts / overrides (a contactless
    /// zone at the base priors). Add contacts / levels with the setters.
    #[must_use]
    pub fn new(zone: usize) -> Self {
        Self {
            zone,
            goc_depth_m: None,
            owc_depth_m: None,
            porosity: None,
            net_to_gross: None,
            water_saturation: None,
            isochore_structural: None,
        }
    }

    /// Attach this zone's **isochore (thickness) perturbation field**
    /// (`decision_structural_uncertainty_isochore`) — the correlated structural
    /// uncertainty of the zone's DEEPER bounding horizon, applied in thickness space
    /// so ordering and exact merges survive every draw by construction. See
    /// [`PerturbationField`].
    #[must_use]
    pub fn with_isochore_structural(mut self, field: PerturbationField) -> Self {
        self.isochore_structural = Some(field);
        self
    }

    /// Set this zone's OWC/GWC (lower) contact depth \[m, positive-down\].
    #[must_use]
    pub fn with_owc(mut self, owc_depth_m: f64) -> Self {
        self.owc_depth_m = Some(owc_depth_m);
        self
    }

    /// Set this zone's GOC (upper) contact depth \[m\] — makes it a two-contact zone
    /// (needs an OWC below).
    #[must_use]
    pub fn with_goc(mut self, goc_depth_m: f64) -> Self {
        self.goc_depth_m = Some(goc_depth_m);
        self
    }

    /// Set this zone's per-zone property levels (fractions) — the per-zone
    /// distribution shift the MC varies.
    #[must_use]
    pub fn with_priors(mut self, porosity: f64, net_to_gross: f64, water_saturation: f64) -> Self {
        self.porosity = Some(porosity);
        self.net_to_gross = Some(net_to_gross);
        self.water_saturation = Some(water_saturation);
        self
    }
}

/// The per-realization sampled scalar set (+ optional structural perturbation).
///
/// `#[non_exhaustive]`: construct via [`RealizationDraw::new`] and set optional
/// fields with the `with_*` builders; new structural fields are additive.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct RealizationDraw {
    /// Areal footprint [m²] (sets cell spacing).
    pub area_m2: f64,
    /// Gross column thickness [m] (top→base offset).
    pub gross_height_m: f64,
    /// Hydrocarbon contact depth [m, positive down] this realization. When
    /// `goc_depth_m` is set this is the lower (OWC/FWL) contact; otherwise it is
    /// the single hydrocarbon contact.
    pub contact_depth_m: f64,
    /// Optional gas–oil contact depth [m, positive down] for a gas-cap + oil-rim
    /// realization. `Some(_)` makes the realized model a two-contact column
    /// (gas cap above `goc_depth_m`, oil leg down to `contact_depth_m`); `None`
    /// = a single-contact column (the MVP). Must be shallower than
    /// `contact_depth_m`.
    pub goc_depth_m: Option<f64>,
    /// Porosity prior (fraction).
    pub porosity: f64,
    /// Net-to-gross prior (fraction).
    pub net_to_gross: f64,
    /// Water saturation prior (fraction).
    pub water_saturation: f64,
    /// The RNG seed index that produced this draw — carried into provenance so a
    /// realization is reproducible/auditable.
    pub seed_index: u64,
    /// Optional structural perturbation of the control lattice (`None` = scalar
    /// draw only, the MVP).
    pub structural: Option<StructuralPerturbation>,
    /// Optional **correlated TOP-surface depth perturbation field**
    /// (`decision_structural_uncertainty_isochore`): a `N(0, sd²)` field added to the
    /// top surface this draw (both the 2-surface and stack paths). The stack's deeper
    /// horizons then ride this perturbed top plus their own (optionally perturbed)
    /// isochores ([`ZoneDraw::isochore_structural`]). `None` = the top is
    /// draw-invariant (aside from the MVP `structural` control shifts, if any).
    pub top_structural: Option<PerturbationField>,
    /// Optional gas-cap connate-water override (R3): with a drawn GOC, gas-zone
    /// cells use this single scalar for `(1 - Sw)` instead of the shared `SW`
    /// cube, so a single cube does not over-state gas-cap OGIP. `None` = the gas
    /// cap uses the cube. No effect on a single-contact (no-GOC) draw.
    pub sw_gas: Option<f64>,
    /// Per-property additive **level shifts** for [`crate::McMode::LevelShift`]
    /// geostatistical properties (`decision_mc_composition`): `(property, delta)`.
    /// The template adds `delta` to every cell of that property's once-propagated
    /// pattern, so the spatial pattern is preserved and only the level moves.
    /// Absent property = no shift. Ignored for [`crate::McMode::Resimulate`]
    /// properties (they redraw the whole field).
    pub property_shifts: Vec<(String, f64)>,
    /// Per-zone draws for a stack-aware
    /// [`crate::StaticModelTemplate::from_horizon_stack`] MC (P8): each zone's
    /// optional contacts + per-zone property levels. Empty on the 2-surface path
    /// (that draw uses `contact_depth_m` / `goc_depth_m` + the base priors).
    pub zones: Vec<ZoneDraw>,
}

impl RealizationDraw {
    /// A scalar draw (no structural perturbation). The seven load-bearing scalars
    /// the MVP MC loop varies; `structural` defaults to `None`.
    #[must_use]
    pub fn new(
        area_m2: f64,
        gross_height_m: f64,
        contact_depth_m: f64,
        porosity: f64,
        net_to_gross: f64,
        water_saturation: f64,
        seed_index: u64,
    ) -> Self {
        Self {
            area_m2,
            gross_height_m,
            contact_depth_m,
            goc_depth_m: None,
            porosity,
            net_to_gross,
            water_saturation,
            seed_index,
            structural: None,
            top_structural: None,
            sw_gas: None,
            property_shifts: Vec::new(),
            zones: Vec::new(),
        }
    }

    /// Attach a per-zone draw for a stack-aware MC (P8). Replaces any prior draw for
    /// the same zone index. See [`ZoneDraw`] and
    /// [`crate::StaticModelTemplate::from_horizon_stack`].
    #[must_use]
    pub fn with_zone_draw(mut self, zone: ZoneDraw) -> Self {
        self.zones.retain(|z| z.zone != zone.zone);
        self.zones.push(zone);
        self
    }

    /// Set an additive level shift for a [`crate::McMode::LevelShift`] property
    /// (`decision_mc_composition`): the template adds `delta` to every cell of that
    /// property's once-propagated pattern. Replaces any prior shift for the same
    /// property. No effect on a [`crate::McMode::Resimulate`] property.
    #[must_use]
    pub fn with_property_shift(mut self, property: impl Into<String>, delta: f64) -> Self {
        let property = property.into();
        self.property_shifts.retain(|(p, _)| p != &property);
        self.property_shifts.push((property, delta));
        self
    }

    /// The additive level shift for `property` (0.0 if none set).
    #[must_use]
    pub fn property_shift(&self, property: &str) -> f64 {
        self.property_shifts
            .iter()
            .find(|(p, _)| p == property)
            .map_or(0.0, |(_, d)| *d)
    }

    /// Attach a structural perturbation (P5 additive path).
    #[must_use]
    pub fn with_structural(mut self, perturbation: StructuralPerturbation) -> Self {
        self.structural = Some(perturbation);
        self
    }

    /// Attach a correlated **TOP-surface depth perturbation field**
    /// (`decision_structural_uncertainty_isochore`). See [`PerturbationField`].
    #[must_use]
    pub fn with_top_structural(mut self, field: PerturbationField) -> Self {
        self.top_structural = Some(field);
        self
    }

    /// Set a gas-cap connate-water override for this draw (R3): with a drawn GOC,
    /// the gas cap uses `sw_gas` instead of the shared `SW` cube. No effect
    /// without a GOC (single-contact column).
    #[must_use]
    pub fn with_sw_gas(mut self, sw_gas: f64) -> Self {
        self.sw_gas = Some(sw_gas);
        self
    }

    /// Set a gas–oil contact depth, making this a gas-cap + oil-rim realization
    /// (two-contact column). `goc_depth_m` must be shallower than
    /// `contact_depth_m` (validated at [`crate::StaticModelTemplate::realize`]).
    #[must_use]
    pub fn with_goc(mut self, goc_depth_m: f64) -> Self {
        self.goc_depth_m = Some(goc_depth_m);
        self
    }

    // --- named scalar setters (V8/R6 ergonomics) ---
    // Fluent overrides so a sampler can build a draw field-by-field from a base
    // rather than positionally through `::new`; each returns a new draw. `::new`
    // stays for compat and the load-bearing seven-scalar path.

    /// Set the areal footprint [m²].
    #[must_use]
    pub fn with_area(mut self, area_m2: f64) -> Self {
        self.area_m2 = area_m2;
        self
    }

    /// Set the gross column thickness / base-relief level [m].
    #[must_use]
    pub fn with_gross(mut self, gross_height_m: f64) -> Self {
        self.gross_height_m = gross_height_m;
        self
    }

    /// Set the (lower, OWC/FWL when a GOC is present) hydrocarbon contact [m].
    #[must_use]
    pub fn with_contact(mut self, contact_depth_m: f64) -> Self {
        self.contact_depth_m = contact_depth_m;
        self
    }

    /// Set the porosity prior (fraction).
    #[must_use]
    pub fn with_porosity(mut self, porosity: f64) -> Self {
        self.porosity = porosity;
        self
    }

    /// Set the net-to-gross prior (fraction).
    #[must_use]
    pub fn with_ntg(mut self, net_to_gross: f64) -> Self {
        self.net_to_gross = net_to_gross;
        self
    }

    /// Set the water-saturation prior (fraction).
    #[must_use]
    pub fn with_sw(mut self, water_saturation: f64) -> Self {
        self.water_saturation = water_saturation;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn named_setters_override_only_their_field() {
        let base = RealizationDraw::new(100.0, 50.0, 5000.0, 0.25, 0.8, 0.3, 7);
        let d = base
            .clone()
            .with_area(200.0)
            .with_gross(60.0)
            .with_contact(5100.0)
            .with_porosity(0.30)
            .with_ntg(0.9)
            .with_sw(0.2);
        assert_eq!(
            (d.area_m2, d.gross_height_m, d.contact_depth_m),
            (200.0, 60.0, 5100.0)
        );
        assert_eq!(
            (d.porosity, d.net_to_gross, d.water_saturation),
            (0.30, 0.9, 0.2)
        );
        // Untouched fields carry through from the base.
        assert_eq!(d.seed_index, 7);
        assert_eq!(d.goc_depth_m, None);
        assert_eq!(d.structural, None);
    }
}
