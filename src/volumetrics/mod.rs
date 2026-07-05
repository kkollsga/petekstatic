//! `srs-volumetrics` — GRV, in-place volumetrics and FVF from the populated grid.
//!
//! This is the static-modelling layer's **volumetrics** half (graph
//! `decision_layer_charters`): GRV/in-place over a `StaticModel`'s grid + cubes,
//! with FVF entering as a validated scalar input (no PVT crosses the seam).
//!
//! [`compute_in_place`] is the deterministic answer: per-cell HCPV summed over
//! the hydrocarbon column → GRV, STOIIP, GIIP (`grv_from_grid_spec`).
//! [`compute_in_place_two_contact`] splits a gas cap + oil rim (GOC + OWC/FWL)
//! into separate gas-zone and oil-zone volumes (geometry only; no PVT).
//! [`populate_constant`] seeds the Day-1 prior properties (the no-data case of
//! `grid_property_population_spec`). [`OilFvf`]/[`GasFvf`] are the self-validating
//! FVF value types the in-place conversion consumes.

mod backing;
mod fvf;
mod grv;
mod names;
mod populate;
mod valid;

pub use backing::{compute_clipped, CellSlab, Clip, GridSource, SlabSource};
pub use fvf::{GasFvf, OilFvf};
pub use grv::{
    compute_in_place, compute_in_place_summary, compute_in_place_two_contact,
    compute_in_place_two_contact_summary, compute_in_place_two_contact_zone, compute_in_place_zone,
    compute_zone_bulk, InPlace, ZoneVolumes,
};
pub use names::{NTG, PORO, SW};
pub use populate::{populate_constant, ConstantPriors};
pub use valid::{validate_fraction, validate_positive};
