//! `srs-petro` — petrophysics. The MVP delivers **log upscaling**
//! (`log_upscaling_spec`): length-weighted power-law means with each property
//! weighted by what it conserves (porosity by length, Sw by pore volume, NTG by
//! net length, k bracketed by Cardwell-Parsons bounds). Facies blocking and
//! flow-based k upscaling are later work.

mod power_mean;
mod upscale;

pub use power_mean::{
    arithmetic_mean, geometric_mean, harmonic_mean, power_law_mean, WeightedSample,
};
pub use upscale::{net_to_gross, perm_bounds, upscale_porosity, upscale_sw, NetSample, SwSample};
