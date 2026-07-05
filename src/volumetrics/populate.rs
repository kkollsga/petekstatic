//! Prior-only (constant) property population — the no-data degenerate of
//! `grid_property_population_spec`. With no wells, kriging-with-trend reduces to
//! the trend, and a flat prior trend is a single constant per property. The
//! geostatistical case (variograms, kriging, SGS) lands once wells/control
//! points exist (`srs-petro` + live-refine).

use crate::error::StaticError;
use crate::grid::{Grid, Property};
use crate::volumetrics::names::{NTG, PORO, SW};
use crate::volumetrics::valid::validate_fraction;
use serde::{Deserialize, Serialize};

/// Day-1 priors: a single representative value per property (fractions).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ConstantPriors {
    pub porosity: f64,
    pub net_to_gross: f64,
    pub water_saturation: f64,
}

/// Populate `PORO`, `NTG`, `SW` as constants across every cell.
///
/// # Errors
/// Returns [`StaticError::InvalidInput`] if any prior is not a fraction in `[0, 1]`.
pub fn populate_constant(grid: &mut Grid, priors: ConstantPriors) -> Result<(), StaticError> {
    validate_fraction("porosity", priors.porosity)?;
    validate_fraction("net_to_gross", priors.net_to_gross)?;
    validate_fraction("water_saturation", priors.water_saturation)?;

    let n = grid.cell_count();
    let props = grid.properties_mut();
    // Recycle each cube's existing value buffer in place (take → refill → reinstall)
    // so the steady-state `realize_into` MC path allocates no fresh cube; on a cold
    // build the take returns an empty Vec and the resize allocates once — identical
    // outcome, one code path.
    for (name, value) in [
        (PORO, priors.porosity),
        (NTG, priors.net_to_gross),
        (SW, priors.water_saturation),
    ] {
        let mut values = props.take_values(name);
        values.clear();
        values.resize(n, value);
        props.set(Property {
            name: name.to_string(),
            values,
        })?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::{build_box, BoxSpec, Dims};

    #[test]
    fn populates_all_three_properties() {
        let mut grid = build_box(BoxSpec::square(40.0, 20.0, Dims::new(2, 2, 2).unwrap())).unwrap();
        populate_constant(
            &mut grid,
            ConstantPriors {
                porosity: 0.2,
                net_to_gross: 0.8,
                water_saturation: 0.3,
            },
        )
        .unwrap();
        assert_eq!(grid.properties().len(), 3);
        assert!((grid.properties().get(PORO).unwrap().values[0] - 0.2).abs() < 1e-12);
    }

    #[test]
    fn rejects_out_of_range_prior() {
        let mut grid = build_box(BoxSpec::square(40.0, 20.0, Dims::new(1, 1, 1).unwrap())).unwrap();
        let bad = ConstantPriors {
            porosity: 1.5,
            net_to_gross: 0.8,
            water_saturation: 0.3,
        };
        assert!(populate_constant(&mut grid, bad).is_err());
    }
}
