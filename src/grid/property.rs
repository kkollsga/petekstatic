//! Cell properties: named per-cell scalar arrays (porosity, Sw, net-to-gross,
//! ...). Each property has one value per cell, indexed by the grid's linear
//! cell index.

use crate::error::StaticError;
use std::collections::HashMap;

/// A named per-cell scalar field. `values.len()` equals the grid cell count.
#[derive(Debug, Clone)]
pub struct Property {
    pub name: String,
    pub values: Vec<f64>,
}

impl Property {
    /// A property filled with a single constant value.
    #[must_use]
    pub fn constant(name: impl Into<String>, value: f64, cell_count: usize) -> Self {
        Self {
            name: name.into(),
            values: vec![value; cell_count],
        }
    }
}

/// The set of properties attached to a grid, keyed by name.
#[derive(Debug, Clone, Default)]
pub struct Properties {
    cell_count: usize,
    map: HashMap<String, Property>,
}

impl Properties {
    /// An empty property set for a grid of `cell_count` cells.
    #[must_use]
    pub fn new(cell_count: usize) -> Self {
        Self {
            cell_count,
            map: HashMap::new(),
        }
    }

    /// Insert or replace a property.
    ///
    /// # Errors
    /// Returns [`StaticError::InvalidInput`] if its length != the grid cell count.
    pub fn set(&mut self, prop: Property) -> Result<(), StaticError> {
        if prop.values.len() != self.cell_count {
            return Err(StaticError::InvalidInput(format!(
                "property '{}' has {} values, expected {}",
                prop.name,
                prop.values.len(),
                self.cell_count
            )));
        }
        self.map.insert(prop.name.clone(), prop);
        Ok(())
    }

    /// Borrow a property by name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Property> {
        self.map.get(name)
    }

    /// Take a property's value buffer **out** for in-place refill, removing the
    /// property and returning its `Vec` (capacity retained), or a fresh empty
    /// `Vec` when the property is absent. Pair with [`Properties::set`] after
    /// overwriting the buffer to recycle the cube allocation across a Monte-Carlo
    /// draw (`StaticModelTemplate::realize_into`) — no fresh cube alloc on the
    /// steady-state path. The caller **must** refill it to the current cell count.
    #[must_use]
    pub fn take_values(&mut self, name: &str) -> Vec<f64> {
        self.map.remove(name).map(|p| p.values).unwrap_or_default()
    }

    /// Reset the cell count every cube must match — the geometry was rebuilt with
    /// a new layer count. Existing cubes are **not** resized here; the following
    /// populate refills them to the new length (see [`Properties::take_values`]).
    /// Crate-internal: only the geometry-recycling path
    /// ([`crate::grid::Grid::install_geometry`]) calls it.
    pub(crate) fn set_cell_count(&mut self, cell_count: usize) {
        self.cell_count = cell_count;
    }

    /// Number of properties stored.
    #[must_use]
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Whether no properties are stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// The stored property names (unordered).
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.map.keys().map(String::as_str)
    }
}
