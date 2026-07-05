//! Grid dimensions and the `(i,j,k)` <-> linear cell index mapping.
//!
//! `i` runs fastest, then `j`, then `k` (k = layer, increasing downward).

use petekstatic_error::StaticError;

/// Grid dimensions: cell counts along i, j, k.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Dims {
    pub ni: usize,
    pub nj: usize,
    pub nk: usize,
}

/// An `(i,j,k)` cell index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ijk {
    pub i: usize,
    pub j: usize,
    pub k: usize,
}

impl Ijk {
    /// Construct an index.
    #[must_use]
    pub const fn new(i: usize, j: usize, k: usize) -> Self {
        Self { i, j, k }
    }
}

impl Dims {
    /// Construct dimensions; all axes must be non-zero.
    ///
    /// # Errors
    /// Returns [`StaticError::Grid`] if any dimension is zero.
    pub fn new(ni: usize, nj: usize, nk: usize) -> Result<Self, StaticError> {
        if ni == 0 || nj == 0 || nk == 0 {
            return Err(StaticError::Grid(format!(
                "dimensions must be non-zero, got ({ni},{nj},{nk})"
            )));
        }
        Ok(Self { ni, nj, nk })
    }

    /// Total number of cells.
    #[must_use]
    pub const fn cell_count(self) -> usize {
        self.ni * self.nj * self.nk
    }

    /// Number of pillars in the areal lattice: `(ni+1) * (nj+1)`.
    #[must_use]
    pub const fn pillar_count(self) -> usize {
        (self.ni + 1) * (self.nj + 1)
    }

    /// Linear pillar index for areal pillar `(ip, jp)`, `ip in 0..=ni`.
    #[must_use]
    #[inline]
    pub const fn pillar_linear(self, ip: usize, jp: usize) -> usize {
        jp * (self.ni + 1) + ip
    }

    /// `(i,j,k)` -> linear cell index, or `None` if out of bounds.
    #[must_use]
    #[inline]
    pub fn linear(self, c: Ijk) -> Option<usize> {
        if c.i >= self.ni || c.j >= self.nj || c.k >= self.nk {
            return None;
        }
        Some((c.k * self.nj + c.j) * self.ni + c.i)
    }

    /// Iterate all cell indices in linear order (`i` fastest).
    pub fn iter(self) -> impl Iterator<Item = Ijk> {
        (0..self.nk).flat_map(move |k| {
            (0..self.nj).flat_map(move |j| (0..self.ni).map(move |i| Ijk::new(i, j, k)))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_dimension_rejected() {
        assert!(Dims::new(0, 1, 1).is_err());
        assert!(Dims::new(2, 3, 4).is_ok());
    }

    #[test]
    fn linear_indexing_is_dense_and_ordered() {
        let d = Dims::new(2, 3, 4).unwrap();
        assert_eq!(d.cell_count(), 24);
        let all: Vec<usize> = d.iter().map(|c| d.linear(c).unwrap()).collect();
        assert_eq!(all, (0..24).collect::<Vec<_>>());
    }

    #[test]
    fn out_of_bounds_is_none() {
        let d = Dims::new(2, 2, 2).unwrap();
        assert_eq!(d.linear(Ijk::new(2, 0, 0)), None);
        assert_eq!(d.linear(Ijk::new(1, 1, 1)), Some(7));
    }

    #[test]
    fn pillar_lattice_sized_correctly() {
        let d = Dims::new(2, 3, 1).unwrap();
        assert_eq!(d.pillar_count(), 3 * 4);
        assert_eq!(d.pillar_linear(2, 3), 3 * 3 + 2);
    }
}
