//! Minimum-curvature surface interpolation.
//!
//! Solve `z(x,y)` on the areal pillar lattice minimizing total squared curvature
//! subject to hard control points: biharmonic `∇⁴z = 0` in the gaps (Briggs
//! 1974), blended toward harmonic `∇²z = 0` by a tension `T in [0,1]` (Smith &
//! Wessel 1990). Iterated by SOR to a constraint tolerance.
//!
//! Center-node update solving `(1-T)∇⁴z - T∇²z = 0` (uniform spacing):
//! `z0 = [ (1-T)(8·E1 - 2·D - W2) + T·E1 ] / [ 20(1-T) + 4T ]`,
//! where `E1` = 4 edge neighbours, `D` = 4 diagonals, `W2` = 4 two-away nodes.
//! Out-of-lattice stencil nodes are synthesized by linear extrapolation (the
//! natural minimum-curvature boundary condition), so a planar regional-dip field
//! is an exact fixed point everywhere — the gridder honours a regional trend
//! rather than flattening it at the edges.
//!
//! Caveats pending validation: the boundary condition, tension discretization
//! and convergence tolerance are implementation choices that may need tuning
//! against structured/faulted cases.

use petekstatic_error::StaticError;
use petektools::{gridding::grid_min_curvature_seeded, Lattice};
use serde::{Deserialize, Serialize};

/// A control point pinned on areal node `(ip, jp)` to depth `z`.
#[derive(Debug, Clone, Copy)]
pub struct Control {
    pub ip: usize,
    pub jp: usize,
    pub z: f64,
}

/// Settings for the relaxation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SolveOpts {
    /// Tension in `[0, 1]` (0 = pure min-curvature, 1 = harmonic).
    pub tension: f64,
    /// SOR over-relaxation factor in `(0, 2)`.
    pub omega: f64,
    /// Stop when the max nodal change in a sweep falls below this.
    pub tol: f64,
    /// Hard iteration cap.
    pub max_iter: usize,
}

impl Default for SolveOpts {
    fn default() -> Self {
        Self {
            tension: 0.25,
            omega: 1.5,
            tol: 1e-6,
            max_iter: 20_000,
        }
    }
}

/// A solved surface: depth `z` at every node of a `(ni+1) x (nj+1)` lattice.
#[derive(Debug, Clone)]
pub struct Surface {
    nx: usize, // ni + 1
    ny: usize, // nj + 1
    z: Vec<f64>,
}

impl Surface {
    /// A flat surface at constant depth `z` over a `nx x ny` node lattice.
    ///
    /// # Panics
    /// Panics if `nx < 2` or `ny < 2`.
    #[must_use]
    pub fn constant(nx: usize, ny: usize, z: f64) -> Self {
        assert!(nx >= 2 && ny >= 2, "surface lattice must be at least 2x2");
        Self {
            nx,
            ny,
            z: vec![z; nx * ny],
        }
    }

    /// A copy shifted down by `dz` (a conformable base at constant thickness).
    #[must_use]
    pub fn offset_by(&self, dz: f64) -> Self {
        Self {
            nx: self.nx,
            ny: self.ny,
            z: self.z.iter().map(|v| v + dz).collect(),
        }
    }

    /// A copy shifted down by a **per-node** field `dz_m` (metres, row-major
    /// `jp * nx + ip`, one value per lattice node) — a base surface that follows
    /// real relief rather than a constant offset
    /// (`decision_template_gross_scaling`).
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if `dz_m.len() != nx * ny` or a value is
    /// not finite.
    pub fn offset_by_field(&self, dz_m: &[f64]) -> Result<Self, StaticError> {
        if dz_m.len() != self.nx * self.ny {
            return Err(StaticError::InvalidInput(format!(
                "offset field has {} values, expected {} ({}x{})",
                dz_m.len(),
                self.nx * self.ny,
                self.nx,
                self.ny
            )));
        }
        if let Some(bad) = dz_m.iter().find(|v| !v.is_finite()) {
            return Err(StaticError::InvalidInput(format!(
                "offset field must be finite, got {bad}"
            )));
        }
        Ok(Self {
            nx: self.nx,
            ny: self.ny,
            z: self.z.iter().zip(dz_m).map(|(v, d)| v + d).collect(),
        })
    }

    /// Areal node count along i (= ni + 1).
    #[must_use]
    pub fn nx(&self) -> usize {
        self.nx
    }
    /// Areal node count along j (= nj + 1).
    #[must_use]
    pub fn ny(&self) -> usize {
        self.ny
    }
    /// Depth at node `(ip, jp)`.
    #[must_use]
    pub fn z(&self, ip: usize, jp: usize) -> f64 {
        self.z[jp * self.nx + ip]
    }

    /// Guard that `self` (a base surface) sits at or below `top` at every node —
    /// depth is positive-down, so a valid base has `base_z >= top_z` (non-negative
    /// gross). A crossing (`base_z < top_z`) is a thin/crossing framework that
    /// would silently collapse GRV.
    ///
    /// With `clamp` the crossings are pulled to the top (zero-thickness at those
    /// nodes only) and the clamped base is returned. Without it, any crossing is a
    /// [`StaticError::CrossedSurfaces`] reporting the offending node count and the
    /// worst (most negative) separation.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the lattices differ in size;
    /// [`StaticError::CrossedSurfaces`] if `!clamp` and the base crosses the top.
    pub fn guard_below(&self, top: &Surface, clamp: bool) -> Result<Surface, StaticError> {
        if self.nx != top.nx || self.ny != top.ny {
            return Err(StaticError::InvalidInput(format!(
                "base lattice {}x{} does not match top {}x{}",
                self.nx, self.ny, top.nx, top.ny
            )));
        }
        let mut nodes = 0usize;
        let mut worst = 0.0_f64;
        for (b, t) in self.z.iter().zip(&top.z) {
            let sep = b - t;
            if sep < 0.0 {
                nodes += 1;
                worst = worst.min(sep);
            }
        }
        if nodes == 0 {
            return Ok(self.clone());
        }
        if !clamp {
            return Err(StaticError::CrossedSurfaces {
                nodes,
                worst_m: worst,
            });
        }
        Ok(Surface {
            nx: self.nx,
            ny: self.ny,
            z: self.z.iter().zip(&top.z).map(|(b, t)| b.max(*t)).collect(),
        })
    }

    /// Post-gridding order-repair: where this base sits less than `min_thickness_m`
    /// below `top` (thickness `base_z - top_z < min_thickness_m`, including a true
    /// crossing), pull the base **down** to exactly `top_z + min_thickness_m`,
    /// **preserving the top** — the top is the better-constrained seismic pick, so
    /// the (softer) base yields to it. Independent gridding of Top and Base can
    /// overshoot at thin margins and re-introduce a crossing a pointwise
    /// pre-repair had removed; this repairs the gridded result.
    ///
    /// Returns the repaired base plus `(repaired, worst_m)`: how many nodes were
    /// pushed and the worst (smallest, most negative → a true crossing) original
    /// `base_z - top_z` separation among them (`0` / `0.0` when nothing crossed).
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the lattices differ in size or
    /// `min_thickness_m` is not finite and `>= 0`.
    pub fn repair_min_thickness(
        &self,
        top: &Surface,
        min_thickness_m: f64,
    ) -> Result<(Surface, usize, f64), StaticError> {
        if self.nx != top.nx || self.ny != top.ny {
            return Err(StaticError::InvalidInput(format!(
                "base lattice {}x{} does not match top {}x{}",
                self.nx, self.ny, top.nx, top.ny
            )));
        }
        if !(min_thickness_m.is_finite() && min_thickness_m >= 0.0) {
            return Err(StaticError::InvalidInput(format!(
                "min_thickness_m must be finite and >= 0, got {min_thickness_m}"
            )));
        }
        let mut repaired = 0usize;
        let mut worst = 0.0_f64;
        let z: Vec<f64> = self
            .z
            .iter()
            .zip(&top.z)
            .map(|(b, t)| {
                let floor = t + min_thickness_m;
                if *b < floor {
                    repaired += 1;
                    worst = worst.min(b - t);
                    floor
                } else {
                    *b
                }
            })
            .collect();
        Ok((
            Surface {
                nx: self.nx,
                ny: self.ny,
                z,
            },
            repaired,
            worst,
        ))
    }

    /// Repair-precedence twin of [`Self::guard_below`] that moves the **upper**
    /// surface instead of the lower — used when the upper is a *derived* horizon
    /// (a tops-only drape) and `lower` is *mapped*: the mapped surface is
    /// authoritative, so the derived one yields. Ensures `self` (the upper) sits at
    /// or above `lower` at every node (`self_z <= lower_z`, positive-down). A
    /// crossing (`self_z > lower_z`) with `clamp` pulls the upper **up** to the
    /// mapped lower (zero-thickness at those nodes); without `clamp` it is a
    /// [`StaticError::CrossedSurfaces`].
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the lattices differ in size;
    /// [`StaticError::CrossedSurfaces`] if `!clamp` and the upper crosses the lower.
    pub fn guard_above(&self, lower: &Surface, clamp: bool) -> Result<Surface, StaticError> {
        if self.nx != lower.nx || self.ny != lower.ny {
            return Err(StaticError::InvalidInput(format!(
                "upper lattice {}x{} does not match lower {}x{}",
                self.nx, self.ny, lower.nx, lower.ny
            )));
        }
        let mut nodes = 0usize;
        let mut worst = 0.0_f64;
        for (u, l) in self.z.iter().zip(&lower.z) {
            let sep = l - u;
            if sep < 0.0 {
                nodes += 1;
                worst = worst.min(sep);
            }
        }
        if nodes == 0 {
            return Ok(self.clone());
        }
        if !clamp {
            return Err(StaticError::CrossedSurfaces {
                nodes,
                worst_m: worst,
            });
        }
        Ok(Surface {
            nx: self.nx,
            ny: self.ny,
            z: self
                .z
                .iter()
                .zip(&lower.z)
                .map(|(u, l)| u.min(*l))
                .collect(),
        })
    }

    /// Repair-precedence twin of [`Self::repair_min_thickness`] that moves the
    /// **upper** surface up instead of the lower down — the derived-yields-to-mapped
    /// case (see [`Self::guard_above`]). Where this upper sits less than
    /// `min_thickness_m` above `lower` (`lower_z - self_z < min_thickness_m`,
    /// including a crossing), pull the upper **up** to exactly `lower_z -
    /// min_thickness_m`, preserving the (mapped) lower. Returns
    /// `(repaired, nodes, worst)` — `worst` the smallest (most negative) original
    /// `lower_z - self_z`.
    ///
    /// # Errors
    /// [`StaticError::InvalidInput`] if the lattices differ in size or
    /// `min_thickness_m` is not finite and `>= 0`.
    pub fn repair_min_thickness_from_below(
        &self,
        lower: &Surface,
        min_thickness_m: f64,
    ) -> Result<(Surface, usize, f64), StaticError> {
        if self.nx != lower.nx || self.ny != lower.ny {
            return Err(StaticError::InvalidInput(format!(
                "upper lattice {}x{} does not match lower {}x{}",
                self.nx, self.ny, lower.nx, lower.ny
            )));
        }
        if !(min_thickness_m.is_finite() && min_thickness_m >= 0.0) {
            return Err(StaticError::InvalidInput(format!(
                "min_thickness_m must be finite and >= 0, got {min_thickness_m}"
            )));
        }
        let mut repaired = 0usize;
        let mut worst = 0.0_f64;
        let z: Vec<f64> = self
            .z
            .iter()
            .zip(&lower.z)
            .map(|(u, l)| {
                let ceil = l - min_thickness_m;
                if *u > ceil {
                    repaired += 1;
                    worst = worst.min(l - u);
                    ceil
                } else {
                    *u
                }
            })
            .collect();
        Ok((
            Surface {
                nx: self.nx,
                ny: self.ny,
                z,
            },
            repaired,
            worst,
        ))
    }

    /// Apply an [`ExtrapolationPolicy`] to this solved surface given the `controls`
    /// (data nodes) that conditioned it: nodes at or within `start_cells` of the
    /// nearest datum keep the solve untouched (data nodes themselves are always
    /// bit-unchanged — they are hard controls); beyond that, the solve blends
    /// linearly toward the **nearest datum's value** over `decay_cells`, so far
    /// extrapolation flattens to nearest-data instead of running the kernel's
    /// natural-dip (linear) extension unbounded into a data void. Distances are in
    /// lattice **cells** (node units). [`ExtrapolationPolicy::NaturalDip`] returns
    /// the surface unchanged (the legacy behaviour). A non-positive `decay_cells`
    /// acts as a hard clamp to nearest-data beyond `start_cells`.
    ///
    /// Nearest-datum assignment uses an 8-neighbour two-pass chamfer sweep
    /// (exact Euclidean distance to the propagated source; near-exact assignment) —
    /// O(nodes), deterministic.
    #[must_use]
    pub fn taper_beyond_data(&self, controls: &[Control], policy: ExtrapolationPolicy) -> Surface {
        let ExtrapolationPolicy::DecayToData {
            start_cells,
            decay_cells,
        } = policy
        else {
            return self.clone();
        };
        if controls.is_empty() {
            return self.clone();
        }
        let (nx, ny) = (self.nx, self.ny);
        // source[node] = index into `controls` of the (approx) nearest datum.
        let mut source: Vec<usize> = vec![usize::MAX; nx * ny];
        let mut dist2: Vec<f64> = vec![f64::INFINITY; nx * ny];
        for (c_idx, c) in controls.iter().enumerate() {
            let idx = c.jp * nx + c.ip;
            source[idx] = c_idx;
            dist2[idx] = 0.0;
        }
        let d2 = |c: &Control, ip: usize, jp: usize| -> f64 {
            let dx = ip as f64 - c.ip as f64;
            let dy = jp as f64 - c.jp as f64;
            dx * dx + dy * dy
        };
        let relax = |ip: usize,
                     jp: usize,
                     nip: usize,
                     njp: usize,
                     src: &mut Vec<usize>,
                     dst: &mut Vec<f64>| {
            let n_idx = njp * nx + nip;
            let s = src[n_idx];
            if s == usize::MAX {
                return;
            }
            let cand = d2(&controls[s], ip, jp);
            let idx = jp * nx + ip;
            if cand < dst[idx] {
                dst[idx] = cand;
                src[idx] = s;
            }
        };
        // Forward pass (W, N, NW, NE), then backward pass (E, S, SE, SW).
        for jp in 0..ny {
            for ip in 0..nx {
                if ip > 0 {
                    relax(ip, jp, ip - 1, jp, &mut source, &mut dist2);
                }
                if jp > 0 {
                    relax(ip, jp, ip, jp - 1, &mut source, &mut dist2);
                    if ip > 0 {
                        relax(ip, jp, ip - 1, jp - 1, &mut source, &mut dist2);
                    }
                    if ip + 1 < nx {
                        relax(ip, jp, ip + 1, jp - 1, &mut source, &mut dist2);
                    }
                }
            }
        }
        for jp in (0..ny).rev() {
            for ip in (0..nx).rev() {
                if ip + 1 < nx {
                    relax(ip, jp, ip + 1, jp, &mut source, &mut dist2);
                }
                if jp + 1 < ny {
                    relax(ip, jp, ip, jp + 1, &mut source, &mut dist2);
                    if ip + 1 < nx {
                        relax(ip, jp, ip + 1, jp + 1, &mut source, &mut dist2);
                    }
                    if ip > 0 {
                        relax(ip, jp, ip - 1, jp + 1, &mut source, &mut dist2);
                    }
                }
            }
        }
        let start = start_cells.max(0.0);
        let z: Vec<f64> = (0..nx * ny)
            .map(|idx| {
                let s = source[idx];
                let d = dist2[idx].sqrt();
                let w = if decay_cells > 0.0 {
                    ((d - start) / decay_cells).clamp(0.0, 1.0)
                } else if d > start {
                    1.0
                } else {
                    0.0
                };
                if w == 0.0 {
                    self.z[idx]
                } else {
                    self.z[idx] * (1.0 - w) + controls[s].z * w
                }
            })
            .collect();
        Surface { nx, ny, z }
    }
}

/// How a solved stack surface (or isochore field) behaves **beyond its data** —
/// the region of the lattice farther from every conditioning datum than the data
/// hull. The kernel's natural-dip boundary linearly extends the local gradient,
/// which is correct near data but runs **unbounded** into a data void (tens of
/// metres over a margin). This policy makes the behaviour explicit and owner-
/// visible (`with_extrapolation` on the builder/template).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ExtrapolationPolicy {
    /// Legacy behaviour: the kernel's natural-dip (linear) extension runs
    /// unbounded beyond the data. Only appropriate when the caller KNOWS the
    /// regional dip continues (e.g. a clipped window of a larger mapped surface).
    NaturalDip,
    /// Conservative default: within `start_cells` (lattice cells) of the nearest
    /// datum the solve is untouched; beyond, it blends linearly toward the nearest
    /// datum's value over `decay_cells`, flattening far extrapolation to
    /// nearest-data. Data nodes themselves are always exact.
    DecayToData {
        /// Distance (cells) from the nearest datum where the taper begins.
        start_cells: f64,
        /// Ramp length (cells) over which the solve blends to nearest-data.
        /// `<= 0` acts as a hard clamp beyond `start_cells`.
        decay_cells: f64,
    },
}

impl Default for ExtrapolationPolicy {
    /// The conservative default: keep the solved dip for 2 cells beyond the data,
    /// then decay to the nearest-data value over the next 4 cells.
    fn default() -> Self {
        ExtrapolationPolicy::DecayToData {
            start_cells: 2.0,
            decay_cells: 4.0,
        }
    }
}

/// Grid a minimum-curvature surface over a `(ni+1) x (nj+1)` node lattice from
/// control points. The surface is seeded at the mean control depth.
///
/// # Errors
/// [`StaticError::InvalidInput`] if the lattice is empty, there are no controls, a
/// control is off-lattice, or options are out of range.
pub fn solve_surface(
    nx: usize,
    ny: usize,
    controls: &[Control],
    opts: SolveOpts,
) -> Result<Surface, StaticError> {
    if nx < 2 || ny < 2 {
        return Err(StaticError::Grid(format!(
            "surface lattice must be at least 2x2 (ni,nj >= 1), got {nx}x{ny}"
        )));
    }
    if controls.is_empty() {
        return Err(StaticError::InvalidInput(
            "minimum-curvature surface needs at least one control point".into(),
        ));
    }
    if !(opts.tension >= 0.0 && opts.tension <= 1.0) {
        return Err(StaticError::OutOfRange(format!(
            "tension must be in [0,1], got {}",
            opts.tension
        )));
    }
    if !(opts.omega > 0.0 && opts.omega < 2.0) {
        return Err(StaticError::OutOfRange(format!(
            "omega must be in (0,2), got {}",
            opts.omega
        )));
    }

    let n = nx * ny;
    let mut fixed = vec![false; n];
    let mut z = vec![0.0; n];
    let mut sum = 0.0;
    for c in controls {
        if c.ip >= nx || c.jp >= ny {
            return Err(StaticError::Grid(format!(
                "control ({},{}) outside lattice {nx}x{ny}",
                c.ip, c.jp
            )));
        }
        let idx = c.jp * nx + c.ip;
        z[idx] = c.z;
        fixed[idx] = true;
        sum += c.z;
    }
    let seed = sum / controls.len() as f64;
    for (idx, zi) in z.iter_mut().enumerate() {
        if !fixed[idx] {
            *zi = seed;
        }
    }

    let t = opts.tension;
    let denom = 20.0 * (1.0 - t) + 4.0 * t;
    for _ in 0..opts.max_iter {
        let mut max_change = 0.0_f64;
        for jp in 0..ny {
            for ip in 0..nx {
                let idx = jp * nx + ip;
                if fixed[idx] {
                    continue;
                }
                let target = update_node(&z, nx, ny, ip, jp, t, denom);
                let old = z[idx];
                let new = old + opts.omega * (target - old);
                z[idx] = new;
                max_change = max_change.max((new - old).abs());
            }
        }
        if max_change < opts.tol {
            break;
        }
    }

    Ok(Surface { nx, ny, z })
}

/// A surface **produced by the petekTools warm-start kernel** — the only thing
/// [`solve_surface_seeded`] accepts as a seed.
///
/// ## Why the newtype survives kernel unification (reassessed 2026-07-04)
/// The newtype's *original* justification was kernel **divergence**: petekTools'
/// kernel lacked the linear-extrapolation natural-dip boundary of the cold
/// [`solve_surface`], so the two solvers converged to *different* interior fixed
/// points (a 12.48 m/ft plane sag) and seeding warm from cold silently violated
/// kernel space. **petekTools f81b6a6 landed the natural-dip boundary**, so the
/// two kernels now share one fixed point (per-node plane parity gated in
/// `tests/adoption_readiness.rs`) — that divergence rationale is retired.
///
/// The newtype is **kept** on its second, still-live justification: it is a
/// **provenance / staleness guard** on the warm chain. A seed must be a genuine
/// same-kernel output — [`KernelSurface::flat`] (a constant depth, a fixed point of
/// the kernel — the safe bootstrap) or a prior [`solve_surface_seeded`] — never an
/// arbitrary [`Surface`]: a cold solve carries caller-chosen [`SolveOpts`]
/// (tension/ω/tol) whose fixed point can differ from the kernel's fixed defaults, and
/// a hand-built or foreign field is simply stale. Keeping the barrier is zero-cost
/// (open/closed, house style) and keeps the warm-start chain sound; it intentionally
/// has no `From<Surface>` (`decision_gridder_kernel_unification`).
#[derive(Debug, Clone)]
pub struct KernelSurface(Surface);

impl KernelSurface {
    /// A flat kernel-space seed at constant depth — the one safe way to bootstrap
    /// the warm-start chain (constant is a fixed point of both kernels).
    ///
    /// # Panics
    /// Panics if `nx < 2` or `ny < 2`.
    #[must_use]
    pub fn flat(nx: usize, ny: usize, z: f64) -> Self {
        Self(Surface::constant(nx, ny, z))
    }

    /// Borrow the underlying [`Surface`] (for layering / offset / node reads).
    #[must_use]
    pub fn surface(&self) -> &Surface {
        &self.0
    }

    /// Areal node count along i (= ni + 1).
    #[must_use]
    pub fn nx(&self) -> usize {
        self.0.nx
    }
    /// Areal node count along j (= nj + 1).
    #[must_use]
    pub fn ny(&self) -> usize {
        self.0.ny
    }
    /// Depth at node `(ip, jp)`.
    #[must_use]
    pub fn z(&self, ip: usize, jp: usize) -> f64 {
        self.0.z(ip, jp)
    }
}

/// Warm-start refine solve: re-grid the same lattice from a prior converged
/// [`KernelSurface`] seed instead of a cold mean-depth seed — the load-bearing
/// per-realization regeneration optimization (SPEC §7a). The warm-start SOR is
/// **delegated to petekTools** (`grid_min_curvature_seeded`, the
/// `ConvergentGridder` kernel): re-seeding from a nearby field converges in far
/// fewer sweeps than a cold solve while reaching the same field to the solver
/// tolerance (`warm == cold` within that kernel).
///
/// Node `(ip, jp)` maps to petekTools lattice node `(ip, jp)` on a unit-spaced
/// [`Lattice`]; the returned [`KernelSurface`] shares the seed's `(nx, ny)`.
///
/// # Kernel constraint (enforced by the type)
/// The seed is a [`KernelSurface`], so it can only be a flat bootstrap or a prior
/// petekTools-kernel output — a cold [`solve_surface`] `Surface` will not type-check
/// here. The kernels now share one fixed point (petekTools f81b6a6 natural-dip), so
/// this is no longer a divergence firewall but a **provenance/staleness guard** on
/// the warm chain (a cold `Surface` may carry different `SolveOpts`, or be foreign/
/// stale). See [`KernelSurface`] for the full reassessment.
///
/// `SolveOpts` is intentionally absent: petekTools owns its solver parameters
/// (ω = 1.5, tol = 1e-6, max_iter = 5000, pure biharmonic).
///
/// # Errors
/// [`StaticError::InvalidInput`] if there are no controls or the seed lattice is
/// smaller than `2x2`; [`StaticError::Grid`] if a control is off the seed lattice.
pub fn solve_surface_seeded(
    seed: &KernelSurface,
    controls: &[Control],
) -> Result<KernelSurface, StaticError> {
    let (nx, ny) = (seed.0.nx, seed.0.ny);
    if nx < 2 || ny < 2 {
        return Err(StaticError::Grid(format!(
            "seed lattice must be at least 2x2, got {nx}x{ny}"
        )));
    }
    if controls.is_empty() {
        return Err(StaticError::InvalidInput(
            "seeded surface needs at least one control point".into(),
        ));
    }
    let mut coords = Vec::with_capacity(controls.len());
    for c in controls {
        if c.ip >= nx || c.jp >= ny {
            return Err(StaticError::Grid(format!(
                "control ({},{}) outside lattice {nx}x{ny}",
                c.ip, c.jp
            )));
        }
        coords.push([c.ip as f64, c.jp as f64, c.z]);
    }

    // Seam conversion: Surface (z[jp*nx+ip]) -> petekTools Array2 shape (nx, ny),
    // indexed [[ip, jp]]; unit-spaced lattice so node (ip,jp) snaps to itself.
    let lattice = Lattice::regular(0.0, 0.0, 1.0, 1.0, nx, ny);
    let mut seed_arr = ndarray::Array2::<f64>::zeros((nx, ny));
    for jp in 0..ny {
        for ip in 0..nx {
            seed_arr[[ip, jp]] = seed.0.z(ip, jp);
        }
    }

    let field = grid_min_curvature_seeded(&coords, &lattice, Some(&seed_arr))
        .map_err(|e| StaticError::Grid(format!("petektools seeded grid failed: {e}")))?;

    let mut z = vec![0.0; nx * ny];
    for jp in 0..ny {
        for ip in 0..nx {
            z[jp * nx + ip] = field[[ip, jp]];
        }
    }
    Ok(KernelSurface(Surface { nx, ny, z }))
}

/// Solve a surface from node controls to the kernel's **true fixed point** —
/// the structure-build solve entry (`structure_fidelity` audit S2). Two measures
/// against the slow-mode stall of a plain seeded solve:
///
/// 1. **Plane detrending.** The kernel's slowest modes are the affine ones: a
///    flat bootstrap seed carries no dip, and under the natural-dip boundary the
///    planar error component decays glacially (metres left after the full sweep
///    budget). A plane is an **exact fixed point** of the tensioned stencil
///    (boundary included), so superposition is exact: fit a least-squares plane
///    to the controls, solve the detrended residual problem (no affine component
///    → converges quickly), and add the plane back.
/// 2. **Fixed-point restarts.** Re-seed the solve from its own output until a
///    whole re-solve moves < 1 mm (a converged re-solve exits in ~1 sweep), so
///    the returned field is the solver's fixed point, not a budget artifact.
///
/// Hard-control contract unchanged: every control node is honoured exactly.
///
/// # Errors
/// [`StaticError`] if there are no controls, the lattice is degenerate, or the
/// underlying seeded solve fails.
pub fn solve_surface_converged(
    nx: usize,
    ny: usize,
    controls: &[Control],
) -> Result<KernelSurface, StaticError> {
    if controls.is_empty() {
        return Err(StaticError::InvalidInput(
            "structural solve needs at least one control point".into(),
        ));
    }
    // Least-squares plane z ≈ a + b·ip + c·jp over the controls (normal
    // equations; falls back to the mean-only plane when degenerate — < 3
    // controls or collinear).
    let n = controls.len() as f64;
    let (mut sx, mut sy, mut sz) = (0.0, 0.0, 0.0);
    let (mut sxx, mut sxy, mut syy, mut sxz, mut syz) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for c in controls {
        let (x, y, z) = (c.ip as f64, c.jp as f64, c.z);
        sx += x;
        sy += y;
        sz += z;
        sxx += x * x;
        sxy += x * y;
        syy += y * y;
        sxz += x * z;
        syz += y * z;
    }
    let (mx, my, mz) = (sx / n, sy / n, sz / n);
    let (cxx, cxy, cyy) = (sxx / n - mx * mx, sxy / n - mx * my, syy / n - my * my);
    let (cxz, cyz) = (sxz / n - mx * mz, syz / n - my * mz);
    let det = cxx * cyy - cxy * cxy;
    let (b, c_coef) = if det.abs() > 1e-9 {
        ((cxz * cyy - cyz * cxy) / det, (cyz * cxx - cxz * cxy) / det)
    } else {
        (0.0, 0.0)
    };
    let a = mz - b * mx - c_coef * my;
    let plane = |ip: usize, jp: usize| a + b * ip as f64 + c_coef * jp as f64;

    // Solve the detrended residual problem (flat bootstrap at its mean — the
    // residuals have no affine component left, so the seed is already close).
    let detrended: Vec<Control> = controls
        .iter()
        .map(|c| Control {
            ip: c.ip,
            jp: c.jp,
            z: c.z - plane(c.ip, c.jp),
        })
        .collect();
    let mean_r = detrended.iter().map(|c| c.z).sum::<f64>() / n;
    let _prof = std::env::var_os("SRS_PROFILE").is_some();
    let _t0 = std::time::Instant::now();
    let mut surf = solve_surface_seeded(&KernelSurface::flat(nx, ny, mean_r), &detrended)?;
    const FIXED_POINT_TOL_M: f64 = 1e-3;
    const MAX_RESTARTS: usize = 16;
    let mut _restarts = 0usize;
    let mut _last_moved = f64::NAN;
    if _prof {
        eprintln!(
            "[SRS_PROFILE] solve_converged START nx={nx} ny={ny} controls={} first_solve_ms={:.1}",
            controls.len(),
            _t0.elapsed().as_secs_f64() * 1e3,
        );
    }
    for _ in 0..MAX_RESTARTS {
        let _tr = std::time::Instant::now();
        let next = solve_surface_seeded(&surf, &detrended)?;
        let mut moved = 0.0_f64;
        for jp in 0..ny {
            for ip in 0..nx {
                moved = moved.max((next.z(ip, jp) - surf.z(ip, jp)).abs());
            }
        }
        surf = next;
        _restarts += 1;
        _last_moved = moved;
        if _prof {
            // Incremental per-restart progress so a timeout-killed run still
            // reveals how far the fixed-point loop got and why it isn't stopping.
            eprintln!(
                "[SRS_PROFILE]   restart {_restarts}/{MAX_RESTARTS} moved={moved:.4}m restart_ms={:.1} cum_ms={:.1}",
                _tr.elapsed().as_secs_f64() * 1e3,
                _t0.elapsed().as_secs_f64() * 1e3,
            );
        }
        if moved < FIXED_POINT_TOL_M {
            break;
        }
    }
    if _prof {
        eprintln!(
            "[SRS_PROFILE] solve_converged DONE controls={} restarts={_restarts} maxed={} last_moved={_last_moved:.4} ms={:.1}",
            controls.len(),
            _restarts == MAX_RESTARTS,
            _t0.elapsed().as_secs_f64() * 1e3,
        );
    }
    // Add the plane back (exact: a plane is a fixed point of the stencil, and
    // the problem is linear, so solution(controls) = solution(residuals) + plane).
    let mut z = vec![0.0; nx * ny];
    for jp in 0..ny {
        for ip in 0..nx {
            z[jp * nx + ip] = surf.z(ip, jp) + plane(ip, jp);
        }
    }
    Ok(KernelSurface(Surface { nx, ny, z }))
}

/// The relaxation target for a free node: the blended biharmonic/harmonic
/// update. Out-of-lattice stencil nodes are synthesized by linear extrapolation
/// ([`z_at`]) — the natural minimum-curvature boundary condition, which keeps
/// linear (planar regional-dip) fields as exact fixed points everywhere.
fn update_node(z: &[f64], nx: usize, ny: usize, ip: usize, jp: usize, t: f64, denom: f64) -> f64 {
    let i = ip as isize;
    let j = jp as isize;
    let at = |di: isize, dj: isize| z_at(z, nx, ny, i + di, j + dj);

    let e1 = at(-1, 0) + at(1, 0) + at(0, -1) + at(0, 1);
    let d = at(-1, -1) + at(1, -1) + at(-1, 1) + at(1, 1);
    let w2 = at(-2, 0) + at(2, 0) + at(0, -2) + at(0, 2);
    ((1.0 - t) * (8.0 * e1 - 2.0 * d - w2) + t * e1) / denom
}

/// Read `z` at a possibly out-of-lattice node, extending the field linearly
/// beyond each edge (so a planar field is reproduced exactly). Bounded recursion
/// — the stencil reaches at most two nodes past an edge.
fn z_at(z: &[f64], nx: usize, ny: usize, i: isize, j: isize) -> f64 {
    let nxi = nx as isize;
    let nyi = ny as isize;
    // Step inward by one and extend linearly (z[t] = 2*z[t+1] - z[t+2]).
    if i < 0 {
        return 2.0 * z_at(z, nx, ny, i + 1, j) - z_at(z, nx, ny, i + 2, j);
    }
    if i >= nxi {
        return 2.0 * z_at(z, nx, ny, i - 1, j) - z_at(z, nx, ny, i - 2, j);
    }
    if j < 0 {
        return 2.0 * z_at(z, nx, ny, i, j + 1) - z_at(z, nx, ny, i, j + 2);
    }
    if j >= nyi {
        return 2.0 * z_at(z, nx, ny, i, j - 1) - z_at(z, nx, ny, i, j - 2);
    }
    z[(j * nxi + i) as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn solve(nx: usize, ny: usize, c: &[Control]) -> Surface {
        solve_surface(nx, ny, c, SolveOpts::default()).unwrap()
    }

    #[test]
    fn constant_controls_give_constant_surface() {
        let c = [
            Control {
                ip: 0,
                jp: 0,
                z: 5000.0,
            },
            Control {
                ip: 9,
                jp: 9,
                z: 5000.0,
            },
        ];
        let s = solve(10, 10, &c);
        for jp in 0..10 {
            for ip in 0..10 {
                assert!(
                    (s.z(ip, jp) - 5000.0).abs() < 1e-4,
                    "node ({ip},{jp})={}",
                    s.z(ip, jp)
                );
            }
        }
    }

    #[test]
    fn repair_min_thickness_floors_thin_columns_and_leaves_the_rest() {
        let top = Surface {
            nx: 2,
            ny: 2,
            z: vec![0.0; 4],
        };
        // Nodes: a crossing (-5), a thin-but-positive (0.5 < 2), and two OK (3, 10).
        let base = Surface {
            nx: 2,
            ny: 2,
            z: vec![-5.0, 0.5, 3.0, 10.0],
        };
        let (rep, columns, worst_m) = base.repair_min_thickness(&top, 2.0).unwrap();
        assert_eq!(columns, 2, "two nodes sit below the 2 m floor");
        assert!(
            (worst_m - -5.0).abs() < 1e-12,
            "worst = most-negative original sep: {worst_m}"
        );
        assert!(
            (rep.z(0, 0) - 2.0).abs() < 1e-12,
            "crossing floored to top+2"
        );
        assert!((rep.z(1, 0) - 2.0).abs() < 1e-12, "thin floored to top+2");
        assert_eq!(rep.z(0, 1), 3.0, "OK node bit-unchanged");
        assert_eq!(rep.z(1, 1), 10.0, "OK node bit-unchanged");
        // A base already at/below the floor repairs nothing.
        let clean = Surface {
            nx: 2,
            ny: 2,
            z: vec![5.0; 4],
        };
        let (_, n, w) = clean.repair_min_thickness(&top, 2.0).unwrap();
        assert_eq!((n, w), (0, 0.0), "nothing repaired -> zero count/worst");
        // Guards: a mismatched lattice and a non-finite / negative floor error.
        let wrong = Surface {
            nx: 3,
            ny: 1,
            z: vec![0.0; 3],
        };
        assert!(base.repair_min_thickness(&wrong, 2.0).is_err());
        assert!(base.repair_min_thickness(&top, -1.0).is_err());
        assert!(base.repair_min_thickness(&top, f64::NAN).is_err());
    }

    #[test]
    fn repair_from_below_moves_the_upper_up_preserving_the_mapped_lower() {
        // Repair-precedence twin: the UPPER (a derived surface) yields to the LOWER
        // (a mapped one). `lower` is the authoritative mapped surface at depth 10.
        let lower = Surface {
            nx: 2,
            ny: 2,
            z: vec![10.0; 4],
        };
        // Upper nodes: a crossing (15 > 10), a thin-but-above (9.5, gap 0.5 < 2), and
        // two OK (7, 3).
        let upper = Surface {
            nx: 2,
            ny: 2,
            z: vec![15.0, 9.5, 7.0, 3.0],
        };
        let (rep, columns, worst_m) = upper.repair_min_thickness_from_below(&lower, 2.0).unwrap();
        assert_eq!(
            columns, 2,
            "two nodes sit within the 2 m floor of the lower"
        );
        assert!(
            (worst_m - -5.0).abs() < 1e-12,
            "worst = most-negative original lower-upper: {worst_m}"
        );
        assert!(
            (rep.z(0, 0) - 8.0).abs() < 1e-12,
            "crossing lifted to lower-2"
        );
        assert!((rep.z(1, 0) - 8.0).abs() < 1e-12, "thin lifted to lower-2");
        assert_eq!(rep.z(0, 1), 7.0, "OK node bit-unchanged");
        assert_eq!(rep.z(1, 1), 3.0, "OK node bit-unchanged");
        // The mapped lower is never touched by this operation.
        assert!(lower.z.iter().all(|&z| z == 10.0));
    }

    #[test]
    fn guard_above_clamps_or_errors_a_derived_upper_crossing_a_mapped_lower() {
        let lower = Surface {
            nx: 2,
            ny: 2,
            z: vec![10.0; 4],
        };
        let upper = Surface {
            nx: 2,
            ny: 2,
            z: vec![15.0, 5.0, 12.0, 3.0], // two nodes cross below the lower
        };
        // clamp: the crossing nodes are pulled UP to the mapped lower (zero gap there).
        let clamped = upper.guard_above(&lower, true).unwrap();
        assert_eq!(clamped.z(0, 0), 10.0, "crossing clamped to lower");
        assert_eq!(clamped.z(0, 1), 10.0, "crossing clamped to lower");
        assert_eq!(clamped.z(1, 0), 5.0, "non-crossing untouched");
        assert_eq!(clamped.z(1, 1), 3.0, "non-crossing untouched");
        // no-clamp: a crossing errors.
        assert!(matches!(
            upper.guard_above(&lower, false),
            Err(StaticError::CrossedSurfaces { nodes: 2, .. })
        ));
        // No crossing -> returned unchanged.
        let ok = Surface {
            nx: 2,
            ny: 2,
            z: vec![1.0, 2.0, 3.0, 4.0],
        };
        assert!(ok.guard_above(&lower, false).is_ok());
    }

    #[test]
    fn reproduces_a_plane() {
        // Plane z = 5000 + 2*ip + 3*jp has zero curvature -> min-curvature exact.
        let plane = |ip: usize, jp: usize| 5000.0 + 2.0 * ip as f64 + 3.0 * jp as f64;
        let mut c = Vec::new();
        // Constrain the four corners + center (enough to pin the plane).
        for &(ip, jp) in &[(0, 0), (12, 0), (0, 12), (12, 12), (6, 6)] {
            c.push(Control {
                ip,
                jp,
                z: plane(ip, jp),
            });
        }
        let s = solve_surface(
            13,
            13,
            &c,
            SolveOpts {
                tol: 1e-9,
                max_iter: 60_000,
                ..SolveOpts::default()
            },
        )
        .unwrap();
        for jp in 0..13 {
            for ip in 0..13 {
                let got = s.z(ip, jp);
                assert!(
                    (got - plane(ip, jp)).abs() < 0.05,
                    "({ip},{jp}) {got} != {}",
                    plane(ip, jp)
                );
            }
        }
    }

    #[test]
    fn honors_control_points_exactly() {
        let c = [
            Control {
                ip: 0,
                jp: 0,
                z: 5000.0,
            },
            Control {
                ip: 5,
                jp: 5,
                z: 5200.0,
            },
            Control {
                ip: 9,
                jp: 9,
                z: 5050.0,
            },
        ];
        let s = solve(10, 10, &c);
        assert!((s.z(0, 0) - 5000.0).abs() < 1e-9);
        assert!((s.z(5, 5) - 5200.0).abs() < 1e-9);
        assert!((s.z(9, 9) - 5050.0).abs() < 1e-9);
        // The bump pulls nearby nodes up between the lows.
        assert!(s.z(5, 5) > s.z(2, 2));
    }

    #[test]
    fn rejects_bad_input() {
        assert!(solve_surface(0, 5, &[], SolveOpts::default()).is_err());
        assert!(solve_surface(5, 5, &[], SolveOpts::default()).is_err());
        let c = [Control {
            ip: 9,
            jp: 0,
            z: 1.0,
        }];
        assert!(solve_surface(5, 5, &c, SolveOpts::default()).is_err()); // off-lattice
    }

    // --- seeded refine path (petekTools ConvergentGridder delegation) ---

    #[test]
    fn seeded_honors_controls_as_hard_constraints() {
        // Warm-solve on a flat seed; the added control must be held exactly.
        let seed = KernelSurface::flat(12, 10, 5000.0);
        let c = [Control {
            ip: 6,
            jp: 7,
            z: 4900.0,
        }];
        let s = solve_surface_seeded(&seed, &c).unwrap();
        assert!(
            (s.z(6, 7) - 4900.0).abs() < 1e-6,
            "control node = {}",
            s.z(6, 7)
        );
        assert_eq!((s.nx(), s.ny()), (12, 10));
    }

    // Re-verified 2026-07-04 after petekTools f81b6a6 (natural-dip boundary) landed
    // (un-ignored; was blocked on the concurrent kernel WIP). The kernel now runs to
    // the family cold cap (20k sweeps, abs TOL 1e-6), so seeding it from a CONVERGED
    // field of the same kernel reproduces that field to solver tolerance — a warm
    // re-solve of the fixed point stops in ~1 sweep. Tolerance tightened to 1e-6.
    //
    // Config note: a **well-determined** control set (edge/corner nodes pinned, as
    // every defined-node Top surface produces in the real regeneration seam) so a
    // single seeded solve actually reaches the fixed point. With only a few interior
    // controls the smooth (level) mode relaxes slowly under the near-Neumann
    // natural-dip boundary, so a lone flat-bootstrap solve is not yet converged and
    // its warm successor keeps drifting toward the fixed point — see the todos
    // follow-up (the build-vs-realize gap on sparse configs rides on R2 kernel
    // unification, not on natural-dip).
    #[test]
    fn seeded_warm_equals_cold_of_same_kernel() {
        // The continuity guarantee the regeneration seam relies on: seeding the
        // petekTools kernel from a converged field of the SAME kernel reproduces it
        // (warm == cold to the solver tolerance). Cold reference = seed a flat field,
        // solve once; warm = seed from that result, solve again.
        let flat = KernelSurface::flat(12, 10, 5000.0);
        // Corners + an interior high — a well-determined framework (the level mode is
        // pinned), so one seeded solve converges.
        let controls: Vec<Control> = [
            (0, 0, 5000.0),
            (11, 0, 5010.0),
            (0, 9, 4990.0),
            (11, 9, 5005.0),
            (5, 5, 5030.0),
        ]
        .iter()
        .map(|&(ip, jp, z)| Control { ip, jp, z })
        .collect();
        let cold = solve_surface_seeded(&flat, &controls).unwrap();
        let warm = solve_surface_seeded(&cold, &controls).unwrap();
        for jp in 0..10 {
            for ip in 0..12 {
                assert!(
                    (warm.z(ip, jp) - cold.z(ip, jp)).abs() < 1e-6,
                    "({ip},{jp}) warm {} vs cold {}",
                    warm.z(ip, jp),
                    cold.z(ip, jp)
                );
            }
        }
    }

    #[test]
    fn seeded_rejects_bad_input() {
        let seed = KernelSurface::flat(5, 5, 0.0);
        assert!(solve_surface_seeded(&seed, &[]).is_err()); // no controls
        let off = [Control {
            ip: 9,
            jp: 0,
            z: 1.0,
        }];
        assert!(solve_surface_seeded(&seed, &off).is_err()); // off-lattice
        let tiny = KernelSurface::flat(2, 2, 0.0);
        let c = [Control {
            ip: 0,
            jp: 0,
            z: 1.0,
        }];
        assert!(solve_surface_seeded(&tiny, &c).is_ok()); // 2x2 is the minimum
    }

    // --- R5 degenerate-input property tests for the order-repair kernel
    // (`repair_min_thickness`): the post-repair base must sit at least
    // `min_thickness` below the top at EVERY node, whatever the input crossing
    // (inverted / zero-thickness / random), and the repair is idempotent. A NaN
    // `min_thickness` is a typed error, never a silent bad surface. ---

    proptest::proptest! {
        #![proptest_config(proptest::prelude::ProptestConfig::with_cases(150))]

        #[test]
        fn prop_repair_min_thickness_enforces_floor_and_is_idempotent(
            // Per-node base offsets relative to a flat top at 5000 — including
            // negatives (base above top = a crossing) and zeros.
            offsets in proptest::collection::vec(-60.0f64..80.0, 9),
            min_t in 0.0f64..30.0,
        ) {
            let top = Surface::constant(3, 3, 5000.0);
            let base = top.offset_by_field(&offsets).unwrap();
            let (repaired, count, worst) = base.repair_min_thickness(&top, min_t).unwrap();
            // (a) the floor holds everywhere post-repair.
            for jp in 0..3 {
                for ip in 0..3 {
                    let sep = repaired.z(ip, jp) - top.z(ip, jp);
                    proptest::prop_assert!(
                        sep >= min_t - 1e-6,
                        "node ({ip},{jp}) thickness {sep} < min {min_t}"
                    );
                }
            }
            // (b) reported count matches the nodes that were below the floor.
            let expected = offsets.iter().filter(|&&o| o < min_t).count();
            proptest::prop_assert_eq!(count, expected);
            proptest::prop_assert!(worst <= 1e-9);
            // (c) idempotent: repairing the repaired base moves nothing.
            let (again, count2, _) = repaired.repair_min_thickness(&top, min_t).unwrap();
            proptest::prop_assert_eq!(count2, 0);
            for jp in 0..3 {
                for ip in 0..3 {
                    proptest::prop_assert!((again.z(ip, jp) - repaired.z(ip, jp)).abs() < 1e-9);
                }
            }
        }
    }

    #[test]
    fn repair_min_thickness_rejects_nonfinite_min() {
        let top = Surface::constant(3, 3, 5000.0);
        let base = top.offset_by(20.0);
        assert!(base.repair_min_thickness(&top, f64::NAN).is_err());
        assert!(base.repair_min_thickness(&top, -1.0).is_err());
    }
}
