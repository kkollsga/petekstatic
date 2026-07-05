//! `srs-gridder` — the convergent gridder.
//!
//! [`solve_surface`] grids a minimum-curvature surface from depth control points
//! (`convergent_gridder_spec`, MEDIUM); [`solve_surface_seeded`] is the
//! warm-start refine path (SPEC §7a) — it delegates the seeded SOR to petekTools'
//! `ConvergentGridder` kernel. [`layer_grid`] turns top/base surfaces into a
//! conformable k-layered corner-point grid (`layer_interpolation_spec`, high).
//! The unfaulted, no-horizon case degenerates to the box grid
//! (`srs_grid::build_box`). Faults / NNCs are deferred — see
//! `question_gridder_spec`.

mod layering;
mod surface;

pub use layering::{
    layer_grid, layer_grid_stack, layer_grid_stack_into, Conformity, LayerScratch, LayeredGrid,
    StackLayering, StackedLayeredGrid, StackedZone, StreamingLayering, ZoneLayerSpec, MAX_NK,
};
pub use surface::{
    solve_surface, solve_surface_converged, solve_surface_seeded, Control, ExtrapolationPolicy,
    KernelSurface, SolveOpts, Surface,
};
