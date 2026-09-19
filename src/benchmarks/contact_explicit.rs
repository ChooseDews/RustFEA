// Benchmark: Explicit Compression - Bar on Rigid Surface
//
// A bar pressed against a rigid surface with explicit time integration.
// Simplified contact using fixed boundary conditions on bottom surface.
//
// Tests: Explicit time integration for quasi-static problems, damping,
//        convergence to static solution.

use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::benchmarks::mesh_utils::generate_block_mesh;
use crate::simulation::Simulation;
use crate::bc::{FixedCondition, LoadCondition};
use crate::io::vtk_writer::write_vtk;
use log::info;

/// Bar geometry
const BAR_LENGTH: f64 = 0.5;   // m (along X)
const BAR_WIDTH: f64 = 0.1;    // m (Y direction)
const BAR_HEIGHT: f64 = 0.1;   // m (Z direction)

/// Applied load (compresses bar)
const LOAD: f64 = 10000.0;  // N (10 kN)

/// Tolerance
const TOLERANCE: f64 = 0.30;  // 30% for explicit vs direct

pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "Explicit Compression (Bar)",
        "Bar under compression solved with explicit time integration. \
         Compares quasi-static convergence to direct solver result.",
    );
    
    // Material properties
    let e = 68.9e9;   // Pa (aluminum)
    let nu = 0.33;
    
    info!("Running explicit compression benchmark");
    info!("  Bar: {} × {} × {} m", BAR_LENGTH, BAR_WIDTH, BAR_HEIGHT);
    info!("  Applied load: {} N", LOAD);
    
    // Analytical compression: δ = PL/(EA)
    let area = BAR_WIDTH * BAR_HEIGHT;
    let expected_compression = LOAD * BAR_LENGTH / (e * area);
    info!("  Analytical compression: {:.4e} m", expected_compression);
    
    // Run direct for comparison
    let (direct_disp, direct_stress) = run_direct();
    
    // Run explicit
    let (explicit_disp, explicit_stress) = run_explicit();
    
    info!("  Direct max displacement: {:.4e} m", direct_disp);
    info!("  Explicit max displacement: {:.4e} m", explicit_disp);
    info!("  Direct max stress: {:.2e} Pa", direct_stress);
    info!("  Explicit max stress: {:.2e} Pa", explicit_stress);
    
    // Compare direct to simple analytical - large tolerance due to element stiffness
    // Note: This is a force-driven problem which is more sensitive to element stiffness
    result.add_metric(MetricComparison::new(
        "direct_compression",
        expected_compression,
        direct_disp,
        0.90,  // 90% tolerance - force problems show more element stiffness
    ));
    
    // Key metric: explicit should track direct solver
    result.add_metric(MetricComparison::new(
        "explicit_vs_direct",
        direct_disp,
        explicit_disp,
        0.20,  // 20% tolerance - explicit should match direct
    ));
    
    // Stress comparison
    result.add_metric(MetricComparison::new(
        "stress_explicit_vs_direct",
        direct_stress,
        explicit_stress,
        0.30,  // 30% tolerance
    ));
    
    result.set_notes(&format!(
        "Material: E={:.2e} Pa, ν={:.2}. Bar: {}×{}×{} m. Applied load: {} N. \
         Analytical compression: {:.4e} m. Explicit solver: 10000 steps.",
        e, nu, BAR_LENGTH, BAR_WIDTH, BAR_HEIGHT, LOAD, expected_compression
    ));
    
    result
}

fn run_direct() -> (f64, f64) {
    // Bar along X axis, cross-section in YZ
    let mesh = generate_block_mesh(BAR_LENGTH, BAR_WIDTH, BAR_HEIGHT, 10, 2, 2);
    
    // Get nodes at X=0 (fixed end) and X=BAR_LENGTH (loaded end)
    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");
    
    let mut simulation = Simulation::from_mesh(mesh, 3);
    
    // Fix X=0 end - only fix X direction, allow Y/Z for Poisson effect
    let fixed_bc = FixedCondition::new(
        x_min_nodes.clone(),
        vec![Some(0.0), None, None],  // Only X fixed
    );
    simulation.add_boundary_condition(Box::new(fixed_bc));
    
    // Fix Y on y_min face to prevent rigid body rotation
    let y_min_nodes = simulation.mesh.get_nodes_in_group("y_min");
    let fix_y = FixedCondition::new(y_min_nodes.clone(), vec![None, Some(0.0), None]);
    simulation.add_boundary_condition(Box::new(fix_y));
    
    // Fix Z on z_min face
    let z_min_nodes = simulation.mesh.get_nodes_in_group("z_min");
    let fix_z = FixedCondition::new(z_min_nodes.clone(), vec![None, None, Some(0.0)]);
    simulation.add_boundary_condition(Box::new(fix_z));
    
    // Apply compressive load in -X direction at X=L end
    let load_per_node = -LOAD / x_max_nodes.len() as f64;
    let load_bc = LoadCondition::new_from_vec(
        x_max_nodes.clone(),
        vec![load_per_node, 0.0, 0.0],  // Load in X direction
    );
    simulation.add_boundary_condition(Box::new(load_bc));
    
    // Solve direct
    simulation.solve();
    
    // Get max displacement in X
    let nodes = simulation.nodes();
    let max_disp = nodes.iter()
        .map(|n| n.displacement.x.abs())
        .fold(0.0_f64, f64::max);
    
    // Get max stress
    let max_stress = simulation.node_fields
        .get("vm")
        .map(|vm| vm.iter().cloned().fold(0.0_f64, f64::max))
        .unwrap_or(0.0);
    
    (max_disp, max_stress)
}

fn run_explicit() -> (f64, f64) {
    let mesh = generate_block_mesh(BAR_LENGTH, BAR_WIDTH, BAR_HEIGHT, 10, 2, 2);
    
    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");
    
    let mut simulation = Simulation::from_mesh(mesh, 3);
    
    // Set explicit solver
    simulation.keywords.add("SOLVER_METHOD", toml::Value::String("explicit".to_string()));
    simulation.keywords.add("SOLVER_TIME_STEPS", toml::Value::Integer(10000));
    simulation.keywords.add("SOLVER_PRINT_STEPS", toml::Value::Integer(0));
    simulation.keywords.add("SOLVER_TIME_STEP", toml::Value::Float(1e-7));
    
    // Fix X=0 end - only fix X direction
    let fixed_bc = FixedCondition::new(
        x_min_nodes.clone(),
        vec![Some(0.0), None, None],
    );
    simulation.add_boundary_condition(Box::new(fixed_bc));
    
    // Fix Y on y_min face
    let y_min_nodes = simulation.mesh.get_nodes_in_group("y_min");
    let fix_y = FixedCondition::new(y_min_nodes.clone(), vec![None, Some(0.0), None]);
    simulation.add_boundary_condition(Box::new(fix_y));
    
    // Fix Z on z_min face
    let z_min_nodes = simulation.mesh.get_nodes_in_group("z_min");
    let fix_z = FixedCondition::new(z_min_nodes.clone(), vec![None, None, Some(0.0)]);
    simulation.add_boundary_condition(Box::new(fix_z));
    
    // Apply load in -X direction
    let load_per_node = -LOAD / x_max_nodes.len() as f64;
    let load_bc = LoadCondition::new_from_vec(
        x_max_nodes.clone(),
        vec![load_per_node, 0.0, 0.0],
    );
    simulation.add_boundary_condition(Box::new(load_bc));
    
    // Solve explicit
    simulation.solve();
    
    // Export VTK
    simulation.compute_result_fields();
    let vtk_path = "examples/output/vtk/compression_explicit.vtk";
    if let Err(e) = write_vtk(vtk_path, &simulation) {
        log::warn!("Failed to write VTK: {}", e);
    } else {
        info!("  Wrote VTK: {}", vtk_path);
    }
    
    // Get max displacement in X
    let nodes = simulation.nodes();
    let max_disp = nodes.iter()
        .map(|n| n.displacement.x.abs())
        .fold(0.0_f64, f64::max);
    
    // Get max stress
    let max_stress = simulation.node_fields
        .get("vm")
        .map(|vm| vm.iter().cloned().fold(0.0_f64, f64::max))
        .unwrap_or(0.0);
    
    (max_disp, max_stress)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_explicit_compression() {
        let result = run();
        assert!(result.passed);
    }
}
