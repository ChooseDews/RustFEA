// Benchmark: C3D20 Cantilever Beam
//
// Tests the 20-node quadratic hexahedral element (C3D20) under bending load.
// A cantilever beam fixed at one end with a tip load at the other.
//
// Analytical solution for tip deflection:
//   δ = P * L³ / (3 * E * I)
// where:
//   P = applied load
//   L = beam length
//   E = Young's modulus
//   I = second moment of area = b * h³ / 12
//
// C3D20 elements should show much better bending behavior than C3D8 due to
// quadratic shape functions which reduce shear locking.

use crate::bc::{FixedCondition, LoadCondition};
use crate::benchmarks::mesh_utils::generate_block_mesh_c3d20;
use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::simulation::Simulation;
use log::info;
use nalgebra::DVector;

/// Material properties (Aluminum 6061-T6)
const E: f64 = 68.9e9; // Young's modulus (Pa)
const NU: f64 = 0.33; // Poisson's ratio

/// Geometry - beam along X axis, cross-section in YZ
const LENGTH: f64 = 1.0; // Beam length in x-direction (m)
const WIDTH: f64 = 0.1; // Width in y-direction (m)
const HEIGHT: f64 = 0.1; // Height in z-direction (m)

/// Loading
const TIP_LOAD: f64 = 1000.0; // Tip load in -Z direction (N)

/// Tolerance - C3D20 should be much more accurate than C3D8 for bending
/// But still allow for some mesh sensitivity
const TOLERANCE: f64 = 0.25; // 25% tolerance

/// Run the C3D20 cantilever beam benchmark
pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "C3D20 Cantilever Beam",
        "20-node quadratic hexahedral elements under bending. \
         Tests tip deflection against analytical Euler-Bernoulli beam theory.",
    );

    // Second moment of area
    let i_zz = WIDTH * HEIGHT.powi(3) / 12.0;

    // Analytical tip deflection (downward, so negative in our coordinate system)
    let delta_analytical = TIP_LOAD * LENGTH.powi(3) / (3.0 * E * i_zz);

    info!("C3D20 Cantilever Beam Benchmark");
    info!("  Analytical tip deflection: {:.6e} m", delta_analytical);
    info!("  I_zz = {:.6e} m⁴", i_zz);

    // Run with multiple mesh refinements
    // Using roughly 20 elements total as requested
    let refinements = vec![
        (5, 2, 2, "20_elem"),  // 5*2*2 = 20 elements
        (10, 2, 2, "40_elem"), // 10*2*2 = 40 elements (refinement check)
    ];

    for (nx, ny, nz, mesh_name) in refinements {
        info!(
            "Running {} mesh ({} x {} x {} = {} elements)",
            mesh_name,
            nx,
            ny,
            nz,
            nx * ny * nz
        );

        let mesh_result = run_single_mesh(nx, ny, nz, delta_analytical);

        // Add metrics with mesh name prefix
        for metric in mesh_result.metrics {
            let mut named_metric = metric.clone();
            named_metric.name = format!("{}_{}", mesh_name, metric.name);
            result.add_metric(named_metric);
        }
    }

    result.set_notes(&format!(
        "Material: E={:.2e} Pa, ν={:.2}. Geometry: L={}m, {}×{}m cross-section. \
         I={:.4e} m⁴. Tip load P={}N in -Z direction. Element: C3D20",
        E, NU, LENGTH, WIDTH, HEIGHT, i_zz, TIP_LOAD
    ));

    result
}

fn run_single_mesh(n_x: usize, n_y: usize, n_z: usize, delta_analytical: f64) -> BenchmarkResult {
    let mut result = BenchmarkResult::new("single_mesh", "");

    // Generate C3D20 mesh
    let mut mesh = generate_block_mesh_c3d20(LENGTH, WIDTH, HEIGHT, n_x, n_y, n_z);

    // Shift to center cross-section at Y=0, Z=0 (beam centered on X axis)
    for (_, node) in mesh.nodes.iter_mut() {
        node.coordinates[1] -= WIDTH / 2.0;
        node.coordinates[2] -= HEIGHT / 2.0;
    }

    // Get node groups before moving mesh into simulation
    let x_min_nodes = mesh.get_nodes_in_group("x_min"); // Fixed end
    let x_max_nodes = mesh.get_nodes_in_group("x_max"); // Load end

    info!(
        "  Mesh: {} nodes, {} elements",
        mesh.nodes.len(),
        mesh.elements.len()
    );
    info!(
        "  Fixed end nodes: {}, Load end nodes: {}",
        x_min_nodes.len(),
        x_max_nodes.len()
    );

    // Create simulation
    let mut sim = Simulation::from_mesh(mesh, 3);

    // Apply boundary conditions:
    // 1. Fix the root (x=0 face) completely - cantilever support
    let fixed_bc = FixedCondition::new(x_min_nodes.clone(), vec![Some(0.0), Some(0.0), Some(0.0)]);
    sim.add_boundary_condition(Box::new(fixed_bc));

    // 2. Apply tip load in -Z direction on x=L face
    // LoadCondition distributes the total force over all nodes in the group
    let load_bc = LoadCondition::new(
        x_max_nodes.clone(),
        DVector::from_vec(vec![0.0, 0.0, -TIP_LOAD]),
    );
    sim.add_boundary_condition(Box::new(load_bc));

    // Solve
    sim.solve();

    // Get nodes for verification
    let nodes = sim.nodes();

    // Compute average tip deflection (should be negative, i.e., downward)
    let tip_deflections: Vec<f64> = x_max_nodes
        .iter()
        .filter_map(|&id| nodes.get(id))
        .map(|node| node.displacement.z)
        .collect();

    let avg_tip_deflection = if !tip_deflections.is_empty() {
        tip_deflections.iter().sum::<f64>() / tip_deflections.len() as f64
    } else {
        0.0
    };

    // Deflection should be negative (downward); analytical is positive magnitude
    let computed_delta = -avg_tip_deflection;

    info!(
        "  Computed tip deflection: {:.6e} m (analytical: {:.6e} m)",
        computed_delta, delta_analytical
    );

    // Compute error
    let relative_error = (computed_delta - delta_analytical).abs() / delta_analytical;
    info!("  Relative error: {:.2}%", relative_error * 100.0);

    // Add metrics
    result.add_metric(MetricComparison::new(
        "tip_deflection",
        delta_analytical,
        computed_delta,
        TOLERANCE,
    ));

    // Stiffness ratio: analytical_deflection / computed_deflection
    // > 1 means FEA is too stiff (less deflection than expected)
    // < 1 means FEA is too soft
    let stiffness_ratio = if computed_delta.abs() > 1e-15 {
        delta_analytical / computed_delta
    } else {
        f64::INFINITY
    };

    info!(
        "  Stiffness ratio: {:.3} (1.0 = perfect, >1 = too stiff)",
        stiffness_ratio
    );

    result.add_metric(MetricComparison::new(
        "stiffness_ratio",
        1.0,
        stiffness_ratio,
        0.5, // Allow up to 50% stiffness ratio deviation
    ));

    // Report mesh info
    let num_nodes = nodes.len();
    let num_elements = sim.elements().len();
    result.set_notes(&format!(
        "Mesh: {} C3D20 elements, {} nodes",
        num_elements, num_nodes
    ));

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c3d20_cantilever_beam() {
        // Initialize logging for test output
        let _ = env_logger::builder()
            .filter_level(log::LevelFilter::Info)
            .is_test(true)
            .try_init();

        let result = run();

        println!("\nC3D20 Cantilever Beam Results:");
        println!("==============================");
        for metric in &result.metrics {
            let status = if metric.passed() { "PASS" } else { "FAIL" };
            println!(
                "  {} [{}]: analytical={:.4e}, computed={:.4e}, error={:.2}%",
                metric.name,
                status,
                metric.analytical,
                metric.computed,
                metric.relative_error.abs() * 100.0
            );
        }

        assert!(
            result.passed,
            "C3D20 cantilever beam benchmark failed. Metrics: {:?}",
            result.metrics
        );
    }

    #[test]
    fn test_c3d20_cantilever_single_element() {
        // Test with minimal mesh to verify basic functionality
        let _ = env_logger::builder()
            .filter_level(log::LevelFilter::Info)
            .is_test(true)
            .try_init();

        // Single element along length - this is a stress test
        let mut mesh = generate_block_mesh_c3d20(LENGTH, WIDTH, HEIGHT, 1, 1, 1);

        // Center the beam
        for (_, node) in mesh.nodes.iter_mut() {
            node.coordinates[1] -= WIDTH / 2.0;
            node.coordinates[2] -= HEIGHT / 2.0;
        }

        let x_min_nodes = mesh.get_nodes_in_group("x_min");
        let x_max_nodes = mesh.get_nodes_in_group("x_max");

        println!(
            "Single element mesh: {} nodes, {} elements",
            mesh.nodes.len(),
            mesh.elements.len()
        );
        println!(
            "Fixed nodes: {}, Load nodes: {}",
            x_min_nodes.len(),
            x_max_nodes.len()
        );

        let mut sim = Simulation::from_mesh(mesh, 3);

        let fixed_bc =
            FixedCondition::new(x_min_nodes.clone(), vec![Some(0.0), Some(0.0), Some(0.0)]);
        sim.add_boundary_condition(Box::new(fixed_bc));

        let load_bc = LoadCondition::new(
            x_max_nodes.clone(),
            DVector::from_vec(vec![0.0, 0.0, -TIP_LOAD]),
        );
        sim.add_boundary_condition(Box::new(load_bc));

        // This should at least not panic
        sim.solve();

        let nodes = sim.nodes();
        let tip_z: Vec<f64> = x_max_nodes
            .iter()
            .filter_map(|&id| nodes.get(id))
            .map(|n| n.displacement.z)
            .collect();

        let avg_tip = tip_z.iter().sum::<f64>() / tip_z.len() as f64;
        println!("Single element tip deflection: {:.6e} m", avg_tip);

        // Just verify we got a reasonable negative deflection
        assert!(
            avg_tip < 0.0,
            "Expected downward (negative) deflection, got {}",
            avg_tip
        );
    }
}
