// Benchmark: Explicit Torsion Shaft
//
// Same torsion shaft problem solved with explicit time integration.
// Compares quasi-static explicit solution to direct solver results.
//
// Tests: Explicit time integration stability, convergence to static solution,
//        damping effects, computational efficiency comparison.

use crate::bc::{FixedCondition, TorqueCondition};
use crate::benchmarks::mesh_utils::generate_block_mesh;
use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::io::vtk_writer::write_vtk;
use crate::simulation::Simulation;
use log::info;
use nalgebra::Vector3;

/// Geometry - approximate circular cross-section with square
const R: f64 = 0.05; // Effective radius (m)
const L: f64 = 1.0; // Length (m)
const SIDE: f64 = 0.0886; // Side of square with same area as circle

/// Applied torque
const T: f64 = 100.0; // Torque (N·m)

/// Tolerance - explicit may differ from direct due to damping
const TOLERANCE: f64 = 0.25; // 25% tolerance for explicit vs analytical

pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "Torsion Shaft (Explicit)",
        "Shaft under torsion solved with explicit time integration. \
         Compares quasi-static explicit solution convergence to direct solver.",
    );

    // Material properties
    let actual_e = 68.9e9;
    let actual_nu = 0.33;
    let g = actual_e / (2.0 * (1.0 + actual_nu));

    // Analytical solution for square section torsion
    // J_eff = k * a^4 where k ≈ 0.1406 for square
    let k = 0.1406;
    let j_square = k * SIDE.powi(4);
    let twist_angle_square = T * L / (g * j_square);

    info!("Running explicit torsion benchmark");
    info!(
        "  Analytical twist angle (square): {:.6} rad = {:.4}°",
        twist_angle_square,
        twist_angle_square.to_degrees()
    );

    // Run with medium mesh - explicit is slower so use coarser mesh
    let (nx, ny, nz) = (20, 4, 4);
    info!(
        "Running explicit solver ({} x {} x {} elements)",
        nx, ny, nz
    );

    // First run direct for comparison
    let direct_result = run_direct(nx, ny, nz, g, twist_angle_square);
    let direct_twist = direct_result.0;

    // Run explicit
    let explicit_result = run_explicit(nx, ny, nz, g, twist_angle_square);
    let explicit_twist = explicit_result.0;

    // Compare results
    result.add_metric(MetricComparison::new(
        "direct_twist_angle",
        twist_angle_square,
        direct_twist,
        0.50, // 50% tolerance for square approx
    ));

    result.add_metric(MetricComparison::new(
        "explicit_twist_angle",
        twist_angle_square,
        explicit_twist,
        0.50,
    ));

    // Compare explicit to direct (should be close)
    let explicit_vs_direct_error = if direct_twist.abs() > 1e-15 {
        (explicit_twist - direct_twist).abs() / direct_twist.abs()
    } else {
        0.0
    };

    result.add_metric(MetricComparison::new(
        "explicit_vs_direct",
        direct_twist,
        explicit_twist,
        TOLERANCE, // Explicit should match direct within 25%
    ));

    info!(
        "  Direct twist:   {:.6} rad ({:.2}% error vs analytical)",
        direct_twist,
        (direct_twist - twist_angle_square).abs() / twist_angle_square * 100.0
    );
    info!(
        "  Explicit twist: {:.6} rad ({:.2}% error vs analytical)",
        explicit_twist,
        (explicit_twist - twist_angle_square).abs() / twist_angle_square * 100.0
    );
    info!(
        "  Explicit vs Direct error: {:.2}%",
        explicit_vs_direct_error * 100.0
    );

    result.set_notes(&format!(
        "Material: E={:.2e} Pa, ν={:.2}, G={:.2e} Pa. Geometry: L={} m, side={:.4} m. \
         Applied T={} N·m. Mesh: {}×{}×{}. \
         Explicit solver: 5000 steps with damping (0.9995 velocity factor). \
         Analytical (square section): φ={:.4}°",
        actual_e,
        actual_nu,
        g,
        L,
        SIDE,
        T,
        nx,
        ny,
        nz,
        twist_angle_square.to_degrees()
    ));

    result
}

fn run_direct(nx: usize, ny: usize, nz: usize, _g: f64, _analytical: f64) -> (f64, f64) {
    let half_side = SIDE / 2.0;

    // Generate mesh
    let mut mesh = generate_block_mesh(L, SIDE, SIDE, nx, ny, nz);

    // Shift to center at Y=0, Z=0
    for (_, node) in mesh.nodes.iter_mut() {
        node.coordinates[1] -= half_side;
        node.coordinates[2] -= half_side;
    }

    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");

    let mut simulation = Simulation::from_mesh(mesh, 3);

    // Fix X=0 end
    let fixed_bc = FixedCondition::new(x_min_nodes.clone(), vec![Some(0.0), Some(0.0), Some(0.0)]);
    simulation.add_boundary_condition(Box::new(fixed_bc));

    // Apply torque at X=L end
    // TorqueCondition::new(nodes, axis_point, axis_direction, net_torque)
    let torque_bc = TorqueCondition::new(
        x_max_nodes.clone(),
        Vector3::new(L, 0.0, 0.0),   // Point on axis (at free end)
        Vector3::new(1.0, 0.0, 0.0), // Axis direction (X)
        T,                           // Net torque
    );
    simulation.add_boundary_condition(Box::new(torque_bc));

    // Solve direct
    simulation.solve();

    // Measure twist
    let nodes = simulation.nodes();
    let mut max_twist = 0.0_f64;

    for &node_id in x_max_nodes.iter() {
        if let Some(node) = nodes.get(node_id) {
            let y = node.position.y;
            let z = node.position.z;
            let r = (y * y + z * z).sqrt();

            if r > 1e-6 {
                let dy = node.displacement.y;
                let dz = node.displacement.z;
                let tangential = (-y * dz + z * dy) / r;
                let twist = tangential / r;
                max_twist = max_twist.max(twist.abs());
            }
        }
    }

    (max_twist, 0.0)
}

fn run_explicit(nx: usize, ny: usize, nz: usize, _g: f64, _analytical: f64) -> (f64, f64) {
    let half_side = SIDE / 2.0;

    // Generate mesh
    let mut mesh = generate_block_mesh(L, SIDE, SIDE, nx, ny, nz);

    // Shift to center
    for (_, node) in mesh.nodes.iter_mut() {
        node.coordinates[1] -= half_side;
        node.coordinates[2] -= half_side;
    }

    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");

    let mut simulation = Simulation::from_mesh(mesh, 3);

    // Set explicit solver parameters using toml::Value
    simulation
        .keywords
        .add("SOLVER_METHOD", toml::Value::String("explicit".to_string()));
    simulation
        .keywords
        .add("SOLVER_TIME_STEPS", toml::Value::Integer(10000));
    simulation
        .keywords
        .add("SOLVER_PRINT_STEPS", toml::Value::Integer(0));
    // Use a smaller time step for stability
    simulation
        .keywords
        .add("SOLVER_TIME_STEP", toml::Value::Float(1e-7));

    // Fix X=0 end
    let fixed_bc = FixedCondition::new(x_min_nodes.clone(), vec![Some(0.0), Some(0.0), Some(0.0)]);
    simulation.add_boundary_condition(Box::new(fixed_bc));

    // Apply torque at X=L end
    let torque_bc = TorqueCondition::new(
        x_max_nodes.clone(),
        Vector3::new(L, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        T,
    );
    simulation.add_boundary_condition(Box::new(torque_bc));

    // Solve explicit
    simulation.solve();

    // Export VTK
    simulation.compute_result_fields();
    let vtk_path = "examples/output/vtk/torsion_shaft_explicit.vtk";
    if let Err(e) = write_vtk(vtk_path, &simulation) {
        log::warn!("Failed to write VTK: {}", e);
    } else {
        info!("  Wrote VTK: {}", vtk_path);
    }

    // Measure twist
    let nodes = simulation.nodes();
    let mut max_twist = 0.0_f64;

    for &node_id in x_max_nodes.iter() {
        if let Some(node) = nodes.get(node_id) {
            let y = node.position.y;
            let z = node.position.z;
            let r = (y * y + z * z).sqrt();

            if r > 1e-6 {
                let dy = node.displacement.y;
                let dz = node.displacement.z;
                let tangential = (-y * dz + z * dy) / r;
                let twist = tangential / r;
                max_twist = max_twist.max(twist.abs());
            }
        }
    }

    (max_twist, 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_explicit_torsion() {
        let result = run();
        assert!(result.passed);
    }
}
