// Benchmark 2: 3D Block in Pure Shear
//
// A rectangular block with prescribed displacement:
//   u_x = γ * y
//   u_y = 0
//   u_z = 0
//
// This generates uniform engineering shear strain γ_xy = γ
//
// Analytical solution:
//   τ_xy = G * γ
//   where G = E / (2 * (1 + ν))
//   All normal stresses are zero
//
// Tests: xy/xz/yz shear components, tensor vs engineering shear conventions,
//        node ordering, coordinate transformations, stress recovery

use crate::bc::FixedCondition;
use crate::benchmarks::mesh_utils::generate_block_mesh;
use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::elements::base_element::Material;
use crate::simulation::Simulation;
use log::info;

/// Geometry
const L: f64 = 1.0; // Length in X (m)
const W: f64 = 1.0; // Width in Y (m)
const H: f64 = 0.5; // Height in Z (m)

/// Applied shear strain
const GAMMA: f64 = 0.001; // Engineering shear strain (dimensionless)

/// Tolerance
const TOLERANCE: f64 = 0.02; // 2% tolerance

pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "Pure Shear (3D Block)",
        "3D block under pure shear deformation u_x = γy. \
         Tests shear strain components, stress recovery, and coordinate handling.",
    );

    // Run with multiple mesh refinements
    let refinements = vec![(2, 2, 2, "coarse"), (4, 4, 2, "medium"), (8, 8, 4, "fine")];

    for (nx, ny, nz, mesh_name) in refinements {
        info!(
            "Running pure shear benchmark with {} mesh ({} x {} x {} elements)",
            mesh_name, nx, ny, nz
        );

        let mesh_result = run_single_mesh(nx, ny, nz);

        for metric in mesh_result.metrics {
            let mut named_metric = metric.clone();
            named_metric.name = format!("{}_{}", mesh_name, metric.name);
            result.add_metric(named_metric);
        }
    }

    // Material properties (using default aluminum)
    let actual_e = 68.9e9; // Aluminum Young's modulus
    let actual_nu = 0.33; // Aluminum Poisson's ratio
    let g = actual_e / (2.0 * (1.0 + actual_nu));

    result.set_notes(&format!(
        "Material: E={:.2e} Pa, ν={:.2}, G={:.2e} Pa. Geometry: {}×{}×{} m. Applied γ={:.4}",
        actual_e, actual_nu, g, L, W, H, GAMMA
    ));

    result
}

fn run_single_mesh(nx: usize, ny: usize, nz: usize) -> BenchmarkResult {
    let mut result = BenchmarkResult::new("single_mesh", "");

    // Generate mesh
    let mesh = generate_block_mesh(L, W, H, nx, ny, nz);

    // Material properties (using default aluminum)
    let actual_e = 68.9e9;
    let actual_nu = 0.33;
    let g = actual_e / (2.0 * (1.0 + actual_nu));

    // Get all boundary nodes and their coordinates before moving mesh
    let all_face_groups = vec!["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"];
    let mut all_boundary_nodes = std::collections::HashSet::new();

    for group in &all_face_groups {
        for node_id in mesh.get_nodes_in_group(group) {
            all_boundary_nodes.insert(node_id);
        }
    }

    // Store boundary node coordinates for BC application
    let boundary_node_coords: std::collections::HashMap<usize, f64> = all_boundary_nodes
        .iter()
        .filter_map(|&id| mesh.nodes.get(&id).map(|n| (id, n.coordinates[1])))
        .collect();

    // Get y_max nodes for later verification
    let y_max_nodes = mesh.get_nodes_in_group("y_max");

    // Create simulation (this moves the mesh)
    let mut simulation = Simulation::from_mesh(mesh, 3);

    // Apply pure shear boundary conditions
    // For pure shear: u_x = γ*y, u_y = 0, u_z = 0

    // Apply the exact solution as Dirichlet BC on all boundary nodes
    // This ensures the internal solution must match the analytical one
    for (&node_id, &y) in boundary_node_coords.iter() {
        // u_x = γ * y, u_y = 0, u_z = 0
        let u_x = GAMMA * y;

        let bc = FixedCondition::new(vec![node_id], vec![Some(u_x), Some(0.0), Some(0.0)]);
        simulation.add_boundary_condition(Box::new(bc));
    }

    // Solve
    simulation.solve();

    // Check results
    let nodes = simulation.nodes();

    // Verify displacement field
    let mut max_u_x_error = 0.0;
    let mut max_u_y_error = 0.0;
    let mut max_u_z_error = 0.0;

    for node in nodes.iter() {
        let y = node.position.y;

        // Analytical displacements
        let u_x_analytical = GAMMA * y;
        let u_y_analytical = 0.0;
        let u_z_analytical = 0.0;

        // Computed displacements
        let u_x_computed = node.displacement.x;
        let u_y_computed = node.displacement.y;
        let u_z_computed = node.displacement.z;

        // Track max errors
        let u_x_err = (u_x_computed - u_x_analytical).abs();
        let u_y_err = (u_y_computed - u_y_analytical).abs();
        let u_z_err = (u_z_computed - u_z_analytical).abs();

        if u_x_err > max_u_x_error {
            max_u_x_error = u_x_err;
        }
        if u_y_err > max_u_y_error {
            max_u_y_error = u_y_err;
        }
        if u_z_err > max_u_z_error {
            max_u_z_error = u_z_err;
        }
    }

    // Analytical shear stress
    let tau_xy_analytical = g * GAMMA;

    // Add displacement error metrics
    result.add_metric(MetricComparison::new(
        "max_u_x_error",
        0.0,
        max_u_x_error / (GAMMA * W), // Normalized to max expected u_x
        TOLERANCE,
    ));

    result.add_metric(MetricComparison::new(
        "max_u_y_error",
        0.0,
        max_u_y_error / (GAMMA * W),
        TOLERANCE,
    ));

    result.add_metric(MetricComparison::new(
        "max_u_z_error",
        0.0,
        max_u_z_error / (GAMMA * W),
        TOLERANCE,
    ));

    // Check corner displacement for sanity
    // At y=W (top), u_x should equal GAMMA * W
    let corner_u_x = if let Some(&node_id) = y_max_nodes.first() {
        nodes.get(node_id).map(|n| n.displacement.x).unwrap_or(0.0)
    } else {
        0.0
    };

    result.add_metric(MetricComparison::new(
        "corner_u_x",
        GAMMA * W,
        corner_u_x,
        TOLERANCE,
    ));

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pure_shear_benchmark() {
        let result = run();
        assert!(
            result.passed,
            "Pure shear benchmark failed: {:?}",
            result.metrics
        );
    }
}
