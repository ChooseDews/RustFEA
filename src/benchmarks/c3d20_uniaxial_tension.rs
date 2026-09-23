// Benchmark: C3D20 Uniaxial Tension
// Tests the 20-node quadratic hexahedral element under uniaxial tension
// Compares displacement to analytical solution

use crate::bc::FixedCondition;
use crate::benchmarks::mesh_utils::generate_block_mesh_c3d20;
use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::elements::Material;
use crate::simulation::Simulation;
use log::info;
use std::time::Instant;

/// Material properties (Aluminum 6061-T6 for consistency with Material::aluminum())
const E: f64 = 68.9e9; // Young's modulus (Pa)
const NU: f64 = 0.33; // Poisson's ratio
const DENSITY: f64 = 2700.0; // Density (kg/m³)

/// Geometry
const LENGTH: f64 = 1.0; // Length in x-direction (m)
const HEIGHT: f64 = 0.1; // Height in y-direction (m)
const WIDTH: f64 = 0.1; // Width in z-direction (m)

/// Loading
const DELTA: f64 = 0.001; // Prescribed displacement at x_max (m)

/// Tolerance for benchmark pass/fail
/// C3D20 elements should be very accurate for this test
const TOLERANCE: f64 = 0.02; // 2% tolerance

/// Run the C3D20 uniaxial tension benchmark
pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "C3D20 Uniaxial Tension",
        "20-node quadratic hexahedral element under uniaxial tension. \
         Tests displacement field against analytical solution.",
    );

    // Run with multiple mesh refinements
    let refinements = vec![(2, 1, 1, "coarse"), (4, 2, 2, "medium")];

    for (nx, ny, nz, mesh_name) in refinements {
        info!(
            "Running C3D20 uniaxial tension benchmark with {} mesh ({} x {} x {} elements)",
            mesh_name, nx, ny, nz
        );

        let mesh_result = run_single_mesh(nx, ny, nz);

        // Add metrics with mesh name prefix
        for metric in mesh_result.metrics {
            let mut named_metric = metric.clone();
            named_metric.name = format!("{}_{}", mesh_name, metric.name);
            result.add_metric(named_metric);
        }
    }

    result.set_notes(&format!(
        "Material: E={:.2e} Pa, ν={:.2}. Geometry: {}×{}×{} m. Applied Δ={:.4} m. Element: C3D20",
        E, NU, LENGTH, HEIGHT, WIDTH, DELTA
    ));

    result
}

fn run_single_mesh(n_x: usize, n_y: usize, n_z: usize) -> BenchmarkResult {
    let mut result = BenchmarkResult::new("single_mesh", "");

    let start = Instant::now();

    // Generate C3D20 mesh
    let mesh = generate_block_mesh_c3d20(LENGTH, HEIGHT, WIDTH, n_x, n_y, n_z);

    // Get node groups before moving mesh into simulation
    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");

    // Store node coordinates for later comparison
    let node_coords: std::collections::HashMap<usize, (f64, f64, f64)> = mesh
        .nodes
        .iter()
        .map(|(id, node)| {
            (
                *id,
                (
                    node.coordinates[0],
                    node.coordinates[1],
                    node.coordinates[2],
                ),
            )
        })
        .collect();

    // Create simulation (this moves the mesh)
    // Note: Material::aluminum() is used by default in mesh -> element conversion
    let mut sim = Simulation::from_mesh(mesh, 3);

    // Apply boundary conditions:
    // 1. Fix x=0 face (all directions) for stability
    let fixed_bc = FixedCondition::new(
        x_min_nodes.clone(),
        vec![Some(0.0), Some(0.0), Some(0.0)], // Fix u_x, u_y, u_z = 0
    );
    sim.add_boundary_condition(Box::new(fixed_bc));

    // 2. Apply prescribed displacement at x=L face
    let prescribed_bc = FixedCondition::new(
        x_max_nodes.clone(),
        vec![Some(DELTA), None, None], // u_x = DELTA, u_y and u_z free
    );
    sim.add_boundary_condition(Box::new(prescribed_bc));

    // Solve
    sim.solve();

    result.set_elapsed(start.elapsed().as_secs_f64() * 1000.0);

    // Analytical solution:
    // Due to fixed BC at x=0 (all directions), Poisson contraction is partially suppressed
    // For a free-to-contract bar: u_y = -ν * ε_xx * y, u_z = -ν * ε_xx * z
    // With fixed y,z at x=0, the contraction is still present but constrained

    let strain_xx = DELTA / LENGTH;

    // Get nodes for verification
    let nodes = sim.nodes();

    // Calculate average displacement error across all nodes for u_x
    // For u_x, analytical is: u_x = (Δ/L) * x
    let mut total_u_x_error = 0.0;
    let mut count = 0;

    for (node_id, &(x, _y, _z)) in &node_coords {
        let node = nodes.get(*node_id).unwrap();
        let u_x_analytical = strain_xx * x;
        total_u_x_error += (node.displacement.x - u_x_analytical).abs();
        count += 1;
    }

    let avg_u_x_error = total_u_x_error / count as f64;

    // Check prescribed displacement at x_max
    let sample_x_max_node = x_max_nodes.first().cloned();
    let sample_u_x = sample_x_max_node
        .and_then(|id| nodes.get(id))
        .map(|n| n.displacement.x)
        .unwrap_or(0.0);

    result.add_metric(MetricComparison::new(
        "u_x_prescribed",
        DELTA,
        sample_u_x,
        TOLERANCE,
    ));

    result.add_metric(MetricComparison::new(
        "avg_u_x_error",
        0.0,
        avg_u_x_error / DELTA, // Normalized
        TOLERANCE,
    ));

    // Report mesh info
    let num_nodes = nodes.len();
    let num_elements = sim.elements().len();
    result.set_notes(&format!(
        "Mesh: {} C3D20 elements, {} nodes. Strain ε_xx = {:.6}",
        num_elements, num_nodes, strain_xx
    ));

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c3d20_uniaxial_tension_benchmark() {
        let result = run();
        assert!(
            result.passed,
            "C3D20 uniaxial tension benchmark failed: {:?}",
            result.metrics
        );
    }

    #[test]
    fn test_c3d20_mesh_generation() {
        // Test basic mesh generation
        let mesh = generate_block_mesh_c3d20(1.0, 0.1, 0.1, 1, 1, 1);

        // For 1x1x1 elements with C3D20:
        // C3D20 has 20 nodes: 8 corners + 12 mid-edge nodes
        // (No face-center or body-center nodes)
        assert_eq!(
            mesh.nodes.len(),
            20,
            "Expected 20 nodes for 1x1x1 C3D20 mesh"
        );
        assert_eq!(
            mesh.elements.len(),
            1,
            "Expected 1 element for 1x1x1 C3D20 mesh"
        );

        // Check that element has 20 nodes
        let elem = mesh.elements.get(&0).unwrap();
        assert_eq!(
            elem.connectivity.len(),
            20,
            "C3D20 element should have 20 nodes"
        );

        // Check all node IDs in connectivity are valid
        for &node_id in &elem.connectivity {
            assert!(
                mesh.nodes.contains_key(&node_id),
                "Node {} in element connectivity not found in mesh",
                node_id
            );
        }

        // Check node groups
        let x_min = mesh.get_nodes_in_group("x_min");
        let x_max = mesh.get_nodes_in_group("x_max");

        println!("x_min nodes: {:?}", x_min);
        println!("x_max nodes: {:?}", x_max);
        println!("Element connectivity: {:?}", elem.connectivity);

        // Verify x_min and x_max nodes are in the element
        for &node_id in &x_min {
            assert!(
                elem.connectivity.contains(&node_id),
                "x_min node {} not in element connectivity",
                node_id
            );
        }
        for &node_id in &x_max {
            assert!(
                elem.connectivity.contains(&node_id),
                "x_max node {} not in element connectivity",
                node_id
            );
        }
    }
}
