use crate::bc::{BodyForce, FixedCondition};
/// Gravity Benchmark: Cantilever beam under self-weight
///
/// Tests the BodyForce boundary condition by comparing FEA results to analytical
/// solution for a cantilever beam deflecting under its own weight.
///
/// # Analytical Solution
/// Maximum deflection at free end: δ_max = qL⁴ / (8EI)
/// where:
/// - q = distributed load per unit length = ρAg
/// - L = beam length
/// - E = Young's modulus
/// - I = second moment of area = bh³/12 for rectangular cross-section
use crate::benchmarks::mesh_utils::generate_block_mesh;
use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::simulation::Simulation;
use log::info;

pub fn run() -> BenchmarkResult {
    info!("Running Gravity Benchmark (Self-Weight Cantilever)");

    let mut result = BenchmarkResult::new(
        "Self-Weight Cantilever",
        "Tests body force BC with gravity on cantilever beam",
    );

    // Beam geometry
    let length: f64 = 1.0; // m
    let width: f64 = 0.1; // m
    let height: f64 = 0.1; // m

    // Material properties (Aluminum - same as other benchmarks)
    let e: f64 = 68.9e9; // Pa
    let nu: f64 = 0.33;
    let density: f64 = 2700.0; // kg/m³

    // Gravitational acceleration
    let g: f64 = 9.81; // m/s²

    // Analytical solution
    // I = bh³/12 for rectangular cross-section
    let i_moment = width * height.powi(3) / 12.0;
    // q = distributed load per unit length = ρ * A * g
    let area = width * height;
    let q = density * area * g;
    // Maximum deflection at free end: δ = qL⁴/(8EI)
    let delta_analytical = q * length.powi(4) / (8.0 * e * i_moment);

    info!("Analytical max deflection: {:.6e} m", delta_analytical);
    info!("Load per unit length q = {:.2} N/m", q);

    // Generate mesh - beam along X axis, 20x4x4 elements
    let mut mesh = generate_block_mesh(length, width, height, 20, 4, 4);

    // Shift to center cross-section at Y=0, Z=0 (like cantilever_beam benchmark)
    for (_, node) in mesh.nodes.iter_mut() {
        node.coordinates[1] -= width / 2.0;
        node.coordinates[2] -= height / 2.0;
    }

    // Get fixed nodes (x = 0 face) - must get before moving mesh to simulation
    let fixed_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");
    info!("Fixed nodes at x=0: {}", fixed_nodes.len());

    // Create simulation
    let mut simulation = Simulation::from_mesh(mesh, 3);

    // Apply fixed BC at one end (x = 0) - use FixedCondition::new pattern
    simulation.add_boundary_condition(Box::new(FixedCondition::new(
        fixed_nodes,
        vec![Some(0.0), Some(0.0), Some(0.0)],
    )));

    // Apply gravity body force (negative Z direction)
    simulation.add_boundary_condition(Box::new(BodyForce::gravity_from_vec(0.0, 0.0, -g)));

    // Solve
    simulation.solve();

    // Find maximum Z-displacement at free end (x = length)
    let tip_deflections: Vec<f64> = x_max_nodes
        .iter()
        .filter_map(|&id| simulation.nodes().get(id))
        .map(|node| node.displacement.z.abs())
        .collect();

    let max_z_disp = tip_deflections.iter().cloned().fold(0.0_f64, f64::max);

    info!("FEA max Z displacement: {:.6e} m", max_z_disp);

    // Calculate error
    let error = ((max_z_disp - delta_analytical) / delta_analytical).abs() * 100.0;
    info!("Error: {:.2}%", error);

    result.add_metric(MetricComparison::new(
        "Max Z Displacement",
        delta_analytical,
        max_z_disp,
        0.25, // 25% tolerance (solid elements in bending are stiff due to shear locking)
    ));

    result.set_notes(&format!(
        "Beam: {}m × {}m × {}m, Aluminum (E={:.1e} Pa, ρ={} kg/m³)\n\
         Load: q = ρAg = {:.2} N/m\n\
         Analytical δ_max = qL⁴/(8EI) = {:.4e} m\n\
         FEA δ_max = {:.4e} m\n\
         Error: {:.2}%",
        length, width, height, e, density, q, delta_analytical, max_z_disp, error
    ));

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gravity_benchmark() {
        let result = run();
        assert!(result.passed, "Gravity benchmark failed");
    }
}
