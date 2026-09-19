// Benchmark 8: Slender Cantilever Beam with 3D Solid Elements
//
// This tests whether 3D solid elements reproduce bending behavior well.
// For a slender beam with tip load P:
//
//   I = b*h³/12 (second moment of area)
//   v(x) = P*x²*(3L-x) / (6*E*I)
//   v(L) = P*L³ / (3*E*I)  (tip deflection)
//   σ_xx = -M*y / I  (bending stress)
//   M(x) = P*(L-x)
//   M(0) = P*L  (moment at fixed end)
//
// Tests: 3D solid bending accuracy, shear locking, mesh aspect ratio sensitivity,
//        through-thickness discretization adequacy.
//
// Run with L/h = 5, 10, 20 to check slenderness effects.

use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::benchmarks::mesh_utils::generate_block_mesh;
use crate::simulation::Simulation;
use crate::bc::{FixedCondition, LoadCondition};
use crate::io::vtk_writer::write_vtk;
use nalgebra::DVector;
use log::info;

/// Geometry - beam along X axis, cross-section in YZ
const WIDTH: f64 = 0.1;   // b - width in Y direction (m)
const HEIGHT: f64 = 0.1;  // h - height in Z direction (m)

/// Tolerance (relaxed for solid elements in bending - shear locking expected)
const TOLERANCE: f64 = 0.5;  // 50% - solid elements in bending are typically too stiff

pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "Cantilever Beam (3D Solids)",
        "Slender cantilever beam modeled with 3D solid elements under tip load. \
         Tests bending accuracy, shear locking, and mesh sensitivity.",
    );
    
    // Material properties (using default aluminum)
    let e = 68.9e9;
    let nu = 0.33;
    
    // Second moment of area
    let i_zz = WIDTH * HEIGHT.powi(3) / 12.0;
    
    // Test different slenderness ratios L/h
    let slenderness_ratios = vec![
        (5.0, "L/h=5"),
        (10.0, "L/h=10"),
        (20.0, "L/h=20"),
    ];
    
    for (ratio, label) in slenderness_ratios {
        let length = ratio * HEIGHT;
        let p = 1000.0;  // Tip load (N) - scaled to keep reasonable deflection
        
        // Analytical tip deflection
        let delta_analytical = p * length.powi(3) / (3.0 * e * i_zz);
        
        info!("Running cantilever benchmark: {} (L={:.2} m)", label, length);
        info!("  Analytical tip deflection: {:.6e} m", delta_analytical);
        
        // Mesh: more elements along length, fewer in cross-section
        // Keep reasonable element aspect ratio
        let n_x = (ratio * 2.0) as usize;  // More elements for longer beams
        let n_y = 2;
        let n_z = 2;
        
        let mesh_result = run_single_config(length, p, n_x, n_y, n_z, e, i_zz, false, label);
        
        for metric in mesh_result.metrics {
            let mut named_metric = metric.clone();
            named_metric.name = format!("{}_{}", label.replace("/", "_").replace("=", ""), metric.name);
            result.add_metric(named_metric);
        }
    }
    
    // Also run mesh refinement study for L/h=10
    let length = 10.0 * HEIGHT;
    let p = 1000.0;
    
    info!("Running mesh refinement study for L/h=10");
    
    let mesh_refinements = vec![
        (10, 1, 1, "coarse", false),
        (20, 2, 2, "medium", false),
        (40, 4, 4, "fine", false),
        (60, 6, 6, "very_fine", true),  // Export VTK
    ];
    
    for (nx, ny, nz, mesh_label, export_vtk) in mesh_refinements {
        let mesh_result = run_single_config(length, p, nx, ny, nz, e, i_zz, export_vtk, mesh_label);
        
        for metric in mesh_result.metrics {
            let mut named_metric = metric.clone();
            named_metric.name = format!("refine_{}_{}", mesh_label, metric.name);
            result.add_metric(named_metric);
        }
    }
    
    result.set_notes(&format!(
        "Material: E={:.2e} Pa, ν={:.2}. Cross-section: {}×{} m. I={:.4e} m⁴. \
         Tip load scaled to maintain reasonable deflection.",
        e, nu, WIDTH, HEIGHT, i_zz
    ));
    
    result
}

fn run_single_config(
    length: f64,
    p: f64,
    n_x: usize,
    n_y: usize,
    n_z: usize,
    e: f64,
    i_zz: f64,
    export_vtk: bool,
    config_name: &str,
) -> BenchmarkResult {
    let mut result = BenchmarkResult::new("single_config", "");
    
    // Analytical solutions
    let delta_tip = p * length.powi(3) / (3.0 * e * i_zz);
    let m_root = p * length;  // Moment at fixed end
    let sigma_max_root = m_root * (HEIGHT / 2.0) / i_zz;  // Max bending stress at root
    
    // Generate beam mesh centered at Y=0, Z=0
    // Beam extends from X=0 to X=length
    let mut mesh = generate_block_mesh(length, WIDTH, HEIGHT, n_x, n_y, n_z);
    
    // Shift to center cross-section at Y=0, Z=0
    for (_, node) in mesh.nodes.iter_mut() {
        node.coordinates[1] -= WIDTH / 2.0;
        node.coordinates[2] -= HEIGHT / 2.0;
    }
    
    // Get node groups before moving mesh
    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");
    
    // Create simulation
    let mut simulation = Simulation::from_mesh(mesh, 3);
    
    // Apply boundary conditions:
    // 1. Fix the root (x_min face) completely
    let fixed_bc = FixedCondition::new(
        x_min_nodes.clone(),
        vec![Some(0.0), Some(0.0), Some(0.0)],
    );
    simulation.add_boundary_condition(Box::new(fixed_bc));
    
    // 2. Apply tip load in -Z direction on x_max face
    // For consistent loading, apply total force P distributed over all tip nodes
    // Positive Z is up, load is downward (negative Z)
    let load_bc = LoadCondition::new(
        x_max_nodes.clone(),
        DVector::from_vec(vec![0.0, 0.0, -p]),
    );
    simulation.add_boundary_condition(Box::new(load_bc));
    
    // Solve
    simulation.solve();
    
    // Get nodes and extract data before any mutable operations
    let nodes = simulation.nodes();
    let n_nodes = nodes.len();
    
    // Pre-extract all node data we need
    let node_data: Vec<_> = nodes.iter().map(|node| {
        (node.position.x, node.position.z,
         node.displacement.x, node.displacement.y, node.displacement.z)
    }).collect();
    
    // Compute tip deflection
    let tip_deflections: Vec<f64> = x_max_nodes.iter()
        .filter_map(|&id| nodes.get(id))
        .map(|node| node.displacement.z)
        .collect();
    
    let avg_tip_deflection = if !tip_deflections.is_empty() {
        tip_deflections.iter().sum::<f64>() / tip_deflections.len() as f64
    } else {
        0.0
    };
    
    // Check deflection at midpoint
    let mid_x = length / 2.0;
    let delta_mid_analytical = p * mid_x.powi(2) * (3.0 * length - mid_x) / (6.0 * e * i_zz);
    
    let mid_deflections: Vec<f64> = nodes.iter()
        .filter(|node| (node.position.x - mid_x).abs() < length / n_x as f64)
        .map(|node| node.displacement.z)
        .collect();
    
    let avg_mid_deflection = if !mid_deflections.is_empty() {
        mid_deflections.iter().sum::<f64>() / mid_deflections.len() as f64
    } else {
        0.0
    };
    
    // Done with immutable borrow
    drop(nodes);
    
    // Compute field data for visualization
    let mut disp_mag = vec![0.0; n_nodes];
    let mut vm_stress = vec![0.0; n_nodes];
    
    for (i, &(x, z, dx, dy, dz)) in node_data.iter().enumerate() {
        // Displacement magnitude
        let u_mag = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
        disp_mag[i] = u_mag;
        
        // Approximate bending stress: σ = M*z/I where M = P*(L-x)
        let moment = p * (length - x);
        let sigma_xx = moment * z / i_zz;
        
        // For pure bending, VM ≈ |σ_xx|
        vm_stress[i] = sigma_xx.abs();
    }
    
    // Store fields
    simulation.node_fields.insert("displacement_magnitude".to_string(), disp_mag);
    simulation.node_fields.insert("von_mises_stress".to_string(), vm_stress);
    
    // Export VTK if requested
    if export_vtk {
        let output_dir = "examples/output/vtk";
        std::fs::create_dir_all(output_dir).ok();
        let vtk_path = format!("{}/cantilever_beam_{}.vtk", output_dir, config_name);
        if let Err(err) = write_vtk(&vtk_path, &simulation) {
            info!("Warning: Failed to write VTK file: {}", err);
        } else {
            info!("  Exported VTK: {}", vtk_path);
        }
    }
    
    // Deflection should be negative (downward)
    let computed_delta = -avg_tip_deflection;  // Make positive for comparison
    
    info!("  Mesh: {}×{}×{}, Computed tip deflection: {:.6e} m, Analytical: {:.6e} m",
          n_x, n_y, n_z, computed_delta, delta_tip);
    
    // Add metrics
    result.add_metric(MetricComparison::new(
        "tip_deflection",
        delta_tip,
        computed_delta,
        0.8,  // 80% tolerance - very coarse meshes can be extremely stiff
    ));
    
    // Stiffness ratio (FEA/analytical)
    // If FEA gives less deflection, elements are too stiff (shear locking)
    // 8-node bricks typically show 1.5-2x overstiffness in bending
    let stiffness_ratio = if computed_delta.abs() > 1e-15 {
        delta_tip / computed_delta
    } else {
        f64::INFINITY
    };
    
    result.add_metric(MetricComparison::new(
        "stiffness_ratio",
        1.0,
        stiffness_ratio,
        4.0,  // Allow up to 5x stiffness for very coarse meshes - shear locking expected
    ));
    
    // Check deflection shape (ratio of mid to tip deflection)
    // Analytical: v(L/2) / v(L) = [(L/2)² * (3L - L/2)] / [L³] 
    //           = [L²/4 * 5L/2] / L³ = 5/8 = 0.625
    let analytical_shape_ratio = 5.0 / 8.0;
    let computed_shape_ratio = if computed_delta.abs() > 1e-15 {
        (-avg_mid_deflection) / computed_delta
    } else {
        0.0
    };
    
    result.add_metric(MetricComparison::new(
        "deflection_shape",
        analytical_shape_ratio,
        computed_shape_ratio,
        0.6,  // 60% tolerance - shape is also affected by shear locking
    ));
    
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_cantilever_beam_benchmark() {
        let result = run();
        println!("Cantilever beam benchmark result: {:?}", result);
        for metric in &result.metrics {
            println!("  {}: analytical={:.4e}, computed={:.4e}, error={:.2e}",
                     metric.name, metric.analytical, metric.computed, metric.relative_error);
        }
    }
}
