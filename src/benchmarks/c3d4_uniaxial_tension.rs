// Benchmark: C3D4 Uniaxial Tension
// Tests the 4-node linear tetrahedral element under uniaxial tension
// Compares displacement to analytical solution

use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::mesh::{MeshAssembly, MeshNode, MeshElement, ElementGroup};
use crate::bc::FixedCondition;
use crate::elements::Material;
use crate::simulation::Simulation;
use std::collections::HashMap;
use std::time::Instant;
use log::info;

/// Material properties (Aluminum 6061-T6 for consistency with Material::aluminum())
const E: f64 = 68.9e9;   // Young's modulus (Pa)
const NU: f64 = 0.33;    // Poisson's ratio
const DENSITY: f64 = 2700.0;  // Density (kg/m³)

/// Geometry
const LENGTH: f64 = 1.0;   // Length in x-direction (m)
const HEIGHT: f64 = 0.1;   // Height in y-direction (m)
const WIDTH: f64 = 0.1;    // Width in z-direction (m)

/// Loading
const DELTA: f64 = 0.001;  // Prescribed displacement at x_max (m)

/// Tolerance for benchmark pass/fail
/// C3D4 elements are known to be stiff due to volumetric locking
/// We allow a larger tolerance than for higher-order elements
const TOLERANCE: f64 = 0.10;  // 10% tolerance

/// Generate a block mesh of C3D4 elements
/// 
/// Creates a rectangular block meshed with tetrahedral elements.
/// Each hexahedral cell is divided into 6 tetrahedra.
/// 
/// # Arguments
/// * `length` - Length in x direction
/// * `width` - Width in y direction  
/// * `height` - Height in z direction
/// * `nx` - Number of hex cells in x direction
/// * `ny` - Number of hex cells in y direction
/// * `nz` - Number of hex cells in z direction
pub fn generate_block_mesh_c3d4(
    length: f64,
    width: f64,
    height: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> MeshAssembly {
    use crate::mesh::NodeGroup;
    
    let mut nodes: HashMap<usize, MeshNode> = HashMap::new();
    let mut elements: HashMap<usize, MeshElement> = HashMap::new();
    let mut node_groups: HashMap<String, NodeGroup> = HashMap::new();
    let mut element_groups: HashMap<String, ElementGroup> = HashMap::new();
    
    // Initialize boundary node groups
    node_groups.insert("x_min".to_string(), NodeGroup { nodes: Vec::new(), name: "x_min".to_string() });
    node_groups.insert("x_max".to_string(), NodeGroup { nodes: Vec::new(), name: "x_max".to_string() });
    node_groups.insert("y_min".to_string(), NodeGroup { nodes: Vec::new(), name: "y_min".to_string() });
    node_groups.insert("y_max".to_string(), NodeGroup { nodes: Vec::new(), name: "y_max".to_string() });
    node_groups.insert("z_min".to_string(), NodeGroup { nodes: Vec::new(), name: "z_min".to_string() });
    node_groups.insert("z_max".to_string(), NodeGroup { nodes: Vec::new(), name: "z_max".to_string() });
    
    // Spacing
    let dx = length / nx as f64;
    let dy = width / ny as f64;
    let dz = height / nz as f64;
    
    // Tolerance for boundary detection
    let tol = 1e-10;
    
    // Create nodes - grid of (nx+1) x (ny+1) x (nz+1)
    let mut node_id = 0;
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                let x = i as f64 * dx;
                let y = j as f64 * dy;
                let z = k as f64 * dz;
                
                nodes.insert(node_id, MeshNode {
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                
                // Add to boundary groups
                if x.abs() < tol {
                    node_groups.get_mut("x_min").unwrap().nodes.push(node_id);
                }
                if (x - length).abs() < tol {
                    node_groups.get_mut("x_max").unwrap().nodes.push(node_id);
                }
                if y.abs() < tol {
                    node_groups.get_mut("y_min").unwrap().nodes.push(node_id);
                }
                if (y - width).abs() < tol {
                    node_groups.get_mut("y_max").unwrap().nodes.push(node_id);
                }
                if z.abs() < tol {
                    node_groups.get_mut("z_min").unwrap().nodes.push(node_id);
                }
                if (z - height).abs() < tol {
                    node_groups.get_mut("z_max").unwrap().nodes.push(node_id);
                }
                
                node_id += 1;
            }
        }
    }
    
    // Create tetrahedral elements
    // Each hexahedral cell is divided into 6 tetrahedra
    let mut elem_id = 0;
    let mut all_elements = Vec::new();
    
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                // Get the 8 corner nodes of the hexahedral cell
                // Numbering:
                //     4-------7
                //    /|      /|
                //   / |     / |
                //  5-------6  |
                //  |  0----|--3
                //  | /     | /
                //  |/      |/
                //  1-------2
                
                let n0 = i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
                let n1 = n0 + 1;
                let n2 = n0 + 1 + (nx + 1);
                let n3 = n0 + (nx + 1);
                let n4 = n0 + (nx + 1) * (ny + 1);
                let n5 = n4 + 1;
                let n6 = n4 + 1 + (nx + 1);
                let n7 = n4 + (nx + 1);
                
                // Standard Freudenthal decomposition divides a cube into 6 tets
                let tet_connectivity = [
                    [n0, n1, n3, n5],  // tet 0
                    [n1, n2, n3, n6],  // tet 1
                    [n3, n5, n6, n7],  // tet 2
                    [n0, n3, n4, n5],  // tet 3
                    [n1, n3, n5, n6],  // tet 4
                    [n3, n4, n5, n7],  // tet 5
                ];
                
                for conn in tet_connectivity.iter() {
                    elements.insert(elem_id, MeshElement {
                        connectivity: conn.to_vec(),
                        name: format!("tet_{}", elem_id),
                        el_type: "C3D4".to_string(),
                        id: elem_id,
                    });
                    all_elements.push(elem_id);
                    elem_id += 1;
                }
            }
        }
    }
    
    // Create element group
    element_groups.insert("all_elements".to_string(), ElementGroup {
        elements: all_elements,
        el_type: "C3D4".to_string(),
        name: "all_elements".to_string(),
    });
    
    MeshAssembly {
        nodes,
        elements,
        node_groups,
        element_groups,
        bodies: Vec::new(),
        name: "c3d4_block".to_string(),
    }
}

/// Run the C3D4 uniaxial tension benchmark
pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "C3D4 Uniaxial Tension",
        "4-node linear tetrahedral element under uniaxial tension. \
         Tests displacement field against analytical solution."
    );
    
    // Run with multiple mesh refinements
    let refinements = vec![
        (2, 1, 1, "coarse"),
        (4, 2, 2, "medium"),
    ];
    
    for (nx, ny, nz, mesh_name) in refinements {
        info!("Running C3D4 uniaxial tension benchmark with {} mesh ({} x {} x {} hex cells)", 
              mesh_name, nx, ny, nz);
        
        let mesh_result = run_single_mesh(nx, ny, nz);
        
        // Add metrics with mesh name prefix
        for metric in mesh_result.metrics {
            let mut named_metric = metric.clone();
            named_metric.name = format!("{}_{}", mesh_name, metric.name);
            result.add_metric(named_metric);
        }
    }
    
    result.set_notes(&format!(
        "Material: E={:.2e} Pa, ν={:.2}. Geometry: {}×{}×{} m. Applied Δ={:.4} m. Element: C3D4",
        E, NU, LENGTH, HEIGHT, WIDTH, DELTA
    ));
    
    result
}

fn run_single_mesh(n_x: usize, n_y: usize, n_z: usize) -> BenchmarkResult {
    let mut result = BenchmarkResult::new("single_mesh", "");
    
    let start = Instant::now();
    
    // Generate C3D4 mesh
    let mesh = generate_block_mesh_c3d4(LENGTH, HEIGHT, WIDTH, n_x, n_y, n_z);
    
    // Get node groups before moving mesh into simulation
    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");
    
    // Store node coordinates for later comparison
    let node_coords: std::collections::HashMap<usize, (f64, f64, f64)> = mesh.nodes.iter()
        .map(|(id, node)| (*id, (node.coordinates[0], node.coordinates[1], node.coordinates[2])))
        .collect();
    
    // Create simulation (this moves the mesh)
    // Note: Material::aluminum() is used by default in mesh -> element conversion
    let mut sim = Simulation::from_mesh(mesh, 3);
    
    // Apply boundary conditions:
    // 1. Fix x=0 face (all directions) for stability
    let fixed_bc = FixedCondition::new(
        x_min_nodes.clone(),
        vec![Some(0.0), Some(0.0), Some(0.0)],  // Fix u_x, u_y, u_z = 0
    );
    sim.add_boundary_condition(Box::new(fixed_bc));
    
    // 2. Apply prescribed displacement at x=L face
    let prescribed_bc = FixedCondition::new(
        x_max_nodes.clone(),
        vec![Some(DELTA), None, None],  // u_x = DELTA, u_y and u_z free
    );
    sim.add_boundary_condition(Box::new(prescribed_bc));
    
    // Solve
    sim.solve();
    
    result.set_elapsed(start.elapsed().as_secs_f64() * 1000.0);
    
    // Analytical solution
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
        avg_u_x_error / DELTA,  // Normalized
        TOLERANCE,
    ));
    
    // Report mesh info
    let num_nodes = nodes.len();
    let num_elements = sim.elements().len();
    result.set_notes(&format!(
        "Mesh: {} C3D4 elements, {} nodes. Strain ε_xx = {:.6}",
        num_elements, num_nodes, strain_xx
    ));
    
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_c3d4_uniaxial_tension_benchmark() {
        let result = run();
        assert!(result.passed, "C3D4 uniaxial tension benchmark failed: {:?}", result.metrics);
    }
    
    #[test]
    fn test_c3d4_mesh_generation() {
        // Test basic mesh generation
        let mesh = generate_block_mesh_c3d4(1.0, 0.1, 0.1, 1, 1, 1);
        
        // For 1x1x1 hex cells with C3D4:
        // Nodes: (1+1)*(1+1)*(1+1) = 8 corner nodes
        // Elements: 6 tetrahedra per hex cell = 6 elements
        assert_eq!(mesh.nodes.len(), 8, "Expected 8 nodes for 1x1x1 C3D4 mesh");
        assert_eq!(mesh.elements.len(), 6, "Expected 6 elements for 1x1x1 C3D4 mesh");
        
        // Check that elements have 4 nodes
        for (_, elem) in &mesh.elements {
            assert_eq!(elem.connectivity.len(), 4, "C3D4 element should have 4 nodes");
        }
        
        // Check all node IDs in connectivity are valid
        for (_, elem) in &mesh.elements {
            for &node_id in &elem.connectivity {
                assert!(mesh.nodes.contains_key(&node_id), 
                    "Node {} in element connectivity not found in mesh", node_id);
            }
        }
        
        // Check node groups
        let x_min = mesh.get_nodes_in_group("x_min");
        let x_max = mesh.get_nodes_in_group("x_max");
        
        println!("x_min nodes: {:?}", x_min);
        println!("x_max nodes: {:?}", x_max);
        
        // x_min should have 4 nodes (the y-z face at x=0)
        assert_eq!(x_min.len(), 4, "Expected 4 nodes on x_min face");
        assert_eq!(x_max.len(), 4, "Expected 4 nodes on x_max face");
    }
}
