// Mesh generation utilities for benchmarks
// Provides functions to create simple analytical meshes for validation

use crate::mesh::{MeshAssembly, MeshNode, MeshElement, NodeGroup, ElementGroup, Body};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Generate a rectangular block mesh with specified divisions
/// Returns a MeshAssembly with node groups for each face
pub fn generate_block_mesh(
    l_x: f64, l_y: f64, l_z: f64,
    n_x: usize, n_y: usize, n_z: usize,
) -> MeshAssembly {
    let dx = l_x / n_x as f64;
    let dy = l_y / n_y as f64;
    let dz = l_z / n_z as f64;
    
    let mut nodes = HashMap::new();
    let mut elements = HashMap::new();
    let mut node_groups: HashMap<String, NodeGroup> = HashMap::new();
    let mut element_groups: HashMap<String, ElementGroup> = HashMap::new();
    
    // Create face node groups
    let mut x_min_nodes = Vec::new();
    let mut x_max_nodes = Vec::new();
    let mut y_min_nodes = Vec::new();
    let mut y_max_nodes = Vec::new();
    let mut z_min_nodes = Vec::new();
    let mut z_max_nodes = Vec::new();
    
    // Generate nodes
    let mut node_id = 0;
    for k in 0..=n_z {
        for j in 0..=n_y {
            for i in 0..=n_x {
                let x = i as f64 * dx;
                let y = j as f64 * dy;
                let z = k as f64 * dz;
                
                nodes.insert(node_id, MeshNode {
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                
                // Categorize by face
                if i == 0 { x_min_nodes.push(node_id); }
                if i == n_x { x_max_nodes.push(node_id); }
                if j == 0 { y_min_nodes.push(node_id); }
                if j == n_y { y_max_nodes.push(node_id); }
                if k == 0 { z_min_nodes.push(node_id); }
                if k == n_z { z_max_nodes.push(node_id); }
                
                node_id += 1;
            }
        }
    }
    
    // Generate elements
    let mut elem_id = 0;
    let mut all_elements = Vec::new();
    
    for k in 0..n_z {
        for j in 0..n_y {
            for i in 0..n_x {
                let n0 = i + j * (n_x + 1) + k * (n_x + 1) * (n_y + 1);
                let n1 = n0 + 1;
                let n2 = n0 + (n_x + 1) + 1;
                let n3 = n0 + (n_x + 1);
                let n4 = n0 + (n_x + 1) * (n_y + 1);
                let n5 = n4 + 1;
                let n6 = n4 + (n_x + 1) + 1;
                let n7 = n4 + (n_x + 1);
                
                elements.insert(elem_id, MeshElement {
                    connectivity: vec![n0, n1, n2, n3, n4, n5, n6, n7],
                    name: "block_body".to_string(),
                    el_type: "C3D8".to_string(),
                    id: elem_id,
                });
                
                all_elements.push(elem_id);
                elem_id += 1;
            }
        }
    }
    
    // Create node groups
    node_groups.insert("x_min".to_string(), NodeGroup { nodes: x_min_nodes.clone(), name: "x_min".to_string() });
    node_groups.insert("x_max".to_string(), NodeGroup { nodes: x_max_nodes.clone(), name: "x_max".to_string() });
    node_groups.insert("y_min".to_string(), NodeGroup { nodes: y_min_nodes.clone(), name: "y_min".to_string() });
    node_groups.insert("y_max".to_string(), NodeGroup { nodes: y_max_nodes.clone(), name: "y_max".to_string() });
    node_groups.insert("z_min".to_string(), NodeGroup { nodes: z_min_nodes.clone(), name: "z_min".to_string() });
    node_groups.insert("z_max".to_string(), NodeGroup { nodes: z_max_nodes.clone(), name: "z_max".to_string() });
    
    // Create element group
    element_groups.insert("block_body".to_string(), ElementGroup {
        elements: all_elements.clone(),
        name: "block_body".to_string(),
        el_type: "C3D8".to_string(),
    });
    
    // Create body
    let all_nodes: Vec<usize> = nodes.keys().cloned().collect();
    let body = Body {
        elements: all_elements,
        nodes: all_nodes,
        name: "block_body".to_string(),
    };
    
    MeshAssembly {
        nodes,
        elements,
        element_groups,
        node_groups,
        bodies: vec![body],
        name: "block_mesh".to_string(),
    }
}

/// Generate a rectangular block mesh with C3D20 (20-node quadratic hexahedral) elements
/// Returns a MeshAssembly with node groups for each face
/// 
/// Note: C3D20 (serendipity) elements only have corner and mid-edge nodes,
/// NOT face-center or body-center nodes. This generator only creates the 
/// nodes that are actually used by the elements.
pub fn generate_block_mesh_c3d20(
    l_x: f64, l_y: f64, l_z: f64,
    n_x: usize, n_y: usize, n_z: usize,
) -> MeshAssembly {
    // For C3D20 elements, we use a 2n+1 grid per direction for reference,
    // but only create nodes at corner and mid-edge positions
    let nn_x = 2 * n_x + 1;
    let nn_y = 2 * n_y + 1;
    let nn_z = 2 * n_z + 1;
    
    let dx = l_x / (2 * n_x) as f64;
    let dy = l_y / (2 * n_y) as f64;
    let dz = l_z / (2 * n_z) as f64;
    
    let mut nodes = HashMap::new();
    let mut elements = HashMap::new();
    let mut node_groups: HashMap<String, NodeGroup> = HashMap::new();
    let mut element_groups: HashMap<String, ElementGroup> = HashMap::new();
    
    // Helper to compute grid index (for reference)
    let grid_index = |i: usize, j: usize, k: usize| -> usize {
        i + j * nn_x + k * nn_x * nn_y
    };
    
    // First pass: collect all unique node grid positions used by C3D20 elements
    let mut grid_positions = std::collections::BTreeSet::new();
    
    for k in 0..n_z {
        for j in 0..n_y {
            for i in 0..n_x {
                let i0 = 2 * i;
                let j0 = 2 * j;
                let k0 = 2 * k;
                
                // Corner nodes
                grid_positions.insert(grid_index(i0, j0, k0));
                grid_positions.insert(grid_index(i0 + 2, j0, k0));
                grid_positions.insert(grid_index(i0 + 2, j0 + 2, k0));
                grid_positions.insert(grid_index(i0, j0 + 2, k0));
                grid_positions.insert(grid_index(i0, j0, k0 + 2));
                grid_positions.insert(grid_index(i0 + 2, j0, k0 + 2));
                grid_positions.insert(grid_index(i0 + 2, j0 + 2, k0 + 2));
                grid_positions.insert(grid_index(i0, j0 + 2, k0 + 2));
                
                // Mid-edge nodes on bottom and top faces
                grid_positions.insert(grid_index(i0 + 1, j0, k0));
                grid_positions.insert(grid_index(i0 + 2, j0 + 1, k0));
                grid_positions.insert(grid_index(i0 + 1, j0 + 2, k0));
                grid_positions.insert(grid_index(i0, j0 + 1, k0));
                grid_positions.insert(grid_index(i0 + 1, j0, k0 + 2));
                grid_positions.insert(grid_index(i0 + 2, j0 + 1, k0 + 2));
                grid_positions.insert(grid_index(i0 + 1, j0 + 2, k0 + 2));
                grid_positions.insert(grid_index(i0, j0 + 1, k0 + 2));
                
                // Mid-edge nodes on vertical edges
                grid_positions.insert(grid_index(i0, j0, k0 + 1));
                grid_positions.insert(grid_index(i0 + 2, j0, k0 + 1));
                grid_positions.insert(grid_index(i0 + 2, j0 + 2, k0 + 1));
                grid_positions.insert(grid_index(i0, j0 + 2, k0 + 1));
            }
        }
    }
    
    // Create mapping from grid index to consecutive node ID
    let mut grid_to_node: HashMap<usize, usize> = HashMap::new();
    let mut node_id = 0;
    for &grid_idx in &grid_positions {
        grid_to_node.insert(grid_idx, node_id);
        node_id += 1;
    }
    
    // Create nodes with consecutive IDs
    let mut x_min_nodes = Vec::new();
    let mut x_max_nodes = Vec::new();
    let mut y_min_nodes = Vec::new();
    let mut y_max_nodes = Vec::new();
    let mut z_min_nodes = Vec::new();
    let mut z_max_nodes = Vec::new();
    
    for &grid_idx in &grid_positions {
        let node_id = grid_to_node[&grid_idx];
        
        // Compute coordinates from grid index
        let i = grid_idx % nn_x;
        let j = (grid_idx / nn_x) % nn_y;
        let k = grid_idx / (nn_x * nn_y);
        
        let x = i as f64 * dx;
        let y = j as f64 * dy;
        let z = k as f64 * dz;
        
        nodes.insert(node_id, MeshNode {
            coordinates: vec![x, y, z],
            id: node_id,
        });
        
        // Categorize by face
        if x.abs() < 1e-10 { x_min_nodes.push(node_id); }
        if (x - l_x).abs() < 1e-10 { x_max_nodes.push(node_id); }
        if y.abs() < 1e-10 { y_min_nodes.push(node_id); }
        if (y - l_y).abs() < 1e-10 { y_max_nodes.push(node_id); }
        if z.abs() < 1e-10 { z_min_nodes.push(node_id); }
        if (z - l_z).abs() < 1e-10 { z_max_nodes.push(node_id); }
    }
    
    // Generate C3D20 elements with remapped node IDs
    let mut elem_id = 0;
    let mut all_elements = Vec::new();
    
    for k in 0..n_z {
        for j in 0..n_y {
            for i in 0..n_x {
                let i0 = 2 * i;
                let j0 = 2 * j;
                let k0 = 2 * k;
                
                // Map grid indices to consecutive node IDs
                let n0 = grid_to_node[&grid_index(i0, j0, k0)];
                let n1 = grid_to_node[&grid_index(i0 + 2, j0, k0)];
                let n2 = grid_to_node[&grid_index(i0 + 2, j0 + 2, k0)];
                let n3 = grid_to_node[&grid_index(i0, j0 + 2, k0)];
                let n4 = grid_to_node[&grid_index(i0, j0, k0 + 2)];
                let n5 = grid_to_node[&grid_index(i0 + 2, j0, k0 + 2)];
                let n6 = grid_to_node[&grid_index(i0 + 2, j0 + 2, k0 + 2)];
                let n7 = grid_to_node[&grid_index(i0, j0 + 2, k0 + 2)];
                
                let n8 = grid_to_node[&grid_index(i0 + 1, j0, k0)];
                let n9 = grid_to_node[&grid_index(i0 + 2, j0 + 1, k0)];
                let n10 = grid_to_node[&grid_index(i0 + 1, j0 + 2, k0)];
                let n11 = grid_to_node[&grid_index(i0, j0 + 1, k0)];
                let n12 = grid_to_node[&grid_index(i0 + 1, j0, k0 + 2)];
                let n13 = grid_to_node[&grid_index(i0 + 2, j0 + 1, k0 + 2)];
                let n14 = grid_to_node[&grid_index(i0 + 1, j0 + 2, k0 + 2)];
                let n15 = grid_to_node[&grid_index(i0, j0 + 1, k0 + 2)];
                
                let n16 = grid_to_node[&grid_index(i0, j0, k0 + 1)];
                let n17 = grid_to_node[&grid_index(i0 + 2, j0, k0 + 1)];
                let n18 = grid_to_node[&grid_index(i0 + 2, j0 + 2, k0 + 1)];
                let n19 = grid_to_node[&grid_index(i0, j0 + 2, k0 + 1)];
                
                elements.insert(elem_id, MeshElement {
                    connectivity: vec![
                        n0, n1, n2, n3, n4, n5, n6, n7,
                        n8, n9, n10, n11, n12, n13, n14, n15,
                        n16, n17, n18, n19
                    ],
                    name: "block_body".to_string(),
                    el_type: "C3D20".to_string(),
                    id: elem_id,
                });
                
                all_elements.push(elem_id);
                elem_id += 1;
            }
        }
    }
    
    // Sort node groups for consistent ordering
    x_min_nodes.sort();
    x_max_nodes.sort();
    y_min_nodes.sort();
    y_max_nodes.sort();
    z_min_nodes.sort();
    z_max_nodes.sort();
    
    node_groups.insert("x_min".to_string(), NodeGroup { nodes: x_min_nodes.clone(), name: "x_min".to_string() });
    node_groups.insert("x_max".to_string(), NodeGroup { nodes: x_max_nodes.clone(), name: "x_max".to_string() });
    node_groups.insert("y_min".to_string(), NodeGroup { nodes: y_min_nodes.clone(), name: "y_min".to_string() });
    node_groups.insert("y_max".to_string(), NodeGroup { nodes: y_max_nodes.clone(), name: "y_max".to_string() });
    node_groups.insert("z_min".to_string(), NodeGroup { nodes: z_min_nodes.clone(), name: "z_min".to_string() });
    node_groups.insert("z_max".to_string(), NodeGroup { nodes: z_max_nodes.clone(), name: "z_max".to_string() });
    
    // Create element group
    element_groups.insert("block_body".to_string(), ElementGroup {
        elements: all_elements.clone(),
        name: "block_body".to_string(),
        el_type: "C3D20".to_string(),
    });
    
    // Create body
    let all_nodes: Vec<usize> = nodes.keys().cloned().collect();
    let body = Body {
        elements: all_elements,
        nodes: all_nodes,
        name: "block_body".to_string(),
    };
    
    MeshAssembly {
        nodes,
        elements,
        element_groups,
        node_groups,
        bodies: vec![body],
        name: "block_mesh_c3d20".to_string(),
    }
}

/// Generate a cylindrical mesh (for torsion shaft)
/// Axis aligned with Z direction, centered at origin
pub fn generate_cylinder_mesh(
    radius: f64,
    length: f64,
    n_radial: usize,    // divisions in radial direction
    n_circumferential: usize, // divisions around circumference
    n_axial: usize,     // divisions along length
) -> MeshAssembly {
    let mut nodes = HashMap::new();
    let mut elements = HashMap::new();
    let mut node_groups: HashMap<String, NodeGroup> = HashMap::new();
    let mut element_groups: HashMap<String, ElementGroup> = HashMap::new();
    
    let dr = radius / n_radial as f64;
    let d_theta = 2.0 * PI / n_circumferential as f64;
    let dz = length / n_axial as f64;
    
    let mut z_min_nodes = Vec::new();
    let mut z_max_nodes = Vec::new();
    let mut outer_nodes = Vec::new();
    
    // Generate nodes
    let mut node_id = 0;
    
    // Center axis nodes (r=0)
    for k in 0..=n_axial {
        let z = k as f64 * dz;
        nodes.insert(node_id, MeshNode {
            coordinates: vec![0.0, 0.0, z],
            id: node_id,
        });
        
        if k == 0 { z_min_nodes.push(node_id); }
        if k == n_axial { z_max_nodes.push(node_id); }
        
        node_id += 1;
    }
    
    // Ring nodes
    for i in 1..=n_radial {
        let r = i as f64 * dr;
        for j in 0..n_circumferential {
            let theta = j as f64 * d_theta;
            let x = r * theta.cos();
            let y = r * theta.sin();
            
            for k in 0..=n_axial {
                let z = k as f64 * dz;
                nodes.insert(node_id, MeshNode {
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                
                if k == 0 { z_min_nodes.push(node_id); }
                if k == n_axial { z_max_nodes.push(node_id); }
                if i == n_radial { outer_nodes.push(node_id); }
                
                node_id += 1;
            }
        }
    }
    
    // For simplicity, we'll use a simpler approach - generate a blocky approximation
    // This is a simplified version; a production version would use proper wedge elements
    
    // Create node groups
    node_groups.insert("z_min".to_string(), NodeGroup { nodes: z_min_nodes.clone(), name: "z_min".to_string() });
    node_groups.insert("z_max".to_string(), NodeGroup { nodes: z_max_nodes.clone(), name: "z_max".to_string() });
    node_groups.insert("outer".to_string(), NodeGroup { nodes: outer_nodes.clone(), name: "outer".to_string() });
    
    // Note: Element generation for cylinders is complex with hex elements
    // For the benchmark, we'll use the block mesh and apply appropriate BCs
    
    let all_nodes: Vec<usize> = nodes.keys().cloned().collect();
    let body = Body {
        elements: vec![],
        nodes: all_nodes,
        name: "cylinder_body".to_string(),
    };
    
    MeshAssembly {
        nodes,
        elements,
        element_groups,
        node_groups,
        bodies: vec![body],
        name: "cylinder_mesh".to_string(),
    }
}

/// Generate a hollow spherical mesh (for pressure vessel)
/// Using a sector approach with 1/8 symmetry
pub fn generate_hollow_sphere_sector_mesh(
    inner_radius: f64,
    outer_radius: f64,
    n_radial: usize,
    n_phi: usize,    // divisions in phi (from z-axis)
    n_theta: usize,  // divisions in theta (around z-axis)
    phi_max: f64,    // max phi angle (PI/2 for quarter)
    theta_max: f64,  // max theta angle (PI/2 for quarter)
) -> MeshAssembly {
    let mut nodes = HashMap::new();
    let mut elements = HashMap::new();
    let mut node_groups: HashMap<String, NodeGroup> = HashMap::new();
    let mut element_groups: HashMap<String, ElementGroup> = HashMap::new();
    
    let dr = (outer_radius - inner_radius) / n_radial as f64;
    let d_phi = phi_max / n_phi as f64;
    let d_theta = theta_max / n_theta as f64;
    
    let mut inner_nodes = Vec::new();
    let mut outer_nodes = Vec::new();
    let mut x_sym_nodes = Vec::new();  // theta = 0 plane
    let mut y_sym_nodes = Vec::new();  // theta = theta_max plane
    let mut z_sym_nodes = Vec::new();  // phi = phi_max plane
    
    // Generate nodes in spherical coordinates
    let mut node_id = 0;
    for i in 0..=n_radial {
        let r = inner_radius + i as f64 * dr;
        for j in 0..=n_phi {
            let phi = j as f64 * d_phi;
            for k in 0..=n_theta {
                let theta = k as f64 * d_theta;
                
                let x = r * phi.sin() * theta.cos();
                let y = r * phi.sin() * theta.sin();
                let z = r * phi.cos();
                
                nodes.insert(node_id, MeshNode {
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                
                // Categorize by surface
                if i == 0 { inner_nodes.push(node_id); }
                if i == n_radial { outer_nodes.push(node_id); }
                if k == 0 { x_sym_nodes.push(node_id); }  // x-z plane
                if k == n_theta { y_sym_nodes.push(node_id); }  // y-z plane (if theta_max = PI/2)
                if j == n_phi { z_sym_nodes.push(node_id); }  // x-y plane (if phi_max = PI/2)
                
                node_id += 1;
            }
        }
    }
    
    // Generate hex elements
    let mut elem_id = 0;
    let mut all_elements = Vec::new();
    
    let nodes_per_ring = n_theta + 1;
    let nodes_per_shell = (n_phi + 1) * nodes_per_ring;
    
    for i in 0..n_radial {
        for j in 0..n_phi {
            for k in 0..n_theta {
                // Node indices for inner shell
                let n0 = i * nodes_per_shell + j * nodes_per_ring + k;
                let n1 = n0 + 1;
                let n2 = n0 + nodes_per_ring + 1;
                let n3 = n0 + nodes_per_ring;
                
                // Node indices for outer shell
                let n4 = n0 + nodes_per_shell;
                let n5 = n1 + nodes_per_shell;
                let n6 = n2 + nodes_per_shell;
                let n7 = n3 + nodes_per_shell;
                
                elements.insert(elem_id, MeshElement {
                    connectivity: vec![n0, n1, n2, n3, n4, n5, n6, n7],
                    name: "sphere_body".to_string(),
                    el_type: "C3D8".to_string(),
                    id: elem_id,
                });
                
                all_elements.push(elem_id);
                elem_id += 1;
            }
        }
    }
    
    // Create node groups
    node_groups.insert("inner".to_string(), NodeGroup { nodes: inner_nodes.clone(), name: "inner".to_string() });
    node_groups.insert("outer".to_string(), NodeGroup { nodes: outer_nodes.clone(), name: "outer".to_string() });
    node_groups.insert("x_sym".to_string(), NodeGroup { nodes: x_sym_nodes.clone(), name: "x_sym".to_string() });
    node_groups.insert("y_sym".to_string(), NodeGroup { nodes: y_sym_nodes.clone(), name: "y_sym".to_string() });
    node_groups.insert("z_sym".to_string(), NodeGroup { nodes: z_sym_nodes.clone(), name: "z_sym".to_string() });
    
    // Create element group
    element_groups.insert("sphere_body".to_string(), ElementGroup {
        elements: all_elements.clone(),
        name: "sphere_body".to_string(),
        el_type: "C3D8".to_string(),
    });
    
    let all_nodes: Vec<usize> = nodes.keys().cloned().collect();
    let body = Body {
        elements: all_elements,
        nodes: all_nodes,
        name: "sphere_body".to_string(),
    };
    
    MeshAssembly {
        nodes,
        elements,
        element_groups,
        node_groups,
        bodies: vec![body],
        name: "hollow_sphere_mesh".to_string(),
    }
}

/// Generate a beam mesh (for cantilever benchmark)
/// Long rectangular block with aspect ratio control
pub fn generate_beam_mesh(
    length: f64,
    width: f64,
    height: f64,
    n_length: usize,
    n_width: usize,
    n_height: usize,
) -> MeshAssembly {
    // Beam oriented along X axis
    generate_block_mesh(length, width, height, n_length, n_width, n_height)
}
