//! Built-in example simulations for demonstration and learning

use rust_fea::mesh::{MeshAssembly, MeshElement, MeshNode, NodeGroup, ElementGroup};
use crate::state::{
    BoundaryConditionConfig, FixedBcConfig, LoadBcConfig, TorqueBcConfig, 
    ContactBcConfig, SolverType,
};

/// Example types available to load
#[derive(Clone, Copy, PartialEq)]
pub enum ExampleType {
    CantileverBeam,
    TorqueShaft,
    ContactBlocks,
}

impl ExampleType {
    pub fn name(&self) -> &'static str {
        match self {
            ExampleType::CantileverBeam => "Cantilever Beam",
            ExampleType::TorqueShaft => "Torque Shaft",
            ExampleType::ContactBlocks => "Contact Blocks",
        }
    }
    
    pub fn description(&self) -> &'static str {
        match self {
            ExampleType::CantileverBeam => 
                "A beam fixed at one end with a load applied at the other end. Classic bending problem.",
            ExampleType::TorqueShaft => 
                "A cylindrical shaft with one end fixed and torque applied to the other end. Torsion problem.",
            ExampleType::ContactBlocks => 
                "Two blocks stacked with contact between them. Demonstrates contact mechanics.",
        }
    }
}

/// Mesh resolution presets
#[derive(Clone, Copy, PartialEq, Default)]
pub enum MeshResolution {
    Coarse,
    #[default]
    Medium,
    Fine,
    VeryFine,
    Custom,
}

impl MeshResolution {
    pub fn name(&self) -> &'static str {
        match self {
            MeshResolution::Coarse => "Coarse",
            MeshResolution::Medium => "Medium",
            MeshResolution::Fine => "Fine",
            MeshResolution::VeryFine => "Very Fine",
            MeshResolution::Custom => "Custom",
        }
    }
    
    /// Get multiplier for base mesh divisions
    pub fn multiplier(&self) -> usize {
        match self {
            MeshResolution::Coarse => 1,
            MeshResolution::Medium => 2,
            MeshResolution::Fine => 4,
            MeshResolution::VeryFine => 8,
            MeshResolution::Custom => 2,
        }
    }
}

/// Configuration for example generation
#[derive(Clone)]
pub struct ExampleConfig {
    pub example_type: ExampleType,
    pub resolution: MeshResolution,
    /// Custom mesh parameters (used when resolution is Custom)
    pub custom_params: ExampleMeshParams,
    /// Load magnitude
    pub load_magnitude: f64,
    /// Geometry scale factor
    pub scale: f64,
}

impl Default for ExampleConfig {
    fn default() -> Self {
        Self {
            example_type: ExampleType::CantileverBeam,
            resolution: MeshResolution::Medium,
            custom_params: ExampleMeshParams::default(),
            load_magnitude: 10000.0,
            scale: 1.0,
        }
    }
}

/// Custom mesh parameters for each example type
#[derive(Clone)]
pub struct ExampleMeshParams {
    // Cantilever beam
    pub beam_length: f64,
    pub beam_height: f64,
    pub beam_width: f64,
    pub beam_nx: usize,
    pub beam_ny: usize,
    pub beam_nz: usize,
    
    // Torque shaft
    pub shaft_radius: f64,
    pub shaft_length: f64,
    pub shaft_n_radial: usize,
    pub shaft_n_height: usize,
    pub shaft_n_layers: usize,
    
    // Contact blocks
    pub block_size: f64,
    pub block_gap: f64,
    pub block_divisions: usize,
}

impl Default for ExampleMeshParams {
    fn default() -> Self {
        Self {
            // Cantilever beam defaults
            beam_length: 10.0,
            beam_height: 1.0,
            beam_width: 1.0,
            beam_nx: 20,
            beam_ny: 4,
            beam_nz: 4,
            
            // Torque shaft defaults
            shaft_radius: 0.5,
            shaft_length: 5.0,
            shaft_n_radial: 12,
            shaft_n_height: 15,
            shaft_n_layers: 3,
            
            // Contact blocks defaults
            block_size: 1.0,
            block_gap: 0.01,
            block_divisions: 5,
        }
    }
}

impl ExampleMeshParams {
    /// Get params scaled by resolution multiplier
    pub fn with_resolution(&self, resolution: MeshResolution) -> Self {
        let m = resolution.multiplier();
        Self {
            beam_nx: (self.beam_nx / 2) * m,
            beam_ny: (self.beam_ny / 2).max(1) * m,
            beam_nz: (self.beam_nz / 2).max(1) * m,
            
            shaft_n_radial: (self.shaft_n_radial / 2).max(4) * m,
            shaft_n_height: (self.shaft_n_height / 2) * m,
            shaft_n_layers: (self.shaft_n_layers).max(1).min(m * 2),
            
            block_divisions: (self.block_divisions / 2).max(1) * m,
            
            ..*self
        }
    }
}

/// Generated example with mesh and boundary conditions
pub struct Example {
    pub name: String,
    pub mesh: MeshAssembly,
    pub boundary_conditions: Vec<BoundaryConditionConfig>,
    pub solver_type: SolverType,
    pub description: String,
}

/// Load a built-in example with default settings
pub fn load_example(example_type: ExampleType) -> Example {
    let config = ExampleConfig {
        example_type,
        resolution: MeshResolution::Medium,
        ..Default::default()
    };
    load_example_with_config(&config)
}

/// Load a built-in example with custom configuration
pub fn load_example_with_config(config: &ExampleConfig) -> Example {
    let params = if config.resolution == MeshResolution::Custom {
        config.custom_params.clone()
    } else {
        ExampleMeshParams::default().with_resolution(config.resolution)
    };
    
    match config.example_type {
        ExampleType::CantileverBeam => create_cantilever_beam(&params, config.load_magnitude, config.scale),
        ExampleType::TorqueShaft => create_torque_shaft(&params, config.load_magnitude, config.scale),
        ExampleType::ContactBlocks => create_contact_blocks(&params, config.load_magnitude, config.scale),
    }
}

/// Create a cantilever beam example
fn create_cantilever_beam(params: &ExampleMeshParams, load: f64, scale: f64) -> Example {
    let mut mesh = MeshAssembly::empty();
    mesh.name = "Cantilever Beam".to_string();
    
    // Scaled dimensions
    let length = params.beam_length * scale;
    let height = params.beam_height * scale;
    let width = params.beam_width * scale;
    
    let nx = params.beam_nx;
    let ny = params.beam_ny;
    let nz = params.beam_nz;
    
    // Create nodes
    let mut node_id = 0;
    let mut node_ids: Vec<Vec<Vec<usize>>> = Vec::new();
    
    for iz in 0..=nz {
        let mut plane = Vec::new();
        for iy in 0..=ny {
            let mut row = Vec::new();
            for ix in 0..=nx {
                let x = (ix as f64 / nx as f64) * length;
                let y = (iy as f64 / ny as f64) * height;
                let z = (iz as f64 / nz as f64) * width;
                
                mesh.nodes.insert(node_id, MeshNode { 
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                row.push(node_id);
                node_id += 1;
            }
            plane.push(row);
        }
        node_ids.push(plane);
    }
    
    // Create brick elements
    let mut elem_id = 0;
    let mut element_ids = Vec::new();
    
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let connectivity = vec![
                    node_ids[iz][iy][ix],
                    node_ids[iz][iy][ix + 1],
                    node_ids[iz][iy + 1][ix + 1],
                    node_ids[iz][iy + 1][ix],
                    node_ids[iz + 1][iy][ix],
                    node_ids[iz + 1][iy][ix + 1],
                    node_ids[iz + 1][iy + 1][ix + 1],
                    node_ids[iz + 1][iy + 1][ix],
                ];
                
                mesh.elements.insert(elem_id, MeshElement {
                    el_type: "C3D8".to_string(),
                    connectivity,
                    name: format!("Element_{}", elem_id),
                    id: elem_id,
                });
                element_ids.push(elem_id);
                elem_id += 1;
            }
        }
    }
    
    mesh.element_groups.insert("all_elements".to_string(), ElementGroup {
        elements: element_ids,
        name: "all_elements".to_string(),
        el_type: "C3D8".to_string(),
    });
    
    // Fixed end node group (x = 0)
    let mut fixed_nodes = Vec::new();
    for iz in 0..=nz {
        for iy in 0..=ny {
            fixed_nodes.push(node_ids[iz][iy][0]);
        }
    }
    mesh.node_groups.insert("fixed_end".to_string(), NodeGroup {
        nodes: fixed_nodes,
        name: "fixed_end".to_string(),
    });
    
    // Load end node group (x = length)
    let mut load_nodes = Vec::new();
    for iz in 0..=nz {
        for iy in 0..=ny {
            load_nodes.push(node_ids[iz][iy][nx]);
        }
    }
    mesh.node_groups.insert("load_end".to_string(), NodeGroup {
        nodes: load_nodes.clone(),
        name: "load_end".to_string(),
    });
    
    // Boundary conditions
    let boundary_conditions = vec![
        BoundaryConditionConfig::Fixed(FixedBcConfig {
            name: "Fixed Support".to_string(),
            node_group: "fixed_end".to_string(),
            constrain_x: Some(0.0),
            constrain_y: Some(0.0),
            constrain_z: Some(0.0),
        }),
        BoundaryConditionConfig::Load(LoadBcConfig {
            name: "Tip Load".to_string(),
            node_group: "load_end".to_string(),
            force_x: 0.0,
            force_y: -load / load_nodes.len() as f64,
            force_z: 0.0,
        }),
    ];
    
    let total_nodes = mesh.nodes.len();
    let total_elements = mesh.elements.len();
    
    Example {
        name: "Cantilever Beam".to_string(),
        mesh,
        boundary_conditions,
        solver_type: SolverType::Direct,
        description: format!(
            "Cantilever beam ({:.1}×{:.1}×{:.1} m) fixed at one end with {:.0} N tip load.\n\
             Mesh: {} nodes, {} elements ({} divisions along length)",
            length, height, width, load, total_nodes, total_elements, nx
        ),
    }
}

/// Create a torque shaft example
fn create_torque_shaft(params: &ExampleMeshParams, load: f64, scale: f64) -> Example {
    let mut mesh = MeshAssembly::empty();
    mesh.name = "Torque Shaft".to_string();
    
    let radius = params.shaft_radius * scale;
    let length = params.shaft_length * scale;
    let n_radial = params.shaft_n_radial;
    let n_height = params.shaft_n_height;
    let n_layers = params.shaft_n_layers;
    
    // Create nodes for solid cylinder
    let mut node_id = 0;
    let mut node_map: std::collections::HashMap<(usize, usize, usize), usize> = std::collections::HashMap::new();
    
    // Center line nodes
    for ih in 0..=n_height {
        let z = (ih as f64 / n_height as f64) * length;
        mesh.nodes.insert(node_id, MeshNode {
            coordinates: vec![0.0, 0.0, z],
            id: node_id,
        });
        node_map.insert((0, 0, ih), node_id);
        node_id += 1;
    }
    
    // Radial layers
    for layer in 1..=n_layers {
        let r = radius * (layer as f64 / n_layers as f64);
        for ih in 0..=n_height {
            let z = (ih as f64 / n_height as f64) * length;
            for ir in 0..n_radial {
                let theta = 2.0 * std::f64::consts::PI * (ir as f64 / n_radial as f64);
                let x = r * theta.cos();
                let y = r * theta.sin();
                
                mesh.nodes.insert(node_id, MeshNode {
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                node_map.insert((layer, ir, ih), node_id);
                node_id += 1;
            }
        }
    }
    
    // Create elements
    let mut elem_id = 0;
    let mut element_ids = Vec::new();
    
    for ih in 0..n_height {
        // Inner wedge elements
        for ir in 0..n_radial {
            let ir_next = (ir + 1) % n_radial;
            
            let center_bottom = node_map[&(0, 0, ih)];
            let center_top = node_map[&(0, 0, ih + 1)];
            let n0 = node_map[&(1, ir, ih)];
            let n1 = node_map[&(1, ir_next, ih)];
            let n4 = node_map[&(1, ir, ih + 1)];
            let n5 = node_map[&(1, ir_next, ih + 1)];
            
            let connectivity = vec![
                center_bottom, n0, n1, center_bottom,
                center_top, n4, n5, center_top,
            ];
            
            mesh.elements.insert(elem_id, MeshElement {
                el_type: "C3D8".to_string(),
                connectivity,
                name: format!("Element_{}", elem_id),
                id: elem_id,
            });
            element_ids.push(elem_id);
            elem_id += 1;
        }
        
        // Outer hex elements
        for layer in 1..n_layers {
            for ir in 0..n_radial {
                let ir_next = (ir + 1) % n_radial;
                
                let connectivity = vec![
                    node_map[&(layer, ir, ih)],
                    node_map[&(layer + 1, ir, ih)],
                    node_map[&(layer + 1, ir_next, ih)],
                    node_map[&(layer, ir_next, ih)],
                    node_map[&(layer, ir, ih + 1)],
                    node_map[&(layer + 1, ir, ih + 1)],
                    node_map[&(layer + 1, ir_next, ih + 1)],
                    node_map[&(layer, ir_next, ih + 1)],
                ];
                
                mesh.elements.insert(elem_id, MeshElement {
                    el_type: "C3D8".to_string(),
                    connectivity,
                    name: format!("Element_{}", elem_id),
                    id: elem_id,
                });
                element_ids.push(elem_id);
                elem_id += 1;
            }
        }
    }
    
    mesh.element_groups.insert("all_elements".to_string(), ElementGroup {
        elements: element_ids,
        name: "all_elements".to_string(),
        el_type: "C3D8".to_string(),
    });
    
    // Node groups
    let mut fixed_nodes = Vec::new();
    fixed_nodes.push(node_map[&(0, 0, 0)]);
    for layer in 1..=n_layers {
        for ir in 0..n_radial {
            fixed_nodes.push(node_map[&(layer, ir, 0)]);
        }
    }
    mesh.node_groups.insert("fixed_end".to_string(), NodeGroup {
        nodes: fixed_nodes,
        name: "fixed_end".to_string(),
    });
    
    let mut torque_nodes = Vec::new();
    torque_nodes.push(node_map[&(0, 0, n_height)]);
    for layer in 1..=n_layers {
        for ir in 0..n_radial {
            torque_nodes.push(node_map[&(layer, ir, n_height)]);
        }
    }
    mesh.node_groups.insert("torque_end".to_string(), NodeGroup {
        nodes: torque_nodes,
        name: "torque_end".to_string(),
    });
    
    let boundary_conditions = vec![
        BoundaryConditionConfig::Fixed(FixedBcConfig {
            name: "Fixed Support".to_string(),
            node_group: "fixed_end".to_string(),
            constrain_x: Some(0.0),
            constrain_y: Some(0.0),
            constrain_z: Some(0.0),
        }),
        BoundaryConditionConfig::Torque(TorqueBcConfig {
            name: "Applied Torque".to_string(),
            node_group: "torque_end".to_string(),
            axis_point: [0.0, 0.0, length],
            axis_direction: [0.0, 0.0, 1.0],
            magnitude: load,
        }),
    ];
    
    let total_nodes = mesh.nodes.len();
    let total_elements = mesh.elements.len();
    
    Example {
        name: "Torque Shaft".to_string(),
        mesh,
        boundary_conditions,
        solver_type: SolverType::Direct,
        description: format!(
            "Cylindrical shaft (R={:.2}m, L={:.1}m) with {:.0} N·m torque.\n\
             Mesh: {} nodes, {} elements ({} radial, {} axial divisions)",
            radius, length, load, total_nodes, total_elements, n_radial, n_height
        ),
    }
}

/// Create contact blocks example
fn create_contact_blocks(params: &ExampleMeshParams, load: f64, scale: f64) -> Example {
    let mut mesh = MeshAssembly::empty();
    mesh.name = "Contact Blocks".to_string();
    
    let size = params.block_size * scale;
    let gap = params.block_gap * scale;
    let n = params.block_divisions;
    
    // Create bottom block
    let mut node_id = 0;
    let mut bottom_node_ids: Vec<Vec<Vec<usize>>> = Vec::new();
    
    for iz in 0..=n {
        let mut plane = Vec::new();
        for iy in 0..=n {
            let mut row = Vec::new();
            for ix in 0..=n {
                let x = (ix as f64 / n as f64) * size;
                let y = (iy as f64 / n as f64) * size;
                let z = (iz as f64 / n as f64) * size;
                
                mesh.nodes.insert(node_id, MeshNode { 
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                row.push(node_id);
                node_id += 1;
            }
            plane.push(row);
        }
        bottom_node_ids.push(plane);
    }
    
    // Create top block
    let mut top_node_ids: Vec<Vec<Vec<usize>>> = Vec::new();
    let top_offset = size + gap;
    
    for iz in 0..=n {
        let mut plane = Vec::new();
        for iy in 0..=n {
            let mut row = Vec::new();
            for ix in 0..=n {
                let x = (ix as f64 / n as f64) * size;
                let y = (iy as f64 / n as f64) * size;
                let z = top_offset + (iz as f64 / n as f64) * size;
                
                mesh.nodes.insert(node_id, MeshNode { 
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                row.push(node_id);
                node_id += 1;
            }
            plane.push(row);
        }
        top_node_ids.push(plane);
    }
    
    // Create elements
    let mut elem_id = 0;
    let mut element_ids = Vec::new();
    
    // Bottom block elements
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let connectivity = vec![
                    bottom_node_ids[iz][iy][ix],
                    bottom_node_ids[iz][iy][ix + 1],
                    bottom_node_ids[iz][iy + 1][ix + 1],
                    bottom_node_ids[iz][iy + 1][ix],
                    bottom_node_ids[iz + 1][iy][ix],
                    bottom_node_ids[iz + 1][iy][ix + 1],
                    bottom_node_ids[iz + 1][iy + 1][ix + 1],
                    bottom_node_ids[iz + 1][iy + 1][ix],
                ];
                
                mesh.elements.insert(elem_id, MeshElement {
                    el_type: "C3D8".to_string(),
                    connectivity,
                    name: format!("Bottom_{}", elem_id),
                    id: elem_id,
                });
                element_ids.push(elem_id);
                elem_id += 1;
            }
        }
    }
    
    // Top block elements
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let connectivity = vec![
                    top_node_ids[iz][iy][ix],
                    top_node_ids[iz][iy][ix + 1],
                    top_node_ids[iz][iy + 1][ix + 1],
                    top_node_ids[iz][iy + 1][ix],
                    top_node_ids[iz + 1][iy][ix],
                    top_node_ids[iz + 1][iy][ix + 1],
                    top_node_ids[iz + 1][iy + 1][ix + 1],
                    top_node_ids[iz + 1][iy + 1][ix],
                ];
                
                mesh.elements.insert(elem_id, MeshElement {
                    el_type: "C3D8".to_string(),
                    connectivity,
                    name: format!("Top_{}", elem_id),
                    id: elem_id,
                });
                element_ids.push(elem_id);
                elem_id += 1;
            }
        }
    }
    
    mesh.element_groups.insert("all_elements".to_string(), ElementGroup {
        elements: element_ids,
        name: "all_elements".to_string(),
        el_type: "C3D8".to_string(),
    });
    
    // Node groups
    let mut fixed_nodes = Vec::new();
    for iy in 0..=n {
        for ix in 0..=n {
            fixed_nodes.push(bottom_node_ids[0][iy][ix]);
        }
    }
    mesh.node_groups.insert("fixed_base".to_string(), NodeGroup {
        nodes: fixed_nodes,
        name: "fixed_base".to_string(),
    });
    
    let mut bottom_contact = Vec::new();
    for iy in 0..=n {
        for ix in 0..=n {
            bottom_contact.push(bottom_node_ids[n][iy][ix]);
        }
    }
    mesh.node_groups.insert("bottom_contact".to_string(), NodeGroup {
        nodes: bottom_contact,
        name: "bottom_contact".to_string(),
    });
    
    let mut top_contact = Vec::new();
    for iy in 0..=n {
        for ix in 0..=n {
            top_contact.push(top_node_ids[0][iy][ix]);
        }
    }
    mesh.node_groups.insert("top_contact".to_string(), NodeGroup {
        nodes: top_contact,
        name: "top_contact".to_string(),
    });
    
    let mut load_nodes = Vec::new();
    for iy in 0..=n {
        for ix in 0..=n {
            load_nodes.push(top_node_ids[n][iy][ix]);
        }
    }
    mesh.node_groups.insert("load_surface".to_string(), NodeGroup {
        nodes: load_nodes.clone(),
        name: "load_surface".to_string(),
    });
    
    let boundary_conditions = vec![
        BoundaryConditionConfig::Fixed(FixedBcConfig {
            name: "Fixed Base".to_string(),
            node_group: "fixed_base".to_string(),
            constrain_x: Some(0.0),
            constrain_y: Some(0.0),
            constrain_z: Some(0.0),
        }),
        BoundaryConditionConfig::Contact(ContactBcConfig {
            name: "Block Contact".to_string(),
            primary_surface: "bottom_contact".to_string(),
            secondary_surface: "top_contact".to_string(),
        }),
        BoundaryConditionConfig::Load(LoadBcConfig {
            name: "Compressive Load".to_string(),
            node_group: "load_surface".to_string(),
            force_x: 0.0,
            force_y: 0.0,
            force_z: -load / load_nodes.len() as f64,
        }),
    ];
    
    let total_nodes = mesh.nodes.len();
    let total_elements = mesh.elements.len();
    
    Example {
        name: "Contact Blocks".to_string(),
        mesh,
        boundary_conditions,
        solver_type: SolverType::Explicit,
        description: format!(
            "Two {:.1}m³ blocks with {:.3}m gap, {:.0} N compressive load.\n\
             Mesh: {} nodes, {} elements ({} divisions per side)",
            size, gap, load, total_nodes, total_elements, n
        ),
    }
}
