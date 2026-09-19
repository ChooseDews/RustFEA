use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::Vector3;
use serde::{Serialize, Deserialize};
use std::fmt;
use crate::bc::BoundaryConditionType;
use log::debug;

/// Represents a pressure boundary condition applied to surface elements
/// 
/// # Mathematical Formulation
/// Pressure is a normal traction applied uniformly over a surface.
/// The equivalent nodal forces are computed by integrating the pressure
/// over the surface: F = ∫∫ p * n * N dA
/// 
/// where:
/// - p is the pressure magnitude (positive = compression into surface)
/// - n is the outward surface normal
/// - N are the shape functions
/// - dA is the differential area element
/// 
/// For a 4-node surface element, we use 2x2 Gauss quadrature.
#[derive(Serialize, Deserialize, Debug)]
pub struct PressureCondition {
    /// Surface element IDs (4-node quad elements)
    elements: Vec<usize>,
    /// Pressure magnitude (positive = compression, negative = tension)
    pressure: f64,
    /// Cached node list for get_nodes()
    #[serde(skip, default = "Vec::new")]
    nodes: Vec<usize>,
    /// Precomputed nodal forces (node_id -> force vector)
    #[serde(skip, default = "Vec::new")]
    nodal_forces: Vec<(usize, Vector3<f64>)>,
}

impl PressureCondition {
    pub fn new(elements: Vec<usize>, pressure: f64) -> Self {
        PressureCondition { 
            elements, 
            pressure, 
            nodes: Vec::new(),
            nodal_forces: Vec::new(),
        }
    }

    /// Get 2x2 Gauss points for surface integration
    /// Returns (xi, eta, weight)
    fn get_surface_gauss_points() -> [(f64, f64, f64); 4] {
        let a = 1.0 / 3.0_f64.sqrt();
        [
            (-a, -a, 1.0),
            ( a, -a, 1.0),
            ( a,  a, 1.0),
            (-a,  a, 1.0),
        ]
    }

    /// Shape functions for 4-node quad element
    fn shape_functions(xi: f64, eta: f64) -> [f64; 4] {
        [
            0.25 * (1.0 - xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 + eta),
            0.25 * (1.0 - xi) * (1.0 + eta),
        ]
    }

    /// Shape function derivatives w.r.t. xi and eta
    /// Returns dN/dxi and dN/deta for each node
    fn shape_derivatives(xi: f64, eta: f64) -> ([f64; 4], [f64; 4]) {
        let dn_dxi = [
            -0.25 * (1.0 - eta),
             0.25 * (1.0 - eta),
             0.25 * (1.0 + eta),
            -0.25 * (1.0 + eta),
        ];
        let dn_deta = [
            -0.25 * (1.0 - xi),
            -0.25 * (1.0 + xi),
             0.25 * (1.0 + xi),
             0.25 * (1.0 - xi),
        ];
        (dn_dxi, dn_deta)
    }

    /// Compute surface normal and Jacobian determinant at a point
    /// Returns (normal, |J|) where normal points outward
    fn compute_surface_jacobian(
        node_positions: &[Vector3<f64>; 4],
        xi: f64,
        eta: f64,
    ) -> (Vector3<f64>, f64) {
        let (dn_dxi, dn_deta) = Self::shape_derivatives(xi, eta);
        
        // Compute tangent vectors
        let mut dx_dxi = Vector3::zeros();
        let mut dx_deta = Vector3::zeros();
        
        for i in 0..4 {
            dx_dxi += dn_dxi[i] * node_positions[i];
            dx_deta += dn_deta[i] * node_positions[i];
        }
        
        // Normal = tangent1 x tangent2 (right-hand rule gives outward normal
        // if nodes are ordered counter-clockwise when viewed from outside)
        let normal_unnormalized = dx_dxi.cross(&dx_deta);
        let jacobian = normal_unnormalized.norm();
        
        let normal = if jacobian > 1e-12 {
            normal_unnormalized / jacobian
        } else {
            Vector3::zeros()
        };
        
        (normal, jacobian)
    }
}

#[typetag::serde]
impl BoundaryCondition for PressureCondition {
    fn initalize(&mut self, simulation: &Simulation) {
        debug!("Initializing pressure condition with {} elements", self.elements.len());
        
        use std::collections::{HashMap, HashSet};
        let mut node_set = HashSet::new();
        let mut nodal_forces_map: HashMap<usize, Vector3<f64>> = HashMap::new();
        
        for &elem_id in &self.elements {
            let element = match simulation.get_element(elem_id) {
                Some(e) => e,
                None => {
                    debug!("Warning: Element {} not found for pressure BC", elem_id);
                    continue;
                }
            };
            
            let connectivity = element.get_connectivity();
            if connectivity.len() != 4 {
                debug!("Warning: Pressure BC requires 4-node surface elements, got {} nodes", 
                       connectivity.len());
                continue;
            }
            
            // Get node positions
            let mut node_positions = [Vector3::zeros(); 4];
            for (i, &node_id) in connectivity.iter().enumerate() {
                if let Some(node) = simulation.get_node(node_id) {
                    node_positions[i] = node.position;
                    node_set.insert(node_id);
                }
            }
            
            // Integrate pressure over surface using Gauss quadrature
            let gauss_points = Self::get_surface_gauss_points();
            
            for (xi, eta, weight) in gauss_points {
                let n = Self::shape_functions(xi, eta);
                let (normal, jacobian) = Self::compute_surface_jacobian(&node_positions, xi, eta);
                
                // Force per unit area: -pressure * normal (negative because pressure acts inward)
                // For positive pressure (compression), force should push INTO the surface
                let traction = -self.pressure * normal;
                
                // Distribute to nodes: F_i = ∫ N_i * t * |J| dξdη
                for (i, &node_id) in connectivity.iter().enumerate() {
                    let force_contribution = n[i] * traction * jacobian * weight;
                    *nodal_forces_map.entry(node_id).or_insert(Vector3::zeros()) += force_contribution;
                }
            }
        }
        
        // Store results
        self.nodes = node_set.into_iter().collect();
        self.nodal_forces = nodal_forces_map.into_iter().collect();
        
        debug!("Pressure BC initialized: {} nodes, total force magnitude: {:.6e}",
               self.nodes.len(),
               self.nodal_forces.iter().map(|(_, f)| f.norm()).sum::<f64>());
    }

    fn apply(&mut self, simulation: &mut Simulation) {
        // Apply precomputed forces
        for (node_id, force) in &self.nodal_forces {
            for i in 0..3 {
                let global_index = simulation.get_global_index(*node_id, i);
                simulation.load_vector[global_index] += force[i];
            }
        }
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }

    fn type_name(&self) -> BoundaryConditionType {
        BoundaryConditionType::Pressure
    }
}

impl fmt::Display for PressureCondition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PressureCondition: {} elements, {} nodes, P={:.3e}", 
               self.elements.len(), self.nodes.len(), self.pressure)
    }
}
