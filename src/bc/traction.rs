use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::Vector3;
use serde::{Serialize, Deserialize};
use std::fmt;
use crate::bc::BoundaryConditionType;
use log::debug;

/// Represents a traction (surface force per unit area) boundary condition
/// 
/// # Mathematical Formulation
/// Traction is a force per unit area applied to a surface. Unlike pressure,
/// traction can have arbitrary direction (not just normal to surface).
/// 
/// The equivalent nodal forces are computed by integrating:
/// F = ∫∫ t * N dA
/// 
/// where:
/// - t is the traction vector (force per unit area)
/// - N are the shape functions
/// - dA is the differential area element
/// 
/// # Traction Types
/// - `Uniform`: Same traction vector everywhere on surface
/// - `Normal`: Traction in surface normal direction (like pressure but can vary)
/// - `Tangential`: Traction tangent to surface (for shear loading)
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum TractionType {
    /// Uniform traction vector in global coordinates
    Uniform(Vector3<f64>),
    /// Normal traction (magnitude, positive = tension outward)
    Normal(f64),
    /// Tangential traction along local xi direction
    TangentialXi(f64),
    /// Tangential traction along local eta direction  
    TangentialEta(f64),
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Traction {
    /// Surface element IDs (4-node quad elements)
    elements: Vec<usize>,
    /// Type of traction to apply
    traction_type: TractionType,
    /// Cached node list
    #[serde(skip, default = "Vec::new")]
    nodes: Vec<usize>,
    /// Precomputed nodal forces
    #[serde(skip, default = "Vec::new")]
    nodal_forces: Vec<(usize, Vector3<f64>)>,
}

impl Traction {
    /// Create a uniform traction BC with a given force per unit area vector
    pub fn new(elements: Vec<usize>, traction: Vector3<f64>) -> Self {
        Traction { 
            elements, 
            traction_type: TractionType::Uniform(traction),
            nodes: Vec::new(),
            nodal_forces: Vec::new(),
        }
    }

    /// Create a normal traction BC (positive = tension outward)
    pub fn new_normal(elements: Vec<usize>, magnitude: f64) -> Self {
        Traction {
            elements,
            traction_type: TractionType::Normal(magnitude),
            nodes: Vec::new(),
            nodal_forces: Vec::new(),
        }
    }

    /// Create a tangential (shear) traction BC along xi direction
    pub fn new_shear_xi(elements: Vec<usize>, magnitude: f64) -> Self {
        Traction {
            elements,
            traction_type: TractionType::TangentialXi(magnitude),
            nodes: Vec::new(),
            nodal_forces: Vec::new(),
        }
    }

    /// Create a tangential (shear) traction BC along eta direction
    pub fn new_shear_eta(elements: Vec<usize>, magnitude: f64) -> Self {
        Traction {
            elements,
            traction_type: TractionType::TangentialEta(magnitude),
            nodes: Vec::new(),
            nodal_forces: Vec::new(),
        }
    }

    /// Get 2x2 Gauss points for surface integration
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

    /// Shape function derivatives
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

    /// Compute surface basis vectors and Jacobian at a point
    /// Returns (normal, tangent_xi, tangent_eta, |J|)
    fn compute_surface_basis(
        node_positions: &[Vector3<f64>; 4],
        xi: f64,
        eta: f64,
    ) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>, f64) {
        let (dn_dxi, dn_deta) = Self::shape_derivatives(xi, eta);
        
        let mut dx_dxi = Vector3::zeros();
        let mut dx_deta = Vector3::zeros();
        
        for i in 0..4 {
            dx_dxi += dn_dxi[i] * node_positions[i];
            dx_deta += dn_deta[i] * node_positions[i];
        }
        
        let normal_unnormalized = dx_dxi.cross(&dx_deta);
        let jacobian = normal_unnormalized.norm();
        
        let (normal, tangent_xi, tangent_eta) = if jacobian > 1e-12 {
            let n = normal_unnormalized / jacobian;
            let t_xi = dx_dxi.normalize();
            let t_eta = n.cross(&t_xi); // Ensure orthogonal
            (n, t_xi, t_eta)
        } else {
            (Vector3::zeros(), Vector3::zeros(), Vector3::zeros())
        };
        
        (normal, tangent_xi, tangent_eta, jacobian)
    }

    /// Compute traction vector at a surface point
    fn compute_traction_at_point(
        &self,
        normal: &Vector3<f64>,
        tangent_xi: &Vector3<f64>,
        tangent_eta: &Vector3<f64>,
    ) -> Vector3<f64> {
        match &self.traction_type {
            TractionType::Uniform(t) => *t,
            TractionType::Normal(mag) => *mag * normal,
            TractionType::TangentialXi(mag) => *mag * tangent_xi,
            TractionType::TangentialEta(mag) => *mag * tangent_eta,
        }
    }
}

#[typetag::serde]
impl BoundaryCondition for Traction {
    fn initalize(&mut self, simulation: &Simulation) {
        debug!("Initializing traction condition with {} elements", self.elements.len());
        
        use std::collections::{HashMap, HashSet};
        let mut node_set = HashSet::new();
        let mut nodal_forces_map: HashMap<usize, Vector3<f64>> = HashMap::new();
        
        for &elem_id in &self.elements {
            let element = match simulation.get_element(elem_id) {
                Some(e) => e,
                None => {
                    debug!("Warning: Element {} not found for traction BC", elem_id);
                    continue;
                }
            };
            
            let connectivity = element.get_connectivity();
            if connectivity.len() != 4 {
                debug!("Warning: Traction BC requires 4-node surface elements, got {} nodes", 
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
            
            // Integrate traction over surface
            let gauss_points = Self::get_surface_gauss_points();
            
            for (xi, eta, weight) in gauss_points {
                let n = Self::shape_functions(xi, eta);
                let (normal, tangent_xi, tangent_eta, jacobian) = 
                    Self::compute_surface_basis(&node_positions, xi, eta);
                
                let traction = self.compute_traction_at_point(&normal, &tangent_xi, &tangent_eta);
                
                // Distribute to nodes: F_i = ∫ N_i * t * |J| dξdη
                for (i, &node_id) in connectivity.iter().enumerate() {
                    let force_contribution = n[i] * traction * jacobian * weight;
                    *nodal_forces_map.entry(node_id).or_insert(Vector3::zeros()) += force_contribution;
                }
            }
        }
        
        self.nodes = node_set.into_iter().collect();
        self.nodal_forces = nodal_forces_map.into_iter().collect();
        
        debug!("Traction BC initialized: {} nodes", self.nodes.len());
    }

    fn apply(&mut self, simulation: &mut Simulation) {
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
        BoundaryConditionType::Traction
    }
}

impl fmt::Display for Traction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let type_str = match &self.traction_type {
            TractionType::Uniform(t) => format!("Uniform({:.3e}, {:.3e}, {:.3e})", t.x, t.y, t.z),
            TractionType::Normal(m) => format!("Normal({:.3e})", m),
            TractionType::TangentialXi(m) => format!("TangentialXi({:.3e})", m),
            TractionType::TangentialEta(m) => format!("TangentialEta({:.3e})", m),
        };
        write!(f, "Traction: {} elements, {} nodes, type={}", 
               self.elements.len(), self.nodes.len(), type_str)
    }
}
