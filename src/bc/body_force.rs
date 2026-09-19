use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::Vector3;
use serde::{Serialize, Deserialize};
use std::fmt;
use crate::bc::BoundaryConditionType;
use log::debug;

/// Represents a body force boundary condition (distributed force per unit volume)
/// 
/// # Mathematical Formulation
/// Body forces act throughout the volume of the material, not just on surfaces.
/// Common examples:
/// - Gravity: b = ρg (density × gravitational acceleration)
/// - Centrifugal: b = ρω²r (density × angular velocity² × radius)
/// 
/// The equivalent nodal forces are computed by integrating over element volumes:
/// F = ∫∫∫ b * N dV
/// 
/// For element-wise constant body forces with lumped mass, this simplifies to:
/// F_i = (m_i / ρ) * b  where m_i is the lumped mass at node i
/// 
/// # Body Force Types
/// - `Gravity`: Uniform acceleration field (e.g., [0, 0, -9.81])
/// - `Centrifugal`: Rotation about an axis (point, axis, angular velocity)
/// - `Uniform`: Arbitrary uniform body force per unit volume
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum BodyForceType {
    /// Gravitational acceleration vector [m/s²]
    Gravity(Vector3<f64>),
    /// Centrifugal force: (axis_point, axis_direction, angular_velocity [rad/s])
    Centrifugal {
        axis_point: Vector3<f64>,
        axis_direction: Vector3<f64>,
        angular_velocity: f64,
    },
    /// Uniform body force per unit volume [N/m³]
    Uniform(Vector3<f64>),
}

#[derive(Serialize, Deserialize, Debug)]
pub struct BodyForce {
    /// Type of body force
    force_type: BodyForceType,
    /// Element IDs to apply body force (empty = all active elements)
    elements: Vec<usize>,
    /// Cached affected nodes
    #[serde(skip, default = "Vec::new")]
    nodes: Vec<usize>,
}

impl BodyForce {
    /// Create a gravity body force with standard Earth gravity in -Z direction
    pub fn gravity() -> Self {
        BodyForce {
            force_type: BodyForceType::Gravity(Vector3::new(0.0, 0.0, -9.81)),
            elements: Vec::new(),
            nodes: Vec::new(),
        }
    }

    /// Create a gravity body force with custom acceleration vector
    pub fn gravity_custom(acceleration: Vector3<f64>) -> Self {
        BodyForce {
            force_type: BodyForceType::Gravity(acceleration),
            elements: Vec::new(),
            nodes: Vec::new(),
        }
    }

    /// Create a gravity body force from vector components
    pub fn gravity_from_vec(gx: f64, gy: f64, gz: f64) -> Self {
        BodyForce {
            force_type: BodyForceType::Gravity(Vector3::new(gx, gy, gz)),
            elements: Vec::new(),
            nodes: Vec::new(),
        }
    }

    /// Create a centrifugal body force
    /// 
    /// # Arguments
    /// * `axis_point` - A point on the rotation axis
    /// * `axis_direction` - Direction of rotation axis (will be normalized)
    /// * `angular_velocity` - Angular velocity in rad/s
    pub fn centrifugal(
        axis_point: Vector3<f64>,
        axis_direction: Vector3<f64>,
        angular_velocity: f64,
    ) -> Self {
        BodyForce {
            force_type: BodyForceType::Centrifugal {
                axis_point,
                axis_direction: axis_direction.normalize(),
                angular_velocity,
            },
            elements: Vec::new(),
            nodes: Vec::new(),
        }
    }

    /// Create a centrifugal body force from vector components
    pub fn centrifugal_from_vec(
        axis_point: Vec<f64>,
        axis_direction: Vec<f64>,
        angular_velocity: f64,
    ) -> Self {
        Self::centrifugal(
            Vector3::new(axis_point[0], axis_point[1], axis_point[2]),
            Vector3::new(axis_direction[0], axis_direction[1], axis_direction[2]),
            angular_velocity,
        )
    }

    /// Create a uniform body force per unit volume
    pub fn uniform(force_per_volume: Vector3<f64>) -> Self {
        BodyForce {
            force_type: BodyForceType::Uniform(force_per_volume),
            elements: Vec::new(),
            nodes: Vec::new(),
        }
    }

    /// Restrict body force to specific elements
    pub fn with_elements(mut self, elements: Vec<usize>) -> Self {
        self.elements = elements;
        self
    }

    /// Compute body force at a given position
    fn compute_force_at_position(
        &self,
        position: &Vector3<f64>,
        density: f64,
    ) -> Vector3<f64> {
        match &self.force_type {
            BodyForceType::Gravity(g) => density * g,
            BodyForceType::Centrifugal { axis_point, axis_direction, angular_velocity } => {
                // Vector from axis to point
                let to_point = position - axis_point;
                // Project onto axis
                let parallel = to_point.dot(axis_direction) * axis_direction;
                // Perpendicular component (radial direction from axis)
                let radial = to_point - parallel;
                let r = radial.norm();
                
                if r < 1e-12 {
                    return Vector3::zeros();
                }
                
                // Centrifugal force: F = ρω²r (pointing outward from axis)
                let omega_squared = angular_velocity * angular_velocity;
                density * omega_squared * radial
            }
            BodyForceType::Uniform(b) => *b,
        }
    }
}

#[typetag::serde]
impl BoundaryCondition for BodyForce {
    fn initalize(&mut self, simulation: &Simulation) {
        use std::collections::HashSet;
        let mut node_set = HashSet::new();
        
        // Determine which elements to process
        let elem_ids: Vec<usize> = if self.elements.is_empty() {
            simulation.active_elements()
        } else {
            self.elements.clone()
        };
        
        for &elem_id in &elem_ids {
            if let Some(element) = simulation.get_element(elem_id) {
                for &node_id in element.get_connectivity() {
                    node_set.insert(node_id);
                }
            }
        }
        
        self.nodes = node_set.into_iter().collect();
        debug!("Body force BC initialized: {} nodes from {} elements", 
               self.nodes.len(), elem_ids.len());
    }

    fn apply(&mut self, simulation: &mut Simulation) {
        // Determine which elements to process
        let elem_ids: Vec<usize> = if self.elements.is_empty() {
            simulation.active_elements()
        } else {
            self.elements.clone()
        };
        
        // For each element, compute nodal forces from body force
        for &elem_id in &elem_ids {
            let element = match simulation.get_element(elem_id) {
                Some(e) => e,
                None => continue,
            };
            
            let connectivity = element.get_connectivity().clone();
            let density = element.get_material().get_density();
            
            // Get lumped mass - if empty, compute it on the fly
            let lumped_mass = element.get_lumped_mass().clone();
            let use_precomputed_mass = !lumped_mass.is_empty();
            
            if use_precomputed_mass {
                // Use pre-computed lumped mass
                for (local_idx, &node_id) in connectivity.iter().enumerate() {
                    let node_mass = if local_idx < lumped_mass.len() {
                        lumped_mass[local_idx]
                    } else {
                        continue;
                    };
                    
                    let position = match simulation.get_node(node_id) {
                        Some(n) => n.position,
                        None => continue,
                    };
                    
                    let body_force_per_volume = self.compute_force_at_position(&position, density);
                    let nodal_force = body_force_per_volume * (node_mass / density);
                    
                    for i in 0..3 {
                        let global_index = simulation.get_global_index(node_id, i);
                        simulation.load_vector[global_index] += nodal_force[i];
                    }
                }
            } else {
                // Compute mass on the fly using element mass matrix
                // The mass matrix is NxN where N = number of nodes (e.g., 8 for brick)
                // Lumped mass per node = row sum
                let mass_matrix = element.compute_mass(simulation);
                let num_nodes = connectivity.len();
                
                for (local_idx, &node_id) in connectivity.iter().enumerate() {
                    // Lumped mass for this node is the row sum
                    let node_mass = if local_idx < mass_matrix.nrows() {
                        mass_matrix.row(local_idx).sum()
                    } else {
                        continue;
                    };
                    
                    if node_mass < 1e-20 {
                        continue;
                    }
                    
                    let position = match simulation.get_node(node_id) {
                        Some(n) => n.position,
                        None => continue,
                    };
                    
                    let body_force_per_volume = self.compute_force_at_position(&position, density);
                    let nodal_force = body_force_per_volume * (node_mass / density);
                    
                    for i in 0..3 {
                        let global_index = simulation.get_global_index(node_id, i);
                        simulation.load_vector[global_index] += nodal_force[i];
                    }
                }
            }
        }
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }

    fn type_name(&self) -> BoundaryConditionType {
        BoundaryConditionType::BodyForce
    }
}

impl fmt::Display for BodyForce {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let type_str = match &self.force_type {
            BodyForceType::Gravity(g) => format!("Gravity({:.3}, {:.3}, {:.3})", g.x, g.y, g.z),
            BodyForceType::Centrifugal { angular_velocity, .. } => {
                format!("Centrifugal(ω={:.3} rad/s)", angular_velocity)
            }
            BodyForceType::Uniform(b) => format!("Uniform({:.3e}, {:.3e}, {:.3e})", b.x, b.y, b.z),
        };
        let elem_str = if self.elements.is_empty() {
            "all".to_string()
        } else {
            format!("{}", self.elements.len())
        };
        write!(f, "BodyForce: {} elements, type={}", elem_str, type_str)
    }
}
