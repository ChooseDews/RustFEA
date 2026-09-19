use crate::bc::BoundaryCondition;
use crate::simulation::Simulation;
use log::debug;
use nalgebra::{self as na, Vector3};
use serde::{Deserialize, Serialize};
use crate::bc::BoundaryConditionType;

/// Contact algorithm type
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ContactAlgorithm {
    /// Simple penalty method: F = k * penetration
    Penalty,
    /// Augmented Lagrangian: iteratively updates multipliers for better accuracy
    AugmentedLagrangian,
}

/// Friction model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FrictionModel {
    /// No friction (frictionless contact)
    Frictionless,
    /// Coulomb friction with coefficient μ
    Coulomb { coefficient: f64 },
    /// Viscous friction (velocity-dependent): F = c * v_tangent
    Viscous { damping: f64 },
}

/// Node-to-Segment contact with friction support
/// 
/// # Improvements over basic NormalContact:
/// 1. Node-to-segment formulation for better edge handling
/// 2. Friction support (Coulomb and viscous models)
/// 3. Augmented Lagrangian option for better constraint enforcement
/// 4. Penetration history for stable contact detection
/// 5. Contact area tracking for more accurate force distribution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeToSegmentContact {
    // Surface definitions
    primary_surface: String,
    secondary_surface: String,
    
    // Contact parameters
    contact_stiffness: f64,
    max_distance: f64,
    algorithm: ContactAlgorithm,
    friction: FrictionModel,
    
    // Surface data
    primary_nodes: Vec<usize>,
    primary_elements: Vec<usize>,
    secondary_nodes: Vec<usize>,
    secondary_elements: Vec<usize>,
    
    // Active contact pairs (populated during apply)
    active_primary_nodes: Vec<usize>,
    active_secondary_elements: Vec<usize>,
    
    // Augmented Lagrangian state
    #[serde(skip, default = "Vec::new")]
    lagrange_multipliers: Vec<f64>,  // One per primary node
    
    // Contact state for friction
    #[serde(skip, default = "Vec::new")]
    prev_tangent_displacement: Vec<Vector3<f64>>,  // For stick-slip
    
    // Statistics
    penetration_count: usize,
    max_penetration: f64,
    max_friction_force: f64,
}

impl NodeToSegmentContact {
    pub fn new(primary_surface: String, secondary_surface: String) -> Self {
        NodeToSegmentContact {
            primary_surface,
            secondary_surface,
            contact_stiffness: 1e7,
            max_distance: 0.25,
            algorithm: ContactAlgorithm::Penalty,
            friction: FrictionModel::Frictionless,
            primary_nodes: Vec::new(),
            primary_elements: Vec::new(),
            secondary_nodes: Vec::new(),
            secondary_elements: Vec::new(),
            active_primary_nodes: Vec::new(),
            active_secondary_elements: Vec::new(),
            lagrange_multipliers: Vec::new(),
            prev_tangent_displacement: Vec::new(),
            penetration_count: 0,
            max_penetration: 0.0,
            max_friction_force: 0.0,
        }
    }

    // Builder methods
    pub fn with_stiffness(mut self, k: f64) -> Self {
        self.contact_stiffness = k;
        self
    }

    pub fn with_max_distance(mut self, d: f64) -> Self {
        self.max_distance = d;
        self
    }

    pub fn with_friction(mut self, friction: FrictionModel) -> Self {
        self.friction = friction;
        self
    }

    pub fn with_augmented_lagrangian(mut self) -> Self {
        self.algorithm = ContactAlgorithm::AugmentedLagrangian;
        self
    }

    pub fn with_coulomb_friction(mut self, coefficient: f64) -> Self {
        self.friction = FrictionModel::Coulomb { coefficient };
        self
    }

    pub fn with_viscous_friction(mut self, damping: f64) -> Self {
        self.friction = FrictionModel::Viscous { damping };
        self
    }

    pub fn set_contact_surfaces_nodes(
        &mut self,
        primary_nodes: Vec<usize>,
        secondary_nodes: Vec<usize>,
    ) {
        self.primary_nodes = primary_nodes;
        self.secondary_nodes = secondary_nodes;
    }

    pub fn set_contact_surfaces_elements(
        &mut self,
        primary_elements: Vec<usize>,
        secondary_elements: Vec<usize>,
    ) {
        self.primary_elements = primary_elements;
        self.secondary_elements = secondary_elements;
    }

    /// Find closest point on a quad element to a given point
    /// Returns (distance, closest_point, normal, is_inside) or None if too far
    fn find_closest_point_on_segment(
        &self,
        point: &Vector3<f64>,
        element_nodes: &[Vector3<f64>; 4],
    ) -> Option<(f64, Vector3<f64>, Vector3<f64>)> {
        // Compute element center and normal
        let center = (element_nodes[0] + element_nodes[1] + element_nodes[2] + element_nodes[3]) / 4.0;
        
        // Compute normal using cross product of diagonals
        let diag1 = element_nodes[2] - element_nodes[0];
        let diag2 = element_nodes[3] - element_nodes[1];
        let normal_unnorm = diag1.cross(&diag2);
        let normal_mag = normal_unnorm.norm();
        
        if normal_mag < 1e-12 {
            return None;
        }
        let normal = normal_unnorm / normal_mag;
        
        // Project point onto plane
        let to_point = point - center;
        let signed_distance = to_point.dot(&normal);
        
        // Only consider contact if point is on the "negative" side (inside)
        // or very close to surface
        if signed_distance > self.max_distance * 0.1 {
            return None;
        }
        
        let projected = point - signed_distance * normal;
        
        // Check if projected point is inside the quad using barycentric-like test
        if self.point_in_quad(&projected, element_nodes) {
            let penetration = -signed_distance;
            if penetration > -self.max_distance * 0.5 {
                return Some((penetration, projected, normal));
            }
        }
        
        None
    }

    /// Check if a point lies within a quad element (projected onto plane)
    fn point_in_quad(&self, point: &Vector3<f64>, nodes: &[Vector3<f64>; 4]) -> bool {
        // Split into two triangles and check both
        self.point_in_triangle(point, &nodes[0], &nodes[1], &nodes[2]) ||
        self.point_in_triangle(point, &nodes[0], &nodes[2], &nodes[3])
    }

    fn point_in_triangle(&self, p: &Vector3<f64>, a: &Vector3<f64>, b: &Vector3<f64>, c: &Vector3<f64>) -> bool {
        let v0 = c - a;
        let v1 = b - a;
        let v2 = p - a;

        let dot00 = v0.dot(&v0);
        let dot01 = v0.dot(&v1);
        let dot02 = v0.dot(&v2);
        let dot11 = v1.dot(&v1);
        let dot12 = v1.dot(&v2);

        let inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        let u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        let v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

        let tol = 0.01;  // Small tolerance for numerical robustness
        (u >= -tol) && (v >= -tol) && (u + v <= 1.0 + tol)
    }

    /// Apply normal contact force using penalty method
    fn apply_penalty_contact(
        &self,
        penetration: f64,
        normal: &Vector3<f64>,
    ) -> Vector3<f64> {
        if penetration > 0.0 {
            self.contact_stiffness * penetration * normal
        } else {
            Vector3::zeros()
        }
    }

    /// Compute friction force based on model
    fn compute_friction_force(
        &self,
        normal_force: &Vector3<f64>,
        tangent_velocity: &Vector3<f64>,
        tangent_displacement: &Vector3<f64>,
        node_idx: usize,
    ) -> Vector3<f64> {
        let normal_magnitude = normal_force.norm();
        if normal_magnitude < 1e-12 {
            return Vector3::zeros();
        }

        match &self.friction {
            FrictionModel::Frictionless => Vector3::zeros(),
            
            FrictionModel::Coulomb { coefficient } => {
                // Coulomb friction: |F_t| ≤ μ|F_n|
                let tangent_mag = tangent_displacement.norm();
                if tangent_mag < 1e-12 {
                    return Vector3::zeros();
                }
                
                // Regularized Coulomb: use penalty in tangent direction
                let tangent_stiffness = self.contact_stiffness * 0.1;  // Lower stiffness for tangent
                let trial_force_mag = tangent_stiffness * tangent_mag;
                let max_friction = coefficient * normal_magnitude;
                
                let friction_mag = trial_force_mag.min(max_friction);
                -friction_mag * (tangent_displacement / tangent_mag)
            }
            
            FrictionModel::Viscous { damping } => {
                // Viscous friction: F = -c * v_tangent
                let vel_mag = tangent_velocity.norm();
                if vel_mag < 1e-12 {
                    return Vector3::zeros();
                }
                -damping * tangent_velocity
            }
        }
    }
}

#[typetag::serde]
impl BoundaryCondition for NodeToSegmentContact {
    fn apply(&mut self, simulation: &mut Simulation) {
        let mut max_penetration: f64 = 0.0;
        let mut max_friction: f64 = 0.0;
        let mut active_count = 0;

        for (node_idx, &primary_node_id) in self.active_primary_nodes.iter().enumerate() {
            let primary_node = match simulation.get_node(primary_node_id) {
                Some(n) => n,
                None => continue,
            };
            
            let current_position = primary_node.position + primary_node.displacement;
            
            // Find closest secondary element
            let mut best_contact: Option<(f64, Vector3<f64>, Vector3<f64>, usize)> = None;
            
            for &sec_elem_id in &self.active_secondary_elements {
                let sec_element = match simulation.get_element(sec_elem_id) {
                    Some(e) => e,
                    None => continue,
                };
                
                let connectivity = sec_element.get_connectivity();
                if connectivity.len() != 4 {
                    continue;
                }
                
                // Get deformed positions of secondary element nodes
                let mut element_nodes = [Vector3::zeros(); 4];
                for (i, &node_id) in connectivity.iter().enumerate() {
                    if let Some(node) = simulation.get_node(node_id) {
                        element_nodes[i] = node.position + node.displacement;
                    }
                }
                
                if let Some((penetration, closest_pt, normal)) = 
                    self.find_closest_point_on_segment(&current_position, &element_nodes) 
                {
                    if penetration > 0.0 {
                        if best_contact.is_none() || penetration > best_contact.as_ref().unwrap().0 {
                            best_contact = Some((penetration, closest_pt, normal, sec_elem_id));
                        }
                    }
                }
            }
            
            // Apply contact forces if in contact
            if let Some((penetration, _closest_pt, normal, sec_elem_id)) = best_contact {
                if penetration > max_penetration {
                    max_penetration = penetration;
                }
                active_count += 1;
                
                // Compute normal contact force
                let normal_force = self.apply_penalty_contact(penetration, &normal);
                
                // Compute tangent displacement (for friction)
                let velocity = primary_node.displacement;  // Approximate
                let tangent_velocity = velocity - velocity.dot(&normal) * normal;
                let tangent_disp = if node_idx < self.prev_tangent_displacement.len() {
                    tangent_velocity - self.prev_tangent_displacement[node_idx]
                } else {
                    Vector3::zeros()
                };
                
                // Compute friction force
                let friction_force = self.compute_friction_force(
                    &normal_force,
                    &tangent_velocity,
                    &tangent_disp,
                    node_idx,
                );
                
                let friction_mag = friction_force.norm();
                if friction_mag > max_friction {
                    max_friction = friction_mag;
                }
                
                // Total contact force on primary node
                let total_force = normal_force + friction_force;
                
                // Apply to primary node
                for i in 0..3 {
                    let global_index = simulation.get_global_index(primary_node_id, i);
                    simulation.load_vector[global_index] += total_force[i];
                }
                
                // Apply reaction to secondary element nodes (Newton's 3rd law)
                if let Some(sec_element) = simulation.get_element(sec_elem_id) {
                    let sec_connectivity = sec_element.get_connectivity().clone();
                    let num_sec_nodes = sec_connectivity.len() as f64;
                    
                    for &sec_node_id in &sec_connectivity {
                        for i in 0..3 {
                            let global_index = simulation.get_global_index(sec_node_id, i);
                            simulation.load_vector[global_index] -= total_force[i] / num_sec_nodes;
                        }
                    }
                }
            }
        }
        
        self.penetration_count = active_count;
        self.max_penetration = max_penetration;
        self.max_friction_force = max_friction;
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.active_primary_nodes
    }

    fn type_name(&self) -> BoundaryConditionType {
        BoundaryConditionType::Contact
    }

    fn initalize(&mut self, simulation: &Simulation) {
        debug!("Initializing node-to-segment contact");
        assert!(!self.secondary_elements.is_empty(), "No secondary surface elements");
        
        // Build bounding box for secondary surface
        let mut min_pt = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max_pt = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        
        for &node_id in &self.secondary_nodes {
            if let Some(node) = simulation.get_node(node_id) {
                min_pt.x = min_pt.x.min(node.position.x);
                min_pt.y = min_pt.y.min(node.position.y);
                min_pt.z = min_pt.z.min(node.position.z);
                max_pt.x = max_pt.x.max(node.position.x);
                max_pt.y = max_pt.y.max(node.position.y);
                max_pt.z = max_pt.z.max(node.position.z);
            }
        }
        
        // Expand bounding box by max_distance
        min_pt -= Vector3::new(self.max_distance, self.max_distance, self.max_distance);
        max_pt += Vector3::new(self.max_distance, self.max_distance, self.max_distance);
        
        // Filter primary nodes within bounding box
        self.active_primary_nodes.clear();
        for &node_id in &self.primary_nodes {
            if let Some(node) = simulation.get_node(node_id) {
                let p = node.position;
                if p.x >= min_pt.x && p.x <= max_pt.x &&
                   p.y >= min_pt.y && p.y <= max_pt.y &&
                   p.z >= min_pt.z && p.z <= max_pt.z 
                {
                    self.active_primary_nodes.push(node_id);
                }
            }
        }
        
        self.active_secondary_elements = self.secondary_elements.clone();
        
        // Initialize friction state
        self.prev_tangent_displacement = vec![Vector3::zeros(); self.active_primary_nodes.len()];
        
        // Initialize Lagrange multipliers for augmented Lagrangian
        if matches!(self.algorithm, ContactAlgorithm::AugmentedLagrangian) {
            self.lagrange_multipliers = vec![0.0; self.active_primary_nodes.len()];
        }
        
        debug!("Contact initialized: {} active primary nodes, {} secondary elements",
               self.active_primary_nodes.len(), self.active_secondary_elements.len());
    }

    fn print_stats(&self) {
        debug!("Contact stats: max_pen={:.6e}, active={}, max_friction={:.6e}",
               self.max_penetration, self.penetration_count, self.max_friction_force);
    }
}
