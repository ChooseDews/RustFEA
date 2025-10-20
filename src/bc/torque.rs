use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::Vector3;
use serde::{Serialize, Deserialize};
use crate::bc::BoundaryConditionType;

/// Represents a torque boundary condition
/// # Mathematical Formulation
/// A torque condition applies forces tangent to an axis following the right-hand rule.
/// For a given axis (defined by a point and direction) and a torque magnitude,
/// forces are applied to each node perpendicular to both the axis and the radial vector
/// from the axis to the node.
/// 
/// The force direction at each point is calculated as: `f_dir = axis_direction × r_perp`
/// where `r_perp` is the perpendicular vector from the axis to the node.
/// 
/// The force at each node is proportional to its distance from the axis (F_i = k × r_i),
/// where k is scaled such that the total torque contribution from all nodes equals exactly
/// the specified net torque magnitude: k = T_total / Σ(r_i²).
/// 
/// # Fields
/// * `nodes`: The nodes to apply the torque to.
/// * `axis_point`: A point on the axis of rotation.
/// * `axis_direction`: The direction vector of the axis (must be normalized).
/// * `net_torque`: The total torque magnitude to apply.
#[derive(Serialize, Deserialize, Debug)] 
pub struct TorqueCondition {
    nodes: Vec<usize>,
    axis_point: Vector3<f64>,
    axis_direction: Vector3<f64>,
    net_torque: f64,
}

impl TorqueCondition {
    /// Creates a new torque condition
    /// 
    /// # Arguments
    /// * `nodes` - List of node IDs to apply the torque to
    /// * `axis_point` - A point on the axis of rotation
    /// * `axis_direction` - Direction vector of the axis (will be normalized)
    /// * `net_torque` - Total torque magnitude to apply
    pub fn new(nodes: Vec<usize>, axis_point: Vector3<f64>, axis_direction: Vector3<f64>, net_torque: f64) -> Self {
        let axis_direction = axis_direction.normalize();
        TorqueCondition { 
            nodes, 
            axis_point, 
            axis_direction, 
            net_torque 
        }
    }

    /// Creates a new torque condition from vectors
    /// 
    /// # Arguments
    /// * `nodes` - List of node IDs to apply the torque to
    /// * `axis_point` - A point on the axis of rotation as [x, y, z]
    /// * `axis_direction` - Direction vector of the axis as [x, y, z] (will be normalized)
    /// * `net_torque` - Total torque magnitude to apply
    pub fn new_from_vec(nodes: Vec<usize>, axis_point: Vec<f64>, axis_direction: Vec<f64>, net_torque: f64) -> Self {
        let axis_point = Vector3::new(axis_point[0], axis_point[1], axis_point[2]);
        let axis_direction = Vector3::new(axis_direction[0], axis_direction[1], axis_direction[2]);
        TorqueCondition::new(nodes, axis_point, axis_direction, net_torque)
    }
}

#[typetag::serde]
impl BoundaryCondition for TorqueCondition {
    fn apply(&mut self, simulation: &mut Simulation) {
        // First pass: calculate radii and sum of r²
        let mut node_data: Vec<(f64, Vector3<f64>)> = Vec::with_capacity(self.nodes.len());
        let mut r_squared_sum = 0.0;
        
        for &node_id in &self.nodes {
            if let Some(node) = simulation.get_node(node_id) {
                // Vector from axis point to node
                let to_node = node.position - self.axis_point;
                
                // Project onto axis to get parallel component
                let parallel_component = to_node.dot(&self.axis_direction) * self.axis_direction;
                
                // Perpendicular component (radial vector from axis to node)
                let r_perp = to_node - parallel_component;
                let radius = r_perp.norm();
                
                r_squared_sum += radius * radius;
                node_data.push((radius, r_perp));
            } else {
                node_data.push((0.0, Vector3::zeros()));
            }
        }
        
        if r_squared_sum < 1e-12 {
            eprintln!("Warning: All nodes in torque condition are on the axis of rotation");
            return;
        }
        
        // Second pass: apply forces
        // Force is distributed such that F_i is proportional to r_i
        // This ensures: Total torque = Σ(r_i × F_i) = Σ(r_i × (k*r_i)) = k*Σ(r_i²) = net_torque
        // Therefore: k = net_torque / Σ(r_i²)
        let force_scale = self.net_torque / r_squared_sum;
        
        for (idx, &node_id) in self.nodes.iter().enumerate() {
            let (radius, r_perp) = node_data[idx];
            
            if radius < 1e-12 {
                continue; // Skip nodes on the axis
            }
            
            // Force direction using right-hand rule: axis × r_perp
            let force_direction = self.axis_direction.cross(&r_perp).normalize();
            
            // Force magnitude: F_i = (net_torque / Σ(r_i²)) * r_i
            // This gives: torque_i = r_i * F_i = r_i * force_scale * r_i = force_scale * r_i²
            // Total torque = Σ(force_scale * r_i²) = force_scale * Σ(r_i²) = net_torque ✓
            let force_magnitude = force_scale * radius;
            let force = force_direction * force_magnitude;
            
            // Apply force to load vector
            for i in 0..3 {
                let global_index = simulation.get_global_index(node_id, i);
                simulation.load_vector[global_index] += force[i];
            }
        }
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }

    fn type_name(&self) -> BoundaryConditionType {
        BoundaryConditionType::Torque
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_torque_perpendicular_force() {
        // Test that force is perpendicular to both axis and radial vector
        let axis_point: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
        let axis_direction: Vector3<f64> = Vector3::new(0.0, 0.0, 1.0);
        let node_position: Vector3<f64> = Vector3::new(1.0, 0.0, 0.0);
        
        let to_node = node_position - axis_point;
        let parallel = to_node.dot(&axis_direction) * axis_direction;
        let r_perp = to_node - parallel;
        let force_dir = axis_direction.cross(&r_perp).normalize();
        
        // Force should be in y-direction for this configuration
        let expected: Vector3<f64> = Vector3::new(0.0, 1.0, 0.0);
        assert!((force_dir - expected).norm() < 1e-10);
    }
}

