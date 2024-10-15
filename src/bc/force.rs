use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::DVector;
use serde::{Serialize, Deserialize};
use crate::bc::BoundaryConditionType;


/// Represents a load boundary condition
/// # Mathematical Formulation
/// A load condition is a natural boundary condition, where the nodes are loaded with a specific force.
/// This feild can be applied to some or all nodes in the mesh.
/// Here we are applying a condition on the derivative of the field, in the case of solid mechanics this would be a force,
/// # Fields
/// * `nodes`: The nodes to apply the load condition to.
/// * `force`: The force to apply to the nodes.
#[derive(Serialize, Deserialize, Debug)] 
pub struct LoadCondition {
    nodes: Vec<u32>,
    force: DVector<f64>
}

impl LoadCondition {
    pub fn new(nodes: Vec<u32>, net_force: DVector<f64>) -> Self {
        let s = nodes.len() as f64;
        let force = net_force / s;
        LoadCondition { nodes, force }
    }
    pub fn new_from_vec(nodes: Vec<u32>, net_force: Vec<f64>) -> Self {
        let force = DVector::from(net_force);
        LoadCondition::new(nodes, force)
    }
}

#[typetag::serde]
impl BoundaryCondition for LoadCondition {
    fn apply(&mut self, simulation: &mut Simulation) {
        for &node_id in &self.nodes {
            for i in 0..self.force.len() {
                let global_index = simulation.get_global_index(node_id, i as u32);
                simulation.load_vector[global_index as usize] += self.force[i];
            }
        }
    }

    fn get_nodes(&self) -> &Vec<u32> {
        &self.nodes
    }

    fn type_name(&self) -> BoundaryConditionType {
        BoundaryConditionType::Force
    }
}