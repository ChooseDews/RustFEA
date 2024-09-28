use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::DVector;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)] 
pub struct LoadCondition {
    nodes: Vec<usize>,
    force: DVector<f64>
}

impl LoadCondition {
    pub fn new(nodes: Vec<usize>, net_force: DVector<f64>) -> Self {
        let s = nodes.len() as f64;
        let force = net_force / s;
        LoadCondition { nodes, force }
    }
    pub fn new_from_vec(nodes: Vec<usize>, net_force: Vec<f64>) -> Self {
        let force = DVector::from(net_force);
        LoadCondition::new(nodes, force)
    }
}

#[typetag::serde]
impl BoundaryCondition for LoadCondition {
    fn apply(&self, simulation: &mut Simulation) {
        for &node_id in &self.nodes {
            for i in 0..self.force.len() {
                let global_index = simulation.get_global_index(node_id, i);
                simulation.load_vector[global_index] += self.force[i];
            }
        }
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }

    fn type_name(&self) -> &str {
        "LoadCondition"
    }
}
