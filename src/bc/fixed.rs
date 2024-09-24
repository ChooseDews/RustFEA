use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::{DMatrix, DVector, Vector3};

pub struct FixedCondition {
    nodes: Vec<usize>,
}

impl FixedCondition {
    pub fn new(nodes: Vec<usize>) -> Self {
        FixedCondition { nodes }
    }
}

impl BoundaryCondition for FixedCondition {
    fn apply(&self, simulation: &mut Simulation) {
        for &node_id in &self.nodes {
            let node = simulation.get_node_mut(node_id).unwrap();
            node.set_displacement(0.0, 0.0, 0.0);
        }
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }
}
