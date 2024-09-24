use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::DVector;

pub struct ForceCondition {
    nodes: Vec<usize>,
    force: DVector<f64>,
}

impl ForceCondition {
    pub fn new(nodes: Vec<usize>, force: DVector<f64>) -> Self {
        ForceCondition { nodes, force }
    }
}

impl BoundaryCondition for ForceCondition {
    fn apply(&self, simulation: &mut Simulation) {
        for &node_id in &self.nodes {
            //panic not implemented
            panic!("Not implemented");
        }
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }
}
