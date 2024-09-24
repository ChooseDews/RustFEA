use crate::simulation::Simulation;

pub trait BoundaryCondition {
    fn apply(&self, simulation: &mut Simulation);
    fn get_nodes(&self) -> &Vec<usize>;
}
