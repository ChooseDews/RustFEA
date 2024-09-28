use crate::simulation::Simulation;


/// Represents a boundary condition
#[typetag::serde]
pub trait BoundaryCondition {
    fn apply(&self, simulation: &mut Simulation);
    fn get_nodes(&self) -> &Vec<usize>;
    fn type_name(&self) -> &str;
}

