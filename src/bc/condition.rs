use crate::simulation::Simulation;
use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq, Hash)]
pub enum BoundaryConditionType {
    Fixed,
    Force,
    Contact
}

/// Represents a boundary condition
#[typetag::serde]
pub trait BoundaryCondition {
    fn apply(&mut self, simulation: &mut Simulation);
    fn initalize(&mut self, simulation: &Simulation){}
    fn get_nodes(&self) -> &Vec<usize>;
    fn type_name(&self) -> BoundaryConditionType;
    fn print_stats(&self){}
}

