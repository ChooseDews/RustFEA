use crate::bc::BoundaryCondition;
use crate::simulation::Simulation;
use serde::{Serialize, Deserialize};


#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ContactType {
    Penalty
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactCondition {
    nodes: Vec<usize>,
    contact_surfaces: Vec<String>,  
    contact_type: ContactType,
    contact_surfaces_nodes: (Vec<f64>, Vec<f64>)
}

#[typetag::serde]
impl BoundaryCondition for ContactCondition {
    fn apply(&self, simulation: &mut Simulation) {
        // Implement the logic to apply the contact condition to the simulation
    }
    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }
    fn type_name(&self) -> &str {
        "ContactCondition"
    }
}

