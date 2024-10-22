use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::{DMatrix, DVector, Vector3};
use nalgebra_sparse::ops::Op;
use serde::{Serialize, Deserialize};
use std::fmt;
use crate::bc::BoundaryConditionType;


#[derive(Serialize, Deserialize, Debug)]
pub struct Traction {
    elements: Vec<usize>,
    traction: Vector3<f64>,
    #[serde(skip, default = "Vec::new")]
    nodes: Vec<usize>,
}

impl Traction  {
    pub fn new(elements: Vec<usize>, traction: Vector3<f64>) -> Self {
        Traction { elements, traction, nodes: Vec::new() }
    }
}

#[typetag::serde]
impl BoundaryCondition for Traction {
    fn apply(&mut self, simulation: &mut Simulation) {
        todo!() //need to loop over surface elements and apply integrate pressure over the surface. We could consider making a general traction BC then pressure is just a special case of normal traction
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }

    fn type_name(&self) -> BoundaryConditionType {
        BoundaryConditionType::Traction
    }
}

impl fmt::Display for Traction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Traction: {} nodes, T={}", self.nodes.len(), self.traction)
    }
}
