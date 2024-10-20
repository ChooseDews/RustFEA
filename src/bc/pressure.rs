use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::{DMatrix, DVector, Vector3};
use nalgebra_sparse::ops::Op;
use serde::{Serialize, Deserialize};
use std::fmt;
use crate::bc::BoundaryConditionType;


#[derive(Serialize, Deserialize, Debug)]
pub struct PressureCondition {
    elements: Vec<usize>,
    pressure: f64,
    #[serde(skip, default = "Vec::new")]
    nodes: Vec<usize>,
}

impl PressureCondition  {
    pub fn new(elements: Vec<usize>, pressure: f64) -> Self {
        PressureCondition { elements, pressure, nodes: Vec::new() }
    }
}

#[typetag::serde]
impl BoundaryCondition for PressureCondition {
    fn apply(&mut self, simulation: &mut Simulation) {
        todo!() //need to loop over surface elements and apply integrate pressure over the surface. We could consider making a general traction BC then pressure is just a special case of normal traction
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }

    fn type_name(&self) -> BoundaryConditionType {
        BoundaryConditionType::Pressure
    }
}

impl fmt::Display for PressureCondition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PressureCondition: {} nodes, P={}", self.nodes.len(), self.pressure)
    }
}
