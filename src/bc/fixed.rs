use crate::simulation::Simulation;
use crate::bc::BoundaryCondition;
use nalgebra as na;
use na::{DMatrix, DVector, Vector3};
use nalgebra_sparse::ops::Op;
use serde::{Serialize, Deserialize};
use std::fmt;


/// Represents a fixed boundary condition
/// # Mathematical Formulation
/// A fixed condition is an essential boundary condition, where the nodes are fixed to a specific value.
/// This feild can be applied to some or all nodes in the mesh.
/// # Fields
/// * `nodes`: The nodes to apply the fixed condition to.
/// * `fixed_values`: The fixed values for each node. This is an optional field as some dof may be free
#[derive(Serialize, Deserialize, Debug)]
pub struct FixedCondition {
    nodes: Vec<usize>,
    fixed_values: Vec<Option<f64>>
}

impl FixedCondition {
    pub fn new(nodes: Vec<usize>, fixed_values: Vec<Option<f64>>) -> Self {
        FixedCondition { nodes, fixed_values }
    }
    pub fn all_3d(nodes: Vec<usize>, fixed_value: f64) -> Self {
        let fixed_values = vec![Some(fixed_value); 3];
        FixedCondition::new(nodes, fixed_values)
    }
    pub fn static_3d(nodes: Vec<usize>) -> Self {
        FixedCondition::all_3d(nodes, 0.0)
    }
}

#[typetag::serde]
impl BoundaryCondition for FixedCondition {
    fn apply(&self, simulation: &mut Simulation) {
        for &node_id in &self.nodes {
            for (i, value) in self.fixed_values.iter().enumerate() {
                if value.is_none() { continue };
                let global_index: usize = simulation.get_global_index(node_id, i);
                simulation.fixed_global_nodal_values.insert(global_index, value.unwrap());
            }
        }
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.nodes
    }

    fn type_name(&self) -> &str {
        "FixedCondition"
    }
}

impl fmt::Display for FixedCondition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FixedCondition: {} nodes, {} fixed values", self.nodes.len(), self.fixed_values.len())
    }
}
