use log::{debug, info, trace};
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;

use crate::bc::condition::BoundaryCondition;
use crate::elements::base_element::{BaseElement, ElementFields};
use crate::io::matrix_writer::{write_hashmap_sparse_matrix, write_vector};
use crate::mesh::Mesh;
use crate::node::Node;
use crate::solver::{direct_choslky, direct_solve};

use crate::io::input_reader::Keywords;

#[derive(Serialize, Deserialize)]
pub struct Simulation {
    //mesh data
    pub nodes: Vec<Node>,
    elements: Vec<Box<dyn BaseElement>>,
    pub node_feilds: HashMap<String, Vec<f64>>,
    //boundary conditions
    pub boundary_conditions: Vec<Box<dyn BoundaryCondition>>,
    //simulation data
    #[serde(skip, default = "Vec::new")]
    pub load_vector: Vec<f64>, //global force vector
    #[serde(skip, default = "HashMap::new")]
    pub fixed_global_nodal_values: HashMap<usize, f64>,
    pub dofs: usize,
    pub keywords: Keywords,
}

impl fmt::Display for Simulation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Simulation [NodeCount: {}, ElementCount: {}, DOFs: {}]",
            self.nodes.len(),
            self.elements.len(),
            self.dofs
        )
    }
}

#[derive(Debug)]
pub struct NodeAvgValue {
    value: f64,
    count: usize,
}

impl NodeAvgValue {
    pub fn new() -> Self {
        NodeAvgValue {
            value: 0.0,
            count: 0,
        }
    }
    pub fn add_value(&mut self, value: f64) {
        self.value += value;
        self.count += 1;
    }
    pub fn get_avg(&self) -> f64 {
        self.value / self.count as f64
    }
}

impl Simulation {
    pub fn new() -> Self {
        Simulation {
            nodes: Vec::new(),
            elements: Vec::new(),
            node_feilds: HashMap::new(),
            boundary_conditions: Vec::new(),
            load_vector: Vec::new(),
            fixed_global_nodal_values: HashMap::new(),
            keywords: Keywords::new(),
            dofs: 3,
        }
    }

    //create from array of nodes and elements
    pub fn from_arrays(nodes: Vec<Node>, elements: Vec<Box<dyn BaseElement>>, dofs: usize) -> Self {
        let mut simulation = Simulation::new();
        simulation.set_dofs(dofs);
        for node in nodes {
            simulation.add_node(node);
        }
        for element in elements {
            simulation.add_element(element);
        }
        simulation
    }
    pub fn from_mesh(mesh: &Mesh, dofs: usize) -> Self {
        let elements = mesh.convert_to_elements();
        let nodes = mesh.convert_to_nodes();
        Simulation::from_arrays(nodes, elements, dofs)
    }
    pub fn get_global_index(&self, node_id: usize, dof: usize) -> usize {
        node_id * self.dofs + dof
    }

    pub fn set_keywords(&mut self, keywords: Keywords) {
        self.keywords = keywords;
    }

    pub fn set_dofs(&mut self, dofs: usize) {
        self.dofs = dofs;
    }

    //initalize function for once nodes+elements are fixed
    pub fn initialize(&mut self) {
        let n = self.nodes.len() * self.dofs;
        self.load_vector = vec![0.0; n];
    }

    pub fn add_node(&mut self, node: Node) {
        self.nodes.push(node);
    }

    pub fn add_boundary_condition(&mut self, bc: Box<dyn BoundaryCondition>) {
        self.boundary_conditions.push(bc);
    }

    pub fn add_element(&mut self, element: Box<dyn BaseElement>) {
        self.elements.push(element);
    }

    pub fn get_node(&self, id: usize) -> Option<&Node> {
        self.nodes.get(id)
    }

    pub fn get_node_mut(&mut self, id: usize) -> Option<&mut Node> {
        self.nodes.get_mut(id)
    }

    pub fn get_nodes(&self, ids: &[usize]) -> Vec<&Node> {
        ids.iter().filter_map(|&id| self.get_node(id)).collect()
    }

    pub fn get_element(&self, id: usize) -> Option<&Box<dyn BaseElement>> {
        self.elements.get(id)
    }

    pub fn get_element_mut(&mut self, id: usize) -> Option<&mut Box<dyn BaseElement>> {
        self.elements.get_mut(id)
    }

    pub fn nodes(&self) -> &Vec<Node> {
        &self.nodes
    }

    pub fn nodes_mut(&mut self) -> &mut Vec<Node> {
        &mut self.nodes
    }

    pub fn elements(&self) -> &Vec<Box<dyn BaseElement>> {
        &self.elements
    }

    pub fn elements_mut(&mut self) -> &mut Vec<Box<dyn BaseElement>> {
        &mut self.elements
    }

    pub fn get_global_force(&self) -> Vec<f64> {
        self.load_vector.clone()
    }

    pub fn handle_bc(&mut self) {
        let mut boundary_conditions = std::mem::take(&mut self.boundary_conditions);
        for bc in &mut boundary_conditions {
            bc.apply(self);
        }
        self.boundary_conditions = boundary_conditions;
    }

    pub fn get_specified_bc(&mut self) -> Vec<(usize, f64)> {
        let mut specified_bc = Vec::new();
        for (global_index, value) in &self.fixed_global_nodal_values {
            specified_bc.push((*global_index, *value));
        }
        specified_bc
    }

    pub fn assemble_global_force(&mut self) {
        self.initialize();
        self.handle_bc();
    }

    pub fn set_element_stiffness(&mut self, id: usize, stiffness: DMatrix<f64>) {
        self.get_element_mut(id).unwrap().set_stiffness(stiffness);
    }

    pub fn compute_element_stiffness(&self, id: usize) -> DMatrix<f64> {
        self.get_element(id).unwrap().compute_stiffness(&self)
    }

    pub fn compute_all_element_stiffness(&mut self) {
        for i in 0..self.elements.len() {
            self.set_element_stiffness(i, self.compute_element_stiffness(i));
        }
    }

    pub fn assemble(&mut self) -> (HashMap<(usize, usize), f64>, Vec<f64>) {
        debug!("Starting assembly process");
        let mut global_stiffness_matrix = HashMap::new();
        let dof = 3;
        self.assemble_global_force(); //populates global force and fixed global nodal values
        let mut global_force = self.get_global_force();
        let specified_bc = self.get_specified_bc();
        self.compute_all_element_stiffness();

        trace!("Assembling global stiffness matrix");
        for i in 0..self.elements.len() {
            let element = self.get_element(i).unwrap();
            let element_stiffness_matrix = element.get_stiffness();
            let element_connectivity = element.get_connectivity();
            let node_count = element_connectivity.len();
            for i in 0..node_count {
                for j in 0..node_count {
                    for k in 0..dof {
                        for l in 0..dof {
                            let global_i = element_connectivity[i] * dof + k;
                            let global_j = element_connectivity[j] * dof + l;
                            let key = (global_i, global_j);
                            let value = element_stiffness_matrix[(i * dof + k, j * dof + l)];
                            *global_stiffness_matrix.entry(key).or_insert(0.0) += value;
                        }
                    }
                }
            }
        }
        let extra_stiffness = self
            .keywords
            .get_single_float("EXTRA_STIFFNESS")
            .unwrap_or(1e12);
        for (g_index, value) in specified_bc {
            let mut v = global_stiffness_matrix[&(g_index, g_index)];
            v += extra_stiffness;
            global_stiffness_matrix.insert((g_index, g_index), v);
            if value.abs() > 0.0 {
                global_force[g_index] += value * extra_stiffness; // F = K*u so K_extra*u_extra = -F_extra
            }
        }

        debug!("Assembly process completed");
        (global_stiffness_matrix, global_force)
    }

    pub fn solve_explicit(&mut self) {
        let dt = self.keywords.get_single_float("TIMESTEP").unwrap_or(1.0);
        let total_time = self.keywords.get_single_float("TOTAL_TIME").unwrap_or(1.0);
        let time_steps = (total_time / dt) as usize;

        for step in 0..time_steps {
            println!("Time Step: {}", step);
        }
    }

    pub fn solve(&mut self) {
        info!("Starting simulation solve process");
        debug!("Assembling system");
        let (mut global_stiffness_matrix, mut global_force) = self.assemble();

        info!("Solving system");
        let u = direct_solve(&global_stiffness_matrix, &global_force);

        info!("Performing post-solve computations");
        //populate node displacements
        for node in self.nodes_mut() {
            let node_id = node.id;
            node.set_displacement(u[node_id * 3], u[node_id * 3 + 1], u[node_id * 3 + 2])
        }

        let element_count = self.elements.len();

        //average element feilds for each node
        let feilds = self.compute_element_properties(0).get_feild_names();

        //create hashmap -> [nodes]- > NodeAvgValue
        let mut node_feilds: HashMap<String, Vec<NodeAvgValue>> = HashMap::new();
        for feild in &feilds {
            let mut node_avg_values: Vec<NodeAvgValue> = Vec::new();
            for _ in 0..self.nodes.len() {
                node_avg_values.push(NodeAvgValue::new());
            }
            node_feilds.insert(feild.to_string(), node_avg_values);
        }

        for i in 0..element_count {
            let element = self.get_element(i).unwrap();
            let element_props = element.compute_element_nodal_properties(&self);
            let connectivity = element.get_connectivity();
            for feild in &feilds {
                let feild_values = element_props.field.get(feild).unwrap();
                for (j, value) in feild_values.iter().enumerate() {
                    let node_id = connectivity[j];
                    let node_avg_value = node_feilds.get_mut(feild).unwrap();
                    node_avg_value[node_id].add_value(*value);
                }
            }
        }

        self.node_feilds = HashMap::new();

        //avg node_feilds
        for feild in &feilds {
            let node_avg_values = node_feilds.get_mut(feild).unwrap();
            let mut node_feilds_value: Vec<f64> = Vec::new();
            for node_avg_value in node_avg_values {
                node_feilds_value.push(node_avg_value.get_avg());
            }
            self.node_feilds
                .insert(feild.to_string(), node_feilds_value);
        }

        // println!("{:?}", self.node_feilds);
    }

    pub fn compute_element_properties(&self, id: usize) -> ElementFields {
        let element = self.get_element(id).unwrap();
        element.compute_element_nodal_properties(&self)
    }

    pub fn print(&self) {
        debug!("{}", self);
        debug!("Boundary Conditions: {}", self.boundary_conditions.len());
        debug!("Node Fields: {}", self.node_feilds.len());

        if !self.node_feilds.is_empty() {
            debug!("Field names:");
            for field_name in self.node_feilds.keys() {
                debug!("  - {}", field_name);
            }
        }

        debug!("Load Vector Size: {}", self.load_vector.len());
        debug!(
            "Fixed Nodal Values: {}",
            self.fixed_global_nodal_values.len()
        );

        // Print element types and counts
        let mut element_types = std::collections::HashMap::new();
        for element in &self.elements {
            let type_name = element.type_name();
            *element_types.entry(type_name).or_insert(0) += 1;
        }
        debug!("Element Types:");
        for (type_name, count) in element_types {
            debug!("  - {}: {}", type_name, count);
        }

        // Print boundary condition types and counts
        let mut bc_types = std::collections::HashMap::new();
        for bc in &self.boundary_conditions {
            let type_name = bc.type_name();
            *bc_types.entry(type_name).or_insert(0) += 1;
        }
        debug!("Boundary Condition Types:");
        for (type_name, count) in bc_types {
            debug!("  - {}: {}", type_name, count);
        }
    }
}
