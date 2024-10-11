use log::{debug, info, trace, warn};
use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::ops::Op;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::path::Path;

use crate::bc::condition::BoundaryCondition;
use crate::elements::base_element::{BaseElement, ElementFields};
use crate::io::matrix_writer::{write_hashmap_sparse_matrix, write_vector};
use crate::mesh::MeshAssembly;
use crate::node::Node;
use crate::solver::{direct_choslky, direct_solve};

use crate::io::vtk_writer::write_vtk;
use crate::utilities::Keywords;
#[derive(Serialize, Deserialize)]
pub struct Simulation {
    pub nodes: Vec<Node>,
    pub node_feilds: HashMap<String, Vec<f64>>,
    pub elements: HashMap<usize, Box<dyn BaseElement>>,

    pub dofs: usize,
    pub keywords: Keywords,
    pub mesh: MeshAssembly,
    pub boundary_conditions: Vec<Box<dyn BoundaryCondition>>,

    #[serde(skip, default = "Vec::new")]
    pub load_vector: Vec<f64>, //global force vector
    #[serde(skip, default = "HashMap::new")]
    pub fixed_global_nodal_values: HashMap<usize, f64>,

    #[serde(skip, default = "Vec::new")]
    pub active_elements: Vec<usize>,
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
            elements: HashMap::new(),
            node_feilds: HashMap::new(),
            boundary_conditions: Vec::new(),
            load_vector: Vec::new(),
            fixed_global_nodal_values: HashMap::new(),
            keywords: Keywords::new(),
            dofs: 3,
            mesh: MeshAssembly::empty(),
            active_elements: Vec::new(),
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

    pub fn check_node_ordering(&self) { //required for correct assembly. Not required for elements as you never need to assemble by index. TODO: Could be replaced by mapping at solve time
        debug!("Checking node ordering");
        let n_count = self.nodes.len();
        for i in 0..n_count {
            let node = self.get_node(i).expect(format!("Node {} out of {} not found", i, n_count).as_str());
            if node.id != i {
                panic!("Internal Node Id {} is out of order from expected {}", node.id, i);
            }
        }
    }

    pub fn set_mesh(&mut self, mesh: MeshAssembly) {
        self.mesh = mesh;
    }

    pub fn from_mesh(mesh: MeshAssembly, dofs: usize) -> Self {
        let elements = mesh.convert_to_elements();
        let nodes = mesh.convert_to_nodes();
        let mut simulation = Simulation::from_arrays(nodes, elements, dofs);
        simulation.set_mesh(mesh);
        simulation
    }

    pub fn get_global_index(&self, node_id: usize, dof: usize) -> usize {
        node_id * self.dofs + dof
    }

    pub fn get_global_index_mut(&mut self, node_id: usize, dof: usize) -> usize {
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

    pub fn one_time_init(&mut self) {
        self.set_active_elements();
        self.init_bc();
        self.check_node_ordering();
    }

    pub fn init_bc(&mut self) {
        //detach bc from simulation then pass in self to
        let mut boundary_conditions = std::mem::take(&mut self.boundary_conditions);
        for bc in &mut boundary_conditions {
            bc.initalize(self);
        }
        self.boundary_conditions = boundary_conditions;
    }

    pub fn add_node(&mut self, node: Node) {
        self.nodes.insert(node.id, node);
    }

    pub fn add_boundary_condition(&mut self, bc: Box<dyn BoundaryCondition>) {
        self.boundary_conditions.push(bc);
    }

    pub fn add_element(&mut self, element: Box<dyn BaseElement>) {
        self.elements.insert(element.get_id(), element);
    }

    pub fn get_node(&self, index: usize) -> Option<&Node> {
        self.nodes.get(index)
    }

    pub fn get_node_mut(&mut self, index: usize) -> Option<&mut Node> {
        self.nodes.get_mut(index)
    }

    pub fn get_nodes(&self, ids: &[usize]) -> Vec<&Node> {
        ids.iter().filter_map(|&id| self.get_node(id)).collect()
    }

    pub fn get_element(&self, index: usize) -> Option<&Box<dyn BaseElement>> {
        self.elements.get(&index)
    }

    pub fn get_element_mut(&mut self, id: usize) -> Option<&mut Box<dyn BaseElement>> {
        self.elements.get_mut(&id)
    }

    pub fn nodes(&self) -> &Vec<Node> {
        &self.nodes
    }

    pub fn nodes_mut(&mut self) -> &mut Vec<Node> {
        &mut self.nodes
    }

    pub fn elements(&self) -> &HashMap<usize, Box<dyn BaseElement>> {
        &self.elements
    }

    pub fn elements_mut(&mut self) -> &mut HashMap<usize, Box<dyn BaseElement>> {
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
    
    pub fn set_active_elements(&mut self) {
        //collect element.id that are active
        self.active_elements = Vec::new();
        for (id, element) in self.elements.iter() {
            if element.is_active() {
                self.active_elements.push(*id);
            }
        }
        debug!("Total Active Elements: {}", self.active_elements.len());
    }

    pub fn active_elements(&self) -> Vec<usize> {
        self.active_elements.clone()
    }


    pub fn compute_all_element_stiffness(&mut self) {
        let mut elements = std::mem::take(&mut self.elements);
        for (_, element) in &mut elements {
            element.compute_stiffness(self);
        }
        self.elements = elements;
    }

    pub fn compute_element_mass(&self, id: usize) -> DMatrix<f64> {
        self.get_element(id).unwrap().compute_mass(&self)
    }

    pub fn set_element_mass(&mut self, id: usize, mass: DMatrix<f64>) {
        self.get_element_mut(id).unwrap().set_mass(mass);
    }

    pub fn compute_all_element_mass(&mut self) {
        let mut total_mass = 0.0;
        for i in self.active_elements().iter() {
            let mass = self.compute_element_mass(*i);
            let element = self.get_element_mut(*i).unwrap();
            total_mass += element.set_lumped_mass(&mass);
            element.set_mass(mass);
        }
        debug!("Total mass: {}", total_mass);
    }

    pub fn assemble(&mut self) -> (HashMap<(usize, usize), f64>, Vec<f64>) {
        debug!("Starting assembly process");
        let mut global_stiffness_matrix = HashMap::new();
        let dof: usize = 3;
        self.assemble_global_force(); //populates global force and fixed global nodal values
        let mut global_force = self.get_global_force();
        let specified_bc = self.get_specified_bc();
        self.compute_all_element_stiffness();
        self.compute_all_element_mass();

        trace!("Assembling global stiffness matrix");
        for i in self.active_elements().iter() {
            let element = self.get_element(*i).unwrap();
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
            .get_float("SOLVER_EXTRA_STIFFNESS")
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

    pub fn compute_global_mass_matrix_diagonal(&self) -> DVector<f64> {
        let total_dofs = self.nodes.len() * self.dofs;
        let mut global_mass_diagonal = DVector::zeros(total_dofs);
        for i in self.active_elements().iter() {
            let element = self.get_element(*i).unwrap();
            let element_mass = element.get_lumped_mass();
            let connectivity = element.get_connectivity();
            for (local_index, &node_id) in connectivity.iter().enumerate() {
                for dof in 0..self.dofs {
                    let global_index = self.get_global_index(node_id, dof);
                    global_mass_diagonal[global_index] += element_mass[local_index];
                }
            }
        }
        global_mass_diagonal
    }

    pub fn compute_force_vector(&self) -> DVector<f64> {
        let mut force_vector = DVector::zeros(self.nodes.len() * self.dofs);
        for i in self.active_elements().iter() {
            let element = self.get_element(*i).unwrap();
            let force = element.compute_force(&self);
            let connectivity = element.get_connectivity();
            for (local_index, &node_id) in connectivity.iter().enumerate() {
                for dof in 0..self.dofs {
                    let global_index = self.get_global_index(node_id, dof);
                    force_vector[global_index] += force[local_index * self.dofs + dof];
                }
            }
        }
        force_vector
    }

    pub fn displacement_vector(&mut self) -> DVector<f64> {
        let dofs = self.dofs;
        let mut displacement_vector = DVector::zeros(self.nodes.len() * dofs);
        for (node_id, node) in self.nodes.iter().enumerate() {
            for dof in 0..dofs {
                let global_index = self.get_global_index(node_id, dof);
                displacement_vector[global_index] = node.displacement[dof];
            }
        }
        displacement_vector
    }

    pub fn prep_output_directory(&self) {
        let clear_directory = self
            .keywords
            .get_bool("OUTPUT_CLEAR_DIRECTORY")
            .unwrap_or(false);
        let output_vtk = self
            .keywords
            .get_string("OUTPUT_VTK")
            .expect("No output vtk specified");
        let output_vtk_dir = Path::new(&output_vtk).parent().unwrap();
        if !output_vtk_dir.exists() {
            info!("Creating output directory: {}", output_vtk_dir.display());
            std::fs::create_dir_all(output_vtk_dir).expect("Failed to create output directory")
        }
        if clear_directory {
            info!("Clearing output directory: {}", output_vtk_dir.display());
            std::fs::remove_dir_all(output_vtk_dir).expect("Failed to clear output directory");
            std::fs::create_dir_all(output_vtk_dir).expect("Failed to create output directory")
        }
    }

    pub fn solve_explicit(&mut self) {
        self.prep_output_directory();
        let output_vtk = self
            .keywords
            .get_string("OUTPUT_VTK")
            .expect("No output vtk specified");

        //compute fundumental dt
        let vel = self.get_element(self.active_elements()[0]).unwrap().get_material().get_wave_speed();
        let fundumental_dt = self.mesh.compute_dt(vel);
        let dt = self
            .keywords
            .get_float("SOLVER_TIME_STEP")
            .unwrap_or(fundumental_dt * 0.75);
        if dt > fundumental_dt {
            warn!(
                "Time step is greater than the fundumental dt: {} > {}",
                dt, fundumental_dt
            );
        }
        let time_steps = self.keywords.get_int("SOLVER_TIME_STEPS").unwrap_or(100);
        let print_steps = self.keywords.get_int("SOLVER_PRINT_STEPS").unwrap_or(0);
        let vtk_save_steps = self.keywords.get_int("SOLVER_VTK_SAVE_STEPS").unwrap_or(0); //0 means no vtk saving

        //we will not assemble the global stiffness matrix
        let dof = 3;
        self.assemble_global_force(); //populates global force and fixed global nodal values
        let specified_bc = self.get_specified_bc();
        self.compute_all_element_stiffness();
        self.compute_all_element_mass();
        let global_mass_matrix_diagonal = self.compute_global_mass_matrix_diagonal();
        for i in 0..self.nodes.len() {
            self.get_node_mut(i)
                .expect("Node not found")
                .set_mass(global_mass_matrix_diagonal[i]);
        }
        let n = self.nodes.len() * self.dofs;
        debug!("Steps: {}; dt: {}; DOF: {}", time_steps, dt, n);
        //start doing the time marching
        let mut u = self.displacement_vector();
        let mut u_dot = DVector::zeros(self.nodes.len() * self.dofs);
        let mut u_half_dot = DVector::zeros(self.nodes.len() * self.dofs);

        let mut i = 0;
        let mut t = 0.0;
        while i < time_steps {
            // Compute forces
            self.assemble_global_force(); //populates global force and fixed global nodal values
            let external_force = DVector::from_vec(self.load_vector.clone());
            let internal_force = self.compute_force_vector();
            let residual_force = &external_force - &internal_force;

            // Update velocity and position
            let u_dotdot = residual_force.component_div(&global_mass_matrix_diagonal);
            u_half_dot = &u_dot + 0.5 * dt * &u_dotdot;
            u += dt * &u_half_dot;

            // Apply boundary conditions
            for (global_index, value) in &specified_bc {
                u[*global_index] = *value;
                u_dot[*global_index] = (value - u[*global_index]) / dt;
                u_half_dot[*global_index] = u_dot[*global_index];
            }

            // Update nodal displacements
            for (node_id, node) in self.nodes_mut().iter_mut().enumerate() {
                node.set_displacement(u[node_id * 3], u[node_id * 3 + 1], u[node_id * 3 + 2]);
            }

            // Compute new forces and update velocity
            let new_internal_force = self.compute_force_vector();
            let new_residual_force = &external_force - &new_internal_force;
            let new_u_dotdot = new_residual_force.component_div(&global_mass_matrix_diagonal);
            u_dot = &u_half_dot + 0.5 * dt * &new_u_dotdot;
            u_dot *= 0.9995;

            if vtk_save_steps > 0 && i % vtk_save_steps == 0 {
                self.compute_result_feilds();
                let output_vtk = output_vtk.replace(".vtk", &format!("_step_{}.vtk", i));
                write_vtk(output_vtk.as_str(), &self);
            }
            if print_steps > 0 && i % print_steps == 0 {
                println!("[{}] Time: {}", i, t);
                let res_vals: Vec<String> = residual_force.iter().map(|x| x.to_string()).collect();
                println!("Residual force: {}", res_vals.join(", "));
                let u_vals: Vec<String> = u.iter().map(|x| x.to_string()).collect();
                println!("\n Displacment: {}", u_vals.join(", "));
            }
            i += 1;
            t += dt;
        }

        self.compute_result_feilds();
    }

    pub fn solve(&mut self) {
        let start_time = std::time::Instant::now();
        let method = self
            .keywords
            .get_string("SOLVER_METHOD")
            .unwrap_or("direct".to_string());

        self.one_time_init();
        match method.as_str() {
            "direct" => self.solve_direct(),
            "explicit" => self.solve_explicit(),
            _ => panic!("Invalid solve method"),
        };
        let duration = start_time.elapsed();
        debug!("Finished Simulation Solve. Took: {:?}", duration);
    }

    pub fn solve_direct(&mut self) {
        //direct solve
        info!("Starting simulation solve process");
        debug!("Assembling system");
        let (mut global_stiffness_matrix, mut global_force) = self.assemble();
        info!("Solving system");
        let u = direct_solve(&global_stiffness_matrix, &global_force);
        info!("Performing post-solve computations");
        //populate node displacements
        for (node_id, node) in self.nodes_mut().iter_mut().enumerate() {
            node.set_displacement(u[node_id * 3], u[node_id * 3 + 1], u[node_id * 3 + 2])
        }
        self.compute_result_feilds();
    }

    pub fn compute_result_feilds(&mut self) {
        let element_count = self.elements.len();

        //average element feilds for each node
        let mut feilds: Vec<String> = Vec::new();
        let mut node_feilds: HashMap<String, Vec<NodeAvgValue>> = HashMap::new();


        for i in self.active_elements().iter() {

            if feilds.is_empty() { //use first element to get & initalize fields
                feilds = self.compute_element_properties(*i).get_feild_names();
                if !feilds.is_empty() {
                    for feild in &feilds {
                        let mut node_avg_values: Vec<NodeAvgValue> = Vec::new();
                        for _ in 0..self.nodes.len() {
                        node_avg_values.push(NodeAvgValue::new());
                        }
                        node_feilds.insert(feild.to_string(), node_avg_values);
                    }
                }else{
                    warn!("No feilds found for first element {}", i);
                    continue;
                }
            }



            let element = self.get_element(*i).unwrap();
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
        for (_index, element) in self.elements.iter() {
            let type_name = element.type_name();
            *element_types.entry(type_name).or_insert(0) += 1;
        }
        debug!("Element Types:");
        for (type_name, count) in element_types {
            debug!("  - {:?}: {}", type_name, count);
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