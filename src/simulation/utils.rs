use crate::elements::base_element::{BaseElement, ElementFields};
use log::{debug, info, warn};
use std::collections::HashMap;
use std::fmt;
use std::path::Path;
use crate::mesh::MeshAssembly;
use crate::node::Node;
use crate::utilities::Keywords;

use super::Simulation;
use crate::simulation::base::SimulationStep;

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

impl Default for NodeAvgValue {
    fn default() -> Self {
        Self::new()
    }
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
        if self.count == 0 {
            return 0.0;
        }
        self.value / self.count as f64
    }
}

impl Default for Simulation {
    fn default() -> Self {
        Self::new()
    }
}

impl Simulation {
    pub fn new() -> Self {
        Simulation {
            nodes: Vec::new(),
            elements: HashMap::new(),
            node_fields: HashMap::new(),
            boundary_conditions: Vec::new(),
            load_vector: Vec::new(),
            fixed_global_nodal_values: HashMap::new(),
            keywords: Keywords::new(),
            dofs: 3,
            mesh: MeshAssembly::empty(),
            active_elements: Vec::new(),
            worker_count: std::env::var("RAYON_NUM_THREADS")
                .map(|s| s.parse::<usize>().unwrap_or(8))
                .unwrap_or(8),
            steps: Vec::new(),
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

    /// Get all timesteps
    pub fn get_steps(&self) -> &Vec<SimulationStep> {
        &self.steps
    }

    /// Get element fields for a specific time
    pub fn get_step_at_time(&self, time: f64, tolerance: f64) -> Option<&HashMap<usize, ElementFields>> {
        self.steps.iter()
            .find(|step| (step.time - time).abs() < tolerance)
            .map(|step| &step.element_fields)
    }

    /// Get element fields for a specific iteration
    pub fn get_step_at_iteration(&self, iteration: u64) -> Option<&HashMap<usize, ElementFields>> {
        self.steps.iter()
            .find(|step| step.iteration == iteration)
            .map(|step| &step.element_fields)
    }

    /// Get the latest timestep
    pub fn get_latest_step(&self) -> Option<&SimulationStep> {
        self.steps.last()
    }

    /// Clear all timesteps
    pub fn clear_steps(&mut self) {
        self.steps.clear();
    }

    pub fn check_node_ordering(&self) {
        //required for correct assembly. Not required for elements as you never need to assemble by index.
        debug!("Checking node ordering");
        let n_count = self.nodes.len();
        for i in 0..n_count {
            let node = self
                .get_node(i)
                .unwrap_or_else(|| panic!("Node {} out of {} not found", i, n_count));
            if node.id != i {
                panic!(
                    "Internal Node Id {} is out of order from expected {}",
                    node.id, i
                );
            }
        }
    }

    pub fn prep_output_directory(&self) {
        let clear_directory = self
            .keywords
            .get_bool("OUTPUT_CLEAR_DIRECTORY")
            .unwrap_or(false);
        let output_vtk = self.keywords.get_string("OUTPUT_VTK");
        if output_vtk.is_none() {
            warn!("No output vtk specified");
            return;
        }
        let output_vtk = output_vtk.unwrap();
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

    pub fn compute_result_fields(&mut self) {
        let element_count = self.elements.len();
        let mut fields: Vec<String> = Vec::new();
        let mut node_fields: HashMap<String, Vec<NodeAvgValue>> = HashMap::new();

        for i in self.active_elements().iter() {
            if fields.is_empty() {
                //use first element to get & initalize fields
                fields = self.compute_element_properties(*i).get_field_names();
                if !fields.is_empty() {
                    for field in &fields {
                        let mut node_avg_values: Vec<NodeAvgValue> = Vec::new();
                        for _ in 0..self.nodes.len() {
                            node_avg_values.push(NodeAvgValue::new());
                        }
                        node_fields.insert(field.to_string(), node_avg_values);
                    }
                } else {
                    warn!("No fields found for first element {}", i);
                    continue;
                }
            }

            let element = self.get_element(*i).unwrap();
            let element_props = element.compute_element_nodal_properties(self);
            let connectivity = element.get_connectivity();
            for field in &fields {
                let field_values = element_props.field.get(field).unwrap();
                for (j, value) in field_values.iter().enumerate() {
                    let node_id = connectivity[j];
                    node_fields.get_mut(field).unwrap()[node_id].add_value(*value);
                }
            }
        }

        self.node_fields = HashMap::new();
        //avg node_fields
        for field in &fields {
            let node_avg_values = node_fields.get_mut(field).unwrap();
            let mut node_fields_value: Vec<f64> = Vec::new();
            for node_avg_value in node_avg_values {
                node_fields_value.push(node_avg_value.get_avg());
            }
            self.node_fields
                .insert(field.to_string(), node_fields_value);
        }
    }

    pub fn compute_element_properties(&self, id: usize) -> ElementFields {
        let element = self.get_element(id).unwrap();
        element.compute_element_nodal_properties(self)
    }

    pub fn print(&self) {
        debug!("{}", self);
        debug!("Boundary Conditions: {}", self.boundary_conditions.len());
        debug!("Node Fields: {}", self.node_fields.len());

        if !self.node_fields.is_empty() {
            debug!("Field names:");
            for field_name in self.node_fields.keys() {
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
            debug!("  - {:?}: {}", type_name, count);
        }
    }
}
