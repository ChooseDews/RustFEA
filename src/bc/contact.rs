use crate::{bc::BoundaryCondition, utilities::check_for_nans};
use crate::simulation::Simulation;
use log::debug;
use nalgebra::{self as na, DVector};
use serde::{Deserialize, Serialize};
use crate::bc::BoundaryConditionType;


#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ContactType {
    Penalty,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NormalContact {
    //surface name ids
    primary_surface: String,
    secondary_surface: String,

    //contact parameters
    contact_stiffness: f64,
    max_distance: f64,
    contact_type: ContactType,

    //info
    primary_nodes: Vec<u32>,
    primary_elements: Vec<u32>,
    secondary_nodes: Vec<u32>,
    secondary_elements: Vec<u32>, //currently expected to be a surface element of sometim

    //used for computation
    active_primary_nodes: Vec<u32>,
    active_secondary_elements: Vec<u32>,
    
    penetration_count: u32,
    max_penetration: f64,
}

impl NormalContact {
    pub fn new(primary_surface: String, secondary_surface: String) -> Self {
        NormalContact {
            primary_surface,
            secondary_surface,
            contact_stiffness: 10000000.0,
            max_distance: 0.25,
            active_primary_nodes: Vec::new(),
            contact_type: ContactType::Penalty,
            active_secondary_elements: Vec::new(),
            primary_nodes: Vec::new(),
            primary_elements: Vec::new(),
            secondary_nodes: Vec::new(),
            secondary_elements: Vec::new(),
            penetration_count: 0,
            max_penetration: 0.0,
        }
    }
    pub fn set_contact_surfaces_nodes(
        &mut self,
        primary_nodes: Vec<u32>,
        secondary_nodes: Vec<u32>,
    ) {
        self.primary_nodes = primary_nodes;
        self.secondary_nodes = secondary_nodes;
    }
    //Mesh id which is not necessary the same as the element id of the simulation
    pub fn set_contact_surfaces_elements(
        &mut self,
        primary_elements: Vec<u32>,
        secondary_elements: Vec<u32>,
    ) {
        self.primary_elements = primary_elements;
        self.secondary_elements = secondary_elements;
    }
    pub fn set_penetration_count(&mut self, count: u32) {
        self.penetration_count = count;
    }
    pub fn set_max_penetration(&mut self, max_penetration: f64) {
        self.max_penetration = max_penetration;
    }
}

#[typetag::serde]
impl BoundaryCondition for NormalContact {
    fn apply(&mut self, simulation: &mut Simulation) {
        // Implement the logic to apply the contact condition to the simulation

        //for each active node check signed distance to each secondary element
        //if the distance is less than the max distance then apply the contact condition
        //to the node and the element
        let mut max_penetration: f64 = 0.0;
        let mut active_nodes = 0;
        for primary_node_index in &self.active_primary_nodes {
            let primary_node = simulation.get_node(*primary_node_index).unwrap();
            let mut load_vector = None;
            for secondary_element in &self.active_secondary_elements {
                let secondary_element = simulation.get_element(*secondary_element).unwrap();
                let point = primary_node.position + primary_node.displacement;
                let distance = secondary_element.get_signed_distance_vector(point, simulation);
                if distance.is_some() {
                    let signed_distance = distance.unwrap();
                    let contact_stiffness = self.contact_stiffness;
                    let contact_force = contact_stiffness * signed_distance;
                    load_vector = Some(contact_force);
                    let magnitude = signed_distance.magnitude();
                    if magnitude > max_penetration {
                        max_penetration = magnitude;
                    }
                    active_nodes += 1;
                    break;
                }
            }
            if load_vector.is_some() {
                for i in 0..3 {
                    let global_index = simulation.get_global_index(*primary_node_index, i);
                    // assert!(!check_for_nans(&load_vector.unwrap()), "#2 contactforce vector contains NaNs, node: {} dof: {}", *primary_node_index, i);
                    simulation.load_vector[global_index as usize] += load_vector.unwrap()[i as usize];
                }
            }

        }
        self.set_penetration_count(active_nodes);
        self.set_max_penetration(max_penetration);
        //  debug!("max penetration: {} active nodes: {}", max_penetration, active_nodes);
    }

    fn get_nodes(&self) -> &Vec<u32> {
        &self.active_primary_nodes
    }
    fn type_name(&self) -> BoundaryConditionType {
        BoundaryConditionType::Contact
    }
    fn initalize(&mut self, simulation: &Simulation) {
        debug!("Initializing contact condition");
        assert!(
            self.secondary_elements.len() > 0,
            "No main 2nd surface elements found"
        );
        //rule out elements which are not close by the max distance from center
        let mut min_point = na::Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max_point: nalgebra::Matrix<
            f64,
            nalgebra::Const<3>,
            nalgebra::Const<1>,
            nalgebra::ArrayStorage<f64, 3, 1>,
        > = na::Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for node in &self.secondary_nodes {
            let node = simulation.get_node(*node).unwrap();
            min_point.x = min_point.x.min(node.position.x);
            min_point.y = min_point.y.min(node.position.y);
            min_point.z = min_point.z.min(node.position.z);
            max_point.x = max_point.x.max(node.position.x);
            max_point.y = max_point.y.max(node.position.y);
            max_point.z = max_point.z.max(node.position.z);
        }
        min_point = min_point.add_scalar(-self.max_distance);
        max_point = max_point.add_scalar(self.max_distance);

        for node in &self.primary_nodes {
            let node = simulation.get_node(*node).unwrap();
            if node.position.x >= min_point.x
                && node.position.x <= max_point.x
                && node.position.y >= min_point.y
                && node.position.y <= max_point.y
                && node.position.z >= min_point.z
                && node.position.z <= max_point.z
            {
                self.active_primary_nodes.push(node.id);
            }
        }

        self.active_secondary_elements = self.secondary_elements.clone();
        debug!("{:?}", self);
    }

    fn print_stats(&self){
        debug!("Contact stats: max penetration: {} active nodes: {}", self.max_penetration, self.penetration_count);
    }
}
