use crate::bc::BoundaryCondition;
use crate::simulation::Simulation;
use log::debug;
use nalgebra as na;
use serde::{Deserialize, Serialize};
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
    primary_nodes: Vec<usize>,
    primary_elements: Vec<usize>,
    secondary_nodes: Vec<usize>,
    secondary_elements: Vec<usize>, //currently expected to be a surface element of sometim

    //used for computation
    active_primary_nodes: Vec<usize>,
    active_secondary_elements: Vec<usize>,
}

impl NormalContact {
    pub fn new(primary_surface: String, secondary_surface: String) -> Self {
        NormalContact {
            primary_surface,
            secondary_surface,
            contact_stiffness: 100000000.0,
            max_distance: 0.25,
            active_primary_nodes: Vec::new(),
            contact_type: ContactType::Penalty,
            active_secondary_elements: Vec::new(),
            primary_nodes: Vec::new(),
            primary_elements: Vec::new(),
            secondary_nodes: Vec::new(),
            secondary_elements: Vec::new(),
        }
    }
    pub fn set_contact_surfaces_nodes(
        &mut self,
        primary_nodes: Vec<usize>,
        secondary_nodes: Vec<usize>,
    ) {
        self.primary_nodes = primary_nodes;
        self.secondary_nodes = secondary_nodes;
    }
    //Mesh id which is not necessary the same as the element id of the simulation
    pub fn set_contact_surfaces_elements(
        &mut self,
        primary_elements: Vec<usize>,
        secondary_elements: Vec<usize>,
    ) {
        self.primary_elements = primary_elements;
        self.secondary_elements = secondary_elements;
    }
}

#[typetag::serde]
impl BoundaryCondition for NormalContact {
    fn apply(&self, simulation: &mut Simulation) {
        // Implement the logic to apply the contact condition to the simulation

        //for each active node check signed distance to each secondary element
        //if the distance is less than the max distance then apply the contact condition
        //to the node and the element
        for primary_node in &self.active_primary_nodes {
            let primary_node = simulation.nodes[*primary_node];
            let mut min_distance = f64::INFINITY;
            let mut vector_distance = na::Vector3::new(0.0, 0.0, 0.0);
            for secondary_element in &self.active_secondary_elements {
                let secondary_element = &simulation.elements[*secondary_element];
                let distance = secondary_element.get_signed_distance_vector(primary_node.position, simulation);
                if distance.is_some() {
                    let distance = distance.unwrap();
                    let distance_norm = distance.norm();
                    if distance_norm < min_distance {
                        min_distance = distance_norm;
                        vector_distance = distance;
                    }
                }
            }
            println!("Min distance: {}", min_distance);
            println!("Vector distance: {}", vector_distance);
        }
    }

    fn get_nodes(&self) -> &Vec<usize> {
        &self.active_primary_nodes
    }
    fn type_name(&self) -> &str {
        "ContactCondition"
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
            let node = simulation.nodes[*node];
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
            let node = simulation.nodes[*node];
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

        println!("{:?}", self);
    }
}
