use std::collections::HashMap;
use crate::node::Node;
use crate::elements::{ BaseElement, BrickElement, Material };
use serde::{Serialize, Deserialize};
use crate::io::file::{seralized_read, seralized_write};
use std::fmt;
use log::{debug, info, trace};
#[derive(Debug, Serialize, Deserialize)]

pub struct MeshNode {
    pub coordinates: Vec<f64>,
    pub id: usize,
}

impl MeshNode {
    pub fn to_node(&self, id: isize) -> Node {
        let mut n_id = self.id;
        if id > -1{
            n_id = id as usize;
        }
        Node::new(n_id, self.coordinates[0], self.coordinates[1], self.coordinates[2])
    }
    pub fn distance(&self, other: &MeshNode) -> f64 {
        let dx = self.coordinates[0] - other.coordinates[0];
        let dy = self.coordinates[1] - other.coordinates[1];
        let dz = self.coordinates[2] - other.coordinates[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MeshElement {
    pub connectivity: Vec<usize>,
    pub name: String,
    pub el_type: String,
    pub id: usize
}


impl MeshElement {
    pub fn to_element(&self, id: isize) -> Box<dyn BaseElement> {
        match self.el_type.as_str() {
            "C3D8" => {
                let mut connectivity = Vec::new();
                for node_id in &self.connectivity {
                    connectivity.push(*node_id);
                }
                let mut el_id = self.id;
                if id > -1{
                    el_id = id as usize;
                }
                Box::new(BrickElement::new(el_id, self.connectivity.clone(), Material::aluminum()))
            }
            _ => {
                panic!("Element type not supported: {}", self.el_type);
            }
        }
    }
}




#[derive(Debug, Serialize, Deserialize)]
pub struct ElementGroup {
    pub elements: Vec<usize>, //indices of elements in the mesh
    pub name: String,
    pub el_type: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct NodeGroup {
    pub nodes: Vec<usize>, //indices of nodes in the mesh
    pub name: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Mesh {
    pub nodes: HashMap<usize, MeshNode>,
    pub elements: HashMap<usize, MeshElement>,
    pub element_groups: HashMap<String, ElementGroup>,
    pub node_groups: HashMap<String, NodeGroup>,
    pub name: String,
}

impl Mesh {
    pub fn empty() -> Self {
        Mesh {
            nodes: HashMap::new(),
            elements: HashMap::new(),
            element_groups: HashMap::new(),
            node_groups: HashMap::new(),
            name: "None".to_string(),
        }
    }

    /// Loads a mesh from a file.
    /// 
    /// This function reads a serialized mesh from the specified file and returns a new `Mesh` instance.
    /// 
    /// # Arguments
    /// * `filename`: The path to the file containing the serialized mesh.
    /// 
    /// # Returns
    /// A new `Mesh` instance.
    pub fn load(filename: &str) -> Self {
        let mesh: Mesh = seralized_read(filename);
        mesh
    }

    pub fn save(&self, filename: &str) {
        seralized_write(filename, &self);
    }

    pub fn get_nodes_in_group(&self, group_name: &str) -> Vec<usize> {
        let group = self.node_groups.get(group_name).expect(&format!("Node group: {} not found in mesh", group_name));
        group.nodes.clone()
    }

    //curve with sin(x) for fun
    pub fn half_sine_func(node: &Vec<f64>) -> Vec<f64> {
        let x = node[0];
        let y = node[1];
        let z = node[2];
        vec![x, y + 0.5 * z.sin(), z]
    }

    pub fn apply_func_to_nodel_positions(&mut self, func: &dyn Fn(&Vec<f64>) -> Vec<f64>) {
        for node in self.nodes.values_mut() {
            let new_pos = func(&node.coordinates);
            node.coordinates = new_pos;
        }
    }

    pub fn apply_func(&mut self) {
        self.apply_func_to_nodel_positions(&Mesh::half_sine_func);
    }

    pub fn print_info(&self) {
        debug!("{}", self);
        for (name, group) in &self.element_groups {
            debug!("  -> {} - {} - {} elements", group.name, group.el_type, group.elements.len());
        }
        for (name, group) in &self.node_groups {
            debug!("  -> {} - {} nodes", group.name, group.nodes.len());
        }
    }

    //get closest pair of nodes
    pub fn get_closest_nodes(&self) -> (usize, usize, f64) {
        let mut min_dist = std::f64::MAX;
        let mut closest_nodes = (0, 0, 0.0);
        let nodes: Vec<(&usize, &MeshNode)> = self.nodes.iter().collect();
        for (i, (&id1, node1)) in nodes.iter().enumerate() {
            for (&id2, node2) in nodes[i+1..].iter() { //avoid duplicate pairs
                let dist = node1.distance(node2);
                if dist < min_dist {
                    min_dist = dist;
                    closest_nodes = (id1, id2, dist);
                }
            }
        }
        closest_nodes
    }

    //compute fundumental dt
    pub fn compute_dt(&self, wave_speed: f64) -> f64 {
        let (node1, node2, dist) = self.get_closest_nodes();
        let dt = dist / wave_speed;
        debug!("Closest nodes: {} and {} with distance {} and dt: {}", node1, node2, dist, dt);
        dt
    }


    /// Converts the mesh nodes to a vector of nodes.
    /// 
    /// This function iterates through the nodes in the mesh, converts each mesh node to a `Node` object,
    /// and adds them to the `nodes` vector.
    /// 
    /// # Returns
    /// A vector of `Node` objects.
    pub fn convert_to_nodes(&self) -> Vec<Node> {
        let mut nodes = Vec::new();
        let n = self.nodes.len();
        for node_id in 0..n {
            let mesh_node = self.nodes.get(&node_id).unwrap();
            nodes.push(mesh_node.to_node(node_id as isize));
        }
        nodes
    }
    /// Converts the mesh elements to a vector of elements.
    /// 
    /// This function iterates through the elements in the mesh, converts each mesh element to a `BaseElement` object,
    /// and adds them to the `elements` vector.
    /// 
    /// # Returns
    /// A vector of `BaseElement` objects.
    pub fn convert_to_elements(&self) -> Vec<Box<dyn BaseElement>> {
        let mut elements = Vec::new();
        let mut i = 0;
        //find element group that contains volume elements
        let mut volume_group = None;
        for group in self.element_groups.values() {
            if group.el_type == "C3D8" {
                volume_group = Some(group);
                break;
            }
        }
        if volume_group.is_none() {
            panic!("No volume elements found in mesh!");
        }
        let volume_group = volume_group.unwrap();
        for element_id in &volume_group.elements {
            let mesh_element = self.elements.get(element_id).unwrap();
            elements.push(mesh_element.to_element(i));
            i += 1;
        }
        elements
    }

}

impl fmt::Display for Mesh {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Mesh: {} [Nodes: {}, Elements: {}, Element Groups: {}, Node Groups: {}]", 
               self.name, self.nodes.len(), self.elements.len(), 
               self.element_groups.len(), self.node_groups.len())
    }
}
