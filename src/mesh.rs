use std::collections::HashMap;
use crate::node::Node;
use crate::elements::{ BaseElement, BrickElement, Material, FourNodeElement};
use serde::{Serialize, Deserialize};
use crate::io::file::{seralized_read, seralized_write};
use std::fmt;
use log::{debug, info, trace};
use std::ops::AddAssign;
use std::collections::HashSet;



#[derive(Debug, Serialize, Deserialize, Clone)]

pub struct MeshNode {
    pub coordinates: Vec<f64>,
    pub id: u32,
}

impl MeshNode {
    pub fn to_node(&self, id: isize) -> Node {
        let mut n_id = self.id;
        if id > -1{
            n_id = id as u32;
        }
        Node::new(n_id, self.coordinates[0], self.coordinates[1], self.coordinates[2])
    }
    pub fn distance(&self, other: &MeshNode) -> f64 {
        let dx = self.coordinates[0] - other.coordinates[0];
        let dy = self.coordinates[1] - other.coordinates[1];
        let dz = self.coordinates[2] - other.coordinates[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    pub fn offset(&self, offset: &HashMap<u32, u32>) -> MeshNode {
        MeshNode {
            coordinates: self.coordinates.clone(),
            id: offset.get(&self.id).unwrap().clone()
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MeshElement {
    pub connectivity: Vec<u32>,
    pub name: String,
    pub el_type: String,
    pub id: u32
}


impl MeshElement {
    pub fn to_element(&self) -> Box<dyn BaseElement> {
        match self.el_type.as_str() {
            "C3D8" => {
                let mut connectivity = Vec::new();
                for node_id in &self.connectivity {
                    connectivity.push(*node_id);
                }
                let mut el_id = self.id;
                Box::new(BrickElement::new(el_id, self.connectivity.clone(), Material::aluminum()))
            }
            "CPS4" => {
                let mut connectivity = Vec::new();
                for node_id in &self.connectivity {
                    connectivity.push(*node_id);
                }
                let mut el_id = self.id;
                Box::new(FourNodeElement::new(el_id, self.connectivity.clone(), Material::empty()))
            }
            _ => {
                panic!("Element type not supported: {}", self.el_type);
            }
        }
    }

    pub fn offset(&self, el_offset: &HashMap<u32, u32>, node_offset: &HashMap<u32, u32>) -> MeshElement {
        let mut connectivity = Vec::new();
        for node_id in &self.connectivity {
            connectivity.push(node_offset.get(node_id).unwrap().clone());
        }
        MeshElement {
            connectivity: connectivity,
            name: self.name.clone(),
            el_type: self.el_type.clone(),
            id: el_offset.get(&self.id).unwrap().clone(),
        }
    }
}




#[derive(Debug, Serialize, Deserialize)]
pub struct ElementGroup {
    pub elements: Vec<u32>, //indices of elements in the mesh
    pub name: String,
    pub el_type: String,
}

impl ElementGroup {
    pub fn offset(&self, el_offset: &HashMap<u32, u32>) -> ElementGroup {
        ElementGroup {
            elements: self.elements.iter().map(|e| el_offset.get(e).unwrap().clone()).collect(),
            name: self.name.clone(),
            el_type: self.el_type.clone(),
        }
    }
}
#[derive(Debug, Serialize, Deserialize)]
pub struct NodeGroup {
    pub nodes: Vec<u32>, //indices of nodes in the mesh
    pub name: String,
}

impl NodeGroup {
    pub fn offset(&self, node_offset: &HashMap<u32, u32>) -> NodeGroup {
        NodeGroup {
            nodes: self.nodes.iter().map(|n| node_offset.get(n).unwrap().clone()).collect(),
            name: self.name.clone(),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Body {
    pub elements: Vec<u32>,
    pub nodes: Vec<u32>,
    pub name: String
}

impl Body {
    pub fn offset(&self, el_offset: &HashMap<u32, u32>, node_offset: &HashMap<u32, u32> ) -> Body {
        Body {
            elements: self.elements.iter().map(|e| el_offset.get(e).unwrap().clone()).collect(),
            nodes: self.nodes.iter().map(|n| node_offset.get(n).unwrap().clone()).collect(),
            name: self.name.clone(),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MeshAssembly {
    pub nodes: HashMap<u32, MeshNode>,
    pub elements: HashMap<u32, MeshElement>,
    pub element_groups: HashMap<String, ElementGroup>,
    pub node_groups: HashMap<String, NodeGroup>,
    pub bodies: Vec<Body>,
    pub name: String,
}


impl MeshAssembly {
    pub fn empty() -> Self {
        MeshAssembly {
            nodes: HashMap::new(),
            elements: HashMap::new(),
            element_groups: HashMap::new(),
            node_groups: HashMap::new(),
            bodies: Vec::new(),
            name: "Unamed".to_string(),
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
        let mesh: MeshAssembly = seralized_read(filename);
        mesh
    }

    pub fn single_body(&mut self) {
        let element_ids = self.elements.keys().cloned().collect::<Vec<u32>>();
        let node_ids = self.nodes.keys().cloned().collect::<Vec<u32>>();
        let mut body = Body {
            elements: element_ids,
            nodes: node_ids,
            name: "Body".to_string(),
        };
        self.bodies.push(body);
    }

    pub fn multiple_bodies(&mut self, bodies: Vec<String>) {
        for body_name in bodies {
            let el_group = self.element_groups.get(&body_name).expect(&format!("Element group: {} not found in mesh", body_name));
            let node_group = self.node_groups.get(&body_name).expect(&format!("Node group: {} not found in mesh", body_name));
            let body = Body {
                elements: el_group.elements.clone(),
                nodes: node_group.nodes.clone(),
                name: body_name.clone(),
            };
            self.bodies.push(body);
        }
    }

    pub fn set_body_name(&mut self, name: &str) {
        self.bodies[0].name = name.to_string();
    }

    pub fn save(&self, filename: &str) {
        seralized_write(filename, &self);
    }

    pub fn get_nodes_in_group(&self, group_name: &str) -> Vec<u32> {
        let group = self.node_groups.get(group_name).expect(&format!("Node group: {} not found in mesh", group_name));
        group.nodes.clone()
    }

    pub fn get_elements_in_group(&self, group_name: &str) -> Vec<u32> {
        let group = self.element_groups.get(group_name).expect(&format!("Element group: {} not found in mesh", group_name));
        group.elements.clone()
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
        self.apply_func_to_nodel_positions(&MeshAssembly::half_sine_func);
    }

    pub fn print_info(&self) {
        debug!("{}", self);
        debug!("Element groups: ");
        for (name, group) in &self.element_groups {
            debug!("  -> {} - {} - {} elements", name, group.el_type, group.elements.len());
        }
        debug!("Node groups: ");
        for (name, group) in &self.node_groups {
            debug!("  -> {} - {} nodes", name, group.nodes.len());
        }
        debug!("Bodies: ");
        for body in &self.bodies {
            debug!("  -> {} - {} elements, {} nodes", body.name, body.elements.len(), body.nodes.len());
        }
    }

    //get closest pair of nodes
    pub fn get_closest_nodes(&self) -> (u32, u32, f64) {
        let overlap_tol = 1e-6;
        let mut min_dist = std::f64::MAX;
        let mut closest_nodes = (0, 0, 0.0);
        let nodes: Vec<(&u32, &MeshNode)> = self.nodes.iter().collect();
        for (i, (&id1, node1)) in nodes.iter().enumerate() {
            for (&id2, node2) in nodes[i+1..].iter() { //avoid duplicate pairs
                let dist = node1.distance(node2);
                if dist < min_dist && dist > overlap_tol {
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
            let node_id = node_id as u32;
            let mesh_node: &MeshNode = self.nodes.get(&node_id).unwrap();
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
    pub fn convert_to_brick_elements(&self) -> Vec<Box<dyn BaseElement>> {
        let mut elements = Vec::new();
        let mut element_index = 0;

        // Collect all C3D8 element groups
        let volume_groups: Vec<&ElementGroup> = self.element_groups.values()
            .filter(|group| group.el_type == "C3D8")
            .collect();

        if volume_groups.is_empty() {
            panic!("No volume elements (C3D8) found in mesh!");
        }

        // Process all C3D8 elements from all relevant groups
        for group in volume_groups {
            for element_id in &group.elements {
                if let Some(mesh_element) = self.elements.get(element_id) {
                    elements.push(mesh_element.to_element());
                    element_index += 1;
                } else {
                    debug!("Element with id {} not found in mesh", element_id);
                }
            }
        }

        debug!("Converted {} C3D8 elements", elements.len());
        elements
    }


    pub fn convert_to_four_node_elements(&self) -> Vec<Box<dyn BaseElement>> {
        let mut elements = Vec::new();
        let mut element_index = 0;

        // Collect all CPS4 element groups
        let volume_groups: Vec<&ElementGroup> = self.element_groups.values()
            .filter(|group| group.el_type == "CPS4")
            .collect();

        if volume_groups.is_empty() {
            return vec![];
        }

        for group in volume_groups {
            for element_id in &group.elements {
                if let Some(mesh_element) = self.elements.get(element_id) {
                    elements.push(mesh_element.to_element());
                } else {
                    debug!("Element with id {} not found in mesh", element_id);
                }
            }
        }

        debug!("Converted {} CPS4 elements", elements.len());
        elements
    }


    pub fn convert_to_elements(&self) -> Vec<Box<dyn BaseElement>> {
        let mut elements = Vec::new();
        elements.extend(self.convert_to_brick_elements());
        elements.extend(self.convert_to_four_node_elements());
        elements
    }


    
    //handle a RHS add of another mesh. Other mesh assumed to start from 0 node id and 0 element id
    pub fn add_mesh(&mut self, other: &MeshAssembly) {
        let min_node_id = other.nodes.keys().min().unwrap();
        let max_node_id = other.nodes.keys().max().unwrap();
        let max_current_id = self.nodes.keys().max().unwrap();
        let node_offset = *max_current_id + 1;

        let min_element_id = other.elements.keys().min().unwrap();
        let max_element_id = other.elements.keys().max().unwrap();
        let max_el_current_id = self.elements.keys().max().unwrap();
        let element_offset = *max_el_current_id + 1;

        //create element_id and node_id maps
        let mut element_id_map = HashMap::new();
        let mut node_id_map = HashMap::new();

        //for nodes and elements its (key) -> (node_offset+i)
        for (i, element) in other.elements.iter() {
            element_id_map.insert(element.id, element_offset + i);
        }
        for (i, node) in other.nodes.iter() {
            node_id_map.insert(node.id, node_offset + i);
        }
        
        for node in other.nodes.values() {
            let new_node: MeshNode = node.offset(&node_id_map);
            self.nodes.insert(new_node.id, new_node);
        }

        for element in other.elements.values() {
            let new_element: MeshElement = element.offset(&element_id_map, &node_id_map);
            self.elements.insert(new_element.id, new_element);
        }

        let mut existing_group_names: HashSet<String> = self.element_groups.keys().cloned().collect();
        existing_group_names.extend(self.node_groups.keys().cloned());

        for group in other.element_groups.values() {
            let new_group = group.offset(&element_id_map);
            let unique_name = self.generate_unique_name(&new_group.name, &existing_group_names);
            self.element_groups.insert(unique_name.clone(), new_group);
            existing_group_names.insert(unique_name);
        }

        for group in other.node_groups.values() {
            let new_group = group.offset(&node_id_map);
            let unique_name = self.generate_unique_name(&new_group.name, &existing_group_names);
            self.node_groups.insert(unique_name.clone(), new_group);
            existing_group_names.insert(unique_name);
        }

        for body in other.bodies.iter() {
            let new_body = body.offset(&element_id_map, &node_id_map);
            self.bodies.push(new_body);
        }
    }

    fn generate_unique_name(&self, base_name: &str, existing_names: &HashSet<String>) -> String {
        if !existing_names.contains(base_name) {
            return base_name.to_string();
        }

        let mut counter = 1;
        loop {
            let new_name = format!("{}_{}", base_name, counter);
            if !existing_names.contains(&new_name) {
                return new_name;
            }
            counter += 1;
        }
    }


}


impl AddAssign for MeshAssembly {
    fn add_assign(&mut self, other: Self) {
        self.add_mesh(&other);
    }
}


impl fmt::Display for MeshAssembly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Mesh: {} [Nodes: {}, Elements: {}, Element Groups: {}, Node Groups: {}]", 
               self.name, self.nodes.len(), self.elements.len(), 
               self.element_groups.len(), self.node_groups.len())
    }
}
