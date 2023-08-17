use std::collections::HashMap;
use crate::node::Node;
use crate::elements::{ BaseElement, BrickElement, Material };


#[derive(Debug)]
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
}

#[derive(Debug)]
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




#[derive(Debug)]
pub struct ElementGroup {
    pub elements: Vec<usize>, //indices of elements in the mesh
    pub name: String,
    pub el_type: String,
}

#[derive(Debug)]
pub struct NodeGroup {
    pub nodes: Vec<usize>, //indices of nodes in the mesh
    pub name: String,
}

#[derive(Debug)]
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

    pub fn print_info(&self) {
        println!("Mesh: {}", self.name);
        println!("Number of nodes: {}", self.nodes.len());
        println!("Number of elements: {}", self.elements.len());
        println!("Number of element groups: {}", self.element_groups.len());
        let group_names = self.element_groups.keys();
        for name in group_names {
            let group = self.element_groups.get(name).unwrap();
            println!("  -> {} - {} - {} elements", group.name, group.el_type, group.elements.len());
        }
        println!("Number of node groups: {}", self.node_groups.len());
        let group_names = self.node_groups.keys();
        for name in group_names {
            let group = self.node_groups.get(name).unwrap();
            println!("  -> {} - {} nodes", group.name, group.nodes.len());
        }
    }

    pub fn convert_to_nodes(&self) -> Vec<Node> {
        let mut nodes = Vec::new();
        let n = self.nodes.len();
        for node_id in 0..n {
            let mesh_node = self.nodes.get(&node_id).unwrap();
            nodes.push(mesh_node.to_node(node_id as isize));
        }
        nodes
    }

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
