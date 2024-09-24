// src/simulation.rs
use super::node::Node;
use crate::elements::base_element::{BaseElement, ElementFields};
use std::collections::HashMap;
use crate::solver::{direct_solve, direct_choslky};
use crate::io::matrix_writer::{write_hashmap_sparse_matrix, write_vector};
use serde::{Serialize, Deserialize};
use std::fmt;

#[derive(Serialize, Deserialize)]
pub struct Simulation {

    //mesh data
    pub nodes: Vec<Node>,
    
    elements: Vec<Box<dyn BaseElement>>,

    //feild output data
    pub node_feilds: HashMap<String, Vec<f64>>,

    //boundary conditions
}

impl fmt::Display for Simulation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Customize the output format as needed
        write!(f, "Simulation [NodeCount: {}, ElementCount: {}]", self.nodes.len(), self.elements.len())
    }
}


#[derive(Debug)]
pub struct NodeAvgValue {
    value: f64,
    count: usize
}

impl NodeAvgValue {
    pub fn new() -> Self {
        NodeAvgValue {
            value: 0.0,
            count: 0
        }
    }
    pub fn add_value(&mut self, value: f64) {
        self.value += value;
        self.count += 1;
    }
    pub fn get_avg(&self) -> f64 {
        self.value/self.count as f64
    }
}
    

impl Simulation {
    pub fn new() -> Self {
        Simulation {
            nodes: Vec::new(),
            elements: Vec::new(),
            node_feilds: HashMap::new()
        }
    }

    pub fn add_node(&mut self, node: Node) {
        self.nodes.push(node);
    }

    pub fn add_element(&mut self, element: Box<dyn BaseElement>) {
        self.elements.push(element);
    }

    // Get a reference to a node by its ID
    pub fn get_node(&self, id: usize) -> Option<&Node> {
        self.nodes.get(id)
    }

    //get mut ref
    pub fn get_node_mut(&mut self, id: usize) -> Option<&mut Node> {
        self.nodes.get_mut(id)
    }

    pub fn get_nodes(&self, ids: &[usize]) -> Vec<&Node> {
        let mut nodes = Vec::new();
        for id in ids {
            nodes.push(self.get_node(*id).unwrap());
        }
        nodes
    }

    pub fn get_element(&self, id: usize) -> Option<&Box<dyn BaseElement>> {
        self.elements.get(id)
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

    pub fn get_specified_bc(&mut self) -> Vec<(usize, usize, usize)> {
        //return a vector of tuples of node id and degree of freedom and value
        //for now, just return an empty vector
        let mut specified_bc: Vec<(usize, usize, usize)> = Vec::new();

        //for now just fix nodes on z=0 plane
        //find min z coordinate
        let mut min_z = 0.0;
        for node in self.nodes(){
            if node.position[2] < min_z{
                min_z = node.position[2];
            }
        }  
        //find nodes with z coordinate closest to min_z

        let tol = 1e-4;
        let _fixed_nodes = 0;
        for node in self.nodes(){
            if (node.position[2] - min_z).abs() < tol{
                //fix node
                specified_bc.push((node.id, 1, 0));
                specified_bc.push((node.id, 0, 0));
                specified_bc.push((node.id, 2, 0));
        
            }
        }

        specified_bc
    }
    

    pub fn assemble_global_force(&self) -> Vec<f64> {
        
        //for now select z most nodes
        //find max z coordinate
        let mut max_z = 0.0;
        for node in self.nodes(){
            if node.position[2] > max_z{
                max_z = node.position[2];
            }
        }

        //find nodes with z coordinate closest to max_z
        let mut z_nodes: Vec<&Node> = Vec::new();
        let tol = 1e-4;
        for node in self.nodes(){
            if (node.position[2] - max_z).abs() < tol{
                z_nodes.push(node);
            }
        }

        //print z_nodes

        let n = self.nodes.len()*3; //number of degrees of freedom
        let mut global_force: Vec<f64> = vec![0.0; n];
        //loop over nodes and add force to global_force
        let f = 1e7/z_nodes.len() as f64;
        let node_force = [f, 0.0, 0.0]; 
        for node in z_nodes{
            let node_id = node.id;
            global_force[node_id*3] = node_force[0];
            global_force[node_id*3 + 1] = node_force[1];
            global_force[node_id*3 + 2] = node_force[2];
        }


    
        global_force

    }

    pub fn assemble(&mut self) -> (HashMap<(usize, usize), f64>, Vec<f64>, Vec<(usize, usize, usize)>) {
        //assemble the global stiffness matrix and force vector
        //return a hashmap of node id to index in the global stiffness matrix
        //and the global stiffness matrix and force vector
        //for now, just return an empty hashmap and two empty vectors
        let mut global_stiffness_matrix: HashMap<(usize, usize), f64> = HashMap::new();
        let DOF = 3; //degrees of freedom per node (x, y, z)
        let _n = self.nodes.len()*DOF;
        let global_force = self.assemble_global_force();
        let specified_bc = self.get_specified_bc();


        //use hashmap as a sparse matrix to assemble the global stiffness matrix
        //loop over elements
        //print how many nodes and elements there are
        println!("{} nodes and {} elements", self.nodes.len(), self.elements.len());
        for element in self.elements(){
            let element_stiffness_matrix = element.compute_stiffness(&self);
            let _element_force_vector = element.compute_force(&self);
            let element_connectivity = element.get_connectivity();
            let node_count = element_connectivity.len();
            //loop over element_stiffness_matrix and append to global_stiffness_matrix according to element_connectivity
            for i in 0..node_count{
                for j in 0..node_count{
                    for k in 0..DOF{
                        for l in 0..DOF{
                            let global_i = element_connectivity[i]*DOF + k;
                            let global_j = element_connectivity[j]*DOF + l;
                            let key = (global_i, global_j);
                            let value = element_stiffness_matrix[(i*DOF + k, j*DOF + l)];
                            //check if key exists in hashmap
                            if global_stiffness_matrix.contains_key(&key){
                                //add value to existing key
                                global_stiffness_matrix.insert(key, global_stiffness_matrix[&key] + value);
                            }
                            else{
                                //create new key
                                global_stiffness_matrix.insert(key, value);
                            }
                       }
                    }
                }
            }

        }
        (global_stiffness_matrix, global_force, specified_bc)
    }

    

    pub fn solve(&mut self) {



        //assemble global stiffness matrix and force vector
        let (mut global_stiffness_matrix, global_force, specified_bc) = self.assemble();

        println!("Solving System...");
        
        //add stiffness on specified bc to diagonal
        let extra_stiffness = 1e12;
        for (node_id, dof, _value) in specified_bc{
            let i = node_id*3 + dof;
            let j = node_id*3 + dof;
            let mut v = global_stiffness_matrix[&(i, j)];
            v += extra_stiffness;
            global_stiffness_matrix.insert((i, j), v);
        }

        let u = direct_solve(&global_stiffness_matrix, &global_force);
        //let u = direct_choslky(&global_stiffness_matrix, global_force);

        //write global stiff matrix
        write_hashmap_sparse_matrix("temp/big.matrix",&global_stiffness_matrix).unwrap();
        write_vector("temp/big.force", &global_force).unwrap();

        // //compate u and u_new
        // let mut max_diff = 0.0;
        // for (i, value) in u.iter().enumerate(){
        //     let diff = (value - u_new[i]).abs();
        //     if diff > max_diff{
        //         max_diff = diff;
        //     }
        // }
        // println!("max diff: {}", max_diff);

        
        println!("Post Solve Computation...");

        //populate node displacements
        for node in self.nodes_mut(){
            let node_id = node.id;
            node.set_displacement(u[node_id*3], u[node_id*3 + 1], u[node_id*3 + 2])
        }
        
        let element_count = self.elements.len();
     
        //average element feilds for each node
        let feilds = self.compute_element_properties(0).get_feild_names();
    
        //create hashmap -> [nodes]- > NodeAvgValue
        let mut node_feilds: HashMap<String, Vec<NodeAvgValue>> = HashMap::new();
        for feild in &feilds{
            let mut node_avg_values: Vec<NodeAvgValue> = Vec::new();
            for _ in 0..self.nodes.len(){
                node_avg_values.push(NodeAvgValue::new());
            }
            node_feilds.insert(feild.to_string(), node_avg_values);
        }

        for i in 0..element_count{
            let element = self.get_element(i).unwrap();
            let element_props = element.compute_element_nodal_properties(&self);
            let connectivity = element.get_connectivity();
            for feild in &feilds{
                let feild_values = element_props.field.get(feild).unwrap();
                for (j, value) in feild_values.iter().enumerate(){
                    let node_id = connectivity[j];
                    let node_avg_value = node_feilds.get_mut(feild).unwrap();
                    node_avg_value[node_id].add_value(*value);
                }
            }
        }

        self.node_feilds = HashMap::new();

        //avg node_feilds
        for feild in &feilds{
            let node_avg_values = node_feilds.get_mut(feild).unwrap();
            let mut node_feilds_value: Vec<f64> = Vec::new();
            for node_avg_value in node_avg_values{
                node_feilds_value.push(node_avg_value.get_avg());
            }
            self.node_feilds.insert(feild.to_string(), node_feilds_value);
        }

       // println!("{:?}", self.node_feilds);

        
    
    }

    //create from array of nodes and elements
    pub fn from_arrays(nodes: Vec<Node>, elements: Vec<Box<dyn BaseElement>>) -> Self {
        let mut simulation = Simulation::new();
        for node in nodes {
            simulation.add_node(node);
        }
        for element in elements {
            simulation.add_element(element);
        }
        simulation
    }

    pub fn compute_element_properties(&self, id: usize) -> ElementFields {
        let element = self.get_element(id).unwrap();
        element.compute_element_nodal_properties(&self)
    }


}
