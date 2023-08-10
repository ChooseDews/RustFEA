// src/simulation.rs
use super::node::Node;
use crate::elements::base_element::BaseElement;
use std::collections::HashMap;
use na::{ComplexField, Vector3};
use nalgebra as na;
use nalgebra_sparse;


pub struct Simulation {
    nodes: Vec<Node>,
    //we want to support more element types in the future so we use a trait object
    elements: Vec<Box<dyn BaseElement>>,
}

impl Simulation {
    pub fn new() -> Self {
        Simulation {
            nodes: Vec::new(),
            elements: Vec::new(),
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

    pub fn nodes(&self) -> &Vec<Node> {
        &self.nodes
    }

    pub fn nodes_mut(&mut self) -> &mut Vec<Node> {
        &mut self.nodes
    }

    pub fn elements(&self) -> &Vec<Box<dyn BaseElement>> {
        &self.elements
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
        for node in self.nodes(){
            if (node.position[2] - min_z).abs() < tol{
                //fix node
                specified_bc.push((node.id, 0, 0));
                specified_bc.push((node.id, 1, 0));
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
        let f = 1e8/z_nodes.len() as f64;
        let node_force = [-f, 0.0, 0.0]; 
        for node in z_nodes{
            let node_id = node.id;
            global_force[node_id*3] = node_force[0];
            global_force[node_id*3 + 1] = node_force[1];
            global_force[node_id*3 + 2] = node_force[2];
        }


        //apply force on min-y and max-y nodes
        //find min y coordinate and max y coordinate
        let mut min_y = 0.0;
        let mut max_y = 0.0;
        for node in self.nodes(){
            if node.position[1] < min_y{
                min_y = node.position[1];
            }
            if node.position[1] > max_y{
                max_y = node.position[1];
            }
        }

        //find nodes with y coordinate closest to min_y and max_y
        let mut min_y_nodes: Vec<&Node> = Vec::new();
        let mut max_y_nodes: Vec<&Node> = Vec::new();

        for node in self.nodes(){
            if (node.position[1] - min_y).abs() < tol{
                min_y_nodes.push(node);
            }
            if (node.position[1] - max_y).abs() < tol{
                max_y_nodes.push(node);
            }
        }
        //print how many nodes in max_y_nodes and min_y_nodes
        println!("{} nodes in min_y_nodes and {} nodes in max_y_nodes", min_y_nodes.len(), max_y_nodes.len());
        

        //apply force to min_y_nodes and max_y_nodes in opposite directions
        let f = 2e8/min_y_nodes.len() as f64;
        let node_force = [-f, 0.0, 0.0];
        for node in min_y_nodes{
            let node_id = node.id;
            global_force[node_id*3] = node_force[0];
            global_force[node_id*3 + 1] = node_force[1];
            global_force[node_id*3 + 2] = node_force[2];
        }

        let f = 2e8/max_y_nodes.len() as f64;
        let node_force = [f, 0.0, 0.0];
        for node in max_y_nodes{
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
        let n = self.nodes.len()*DOF;



        let global_force = self.assemble_global_force();

        let specified_bc = self.get_specified_bc();


        //use hashmap as a sparse matrix to assemble the global stiffness matrix
        //loop over elements
        //print how many nodes and elements there are
        println!("{} nodes and {} elements", self.nodes.len(), self.elements.len());
        for element in self.elements(){
            let element_stiffness_matrix = element.compute_stiffness(&self);
            let element_force_vector = element.compute_force(&self);
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
        let max_i = self.nodes.len()*3;
        //load into nalgebra sparse matrix
        let mut sparse_matrix: nalgebra_sparse::CooMatrix<f64> = nalgebra_sparse::CooMatrix::new(max_i, max_i);
        let mut max_stiffness_value = 0.0;


        //add stiffness on specified bc to diagonal
        let extra_stiffness = 1e12;
        for (node_id, dof, value) in specified_bc{
            let i = node_id*3 + dof;
            let j = node_id*3 + dof;
            let mut v = global_stiffness_matrix[&(i, j)];
            v += extra_stiffness;
            global_stiffness_matrix.insert((i, j), v);
        }


        for ((i, j), value) in global_stiffness_matrix.iter(){
            let mut v = *value;
            if v.abs() > max_stiffness_value{
                max_stiffness_value = v.abs();
            }
            //add to diagonal if i == j
            if i == j{
                v += 0.0001; //add small value to diagonal to avoid singularity
            }
            sparse_matrix.push(*i, *j, v);
        }
        println!("max stiffness value: {}", max_stiffness_value);

        let csc = nalgebra_sparse::CscMatrix::from(&sparse_matrix);

        println!("solving system...");
        //use cholesky decomposition to solve system
        let cholesky = nalgebra_sparse::factorization::CscCholesky::factor(&csc).unwrap();
        
        //solve system
        let g_force = nalgebra::DVector::from_vec(global_force);
        let u: na::Matrix<f64, na::Dyn, na::Dyn, na::VecStorage<f64, na::Dyn, na::Dyn>> = cholesky.solve(&g_force);

        //populate node displacements
        for node in self.nodes_mut(){
            let node_id = node.id;
            node.set_displacement(u[node_id*3], u[node_id*3 + 1], u[node_id*3 + 2])
        }
        
    
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
}
