extern crate nalgebra as na;
mod node;
mod elements;
mod mesh_generation;
mod io;
mod simulation;
mod utilities;
mod mesh;
mod solver;
mod bc;

use std::time::{SystemTime, UNIX_EPOCH};
use io::file::seralized_write;
use env_logger::Env;
use log::{info, debug, error};

fn main() {
    env_logger::init_from_env(Env::default().default_filter_or("info"));
    println!("Hello, world!");  
    // Example usage
    // let mut node_instance = node::Node::new(1, 0.0, 0.0, 0.0);
    // println!("Node before displacement: {:?}", node_instance);
    
    // node_instance.set_displacement(0.1, 0.2, 0.3);
    // println!("Node after displacement: {:?}", node_instance);

    // let total_disp = node_instance.total_displacement();
    // println!("Total displacement magnitude: {}", total_disp);

    let material = elements::base_element::Material::aluminum();
    let b = 4;

    let base = 2.0;
    let height = 2.0;
    let length = 15.0;

    let (nodes, elements) = mesh_generation::generate_mesh(b, 2*b, 10*b, height,  base, length);
    let mut simulation = simulation::Simulation::from_arrays(nodes, elements.into_iter().map(|e| Box::new(e) as Box<dyn elements::base_element::BaseElement>).collect(), 3);
    simulation.solve();

    let load = 1e7;
    let modulus = material.youngs_modulus;
    let I = (1.0/12.0) * base * height.powi(3);
    let disp_max = length.powi(3) * load / (3.0 * modulus * I);
    info!("Analytical displacement: {}", disp_max);

    //find node with max displacement and print it
    let mut max_disp = 0.0;
    let mut max_node = 0;
    for (id, node) in simulation.nodes.iter().enumerate() {
        let disp = node.displacement.x;
        if disp > max_disp {
            max_disp = disp;
            max_node = node.id;
        }
    }

    info!("Max displacement: {} at node {}", max_disp, max_node);
    let error = (disp_max - max_disp).abs() / disp_max;
    info!("Error: {}", error);

    //println!("K_global: {:?}", K_global);
    // let (global_stiffness_matrix, global_force, specified_bc) = simulation.assemble();
    // println!("Done assembling global stiffness matrix!");
    // //write matrix to file
    // match io::matrix_writer::write_hashmap_sparse_matrix("matrix.txt", &global_stiffness_matrix) {
    //     Ok(()) => println!("Matrix written successfully!"),
    //     Err(e) => println!("Error writing matrix: {}", e),
    // }

    let current_epoch = SystemTime::now().duration_since(UNIX_EPOCH).unwrap();
    let filename = format!("temp/{:?}_mesh.vtk", current_epoch);
    match io::vtk_writer::write_vtk(filename.as_str(), &simulation) {
        Ok(()) => info!("VTK file written successfully!"),
        Err(e) => error!("Error writing VTK file: {}", e),
    }    

    //save project
    debug!("Saving simulation data");
    seralized_write("temp/simulation_run.json", &simulation);
    seralized_write("temp/simulation_run.json.xz", &simulation);
    seralized_write("temp/simulation_run.bin", &simulation);
    seralized_write("temp/simulation_run.bin.xz", &simulation);
    info!("Simulation data saved successfully");
}


