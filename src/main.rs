extern crate nalgebra as na;
mod node;
mod elements;
mod mesh_generation;
mod io;
mod simulation;


fn main() {
    // Example usage
    // let mut node_instance = node::Node::new(1, 0.0, 0.0, 0.0);
    // println!("Node before displacement: {:?}", node_instance);
    
    // node_instance.set_displacement(0.1, 0.2, 0.3);
    // println!("Node after displacement: {:?}", node_instance);

    // let total_disp = node_instance.total_displacement();
    // println!("Total displacement magnitude: {}", total_disp);

    let b = 5;
    let (nodes, elements) = mesh_generation::generate_mesh(4*b, 4*b, 10*b, 1.0, 1.0, 5.0);
    let mut simulation = simulation::Simulation::from_arrays(nodes, elements.into_iter().map(|e| Box::new(e) as Box<dyn elements::base_element::BaseElement>).collect());


    simulation.solve();

    //println!("K_global: {:?}", K_global);
    // let (global_stiffness_matrix, global_force, specified_bc) = simulation.assemble();
    // println!("Done assembling global stiffness matrix!");
    // //write matrix to file
    // match io::matrix_writer::write_hashmap_sparse_matrix("matrix.txt", &global_stiffness_matrix) {
    //     Ok(()) => println!("Matrix written successfully!"),
    //     Err(e) => println!("Error writing matrix: {}", e),
    // }


    match io::vtk_writer::write_vtk("mesh.vtk", &simulation) {
        Ok(()) => println!("VTK file written successfully!"),
        Err(e) => println!("Error writing VTK file: {}", e),
    }    
    

}


