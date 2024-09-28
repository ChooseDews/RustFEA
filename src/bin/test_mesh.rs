use clap::Parser;
use rust_fea::io::mesh_reader::read_file;
use rust_fea::simulation::Simulation;
use rust_fea::io::vtk_writer::write_vtk;
use serde_json;
use std::fs::File;
use std::io::Write;
use rust_fea::io::file::{seralized_read, seralized_write};
use rust_fea::io::project::Project;


/// Mesh Reader CLI
#[derive(Parser, Debug)]
#[command(about="Import a .inp mesh", long_about=None)]
struct Args {
    #[arg(short, long, help="Input file name (.inp file)", required=true)]
    input: String,
}


fn main() {
    let args = Args::parse();
    let mut mesh: rust_fea::mesh::Mesh = read_file(&args.input);
    mesh.print_info();
    //mesh.apply_func();
    let nodes = mesh.convert_to_nodes();
    println!("Number of nodes: {}", nodes.len());
    let elements = mesh.convert_to_elements();
    println!("Number of elements: {}", elements.len());


    let mut simulation = Simulation::from_arrays(nodes, elements, 3);



    write_vtk("mesh_2.vtk", &simulation).unwrap();
    seralized_write("temp/output.mesh.json.xz", &mesh);

    // println!("Assembling problem...");
    // let problem = simulation.assemble();
    // println!("Problem assembled!");
    // println!("Stiffness [K] value count: {}", problem.0.len());

    // seralized_write("temp/mesh.json.xz", &mesh);
    // seralized_write("temp/sim.json.xz", &simulation);
    // let sim_read: Simulation = seralized_read("temp/sim.json.xz");
    // println!("Simulation read from json: {}", sim_read);
    // let project = Project::from(mesh, simulation);
    // seralized_write("temp/project.json.xz", &project);

}
