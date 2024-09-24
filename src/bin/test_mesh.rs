use clap::Parser;
use rust_fea::io::mesh_reader::read_inp_file;
use rust_fea::simulation::Simulation;
use rust_fea::io::vtk_writer::write_vtk;
use serde_json;
use std::fs::File;
use std::io::Write;
use rust_fea::io::file::{write_json_file, read_json_file, seralized_write_json};
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
    let mesh: rust_fea::mesh::Mesh = read_inp_file(&args.input);
    mesh.print_info();

    let nodes = mesh.convert_to_nodes();
    println!("Number of nodes: {}", nodes.len());

    let elements = mesh.convert_to_elements();
    println!("Number of elements: {}", elements.len());

    let mut simulation = Simulation::from_arrays(nodes, elements);
    write_vtk("mesh_2.vtk", &simulation).unwrap();

    println!("Assembling problem...");
    let problem = simulation.assemble();
    println!("Problem assembled!");
    println!("Stiffness [K] value count: {}", problem.0.len());


    //save mesh to json file
    seralized_write_json("temp/mesh.json.xz", &mesh);
    seralized_write_json("temp/mesh.json", &mesh);
    seralized_write_json("temp/sim.json", &simulation);
    

    //read simulation from json file
    let sim_value = read_json_file("temp/sim.json");
    let sim_read: Simulation = serde_json::from_value(sim_value).unwrap();
    println!("Simulation read from json: {}", sim_read);


    //create project
    let project = Project::from(mesh, simulation);
    //save the project
    seralized_write_json("temp/project.json", &project);
    seralized_write_json("temp/project.json.xz", &project);

}
