use clap::Parser;
use rust_fea::io::mesh_reader::read_file;
use rust_fea::mesh::Mesh;
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
    #[arg(value_name = "Path to .inp file")]
    input: String,
    #[arg(short, long, help="Run function on mesh", required=false)]
    run: Option<String>,
    #[arg(short, long, help="Output file name (.mesh.bin.xz file)", required=false)]
    output: Option<String>,
}


fn main() {
    let args = Args::parse();
    let mut mesh: rust_fea::mesh::Mesh = read_file(&args.input);
    mesh.print_info();
    let apply_func = args.run.unwrap_or("none".to_string());
    let output_file = args.input.replace(".inp", ".mesh.bin.xz");
    if apply_func == "sin" {
        mesh.apply_func_to_nodel_positions(&Mesh::half_sine_func);
    } else if apply_func == "none" {
        println!("No function applied to mesh");
    } else {
        println!("Invalid function applied to mesh. So no function applied.");
    }
    mesh.save(&output_file);



    // let nodes = mesh.convert_to_nodes();
    // println!("Number of nodes: {}", nodes.len());
    // let elements = mesh.convert_to_elements();
    // println!("Number of elements: {}", elements.len());
    // let mut simulation = Simulation::from_arrays(nodes, elements, 3);
    // write_vtk("mesh_2.vtk", &simulation).unwrap();
    // seralized_write("temp/output.mesh.json.xz", &mesh);

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
