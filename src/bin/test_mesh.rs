use clap::Parser;
use rust_fea::io::mesh_reader::read_file;
use rust_fea::mesh::MeshAssembly;
use rust_fea::simulation::Simulation;
use rust_fea::io::vtk_writer::write_vtk;
use serde_json;
use std::fs::File;
use std::io::Write;
use rust_fea::io::file::{seralized_read, seralized_write};
use rust_fea::io::project::Project;
use log::LevelFilter;
use env_logger::Builder;


/// Mesh Reader CLI
#[derive(Parser, Debug)]
#[command(about="Import a .inp mesh", long_about=None)]
struct Args {
    #[arg(value_name = "Path to .inp file")]
    input: String,
    #[arg(short, long, help="Path to a second inp file", required=false)] //-a
    add_mesh: Option<String>,
    #[arg(short, long, help="Run function on mesh", required=false)] //-r
    run: Option<String>,
    #[arg(short, long, help="Output file name (.mesh.bin.xz file)", required=false)] //-o
    output: Option<String>,
    #[arg(short, long, help="Output file name (.vtk file)", required=false)] //-v
    paraview: Option<String>,
    #[arg(short, long, help="Verbose output", action = clap::ArgAction::Count)]
    verbose: u8,    
}


fn main() {

    let args = Args::parse();
    // Initialize the logger
    let log_level = match args.verbose {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    Builder::new().filter_level(log_level).init();

    let mut mesh: rust_fea::mesh::MeshAssembly = read_file(&args.input);
    mesh.print_info();
    let apply_func = args.run.unwrap_or("none".to_string());
    let output_file = args.input.replace(".inp", ".mesh.bin.xz");
    let vtk_file = args.paraview.unwrap_or("mesh.vtk".to_string());

    if let Some(add_mesh) = args.add_mesh {

        let add_mesh: MeshAssembly = read_file(&add_mesh);
        println!("Adding mesh! Info: ********************************");
        add_mesh.print_info();
        mesh += add_mesh;
        println!("Mesh after addition ###############################");
        mesh.print_info();

    }




    if apply_func == "sin" {
        mesh.apply_func_to_nodel_positions(&MeshAssembly::half_sine_func);
    } else if apply_func == "none" {
        println!("No function applied to mesh");
    } else {
        println!("Invalid function applied to mesh. So no function applied.");
    }
    mesh.save(&output_file);

    let simulation = Simulation::from_mesh(mesh, 3);
    write_vtk(&vtk_file, &simulation).unwrap();



}
