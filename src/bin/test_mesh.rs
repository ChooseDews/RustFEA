use clap::Parser;
use rust_fea::io::mesh_reader::read_inp_file;
use rust_fea::simulation::Simulation;
use rust_fea::io::vtk_writer::write_vtk;


/// Mesh Reader CLI
#[derive(Parser, Debug)]
#[command(about="Import a .inp mesh", long_about=None)]
struct Args {
    #[arg(short, long, help="Input file name (.inp file)", required=true)]
    input: String,
}


fn main() {
    let args = Args::parse();
    let mesh = read_inp_file(&args.input);
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



}