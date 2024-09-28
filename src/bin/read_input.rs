use clap::Parser;
use rust_fea::io::input_reader::read_simulation_file;
use rust_fea::io::project::Project;

/// Mesh Reader CLI
#[derive(Parser, Debug)]
#[command(about="Import a .sim file", long_about=None)]
struct Args {
    #[arg(short, long, help="Input file name (.sim file)", required=true)]
    input: String,
}


fn main() {
    let args = Args::parse();
    let mut project = Project::from_input_file(&args.input);
    project.get_simulation().solve();
    project.export_vtk();
    project.save();

}
