use clap::Parser;
use rust_fea::io::input_reader::read_simulation_file;
use rust_fea::io::project::Project;
use log::LevelFilter;
use env_logger::Builder;

/// Mesh Reader CLI
#[derive(Parser, Debug)]
#[command(about="Import a .sim file", long_about=None)]
struct Args {
    #[arg(value_name = "Input file name (.sim file)")]
    input: String,
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

    let mut project = Project::from_input_file(&args.input);
    project.solve();
    project.export_vtk();
    project.save();



    
}
