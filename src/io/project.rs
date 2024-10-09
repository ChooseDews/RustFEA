//Main representation of a project which holds all relevant data including potentially multiple simulations and meshes. 
use serde_json::Value;
use std::fs::File;
use std::io::BufReader;
use crate::mesh::MeshAssembly;
use crate::simulation::Simulation; 
use serde::{Serialize, Deserialize};
use crate::io::input_reader::read_simulation_file;
use crate::io::file::{seralized_write, seralized_read};
use crate::utilities::Keywords;
use super::vtk_writer::write_vtk;
use log::{info, debug};

#[derive(Serialize, Deserialize)]
pub struct Project {
    simulations: Vec<Simulation>,
    keywords: Keywords,
}

impl Project {
    /// Creates a new project from an input file (.sim file)
    pub fn from_input_file(file_path: &str) -> Self {
        read_simulation_file(file_path).unwrap()
    }

    pub fn new(simulations: Vec<Simulation>, keywords: Keywords) -> Self {
        Project { keywords, simulations }
    }
    /// Saves the project to a file via serialization.
    pub fn save(&self) -> String {
        let output = self.keywords.get_string("OUTPUT").expect("No output file specified");
        self.save_to_file(output.as_str());
        output
    }

    pub fn save_to_file(&self, file_path: &str) {
        info!("Saving project to {}", file_path);
        seralized_write(file_path, self);
    }
    /// Loads a project from a file and deserializes it to a Project struct.
    pub fn load(file_path: &str) -> Self {
        info!("Loading project from {}", file_path);
        seralized_read(file_path)
    }
    pub fn print(&self) {
        debug!("Printing project information");
        self.keywords.print();
        for simulation in &self.simulations {
            simulation.print();
        }
    }

    pub fn get_simulation(&mut self) -> &mut Simulation {
        &mut self.simulations[0]
    }

    pub fn export_vtk(&self) {
        for (index, simulation) in self.simulations.iter().enumerate() {
            let output_vtk: String = simulation.keywords.get_string("OUTPUT_VTK").expect("No output vtk specified");
            info!("Exporting VTK for simulation {} to {}", index, output_vtk);
            write_vtk(output_vtk.as_str(), simulation);
        }
    }

    pub fn solve(&mut self) {
        for simulation in &mut self.simulations {
            simulation.solve();
        }
    }

}

