//Main representation of a project which holds all relevant data including potentially multiple simulations and meshes. 
use serde_json::Value;
use std::fs::File;
use std::io::BufReader;
use crate::mesh::Mesh;
use crate::simulation::Simulation; 
use serde::{Serialize, Deserialize};
use crate::io::input_reader::{Keywords, read_simulation_file};
use crate::io::file::{seralized_write, seralized_read};

use super::vtk_writer::write_vtk;
use log::{info, debug};

#[derive(Serialize, Deserialize)]
pub struct Project {
    meshes: Vec<Mesh>,
    simulations: Vec<Simulation>,
    keywords: Keywords,
}

impl Project {
    pub fn from(mesh: Mesh, simulation: Simulation) -> Self {
        Project {
            meshes: vec![mesh],
            simulations: vec![simulation],
            keywords: Keywords::new(),
        }
    }
    /// Creates a new project from an input file (.sim file)
    pub fn from_input_file(file_path: &str) -> Self {
        let (keywords, simulations, meshes) = read_simulation_file(file_path).unwrap();
        Project { keywords, simulations, meshes }
    }
    /// Saves the project to a file via serialization.
    pub fn save(&self) -> String {
        let output = self.keywords.get_single_value("OUTPUT").expect("No output file specified");
        info!("Saving project to {}", output);
        seralized_write(&output, self);
        output
    }
    /// Loads a project from a file and deserializes it to a Project struct.
    pub fn load(file_path: &str) -> Self {
        info!("Loading project from {}", file_path);
        seralized_read(file_path)
    }
    pub fn print(&self) {
        debug!("Printing project information");
        self.keywords.print();
        for mesh in &self.meshes {
            debug!("{}", mesh);
        }
        for simulation in &self.simulations {
            debug!("{}", simulation);
        }
    }

    pub fn get_simulation(&mut self) -> &mut Simulation {
        &mut self.simulations[0]
    }

    pub fn export_vtk(&self) {
       let simulation = &self.simulations[0];
       let output_vtk = self.keywords.get_single_value("OUTPUT_VTK").expect("No output vtk specified");
       info!("Exporting VTK to {}", output_vtk);
       write_vtk(output_vtk.as_str(), &simulation);
    }
}

