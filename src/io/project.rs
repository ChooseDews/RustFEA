//Main representation of a project which holds all relevant data including potentially multiple simulations and meshes. 
use serde_json::Value;
use std::fs::File;
use std::io::BufReader;
use crate::mesh::Mesh;
use crate::simulation::Simulation; 
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub struct Project {
    meshes: Vec<Mesh>,
    simulations: Vec<Simulation>,
}

impl Project {
    pub fn from(mesh: Mesh, simulation: Simulation) -> Self {
        Project {
            meshes: vec![mesh],
            simulations: vec![simulation],
        }
    }
}

