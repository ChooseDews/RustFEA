use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::collections::HashMap;
use crate::simulation::Simulation;
use crate::bc::{BoundaryCondition, FixedCondition, LoadCondition};
use crate::io::mesh_reader::read_file;
use std::error::Error;
use serde::{Serialize, Deserialize};
use crate::mesh::Mesh;
use log::{info, debug, warn};

fn remove_comments(line: &str) -> String {
    match line.find('#') {
        Some(index) => line[..index].trim().to_string(),
        None => line.trim().to_string(),
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Keyword {
    keyword: String,
    values: Vec<String>,
}

impl Keyword {
    pub fn print(&self) {
        debug!("{:}={:?}", self.keyword, self.values);
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Keywords {
    values: Vec<Keyword>,
}

impl Keywords {
    pub fn new() -> Self {
        Keywords {
            values: Vec::new()
        }
    }
    /// add a keyword to the keywords
    pub fn add_keyword(&mut self, keyword: String, values: Vec<String>) {
        self.values.push(Keyword { keyword: keyword.to_uppercase(), values });
    }
    /// get the values of the first instance of a keyword
    pub fn get_keyword(&self, keyword: &str) -> Option<&Keyword> {
        self.values.iter().find(|k| k.keyword == keyword.to_uppercase())
    }
    /// Checks if a keyword exists in the keywords.
    pub fn keyword_exists(&self, keyword: &str) -> bool {
        self.values.iter().any(|k| k.keyword == keyword.to_uppercase())
    }   
    /// Retrieves a single value from the keywords. Will join the values with a space.
    pub fn get_single_value(&self, keyword: &str) -> Option<String> {
        if let Some(keyword) = self.values.iter().find(|k| k.keyword == keyword.to_uppercase()) {
            return Some(keyword.values.join(" "));
        }
        None
    }
    /// Retrieves all keywords with the specified keyword.
    pub fn get_keywords(&self, keyword: &str) -> Vec<&Keyword> { //some keywords may have multiple values or instances
        self.values.iter().filter(|k| k.keyword == keyword.to_uppercase()).collect()
    }
    /// Retrieves a single float value from the keywords. Will get the first value and try to parse it as a float.'
    pub fn get_single_float(&self, keyword: &str) -> Option<f64> {
        if let Some(keyword) = self.values.iter().find(|k| k.keyword == keyword.to_uppercase()) {
            if keyword.values[0].parse::<f64>().is_ok() {
                return Some(keyword.values[0].parse::<f64>().unwrap());
            }
        }
        None
    }

    pub fn get_keyword_values_as_floats(&self, keyword: &str) -> Option<Vec<f64>> {
        if let Some(keyword) = self.values.iter().find(|k| k.keyword == keyword.to_uppercase()) {
            return Some(keyword.values.iter().map(|s| s.parse::<f64>().unwrap()).collect());
        }
        None
    }

    pub fn add(&mut self, keyword: &str, values: Vec<&str>) {
        let keyword = Keyword { 
            keyword: keyword.to_uppercase(), 
            values: values.iter().map(|s| s.to_string()).collect() 
        };
        self.values.push(keyword);
    }

    pub fn print(&self) {
        debug!("Keywords ({}):", self.values.len());
        for keyword in &self.values {
            keyword.print();
        }
    }
}

/// Reads a simulation file and returns the keywords, simulations, and meshes.
/// 
/// This function reads a simulation file, processes its content, and extracts relevant information such as material, solver, output, mesh file, name, version, fixed boundary conditions, load boundary conditions, degrees of freedom, and output VTK settings.
/// It also handles the mesh file, reading it from the specified path and converting it into nodes and elements.
/// The function then creates a simulation object with the extracted data and adds the necessary boundary conditions.
/// Finally, it returns the keywords, a single simulation, and the mesh.
///
/// # Arguments
/// * `file_path`: The path to the simulation file. .sim
/// 
/// # Returns
/// A tuple containing the keywords, a single simulation, and the mesh.
pub fn read_simulation_file(file_path: &str) -> Result<(Keywords, Vec<Simulation>, Vec<Mesh>), Box<dyn Error>> {
    info!("Reading simulation file: {}", file_path);
    if !file_path.ends_with(".sim") {
        warn!("File does not have a .sim extension");
        return Err("File must have a .sim extension".into());
    }
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut simulation: Simulation = Simulation::new();
    let mut boundary_condition_data: Vec<(String, String, Vec<String>)> = Vec::new();
    let mut keywords = Keywords::new(); 

    let general_keywords = vec!["MATERIAL", "SOLVER", "OUTPUT", "MESH_FILE", "NAME", "VERSION", "FIXED_BC", "LOAD_BC", "DOF", "OUTPUT_VTK"];

    for line in reader.lines() {
        let line = remove_comments(&line?).trim().to_string();
        if line.is_empty() {continue;}
        let parts: Vec<&str> = line.split_whitespace().collect();
        match parts[0].to_uppercase().as_str() {
            keyword if general_keywords.contains(&keyword) => {
                keywords.add(keyword, parts[1..].to_vec());
            },
            _ => {
                debug!("Found Unknown keyword: {}", parts[0]);
                keywords.add(parts[0].to_uppercase().as_str(), parts[1..].to_vec());
            }
        }
    }
    keywords.print();
    let mesh_file = keywords.get_single_value("MESH_FILE").unwrap();
    let mesh_path = Path::new(file_path).parent().unwrap().join(mesh_file);
    let mesh_path_str = mesh_path.to_str().unwrap();
    let mesh: Mesh = read_file(mesh_path_str);
    // mesh.print_info();
    let nodes = mesh.convert_to_nodes();
    let elements = mesh.convert_to_elements();
    let dofs = keywords.get_single_float("DOF").unwrap_or(3.0) as usize;
    let mut simulation = Simulation::from_arrays(nodes, elements, dofs);
    

    //handle FIXED_BC
    let fixed_bc_values = keywords.get_keywords("FIXED_BC"); //FIXED_BC {group_name} 0.0 _ 0.1
    for fixed_bc_value in fixed_bc_values {
        let values: Vec<String> = fixed_bc_value.values.iter().map(|s| s.to_string()).collect();
        let node_group_name = values[0].clone();
        let node_ids = mesh.get_nodes_in_group(&node_group_name);
        let fixed_values: Vec<Option<f64>> = fixed_bc_value.values[1..].iter().map(|s| {
            if s == "_" {
                None
            } else {
                Some(s.parse::<f64>().unwrap())
            }
        }).collect();
        let fixed_condition: FixedCondition = FixedCondition::new(node_ids, fixed_values);
        // println!("Fixed condition: {:?}", fixed_condition);
        simulation.add_boundary_condition(Box::new(fixed_condition));
    }

    let load_bc_values = keywords.get_keywords("LOAD_BC");
    for load_bc_value in load_bc_values {
        let values: Vec<String> = load_bc_value.values.iter().map(|s| s.to_string()).collect();
        let node_group_name = values[0].clone();
        let node_ids = mesh.get_nodes_in_group(&node_group_name);
        let load_values: Vec<f64> = load_bc_value.values[1..].iter().map(|s| s.parse::<f64>().unwrap()).collect();
        let load_condition: LoadCondition = LoadCondition::new_from_vec(node_ids, load_values);
        // println!("Load condition: {:?}", load_condition);
        simulation.add_boundary_condition(Box::new(load_condition));
    }
    simulation.set_keywords(keywords.clone()); //duplicate keywords for the simulation okay for now. Small and easily compressed
    debug!("Finished reading simulation file");
    Ok((keywords, vec![simulation], vec![mesh]))
}


