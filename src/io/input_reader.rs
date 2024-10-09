use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::collections::HashMap;
use crate::io::file;
use crate::simulation::Simulation;
use crate::bc::{BoundaryCondition, NormalContact, FixedCondition, LoadCondition};
use crate::io::mesh_reader::{read_file_single_body, read_file_multiple_bodies};
use std::error::Error;
use nalgebra::DVector;
use serde::{Serialize, Deserialize};
use crate::mesh::MeshAssembly;
use crate::io::project::Project;
use crate::utilities::Keywords;

use log::{info, debug, warn};

use super::project;

fn remove_comments(line: &str) -> String {
    match line.find('#') {
        Some(index) => line[..index].trim().to_string(),
        None => line.trim().to_string(),
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
pub fn read_simulation_file(file_path: &str) -> Result<Project, Box<dyn Error>> {
    if file_path.ends_with(".toml"){
        return Ok(read_toml_file(file_path));
    }
    Err("File format not supported".into())
}



const RESERVED_KEYS: [&str; 3] = ["SIM", "MATERIAL", "BOUNDARY_CONDITIONS"];
pub fn table_to_keywords(table: toml::Table) -> Keywords {
    let mut keywords = Keywords::new();
    fn recursive_table_to_keywords(table: toml::Table, keywords: &mut Keywords, parent_key: String) {
        for (key, value) in table.iter() {
            if RESERVED_KEYS.contains(&key.to_uppercase().as_str()) {
                continue;
            }
            let key_name: String = if parent_key.is_empty() {key.to_string()} else {format!("{}_{}", &parent_key, key)};
            if value.is_table() {
                let value_table = value.as_table().unwrap();
                recursive_table_to_keywords(value_table.clone(), keywords, key_name.clone());
            }
            else {
                keywords.add(&key_name, value.clone());
            }
        }
    }
    recursive_table_to_keywords(table, &mut keywords, "".to_string());
    keywords
}

pub fn format_path(file_path: &str, input_file_path: &str) -> String {
    let file_path = Path::new(file_path);
    let input_file_path = Path::new(input_file_path);
    let project_dir = input_file_path.parent().unwrap();
    let formatted_path = if file_path.is_absolute() {
        file_path.to_path_buf()
    } else {
        project_dir.join(file_path)
    };
    //return cononical path
    formatted_path.canonicalize().unwrap().to_str().unwrap().to_string()
}


pub fn read_toml_file(file_path: &str) -> Project {
    info!("Reading TOML simulation input file: {}", file_path);
    let file = File::open(file_path).unwrap();
    let mut reader = BufReader::new(file);
    let mut contents = String::new();
    reader.read_to_string(&mut contents).unwrap();
    let value: toml::Table = toml::from_str(&contents).unwrap();
    let mut project_keywords = table_to_keywords(value.clone());
    project_keywords.print();
    let sim_array = value.get("sim").unwrap().as_array().unwrap();
    let mut simulations: Vec<Simulation> = Vec::new();
    for sim in sim_array.iter() {
        let sim_table = sim.as_table().unwrap();
        let mut simulation_keywords = table_to_keywords(sim_table.clone());
        simulation_keywords.print();
        //handle single or multiple meshes
        let mut mesh_list: Vec<MeshAssembly> = Vec::new();
        //check if meshes is a array in the sim table
        if sim_table.get("meshes").is_some() {
            let meshes = sim_table.get("meshes").unwrap().as_array().unwrap();
            for mesh in meshes.iter() {
                let mesh_table = mesh.as_table().unwrap();
                let mesh_file = mesh_table.get("file").unwrap().as_str().unwrap();
                let mesh_path = format_path(&mesh_file, file_path);

                //get array of bodies corresponding to the mesh node_sets
                let volumes = {
                    if mesh_table.get("volumes").is_some() {    
                        mesh_table.get("volumes").unwrap().as_array().unwrap().clone()
                    } else {
                        vec![]
                    }
                };
                let mut bodies: Vec<String> = Vec::new();
                for volume in volumes.iter() {
                    let volume_name = volume.as_str().unwrap();
                    bodies.push(volume_name.to_string());
                }
                let mut mesh = if bodies.is_empty() {  
                    let mesh_name = {
                        if mesh_table.get("name").is_some() {
                            mesh_table.get("name").unwrap().as_str().unwrap().to_string()
                        } else {
                            let name = Path::new(&mesh_file).file_name().unwrap().to_str().unwrap();
                            let name = name.split('.').next().unwrap();
                            format!("{}_body", name)
                        }
                    };
                    let mut mesh = read_file_single_body(&mesh_path);
                    mesh.set_body_name(mesh_name.as_str());
                    mesh
                }else{
                    read_file_multiple_bodies(&mesh_path, bodies)
                };
                mesh_list.push(mesh);
            }
        }else{
            let mesh_file = simulation_keywords.get_string("MESH").unwrap();
            let mesh_path = format_path(&mesh_file, file_path);
            let mesh = read_file_single_body(&mesh_path);
            mesh_list.push(mesh);
        }

        let mut mesh = mesh_list.remove(0);
        for m in mesh_list {
            mesh += m;
        }

        mesh.print_info();





        let mesh_nodes = mesh.convert_to_nodes();
        let mesh_elements = mesh.convert_to_elements();
        let dofs = simulation_keywords.get_int("DOF").unwrap_or(3) as usize;
        let mut simulation = Simulation::from_mesh(mesh, dofs);
        simulation.set_keywords(simulation_keywords);
        let boundary_conditions = sim_table.get("boundary_conditions").unwrap().as_array().unwrap();
        for bc in boundary_conditions.iter() {
            let bc_table = bc.as_table().unwrap();
            let bc_type = bc_table.get("type").unwrap().as_str().unwrap();
            let bc_name = bc_table.get("name").unwrap().as_str().unwrap();
            let emtpy_values = toml::Value::Array(vec![]);
            let bc_values = bc_table.get("values").unwrap_or(&emtpy_values).as_array().unwrap(); //contains a float or false thus option is none
            let bc_values_vec: Vec<Option<f64>> = bc_values.iter().map(|v| {
                if v.is_bool(){
                    if v.as_bool().unwrap() {
                        Some(0.0)
                    } else {
                        None
                    }
                } else if v.is_float() {
                    Some(v.as_float().unwrap())
                } else if v.is_integer() {
                    Some(v.as_integer().unwrap() as f64)
                } else {
                    None
                }
            }).collect();
            let bc_node_ids: Vec<usize> = simulation.mesh.get_nodes_in_group(&bc_name);
            match bc_type {
                "fixed" => {
                    let bc = FixedCondition::new(bc_node_ids, bc_values_vec);
                    simulation.add_boundary_condition(Box::new(bc));
                }
                "load" => {
                    let force = DVector::from(bc_values_vec.iter().map(|v| v.unwrap()).collect::<Vec<f64>>());
                    let bc = LoadCondition::new(bc_node_ids, force);
                    simulation.add_boundary_condition(Box::new(bc));
                }
                "contact" => {
                    let secondary_surface = bc_table.get("secondary").unwrap().as_str().unwrap();
                    let mut bc = NormalContact::new(bc_name.to_string(), secondary_surface.to_string());
                    debug!("Setting contact surfaces nodes for contact condition: {} total nodes: {}", bc_name, bc_node_ids.len());
                    let secondary_surface_nodes = simulation.mesh.get_nodes_in_group(secondary_surface);
                    debug!("Secondary surface nodes: {}", secondary_surface_nodes.len());
                    bc.set_contact_surfaces_nodes(bc_node_ids, secondary_surface_nodes);
                    let primary_surface_elements = simulation.mesh.get_elements_in_group(bc_name);
                    let secondary_surface_elements = simulation.mesh.get_elements_in_group(secondary_surface);
                    debug!("Primary surface elements: {}", primary_surface_elements.len());
                    bc.set_contact_surfaces_elements(primary_surface_elements, secondary_surface_elements);
                    simulation.add_boundary_condition(Box::new(bc));
                }
                _ => {
                    warn!("Unknown boundary condition type: {}", bc_type);
                }
            }
        }
        simulations.push(simulation);
    }
    Project::new(simulations, project_keywords)
}

