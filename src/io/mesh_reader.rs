// written with gmsh with an export to .inp in mind
use std::fs::{read_to_string, canonicalize};
use std::collections::HashMap;
use crate::mesh::{ElementGroup, MeshAssembly, MeshElement, MeshNode, NodeGroup};
use log::{info, debug, trace, error};
use std::io::{BufReader, Read, BufRead};
use std::path::Path;
use std::fs::File;
use xz2::read::XzDecoder;

fn get_reader(filename: &str) -> BufReader<Box<dyn Read>>{
    let file_extension = Path::new(filename).extension().expect("Issue parsing file extension").to_str().unwrap();
    let file = File::open(filename).unwrap();
    let reader: Box<dyn Read> = match file_extension {
        "inp" => Box::new(file),
        "xz" => Box::new(XzDecoder::new(file)),
        _ => panic!("Unsupported file extension: {}", file_extension),
    };
    BufReader::new(reader)
}

//naive method but okay for the size of meshes
fn read_lines(filename: &str) -> Vec<String> {
    let abs_path = canonicalize(filename).unwrap();
    let abs_path = abs_path.to_str().unwrap();
    let mut result = Vec::new();
    let reader = get_reader(filename);
    for line in reader.lines() {
        result.push(line.unwrap());
    }
    result
}

/// Reads a mesh from a file.
/// 
/// This function reads a mesh from a file and returns a new `Mesh` instance.
/// It supports reading from both serialized mesh files (with .mesh extension) and INP files.
/// 
/// # Arguments
/// * `filename`: The path to the file containing the mesh data. .mesh or .inp
/// 
/// # Returns
/// A new `Mesh` instance.
pub fn read_file(filename: &str) -> MeshAssembly {

    if filename.contains(".mesh") {
        info!("Reading mesh from serialized file: {}", filename);
        return MeshAssembly::load(filename);
    }

    //panic of .inp or empty
    if filename.is_empty() { 
        error!("No mesh file specified");
        panic!("No mesh file specified"); 
    }

    info!("Reading mesh from file: {}", filename);

    let mut mesh = MeshAssembly::empty();
    let lines = read_lines(filename);

    let mut current_block = -1; //0 = *NODE, 1 = *ELEMENT, 2 = *ELSET, 3 = *NSET
    let headings = ["*NODE", "*ELEMENT", "*ELSET", "*NSET"];
    let mut params: HashMap<String, String> = HashMap::new();
    for line in lines {
        //check if a heading is contained in the line, then set current_block to the index of that heading
        if line.contains('*') {
            for (i, heading) in headings.iter().enumerate() {
                if line.contains(heading) {
                    current_block = i as i32;
                    //handle params
                    params = HashMap::new();
                    let line_parts = line.split(',').map(|s| s.trim()).collect::<Vec<&str>>();
                    for part in line_parts {
                        if part.contains('*') { continue }
                        let param_parts = part.split('=').map(|s| s.trim()).collect::<Vec<&str>>();
                        if param_parts.len() != 2 { continue }
                        params.insert(param_parts[0].to_string(), param_parts[1].to_string());
                    }
                    break;
                }else{
                    current_block = -1;
                }
            }
            continue;
        }
        
        match current_block {
            0 => {
                let line_parts = line.split(',').map(|s| s.trim()).collect::<Vec<&str>>();
                if line_parts.len() != 4 {
                    println!("Invalid line in *NODE block: {}", line);
                    continue;
                }
                let id: usize = line_parts[0].parse::<usize>().unwrap() - 1;
                let x = line_parts[1].parse::<f64>().unwrap();
                let y = line_parts[2].parse::<f64>().unwrap();
                let z = line_parts[3].parse::<f64>().unwrap();
                let node = MeshNode {
                    coordinates: vec![x, y, z],
                    id
                };
                mesh.nodes.insert(id, node);
            }
            1 => {
                let line_parts = line.split(',').map(|s| s.trim()).collect::<Vec<&str>>();
                let n = line_parts.len();
                let el_type = params.get("type").unwrap();
                let el_set = params.get("ELSET").unwrap();
                if n < 3 {
                    println!("Invalid line in *ELEMENT block: {}", line);
                    continue;
                }
                let id: usize = line_parts[0].parse::<usize>().unwrap() - 1;
                let mut connectivity = Vec::new();
                for i in 1..n {
                    let node_id: usize = line_parts[i].parse::<usize>().unwrap() - 1;
                    connectivity.push(node_id);
                }
                let element = MeshElement {
                    connectivity,
                    name: el_set.to_string(),
                    el_type: el_type.to_string(),
                    id
                };
                mesh.elements.insert(id, element);
                //add element to element group
                if mesh.element_groups.contains_key(el_set) {
                    let group = mesh.element_groups.get_mut(el_set).unwrap();
                    group.elements.push(id);
                } else {
                    let group = ElementGroup {
                        elements: vec![id],
                        name: el_set.to_string(),
                        el_type: el_type.to_string()
                    };
                    mesh.element_groups.insert(el_set.to_string(), group);
                }
            }
            2 => {
                //handle *ELSET block
                let el_set = params.get("ELSET").unwrap();
                let el_nums_parts = line.split(',').map(|s| s.trim());
                let el_nums = el_nums_parts.filter(|s| !s.is_empty()).map(|s| s.parse::<usize>().unwrap()).collect::<Vec<usize>>();
                //index from zero
                let el_nums = el_nums.iter().map(|n| n - 1).collect::<Vec<usize>>();
                if mesh.element_groups.contains_key(el_set) {
                    let group = mesh.element_groups.get_mut(el_set).unwrap();
                    group.elements.extend(el_nums);
                } else {
                    let group = ElementGroup {
                        elements: el_nums,
                        name: el_set.to_string(),
                        el_type: String::new()
                    };
                    mesh.element_groups.insert(el_set.to_string(), group);
                }
            }
            3 => {
                //handle *NSET block
                let n_set = params.get("NSET").unwrap();
                let node_nums_parts = line.split(',').map(|s| s.trim());
                let node_nums = node_nums_parts.filter(|s| !s.is_empty()).map(|s| s.parse::<usize>().unwrap()).collect::<Vec<usize>>();
                //index from zero
                let node_nums = node_nums.iter().map(|n| n - 1).collect::<Vec<usize>>();
                if mesh.node_groups.contains_key(n_set) {
                    let group = mesh.node_groups.get_mut(n_set).unwrap();
                    group.nodes.extend(node_nums);
                } else {
                    let group = NodeGroup {
                        nodes: node_nums,
                        name: n_set.to_string(),
                    };
                    mesh.node_groups.insert(n_set.to_string(), group);
                }
            }
            _ => continue,
        }
    }
    debug!("Finished reading mesh file");
    mesh
}


pub fn read_file_single_body(filename: &str) -> MeshAssembly {
    let mut mesh = read_file(filename);
    mesh.single_body();
    mesh
}

pub fn read_file_multiple_bodies(filename: &str, bodies: Vec<String>) -> MeshAssembly {
    let mut mesh = read_file(filename);
    mesh.multiple_bodies(bodies);
    mesh
}

