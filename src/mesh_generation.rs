// src/mesh_generation.rs
extern crate nalgebra as na;
use crate::elements::brick_element::BrickElement;
use crate::elements::base_element::Material;
use crate::node::Node;

pub fn generate_mesh(n_x: usize, n_y: usize, n_z: usize, s_x: f64, s_y: f64, s_z: f64) -> (Vec<Node>, Vec<BrickElement>) {
    let mut nodes = Vec::new();
    let mut elements = Vec::new();
    let dx = s_x / n_x as f64;
    let dy = s_y / n_y as f64;
    let dz = s_z / n_z as f64;

    // Generating node positions
    let mut id = 0;
    for k in 0..=n_z {
        for j in 0..=n_y {
            for i in 0..=n_x {
                nodes.push(Node::new(id, i as f64 * dx, j as f64 * dy, k as f64 * dz));
                id += 1;
            }
        }
    }
    // Generating element connectivity
    let mut elem_id = 0;
    for k in 0..n_z {
        for j in 0..n_y {
            for i in 0..n_x {
                let bottom_front_left = i + j * (n_x + 1) + k * (n_x + 1) * (n_y + 1) + 1;
                let bottom_front_right = bottom_front_left + 1;
                let bottom_back_left = bottom_front_left + n_x + 1;
                let bottom_back_right = bottom_back_left + 1;
                let top_front_left = bottom_front_left + (n_x + 1) * (n_y + 1);
                let top_front_right = top_front_left + 1;
                let top_back_left = top_front_left + n_x + 1;
                let top_back_right = top_back_left + 1;

                // Using our previously defined BrickElement constructor
                elements.push(BrickElement::new(
                    elem_id,
                    vec![
                        bottom_front_left-1,
                        bottom_front_right-1,
                        bottom_back_right-1,
                        bottom_back_left-1,
                        top_front_left-1,
                        top_front_right-1,
                        top_back_right-1,
                        top_back_left-1,
                    ],
                    Material::aluminum()
                ));
                elem_id += 1;
            }
        }
    }

    (nodes, elements)
}

