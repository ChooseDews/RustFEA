// src/node.rs

extern crate nalgebra as na;

use na::Vector3;

#[derive(Debug, Clone, Copy)]
pub struct Node {
    pub id: usize,
    pub position: Vector3<f64>,
    pub displacement: Vector3<f64>,
}

impl Node {
    // Constructor
    pub fn new(id: usize, x: f64, y: f64, z: f64) -> Self {
        Node {
            id,
            position: Vector3::new(x, y, z),
            displacement: Vector3::new(0.0, 0.0, 0.0),
        }
    }

    // Sets the displacement for the node
    pub fn set_displacement(&mut self, ux: f64, uy: f64, uz: f64) {
        self.displacement = Vector3::new(ux, uy, uz);
    }

    // Returns the total displacement (as a vector magnitude)
    pub fn total_displacement(&self) -> f64 {
        self.displacement.norm()
    }
}
