// src/node.rs
use serde::{Serialize, Deserialize};
extern crate nalgebra as na;

use na::Vector3;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Node {
    pub id: u32,
    pub position: Vector3<f64>,
    pub displacement: Vector3<f64>,
    #[serde(skip, default = "Vector3::zeros")]
    pub velocity: Vector3<f64>,
    #[serde(skip, default = "Vector3::zeros")]
    pub acceleration: Vector3<f64>,
    pub mass: f64,
}

/// Implementation of the `Node` struct.
impl Node {
    /// Constructor for creating a new node.
    /// 
    /// This function initializes a new `Node` object with the specified ID and coordinates.
    /// 
    /// # Arguments
    /// * `id`: The unique identifier for the node.
    /// * `x`: The x-coordinate of the node.
    /// * `y`: The y-coordinate of the node.
    // Constructor
    pub fn new(id: u32, x: f64, y: f64, z: f64) -> Self {
        Node {
            id,
            position: Vector3::new(x, y, z),
            displacement: Vector3::new(0.0, 0.0, 0.0),
            velocity: Vector3::new(0.0, 0.0, 0.0),
            acceleration: Vector3::new(0.0, 0.0, 0.0),
            mass: 1.0,
        }
    }

    /// Sets the displacement for the node
    /// 
    /// This function updates the displacement of the node with the provided values.
    /// 
    /// # Arguments
    /// * `ux`: The displacement in the x-direction.
    /// * `uy`: The displacement in the y-direction.
    /// * `uz`: The displacement in the z-direction.
    pub fn set_displacement(&mut self, ux: f64, uy: f64, uz: f64) {
        self.displacement = Vector3::new(ux, uy, uz);
    }

    pub fn set_mass(&mut self, mass: f64) {
        self.mass = mass;
    }

    /// This function calculates the norm of the displacement vector to get the total displacement.
    /// 
    /// # Returns
    /// The total displacement as a scalar value.
    pub fn total_displacement(&self) -> f64 {
        self.displacement.norm()
    }
}
