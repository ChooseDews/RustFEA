// src/elements/brick20_element.rs
//! 20-node quadratic hexahedral (serendipity) element (C3D20)
//!
//! Node numbering follows Abaqus/GMSH convention:
//!
//! ```text
//!        7----14----6
//!       /|         /|
//!     15 |       13 |
//!     /  19      /  18
//!    4----12----5   |
//!    |   |      |   |
//!    |   3---10-|---2
//!   16  /      17  /
//!    | 11       | 9
//!    |/         |/
//!    0----8-----1
//!
//! Corner nodes: 0-7 (at ξ,η,ζ = ±1)
//! Mid-edge nodes: 8-19
//!   8:  (0,-1,-1)  between 0-1
//!   9:  (1,0,-1)   between 1-2
//!  10:  (0,1,-1)   between 2-3
//!  11:  (-1,0,-1)  between 3-0
//!  12:  (0,-1,1)   between 4-5
//!  13:  (1,0,1)    between 5-6
//!  14:  (0,1,1)    between 6-7
//!  15:  (-1,0,1)   between 7-4
//!  16:  (-1,-1,0)  between 0-4
//!  17:  (1,-1,0)   between 1-5
//!  18:  (1,1,0)    between 2-6
//!  19:  (-1,1,0)   between 3-7
//! ```

use crate::{simulation::Simulation, utilities::check_for_nans};
use super::base_element::{BaseElement, Material, ElementFields, ElementType};
use nalgebra as na;
use na::{DMatrix, DVector, Matrix3, Vector6};
use crate::utilities::compute_von_mises;
use serde::{Serialize, Deserialize};
use log::{debug, trace};

/// Number of nodes in C3D20 element
const NUM_NODES: usize = 20;
/// Number of DOFs (3 per node)
const NUM_DOFS: usize = NUM_NODES * 3;

#[derive(Serialize, Deserialize, Debug)]
pub struct Brick20Element {
    id: usize,
    connectivity: Vec<usize>,
    material: Material,
    #[serde(skip, default = "default_deformation_gradient")]
    deformation_gradient: DMatrix<f64>,
    #[serde(skip, default = "default_stiffness")]
    stiffness: DMatrix<f64>,
    #[serde(skip, default = "default_mass")]
    mass: DMatrix<f64>,
    #[serde(skip, default = "Vec::new")]
    lumped_mass: Vec<f64>,
    active: bool,
    #[serde(skip, default = "nodal_positions_default")]
    nodal_positions: Option<DMatrix<f64>>,
    #[serde(skip, default = "default_internal_force")]
    internal_force: DVector<f64>,
}

fn default_internal_force() -> DVector<f64> {
    DVector::zeros(NUM_DOFS)
}

fn nodal_positions_default() -> Option<DMatrix<f64>> {
    None
}

fn default_deformation_gradient() -> DMatrix<f64> {
    DMatrix::identity(3, 3)
}

fn default_mass() -> DMatrix<f64> {
    DMatrix::zeros(NUM_NODES, NUM_NODES)
}

fn default_stiffness() -> DMatrix<f64> {
    DMatrix::zeros(NUM_DOFS, NUM_DOFS)
}

impl Brick20Element {
    pub fn new(id: usize, connectivity: Vec<usize>, material: Material) -> Self {
        assert_eq!(
            connectivity.len(),
            NUM_NODES,
            "20 nodes required for a C3D20 element"
        );
        Brick20Element {
            id,
            connectivity,
            material,
            deformation_gradient: DMatrix::identity(3, 3),
            stiffness: default_stiffness(),
            mass: default_mass(),
            lumped_mass: Vec::new(),
            active: true,
            nodal_positions: None,
            internal_force: default_internal_force(),
        }
    }

    /// Natural coordinates of the 20 nodes
    fn get_node_coordinates() -> &'static [(f64, f64, f64)] {
        static COORDS: [(f64, f64, f64); 20] = [
            // Corner nodes 0-7
            (-1.0, -1.0, -1.0), // 0
            (1.0, -1.0, -1.0),  // 1
            (1.0, 1.0, -1.0),   // 2
            (-1.0, 1.0, -1.0),  // 3
            (-1.0, -1.0, 1.0),  // 4
            (1.0, -1.0, 1.0),   // 5
            (1.0, 1.0, 1.0),    // 6
            (-1.0, 1.0, 1.0),   // 7
            // Mid-edge nodes 8-19
            (0.0, -1.0, -1.0),  // 8  (between 0-1)
            (1.0, 0.0, -1.0),   // 9  (between 1-2)
            (0.0, 1.0, -1.0),   // 10 (between 2-3)
            (-1.0, 0.0, -1.0),  // 11 (between 3-0)
            (0.0, -1.0, 1.0),   // 12 (between 4-5)
            (1.0, 0.0, 1.0),    // 13 (between 5-6)
            (0.0, 1.0, 1.0),    // 14 (between 6-7)
            (-1.0, 0.0, 1.0),   // 15 (between 7-4)
            (-1.0, -1.0, 0.0),  // 16 (between 0-4)
            (1.0, -1.0, 0.0),   // 17 (between 1-5)
            (1.0, 1.0, 0.0),    // 18 (between 2-6)
            (-1.0, 1.0, 0.0),   // 19 (between 3-7)
        ];
        &COORDS
    }

    /// 3x3x3 Gauss quadrature points for 20-node element
    fn get_gauss_points() -> &'static [(f64, f64, f64, f64)] {
        // 3x3x3 Gauss quadrature (27 points)
        static GP: f64 = 0.7745966692414834; // sqrt(3/5)
        static W1: f64 = 0.5555555555555556; // 5/9
        static W2: f64 = 0.8888888888888889; // 8/9
        
        static GAUSS_POINTS: [(f64, f64, f64, f64); 27] = [
            (-GP, -GP, -GP, W1 * W1 * W1),
            (0.0, -GP, -GP, W2 * W1 * W1),
            (GP, -GP, -GP, W1 * W1 * W1),
            (-GP, 0.0, -GP, W1 * W2 * W1),
            (0.0, 0.0, -GP, W2 * W2 * W1),
            (GP, 0.0, -GP, W1 * W2 * W1),
            (-GP, GP, -GP, W1 * W1 * W1),
            (0.0, GP, -GP, W2 * W1 * W1),
            (GP, GP, -GP, W1 * W1 * W1),
            
            (-GP, -GP, 0.0, W1 * W1 * W2),
            (0.0, -GP, 0.0, W2 * W1 * W2),
            (GP, -GP, 0.0, W1 * W1 * W2),
            (-GP, 0.0, 0.0, W1 * W2 * W2),
            (0.0, 0.0, 0.0, W2 * W2 * W2),
            (GP, 0.0, 0.0, W1 * W2 * W2),
            (-GP, GP, 0.0, W1 * W1 * W2),
            (0.0, GP, 0.0, W2 * W1 * W2),
            (GP, GP, 0.0, W1 * W1 * W2),
            
            (-GP, -GP, GP, W1 * W1 * W1),
            (0.0, -GP, GP, W2 * W1 * W1),
            (GP, -GP, GP, W1 * W1 * W1),
            (-GP, 0.0, GP, W1 * W2 * W1),
            (0.0, 0.0, GP, W2 * W2 * W1),
            (GP, 0.0, GP, W1 * W2 * W1),
            (-GP, GP, GP, W1 * W1 * W1),
            (0.0, GP, GP, W2 * W1 * W1),
            (GP, GP, GP, W1 * W1 * W1),
        ];
        &GAUSS_POINTS
    }

    /// Compute shape functions for 20-node serendipity element
    /// 
    /// For corner nodes (i = 0..7):
    ///   N_i = (1/8)(1 + ξ_i*ξ)(1 + η_i*η)(1 + ζ_i*ζ)(ξ_i*ξ + η_i*η + ζ_i*ζ - 2)
    /// 
    /// For mid-edge nodes (i = 8..19):
    ///   - Nodes 8,10,12,14 (ξ_i = 0): N_i = (1/4)(1 - ξ²)(1 + η_i*η)(1 + ζ_i*ζ)
    ///   - Nodes 9,11,13,15 (η_i = 0): N_i = (1/4)(1 + ξ_i*ξ)(1 - η²)(1 + ζ_i*ζ)
    ///   - Nodes 16,17,18,19 (ζ_i = 0): N_i = (1/4)(1 + ξ_i*ξ)(1 + η_i*η)(1 - ζ²)
    fn get_shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> DVector<f64> {
        let mut n = DVector::zeros(NUM_NODES);
        let node_coords = Self::get_node_coordinates();
        
        // Precompute terms
        let xi2 = xi * xi;
        let eta2 = eta * eta;
        let zeta2 = zeta * zeta;
        
        // Corner nodes (0-7)
        for i in 0..8 {
            let (xi_i, eta_i, zeta_i) = node_coords[i];
            let xi_term = 1.0 + xi_i * xi;
            let eta_term = 1.0 + eta_i * eta;
            let zeta_term = 1.0 + zeta_i * zeta;
            n[i] = 0.125 * xi_term * eta_term * zeta_term * (xi_i * xi + eta_i * eta + zeta_i * zeta - 2.0);
        }
        
        // Mid-edge nodes on bottom face (z = -1): 8, 9, 10, 11
        // Node 8: between 0-1, at (0, -1, -1)
        n[8] = 0.25 * (1.0 - xi2) * (1.0 - eta) * (1.0 - zeta);
        // Node 9: between 1-2, at (1, 0, -1)
        n[9] = 0.25 * (1.0 + xi) * (1.0 - eta2) * (1.0 - zeta);
        // Node 10: between 2-3, at (0, 1, -1)
        n[10] = 0.25 * (1.0 - xi2) * (1.0 + eta) * (1.0 - zeta);
        // Node 11: between 3-0, at (-1, 0, -1)
        n[11] = 0.25 * (1.0 - xi) * (1.0 - eta2) * (1.0 - zeta);
        
        // Mid-edge nodes on top face (z = +1): 12, 13, 14, 15
        // Node 12: between 4-5, at (0, -1, 1)
        n[12] = 0.25 * (1.0 - xi2) * (1.0 - eta) * (1.0 + zeta);
        // Node 13: between 5-6, at (1, 0, 1)
        n[13] = 0.25 * (1.0 + xi) * (1.0 - eta2) * (1.0 + zeta);
        // Node 14: between 6-7, at (0, 1, 1)
        n[14] = 0.25 * (1.0 - xi2) * (1.0 + eta) * (1.0 + zeta);
        // Node 15: between 7-4, at (-1, 0, 1)
        n[15] = 0.25 * (1.0 - xi) * (1.0 - eta2) * (1.0 + zeta);
        
        // Mid-edge nodes on vertical edges: 16, 17, 18, 19
        // Node 16: between 0-4, at (-1, -1, 0)
        n[16] = 0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta2);
        // Node 17: between 1-5, at (1, -1, 0)
        n[17] = 0.25 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta2);
        // Node 18: between 2-6, at (1, 1, 0)
        n[18] = 0.25 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta2);
        // Node 19: between 3-7, at (-1, 1, 0)
        n[19] = 0.25 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta2);
        
        n
    }

    /// Compute shape function derivatives with respect to natural coordinates (ξ, η, ζ)
    /// Returns a 20x3 matrix: [∂N_i/∂ξ, ∂N_i/∂η, ∂N_i/∂ζ]
    fn get_shape_derivatives_local(&self, xi: f64, eta: f64, zeta: f64) -> DMatrix<f64> {
        let mut dn = DMatrix::zeros(NUM_NODES, 3);
        let node_coords = Self::get_node_coordinates();
        
        // Precompute terms
        let xi2 = xi * xi;
        let eta2 = eta * eta;
        let zeta2 = zeta * zeta;
        
        // Corner nodes (0-7)
        // N_i = (1/8)(1 + ξ_i*ξ)(1 + η_i*η)(1 + ζ_i*ζ)(ξ_i*ξ + η_i*η + ζ_i*ζ - 2)
        for i in 0..8 {
            let (xi_i, eta_i, zeta_i) = node_coords[i];
            
            let a = 1.0 + xi_i * xi;
            let b = 1.0 + eta_i * eta;
            let c = 1.0 + zeta_i * zeta;
            let d = xi_i * xi + eta_i * eta + zeta_i * zeta - 2.0;
            
            // ∂N_i/∂ξ = (1/8) * [ξ_i * b * c * d + a * b * c * ξ_i]
            //         = (1/8) * ξ_i * b * c * (d + a)
            //         = (1/8) * ξ_i * b * c * (2*ξ_i*ξ + η_i*η + ζ_i*ζ - 1)
            dn[(i, 0)] = 0.125 * xi_i * b * c * (2.0 * xi_i * xi + eta_i * eta + zeta_i * zeta - 1.0);
            
            // ∂N_i/∂η
            dn[(i, 1)] = 0.125 * eta_i * a * c * (xi_i * xi + 2.0 * eta_i * eta + zeta_i * zeta - 1.0);
            
            // ∂N_i/∂ζ
            dn[(i, 2)] = 0.125 * zeta_i * a * b * (xi_i * xi + eta_i * eta + 2.0 * zeta_i * zeta - 1.0);
        }
        
        // Mid-edge nodes on bottom face (z = -1): 8, 9, 10, 11
        // Node 8: N_8 = (1/4)(1 - ξ²)(1 - η)(1 - ζ)
        dn[(8, 0)] = -0.5 * xi * (1.0 - eta) * (1.0 - zeta);
        dn[(8, 1)] = -0.25 * (1.0 - xi2) * (1.0 - zeta);
        dn[(8, 2)] = -0.25 * (1.0 - xi2) * (1.0 - eta);
        
        // Node 9: N_9 = (1/4)(1 + ξ)(1 - η²)(1 - ζ)
        dn[(9, 0)] = 0.25 * (1.0 - eta2) * (1.0 - zeta);
        dn[(9, 1)] = -0.5 * eta * (1.0 + xi) * (1.0 - zeta);
        dn[(9, 2)] = -0.25 * (1.0 + xi) * (1.0 - eta2);
        
        // Node 10: N_10 = (1/4)(1 - ξ²)(1 + η)(1 - ζ)
        dn[(10, 0)] = -0.5 * xi * (1.0 + eta) * (1.0 - zeta);
        dn[(10, 1)] = 0.25 * (1.0 - xi2) * (1.0 - zeta);
        dn[(10, 2)] = -0.25 * (1.0 - xi2) * (1.0 + eta);
        
        // Node 11: N_11 = (1/4)(1 - ξ)(1 - η²)(1 - ζ)
        dn[(11, 0)] = -0.25 * (1.0 - eta2) * (1.0 - zeta);
        dn[(11, 1)] = -0.5 * eta * (1.0 - xi) * (1.0 - zeta);
        dn[(11, 2)] = -0.25 * (1.0 - xi) * (1.0 - eta2);
        
        // Mid-edge nodes on top face (z = +1): 12, 13, 14, 15
        // Node 12: N_12 = (1/4)(1 - ξ²)(1 - η)(1 + ζ)
        dn[(12, 0)] = -0.5 * xi * (1.0 - eta) * (1.0 + zeta);
        dn[(12, 1)] = -0.25 * (1.0 - xi2) * (1.0 + zeta);
        dn[(12, 2)] = 0.25 * (1.0 - xi2) * (1.0 - eta);
        
        // Node 13: N_13 = (1/4)(1 + ξ)(1 - η²)(1 + ζ)
        dn[(13, 0)] = 0.25 * (1.0 - eta2) * (1.0 + zeta);
        dn[(13, 1)] = -0.5 * eta * (1.0 + xi) * (1.0 + zeta);
        dn[(13, 2)] = 0.25 * (1.0 + xi) * (1.0 - eta2);
        
        // Node 14: N_14 = (1/4)(1 - ξ²)(1 + η)(1 + ζ)
        dn[(14, 0)] = -0.5 * xi * (1.0 + eta) * (1.0 + zeta);
        dn[(14, 1)] = 0.25 * (1.0 - xi2) * (1.0 + zeta);
        dn[(14, 2)] = 0.25 * (1.0 - xi2) * (1.0 + eta);
        
        // Node 15: N_15 = (1/4)(1 - ξ)(1 - η²)(1 + ζ)
        dn[(15, 0)] = -0.25 * (1.0 - eta2) * (1.0 + zeta);
        dn[(15, 1)] = -0.5 * eta * (1.0 - xi) * (1.0 + zeta);
        dn[(15, 2)] = 0.25 * (1.0 - xi) * (1.0 - eta2);
        
        // Mid-edge nodes on vertical edges: 16, 17, 18, 19
        // Node 16: N_16 = (1/4)(1 - ξ)(1 - η)(1 - ζ²)
        dn[(16, 0)] = -0.25 * (1.0 - eta) * (1.0 - zeta2);
        dn[(16, 1)] = -0.25 * (1.0 - xi) * (1.0 - zeta2);
        dn[(16, 2)] = -0.5 * zeta * (1.0 - xi) * (1.0 - eta);
        
        // Node 17: N_17 = (1/4)(1 + ξ)(1 - η)(1 - ζ²)
        dn[(17, 0)] = 0.25 * (1.0 - eta) * (1.0 - zeta2);
        dn[(17, 1)] = -0.25 * (1.0 + xi) * (1.0 - zeta2);
        dn[(17, 2)] = -0.5 * zeta * (1.0 + xi) * (1.0 - eta);
        
        // Node 18: N_18 = (1/4)(1 + ξ)(1 + η)(1 - ζ²)
        dn[(18, 0)] = 0.25 * (1.0 + eta) * (1.0 - zeta2);
        dn[(18, 1)] = 0.25 * (1.0 + xi) * (1.0 - zeta2);
        dn[(18, 2)] = -0.5 * zeta * (1.0 + xi) * (1.0 + eta);
        
        // Node 19: N_19 = (1/4)(1 - ξ)(1 + η)(1 - ζ²)
        dn[(19, 0)] = -0.25 * (1.0 + eta) * (1.0 - zeta2);
        dn[(19, 1)] = 0.25 * (1.0 - xi) * (1.0 - zeta2);
        dn[(19, 2)] = -0.5 * zeta * (1.0 - xi) * (1.0 + eta);
        
        dn
    }

    fn get_x_local(&self, simulation: &Simulation) -> &DMatrix<f64> {
        self.nodal_positions.as_ref().expect("Nodal positions not initialized")
    }

    fn compute_x_local(&self, simulation: &Simulation) -> DMatrix<f64> {
        let mut x = DMatrix::zeros(NUM_NODES, 3);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            x[(i, 0)] = node.position[0];
            x[(i, 1)] = node.position[1];
            x[(i, 2)] = node.position[2];
        }
        x
    }

    fn get_u_local(&self, simulation: &Simulation) -> DVector<f64> {
        let mut u = DVector::zeros(NUM_DOFS);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            u[3 * i] = node.displacement[0];
            u[3 * i + 1] = node.displacement[1];
            u[3 * i + 2] = node.displacement[2];
        }
        u
    }

    fn get_u_local_from_displacement(&self, displacement: &DVector<f64>) -> DVector<f64> {
        let mut u = DVector::zeros(NUM_DOFS);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            u[3 * i] = displacement[3 * node_id];
            u[3 * i + 1] = displacement[3 * node_id + 1];
            u[3 * i + 2] = displacement[3 * node_id + 2];
        }
        u
    }

    /// Compute Jacobian matrix J = dN/dξ * X
    fn compute_jacobian_matrix(&self, x: &DMatrix<f64>, d_n: &DMatrix<f64>) -> Matrix3<f64> {
        let j_mat = d_n.transpose() * x;
        Matrix3::new(
            j_mat[(0, 0)], j_mat[(0, 1)], j_mat[(0, 2)],
            j_mat[(1, 0)], j_mat[(1, 1)], j_mat[(1, 2)],
            j_mat[(2, 0)], j_mat[(2, 1)], j_mat[(2, 2)],
        )
    }

    /// Compute strain-displacement matrix B (6 x 60)
    /// B relates nodal displacements to strains: ε = B * u
    fn compute_b(&self, x: &DMatrix<f64>, j: &Matrix3<f64>, d_n: &DMatrix<f64>) -> DMatrix<f64> {
        let mut b = DMatrix::zeros(6, NUM_DOFS);
        let j_inv = j.try_inverse().expect("Jacobian matrix is singular");
        
        for i in 0..NUM_NODES {
            // Compute global shape function derivatives: ∂N_i/∂x_j = (J^-1)_jk * ∂N_i/∂ξ_k
            // d_n is 20x3 where row i is [∂N_i/∂ξ, ∂N_i/∂η, ∂N_i/∂ζ]
            let dn_dxi = na::Vector3::new(d_n[(i, 0)], d_n[(i, 1)], d_n[(i, 2)]);
            let n_i = j_inv * dn_dxi;  // [∂N_i/∂x, ∂N_i/∂y, ∂N_i/∂z]
            
            // Build B_i block (6x3) for node i
            // Strain vector: [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]
            b[(0, 3 * i)] = n_i[0];     // ε_xx = ∂u/∂x
            b[(1, 3 * i + 1)] = n_i[1]; // ε_yy = ∂v/∂y
            b[(2, 3 * i + 2)] = n_i[2]; // ε_zz = ∂w/∂z
            b[(3, 3 * i)] = n_i[1];     // γ_xy: ∂u/∂y
            b[(3, 3 * i + 1)] = n_i[0]; // γ_xy: ∂v/∂x
            b[(4, 3 * i + 1)] = n_i[2]; // γ_yz: ∂v/∂z
            b[(4, 3 * i + 2)] = n_i[1]; // γ_yz: ∂w/∂y
            b[(5, 3 * i)] = n_i[2];     // γ_xz: ∂u/∂z
            b[(5, 3 * i + 2)] = n_i[0]; // γ_xz: ∂w/∂x
        }
        b
    }

    fn compute_stress(&self, x: &DMatrix<f64>, u: &DVector<f64>, d_n: &DMatrix<f64>) -> DVector<f64> {
        let j = self.compute_jacobian_matrix(x, d_n);
        let b = self.compute_b(x, &j, d_n);
        let c = self.material.get_3d_matrix();
        DVector::from_column_slice((c * b * u).as_slice())
    }

    fn compute_strain(&self, x: &DMatrix<f64>, u: &DVector<f64>, d_n: &DMatrix<f64>) -> DVector<f64> {
        let j = self.compute_jacobian_matrix(x, d_n);
        let b = self.compute_b(x, &j, d_n);
        b * u
    }
}

#[typetag::serde]
impl BaseElement for Brick20Element {
    fn get_id(&self) -> usize {
        self.id
    }

    fn initialize(&mut self, simulation: &Simulation) {
        self.nodal_positions = Some(self.compute_x_local(simulation));
    }

    fn get_connectivity(&self) -> &Vec<usize> {
        &self.connectivity
    }

    fn get_material(&self) -> &Material {
        &self.material
    }

    fn get_deformation_gradient(&self) -> &DMatrix<f64> {
        &self.deformation_gradient
    }

    fn get_global_position(&self, n: &DVector<f64>, simulation: &Simulation) -> na::Vector3<f64> {
        let mut global_position = na::Vector3::zeros();
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            global_position += n[i] * node.position;
        }
        global_position
    }

    fn get_x(&self, simulation: &Simulation) -> DMatrix<f64> {
        let mut x_global = DMatrix::zeros(3, NUM_NODES);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            x_global[(0, i)] = node.position[0];
            x_global[(1, i)] = node.position[1];
            x_global[(2, i)] = node.position[2];
        }
        x_global
    }

    fn get_u(&self, simulation: &Simulation) -> DVector<f64> {
        self.get_u_local(simulation)
    }

    fn get_shape_derivatives(&self, xi: f64, eta: f64, zeta: f64) -> DMatrix<f64> {
        self.get_shape_derivatives_local(xi, eta, zeta)
    }

    fn get_b(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DMatrix<f64> {
        let x = self.get_x_local(simulation);
        let d_n = self.get_shape_derivatives_local(xi, eta, zeta);
        let j = self.compute_jacobian_matrix(x, &d_n);
        self.compute_b(x, &j, &d_n)
    }

    fn compute_stiffness(&mut self, simulation: &Simulation) {
        trace!("Computing stiffness matrix for C3D20 element");
        let mut k = DMatrix::zeros(NUM_DOFS, NUM_DOFS);
        let gauss_points = Self::get_gauss_points();
        let x = self.get_x_local(simulation);
        let c = self.material.get_3d_matrix();

        for &(xi, eta, zeta, weight) in gauss_points {
            let d_n = self.get_shape_derivatives_local(xi, eta, zeta);
            let j = self.compute_jacobian_matrix(x, &d_n);
            let det_j = j.determinant();
            
            if det_j <= 0.0 {
                panic!("Negative or zero Jacobian determinant in C3D20 element {}", self.id);
            }
            
            let b = self.compute_b(x, &j, &d_n);
            k += &b.transpose() * c * &b * det_j * weight;
        }
        self.stiffness = k;
    }

    fn get_stiffness(&self) -> DMatrix<f64> {
        self.stiffness.clone()
    }

    fn compute_mass(&self, simulation: &Simulation) -> DMatrix<f64> {
        let mut m = DMatrix::zeros(NUM_NODES, NUM_NODES);
        let gauss_points = Self::get_gauss_points();
        let density = self.material.density;
        let x = self.get_x_local(simulation);

        for &(xi, eta, zeta, weight) in gauss_points {
            let d_n = self.get_shape_derivatives_local(xi, eta, zeta);
            let j = self.compute_jacobian_matrix(x, &d_n);
            let det_j = j.determinant();
            let n = self.get_shape_functions(xi, eta, zeta);
            m += density * &n * &n.transpose() * det_j * weight;
        }
        m
    }

    fn set_lumped_mass(&mut self, mass: &DMatrix<f64>) -> f64 {
        let mut total_mass = 0.0;
        self.lumped_mass.clear();
        for i in 0..NUM_NODES {
            let row_sum: f64 = mass.row(i).sum();
            self.lumped_mass.push(row_sum);
            total_mass += row_sum;
        }
        total_mass
    }

    fn get_lumped_mass(&self) -> &Vec<f64> {
        &self.lumped_mass
    }

    fn get_mass(&self) -> &DMatrix<f64> {
        &self.mass
    }

    fn set_mass(&mut self, mass: DMatrix<f64>) {
        self.mass = mass;
    }

    fn is_active(&self) -> bool {
        self.active
    }

    fn set_active(&mut self, active: bool) {
        self.active = active;
    }

    fn compute_force(&self, displacement: &DVector<f64>) -> DVector<f64> {
        let f_e = &self.stiffness * self.get_u_local_from_displacement(displacement);
        f_e
    }

    fn add_force(&self, simulation: &Simulation, global_force_vector: &mut DVector<f64>) {
        let f_e = &self.stiffness * self.get_u_local(simulation);
        for (local_index, global_index) in self.connectivity.iter().enumerate() {
            for dof in 0..3 {
                global_force_vector[3 * global_index + dof] += f_e[3 * local_index + dof];
            }
        }
    }

    fn compute_element_nodal_properties(&self, simulation: &Simulation) -> ElementFields {
        trace!("Computing element nodal properties for C3D20 element");
        let mut element_fields = ElementFields::new(self.get_connectivity().to_vec());
        
        // Use 20 evaluation points (at node locations) for output
        let node_coords = Self::get_node_coordinates();
        let u_e = self.get_u_local(simulation);
        let x = self.get_x_local(simulation);

        for (nn, &(xi, eta, zeta)) in node_coords.iter().enumerate() {
            let d_n = self.get_shape_derivatives_local(xi, eta, zeta);
            let strain = self.compute_strain(x, &u_e, &d_n);
            let stress = self.compute_stress(x, &u_e, &d_n);
            let stress_vec = Vector6::new(stress[0], stress[1], stress[2], stress[3], stress[4], stress[5]);
            let vm = compute_von_mises(stress_vec);
            
            element_fields.append_to_field("e_xx", nn, strain[0]);
            element_fields.append_to_field("e_yy", nn, strain[1]);
            element_fields.append_to_field("e_zz", nn, strain[2]);
            element_fields.append_to_field("e_xy", nn, strain[3]);
            element_fields.append_to_field("e_yz", nn, strain[4]);
            element_fields.append_to_field("e_xz", nn, strain[5]);
            element_fields.append_to_field("s_xx", nn, stress[0]);
            element_fields.append_to_field("s_yy", nn, stress[1]);
            element_fields.append_to_field("s_zz", nn, stress[2]);
            element_fields.append_to_field("s_xy", nn, stress[3]);
            element_fields.append_to_field("s_yz", nn, stress[4]);
            element_fields.append_to_field("s_xz", nn, stress[5]);
            element_fields.append_to_field("vm", nn, vm);
        }

        element_fields
    }

    fn type_name(&self) -> ElementType {
        ElementType::Brick20
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_shape_functions_partition_of_unity() {
        let elem = Brick20Element::new(0, (0..20).collect(), Material::aluminum());
        
        // Test at several points that shape functions sum to 1
        let test_points = [
            (0.0, 0.0, 0.0),
            (0.5, 0.5, 0.5),
            (-0.5, 0.3, -0.7),
            (0.8, -0.6, 0.2),
        ];
        
        for (xi, eta, zeta) in test_points {
            let n = elem.get_shape_functions(xi, eta, zeta);
            let sum: f64 = n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, 
                "Shape functions don't sum to 1 at ({}, {}, {}): sum = {}", xi, eta, zeta, sum);
        }
    }
    
    #[test]
    fn test_shape_functions_at_nodes() {
        let elem = Brick20Element::new(0, (0..20).collect(), Material::aluminum());
        let node_coords = Brick20Element::get_node_coordinates();
        
        // At each node, only that node's shape function should be 1
        for (i, &(xi, eta, zeta)) in node_coords.iter().enumerate() {
            let n = elem.get_shape_functions(xi, eta, zeta);
            for j in 0..20 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((n[j] - expected).abs() < 1e-12,
                    "N[{}] at node {} = {}, expected {}", j, i, n[j], expected);
            }
        }
    }
    
    #[test]
    fn test_shape_derivative_consistency() {
        let elem = Brick20Element::new(0, (0..20).collect(), Material::aluminum());
        let eps = 1e-6;
        
        // Test at center
        let (xi, eta, zeta) = (0.3, -0.2, 0.4);
        let dn = elem.get_shape_derivatives_local(xi, eta, zeta);
        
        // Numerical derivatives
        let n_xi_p = elem.get_shape_functions(xi + eps, eta, zeta);
        let n_xi_m = elem.get_shape_functions(xi - eps, eta, zeta);
        let n_eta_p = elem.get_shape_functions(xi, eta + eps, zeta);
        let n_eta_m = elem.get_shape_functions(xi, eta - eps, zeta);
        let n_zeta_p = elem.get_shape_functions(xi, eta, zeta + eps);
        let n_zeta_m = elem.get_shape_functions(xi, eta, zeta - eps);
        
        for i in 0..20 {
            let dn_dxi_num = (n_xi_p[i] - n_xi_m[i]) / (2.0 * eps);
            let dn_deta_num = (n_eta_p[i] - n_eta_m[i]) / (2.0 * eps);
            let dn_dzeta_num = (n_zeta_p[i] - n_zeta_m[i]) / (2.0 * eps);
            
            assert!((dn[(i, 0)] - dn_dxi_num).abs() < 1e-5,
                "dN[{}]/dξ: analytical = {}, numerical = {}", i, dn[(i, 0)], dn_dxi_num);
            assert!((dn[(i, 1)] - dn_deta_num).abs() < 1e-5,
                "dN[{}]/dη: analytical = {}, numerical = {}", i, dn[(i, 1)], dn_deta_num);
            assert!((dn[(i, 2)] - dn_dzeta_num).abs() < 1e-5,
                "dN[{}]/dζ: analytical = {}, numerical = {}", i, dn[(i, 2)], dn_dzeta_num);
        }
    }
}
