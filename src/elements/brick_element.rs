// src/elements/brick_element.rs

use crate::simulation::{Simulation, self};
use super::base_element::{BaseElement, Material};
use na::{DMatrix, DVector};
use nalgebra as na;

pub struct BrickElement {
    id: usize,
    connectivity: Vec<usize>,
    material: Material,
    deformation_gradient: DMatrix<f64>, // ... other properties specific to the 8-node brick element
}

impl BrickElement {
    pub fn new(id: usize, connectivity: Vec<usize>, material: Material) -> Self {
        assert_eq!(
            connectivity.len(),
            8,
            "8 nodes required for a brick element"
        );
        BrickElement {
            id,
            connectivity,
            material,
            deformation_gradient: DMatrix::<f64>::identity(3, 3),
            // ... initialize other properties here
        }
    }

    fn get_gauss_points() -> Vec<(f64, f64, f64, f64)> { //xi, eta, zeta, weight
        //compute and return the gauss points here
        //assume 2x2x2 gauss points for now
        let mut gauss_points: Vec<(f64, f64, f64, f64)> = Vec::new();
        let a = 1.0 / (3.0 as f64).sqrt();
        gauss_points.push((-a, -a, -a, 1.0));
        gauss_points.push((a, -a, -a, 1.0));
        gauss_points.push((a, a, -a, 1.0));
        gauss_points.push((-a, a, -a, 1.0));
        gauss_points.push((-a, -a, a, 1.0));
        gauss_points.push((a, -a, a, 1.0));
        gauss_points.push((a, a, a, 1.0));
        gauss_points.push((-a, a, a, 1.0));
        gauss_points
    }

}

impl BaseElement for BrickElement {
    fn get_id(&self) -> usize {
        self.id
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

    fn get_shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> DVector<f64> {
        //compute N_i in local coordinates (xi, eta, zeta)
        //return as a DVector of length 8
        let mut shape_functions: na::Matrix<
            f64,
            na::Dyn,
            na::Const<1>,
            na::VecStorage<f64, na::Dyn, na::Const<1>>,
        > = DVector::<f64>::zeros(8);
        shape_functions[0] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta);
        shape_functions[1] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta);
        shape_functions[2] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta);
        shape_functions[3] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta);
        shape_functions[4] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta);
        shape_functions[5] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta);
        shape_functions[6] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta);
        shape_functions[7] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta);
        shape_functions
    }

    fn get_global_position(
        &self,
        N: &DVector<f64>,
        simulation: &Simulation,
    ) -> na::Vector3<f64> {
        //compute global position of the element given shape functions N
        //return as a na::Vector3<f64>
        let mut global_position = na::Vector3::<f64>::zeros();
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            global_position += N[i] * node.position;
        }
        global_position
    }

    fn get_X(&self, simulation: &Simulation) -> DMatrix<f64> { //get global position of each node
        //compute X matrix
        //return as a DMatrix<f64>
        let mut X = DMatrix::<f64>::zeros(3, 8);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            X[(0, i)] = node.position[0];
            X[(1, i)] = node.position[1];
            X[(2, i)] = node.position[2];
        }
        X
    }

    fn get_shape_derivatives(&self, xi: f64, eta: f64, zeta: f64) -> DMatrix<f64> {
        //we expect to return a 3x8 matrix of shape derivatives. ie. the derivative of N_i wrt. xi, eta, and zeta for each node i
        //derived with sympy
        let matrix_data = vec![
            [-0.125 * (eta - 1.0) * (zeta - 1.0), -0.125 * (xi - 1.0) * (zeta - 1.0), -0.125 * (eta - 1.0) * (xi - 1.0)],
            [0.125 * (eta - 1.0) * (zeta - 1.0), 0.125 * (xi + 1.0) * (zeta - 1.0), 0.125 * (eta - 1.0) * (xi + 1.0)],
            [-0.125 * (eta + 1.0) * (zeta - 1.0), -0.125 * (xi + 1.0) * (zeta - 1.0), -0.125 * (eta + 1.0) * (xi + 1.0)],
            [0.125 * (eta + 1.0) * (zeta - 1.0), 0.125 * (xi - 1.0) * (zeta - 1.0), 0.125 * (eta + 1.0) * (xi - 1.0)],
            [0.125 * (eta - 1.0) * (zeta + 1.0), 0.125 * (xi - 1.0) * (zeta + 1.0), 0.125 * (eta - 1.0) * (xi - 1.0)],
            [-0.125 * (eta - 1.0) * (zeta + 1.0), -0.125 * (xi + 1.0) * (zeta + 1.0), -0.125 * (eta - 1.0) * (xi + 1.0)],
            [0.125 * (eta + 1.0) * (zeta + 1.0), 0.125 * (xi + 1.0) * (zeta + 1.0), 0.125 * (eta + 1.0) * (xi + 1.0)],
            [-0.125 * (eta + 1.0) * (zeta + 1.0), -0.125 * (xi - 1.0) * (zeta + 1.0), -0.125 * (eta + 1.0) * (xi - 1.0)]
        ];
        let dmatrix = na::DMatrix::from_row_slice(8, 3, &matrix_data.concat());
        dmatrix
    }

    fn compute_jacobian_matrix(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DMatrix<f64> {
        //compute the jacobian matrix
        //return as a DMatrix<f64>
        let X = self.get_X(simulation);
        let dN = self.get_shape_derivatives(xi, eta, zeta);
        let jacobian_matrix = dN.transpose() * X.transpose();
        jacobian_matrix
    }

    fn compute_B(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DMatrix<f64> {
        let mut B = DMatrix::<f64>::zeros(6, 24);
        let J = self.compute_jacobian_matrix(xi, eta, zeta, &simulation);
        let J_inv = J.try_inverse().unwrap();
        let N = self.get_shape_functions(xi, eta, zeta);
        let dN = self.get_shape_derivatives(xi, eta, zeta);
        //compute B_i for each node [6x3] 6 strains and 3 displacements per node
        for i in 0..8 {
            //Compute N_I for each node N_I,M = (delta_N_I/delta_local_J)*J_inv(j,M)
            let mut N_I = DVector::<f64>::zeros(3);
            for m in 0..3 {
                let dN_j = dN.row(i); 
                let J_inv_j = J_inv.row(m); // TODO: HIGH BUG POTENTIAL
                 N_I[m] =  J_inv_j.dot(&dN_j);
            }
            //define B_I for each node
            let mut B_I = DMatrix::<f64>::zeros(6, 3);
            B_I[(0, 0)] = N_I[0];
            B_I[(1, 1)] = N_I[1];
            B_I[(2, 2)] = N_I[2];
            B_I[(3, 1)] = N_I[2];
            B_I[(3, 2)] = N_I[1];
            B_I[(4, 0)] = N_I[2];
            B_I[(4, 2)] = N_I[0];
            B_I[(5, 0)] = N_I[1];
            B_I[(5, 1)] = N_I[0];
            
            //populate B matrix
            for j in 0..6 {
                for k in 0..3 {
                    B[(j, 3 * i + k)] = B_I[(j, k)];
                }
            }
        }
        B
    }


    

    fn compute_stiffness(&self, simulation: &Simulation) -> DMatrix<f64> {
        
        let mut K = DMatrix::<f64>::zeros(24, 24);
        let gauss_points = BrickElement::get_gauss_points();
        for (xi, eta, zeta, weight) in gauss_points {
            let B = self.compute_B(xi, eta, zeta, &simulation);
            let J = self.compute_jacobian_matrix(xi, eta, zeta, &simulation);
            let detJ: f64 = J.determinant();
            let B_T: na::Matrix<f64, na::Dyn, na::Dyn, na::VecStorage<f64, na::Dyn, na::Dyn>> = B.transpose();
            let C = self.material.get_3d_matrix();
            let K_point: na::Matrix<f64, na::Dyn, na::Dyn, na::VecStorage<f64, na::Dyn, na::Dyn>> = B_T * C * B * detJ * weight; // HIGH BUG POTENTIAL
            K += K_point;
        }
        //use the gauss quadrature to compute the stiffness matrix
        K
    }

    fn compute_force(&self, simulation: &Simulation) -> DVector<f64> {
        //return 8*3 zero vector
        let mut f = DVector::<f64>::zeros(24);
        return f;
    }

    
    

    // If you wish to override the default methods from BaseElement, you can do so here.
}
