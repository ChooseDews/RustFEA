// src/elements/brick_element.rs
use crate::{simulation::Simulation, utilities::check_for_nans};
use super::base_element::{BaseElement, Material, ElementFields, ElementType};
use nalgebra as na;
use na::{DMatrix, DVector, SMatrix, SVector, Matrix3};
use crate::utilities::{compute_von_mises};
use serde::{Serialize, Deserialize};
use log::{debug, trace};


#[derive(Serialize, Deserialize, Debug)] 
pub struct BrickElement {
    id: usize,
    connectivity: Vec<usize>,
    material: Material,
    #[serde(skip, default = "default_deformation_gradient")]
    deformation_gradient: DMatrix<f64>, // ... other properties specific to the 8-node brick element
    #[serde(skip, default = "empty_element_matrix")]
    stiffness: SMatrix<f64, 24, 24>,
    #[serde(skip, default = "default_zero_matrix")]
    mass: DMatrix<f64>,
    #[serde(skip, default = "Vec::new")]//lumped mass matrix
    lumped_mass: Vec<f64>,
    active: bool,
}

fn default_deformation_gradient() -> DMatrix<f64> {
    DMatrix::<f64>::identity(3, 3)
}
// Add this function outside the struct
fn default_zero_matrix() -> DMatrix<f64> {
    DMatrix::zeros(24, 24)
}

fn empty_element_matrix() -> na::SMatrix<f64, 24, 24> {
    na::SMatrix::zeros()
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
            stiffness: empty_element_matrix(),
            mass: default_zero_matrix(),
            lumped_mass: Vec::new(),
            active: true,
        }
    }

    fn get_x_local(&self, simulation: &Simulation) -> SMatrix<f64, 8, 3> { //get global position of each node
        let mut X = SMatrix::<f64, 8, 3>::zeros();
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            X[(i, 0)] = node.position[0];
            X[(i, 1)] = node.position[1];
            X[(i, 2)] = node.position[2];
        }
        X
    }

    fn get_gauss_points() -> Vec<(f64, f64, f64, f64)> { //xi, eta, zeta, weight
        trace!("Generating Gauss points for brick element");
        //compute and return the gauss points here
        //assume 2x2x2 gauss points for now
        let mut gauss_points: Vec<(f64, f64, f64, f64)> = Vec::new();
        let a = 1.0 / (3.0 as f64).sqrt();
        //same order as node numbering (useful for computing properties at nodes)
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

    fn get_corner_points() -> Vec<(f64, f64, f64, f64)> {
        //in local coordinates ie. -1 to 1
        let mut corner_points: Vec<(f64, f64, f64, f64)> = Vec::new();
        corner_points.push((-1.0, -1.0, -1.0, 1.0));
        corner_points.push((1.0, -1.0, -1.0, 1.0));
        corner_points.push((1.0, 1.0, -1.0, 1.0));
        corner_points.push((-1.0, 1.0, -1.0, 1.0));
        corner_points.push((-1.0, -1.0, 1.0, 1.0));
        corner_points.push((1.0, -1.0, 1.0, 1.0));
        corner_points.push((1.0, 1.0, 1.0, 1.0));
        corner_points.push((-1.0, 1.0, 1.0, 1.0));
        corner_points
    }

    fn get_u_local(&self, simulation: &Simulation) -> SVector<f64, 24> {
        let mut u = SVector::<f64, 24>::zeros();
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            u[3 * i] = node.displacement[0];
            u[3 * i + 1] = node.displacement[1];
            u[3 * i + 2] = node.displacement[2];
        }
        u
    }

    fn compute_b(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> SMatrix<f64, 6, 24> {
        let mut B = SMatrix::<f64, 6, 24>::zeros();
        let J = self.compute_jacobian_matrix(xi, eta, zeta, &simulation);
        let J_inv = J.try_inverse().unwrap();
        let dN = self.get_shape_derivatives_local(xi, eta, zeta);
        //compute B_i for each node [6x3] 6 strains and 3 displacements per node
        for i in 0..8 {
            //Compute N_I for each node N_I,M = (delta_N_I/delta_local_J)*J_inv(j,M)
            let mut N_I = SVector::<f64, 3>::zeros();
            for m in 0..3 {
                 N_I[m] =  J_inv.row(m).dot(&dN.row(i));
            }
            //define B_I for each node
            let mut B_I = SMatrix::<f64, 6, 3>::zeros();
            B_I[(0, 0)] = N_I[0];
            B_I[(1, 1)] = N_I[1];
            B_I[(2, 2)] = N_I[2];
            B_I[(3, 1)] = N_I[2];
            B_I[(3, 2)] = N_I[1];
            B_I[(4, 0)] = N_I[2];
            B_I[(4, 2)] = N_I[0];
            B_I[(5, 0)] = N_I[1];
            B_I[(5, 1)] = N_I[0];
            
            for j in 0..6 {
                for k in 0..3 {
                    B[(j, 3 * i + k)] = B_I[(j, k)];
                }
            }
        }
        B
    }


    fn compute_jacobian_matrix(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> Matrix3<f64>  {
        let x = self.get_x_local(simulation); //[x] 8x3
        let dn = self.get_shape_derivatives_local(xi, eta, zeta); //[dN] 8x3
        dn.transpose() * x
    }



    fn get_shape_derivatives_local(&self, xi: f64, eta: f64, zeta: f64) -> SMatrix<f64, 8, 3> {
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
        na::SMatrix::from_row_slice(&matrix_data.concat())
    }   


    fn compute_stress(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> SVector<f64, 6> {
        //compute the stress at a given point
        let B = self.compute_b(xi, eta, zeta, &simulation);
        let C = self.material.get_3d_matrix();
        let u = self.get_u(simulation);
        let stress = C * B * u;
        stress
    }

    fn compute_strain(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> SVector<f64, 6> {
        //compute the strain at a given point
        let B = self.compute_b(xi, eta, zeta, &simulation);
        let u = self.get_u_local(simulation);
        let strain = B * u;
        strain
    }



}

#[typetag::serde]
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

    fn get_x(&self, simulation: &Simulation) -> DMatrix<f64> { //get global position of each node
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

    fn get_u(&self, simulation: &Simulation) -> DVector<f64> {
        let mut u = DVector::<f64>::zeros(24);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            u[3 * i] = node.displacement[0];
            u[3 * i + 1] = node.displacement[1];
            u[3 * i + 2] = node.displacement[2];
        }
        u
    }

    fn get_shape_derivatives(&self, xi: f64, eta: f64, zeta: f64) -> DMatrix<f64> {
        //we expect to return a 8x3 matrix of shape derivatives. ie. the derivative of N_i wrt. xi, eta, and zeta for each node i
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



    fn get_b(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DMatrix<f64> {
        todo!()
    }

    
    fn compute_stiffness(&mut self, simulation: &Simulation){
        trace!("Computing stiffness matrix for brick element");
        let mut K = empty_element_matrix();
        let gauss_points = BrickElement::get_corner_points();
        for (xi, eta, zeta, weight) in gauss_points {
            let B = self.compute_b(xi, eta, zeta, &simulation);
            let J = self.compute_jacobian_matrix(xi, eta, zeta, &simulation);
            let C = self.material.get_3d_matrix();
            K +=  B.transpose() * C * B * J.determinant() * weight;
        }
        //use the gauss quadrature to compute the stiffness matrix
        self.stiffness = K;
    }


    fn get_stiffness(&self) -> DMatrix<f64> {
        //construct dynamic matrix from fixed 24x24 matrix
        let mut k =  DMatrix::zeros(24, 24);
        for i in 0..24 {
            for j in 0..24 {
                k[(i, j)] = self.stiffness[(i, j)];
            }
        }
        k
    }


    fn compute_mass(&self, simulation: &Simulation) -> DMatrix<f64> {
        //compute the mass matrix [M] = intergrate(density * N^T * N * detJ * weight)
        let mut M = DMatrix::<f64>::zeros(8, 8);
        let gauss_points = BrickElement::get_gauss_points();
        let density = self.material.density;
        for (xi, eta, zeta, weight) in gauss_points {
            let J = self.compute_jacobian_matrix(xi, eta, zeta, &simulation);
            let detJ: f64 = J.determinant();
            let N = self.get_shape_functions(xi, eta, zeta);
            let N_T = N.transpose();
            let M_point = density * N * N_T * detJ * weight;
            M += M_point;
        }
        M
    }

    fn set_lumped_mass(&mut self, mass: &DMatrix<f64>) -> f64 {
        //sum rows and save to lumped_mass
        let mut total_mass = 0.0;
        for i in 0..8 {
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


    fn compute_force(&self, simulation: &Simulation) -> DVector<f64> {
        let u_e = self.get_u_local(simulation);
        let f_e = self.stiffness * u_e;
        DVector::<f64>::from_column_slice(f_e.as_slice())
    }

    fn compute_element_nodal_properties(&self, simulation: &Simulation) -> ElementFields {
        trace!("Computing element nodal properties for brick element");
        //compute the nodal properties for each node
        let mut element_feilds = ElementFields::new(self.get_connectivity().to_vec());
        let gauss_points = BrickElement::get_gauss_points();

        let mut nn = 0; //node number
        for (xi, eta, zeta, _) in gauss_points {
            let strain = self.compute_strain(xi, eta, zeta, &simulation);
            let stress = self.compute_stress(xi, eta, zeta, &simulation);
            element_feilds.append_to_feild("e_xx", nn, strain[0]);
            element_feilds.append_to_feild("e_yy", nn, strain[1]);
            element_feilds.append_to_feild("e_zz", nn, strain[2]);
            element_feilds.append_to_feild("e_xy", nn, strain[3]);
            element_feilds.append_to_feild("e_yz", nn, strain[4]);
            element_feilds.append_to_feild("e_xz", nn, strain[5]);
            element_feilds.append_to_feild("s_xx", nn, stress[0]);
            element_feilds.append_to_feild("s_yy", nn, stress[1]);
            element_feilds.append_to_feild("s_zz", nn, stress[2]);
            element_feilds.append_to_feild("s_xy", nn, stress[3]);
            element_feilds.append_to_feild("s_yz", nn, stress[4]);
            element_feilds.append_to_feild("s_xz", nn, stress[5]);
            let vm = compute_von_mises(stress);
            element_feilds.append_to_feild("vm", nn, vm);
            nn += 1;
        }

        return element_feilds;
    }

    fn type_name(&self) -> ElementType {
        ElementType::Brick
    }
}
