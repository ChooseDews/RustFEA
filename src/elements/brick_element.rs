// src/elements/brick_element.rs
use crate::{simulation::Simulation, utilities::check_for_nans};
use super::base_element::{BaseElement, Material, ElementFields, ElementType};
use nalgebra as na;
use na::{DMatrix, DVector, SMatrix, SVector, Matrix3};
use crate::utilities::compute_von_mises;
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
    #[serde(skip, default = "nodal_positions")]
    nodal_positions: Option<SMatrix<f64, 8, 3>>,
    #[serde(skip, default = "default_internal_force")]
    internal_force: SVector<f64, 24>,
}

fn default_internal_force() -> SVector<f64, 24> {
    SVector::zeros()
}

fn nodal_positions() -> Option<SMatrix<f64, 8, 3>> {
    Option::None
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
            nodal_positions: Option::None,
            internal_force: default_internal_force(),
        }
    }

    fn get_x_local(&self, simulation: &Simulation) -> &SMatrix<f64, 8, 3> { //get global position of each node
        self.nodal_positions.as_ref().expect("Nodal positions not initialized")
    }

    fn compute_x_local(&mut self, simulation: &Simulation) -> SMatrix<f64, 8, 3> {
        let mut x = SMatrix::<f64, 8, 3>::zeros();
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            x[(i, 0)] = node.position[0];
            x[(i, 1)] = node.position[1];
            x[(i, 2)] = node.position[2];
        }
        x
    }


    fn get_gauss_points() -> &'static [(f64, f64, f64, f64)] {
        static A: f64 = 0.5773502691896257; // 1/sqrt(3)
        static GAUSS_POINTS: [(f64, f64, f64, f64); 8] = [
            (-A, -A, -A, 1.0),
            (A, -A, -A, 1.0),
            (A, A, -A, 1.0),
            (-A, A, -A, 1.0),
            (-A, -A, A, 1.0),
            (A, -A, A, 1.0),
            (A, A, A, 1.0),
            (-A, A, A, 1.0),
        ];
        &GAUSS_POINTS
    }

    fn get_corner_points() -> &'static [(f64, f64, f64, f64)] {
        static CORNER_POINTS: [(f64, f64, f64, f64); 8] = [
            (-1.0, -1.0, -1.0, 1.0),
            (1.0, -1.0, -1.0, 1.0),
            (1.0, 1.0, -1.0, 1.0),
            (-1.0, 1.0, -1.0, 1.0),
            (-1.0, -1.0, 1.0, 1.0),
            (1.0, -1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0, 1.0),
            (-1.0, 1.0, 1.0, 1.0),
        ];
        &CORNER_POINTS
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

    fn compute_b(&self, x: &SMatrix<f64, 8, 3>, j: &Matrix3<f64>, d_n: &SMatrix<f64, 8, 3>) -> SMatrix<f64, 6, 24> {
        let mut b = SMatrix::<f64, 6, 24>::zeros();
        let j_inv = j.try_inverse().unwrap();
        let mut n_i = [0.0; 3];
        //compute B_i for each node [6x3] 6 strains and 3 displacements per node
        for i in 0..8 {
            //Compute N_I for each node N_I,M = (delta_N_I/delta_local_J)*J_inv(j,M)
            for m in 0..3 {
                 n_i[m] =  j_inv.row(m).dot(&d_n.row(i));
            }
            //new method to handle B_I
            let b_i = SMatrix::<f64, 6, 3>::new(
                n_i[0], 0.0, 0.0,
                0.0, n_i[1], 0.0,
                0.0, 0.0, n_i[2],

                0.0, n_i[2], n_i[1],
                n_i[2], 0.0, n_i[0],
                n_i[1], n_i[0], 0.0
            );

            for j in 0..6 {
                for k in 0..3 {
                    b[(j, 3 * i + k)] = b_i[(j, k)];
                }
            }
        }
        b
    }


    fn compute_jacobian_matrix(&self, x: &SMatrix<f64, 8, 3>, d_n: &SMatrix<f64, 8, 3>) -> Matrix3<f64>  {
        d_n.transpose() * x
    }



    fn get_shape_derivatives_local(&self, xi: f64, eta: f64, zeta: f64) -> SMatrix<f64, 8, 3> {
        // Pre-compute common terms
        let xi_m = xi - 1.0;
        let xi_p = xi + 1.0;
        let eta_m = eta - 1.0;
        let eta_p = eta + 1.0;
        let zeta_m = zeta - 1.0;
        let zeta_p = zeta + 1.0;

        SMatrix::from_row_slice(&[
            -0.125 * eta_m * zeta_m, -0.125 * xi_m * zeta_m, -0.125 * eta_m * xi_m,
            0.125 * eta_m * zeta_m, 0.125 * xi_p * zeta_m, 0.125 * eta_m * xi_p,
            -0.125 * eta_p * zeta_m, -0.125 * xi_p * zeta_m, -0.125 * eta_p * xi_p,
            0.125 * eta_p * zeta_m, 0.125 * xi_m * zeta_m, 0.125 * eta_p * xi_m,
            0.125 * eta_m * zeta_p, 0.125 * xi_m * zeta_p, 0.125 * eta_m * xi_m,
            -0.125 * eta_m * zeta_p, -0.125 * xi_p * zeta_p, -0.125 * eta_m * xi_p,
            0.125 * eta_p * zeta_p, 0.125 * xi_p * zeta_p, 0.125 * eta_p * xi_p,
            -0.125 * eta_p * zeta_p, -0.125 * xi_m * zeta_p, -0.125 * eta_p * xi_m
        ])
    }   


    fn compute_stress(&self, x: &SMatrix<f64, 8, 3>, u: &SVector<f64, 24>, d_n: &SMatrix<f64, 8, 3>) -> SVector<f64, 6> {
        let j = self.compute_jacobian_matrix(x, d_n);
        let b = self.compute_b(x, &j, d_n);
        let c = self.material.get_3d_matrix();
        c * b * u
    }

    fn compute_strain(&self, x: &SMatrix<f64, 8, 3>, u: &SVector<f64, 24>, d_n: &SMatrix<f64, 8, 3>) -> SVector<f64, 6> {
        let J = self.compute_jacobian_matrix(x, d_n);
        self.compute_b(x, &J, d_n) * u
    }

    fn get_shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> SVector<f64, 8> {
        // Pre-compute common terms
        let xi_m = 1.0 - xi;
        let xi_p = 1.0 + xi;
        let eta_m = 1.0 - eta;
        let eta_p = 1.0 + eta;
        let zeta_m = 1.0 - zeta;
        let zeta_p = 1.0 + zeta;

        SVector::from([
            0.125 * xi_m * eta_m * zeta_m,
            0.125 * xi_p * eta_m * zeta_m,
            0.125 * xi_p * eta_p * zeta_m,
            0.125 * xi_m * eta_p * zeta_m,
            0.125 * xi_m * eta_m * zeta_p,
            0.125 * xi_p * eta_m * zeta_p,
            0.125 * xi_p * eta_p * zeta_p,
            0.125 * xi_m * eta_p * zeta_p
        ])
    }



}

#[typetag::serde]
impl BaseElement for BrickElement {
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
        // Pre-compute common terms
        let xi_m = xi - 1.0;
        let xi_p = xi + 1.0;
        let eta_m = eta - 1.0;
        let eta_p = eta + 1.0;
        let zeta_m = zeta - 1.0;
        let zeta_p = zeta + 1.0;

        let eta_m_zeta_m = eta_m * zeta_m;
        let eta_p_zeta_m = eta_p * zeta_m;
        let eta_m_zeta_p = eta_m * zeta_p;
        let eta_p_zeta_p = eta_p * zeta_p;

        let xi_m_zeta_m = xi_m * zeta_m;
        let xi_p_zeta_m = xi_p * zeta_m;
        let xi_m_zeta_p = xi_m * zeta_p;
        let xi_p_zeta_p = xi_p * zeta_p;

        let xi_m_eta_m = xi_m * eta_m;
        let xi_p_eta_m = xi_p * eta_m;
        let xi_m_eta_p = xi_m * eta_p;
        let xi_p_eta_p = xi_p * eta_p;

        // Create the matrix directly with computed values
        DMatrix::from_row_slice(8, 3, &[
            -0.125 * eta_m_zeta_m, -0.125 * xi_m_zeta_m, -0.125 * xi_m_eta_m,
             0.125 * eta_m_zeta_m,  0.125 * xi_p_zeta_m,  0.125 * xi_p_eta_m,
            -0.125 * eta_p_zeta_m, -0.125 * xi_p_zeta_m, -0.125 * xi_p_eta_p,
             0.125 * eta_p_zeta_m,  0.125 * xi_m_zeta_m,  0.125 * xi_m_eta_p,
             0.125 * eta_m_zeta_p,  0.125 * xi_m_zeta_p,  0.125 * xi_m_eta_m,
            -0.125 * eta_m_zeta_p, -0.125 * xi_p_zeta_p, -0.125 * xi_p_eta_m,
             0.125 * eta_p_zeta_p,  0.125 * xi_p_zeta_p,  0.125 * xi_p_eta_p,
            -0.125 * eta_p_zeta_p, -0.125 * xi_m_zeta_p, -0.125 * xi_m_eta_p,
        ])
    }



    fn get_b(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DMatrix<f64> {
        unimplemented!()
    }

    
    fn compute_stiffness(&mut self, simulation: &Simulation){
        trace!("Computing stiffness matrix for brick element");
        let mut k = empty_element_matrix();
        let gauss_points = BrickElement::get_corner_points();
        let x = self.get_x_local(simulation);
        let c = self.material.get_3d_matrix();

        for &(xi, eta, zeta, weight) in gauss_points {
            let d_n = self.get_shape_derivatives_local(xi, eta, zeta);
            let j = self.compute_jacobian_matrix(x, &d_n);
            let b = self.compute_b(x, &j, &d_n);
            k += b.transpose() * c * b * j.determinant() * weight;
        }
        self.stiffness = k;
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
        let mut M = DMatrix::<f64>::zeros(8, 8);
        let gauss_points = BrickElement::get_gauss_points();
        let density = self.material.density;
        let x = self.get_x_local(simulation);

        for &(xi, eta, zeta, weight) in gauss_points {
            let d_n = self.get_shape_derivatives_local(xi, eta, zeta);
            let J = self.compute_jacobian_matrix(x, &d_n);
            let N = self.get_shape_functions(xi, eta, zeta);
            M += density * N * N.transpose() * J.determinant() * weight;
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

    fn add_force(&self, simulation: &Simulation, global_force_vector: &mut DVector<f64>) {
        let f_e = self.stiffness * self.get_u_local(simulation);
        for (local_index, global_index) in self.connectivity.iter().enumerate() {
            for dof in 0..3 {
                global_force_vector[3 * global_index + dof] += f_e[3 * local_index + dof];
            }
        }
    }

    fn compute_element_nodal_properties(&self, simulation: &Simulation) -> ElementFields {
        trace!("Computing element nodal properties for brick element");
        let mut element_fields = ElementFields::new(self.get_connectivity().to_vec());
        let gauss_points = BrickElement::get_gauss_points();
        let u_e = self.get_u_local(simulation);
        let x = self.get_x_local(simulation);

        for (nn, &(xi, eta, zeta, _)) in gauss_points.iter().enumerate() {
            let d_n = self.get_shape_derivatives_local(xi, eta, zeta);
            let strain = self.compute_strain(x, &u_e, &d_n);
            let stress = self.compute_stress(x, &u_e, &d_n);
            let vm = compute_von_mises(stress);
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
        ElementType::Brick
    }
}
