// src/elements/tetrahedral_element.rs
//! 4-node linear tetrahedral element (C3D4)
//!
//! This is the simplest 3D solid element. It uses natural/volume coordinates
//! (L1, L2, L3, L4) where L1 + L2 + L3 + L4 = 1.
//!
//! Node ordering (standard convention):
//!           3
//!          /|\
//!         / | \
//!        /  |  \
//!       /   |   \
//!      /    |    \
//!     /     |     \
//!    0------|------2
//!     \     |     /
//!      \    |    /
//!       \   |   /
//!        \  |  /
//!         \ | /
//!          \|/
//!           1
//!
//! The element uses linear shape functions:
//!   N1 = L1 = 1 - L2 - L3 - L4  (node 0)
//!   N2 = L2                      (node 1)
//!   N3 = L3                      (node 2)
//!   N4 = L4                      (node 3)
//!
//! For a linear tetrahedron, the shape function derivatives are CONSTANT
//! throughout the element, so only 1 Gauss point is needed (at the centroid).

use super::base_element::{BaseElement, ElementFields, ElementType, Material};
use crate::utilities::compute_von_mises;
use crate::{simulation::Simulation, utilities::check_for_nans};
use log::{debug, trace};
use na::{DMatrix, DVector, Matrix3, Matrix4, SMatrix, SVector, Vector3, Vector4};
use nalgebra as na;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct TetElement {
    id: usize,
    connectivity: Vec<usize>,
    material: Material,
    #[serde(skip, default = "default_deformation_gradient")]
    deformation_gradient: DMatrix<f64>,
    #[serde(skip, default = "empty_element_matrix")]
    stiffness: SMatrix<f64, 12, 12>,
    #[serde(skip, default = "default_zero_matrix")]
    stiffness_dmatrix: DMatrix<f64>,
    #[serde(skip, default = "default_zero_matrix")]
    mass: DMatrix<f64>,
    #[serde(skip, default = "Vec::new")]
    lumped_mass: Vec<f64>,
    active: bool,
    #[serde(skip, default = "nodal_positions")]
    nodal_positions: Option<SMatrix<f64, 4, 3>>,
    #[serde(skip, default = "default_internal_force")]
    internal_force: SVector<f64, 12>,
    #[serde(skip, default)]
    volume: f64,
    #[serde(skip, default = "default_b_matrix")]
    b_matrix: Option<SMatrix<f64, 6, 12>>,
}

fn default_internal_force() -> SVector<f64, 12> {
    SVector::zeros()
}

fn nodal_positions() -> Option<SMatrix<f64, 4, 3>> {
    None
}

fn default_deformation_gradient() -> DMatrix<f64> {
    DMatrix::<f64>::identity(3, 3)
}

fn default_zero_matrix() -> DMatrix<f64> {
    DMatrix::zeros(12, 12)
}

fn empty_element_matrix() -> SMatrix<f64, 12, 12> {
    SMatrix::zeros()
}

fn default_b_matrix() -> Option<SMatrix<f64, 6, 12>> {
    None
}

impl TetElement {
    pub fn new(id: usize, connectivity: Vec<usize>, material: Material) -> Self {
        assert_eq!(
            connectivity.len(),
            4,
            "4 nodes required for a tetrahedral element"
        );
        TetElement {
            id,
            connectivity,
            material,
            deformation_gradient: DMatrix::<f64>::identity(3, 3),
            stiffness: empty_element_matrix(),
            stiffness_dmatrix: default_zero_matrix(),
            mass: default_zero_matrix(),
            lumped_mass: Vec::new(),
            active: true,
            nodal_positions: None,
            internal_force: default_internal_force(),
            volume: 0.0,
            b_matrix: None,
        }
    }

    fn get_x_local(&self, simulation: &Simulation) -> &SMatrix<f64, 4, 3> {
        self.nodal_positions
            .as_ref()
            .expect("Nodal positions not initialized")
    }

    fn compute_x_local(&mut self, simulation: &Simulation) -> SMatrix<f64, 4, 3> {
        let mut x = SMatrix::<f64, 4, 3>::zeros();
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            x[(i, 0)] = node.position[0];
            x[(i, 1)] = node.position[1];
            x[(i, 2)] = node.position[2];
        }
        x
    }

    /// V = (1/6) * |det([x2-x1, x3-x1, x4-x1])|
    fn compute_volume(x: &SMatrix<f64, 4, 3>) -> f64 {
        let x1 = Vector3::new(x[(0, 0)], x[(0, 1)], x[(0, 2)]);
        let x2 = Vector3::new(x[(1, 0)], x[(1, 1)], x[(1, 2)]);
        let x3 = Vector3::new(x[(2, 0)], x[(2, 1)], x[(2, 2)]);
        let x4 = Vector3::new(x[(3, 0)], x[(3, 1)], x[(3, 2)]);

        let v1 = x2 - x1;
        let v2 = x3 - x1;
        let v3 = x4 - x1;

        let cross = v1.cross(&v2);
        (cross.dot(&v3) / 6.0).abs()
    }

    /// J = [x2-x1, x3-x1, x4-x1]^T
    fn compute_jacobian(x: &SMatrix<f64, 4, 3>) -> Matrix3<f64> {
        let x1 = Vector3::new(x[(0, 0)], x[(0, 1)], x[(0, 2)]);
        let x2 = Vector3::new(x[(1, 0)], x[(1, 1)], x[(1, 2)]);
        let x3 = Vector3::new(x[(2, 0)], x[(2, 1)], x[(2, 2)]);
        let x4 = Vector3::new(x[(3, 0)], x[(3, 1)], x[(3, 2)]);

        let v1 = x2 - x1;
        let v2 = x3 - x1;
        let v3 = x4 - x1;

        Matrix3::from_columns(&[v1, v2, v3])
    }

    /// Compute strain-displacement (B) matrix. For linear tet, B is constant.
    ///
    /// Shape function derivatives in natural coordinates:
    /// dN/dL = [-1, -1, -1; 1, 0, 0; 0, 1, 0; 0, 0, 1]
    /// Transform to global: dN/dx = dN/dL * J^(-1)
    fn compute_b_matrix(x: &SMatrix<f64, 4, 3>) -> SMatrix<f64, 6, 12> {
        let j = Self::compute_jacobian(x);
        let j_inv = j
            .try_inverse()
            .expect("Jacobian is singular - degenerate element");

        let dn_dl = SMatrix::<f64, 4, 3>::from_row_slice(&[
            -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
        ]);

        let dn_dx = dn_dl * j_inv;

        // B_i = [dNi/dx, 0, 0; 0, dNi/dy, 0; 0, 0, dNi/dz;
        //        dNi/dy, dNi/dx, 0; 0, dNi/dz, dNi/dy; dNi/dz, 0, dNi/dx]
        let mut b = SMatrix::<f64, 6, 12>::zeros();

        for i in 0..4 {
            let dnx = dn_dx[(i, 0)];
            let dny = dn_dx[(i, 1)];
            let dnz = dn_dx[(i, 2)];
            let col = i * 3;

            b[(0, col)] = dnx;
            b[(1, col + 1)] = dny;
            b[(2, col + 2)] = dnz;
            b[(3, col)] = dny;
            b[(3, col + 1)] = dnx;
            b[(4, col + 1)] = dnz;
            b[(4, col + 2)] = dny;
            b[(5, col)] = dnz;
            b[(5, col + 2)] = dnx;
        }

        b
    }

    fn get_u_local(&self, simulation: &Simulation) -> SVector<f64, 12> {
        let mut u = SVector::<f64, 12>::zeros();
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            u[3 * i] = node.displacement[0];
            u[3 * i + 1] = node.displacement[1];
            u[3 * i + 2] = node.displacement[2];
        }
        u
    }

    /// Linear shape functions at natural coordinates (L2, L3, L4). L1 = 1 - L2 - L3 - L4.
    fn get_shape_functions(l2: f64, l3: f64, l4: f64) -> SVector<f64, 4> {
        let l1 = 1.0 - l2 - l3 - l4;
        SVector::from([l1, l2, l3, l4])
    }

    /// Single Gauss point at centroid (L1=L2=L3=L4=1/4), weight = 1/6
    fn get_gauss_points() -> &'static [(f64, f64, f64, f64)] {
        static GAUSS_POINTS: [(f64, f64, f64, f64); 1] = [(0.25, 0.25, 0.25, 1.0 / 6.0)];
        &GAUSS_POINTS
    }

    fn compute_stress(&self, u: &SVector<f64, 12>) -> SVector<f64, 6> {
        let b = self.b_matrix.as_ref().expect("B matrix not initialized");
        let c = self.material.get_3d_matrix();
        c * (b * u)
    }

    fn compute_strain(&self, u: &SVector<f64, 12>) -> SVector<f64, 6> {
        let b = self.b_matrix.as_ref().expect("B matrix not initialized");
        b * u
    }
}

#[typetag::serde]
impl BaseElement for TetElement {
    fn get_id(&self) -> usize {
        self.id
    }

    fn initialize(&mut self, simulation: &Simulation) {
        let x = self.compute_x_local(simulation);
        self.nodal_positions = Some(x);
        self.volume = Self::compute_volume(&x);
        self.b_matrix = Some(Self::compute_b_matrix(&x));
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
        let mut global_position = na::Vector3::<f64>::zeros();
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            global_position += n[i] * node.position;
        }
        global_position
    }

    fn get_x(&self, simulation: &Simulation) -> DMatrix<f64> {
        let mut x = DMatrix::<f64>::zeros(3, 4);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            x[(0, i)] = node.position[0];
            x[(1, i)] = node.position[1];
            x[(2, i)] = node.position[2];
        }
        x
    }

    fn get_u(&self, simulation: &Simulation) -> DVector<f64> {
        let mut u = DVector::<f64>::zeros(12);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            u[3 * i] = node.displacement[0];
            u[3 * i + 1] = node.displacement[1];
            u[3 * i + 2] = node.displacement[2];
        }
        u
    }

    fn get_shape_derivatives(&self, _l2: f64, _l3: f64, _l4: f64) -> DMatrix<f64> {
        DMatrix::from_row_slice(
            4,
            3,
            &[
                -1.0, -1.0, -1.0, // dN1/dL
                1.0, 0.0, 0.0, // dN2/dL
                0.0, 1.0, 0.0, // dN3/dL
                0.0, 0.0, 1.0, // dN4/dL
            ],
        )
    }

    fn get_b(&self, _l2: f64, _l3: f64, _l4: f64, _simulation: &Simulation) -> DMatrix<f64> {
        let b = self.b_matrix.as_ref().expect("B matrix not initialized");
        DMatrix::from_row_slice(6, 12, b.as_slice())
    }

    fn compute_stiffness(&mut self, simulation: &Simulation) {
        trace!("Computing stiffness matrix for tetrahedral element");

        // K = V * B^T * C * B
        let b = self.b_matrix.as_ref().expect("B matrix not initialized");
        let c = self.material.get_3d_matrix();
        let v = self.volume;

        self.stiffness = v * b.transpose() * c * b;
        self.stiffness_dmatrix = DMatrix::from_fn(12, 12, |i, j| self.stiffness[(i, j)]);
    }

    fn get_stiffness(&self) -> &DMatrix<f64> {
        &self.stiffness_dmatrix
    }

    fn compute_mass(&self, simulation: &Simulation) -> DMatrix<f64> {
        // M_ii = ρV/10, M_ij = ρV/20 (i≠j)
        let density = self.material.density;
        let v = self.volume;

        let mut m = DMatrix::<f64>::zeros(4, 4);
        let diag = density * v / 10.0;
        let off_diag = density * v / 20.0;

        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    m[(i, j)] = diag;
                } else {
                    m[(i, j)] = off_diag;
                }
            }
        }
        m
    }

    fn set_lumped_mass(&mut self, mass: &DMatrix<f64>) -> f64 {
        let mut total_mass = 0.0;
        for i in 0..4 {
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

    fn add_force(&self, simulation: &Simulation, global_force_vector: &mut DVector<f64>) {
        let f_e = self.stiffness * self.get_u_local(simulation);
        for (local_index, global_index) in self.connectivity.iter().enumerate() {
            for dof in 0..3 {
                global_force_vector[3 * global_index + dof] += f_e[3 * local_index + dof];
            }
        }
    }

    fn compute_element_nodal_properties(&self, simulation: &Simulation) -> ElementFields {
        trace!("Computing element nodal properties for tetrahedral element");
        let mut element_fields = ElementFields::new(self.get_connectivity().to_vec());
        let u_e = self.get_u_local(simulation);

        // Stress/strain constant throughout linear tet - assign same to all nodes
        let strain = self.compute_strain(&u_e);
        let stress = self.compute_stress(&u_e);
        let vm = compute_von_mises(stress);

        for nn in 0..4 {
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
        ElementType::Tetrahedral
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference tet: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    fn create_reference_tet() -> SMatrix<f64, 4, 3> {
        SMatrix::from_row_slice(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
    }

    #[test]
    fn test_volume_computation() {
        let x = create_reference_tet();
        let volume = TetElement::compute_volume(&x);
        assert!(
            (volume - 1.0 / 6.0).abs() < 1e-10,
            "Volume should be 1/6, got {}",
            volume
        );
    }

    #[test]
    fn test_volume_scaled() {
        let x =
            SMatrix::from_row_slice(&[0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0]);
        let volume = TetElement::compute_volume(&x);
        assert!(
            (volume - 8.0 / 6.0).abs() < 1e-10,
            "Volume should be 8/6, got {}",
            volume
        );
    }

    #[test]
    fn test_shape_functions_partition_of_unity() {
        let test_points = [
            (0.25, 0.25, 0.25),
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.5, 0.25, 0.25),
        ];

        for (l2, l3, l4) in test_points {
            let n = TetElement::get_shape_functions(l2, l3, l4);
            let sum: f64 = n.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "Shape functions should sum to 1 at ({}, {}, {}), got {}",
                l2,
                l3,
                l4,
                sum
            );
        }
    }

    #[test]
    fn test_shape_functions_kronecker_delta() {
        // Node natural coords (L2, L3, L4):
        // Node 0: (0, 0, 0), Node 1: (1, 0, 0), Node 2: (0, 1, 0), Node 3: (0, 0, 1)
        let node_coords = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ];

        for (i, &(l2, l3, l4)) in node_coords.iter().enumerate() {
            let n = TetElement::get_shape_functions(l2, l3, l4);
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (n[j] - expected).abs() < 1e-10,
                    "N{}({}) should be {}, got {}",
                    j,
                    i,
                    expected,
                    n[j]
                );
            }
        }
    }

    #[test]
    fn test_jacobian_matrix() {
        let x = create_reference_tet();
        let j = TetElement::compute_jacobian(&x);

        let expected = Matrix3::<f64>::identity();
        for i in 0..3 {
            for k in 0..3 {
                let diff: f64 = j[(i, k)] - expected[(i, k)];
                assert!(
                    diff.abs() < 1e-10,
                    "Jacobian[{},{}] should be {}, got {}",
                    i,
                    k,
                    expected[(i, k)],
                    j[(i, k)]
                );
            }
        }
    }

    #[test]
    fn test_b_matrix_strain_computation() {
        let x = create_reference_tet();
        let b = TetElement::compute_b_matrix(&x);

        // Uniform tension in x: u = [0.01*x, 0, 0]
        let u = SVector::<f64, 12>::from_row_slice(&[
            0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]);

        let strain = b * u;

        assert!(
            (strain[0] - 0.01).abs() < 1e-10,
            "ε_xx should be 0.01, got {}",
            strain[0]
        );
        assert!(
            strain[1].abs() < 1e-10,
            "ε_yy should be 0, got {}",
            strain[1]
        );
        assert!(
            strain[2].abs() < 1e-10,
            "ε_zz should be 0, got {}",
            strain[2]
        );
        assert!(
            strain[3].abs() < 1e-10,
            "γ_xy should be 0, got {}",
            strain[3]
        );
        assert!(
            strain[4].abs() < 1e-10,
            "γ_yz should be 0, got {}",
            strain[4]
        );
        assert!(
            strain[5].abs() < 1e-10,
            "γ_xz should be 0, got {}",
            strain[5]
        );
    }

    #[test]
    fn test_rigid_body_translation() {
        let x = create_reference_tet();
        let b = TetElement::compute_b_matrix(&x);

        let u = SVector::<f64, 12>::from_row_slice(&[
            1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0,
        ]);

        let strain = b * u;

        for i in 0..6 {
            assert!(
                strain[i].abs() < 1e-10,
                "Strain component {} should be 0 for rigid translation, got {}",
                i,
                strain[i]
            );
        }
    }

    #[test]
    fn test_stiffness_matrix_symmetry() {
        let material = Material::new(200e9, 0.3, 7800.0);

        let x = create_reference_tet();
        let b = TetElement::compute_b_matrix(&x);
        let c = material.get_3d_matrix();
        let v = TetElement::compute_volume(&x);

        let k = v * b.transpose() * c * b;

        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (k[(i, j)] - k[(j, i)]).abs() < 1e-6,
                    "K should be symmetric: K[{},{}]={} != K[{},{}]={}",
                    i,
                    j,
                    k[(i, j)],
                    j,
                    i,
                    k[(j, i)]
                );
            }
        }
    }

    #[test]
    fn test_stiffness_matrix_positive_semidefinite() {
        let material = Material::new(200e9, 0.3, 7800.0);

        let x = create_reference_tet();
        let b = TetElement::compute_b_matrix(&x);
        let c = material.get_3d_matrix();
        let v = TetElement::compute_volume(&x);

        let k = v * b.transpose() * c * b;

        let k_dyn = DMatrix::from_fn(12, 12, |i, j| k[(i, j)]);
        let eigen = k_dyn.symmetric_eigen();

        let max_ev = eigen
            .eigenvalues
            .iter()
            .fold(0.0_f64, |a, &b| a.max(b.abs()));
        let tol = max_ev * 1e-10;

        let mut zero_count = 0;
        let mut positive_count = 0;

        for &ev in eigen.eigenvalues.iter() {
            if ev.abs() < tol {
                zero_count += 1;
            } else if ev > -tol {
                positive_count += 1;
            } else {
                panic!("Found negative eigenvalue: {} (tolerance: {})", ev, tol);
            }
        }

        assert_eq!(
            zero_count, 6,
            "Expected 6 zero eigenvalues (rigid body modes), got {}",
            zero_count
        );
        assert_eq!(
            positive_count, 6,
            "Expected 6 positive eigenvalues, got {}",
            positive_count
        );
    }

    #[test]
    fn test_mass_matrix() {
        let material = Material::new(200e9, 0.3, 1000.0);

        let x = create_reference_tet();
        let v = TetElement::compute_volume(&x);

        let density = material.density;
        let diag = density * v / 10.0;
        let off_diag = density * v / 20.0;

        let expected_mass = density * v;
        let computed_mass = 4.0 * diag + 12.0 * off_diag;

        assert!(
            (computed_mass - expected_mass).abs() < 1e-6,
            "Total mass should be {}, got {}",
            expected_mass,
            computed_mass
        );
    }
}
