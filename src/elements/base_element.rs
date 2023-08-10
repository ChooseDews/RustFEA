extern crate nalgebra as na;
use na::{DVector, DMatrix};
use crate::simulation::Simulation;
pub struct Material {
    //linear elastic material properties
    youngs_modulus: f64,
    poisson_ratio: f64,
    density: f64
}

//3D [C] matrix
impl Material {
    pub fn new(youngs_modulus: f64, poisson_ratio: f64, density: f64) -> Self {
        Material {
            youngs_modulus,
            poisson_ratio,
            density
        }
    }

    pub fn get_3d_matrix(&self) -> DMatrix<f64>{
        let mut C = DMatrix::<f64>::zeros(6,6);
        let E = self.youngs_modulus;
        let v = self.poisson_ratio;
        let G = E / (2.0 * (1.0 + v));
        let lambda = (E * v) / ((1.0 + v) * (1.0 - 2.0 * v));
        let mu = E / (2.0 * (1.0 + v));
        C[(0,0)] = lambda + 2.0 * mu;
        C[(0,1)] = lambda;
        C[(0,2)] = lambda;
        C[(1,0)] = lambda;
        C[(1,1)] = lambda + 2.0 * mu;
        C[(1,2)] = lambda;
        C[(2,0)] = lambda;
        C[(2,1)] = lambda;
        C[(2,2)] = lambda + 2.0 * mu;
        C[(3,3)] = mu;
        C[(4,4)] = mu;
        C[(5,5)] = mu;
        C
    }

    //Aluminum 6061-T6
    pub fn aluminum() -> Self {
        Material {
            youngs_modulus: 68.9e9,
            poisson_ratio: 0.33,
            density: 2700.0
        }
    }

}

pub trait BaseElement {
    fn get_id(&self) -> usize;
    fn get_connectivity(&self) -> &Vec<usize>;
    fn get_material(&self) -> &Material;
    fn get_deformation_gradient(&self) -> &DMatrix<f64>;
    fn get_shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> DVector<f64>;
    fn get_shape_derivatives(&self, xi: f64, eta: f64, zeta: f64) -> DMatrix<f64>;
    fn get_global_position(&self, N: &DVector<f64>, simulation: &Simulation) -> na::Vector3<f64>;
    fn compute_stiffness(&self, simulation: &Simulation) -> DMatrix<f64>;
    fn compute_force(&self, simulation: &Simulation) -> DVector<f64>;
    fn compute_jacobian_matrix(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DMatrix<f64>;
    fn get_X(&self, simulation: &Simulation) -> DMatrix<f64>;
    fn compute_B(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DMatrix<f64>;
}