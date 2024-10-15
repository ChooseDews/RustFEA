extern crate nalgebra as na;
use na::{DVector, DMatrix};
use nalgebra::{Matrix3, Matrix6};
use crate::simulation::Simulation;
use std::collections::HashMap;
use serde::{Serialize, Deserialize};
use typetag;
use log::{debug, trace};


#[derive(Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum ElementType {
    Brick,
    Quad
}


#[derive(Debug, Serialize, Deserialize)]
pub struct ElementFields { //hashmap containing
    pub field: HashMap<String, Vec<f64>>,
    size: usize,
    pub connectivity: Vec<usize>
}


impl ElementFields {
    pub fn new(connectivity: Vec<usize>) -> Self {
        let size = connectivity.len();
        ElementFields {
            field: HashMap::new(),
            size,
            connectivity
        }
    }

    pub fn add_field(&mut self, name: &str, field: Vec<f64>) {
        self.field.insert(name.to_string(), field);
    }

    pub fn find_field(&mut self, name: &str) -> Option<&mut Vec<f64>> {
        //check if field exists if not create it
        if !self.field.contains_key(name) {
            self.field.insert(name.to_string(), vec![0.0; self.size]);
        }
        self.field.get_mut(name)
    }

    pub fn append_to_feild(&mut self, name: &str, nn: usize, value: f64) {
        let current_field = self.find_field(name).unwrap();
        current_field[nn] = value;
    }

    pub fn get_feild_names(&self) -> Vec<String> {
        let mut names = Vec::new();
        for key in self.field.keys() {
            names.push(key.to_string());
        }
        names
    }


}

#[derive(Serialize, Deserialize, Debug)]
pub struct Material {
    //linear elastic material properties
    pub youngs_modulus: f64,
    pub poisson_ratio: f64,
    pub density: f64
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

    pub fn get_3d_matrix(&self) -> Matrix6<f64>{
        let mut C = Matrix6::<f64>::zeros();
        let E = self.youngs_modulus;
        let v = self.poisson_ratio;
        let _g = E / (2.0 * (1.0 + v));
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

    pub fn get_density(&self) -> f64 {
        self.density
    }

    pub fn get_wave_speed(&self) -> f64 {
        (self.youngs_modulus / self.density).sqrt()
    }

    //Aluminum 6061-T6
    pub fn aluminum() -> Self {
        Material {
            youngs_modulus: 68.9e9,
            poisson_ratio: 0.33,
            density: 2700.0 //mass scaling
        }
    }

    pub fn empty() -> Self {
        Material {
            youngs_modulus: 1.0,
            poisson_ratio: 0.2,
            density: 1.0
        }
    }

}


#[typetag::serde]
pub trait BaseElement {

    fn dofs(&self) -> usize { 3 }
    fn get_id(&self) -> usize;
    fn get_connectivity(&self) -> &Vec<usize>;
    fn get_material(&self) -> &Material;
    fn get_deformation_gradient(&self) -> &DMatrix<f64>;
    fn get_shape_functions(&self, xi: f64, eta: f64, zeta: f64) -> DVector<f64>;
    fn get_shape_derivatives(&self, xi: f64, eta: f64, zeta: f64) -> DMatrix<f64>;
    fn get_global_position(&self, n: &DVector<f64>, simulation: &Simulation) -> na::Vector3<f64> {
        na::Vector3::zeros()
    }

    fn get_signed_distance_vector(&self, point: na::Vector3<f64>, simulation: &Simulation) -> Option<na::Vector3<f64>> {
        None
    }

    fn is_active(&self) -> bool {false}
    fn set_active(&mut self, active: bool){}

    fn get_stiffness(&self) -> DMatrix<f64>;
    fn compute_stiffness(&mut self, simulation: &Simulation){}


    fn compute_mass(&self, simulation: &Simulation) -> DMatrix<f64>{
        DMatrix::<f64>::zeros(self.get_connectivity().len(), self.get_connectivity().len())
    }
    fn get_mass(&self) -> &DMatrix<f64>;
    fn set_mass(&mut self, mass: DMatrix<f64>){}
    fn set_lumped_mass(&mut self, mass: &DMatrix<f64>) -> f64;
    fn get_lumped_mass(&self) -> &Vec<f64>;

    fn get_x(&self, simulation: &Simulation) -> DMatrix<f64>;
    fn get_u(&self, simulation: &Simulation) -> DVector<f64>;
    fn get_b(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DMatrix<f64>; //get the strain displacement matrix
    // fn compute_volume(&self, simulation: &Simulation) -> f64;
    fn compute_element_nodal_properties(&self, simulation: &Simulation) -> ElementFields {
        ElementFields::new(self.get_connectivity().clone()) 
    }

    fn type_name(&self) -> ElementType;

    //explicit stuff
    fn compute_force(&self, simulation: &Simulation) -> DVector<f64> {
        DVector::zeros(self.get_connectivity().len() * self.dofs())
    }
    // fn step_explicit(&mut self, simulation: &Simulation, force: DVector<f64>, dt: f64); //self modification


}