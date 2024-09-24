use nalgebra as na;
use na::{DVector};


pub fn compute_von_mises(stress_vector: DVector<f64>) -> f64 {
    //compute and return the von mises stress here
    let s11 = stress_vector[0];
    let s22 = stress_vector[1];
    let s33 = stress_vector[2];
    let s12 = stress_vector[3];
    let s23 = stress_vector[4];
    let s13 = stress_vector[5];
    let s_vm = ((s11 - s22).powi(2) + (s22 - s33).powi(2) + (s33 - s11).powi(2) + 6.0 * (s12.powi(2) + s23.powi(2) + s13.powi(2))) / 2.0;
    s_vm.sqrt()
}

//compute the y displacement of a cantilever beam with a point load at the end
pub fn eular_beam_displacement(x: f64, length: f64, load: f64, modulus: f64, moment_of_inertia: f64) -> f64 {
    let y = load * x.powi(2) * (3.0 * length - x) / (6.0 * modulus * moment_of_inertia);
    y
}