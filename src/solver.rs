//takes in let mut global_stiffness_matrix: HashMap<(usize, usize), f64> = HashMap::new();
use russell_lab::{Matrix, Vector};
use russell_sparse::prelude::*;
use std::collections::HashMap;
use nalgebra as na;
use std::time::Instant;
use log::{info, debug, trace};
use crate::simulation::Simulation;

pub fn get_max_row_col(global_stiffness_matrix: &HashMap<(usize, usize), f64>) -> (usize, usize) {
    let mut max_row = 0;
    let mut max_col = 0;
    for key in global_stiffness_matrix.keys() {
        if key.0 > max_row {
            max_row = key.0;
        }
        if key.1 > max_col {
            max_col = key.1;
        }
    }
    assert!(max_row == max_col, "Matrix is not square");
    (max_row, max_col)
}

pub fn direct_solve(simulation: &Simulation, global_stiffness_matrix: &HashMap<(usize, usize), f64>, global_force_vector: &Vec<f64>) -> Vec<f64> {
    let (max_row, max_col) = get_max_row_col(global_stiffness_matrix);
    debug!("Allocating matrix of size: {}x{}", max_row+1, max_col+1);
    let neq = max_row + 1; // number of equations
    let nnz = global_stiffness_matrix.len();
    let mut umfpack = SolverUMFPACK::new().unwrap();
    let mut coo = SparseMatrix::new_coo(neq, neq, nnz, Sym::YesFull).unwrap();
    for ((row, col), value) in global_stiffness_matrix.iter() {
        coo.put(*row, *col, *value).unwrap();
    } 
    let save_matrix = simulation.keywords.get_keyword("SOLVER_SAVE_STIFFNESS_MATRIX_PATH");
    if save_matrix.is_some() {
        let start = Instant::now();
        let coo = coo.get_coo().unwrap();
        let path = save_matrix.unwrap().value.as_str().unwrap();
        let csc = CscMatrix::from_coo(&coo).unwrap();
        csc.write_matrix_market(path, true).unwrap();
        let duration = start.elapsed();
        info!("Saved matrix to {} in {:?}", path, duration);
    }

    let rhs = Vector::from(global_force_vector);
    let mut x = Vector::new(neq);
    info!("Solving system of size: {}", neq);
    let start = Instant::now();
    umfpack.factorize(&mut coo, None).unwrap();
    umfpack.solve(&mut x, &coo, &rhs, false).unwrap();
    let duration = start.elapsed();
    info!("Solved System: Time elapsed in direct_solve() is: {:?}", duration);
    x.as_data().to_vec()
}

pub fn direct_choslky(global_stiffness_matrix: &HashMap<(usize, usize), f64>, global_force: Vec<f64>) -> Vec<f64> {
    let (max_row, max_col) = get_max_row_col(global_stiffness_matrix);
    debug!("Setting up sparse matrix for Cholesky decomposition");
    
    let mut sparse_matrix: nalgebra_sparse::CooMatrix<f64> = nalgebra_sparse::CooMatrix::new(max_row  +1, max_col  +1);
    for ((i, j), value) in global_stiffness_matrix.iter(){
        let mut v = *value;
        if i == j{
            v += 0.0001; //add small value to diagonal to avoid singularity
        }
        sparse_matrix.push(*i , *j , v);
    }
    let csc = nalgebra_sparse::CscMatrix::from(&sparse_matrix);

    info!("Solving system of size: {} using Cholesky decomposition", max_row+1);
    let b = nalgebra::DVector::from_vec(global_force);
    let start = Instant::now();
    let cholesky = nalgebra_sparse::factorization::CscCholesky::factor(&csc).unwrap();
    let u: na::Matrix<f64, na::Dyn, na::Dyn, na::VecStorage<f64, na::Dyn, na::Dyn>> = cholesky.solve(&b);
    let duration = start.elapsed();
    info!("Cholesky decomposition completed in {:?}", duration);

    return u.data.as_vec().clone();
}