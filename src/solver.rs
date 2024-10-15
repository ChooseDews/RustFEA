//takes in let mut global_stiffness_matrix: HashMap<(u32, u32), f64> = HashMap::new();
use russell_lab::{Matrix, Vector};
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, StrError, LinSolKind};
use std::collections::HashMap;
use nalgebra as na;
use nalgebra_sparse;
use std::time::Instant;
use log::{info, debug, trace};

pub fn get_max_row_col(global_stiffness_matrix: &HashMap<(u32, u32), f64>) -> (u32, u32) {
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

pub fn direct_solve(global_stiffness_matrix: &HashMap<(u32, u32), f64>, global_force_vector: &Vec<f64>) -> Vec<f64> {
    let (max_row, max_col) = get_max_row_col(global_stiffness_matrix);
    debug!("Allocating matrix of size: {}x{}", max_row+1, max_col+1);
    
    //find max row and col
    let neq = max_row + 1; // number of equations
    let nnz = global_stiffness_matrix.len();
    let mut trip = SparseTriplet::new(neq as usize, nnz as usize).unwrap();
    for key in global_stiffness_matrix.keys() {
        let row = key.0;
        let col = key.1;
        let value = global_stiffness_matrix[&key];
        trip.put(row as usize, col as usize, value).unwrap();
    }    
    //allocate rhs
    let rhs = Vector::from(global_force_vector);
    //calculate solution
    let mut config = ConfigSolver::new();
    config.lin_sol_kind(LinSolKind::Mmp);
    info!("Solving system of size: {}", neq);
    let start = Instant::now();
    let (mut solver, x) = Solver::compute(config, &trip, &rhs).unwrap();
    let duration = start.elapsed();
    info!("Solved System: Time elapsed in direct_solve() is: {:?}", duration);
    x.as_data().to_vec()
}

pub fn direct_choslky(global_stiffness_matrix: &HashMap<(u32, u32), f64>, global_force: Vec<f64>) -> Vec<f64> {
    let (max_row, max_col) = get_max_row_col(global_stiffness_matrix);
    debug!("Setting up sparse matrix for Cholesky decomposition");
    
    let mut sparse_matrix: nalgebra_sparse::CooMatrix<f64> = nalgebra_sparse::CooMatrix::new(max_row as usize +1, max_col as usize +1);
    for ((i, j), value) in global_stiffness_matrix.iter(){
        let mut v = *value;
        if i == j{
            v += 0.0001; //add small value to diagonal to avoid singularity
        }
        sparse_matrix.push(*i as usize, *j as usize, v);
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