//takes in let mut global_stiffness_matrix: HashMap<(usize, usize), f64> = HashMap::new();
use russell_lab::{Matrix, Vector};
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, StrError, LinSolKind};
use std::collections::HashMap;
use nalgebra as na;
use nalgebra_sparse;
use std::time::Instant;


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

pub fn direct_solve(global_stiffness_matrix: &HashMap<(usize, usize), f64>, global_force_vector: &Vec<f64>) -> Vec<f64> {
    //find max row and col
    let (max_row, max_col) = get_max_row_col(global_stiffness_matrix);
    println!("Allocating matrix of size: {}x{}", max_row+1, max_col+1);
    //allocate a square matrix
    let neq = max_row + 1; // number of equations
    let nnz = global_stiffness_matrix.len();
    let mut trip = SparseTriplet::new(neq, nnz).unwrap();
    for key in global_stiffness_matrix.keys() {
        let row = key.0;
        let col = key.1;
        let value = global_stiffness_matrix[&key];
        trip.put(row, col, value).unwrap();
    }    
    //allocate rhs
    let rhs = Vector::from(global_force_vector);
    //calculate solution
    let mut config = ConfigSolver::new();
    config.lin_sol_kind(LinSolKind::Mmp);
    println!("direct_solve() Solving system of size: {}", neq);
    let start = Instant::now();
    let (mut solver, x) = Solver::compute(config, &trip, &rhs).unwrap();
    let duration = start.elapsed();
    println!("Solved System: Time elapsed in direct_solve() is: {:?}", duration);
    x.as_data().to_vec()
}


pub fn direct_choslky(global_stiffness_matrix: &HashMap<(usize, usize), f64>, global_force: Vec<f64>) -> Vec<f64> {
    let (max_row, max_col) = get_max_row_col(global_stiffness_matrix);
    let mut sparse_matrix: nalgebra_sparse::CooMatrix<f64> = nalgebra_sparse::CooMatrix::new(max_row+1, max_col+1);
    for ((i, j), value) in global_stiffness_matrix.iter(){
        let mut v = *value;
        if i == j{
            v += 0.0001; //add small value to diagonal to avoid singularity
        }
        sparse_matrix.push(*i, *j, v);
    }
    let csc = nalgebra_sparse::CscMatrix::from(&sparse_matrix);

    println!("direct_choslky() Solving system of size: {}", max_row+1);
    let b = nalgebra::DVector::from_vec(global_force);
    let start = Instant::now();
    let cholesky = nalgebra_sparse::factorization::CscCholesky::factor(&csc).unwrap();
    let u: na::Matrix<f64, na::Dyn, na::Dyn, na::VecStorage<f64, na::Dyn, na::Dyn>> = cholesky.solve(&b);
    let duration = start.elapsed();
    println!("Solved System: Time elapsed in direct_choslky() is: {:?}", duration);

    return u.data.as_vec().clone();
}