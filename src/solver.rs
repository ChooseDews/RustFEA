//! Sparse direct solvers for Kx = F. Primary: faer (pure Rust, SIMD-optimized, Cholesky/LU).

use std::collections::HashMap;
use nalgebra as na;
use std::time::Instant;
use log::{info, debug, trace};

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

/// Solve using faer sparse Cholesky (LLT). Falls back to LU if Cholesky fails.
pub fn direct_solve(global_stiffness_matrix: &HashMap<(usize, usize), f64>, global_force_vector: &Vec<f64>) -> Vec<f64> {
    use faer::sparse::{SparseColMat, Triplet};
    use faer::prelude::*;
    use faer::linalg::solvers::Solve;
    
    let (max_row, _max_col) = get_max_row_col(global_stiffness_matrix);
    let neq = max_row + 1;
    
    info!("Solving system of size: {} using faer sparse solver", neq);
    let start = Instant::now();
    
    // Convert HashMap to triplet format using faer's Triplet struct
    let triplets: Vec<Triplet<usize, usize, f64>> = global_stiffness_matrix
        .iter()
        .map(|((r, c), v)| Triplet::new(*r, *c, *v))
        .collect();
    
    // Create sparse matrix from triplets
    let mat = SparseColMat::<usize, f64>::try_new_from_triplets(
        neq, neq, 
        &triplets
    ).expect("Failed to create sparse matrix");
    
    // Create right-hand side as a column matrix
    let b = faer::Mat::<f64>::from_fn(neq, 1, |i, _j| global_force_vector[i]);
    
    // Try Cholesky first (faster for SPD matrices)
    let x = match mat.sp_cholesky(faer::Side::Lower) {
        Ok(llt) => {
            debug!("Using Cholesky factorization");
            llt.solve(&b)
        }
        Err(e) => {
            // Fall back to LU for non-SPD matrices
            info!("Cholesky failed ({:?}), falling back to LU factorization", e);
            let lu = mat.sp_lu().expect("LU factorization failed");
            lu.solve(&b)
        }
    };
    
    let duration = start.elapsed();
    info!("faer solver completed in {:?}", duration);
    
    // Convert result to Vec<f64>
    (0..neq).map(|i| x[(i, 0)]).collect()
}

/// Fallback solver using nalgebra-sparse Cholesky
pub fn direct_cholesky_nalgebra(global_stiffness_matrix: &HashMap<(usize, usize), f64>, global_force: Vec<f64>) -> Vec<f64> {
    let (max_row, max_col) = get_max_row_col(global_stiffness_matrix);
    
    let mut sparse_matrix: nalgebra_sparse::CooMatrix<f64> = nalgebra_sparse::CooMatrix::new(max_row + 1, max_col + 1);
    for ((i, j), value) in global_stiffness_matrix.iter() {
        let mut v = *value;
        if i == j {
            v += 0.0001;
        }
        sparse_matrix.push(*i, *j, v);
    }
    let csc = nalgebra_sparse::CscMatrix::from(&sparse_matrix);

    info!("Solving system of size: {} using nalgebra Cholesky decomposition", max_row + 1);
    let b = nalgebra::DVector::from_vec(global_force);
    let start = Instant::now();
    let cholesky = nalgebra_sparse::factorization::CscCholesky::factor(&csc).unwrap();
    let u: na::Matrix<f64, na::Dyn, na::Dyn, na::VecStorage<f64, na::Dyn, na::Dyn>> = cholesky.solve(&b);
    let duration = start.elapsed();
    info!("nalgebra Cholesky decomposition completed in {:?}", duration);

    u.data.as_vec().clone()
}
