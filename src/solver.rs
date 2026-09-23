//! Sparse direct solvers for Kx = F. Primary: faer (pure Rust, SIMD-optimized, Cholesky/LU).
//! Falls back to rsparse on WASM where faer is unavailable (~10% faster than nalgebra-sparse).

use std::collections::HashMap;
use nalgebra as na;
use web_time::Instant;
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
/// Only available on native builds. WASM uses nalgebra fallback.
#[cfg(feature = "native")]
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

/// Solve using rsparse Cholesky. Used on WASM where faer is unavailable.
/// rsparse is ~10% faster than nalgebra-sparse.
#[cfg(not(feature = "native"))]
pub fn direct_solve(global_stiffness_matrix: &HashMap<(usize, usize), f64>, global_force_vector: &Vec<f64>) -> Vec<f64> {
    direct_cholesky_rsparse(global_stiffness_matrix, global_force_vector.clone())
}

/// Fast solver using rsparse Cholesky (pure Rust, WASM-compatible)
pub fn direct_cholesky_rsparse(global_stiffness_matrix: &HashMap<(usize, usize), f64>, global_force: Vec<f64>) -> Vec<f64> {
    use rsparse::data::Trpl;
    
    let (max_row, _max_col) = get_max_row_col(global_stiffness_matrix);
    let neq = max_row + 1;
    
    info!("Solving system of size: {} using rsparse Cholesky", neq);
    let start = Instant::now();
    
    // Build triplet matrix
    let mut trpl: Trpl<f64> = Trpl::new();
    for ((i, j), value) in global_stiffness_matrix.iter() {
        let v = if i == j { *value + 0.0001 } else { *value };
        trpl.append(*i, *j, v);
    }
    
    // Convert to CSC and solve
    let mat = trpl.to_sprs();
    let mut b = global_force;
    rsparse::cholsol(&mat, &mut b, 0).expect("rsparse cholsol failed");
    
    let duration = start.elapsed();
    info!("rsparse Cholesky completed in {:?}", duration);
    
    b
}

/// Legacy fallback solver using nalgebra-sparse Cholesky
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


/// Error type for solver failures
#[derive(Debug)]
pub enum SolverError {
    FactorizationFailed(String),
    InvalidMatrix(String),
}

/// Solve linear system from COO triplet format (row, col, val arrays).
/// This is a convenience function for GUI and external callers.
#[cfg(feature = "native")]
pub fn direct_solve_triplet(
    n: usize,
    rows: &[usize],
    cols: &[usize],
    vals: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, SolverError> {
    use faer::sparse::{SparseColMat, Triplet};
    use faer::prelude::*;
    use faer::linalg::solvers::Solve;
    
    info!("Solving system of size: {} using faer sparse solver (triplet input)", n);
    let start = Instant::now();
    
    // Convert to faer triplets
    let triplets: Vec<Triplet<usize, usize, f64>> = rows.iter()
        .zip(cols.iter())
        .zip(vals.iter())
        .map(|((&r, &c), &v)| Triplet::new(r, c, v))
        .collect();
    
    // Create sparse matrix from triplets
    let mat = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets)
        .map_err(|e| SolverError::InvalidMatrix(format!("{:?}", e)))?;
    
    // Create right-hand side as a column matrix
    let b = faer::Mat::<f64>::from_fn(n, 1, |i, _j| rhs[i]);
    
    // Try Cholesky first (faster for SPD matrices)
    let x = match mat.sp_cholesky(faer::Side::Lower) {
        Ok(llt) => {
            debug!("Using Cholesky factorization");
            llt.solve(&b)
        }
        Err(e) => {
            // Fall back to LU for non-SPD matrices
            info!("Cholesky failed ({:?}), falling back to LU factorization", e);
            let lu = mat.sp_lu().map_err(|e| SolverError::FactorizationFailed(format!("{:?}", e)))?;
            lu.solve(&b)
        }
    };
    
    let duration = start.elapsed();
    info!("faer solver completed in {:?}", duration);
    
    // Convert result to Vec<f64>
    Ok((0..n).map(|i| x[(i, 0)]).collect())
}

/// Solve linear system from COO triplet format - WASM version using rsparse.
#[cfg(not(feature = "native"))]
pub fn direct_solve_triplet(
    n: usize,
    rows: &[usize],
    cols: &[usize],
    vals: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, SolverError> {
    use rsparse::data::Trpl;
    
    info!("Solving system of size: {} using rsparse Cholesky (triplet input)", n);
    let start = Instant::now();
    
    // Build triplet matrix
    let mut trpl: Trpl<f64> = Trpl::new();
    for ((&r, &c), &v) in rows.iter().zip(cols.iter()).zip(vals.iter()) {
        let value = if r == c { v + 0.0001 } else { v };
        trpl.append(r, c, value);
    }
    
    // Convert to CSC and solve
    let mat = trpl.to_sprs();
    let mut b = rhs.to_vec();
    rsparse::cholsol(&mat, &mut b, 0)
        .map_err(|e| SolverError::FactorizationFailed(format!("{:?}", e)))?;
    
    let duration = start.elapsed();
    info!("rsparse Cholesky completed in {:?}", duration);
    
    Ok(b)
}
