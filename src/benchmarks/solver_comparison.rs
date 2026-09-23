//! Benchmark comparing faer vs nalgebra-sparse vs rsparse direct solvers.
//!
//! Run with: cargo run --bin solver_bench --release

use web_time::Instant;

/// Get current process memory usage in MB (macOS specific)
#[cfg(target_os = "macos")]
fn get_memory_mb() -> f64 {
    use std::mem::MaybeUninit;

    // mach_task_basic_info struct
    #[repr(C)]
    struct MachTaskBasicInfo {
        virtual_size: u64,
        resident_size: u64,
        resident_size_max: u64,
        user_time: [u32; 2],
        system_time: [u32; 2],
        policy: i32,
        suspend_count: i32,
    }

    extern "C" {
        fn mach_task_self() -> u32;
        fn task_info(
            target_task: u32,
            flavor: i32,
            task_info_out: *mut MachTaskBasicInfo,
            task_info_count: *mut u32,
        ) -> i32;
    }

    const MACH_TASK_BASIC_INFO: i32 = 20;
    const MACH_TASK_BASIC_INFO_COUNT: u32 = 10;

    unsafe {
        let mut info = MaybeUninit::<MachTaskBasicInfo>::uninit();
        let mut count = MACH_TASK_BASIC_INFO_COUNT;
        let result = task_info(
            mach_task_self(),
            MACH_TASK_BASIC_INFO,
            info.as_mut_ptr(),
            &mut count,
        );
        if result == 0 {
            let info = info.assume_init();
            info.resident_size as f64 / (1024.0 * 1024.0)
        } else {
            0.0
        }
    }
}

#[cfg(not(target_os = "macos"))]
fn get_memory_mb() -> f64 {
    0.0 // Not implemented for other platforms
}

/// Result of a solver benchmark run
#[derive(Debug, Clone)]
pub struct SolverBenchmarkResult {
    pub n_dofs: usize,
    pub n_nonzeros: usize,
    pub faer_time_ms: f64,
    pub nalgebra_time_ms: f64,
    pub rsparse_time_ms: f64,
    pub faer_mem_mb: f64,
    pub nalgebra_mem_mb: f64,
    pub rsparse_mem_mb: f64,
    pub faer_vs_nalgebra: f64,
    pub rsparse_vs_nalgebra: f64,
    pub solution_diff_norm: f64,
}

/// Solve using faer sparse Cholesky/LU
#[cfg(feature = "native")]
fn solve_faer(n: usize, triplets: &[(usize, usize, f64)], rhs: &[f64]) -> (Vec<f64>, f64, f64) {
    use faer::linalg::solvers::Solve;
    use faer::prelude::*;
    use faer::sparse::{SparseColMat, Triplet};

    let mem_before = get_memory_mb();
    let start = Instant::now();

    let faer_triplets: Vec<Triplet<usize, usize, f64>> = triplets
        .iter()
        .map(|(r, c, v)| Triplet::new(*r, *c, *v))
        .collect();

    let mat = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &faer_triplets)
        .expect("Failed to create sparse matrix");

    let b = faer::Mat::<f64>::from_fn(n, 1, |i, _j| rhs[i]);

    let x = match mat.sp_cholesky(faer::Side::Lower) {
        Ok(llt) => llt.solve(&b),
        Err(_) => {
            let lu = mat.sp_lu().expect("LU factorization failed");
            lu.solve(&b)
        }
    };

    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
    let mem_after = get_memory_mb();
    let result: Vec<f64> = (0..n).map(|i| x[(i, 0)]).collect();

    (result, elapsed_ms, mem_after - mem_before)
}

/// Solve using nalgebra-sparse Cholesky
fn solve_nalgebra(n: usize, triplets: &[(usize, usize, f64)], rhs: &[f64]) -> (Vec<f64>, f64, f64) {
    use nalgebra::DVector;
    use nalgebra_sparse::factorization::CscCholesky;
    use nalgebra_sparse::CooMatrix;
    use nalgebra_sparse::CscMatrix;

    let mem_before = get_memory_mb();
    let start = Instant::now();

    let mut coo: CooMatrix<f64> = CooMatrix::new(n, n);
    for &(r, c, v) in triplets {
        let value = if r == c { v + 0.0001 } else { v };
        coo.push(r, c, value);
    }
    let csc = CscMatrix::from(&coo);

    let b = DVector::from_vec(rhs.to_vec());
    let cholesky = CscCholesky::factor(&csc).expect("Cholesky factorization failed");
    let u = cholesky.solve(&b);

    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
    let mem_after = get_memory_mb();

    (u.data.as_vec().clone(), elapsed_ms, mem_after - mem_before)
}

/// Solve using rsparse Cholesky
fn solve_rsparse(n: usize, triplets: &[(usize, usize, f64)], rhs: &[f64]) -> (Vec<f64>, f64, f64) {
    use rsparse::data::Trpl;

    let mem_before = get_memory_mb();
    let start = Instant::now();

    let mut trpl: Trpl<f64> = Trpl::new();
    for &(r, c, v) in triplets {
        let value = if r == c { v + 0.0001 } else { v };
        trpl.append(r, c, value);
    }

    let mat = trpl.to_sprs();
    let mut b: Vec<f64> = rhs.to_vec();
    rsparse::cholsol(&mat, &mut b, 0).expect("rsparse cholsol failed");

    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
    let mem_after = get_memory_mb();

    (b, elapsed_ms, mem_after - mem_before)
}

/// Generate a test stiffness matrix (3D elasticity-like sparsity pattern)
fn generate_test_system(n_nodes: usize) -> (usize, Vec<(usize, usize, f64)>, Vec<f64>) {
    let n_dofs = n_nodes * 3;
    let mut triplets = Vec::new();
    let bandwidth = 30.min(n_nodes / 2);

    for i in 0..n_nodes {
        for dof_i in 0..3 {
            let row = i * 3 + dof_i;
            triplets.push((row, row, 1000.0 + (i as f64) * 0.1));

            for j in (i.saturating_sub(bandwidth))..=(i + bandwidth).min(n_nodes - 1) {
                if i != j {
                    for dof_j in 0..3 {
                        let col = j * 3 + dof_j;
                        if col > row {
                            let coupling = -10.0 / ((1 + (i as i64 - j as i64).abs()) as f64);
                            triplets.push((row, col, coupling));
                            triplets.push((col, row, coupling));
                        }
                    }
                }
            }
        }
    }

    let rhs: Vec<f64> = (0..n_dofs)
        .map(|i| (i as f64 * 0.01).sin() * 100.0)
        .collect();
    (n_dofs, triplets, rhs)
}

/// Run the solver comparison benchmark
#[cfg(feature = "native")]
pub fn run_solver_comparison() -> Vec<SolverBenchmarkResult> {
    println!("\n{}", "=".repeat(95));
    println!("SOLVER COMPARISON: faer vs nalgebra-sparse vs rsparse (time + memory)");
    println!("{}\n", "=".repeat(95));

    let node_counts = vec![100, 500, 1000, 2000, 5000, 10000, 20000];
    let mut results = Vec::new();

    for n_nodes in node_counts {
        let (n_dofs, triplets, rhs) = generate_test_system(n_nodes);
        let n_nonzeros = triplets.len();

        println!(
            "Testing system: {} nodes, {} DOFs, {} nonzeros",
            n_nodes, n_dofs, n_nonzeros
        );

        // Warmup
        let _ = solve_faer(n_dofs, &triplets, &rhs);
        let _ = solve_nalgebra(n_dofs, &triplets, &rhs);
        let _ = solve_rsparse(n_dofs, &triplets, &rhs);

        // Measured run
        let (faer_solution, faer_time, faer_mem) = solve_faer(n_dofs, &triplets, &rhs);
        let (_, nalgebra_time, nalgebra_mem) = solve_nalgebra(n_dofs, &triplets, &rhs);
        let (rsparse_solution, rsparse_time, rsparse_mem) = solve_rsparse(n_dofs, &triplets, &rhs);

        let faer_vs_nalgebra = nalgebra_time / faer_time;
        let rsparse_vs_nalgebra = nalgebra_time / rsparse_time;

        let diff_norm: f64 = faer_solution
            .iter()
            .zip(rsparse_solution.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let rel_diff = diff_norm / faer_solution.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();

        println!("  faer:     {:7.2} ms, {:+6.1} MB", faer_time, faer_mem);
        println!(
            "  nalgebra: {:7.2} ms, {:+6.1} MB",
            nalgebra_time, nalgebra_mem
        );
        println!(
            "  rsparse:  {:7.2} ms, {:+6.1} MB",
            rsparse_time, rsparse_mem
        );
        println!(
            "  speedup: faer {:.2}x, rsparse {:.2}x vs nalgebra\n",
            faer_vs_nalgebra, rsparse_vs_nalgebra
        );

        results.push(SolverBenchmarkResult {
            n_dofs,
            n_nonzeros,
            faer_time_ms: faer_time,
            nalgebra_time_ms: nalgebra_time,
            rsparse_time_ms: rsparse_time,
            faer_mem_mb: faer_mem,
            nalgebra_mem_mb: nalgebra_mem,
            rsparse_mem_mb: rsparse_mem,
            faer_vs_nalgebra,
            rsparse_vs_nalgebra,
            solution_diff_norm: rel_diff,
        });
    }

    // Summary tables
    println!("\n{}", "=".repeat(95));
    println!("TIME COMPARISON (ms)");
    println!("{}", "=".repeat(95));
    println!(
        "{:>8} {:>12} {:>12} {:>12} {:>12} {:>14}",
        "DOFs", "faer", "nalgebra", "rsparse", "faer/na", "rsparse/na"
    );
    println!("{}", "-".repeat(75));

    for r in &results {
        println!(
            "{:>8} {:>11.2}ms {:>11.2}ms {:>11.2}ms {:>11.2}x {:>13.2}x",
            r.n_dofs,
            r.faer_time_ms,
            r.nalgebra_time_ms,
            r.rsparse_time_ms,
            r.faer_vs_nalgebra,
            r.rsparse_vs_nalgebra
        );
    }

    println!("\n{}", "=".repeat(95));
    println!("MEMORY DELTA (MB) - change in RSS during solve");
    println!("{}", "=".repeat(95));
    println!(
        "{:>8} {:>12} {:>12} {:>12}",
        "DOFs", "faer", "nalgebra", "rsparse"
    );
    println!("{}", "-".repeat(50));

    for r in &results {
        println!(
            "{:>8} {:>+11.1}MB {:>+11.1}MB {:>+11.1}MB",
            r.n_dofs, r.faer_mem_mb, r.nalgebra_mem_mb, r.rsparse_mem_mb
        );
    }

    // Averages
    let avg_faer_speed: f64 =
        results.iter().map(|r| r.faer_vs_nalgebra).sum::<f64>() / results.len() as f64;
    let avg_rsparse_speed: f64 =
        results.iter().map(|r| r.rsparse_vs_nalgebra).sum::<f64>() / results.len() as f64;
    let total_faer_mem: f64 = results.iter().map(|r| r.faer_mem_mb).sum();
    let total_nalgebra_mem: f64 = results.iter().map(|r| r.nalgebra_mem_mb).sum();
    let total_rsparse_mem: f64 = results.iter().map(|r| r.rsparse_mem_mb).sum();

    println!("\n{}", "=".repeat(95));
    println!("SUMMARY");
    println!("{}", "=".repeat(95));
    println!("Average speedup vs nalgebra:");
    println!("  faer:    {:.2}x faster", avg_faer_speed);
    println!("  rsparse: {:.2}x faster", avg_rsparse_speed);
    println!("\nTotal memory delta across all tests:");
    println!("  faer:    {:+.1} MB", total_faer_mem);
    println!("  nalgebra:{:+.1} MB", total_nalgebra_mem);
    println!("  rsparse: {:+.1} MB", total_rsparse_mem);

    println!("\nFor WASM (rsparse vs nalgebra):");
    if avg_rsparse_speed > 1.0 {
        println!(
            "  ✅ rsparse is {:.0}% FASTER",
            (avg_rsparse_speed - 1.0) * 100.0
        );
    } else {
        println!(
            "  ❌ rsparse is {:.0}% slower",
            (1.0 - avg_rsparse_speed) * 100.0
        );
    }

    results
}

#[cfg(not(feature = "native"))]
pub fn run_solver_comparison() -> Vec<SolverBenchmarkResult> {
    println!("Solver comparison requires native feature (faer not available on WASM)");
    Vec::new()
}
