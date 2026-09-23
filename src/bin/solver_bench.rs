//! Simple benchmark comparing faer vs nalgebra-sparse vs rsparse direct solvers
//!
//! Run with: cargo run --bin solver_bench --release

use rust_fea::benchmarks::solver_comparison;

fn main() {
    env_logger::init();
    solver_comparison::run_solver_comparison();
}
