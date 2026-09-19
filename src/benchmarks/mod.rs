// Benchmarks module for FEA accuracy validation
// Contains analytical benchmarks comparing FEA results to exact solutions

pub mod uniaxial_tension;
pub mod pure_shear;
pub mod hydrostatic_compression;
pub mod torsion_shaft;
pub mod torsion_explicit;
pub mod hollow_sphere;
pub mod cantilever_beam;
pub mod spherical_cavity;
pub mod boussinesq;
pub mod hertz_sphere_flat;
pub mod hertz_sphere_sphere;
pub mod contact_explicit;
pub mod gravity;
pub mod mesh_utils;
pub mod plotting;

use serde::{Deserialize, Serialize};
use std::time::Instant;

/// A single metric comparison between analytical and computed values
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetricComparison {
    pub name: String,
    pub analytical: f64,
    pub computed: f64,
    pub relative_error: f64,
    pub tolerance: f64,
}

impl MetricComparison {
    pub fn new(name: &str, analytical: f64, computed: f64, tolerance: f64) -> Self {
        let relative_error = if analytical.abs() > 1e-15 {
            (computed - analytical) / analytical.abs()
        } else {
            computed - analytical
        };
        
        MetricComparison {
            name: name.to_string(),
            analytical,
            computed,
            relative_error,
            tolerance,
        }
    }
    
    pub fn passed(&self) -> bool {
        self.relative_error.abs() <= self.tolerance
    }
}

/// Result of a single benchmark
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult {
    pub name: String,
    pub description: String,
    pub passed: bool,
    pub metrics: Vec<MetricComparison>,
    pub elapsed_ms: f64,
    pub notes: Option<String>,
}

impl BenchmarkResult {
    pub fn new(name: &str, description: &str) -> Self {
        BenchmarkResult {
            name: name.to_string(),
            description: description.to_string(),
            passed: true,
            metrics: Vec::new(),
            elapsed_ms: 0.0,
            notes: None,
        }
    }
    
    pub fn add_metric(&mut self, metric: MetricComparison) {
        if !metric.passed() {
            self.passed = false;
        }
        self.metrics.push(metric);
    }
    
    pub fn set_notes(&mut self, notes: &str) {
        self.notes = Some(notes.to_string());
    }
    
    pub fn set_elapsed(&mut self, elapsed_ms: f64) {
        self.elapsed_ms = elapsed_ms;
    }
}

/// A collection of benchmarks to run
pub struct BenchmarkSuite {
    benchmarks: Vec<(String, fn() -> BenchmarkResult)>,
}

impl BenchmarkSuite {
    pub fn new() -> Self {
        BenchmarkSuite {
            benchmarks: Vec::new(),
        }
    }
    
    pub fn add_benchmark(&mut self, name: &str, run_fn: fn() -> BenchmarkResult) {
        self.benchmarks.push((name.to_string(), run_fn));
    }
    
    pub fn benchmark_count(&self) -> usize {
        self.benchmarks.len()
    }
    
    pub fn run_all(&self) -> Vec<BenchmarkResult> {
        let mut results = Vec::new();
        
        for (name, run_fn) in &self.benchmarks {
            println!("Running benchmark: {}...", name);
            let start = Instant::now();
            let mut result = run_fn();
            result.set_elapsed(start.elapsed().as_secs_f64() * 1000.0);
            
            let status = if result.passed { "PASS" } else { "FAIL" };
            println!("  {} ({}ms)", status, result.elapsed_ms as u64);
            
            results.push(result);
        }
        
        results
    }
}

impl Default for BenchmarkSuite {
    fn default() -> Self {
        Self::new()
    }
}
