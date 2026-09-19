// Benchmark runner for FEA accuracy validation
// Runs a suite of analytical benchmarks and compares FEA results to exact solutions

use rust_fea::benchmarks::{
    BenchmarkResult, BenchmarkSuite,
    uniaxial_tension, pure_shear, hydrostatic_compression,
    torsion_shaft, cantilever_beam, spherical_cavity,
    boussinesq, hertz_sphere_flat, hertz_sphere_sphere,
    plotting::{format_with_units, ConvergenceStudy, generate_convergence_table},
};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use clap::{Parser, ValueEnum};
use log::info;

#[derive(Parser)]
#[command(name = "run_benchmarks")]
#[command(about = "Run FEA accuracy benchmarks and generate a report")]
struct Args {
    /// Output format for the report
    #[arg(short, long, default_value = "markdown")]
    format: OutputFormat,
    
    /// Output file path (stdout if not specified)
    #[arg(short, long)]
    output: Option<String>,
    
    /// Run only specific benchmarks (comma-separated)
    #[arg(short, long)]
    benchmarks: Option<String>,
    
    /// Verbose output
    #[arg(short, long, default_value = "false")]
    verbose: bool,
}

#[derive(Clone, ValueEnum)]
enum OutputFormat {
    Markdown,
    Json,
    Csv,
}

fn main() {
    env_logger::init();
    
    let args = Args::parse();
    
    let mut suite = BenchmarkSuite::new();
    
    // Register all benchmarks
    let filter = args.benchmarks.as_deref();
    
    let benchmarks_to_run: Vec<(&str, fn() -> BenchmarkResult)> = vec![
        // Level 1: Fundamental Element Verification
        ("uniaxial_tension", uniaxial_tension::run),
        ("pure_shear", pure_shear::run),
        ("hydrostatic_compression", hydrostatic_compression::run),
        // Level 2: 3D Continuum Verification  
        ("torsion_shaft", torsion_shaft::run),
        // ("hollow_sphere", hollow_sphere::run),  // TODO: Fix mesh generation
        ("cantilever_beam", cantilever_beam::run),
        ("spherical_cavity", spherical_cavity::run),
        ("boussinesq", boussinesq::run),
        // Level 3: Contact Verification (analytical validation)
        ("hertz_sphere_flat", hertz_sphere_flat::run),
        ("hertz_sphere_sphere", hertz_sphere_sphere::run),
    ];
    
    for (name, run_fn) in benchmarks_to_run.iter() {
        if let Some(filter_str) = filter {
            let filters: Vec<&str> = filter_str.split(',').collect();
            if !filters.iter().any(|f| name.contains(f)) {
                continue;
            }
        }
        suite.add_benchmark(name, *run_fn);
    }
    
    // Run all benchmarks
    info!("Running {} benchmarks...", suite.benchmark_count());
    let results = suite.run_all();
    
    // Generate report
    let report = match args.format {
        OutputFormat::Markdown => generate_markdown_report(&results),
        OutputFormat::Json => generate_json_report(&results),
        OutputFormat::Csv => generate_csv_report(&results),
    };
    
    // Output report
    match args.output {
        Some(path) => {
            // Create parent directories if needed
            if let Some(parent) = Path::new(&path).parent() {
                std::fs::create_dir_all(parent).ok();
            }
            let mut file = File::create(&path).expect("Failed to create output file");
            file.write_all(report.as_bytes()).expect("Failed to write report");
            println!("Report written to: {}", path);
        }
        None => {
            println!("{}", report);
        }
    }
    
    // Summary
    let passed = results.iter().filter(|r| r.passed).count();
    let total = results.len();
    println!("\n=== Summary ===");
    println!("Passed: {}/{}", passed, total);
    
    if passed < total {
        std::process::exit(1);
    }
}

fn generate_markdown_report(results: &[BenchmarkResult]) -> String {
    let mut report = String::new();
    
    report.push_str("# FEA Benchmark Report\n\n");
    report.push_str(&format!("**Generated:** {}\n\n", chrono::Local::now().format("%Y-%m-%d %H:%M:%S")));
    
    // Quick summary table
    report.push_str("## Summary\n\n");
    report.push_str("| # | Benchmark | Status | Max Error | Time |\n");
    report.push_str("|---|-----------|--------|-----------|------|\n");
    
    for (i, result) in results.iter().enumerate() {
        let status = if result.passed { "✅ PASS" } else { "❌ FAIL" };
        let max_error = result.metrics.iter()
            .map(|m| m.relative_error.abs())
            .fold(0.0_f64, f64::max);
        let time_str = format_time(result.elapsed_ms);
        let error_str = if max_error < 0.001 {
            format!("{:.2e}", max_error)
        } else {
            format!("{:.1}%", max_error * 100.0)
        };
        report.push_str(&format!(
            "| {} | {} | {} | {} | {} |\n",
            i + 1, result.name, status, error_str, time_str
        ));
    }
    
    // Legend
    report.push_str("\n**Legend:** Max Error shown as percentage (>0.1%) or scientific notation (<0.1%)\n\n");
    
    report.push_str("---\n\n");
    report.push_str("## Detailed Results\n\n");
    
    for (i, result) in results.iter().enumerate() {
        report.push_str(&format!("### {}. {}\n\n", i + 1, result.name));
        report.push_str(&format!("**Description:** {}\n\n", result.description));
        
        let status_emoji = if result.passed { "✅" } else { "❌" };
        report.push_str(&format!("**Status:** {} {}\n\n", status_emoji, 
            if result.passed { "PASS" } else { "FAIL" }));
        
        if !result.metrics.is_empty() {
            // Group metrics by mesh level if present
            let has_mesh_levels = result.metrics.iter()
                .any(|m| m.name.contains("coarse") || m.name.contains("medium") || m.name.contains("fine"));
            
            if has_mesh_levels {
                report.push_str("#### Results by Mesh Refinement\n\n");
                report.push_str(generate_grouped_metrics_table(&result.metrics).as_str());
            } else {
                report.push_str("#### Metrics\n\n");
                report.push_str(generate_metrics_table(&result.metrics).as_str());
            }
            
            // Add convergence plot if multiple mesh levels
            if has_mesh_levels {
                report.push_str(&generate_convergence_section(&result.metrics));
            }
        }
        
        if let Some(notes) = &result.notes {
            report.push_str(&format!("\n**Parameters:**\n{}\n", notes));
        }
        
        report.push_str("\n---\n\n");
    }
    
    // Appendix: Theory
    report.push_str("## Appendix: Analytical Solutions\n\n");
    report.push_str(THEORY_APPENDIX);
    
    report
}

fn generate_metrics_table(metrics: &[rust_fea::benchmarks::MetricComparison]) -> String {
    let mut table = String::new();
    
    table.push_str("| Metric | Analytical | FEA | Error | Status |\n");
    table.push_str("|--------|------------|-----|-------|--------|\n");
    
    for metric in metrics {
        let ana_str = format_value_smart(metric.analytical);
        let comp_str = format_value_smart(metric.computed);
        let error_str = format_error(metric.relative_error, metric.tolerance);
        let status = if metric.passed() { "✓" } else { "✗" };
        
        table.push_str(&format!(
            "| {} | {} | {} | {} | {} |\n",
            clean_metric_name(&metric.name),
            ana_str,
            comp_str,
            error_str,
            status
        ));
    }
    
    table
}

fn generate_grouped_metrics_table(metrics: &[rust_fea::benchmarks::MetricComparison]) -> String {
    let mut table = String::new();
    
    // Extract mesh levels
    let levels: Vec<&str> = vec!["coarse", "medium", "fine", "refine_coarse", "refine_medium", "refine_fine"];
    
    // Find unique base metric names
    let mut base_names: Vec<String> = Vec::new();
    for metric in metrics {
        let base = levels.iter()
            .fold(metric.name.clone(), |s, l| s.replace(&format!("{}_", l), ""));
        if !base_names.contains(&base) {
            base_names.push(base);
        }
    }
    
    table.push_str("| Mesh | Metric | Analytical | FEA | Error | Status |\n");
    table.push_str("|------|--------|------------|-----|-------|--------|\n");
    
    for metric in metrics {
        let mesh_level = levels.iter()
            .find(|l| metric.name.starts_with(*l))
            .map(|s| capitalize(s))
            .unwrap_or_else(|| "-".to_string());
        
        let base_name = levels.iter()
            .fold(metric.name.clone(), |s, l| s.replace(&format!("{}_", l), ""));
        
        let ana_str = format_value_smart(metric.analytical);
        let comp_str = format_value_smart(metric.computed);
        let error_str = format_error(metric.relative_error, metric.tolerance);
        let status = if metric.passed() { "✓" } else { "✗" };
        
        table.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} |\n",
            mesh_level,
            clean_metric_name(&base_name),
            ana_str,
            comp_str,
            error_str,
            status
        ));
    }
    
    table
}

fn generate_convergence_section(metrics: &[rust_fea::benchmarks::MetricComparison]) -> String {
    let mut section = String::new();
    
    // Build convergence data
    let levels = vec![
        ("coarse", 1.0),
        ("medium", 0.5),
        ("fine", 0.25),
        ("refine_coarse", 1.0),
        ("refine_medium", 0.5),
        ("refine_fine", 0.25),
    ];
    
    // Find primary error metric
    let error_metrics: Vec<_> = metrics.iter()
        .filter(|m| m.name.contains("error") || m.name.contains("deflection") || m.name.contains("u_r"))
        .collect();
    
    if error_metrics.len() >= 2 {
        section.push_str("\n#### Convergence Analysis\n\n");
        section.push_str("```\n");
        section.push_str("Error vs Mesh Refinement (schematic):\n\n");
        section.push_str("    Error\n");
        section.push_str("      │\n");
        
        // Simple ASCII convergence plot
        let max_err = error_metrics.iter()
            .map(|m| m.relative_error.abs())
            .fold(0.0_f64, f64::max);
        
        if max_err > 1e-10 {
            for metric in error_metrics.iter().take(6) {
                let bar_len = ((metric.relative_error.abs() / max_err) * 20.0) as usize;
                let bar: String = "█".repeat(bar_len.max(1));
                let level = levels.iter()
                    .find(|(l, _)| metric.name.starts_with(*l))
                    .map(|(l, _)| *l)
                    .unwrap_or("?");
                section.push_str(&format!("  {:>8} │{} {:.1}%\n", 
                    level, bar, metric.relative_error.abs() * 100.0));
            }
        }
        
        section.push_str("           └───────────────────────\n");
        section.push_str("              Mesh refinement →\n");
        section.push_str("```\n");
    }
    
    section
}

fn format_value_smart(value: f64) -> String {
    let abs_val = value.abs();
    
    if abs_val == 0.0 {
        return "0".to_string();
    }
    
    // Use SI prefix for engineering values
    if abs_val >= 1e6 || abs_val < 1e-3 {
        format_with_units(value, "")
    } else if abs_val >= 1000.0 {
        format!("{:.1}", value)
    } else if abs_val >= 1.0 {
        format!("{:.3}", value)
    } else {
        format!("{:.4}", value)
    }
}

fn format_error(error: f64, tolerance: f64) -> String {
    let abs_err = error.abs();
    let pct = abs_err * 100.0;
    
    if abs_err < 0.0001 {
        format!("{:.2e}", error)
    } else if pct < 1.0 {
        format!("{:.2}%", pct)
    } else {
        format!("{:.1}% (tol: {:.0}%)", pct, tolerance * 100.0)
    }
}

fn format_time(ms: f64) -> String {
    if ms < 1.0 {
        format!("{:.0} μs", ms * 1000.0)
    } else if ms < 1000.0 {
        format!("{:.0} ms", ms)
    } else {
        format!("{:.1} s", ms / 1000.0)
    }
}

fn clean_metric_name(name: &str) -> String {
    name.replace("_", " ")
        .replace("  ", " ")
}

fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}

fn generate_json_report(results: &[BenchmarkResult]) -> String {
    serde_json::to_string_pretty(results).unwrap_or_else(|_| "{}".to_string())
}

fn generate_csv_report(results: &[BenchmarkResult]) -> String {
    let mut csv = String::new();
    csv.push_str("benchmark,metric,analytical,computed,relative_error,tolerance,passed,error_pct\n");
    
    for result in results {
        for metric in &result.metrics {
            csv.push_str(&format!(
                "{},{},{:.10e},{:.10e},{:.10e},{:.10e},{},{:.4}\n",
                result.name,
                metric.name,
                metric.analytical,
                metric.computed,
                metric.relative_error,
                metric.tolerance,
                metric.relative_error.abs() <= metric.tolerance,
                metric.relative_error.abs() * 100.0
            ));
        }
    }
    
    csv
}

const THEORY_APPENDIX: &str = r#"
### Benchmark 1-3: Fundamental Element Tests

These test uniform stress/strain states that 8-node brick elements should reproduce exactly:

- **Uniaxial Tension:** σ_xx = F/A, ε_xx = σ/E, u_x = ε·x
- **Pure Shear:** τ_xy = Gγ, G = E/[2(1+ν)]
- **Hydrostatic:** σ_kk = -3p, ε_v = -p/K, K = E/[3(1-2ν)]

### Benchmark 4: Torsion Shaft

For circular shaft: φ = TL/(GJ), τ_max = TR/J, J = πR⁴/2

Square sections use torsion constant k ≈ 0.1406 for J_eff = k·a⁴

### Benchmark 5: Hollow Sphere (Lamé Solution)

σ_r = A - B/r³, σ_θ = A + B/(2r³)

where A, B satisfy pressure BCs at inner/outer radii.

### Benchmark 6: Spherical Cavity

For cavity radius a in infinite medium with remote stress p∞:
- σ_r = -p∞(1 - a³/r³)
- σ_θ = -p∞(1 + a³/2r³)
- Stress concentration: 1.5× at surface

### Benchmark 7: Boussinesq Half-Space

Point load P on half-space: σ_z = -3Pz³/(2πR⁵)

Surface displacement: u_z(r,0) = P(1-ν²)/(πEr)

### Benchmarks 9-10: Hertzian Contact

Contact radius: a = (3FR*/4E*)^(1/3)
Max pressure: p₀ = 3F/(2πa²)  
Approach: δ = a²/R*

where:
- 1/R* = 1/R₁ + 1/R₂ (effective radius)
- 1/E* = (1-ν₁²)/E₁ + (1-ν₂²)/E₂ (effective modulus)
"#;
