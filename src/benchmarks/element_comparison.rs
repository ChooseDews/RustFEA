// Benchmark: C3D8 vs C3D20 Element Comparison
// Compares linear (8-node) vs quadratic (20-node) brick elements
// on a cantilever beam bending problem across multiple mesh densities
//
// Key metrics:
// - Accuracy vs analytical solution
// - DOF count vs error (efficiency)
// - Convergence rates

use crate::benchmarks::{BenchmarkResult, MetricComparison};
use crate::benchmarks::mesh_utils::{generate_block_mesh, generate_block_mesh_c3d20};
use crate::bc::{FixedCondition, LoadCondition};
use crate::simulation::Simulation;
use nalgebra::DVector;
use log::info;
use std::fs::File;
use std::io::Write;

/// Material properties (Aluminum 6061-T6 - matches Material::aluminum() default)
const E: f64 = 68.9e9;    // Young's modulus (Pa)
const NU: f64 = 0.33;     // Poisson's ratio

/// Geometry - beam along X axis, cross-section in YZ
const LENGTH: f64 = 1.0;    // Beam length in x-direction (m)
const WIDTH: f64 = 0.1;     // Width in y-direction (m)  
const HEIGHT: f64 = 0.1;    // Height in z-direction (m)

/// Loading
const TIP_LOAD: f64 = 1000.0;  // Tip load in -Z direction (N)

/// Result from a single mesh configuration
#[derive(Debug, Clone)]
pub struct MeshResult {
    pub element_type: String,
    pub mesh_label: String,
    pub num_elements: usize,
    pub num_nodes: usize,
    pub dof_count: usize,
    pub tip_deflection: f64,
    pub error_percent: f64,
    pub solve_time_ms: f64,
}

/// Run the comparison benchmark and return all results
pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "C3D8 vs C3D20 Element Comparison",
        "Compares linear (8-node) and quadratic (20-node) hexahedral elements \
         on a cantilever beam bending problem. Evaluates accuracy vs DOF count."
    );
    
    // Second moment of area
    let i_zz = WIDTH * HEIGHT.powi(3) / 12.0;
    
    // Analytical tip deflection (Euler-Bernoulli beam theory)
    // For load in -Z: δ = P * L³ / (3 * E * I)
    let delta_analytical = TIP_LOAD * LENGTH.powi(3) / (3.0 * E * i_zz);
    
    info!("C3D8 vs C3D20 Element Comparison");
    info!("  Problem: Cantilever Beam Bending");
    info!("  L = {} m, W = {} m, H = {} m", LENGTH, WIDTH, HEIGHT);
    info!("  E = {:.2e} Pa, ν = {}", E, NU);
    info!("  Tip load = {} N", TIP_LOAD);
    info!("  Analytical tip deflection: {:.6e} m", delta_analytical);
    info!("");
    
    // Mesh refinements: (nx, ny, nz)
    // Prioritize elements in the long (X) direction for beam problems
    let c3d8_refinements = vec![
        (1, 1, 1),    // 1 element
        (2, 1, 1),    // 2 elements
        (4, 1, 1),    // 4 elements
        (8, 1, 1),    // 8 elements
        (16, 1, 1),   // 16 elements
        (32, 1, 1),   // 32 elements
        (16, 2, 2),   // 64 elements
        (32, 2, 2),   // 128 elements
        (64, 2, 2),   // 256 elements
        (64, 4, 4),   // 1024 elements
        (128, 4, 4),  // 2048 elements
        (128, 8, 8),  // 8192 elements - very fine
    ];
    
    let c3d20_refinements = vec![
        (1, 1, 1),    // 1 element
        (2, 1, 1),    // 2 elements
        (4, 1, 1),    // 4 elements
        (8, 1, 1),    // 8 elements
        (16, 1, 1),   // 16 elements
        (16, 2, 2),   // 64 elements
        (32, 2, 2),   // 128 elements
        (64, 2, 2),   // 256 elements
        (32, 4, 4),   // 512 elements
        (64, 4, 4),   // 1024 elements
    ];
    
    let mut c3d8_results: Vec<MeshResult> = Vec::new();
    let mut c3d20_results: Vec<MeshResult> = Vec::new();
    
    // Run C3D8 benchmarks
    info!("--- C3D8 (8-node linear hexahedral) ---");
    for &(nx, ny, nz) in &c3d8_refinements {
        match run_c3d8(nx, ny, nz, delta_analytical) {
            Ok(mr) => {
                info!("  {}x{}x{}: {} elements, {} DOF, error = {:.2}%",
                    nx, ny, nz, mr.num_elements, mr.dof_count, mr.error_percent);
                c3d8_results.push(mr);
            }
            Err(e) => {
                info!("  {}x{}x{}: FAILED - {}", nx, ny, nz, e);
            }
        }
    }
    
    info!("");
    info!("--- C3D20 (20-node quadratic hexahedral) ---");
    for &(nx, ny, nz) in &c3d20_refinements {
        match run_c3d20(nx, ny, nz, delta_analytical) {
            Ok(mr) => {
                info!("  {}x{}x{}: {} elements, {} DOF, error = {:.2}%",
                    nx, ny, nz, mr.num_elements, mr.dof_count, mr.error_percent);
                c3d20_results.push(mr);
            }
            Err(e) => {
                info!("  {}x{}x{}: FAILED - {}", nx, ny, nz, e);
            }
        }
    }
    
    // Record metrics
    result.add_metric(MetricComparison::new(
        "analytical_tip_deflection",
        delta_analytical,
        delta_analytical,
        0.0,
    ));
    
    // Record best results for each element type
    if let Some(best_c3d8) = c3d8_results.last() {
        result.add_metric(MetricComparison::new(
            "c3d8_best_tip_deflection",
            delta_analytical,
            best_c3d8.tip_deflection,
            0.30,  // 30% tolerance
        ));
    }
    
    if let Some(best_c3d20) = c3d20_results.last() {
        result.add_metric(MetricComparison::new(
            "c3d20_best_tip_deflection",
            delta_analytical,
            best_c3d20.tip_deflection,
            0.30,  // 30% tolerance
        ));
    }
    
    // Generate comparison table and plot
    let table = generate_comparison_table(&c3d8_results, &c3d20_results, delta_analytical);
    info!("\n{}", table);
    
    let svg = generate_comparison_svg(&c3d8_results, &c3d20_results, delta_analytical);
    
    // Save SVG to file
    if let Ok(mut file) = File::create("element_comparison.svg") {
        let _ = file.write_all(svg.as_bytes());
        info!("\nPlot saved to: element_comparison.svg");
    }
    
    // Generate timing SVG
    let timing_svg = generate_timing_svg(&c3d8_results, &c3d20_results);
    if let Ok(mut file) = File::create("element_timing.svg") {
        let _ = file.write_all(timing_svg.as_bytes());
        info!("Timing plot saved to: element_timing.svg");
    }
    
    // Print ASCII summary
    info!("\n{}", generate_ascii_summary(&c3d8_results, &c3d20_results));
    
    // Overall pass/fail
    result.passed = !c3d8_results.is_empty() && !c3d20_results.is_empty();
    
    result
}

fn run_c3d8(nx: usize, ny: usize, nz: usize, analytical: f64) -> Result<MeshResult, String> {
    use std::time::Instant;
    
    // Generate C3D8 mesh
    let mut mesh = generate_block_mesh(LENGTH, WIDTH, HEIGHT, nx, ny, nz);
    
    let num_elements = nx * ny * nz;
    let num_nodes = mesh.nodes.len();
    let dof_count = num_nodes * 3;
    
    // Shift to center cross-section at Y=0, Z=0
    for (_, node) in mesh.nodes.iter_mut() {
        node.coordinates[1] -= WIDTH / 2.0;
        node.coordinates[2] -= HEIGHT / 2.0;
    }
    
    // Get node groups
    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");
    
    // Create simulation
    let mut sim = Simulation::from_mesh(mesh, 3);
    
    // Apply BCs
    let fixed_bc = FixedCondition::new(
        x_min_nodes.clone(),
        vec![Some(0.0), Some(0.0), Some(0.0)],
    );
    sim.add_boundary_condition(Box::new(fixed_bc));
    
    let load_bc = LoadCondition::new(
        x_max_nodes.clone(),
        DVector::from_vec(vec![0.0, 0.0, -TIP_LOAD]),
    );
    sim.add_boundary_condition(Box::new(load_bc));
    
    // Solve with timing
    let start = Instant::now();
    sim.solve();
    let solve_time_ms = start.elapsed().as_secs_f64() * 1000.0;
    
    // Get tip deflection
    let nodes = sim.nodes();
    let tip_deflections: Vec<f64> = x_max_nodes.iter()
        .filter_map(|&id| nodes.get(id))
        .map(|node| node.displacement.z)
        .collect();
    
    let avg_tip_deflection = if !tip_deflections.is_empty() {
        tip_deflections.iter().sum::<f64>() / tip_deflections.len() as f64
    } else {
        return Err("No tip nodes found".to_string());
    };
    
    // Deflection is negative (downward), analytical is positive magnitude
    let computed_delta = -avg_tip_deflection;
    let error = ((computed_delta - analytical) / analytical).abs() * 100.0;
    
    Ok(MeshResult {
        element_type: "C3D8".to_string(),
        mesh_label: format!("{}x{}x{}", nx, ny, nz),
        num_elements,
        num_nodes,
        dof_count,
        tip_deflection: computed_delta,
        error_percent: error,
        solve_time_ms,
    })
}

fn run_c3d20(nx: usize, ny: usize, nz: usize, analytical: f64) -> Result<MeshResult, String> {
    use std::time::Instant;
    
    // Generate C3D20 mesh
    let mut mesh = generate_block_mesh_c3d20(LENGTH, WIDTH, HEIGHT, nx, ny, nz);
    
    let num_elements = nx * ny * nz;
    let num_nodes = mesh.nodes.len();
    let dof_count = num_nodes * 3;
    
    // Shift to center cross-section at Y=0, Z=0
    for (_, node) in mesh.nodes.iter_mut() {
        node.coordinates[1] -= WIDTH / 2.0;
        node.coordinates[2] -= HEIGHT / 2.0;
    }
    
    // Get node groups
    let x_min_nodes = mesh.get_nodes_in_group("x_min");
    let x_max_nodes = mesh.get_nodes_in_group("x_max");
    
    // Create simulation
    let mut sim = Simulation::from_mesh(mesh, 3);
    
    // Apply BCs
    let fixed_bc = FixedCondition::new(
        x_min_nodes.clone(),
        vec![Some(0.0), Some(0.0), Some(0.0)],
    );
    sim.add_boundary_condition(Box::new(fixed_bc));
    
    let load_bc = LoadCondition::new(
        x_max_nodes.clone(),
        DVector::from_vec(vec![0.0, 0.0, -TIP_LOAD]),
    );
    sim.add_boundary_condition(Box::new(load_bc));
    
    // Solve with timing
    let start = Instant::now();
    sim.solve();
    let solve_time_ms = start.elapsed().as_secs_f64() * 1000.0;
    
    // Get tip deflection
    let nodes = sim.nodes();
    let tip_deflections: Vec<f64> = x_max_nodes.iter()
        .filter_map(|&id| nodes.get(id))
        .map(|node| node.displacement.z)
        .collect();
    
    let avg_tip_deflection = if !tip_deflections.is_empty() {
        tip_deflections.iter().sum::<f64>() / tip_deflections.len() as f64
    } else {
        return Err("No tip nodes found".to_string());
    };
    
    // Deflection is negative (downward), analytical is positive magnitude
    let computed_delta = -avg_tip_deflection;
    let error = ((computed_delta - analytical) / analytical).abs() * 100.0;
    
    Ok(MeshResult {
        element_type: "C3D20".to_string(),
        mesh_label: format!("{}x{}x{}", nx, ny, nz),
        num_elements,
        num_nodes,
        dof_count,
        tip_deflection: computed_delta,
        error_percent: error,
        solve_time_ms,
    })
}

fn generate_comparison_table(c3d8: &[MeshResult], c3d20: &[MeshResult], analytical: f64) -> String {
    let mut output = String::new();
    output.push_str("### Element Comparison Results\n\n");
    output.push_str(&format!("Analytical tip deflection: {:.6e} m\n\n", analytical));
    output.push_str("| Element | Mesh | Elements | Nodes | DOF | Deflection (m) | Error (%) |\n");
    output.push_str("|---------|------|----------|-------|-----|----------------|----------|\n");
    
    for r in c3d8 {
        output.push_str(&format!(
            "| {} | {} | {} | {} | {} | {:.4e} | {:.2} |\n",
            r.element_type, r.mesh_label, r.num_elements, r.num_nodes,
            r.dof_count, r.tip_deflection, r.error_percent
        ));
    }
    
    output.push_str("|---------|------|----------|-------|-----|----------------|----------|\n");
    
    for r in c3d20 {
        output.push_str(&format!(
            "| {} | {} | {} | {} | {} | {:.4e} | {:.2} |\n",
            r.element_type, r.mesh_label, r.num_elements, r.num_nodes,
            r.dof_count, r.tip_deflection, r.error_percent
        ));
    }
    
    output
}

fn generate_ascii_summary(c3d8: &[MeshResult], c3d20: &[MeshResult]) -> String {
    let mut output = String::new();
    output.push_str("### Summary: DOF vs Error\n\n");
    
    output.push_str("C3D8 (linear):\n");
    for r in c3d8 {
        output.push_str(&format!("  {:>6} DOF -> {:>6.2}% error\n", r.dof_count, r.error_percent));
    }
    
    output.push_str("\nC3D20 (quadratic):\n");
    for r in c3d20 {
        output.push_str(&format!("  {:>6} DOF -> {:>6.2}% error\n", r.dof_count, r.error_percent));
    }
    
    // Efficiency comparison
    if let (Some(c3d8_fine), Some(c3d20_coarse)) = (c3d8.last(), c3d20.first()) {
        if c3d20_coarse.error_percent < c3d8_fine.error_percent {
            let dof_ratio = c3d8_fine.dof_count as f64 / c3d20_coarse.dof_count as f64;
            output.push_str(&format!(
                "\n**Key finding:** C3D20 with {} DOF ({:.2}% error) beats \
                 C3D8 with {} DOF ({:.2}% error) - {:.1}x more efficient!\n",
                c3d20_coarse.dof_count, c3d20_coarse.error_percent,
                c3d8_fine.dof_count, c3d8_fine.error_percent, dof_ratio
            ));
        }
    }
    
    output
}

/// Generate SVG plot comparing C3D8 vs C3D20
fn generate_comparison_svg(c3d8: &[MeshResult], c3d20: &[MeshResult], analytical: f64) -> String {
    let width = 800;
    let height = 500;
    let margin = 70;
    let plot_width = width - 2 * margin;
    let plot_height = height - 2 * margin;
    
    // Collect all data for axis scaling
    let all_results: Vec<&MeshResult> = c3d8.iter().chain(c3d20.iter()).collect();
    if all_results.is_empty() {
        return String::from("<svg></svg>");
    }
    
    // Get data ranges (log scale for both axes)
    let all_dof: Vec<f64> = all_results.iter().map(|r| r.dof_count as f64).collect();
    let all_err: Vec<f64> = all_results.iter().map(|r| r.error_percent.max(0.1)).collect();
    
    let log_dof_min = all_dof.iter().map(|x| x.log10()).fold(f64::INFINITY, f64::min);
    let log_dof_max = all_dof.iter().map(|x| x.log10()).fold(f64::NEG_INFINITY, f64::max);
    let log_err_min = all_err.iter().map(|x| x.log10()).fold(f64::INFINITY, f64::min);
    let log_err_max = all_err.iter().map(|x| x.log10()).fold(f64::NEG_INFINITY, f64::max);
    
    // Add padding
    let dof_range = log_dof_max - log_dof_min;
    let err_range = log_err_max - log_err_min;
    let log_dof_min = log_dof_min - dof_range * 0.15;
    let log_dof_max = log_dof_max + dof_range * 0.15;
    let log_err_min = (log_err_min - err_range * 0.2).max(-1.0);
    let log_err_max = log_err_max + err_range * 0.15;
    
    let scale_x = |dof: f64| -> f64 {
        margin as f64 + (dof.log10() - log_dof_min) / (log_dof_max - log_dof_min) * plot_width as f64
    };
    let scale_y = |err: f64| -> f64 {
        margin as f64 + (log_err_max - err.max(0.1).log10()) / (log_err_max - log_err_min) * plot_height as f64
    };
    
    let mut svg = String::new();
    
    // SVG header with embedded styles
    svg.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif">
  <defs>
    <style>
      .title {{ font-size: 16px; font-weight: bold; }}
      .axis-label {{ font-size: 12px; }}
      .tick-label {{ font-size: 10px; }}
      .legend-text {{ font-size: 11px; }}
      .note {{ font-size: 10px; fill: #666; }}
      .grid {{ stroke: #e0e0e0; stroke-width: 1; }}
      .axis {{ stroke: black; stroke-width: 2; }}
      .c3d8-line {{ fill: none; stroke: #1976D2; stroke-width: 2.5; }}
      .c3d8-point {{ fill: #1976D2; }}
      .c3d20-line {{ fill: none; stroke: #D32F2F; stroke-width: 2.5; }}
      .c3d20-point {{ fill: #D32F2F; }}
    </style>
  </defs>
"#, width, height));
    
    // Background
    svg.push_str(&format!(
        r#"  <rect width="{}" height="{}" fill="white"/>
"#, width, height));
    
    // Title
    svg.push_str(&format!(
        r#"  <text x="{}" y="30" text-anchor="middle" class="title">C3D8 vs C3D20: DOF vs Error (Cantilever Beam)</text>
"#, width / 2));
    
    // Grid lines
    for i in 0..=5 {
        let x = margin + i * plot_width / 5;
        svg.push_str(&format!(
            r#"  <line x1="{}" y1="{}" x2="{}" y2="{}" class="grid"/>
"#, x, margin, x, margin + plot_height));
    }
    for i in 0..=4 {
        let y = margin + i * plot_height / 4;
        svg.push_str(&format!(
            r#"  <line x1="{}" y1="{}" x2="{}" y2="{}" class="grid"/>
"#, margin, y, margin + plot_width, y));
    }
    
    // Axes
    svg.push_str(&format!(
        r#"  <line x1="{}" y1="{}" x2="{}" y2="{}" class="axis"/>
  <line x1="{}" y1="{}" x2="{}" y2="{}" class="axis"/>
"#, margin, margin + plot_height, margin + plot_width, margin + plot_height,
    margin, margin, margin, margin + plot_height));
    
    // X-axis label
    svg.push_str(&format!(
        r#"  <text x="{}" y="{}" text-anchor="middle" class="axis-label">Degrees of Freedom (log scale)</text>
"#, width / 2, height - 15));
    
    // Y-axis label
    svg.push_str(&format!(
        r#"  <text x="18" y="{}" text-anchor="middle" class="axis-label" transform="rotate(-90, 18, {})">Error % (log scale)</text>
"#, height / 2, height / 2));
    
    // X-axis tick labels
    let dof_ticks: [f64; 7] = [100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0, 100000.0];
    for dof in dof_ticks {
        if dof.log10() >= log_dof_min && dof.log10() <= log_dof_max {
            let x = scale_x(dof);
            let label = if dof >= 1000.0 {
                format!("{}k", (dof / 1000.0) as i32)
            } else {
                format!("{}", dof as i32)
            };
            svg.push_str(&format!(
                r#"  <text x="{:.0}" y="{}" text-anchor="middle" class="tick-label">{}</text>
"#, x, margin + plot_height + 15, label));
        }
    }
    
    // Y-axis tick labels - all use scale_y() for consistent log positioning
    let err_ticks: [f64; 6] = [0.2, 1.0, 2.0, 5.0, 10.0, 100.0];
    for err in err_ticks {
        if err.log10() >= log_err_min && err.log10() <= log_err_max {
            let y = scale_y(err);
            let label = if err < 1.0 {
                format!("{:.1}%", err)
            } else {
                format!("{:.0}%", err)
            };
            svg.push_str(&format!(
                r#"  <text x="{}" y="{:.0}" text-anchor="end" class="tick-label">{}</text>
"#, margin - 5, y + 4.0, label));
        }
    }
    // Note: 0% cannot be shown on log scale (log(0) = -∞)
    
    // C3D8 line and points (blue)
    if !c3d8.is_empty() {
        let mut path = String::new();
        for (i, r) in c3d8.iter().enumerate() {
            let x = scale_x(r.dof_count as f64);
            let y = scale_y(r.error_percent);
            if i == 0 {
                path.push_str(&format!("M {:.1} {:.1}", x, y));
            } else {
                path.push_str(&format!(" L {:.1} {:.1}", x, y));
            }
        }
        svg.push_str(&format!(r#"  <path d="{}" class="c3d8-line"/>
"#, path));
        
        for r in c3d8 {
            let x = scale_x(r.dof_count as f64);
            let y = scale_y(r.error_percent);
            if r.num_elements == 1 {
                // Star marker for single element
                svg.push_str(&format!(
                    "  <polygon points=\"{:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1}\" class=\"c3d8-point\" stroke=\"#1976D2\" stroke-width=\"1.5\"/>\n", 
                    x, y - 10.0,           // top
                    x + 2.5, y - 3.0,      // inner right-top
                    x + 9.5, y - 3.0,      // outer right-top
                    x + 4.0, y + 2.0,      // inner right
                    x + 6.0, y + 9.0,      // outer right-bottom
                    x, y + 5.0,            // bottom center
                    x - 6.0, y + 9.0,      // outer left-bottom
                    x - 4.0, y + 2.0,      // inner left
                    x - 9.5, y - 3.0,      // outer left-top
                    x - 2.5, y - 3.0       // inner left-top
                ));
                // Label with error % to the left of star
                svg.push_str(&format!(
                    "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" class=\"tick-label\" fill=\"#1976D2\">{:.0}%</text>\n",
                    x - 12.0, y + 4.0, r.error_percent
                ));
            } else {
                svg.push_str(&format!(
                    r#"  <circle cx="{:.1}" cy="{:.1}" r="6" class="c3d8-point"/>
"#, x, y));
            }
        }
    }
    
    // C3D20 line and points (red)
    if !c3d20.is_empty() {
        let mut path = String::new();
        for (i, r) in c3d20.iter().enumerate() {
            let x = scale_x(r.dof_count as f64);
            let y = scale_y(r.error_percent);
            if i == 0 {
                path.push_str(&format!("M {:.1} {:.1}", x, y));
            } else {
                path.push_str(&format!(" L {:.1} {:.1}", x, y));
            }
        }
        svg.push_str(&format!(r#"  <path d="{}" class="c3d20-line"/>
"#, path));
        
        for r in c3d20 {
            let x = scale_x(r.dof_count as f64);
            let y = scale_y(r.error_percent);
            if r.num_elements == 1 {
                // Star marker for single element
                svg.push_str(&format!(
                    "  <polygon points=\"{:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1}\" class=\"c3d20-point\" stroke=\"#D32F2F\" stroke-width=\"1.5\"/>\n", 
                    x, y - 10.0,           // top
                    x + 2.5, y - 3.0,      // inner right-top
                    x + 9.5, y - 3.0,      // outer right-top
                    x + 4.0, y + 2.0,      // inner right
                    x + 6.0, y + 9.0,      // outer right-bottom
                    x, y + 5.0,            // bottom center
                    x - 6.0, y + 9.0,      // outer left-bottom
                    x - 4.0, y + 2.0,      // inner left
                    x - 9.5, y - 3.0,      // outer left-top
                    x - 2.5, y - 3.0       // inner left-top
                ));
                // Label with error % to the left of star
                svg.push_str(&format!(
                    "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" class=\"tick-label\" fill=\"#D32F2F\">{:.0}%</text>\n",
                    x - 12.0, y + 4.0, r.error_percent
                ));
            } else {
                svg.push_str(&format!(
                    r#"  <rect x="{:.1}" y="{:.1}" width="10" height="10" class="c3d20-point"/>
"#, x - 5.0, y - 5.0));
            }
        }
    }
    
    // Legend
    let legend_x = margin + plot_width - 180;
    let legend_y = margin + 15;
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"170\" height=\"80\" fill=\"white\" stroke=\"#ccc\" rx=\"3\"/>\n",
        legend_x, legend_y));
    svg.push_str(&format!(
        "  <circle cx=\"{}\" cy=\"{}\" r=\"5\" class=\"c3d8-point\"/>\n",
        legend_x + 15, legend_y + 18));
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"legend-text\">C3D8 (8-node linear)</text>\n",
        legend_x + 28, legend_y + 22));
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"8\" height=\"8\" class=\"c3d20-point\"/>\n",
        legend_x + 11, legend_y + 36));
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"legend-text\">C3D20 (20-node quad)</text>\n",
        legend_x + 28, legend_y + 44));
    // Single element star legend entry
    svg.push_str(&format!(
        "  <polygon points=\"{},{} {},{} {},{} {},{} {},{} {},{} {},{} {},{} {},{} {},{}\" fill=\"#666\" stroke=\"#666\" stroke-width=\"1\"/>\n",
        legend_x + 15, legend_y + 54,      // top
        legend_x + 17, legend_y + 58,
        legend_x + 22, legend_y + 58,
        legend_x + 18, legend_y + 61,
        legend_x + 20, legend_y + 66,
        legend_x + 15, legend_y + 63,
        legend_x + 10, legend_y + 66,
        legend_x + 12, legend_y + 61,
        legend_x + 8, legend_y + 58,
        legend_x + 13, legend_y + 58
    ));
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"legend-text\">★ = Single Element</text>\n",
        legend_x + 28, legend_y + 66));
    
    // Note
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"note\">Analytical delta = {:.4e} m (Euler-Bernoulli)</text>\n",
        margin + 5, height - 40, analytical));
    
    svg.push_str("</svg>\n");
    svg
}

/// Generate SVG plot of solve time vs DOF
fn generate_timing_svg(c3d8: &[MeshResult], c3d20: &[MeshResult]) -> String {
    let width = 800;
    let height = 500;
    let margin = 70;
    let plot_width = width - 2 * margin;
    let plot_height = height - 2 * margin;
    
    // Collect all data for axis scaling
    let all_results: Vec<&MeshResult> = c3d8.iter().chain(c3d20.iter()).collect();
    if all_results.is_empty() {
        return String::from("<svg></svg>");
    }
    
    // Get data ranges (log scale for both axes)
    let all_dof: Vec<f64> = all_results.iter().map(|r| r.dof_count as f64).collect();
    let all_time: Vec<f64> = all_results.iter().map(|r| r.solve_time_ms.max(0.01)).collect();
    
    let log_dof_min = all_dof.iter().map(|x| x.log10()).fold(f64::INFINITY, f64::min);
    let log_dof_max = all_dof.iter().map(|x| x.log10()).fold(f64::NEG_INFINITY, f64::max);
    let log_time_min = all_time.iter().map(|x| x.log10()).fold(f64::INFINITY, f64::min);
    let log_time_max = all_time.iter().map(|x| x.log10()).fold(f64::NEG_INFINITY, f64::max);
    
    // Add padding
    let dof_range = log_dof_max - log_dof_min;
    let time_range = log_time_max - log_time_min;
    let log_dof_min = log_dof_min - dof_range * 0.15;
    let log_dof_max = log_dof_max + dof_range * 0.15;
    let log_time_min = log_time_min - time_range * 0.2;
    let log_time_max = log_time_max + time_range * 0.15;
    
    let scale_x = |dof: f64| -> f64 {
        margin as f64 + (dof.log10() - log_dof_min) / (log_dof_max - log_dof_min) * plot_width as f64
    };
    let scale_y = |time: f64| -> f64 {
        margin as f64 + (log_time_max - time.max(0.01).log10()) / (log_time_max - log_time_min) * plot_height as f64
    };
    
    let mut svg = String::new();
    
    // SVG header with embedded styles
    svg.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif">
  <defs>
    <style>
      .title {{ font-size: 16px; font-weight: bold; }}
      .axis-label {{ font-size: 12px; }}
      .tick-label {{ font-size: 10px; }}
      .legend-text {{ font-size: 11px; }}
      .note {{ font-size: 10px; fill: #666; }}
      .grid {{ stroke: #e0e0e0; stroke-width: 1; }}
      .axis {{ stroke: black; stroke-width: 2; }}
      .c3d8-line {{ fill: none; stroke: #1976D2; stroke-width: 2.5; }}
      .c3d8-point {{ fill: #1976D2; }}
      .c3d20-line {{ fill: none; stroke: #D32F2F; stroke-width: 2.5; }}
      .c3d20-point {{ fill: #D32F2F; }}
    </style>
  </defs>
"#, width, height));
    
    // Background
    svg.push_str(&format!(
        r#"  <rect width="{}" height="{}" fill="white"/>
"#, width, height));
    
    // Title
    svg.push_str(&format!(
        r#"  <text x="{}" y="30" text-anchor="middle" class="title">C3D8 vs C3D20: DOF vs Solve Time</text>
"#, width / 2));
    
    // Grid lines
    for i in 0..=5 {
        let x = margin + i * plot_width / 5;
        svg.push_str(&format!(
            r#"  <line x1="{}" y1="{}" x2="{}" y2="{}" class="grid"/>
"#, x, margin, x, margin + plot_height));
    }
    for i in 0..=4 {
        let y = margin + i * plot_height / 4;
        svg.push_str(&format!(
            r#"  <line x1="{}" y1="{}" x2="{}" y2="{}" class="grid"/>
"#, margin, y, margin + plot_width, y));
    }
    
    // Axes
    svg.push_str(&format!(
        r#"  <line x1="{}" y1="{}" x2="{}" y2="{}" class="axis"/>
  <line x1="{}" y1="{}" x2="{}" y2="{}" class="axis"/>
"#, margin, margin + plot_height, margin + plot_width, margin + plot_height,
    margin, margin, margin, margin + plot_height));
    
    // X-axis label
    svg.push_str(&format!(
        r#"  <text x="{}" y="{}" text-anchor="middle" class="axis-label">Degrees of Freedom (log scale)</text>
"#, width / 2, height - 15));
    
    // Y-axis label
    svg.push_str(&format!(
        r#"  <text x="18" y="{}" text-anchor="middle" class="axis-label" transform="rotate(-90, 18, {})">Solve Time (ms, log scale)</text>
"#, height / 2, height / 2));
    
    // X-axis tick labels
    let dof_ticks: [f64; 7] = [100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0, 100000.0];
    for dof in dof_ticks {
        if dof.log10() >= log_dof_min && dof.log10() <= log_dof_max {
            let x = scale_x(dof);
            let label = if dof >= 1000.0 {
                format!("{}k", (dof / 1000.0) as i32)
            } else {
                format!("{}", dof as i32)
            };
            svg.push_str(&format!(
                r#"  <text x="{:.0}" y="{}" text-anchor="middle" class="tick-label">{}</text>
"#, x, margin + plot_height + 15, label));
        }
    }
    
    // Y-axis tick labels (time in ms)
    let time_ticks: [f64; 9] = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0];
    for time in time_ticks {
        if time.log10() >= log_time_min && time.log10() <= log_time_max {
            let y = scale_y(time);
            let label = if time >= 1000.0 {
                format!("{:.0}s", time / 1000.0)
            } else if time >= 1.0 {
                format!("{:.0}ms", time)
            } else {
                format!("{:.1}ms", time)
            };
            svg.push_str(&format!(
                r#"  <text x="{}" y="{:.0}" text-anchor="end" class="tick-label">{}</text>
"#, margin - 5, y + 4.0, label));
        }
    }
    
    // C3D8 line and points (blue)
    if !c3d8.is_empty() {
        let mut path = String::new();
        for (i, r) in c3d8.iter().enumerate() {
            let x = scale_x(r.dof_count as f64);
            let y = scale_y(r.solve_time_ms);
            if i == 0 {
                path.push_str(&format!("M {:.1} {:.1}", x, y));
            } else {
                path.push_str(&format!(" L {:.1} {:.1}", x, y));
            }
        }
        svg.push_str(&format!(r#"  <path d="{}" class="c3d8-line"/>
"#, path));
        
        for r in c3d8 {
            let x = scale_x(r.dof_count as f64);
            let y = scale_y(r.solve_time_ms);
            svg.push_str(&format!(
                r#"  <circle cx="{:.1}" cy="{:.1}" r="6" class="c3d8-point"/>
"#, x, y));
        }
    }
    
    // C3D20 line and points (red)
    if !c3d20.is_empty() {
        let mut path = String::new();
        for (i, r) in c3d20.iter().enumerate() {
            let x = scale_x(r.dof_count as f64);
            let y = scale_y(r.solve_time_ms);
            if i == 0 {
                path.push_str(&format!("M {:.1} {:.1}", x, y));
            } else {
                path.push_str(&format!(" L {:.1} {:.1}", x, y));
            }
        }
        svg.push_str(&format!(r#"  <path d="{}" class="c3d20-line"/>
"#, path));
        
        for r in c3d20 {
            let x = scale_x(r.dof_count as f64);
            let y = scale_y(r.solve_time_ms);
            svg.push_str(&format!(
                r#"  <rect x="{:.1}" y="{:.1}" width="10" height="10" class="c3d20-point"/>
"#, x - 5.0, y - 5.0));
        }
    }
    
    // Legend
    let legend_x = margin + 15;
    let legend_y = margin + 15;
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"150\" height=\"55\" fill=\"white\" stroke=\"#ccc\" rx=\"3\"/>\n",
        legend_x, legend_y));
    svg.push_str(&format!(
        "  <circle cx=\"{}\" cy=\"{}\" r=\"5\" class=\"c3d8-point\"/>\n",
        legend_x + 15, legend_y + 18));
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"legend-text\">C3D8 (8-node linear)</text>\n",
        legend_x + 28, legend_y + 22));
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"8\" height=\"8\" class=\"c3d20-point\"/>\n",
        legend_x + 11, legend_y + 36));
    svg.push_str(&format!(
        "  <text x=\"{}\" y=\"{}\" class=\"legend-text\">C3D20 (20-node quad)</text>\n",
        legend_x + 28, legend_y + 44));
    
    svg.push_str("</svg>\n");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_element_comparison() {
        let _ = env_logger::builder()
            .filter_level(log::LevelFilter::Info)
            .is_test(true)
            .try_init();
        
        let result = run();
        
        println!("\n=== Element Comparison Results ===");
        for metric in &result.metrics {
            println!("  {}: analytical={:.4e}, computed={:.4e}, error={:.2}%",
                metric.name, metric.analytical, metric.computed,
                metric.relative_error.abs() * 100.0);
        }
        
        assert!(result.passed, "Benchmark should complete successfully");
    }
    
    #[test]
    fn test_c3d8_single() {
        let i_zz = WIDTH * HEIGHT.powi(3) / 12.0;
        let analytical = TIP_LOAD * LENGTH.powi(3) / (3.0 * E * i_zz);
        
        let result = run_c3d8(4, 2, 2, analytical);
        assert!(result.is_ok(), "C3D8 should solve: {:?}", result.err());
        
        let mr = result.unwrap();
        assert_eq!(mr.num_elements, 16, "Should have 16 elements");
        assert!(mr.error_percent < 100.0, "Error should be reasonable");
    }
    
    #[test]
    fn test_c3d20_single() {
        let i_zz = WIDTH * HEIGHT.powi(3) / 12.0;
        let analytical = TIP_LOAD * LENGTH.powi(3) / (3.0 * E * i_zz);
        
        let result = run_c3d20(4, 2, 2, analytical);
        assert!(result.is_ok(), "C3D20 should solve: {:?}", result.err());
        
        let mr = result.unwrap();
        assert_eq!(mr.num_elements, 16, "Should have 16 elements");
        assert!(mr.error_percent < 100.0, "Error should be reasonable");
    }
    
    #[test]
    fn test_c3d20_more_accurate() {
        let _ = env_logger::builder()
            .filter_level(log::LevelFilter::Info)
            .is_test(true)
            .try_init();
        
        let i_zz = WIDTH * HEIGHT.powi(3) / 12.0;
        let analytical = TIP_LOAD * LENGTH.powi(3) / (3.0 * E * i_zz);
        
        let c3d8 = run_c3d8(10, 2, 2, analytical).expect("C3D8 should work");
        let c3d20 = run_c3d20(10, 2, 2, analytical).expect("C3D20 should work");
        
        println!("\nSame mesh comparison (10x2x2):");
        println!("  C3D8:  {} DOF, {:.2}% error", c3d8.dof_count, c3d8.error_percent);
        println!("  C3D20: {} DOF, {:.2}% error", c3d20.dof_count, c3d20.error_percent);
        
        // C3D20 should be significantly more accurate for bending
        assert!(c3d20.error_percent < c3d8.error_percent,
            "C3D20 ({:.2}%) should be more accurate than C3D8 ({:.2}%) for bending",
            c3d20.error_percent, c3d8.error_percent);
    }
}
