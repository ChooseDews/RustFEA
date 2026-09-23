// Plotting utilities for benchmark reports
// Generates ASCII/Unicode plots and convergence analysis for benchmark results

use std::collections::HashMap;

/// A data point with mesh refinement level and error
#[derive(Debug, Clone)]
pub struct ConvergencePoint {
    pub mesh_size: f64,   // Characteristic element size (m)
    pub dof_count: usize, // Total degrees of freedom
    pub error: f64,       // Relative error
    pub label: String,    // e.g., "coarse", "medium", "fine"
}

/// Convergence study results for a metric
#[derive(Debug, Clone)]
pub struct ConvergenceStudy {
    pub metric_name: String,
    pub points: Vec<ConvergencePoint>,
    pub convergence_rate: Option<f64>, // Estimated order of convergence
    pub unit: String,
}

impl ConvergenceStudy {
    pub fn new(metric_name: &str, unit: &str) -> Self {
        ConvergenceStudy {
            metric_name: metric_name.to_string(),
            points: Vec::new(),
            convergence_rate: None,
            unit: unit.to_string(),
        }
    }

    pub fn add_point(&mut self, mesh_size: f64, dof: usize, error: f64, label: &str) {
        self.points.push(ConvergencePoint {
            mesh_size,
            dof_count: dof,
            error,
            label: label.to_string(),
        });
    }

    /// Calculate convergence rate using least-squares fit
    /// log(error) ≈ log(C) + p * log(h)
    /// Returns convergence rate p
    pub fn calculate_convergence_rate(&mut self) -> Option<f64> {
        if self.points.len() < 2 {
            return None;
        }

        // Filter out zero or negative errors
        let valid_points: Vec<_> = self
            .points
            .iter()
            .filter(|p| p.error > 1e-15 && p.mesh_size > 1e-15)
            .collect();

        if valid_points.len() < 2 {
            return None;
        }

        // Simple linear regression on log-log data
        let n = valid_points.len() as f64;
        let sum_log_h: f64 = valid_points.iter().map(|p| p.mesh_size.ln()).sum();
        let sum_log_e: f64 = valid_points.iter().map(|p| p.error.ln()).sum();
        let sum_log_h_sq: f64 = valid_points.iter().map(|p| p.mesh_size.ln().powi(2)).sum();
        let sum_log_h_log_e: f64 = valid_points
            .iter()
            .map(|p| p.mesh_size.ln() * p.error.ln())
            .sum();

        let denom = n * sum_log_h_sq - sum_log_h.powi(2);
        if denom.abs() < 1e-15 {
            return None;
        }

        let rate = (n * sum_log_h_log_e - sum_log_h * sum_log_e) / denom;
        self.convergence_rate = Some(rate);
        Some(rate)
    }
}

/// Format a value with appropriate SI prefix and units
pub fn format_with_units(value: f64, unit: &str) -> String {
    let (scaled, prefix) = si_prefix(value);
    if unit.is_empty() {
        format!("{:.3} {}", scaled, prefix).trim().to_string()
    } else {
        format!("{:.3} {}{}", scaled, prefix, unit)
    }
}

/// Get SI prefix for a value
fn si_prefix(value: f64) -> (f64, &'static str) {
    let abs_val = value.abs();
    if abs_val == 0.0 {
        return (0.0, "");
    }

    let prefixes: [(f64, &str); 17] = [
        (1e24, "Y"),
        (1e21, "Z"),
        (1e18, "E"),
        (1e15, "P"),
        (1e12, "T"),
        (1e9, "G"),
        (1e6, "M"),
        (1e3, "k"),
        (1e0, ""),
        (1e-3, "m"),
        (1e-6, "μ"),
        (1e-9, "n"),
        (1e-12, "p"),
        (1e-15, "f"),
        (1e-18, "a"),
        (1e-21, "z"),
        (1e-24, "y"),
    ];

    for &(threshold, prefix) in &prefixes {
        if abs_val >= threshold * 0.9999 {
            return (value / threshold, prefix);
        }
    }

    (value, "")
}

/// Generate an ASCII convergence plot
pub fn generate_ascii_convergence_plot(
    study: &ConvergenceStudy,
    width: usize,
    height: usize,
) -> String {
    if study.points.is_empty() {
        return "No data points for plot".to_string();
    }

    let mut output = String::new();

    // Get data range (use log scale for both axes)
    let log_h: Vec<f64> = study
        .points
        .iter()
        .filter(|p| p.mesh_size > 0.0)
        .map(|p| p.mesh_size.log10())
        .collect();
    let log_e: Vec<f64> = study
        .points
        .iter()
        .filter(|p| p.error > 0.0)
        .map(|p| p.error.log10())
        .collect();

    if log_h.is_empty() || log_e.is_empty() {
        return "Insufficient valid data for plot".to_string();
    }

    let h_min = log_h.iter().cloned().fold(f64::INFINITY, f64::min);
    let h_max = log_h.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let e_min = log_e.iter().cloned().fold(f64::INFINITY, f64::min);
    let e_max = log_e.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Add padding
    let h_range = (h_max - h_min).max(1.0);
    let e_range = (e_max - e_min).max(1.0);
    let h_pad = h_range * 0.1;
    let e_pad = e_range * 0.1;

    let h_min = h_min - h_pad;
    let h_max = h_max + h_pad;
    let e_min = e_min - e_pad;
    let e_max = e_max + e_pad;

    // Create plot grid
    let mut grid: Vec<Vec<char>> = vec![vec![' '; width]; height];

    // Draw axes
    for x in 0..width {
        grid[height - 1][x] = '─';
    }
    for y in 0..height {
        grid[y][0] = '│';
    }
    grid[height - 1][0] = '└';

    // Plot data points
    for (lh, le) in log_h.iter().zip(log_e.iter()) {
        let x = ((lh - h_min) / (h_max - h_min) * (width - 2) as f64) as usize + 1;
        let y = height - 2 - ((le - e_min) / (e_max - e_min) * (height - 2) as f64) as usize;
        if x < width && y < height {
            grid[y][x] = '●';
        }
    }

    // Draw convergence rate line if available
    if let Some(rate) = study.convergence_rate {
        // Draw from first to last point
        let x1 = 2;
        let x2 = width - 2;
        for x in x1..=x2 {
            let lh = h_min + (x - 1) as f64 / (width - 2) as f64 * (h_max - h_min);
            // y = rate * x + intercept
            // Use first point to find intercept
            let first_lh = log_h.first().unwrap_or(&0.0);
            let first_le = log_e.first().unwrap_or(&0.0);
            let intercept = first_le - rate * first_lh;
            let le = rate * lh + intercept;

            let y = ((e_max - le) / (e_max - e_min) * (height - 2) as f64) as usize;
            if y < height - 1 && grid[y][x] == ' ' {
                grid[y][x] = '·';
            }
        }
    }

    // Build output
    output.push_str(&format!(
        "```\n{}: Convergence Plot (log-log)\n",
        study.metric_name
    ));
    if let Some(rate) = study.convergence_rate {
        output.push_str(&format!("Convergence rate: p = {:.2}\n", rate));
    }
    output.push_str(&format!("log₁₀(error)\n"));

    for row in &grid {
        output.push_str(&format!(
            "{:.2e} ",
            10_f64.powf(
                e_max
                    - (grid.iter().position(|r| std::ptr::eq(r, row)).unwrap_or(0) as f64
                        / height as f64)
                        * (e_max - e_min)
            )
        ));
        for &c in row {
            output.push(c);
        }
        output.push('\n');
    }

    output.push_str(&format!(
        "          {:>width$}\n",
        "log₁₀(h)",
        width = width - 10
    ));
    output.push_str("```\n");

    output
}

/// Generate a simple text-based bar chart
pub fn generate_bar_chart(data: &[(String, f64)], label: &str, max_width: usize) -> String {
    if data.is_empty() {
        return String::new();
    }

    let max_val = data.iter().map(|(_, v)| v.abs()).fold(0.0_f64, f64::max);
    if max_val < 1e-15 {
        return String::new();
    }

    let mut output = String::new();
    output.push_str(&format!("{} (relative scale)\n", label));
    output.push_str("```\n");

    let label_width = data.iter().map(|(l, _)| l.len()).max().unwrap_or(10);

    for (name, value) in data {
        let bar_len = (value.abs() / max_val * max_width as f64) as usize;
        let bar: String = "█".repeat(bar_len);
        let sign = if *value < 0.0 { "-" } else { " " };
        output.push_str(&format!(
            "{:>width$} │{}{} {:.2e}\n",
            name,
            sign,
            bar,
            value.abs(),
            width = label_width
        ));
    }
    output.push_str("```\n");

    output
}

/// Generate convergence table in markdown
pub fn generate_convergence_table(study: &ConvergenceStudy) -> String {
    let mut output = String::new();

    output.push_str(&format!("#### {} Convergence\n\n", study.metric_name));
    output.push_str("| Mesh | Element Size | DOF | Error | Error (%) |\n");
    output.push_str("|------|--------------|-----|-------|----------|\n");

    for point in &study.points {
        let size_str = format_with_units(point.mesh_size, "m");
        let error_pct = point.error * 100.0;
        output.push_str(&format!(
            "| {} | {} | {} | {:.2e} | {:.2}% |\n",
            point.label, size_str, point.dof_count, point.error, error_pct
        ));
    }

    if let Some(rate) = study.convergence_rate {
        output.push_str(&format!("\n**Convergence rate:** p ≈ {:.2}", rate));
        if rate >= 1.8 {
            output.push_str(" (quadratic ✓)");
        } else if rate >= 0.8 {
            output.push_str(" (linear)");
        } else {
            output.push_str(" (sublinear ⚠)");
        }
        output.push_str("\n");
    }

    output
}

/// Generate a summary of units for a benchmark
pub fn format_material_properties(e: f64, nu: f64) -> String {
    let e_str = format_with_units(e, "Pa");
    format!(
        "E = {} (Young's modulus), ν = {:.2} (Poisson's ratio)",
        e_str, nu
    )
}

pub fn format_geometry_block(l_x: f64, l_y: f64, l_z: f64) -> String {
    format!(
        "Geometry: {} × {} × {} (L × W × H)",
        format_with_units(l_x, "m"),
        format_with_units(l_y, "m"),
        format_with_units(l_z, "m")
    )
}

pub fn format_pressure(p: f64) -> String {
    format_with_units(p, "Pa")
}

pub fn format_stress(sigma: f64) -> String {
    format_with_units(sigma, "Pa")
}

pub fn format_displacement(u: f64) -> String {
    format_with_units(u, "m")
}

pub fn format_force(f: f64) -> String {
    format_with_units(f, "N")
}

pub fn format_length(l: f64) -> String {
    format_with_units(l, "m")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_si_prefix() {
        assert_eq!(format_with_units(1e9, "Pa"), "1.000 GPa");
        assert_eq!(format_with_units(1e6, "Pa"), "1.000 MPa");
        assert_eq!(format_with_units(1e-3, "m"), "1.000 mm");
        assert_eq!(format_with_units(1e-6, "m"), "1.000 μm");
    }

    #[test]
    fn test_convergence_rate() {
        let mut study = ConvergenceStudy::new("test", "");
        // Perfect quadratic convergence: error = h^2
        study.add_point(0.1, 100, 0.01, "coarse");
        study.add_point(0.05, 400, 0.0025, "medium");
        study.add_point(0.025, 1600, 0.000625, "fine");

        let rate = study.calculate_convergence_rate();
        assert!(rate.is_some());
        assert!((rate.unwrap() - 2.0).abs() < 0.1);
    }
}
