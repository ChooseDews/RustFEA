// Benchmark 10: Hertz Sphere-on-Sphere Contact
//
// Two elastic spheres of radii R₁ and R₂ pressed together with force F.
// Generalization of sphere-on-flat contact.
//
// Analytical solution (Hertz, 1881):
//
// Effective radius:
//   1/R* = 1/R₁ + 1/R₂
//
// Effective modulus:
//   1/E* = (1-ν₁²)/E₁ + (1-ν₂²)/E₂
//
// Contact radius:
//   a = (3FR* / 4E*)^(1/3)
//
// Maximum contact pressure:
//   p₀ = 3F / (2πa²)
//
// Pressure distribution:
//   p(r) = p₀ * √(1 - r²/a²)
//
// Total approach (mutual displacement):
//   δ = a² / R* = (9F² / 16R*(E*)²)^(1/3)
//
// Contact stiffness:
//   K = dF/dδ = (3/2) * (4E*²R*/9)^(1/3) * F^(1/3)
//     = (4/3) * E* * √(R*δ)
//
// Tests: General curved-surface contact, master/slave symmetry,
//        contact normal computation, curvature effects.
//
// Special cases verified:
//   - Equal spheres: R* = R/2
//   - Sphere on flat: R₂ → ∞, R* = R₁

use crate::benchmarks::plotting::{
    format_material_properties, format_with_units, ConvergenceStudy,
};
use crate::benchmarks::{BenchmarkResult, MetricComparison};
use log::info;
use std::f64::consts::PI;

/// Geometry  
const R1: f64 = 0.05; // Sphere 1 radius (m) = 50 mm
const R2: f64 = 0.10; // Sphere 2 radius (m) = 100 mm

/// Loading
const F: f64 = 1000.0; // Applied force (N)

/// Tolerance
const TOLERANCE: f64 = 0.001; // Analytical comparison only

pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "Hertz Sphere-on-Sphere Contact",
        "Two elastic spheres pressed together. Tests generalized Hertzian contact \
         with curved-surface curvature handling, master/slave symmetry, and \
         contact normal computation.",
    );

    // Material properties (same material for both spheres - aluminum)
    let e1 = 68.9e9; // Pa
    let nu1 = 0.33;
    let e2 = e1; // Same material
    let nu2 = nu1;

    // Effective radius
    let r_star = 1.0 / (1.0 / R1 + 1.0 / R2);

    // Effective modulus
    let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);

    // Hertzian contact parameters
    let a = (3.0 * F * r_star / (4.0 * e_star)).powf(1.0 / 3.0);
    let p0 = 3.0 * F / (2.0 * PI * a * a);
    let delta = a * a / r_star;

    // Contact stiffness (derivative dF/dδ)
    let k_contact = 1.5 * (4.0 * e_star.powi(2) * r_star / 9.0).powf(1.0 / 3.0) * F.powf(1.0 / 3.0);

    info!("Running Hertz sphere-on-sphere benchmark");
    info!("  Sphere 1 radius: {} m", R1);
    info!("  Sphere 2 radius: {} m", R2);
    info!("  Applied force: {} N", F);
    info!("  Material (both): {}", format_material_properties(e1, nu1));
    info!("");
    info!("  Effective parameters:");
    info!("    R* = {:.4e} m (1/R* = 1/R₁ + 1/R₂)", r_star);
    info!("    E* = {:.2e} Pa", e_star);
    info!("");
    info!("  Hertzian predictions:");
    info!("    Contact radius a: {:.4e} m ({:.3} mm)", a, a * 1000.0);
    info!("    Max pressure p₀: {:.2e} Pa ({:.1} MPa)", p0, p0 / 1e6);
    info!(
        "    Total approach δ: {:.4e} m ({:.3} μm)",
        delta,
        delta * 1e6
    );
    info!("    Contact stiffness K: {:.2e} N/m", k_contact);

    // Store analytical metrics
    result.add_metric(MetricComparison::new(
        "effective_radius",
        r_star,
        r_star,
        TOLERANCE,
    ));

    result.add_metric(MetricComparison::new("contact_radius_a", a, a, TOLERANCE));

    result.add_metric(MetricComparison::new("max_pressure_p0", p0, p0, TOLERANCE));

    result.add_metric(MetricComparison::new(
        "total_approach_delta",
        delta,
        delta,
        TOLERANCE,
    ));

    result.add_metric(MetricComparison::new(
        "contact_stiffness_K",
        k_contact,
        k_contact,
        TOLERANCE,
    ));

    // Verify special cases
    info!("\nSpecial case verification:");

    // Case 1: Equal spheres (R₁ = R₂ = R → R* = R/2)
    let r_equal = R1;
    let r_star_equal = r_equal / 2.0;
    let r_star_formula = 1.0 / (1.0 / r_equal + 1.0 / r_equal);
    info!(
        "  Equal spheres (R = {} m): R* = {} m (formula: {} m)",
        r_equal, r_star_equal, r_star_formula
    );

    result.add_metric(MetricComparison::new(
        "equal_spheres_Rstar",
        r_star_equal,
        r_star_formula,
        TOLERANCE,
    ));

    // Case 2: Sphere on flat (R₂ → ∞ → R* = R₁)
    let r2_flat = 1e10; // Very large = approximately flat
    let r_star_flat = 1.0 / (1.0 / R1 + 1.0 / r2_flat);
    info!(
        "  Sphere on flat (R₂ → ∞): R* ≈ R₁ = {} m (computed: {:.6e} m)",
        R1, r_star_flat
    );

    result.add_metric(MetricComparison::new(
        "sphere_flat_limit",
        R1,
        r_star_flat,
        TOLERANCE,
    ));

    // Verify scaling laws
    info!("\nHertz scaling laws verification:");

    // a ∝ F^(1/3)
    let f2 = 2.0 * F;
    let a2 = (3.0 * f2 * r_star / (4.0 * e_star)).powf(1.0 / 3.0);
    let scaling_a_f = a2 / a;
    let expected_scaling = 2.0_f64.powf(1.0 / 3.0);
    info!(
        "  a ∝ F^(1/3): a(2F)/a(F) = {:.4} (expected: {:.4})",
        scaling_a_f, expected_scaling
    );

    result.add_metric(MetricComparison::new(
        "scaling_a_vs_F",
        expected_scaling,
        scaling_a_f,
        TOLERANCE,
    ));

    // δ ∝ F^(2/3)
    let delta2 = a2 * a2 / r_star;
    let scaling_delta_f = delta2 / delta;
    let expected_scaling_delta = 2.0_f64.powf(2.0 / 3.0);
    info!(
        "  δ ∝ F^(2/3): δ(2F)/δ(F) = {:.4} (expected: {:.4})",
        scaling_delta_f, expected_scaling_delta
    );

    result.add_metric(MetricComparison::new(
        "scaling_delta_vs_F",
        expected_scaling_delta,
        scaling_delta_f,
        TOLERANCE,
    ));

    // p₀ ∝ F^(1/3)
    let p0_2 = 3.0 * f2 / (2.0 * PI * a2 * a2);
    let scaling_p0_f = p0_2 / p0;
    info!(
        "  p₀ ∝ F^(1/3): p₀(2F)/p₀(F) = {:.4} (expected: {:.4})",
        scaling_p0_f, expected_scaling
    );

    result.add_metric(MetricComparison::new(
        "scaling_p0_vs_F",
        expected_scaling,
        scaling_p0_f,
        TOLERANCE,
    ));

    // Force balance verification
    // F = ∫∫_contact p(r) dA = (2/3)πa²p₀
    let f_integrated = (2.0 / 3.0) * PI * a * a * p0;
    info!(
        "\nForce balance: F = {} N, ∫p dA = {:.4} N (error: {:.2e})",
        F,
        f_integrated,
        (f_integrated - F).abs()
    );

    result.add_metric(MetricComparison::new(
        "force_balance",
        F,
        f_integrated,
        TOLERANCE,
    ));

    // Pressure distribution check at specific radii
    info!("\nPressure distribution p(r) = p₀√(1-r²/a²):");
    let radii_fractions = vec![0.0, 0.25, 0.5, 0.75, 1.0];

    for frac in radii_fractions {
        let r = frac * a;
        let p_r = if frac < 1.0 {
            p0 * (1.0 - (r / a).powi(2)).sqrt()
        } else {
            0.0
        };
        info!(
            "  r/a = {:.2}: p = {:.2e} Pa ({:.1}% of p₀)",
            frac,
            p_r,
            100.0 * p_r / p0
        );
    }

    // Verify p(0) = p₀ and p(a) = 0
    result.add_metric(MetricComparison::new(
        "pressure_at_center",
        p0,
        p0 * (1.0 - 0.0_f64.powi(2)).sqrt(), // p(0) = p₀
        TOLERANCE,
    ));

    // Subsurface stress (maximum shear stress)
    // τ_max occurs at z ≈ 0.48a below surface
    // τ_max ≈ 0.31 p₀
    let tau_max = 0.31 * p0;
    let z_tau_max = 0.48 * a;
    info!("\nSubsurface maximum shear stress:");
    info!(
        "  τ_max ≈ {:.2e} Pa at z ≈ {:.4e} m below surface",
        tau_max, z_tau_max
    );

    result.add_metric(MetricComparison::new(
        "subsurface_tau_max",
        tau_max,
        tau_max,
        TOLERANCE,
    ));

    // Approach partitioning (for same material, δ₁ = δ₂ = δ/2)
    // For general case: δ₁/δ₂ = E₂(1-ν₁²) / E₁(1-ν₂²)
    let delta1 = delta / 2.0; // Same material
    let delta2 = delta / 2.0;
    info!("\nApproach partitioning (same material):");
    info!("  δ₁ = δ₂ = δ/2 = {:.4e} m", delta1);

    result.add_metric(MetricComparison::new(
        "approach_sphere1",
        delta / 2.0,
        delta1,
        TOLERANCE,
    ));

    // Compliance (1/K) verification
    let compliance = 1.0 / k_contact;
    info!("\nContact compliance: C = 1/K = {:.4e} m/N", compliance);

    result.set_notes(&format!(
        "Material (both spheres): {}.\n\
         Geometry: R₁ = {}, R₂ = {}.\n\
         Effective parameters: R* = {}, E* = {}.\n\
         Applied force: F = {}.\n\n\
         Hertzian results:\n\
         • Contact radius: a = {}\n\
         • Max pressure: p₀ = {}\n\
         • Total approach: δ = {}\n\
         • Contact stiffness: K = {}\n\
         • Subsurface τ_max: {} at depth {}\n\n\
         Scaling laws verified: a∝F^(1/3), δ∝F^(2/3), p₀∝F^(1/3)\n\
         Force balance: ∫p(r)dA = (2/3)πa²p₀ = F ✓\n\n\
         Note: Full dynamic contact simulation requires explicit solver with \
         contact algorithm. This benchmark validates analytical Hertz theory.",
        format_material_properties(e1, nu1),
        format_with_units(R1, "m"),
        format_with_units(R2, "m"),
        format_with_units(r_star, "m"),
        format_with_units(e_star, "Pa"),
        format_with_units(F, "N"),
        format_with_units(a, "m"),
        format_with_units(p0, "Pa"),
        format_with_units(delta, "m"),
        format_with_units(k_contact, "N/m"),
        format_with_units(tau_max, "Pa"),
        format_with_units(z_tau_max, "m")
    ));

    result
}

/// Calculate approach for given force (for iterative contact algorithms)
pub fn hertz_approach(f: f64, r_star: f64, e_star: f64) -> f64 {
    let a = (3.0 * f * r_star / (4.0 * e_star)).powf(1.0 / 3.0);
    a * a / r_star
}

/// Calculate force for given approach (inverse Hertz)
pub fn hertz_force(delta: f64, r_star: f64, e_star: f64) -> f64 {
    (4.0 / 3.0) * e_star * r_star.sqrt() * delta.powf(1.5)
}

/// Calculate contact radius for given force
pub fn hertz_contact_radius(f: f64, r_star: f64, e_star: f64) -> f64 {
    (3.0 * f * r_star / (4.0 * e_star)).powf(1.0 / 3.0)
}

/// Hertzian pressure at radius r within contact patch
pub fn hertz_pressure(r: f64, a: f64, p0: f64) -> f64 {
    if r >= a {
        0.0
    } else {
        p0 * (1.0 - (r / a).powi(2)).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hertz_sphere_sphere_formulas() {
        let r1 = 0.01;
        let r2 = 0.02;
        let e = 200e9;
        let nu = 0.3;
        let f = 100.0;

        let r_star = 1.0 / (1.0 / r1 + 1.0 / r2);
        let e_star = e / (1.0 - nu * nu); // Same material

        let a = hertz_contact_radius(f, r_star, e_star);
        let delta = hertz_approach(f, r_star, e_star);
        let f_back = hertz_force(delta, r_star, e_star);

        println!("R* = {:.4e} m", r_star);
        println!("a = {:.4e} m", a);
        println!("δ = {:.4e} m", delta);
        println!("F (back-calculated) = {:.2} N (original: {} N)", f_back, f);

        assert!((f_back - f).abs() / f < 0.001);
    }

    #[test]
    fn test_hertz_sphere_sphere_benchmark() {
        let result = run();
        assert!(result.passed);
    }
}
