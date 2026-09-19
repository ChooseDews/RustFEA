# RustFEA Benchmark Suite

Comprehensive analytical benchmarks for validating FEA solver accuracy against exact solutions.

## Quick Start

```bash
# Run all benchmarks
./run_benchmarks.sh

# Or with cargo directly
cargo run --bin run_benchmarks -- --output examples/output/benchmark_report.md
```

## Benchmark Suite Overview

### Level 1: Fundamental Element Verification
These test uniform stress/strain states that 8-node brick elements should reproduce nearly exactly.

| # | Benchmark | Analytical Solution | What it Proves | Expected Error |
|---|-----------|---------------------|----------------|----------------|
| 1 | **Uniaxial Tension** | σ = F/A, ε = σ/E, u = ε·x | Basic 3D elasticity, Poisson effect | < 2% |
| 2 | **Pure Shear** | τ = Gγ, G = E/[2(1+ν)] | Shear terms, off-diagonal strain | < 1% |
| 3 | **Hydrostatic** | σ_kk = -3p, ε_v = -p/K | Volumetric response, bulk modulus | < 1% |

### Level 2: 3D Continuum Verification
These test convergence with nonuniform stress fields.

| # | Benchmark | Analytical Solution | What it Proves | Expected Error |
|---|-----------|---------------------|----------------|----------------|
| 4 | **Torsion Shaft** | φ = TL/(GJ), τ = TR/J | 3D torsion, shear stress gradient | 10-45% (converges) |
| 5 | **Hollow Sphere** | σ_r = A - B/r³ (Lamé) | Curved geometry, pressure loading | (mesh bug - disabled) |
| 6 | **Spherical Cavity** | σ_θ = -1.5p∞ at surface | 3D stress concentration (1.5×) | ~3% avg, ~60% max |
| 7 | **Boussinesq** | σ_z = -3Pz³/(2πR⁵) | Point load singularity, decay | High (singular) |
| 8 | **Cantilever Beam** | δ = PL³/(3EI) | Bending accuracy, shear locking | 20-80% (shear locking) |

### Level 3: Contact Verification
Validate Hertzian contact theory predictions.

| # | Benchmark | Analytical Solution | What it Proves | Expected Error |
|---|-----------|---------------------|----------------|----------------|
| 9 | **Hertz Sphere-Flat** | a = (3FR/4E*)^(1/3) | Contact mechanics, approach | < 10% |
| 10 | **Hertz Sphere-Sphere** | 1/R* = 1/R₁ + 1/R₂ | Curved surface contact | Analytical (exact) |

## Running Benchmarks

### Command Line Options

```bash
# Run all benchmarks
cargo run --bin run_benchmarks

# Output to file (markdown, json, or csv)
cargo run --bin run_benchmarks -- --output report.md
cargo run --bin run_benchmarks -- --format json --output report.json
cargo run --bin run_benchmarks -- --format csv --output results.csv

# Run specific benchmarks (comma-separated filter)
cargo run --bin run_benchmarks -- --benchmarks uniaxial,shear,hertz

# With verbose logging
RUST_LOG=info cargo run --bin run_benchmarks
```

### Shell Script

```bash
./run_benchmarks.sh                           # Run all
./run_benchmarks.sh -o report.md              # Save to file
./run_benchmarks.sh -b uniaxial,torsion       # Filter
./run_benchmarks.sh -f csv -o data.csv        # CSV export
```

## Understanding the Report

### Output Format

The markdown report includes:
1. **Summary table** - Pass/fail status, max error, execution time
2. **Detailed results** - Per-metric comparison with analytical values
3. **Convergence plots** - ASCII visualization of mesh refinement
4. **Theory appendix** - Analytical solution formulas

### Reading the Metrics

| Column | Meaning |
|--------|---------|
| Analytical | Exact solution value |
| FEA | Computed finite element result |
| Error | Relative error: (FEA - Analytical) / |Analytical| |
| Tolerance | Maximum acceptable error for PASS |
| Status | ✓ = within tolerance, ✗ = exceeds tolerance |

### Units

Values are displayed with SI prefixes for readability:
- **GPa, MPa, kPa** - Pressure/stress
- **mm, μm, nm** - Displacement
- **kN, N, mN** - Force

## Understanding Errors

### Shear Locking (Cantilever Beam)

8-node brick elements are known to exhibit "shear locking" in bending-dominated problems:
- The element is too stiff, underpredicting deflection
- Stiffness ratio > 1.0 indicates overstiffness
- Fine meshes reduce this effect but don't eliminate it
- Proper fix requires incompatible modes or reduced integration

### Stress Singularities (Boussinesq)

Point loads create 1/r stress singularities that FEA cannot capture:
- Near-field stresses are always inaccurate
- Mesh refinement helps but never converges
- Focus on far-field displacements for validation

### Domain Truncation (Spherical Cavity)

Infinite domain problems require finite approximations:
- Larger domain ratios (b/a = 10+) improve accuracy
- Boundary effects appear at corners of block mesh
- Average errors are typically much smaller than max errors

## Adding New Benchmarks

1. Create `src/benchmarks/your_benchmark.rs`:

```rust
use crate::benchmarks::{BenchmarkResult, MetricComparison};

pub fn run() -> BenchmarkResult {
    let mut result = BenchmarkResult::new(
        "Benchmark Name",
        "Description of what this tests",
    );
    
    // Run simulation(s)...
    
    // Add metrics
    result.add_metric(MetricComparison::new(
        "metric_name",
        analytical_value,
        computed_value,
        tolerance,  // e.g., 0.05 for 5%
    ));
    
    result.set_notes("Material and geometry parameters");
    result
}
```

2. Add to `src/benchmarks/mod.rs`:
```rust
pub mod your_benchmark;
```

3. Register in `src/bin/run_benchmarks.rs`:
```rust
("your_benchmark", your_benchmark::run),
```

## Mesh Utilities

`src/benchmarks/mesh_utils.rs` provides:

- `generate_block_mesh(lx, ly, lz, nx, ny, nz)` - Rectangular block
- `generate_beam_mesh(...)` - Oriented beam
- `generate_hollow_sphere_sector_mesh(...)` - Spherical sector (has issues)

All meshes include named node groups for boundary conditions:
- `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max` - Face nodes
- `inner`, `outer` - For hollow geometries

## Known Limitations

1. **Hollow sphere mesh** - Element connectivity is incorrect; needs fixing
2. **Pressure BC** - Not implemented; use prescribed displacement instead
3. **Contact benchmarks** - Full contact requires explicit solver
4. **Stress recovery** - Not implemented; only displacement comparisons

## References

1. Timoshenko & Goodier, *Theory of Elasticity*, 3rd ed.
2. Hertz (1881), "On the contact of elastic solids"
3. Boussinesq (1885), *Applications des potentiels*
4. Zienkiewicz & Taylor, *The Finite Element Method*, Vol. 1
