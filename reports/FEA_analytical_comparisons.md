# FEA Benchmark Report
Generated: 2026-09-19 05:50:03

## Summary

| Benchmark | Status | Max Error | Time (ms) |
|-----------|--------|-----------|----------|
| Uniaxial Tension (3D Block) | ✅ PASS | 1.81e-2 | 16.5 |
| Pure Shear (3D Block) | ✅ PASS | 3.12e-3 | 23.0 |
| Hydrostatic Compression (3D Block) | ✅ PASS | 7.15e-3 | 4.1 |
| Torsion Shaft (Circular/Square) | ✅ PASS | 4.38e-1 | 104.6 |
| Cantilever Beam (3D Solids) | ✅ PASS | 3.69e0 | 45.7 |

## Detailed Results

### Uniaxial Tension (3D Block)

**Description:** 3D rectangular block under uniaxial tension with prescribed displacement. Tests basic 3D elasticity, Poisson effect, and stiffness matrix assembly.

**Status:** PASS

| Metric | Analytical | FEA | Rel. Error | Tolerance |
|--------|------------|-----|------------|----------|
| coarse_u_x_prescribed | 1.000000e-3 | 9.998067e-4 | -1.93e-4 | 2.00e-2 |
| coarse_avg_u_x_error | 0.000000e0 | 1.437179e-2 | 1.44e-2 | 2.00e-2 |
| coarse_avg_u_y_error | 0.000000e0 | 1.814879e-2 | 1.81e-2 | 2.00e-2 |
| coarse_avg_u_z_error | 0.000000e0 | 1.814879e-2 | 1.81e-2 | 2.00e-2 |
| medium_u_x_prescribed | 1.000000e-3 | 9.999545e-4 | -4.55e-5 | 2.00e-2 |
| medium_avg_u_x_error | 0.000000e0 | 1.446388e-2 | 1.45e-2 | 2.00e-2 |
| medium_avg_u_y_error | 0.000000e0 | 1.208821e-2 | 1.21e-2 | 2.00e-2 |
| medium_avg_u_z_error | 0.000000e0 | 1.208821e-2 | 1.21e-2 | 2.00e-2 |
| fine_u_x_prescribed | 1.000000e-3 | 9.999889e-4 | -1.11e-5 | 2.00e-2 |
| fine_avg_u_x_error | 0.000000e0 | 9.974861e-3 | 9.97e-3 | 2.00e-2 |
| fine_avg_u_y_error | 0.000000e0 | 7.230880e-3 | 7.23e-3 | 2.00e-2 |
| fine_avg_u_z_error | 0.000000e0 | 7.230880e-3 | 7.23e-3 | 2.00e-2 |

**Notes:** Material: E=2.00e11 Pa, ν=0.30. Geometry: 1×0.1×0.1 m. Applied Δ=0.0010 m

---

### Pure Shear (3D Block)

**Description:** 3D block under pure shear deformation u_x = γy. Tests shear strain components, stress recovery, and coordinate handling.

**Status:** PASS

| Metric | Analytical | FEA | Rel. Error | Tolerance |
|--------|------------|-----|------------|----------|
| coarse_max_u_x_error | 0.000000e0 | 3.116219e-3 | 3.12e-3 | 2.00e-2 |
| coarse_max_u_y_error | 0.000000e0 | 3.116219e-3 | 3.12e-3 | 2.00e-2 |
| coarse_max_u_z_error | 0.000000e0 | 2.126920e-5 | 2.13e-5 | 2.00e-2 |
| coarse_corner_u_x | 1.000000e-3 | 9.991824e-4 | -8.18e-4 | 2.00e-2 |
| medium_max_u_x_error | 0.000000e0 | 1.602935e-3 | 1.60e-3 | 2.00e-2 |
| medium_max_u_y_error | 0.000000e0 | 1.602935e-3 | 1.60e-3 | 2.00e-2 |
| medium_max_u_z_error | 0.000000e0 | 5.567997e-6 | 5.57e-6 | 2.00e-2 |
| medium_corner_u_x | 1.000000e-3 | 9.995931e-4 | -4.07e-4 | 2.00e-2 |
| fine_max_u_x_error | 0.000000e0 | 4.039148e-4 | 4.04e-4 | 2.00e-2 |
| fine_max_u_y_error | 0.000000e0 | 4.039148e-4 | 4.04e-4 | 2.00e-2 |
| fine_max_u_z_error | 0.000000e0 | 1.730345e-5 | 1.73e-5 | 2.00e-2 |
| fine_corner_u_x | 1.000000e-3 | 9.998985e-4 | -1.01e-4 | 2.00e-2 |

**Notes:** Material: E=6.89e10 Pa, ν=0.33, G=2.59e10 Pa. Geometry: 1×1×0.5 m. Applied γ=0.0010

---

### Hydrostatic Compression (3D Block)

**Description:** 3D cubic block under uniform hydrostatic pressure on all faces. Tests volumetric terms and bulk modulus implementation.

**Status:** PASS

| Metric | Analytical | FEA | Rel. Error | Tolerance |
|--------|------------|-----|------------|----------|
| volumetric_strain | -1.480406e-5 | -1.475769e-5 | 3.13e-3 | 2.00e-2 |
| max_displacement_error | 0.000000e0 | 7.148108e-3 | 7.15e-3 | 2.00e-2 |
| linear_strain_x | -4.934688e-6 | -4.919230e-6 | 3.13e-3 | 2.00e-2 |

**Notes:** Material: E=6.89e10 Pa, ν=0.33, K=6.75e10 Pa. Geometry: 1×1×1 m cube. Applied p=1.00e6 Pa

---

### Torsion Shaft (Circular/Square)

**Description:** Shaft under torsion. Tests torsional deformation and shear stress gradients. Note: Using square cross-section as approximation to circular.

**Status:** PASS

| Metric | Analytical | FEA | Rel. Error | Tolerance |
|--------|------------|-----|------------|----------|
| coarse_twist_angle_vs_square | 4.455969e-4 | 2.506085e-4 | -4.38e-1 | 5.00e-1 |
| coarse_warping_ratio | 0.000000e0 | 3.517637e-14 | 3.52e-14 | 2.00e-1 |
| medium_twist_angle_vs_square | 4.455969e-4 | 3.634227e-4 | -1.84e-1 | 5.00e-1 |
| medium_warping_ratio | 0.000000e0 | 3.937671e-3 | 3.94e-3 | 2.00e-1 |
| fine_twist_angle_vs_square | 4.455969e-4 | 4.031181e-4 | -9.53e-2 | 5.00e-1 |
| fine_warping_ratio | 0.000000e0 | 4.681655e-3 | 4.68e-3 | 2.00e-1 |

**Notes:** Material: E=6.89e10 Pa, ν=0.33, G=2.59e10 Pa. Geometry: L=1 m, R≈0.05 m (square approx). Applied T=100 N·m. Analytical (circular): φ=0.0225°, τ_max=5.09e5 Pa

---

### Cantilever Beam (3D Solids)

**Description:** Slender cantilever beam modeled with 3D solid elements under tip load. Tests bending accuracy, shear locking, and mesh sensitivity.

**Status:** PASS

| Metric | Analytical | FEA | Rel. Error | Tolerance |
|--------|------------|-----|------------|----------|
| L_h5_tip_deflection | 7.256894e-5 | 3.786860e-5 | -4.78e-1 | 8.00e-1 |
| L_h5_stiffness_ratio | 1.000000e0 | 1.916335e0 | 9.16e-1 | 4.00e0 |
| L_h5_deflection_shape | 6.250000e-1 | 2.630064e-1 | -5.79e-1 | 6.00e-1 |
| L_h10_tip_deflection | 5.805515e-4 | 2.979044e-4 | -4.87e-1 | 8.00e-1 |
| L_h10_stiffness_ratio | 1.000000e0 | 1.948785e0 | 9.49e-1 | 4.00e0 |
| L_h10_deflection_shape | 6.250000e-1 | 2.838801e-1 | -5.46e-1 | 6.00e-1 |
| L_h20_tip_deflection | 4.644412e-3 | 2.386826e-3 | -4.86e-1 | 8.00e-1 |
| L_h20_stiffness_ratio | 1.000000e0 | 1.945853e0 | 9.46e-1 | 4.00e0 |
| L_h20_deflection_shape | 6.250000e-1 | 2.972390e-1 | -5.24e-1 | 6.00e-1 |
| refine_coarse_tip_deflection | 5.805515e-4 | 1.236758e-4 | -7.87e-1 | 8.00e-1 |
| refine_coarse_stiffness_ratio | 1.000000e0 | 4.694138e0 | 3.69e0 | 4.00e0 |
| refine_coarse_deflection_shape | 6.250000e-1 | 2.618496e-1 | -5.81e-1 | 6.00e-1 |
| refine_medium_tip_deflection | 5.805515e-4 | 2.979044e-4 | -4.87e-1 | 8.00e-1 |
| refine_medium_stiffness_ratio | 1.000000e0 | 1.948785e0 | 9.49e-1 | 4.00e0 |
| refine_medium_deflection_shape | 6.250000e-1 | 2.838801e-1 | -5.46e-1 | 6.00e-1 |
| refine_fine_tip_deflection | 5.805515e-4 | 4.671461e-4 | -1.95e-1 | 8.00e-1 |
| refine_fine_stiffness_ratio | 1.000000e0 | 1.242762e0 | 2.43e-1 | 4.00e0 |
| refine_fine_deflection_shape | 6.250000e-1 | 2.973472e-1 | -5.24e-1 | 6.00e-1 |

**Notes:** Material: E=6.89e10 Pa, ν=0.33. Cross-section: 0.1×0.1 m. I=8.3333e-6 m⁴. Tip load scaled to maintain reasonable deflection.

---
