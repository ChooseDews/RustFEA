#!/usr/bin/env python3
"""
Plot solver scaling benchmark results comparing UMFPACK vs faer.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load results
results_path = Path(__file__).parent / "output" / "solver_scaling_results.json"
with open(results_path) as f:
    results = json.load(f)

# Extract data
dofs = [r['num_dofs'] for r in results]
umfpack_times = [r['umfpack_time_ms'] for r in results]
faer_times = [r['faer_time_ms'] for r in results]
speedups = [r['speedup'] for r in results]

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Solve time vs DOFs (log-log)
ax1.loglog(dofs, umfpack_times, 'o-', color='#1f77b4', linewidth=2, markersize=8, label='UMFPACK')
ax1.loglog(dofs, faer_times, 's-', color='#ff7f0e', linewidth=2, markersize=8, label='faer (Cholesky)')
ax1.set_xlabel('Degrees of Freedom (DOFs)', fontsize=12)
ax1.set_ylabel('Solve Time (ms)', fontsize=12)
ax1.set_title('Solver Performance: UMFPACK vs faer', fontsize=14, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(True, which='both', alpha=0.3)

# Add reference O(n^2) and O(n^1.5) lines
dofs_ref = np.array(dofs)
base_time = faer_times[len(faer_times)//2]
base_dof = dofs[len(dofs)//2]
# O(n^1.5) reference
scale_15 = base_time / (base_dof ** 1.5)
ref_15 = scale_15 * (dofs_ref ** 1.5)
ax1.loglog(dofs_ref, ref_15, '--', color='gray', alpha=0.5, label='O(n^1.5)')
# O(n^2) reference
scale_2 = base_time / (base_dof ** 2)
ref_2 = scale_2 * (dofs_ref ** 2)
ax1.loglog(dofs_ref, ref_2, ':', color='gray', alpha=0.5, label='O(n^2)')
ax1.legend(fontsize=10, loc='upper left')

# Plot 2: Speedup vs DOFs
ax2.semilogx(dofs, speedups, 'D-', color='#2ca02c', linewidth=2, markersize=8)
ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Equal performance')
ax2.fill_between(dofs, 1.0, speedups, alpha=0.3, color='#2ca02c', where=[s > 1 for s in speedups])
ax2.set_xlabel('Degrees of Freedom (DOFs)', fontsize=12)
ax2.set_ylabel('Speedup (UMFPACK / faer)', fontsize=12)
ax2.set_title('faer Speedup over UMFPACK', fontsize=14, fontweight='bold')
ax2.grid(True, which='both', alpha=0.3)
ax2.set_ylim(0, max(speedups) * 1.1)

# Add annotations for key points
max_speedup_idx = speedups.index(max(speedups))
ax2.annotate(f'{speedups[max_speedup_idx]:.2f}x faster', 
             xy=(dofs[max_speedup_idx], speedups[max_speedup_idx]),
             xytext=(dofs[max_speedup_idx]*1.5, speedups[max_speedup_idx]*0.9),
             fontsize=10, fontweight='bold',
             arrowprops=dict(arrowstyle='->', color='#2ca02c'))

plt.tight_layout()

# Save plot
output_path = Path(__file__).parent / "output" / "plots" / "solver_scaling_comparison.png"
output_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"Plot saved to {output_path}")

# Also save as PDF for high quality
pdf_path = output_path.with_suffix('.pdf')
plt.savefig(pdf_path, bbox_inches='tight')
print(f"PDF saved to {pdf_path}")

plt.show()

# Print summary statistics
print("\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)
print(f"DOF Range: {min(dofs):,} - {max(dofs):,}")
print(f"UMFPACK Time Range: {min(umfpack_times):.1f} - {max(umfpack_times):.1f} ms")
print(f"faer Time Range: {min(faer_times):.1f} - {max(faer_times):.1f} ms")
print(f"Speedup Range: {min(speedups):.2f}x - {max(speedups):.2f}x")
print(f"Average Speedup: {np.mean(speedups):.2f}x")
print(f"Max Speedup at: {dofs[max_speedup_idx]:,} DOFs ({results[max_speedup_idx]['mesh_size']})")
