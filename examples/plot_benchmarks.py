#!/usr/bin/env python3
"""
FEA Benchmark Report Generator
Generates plots and a clean HTML report for benchmark results.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
from datetime import datetime
import re

# Configure matplotlib
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
})


def load_results(json_path):
    with open(json_path) as f:
        return json.load(f)


def get_best_error(metrics):
    """Get finest mesh error."""
    mesh_priority = ['ultra_fine', 'very_fine', 'fine', 'refine_ultra_fine', 'refine_very_fine', 
                     'refine_fine', 'medium', 'refine_medium', 'coarse', 'refine_coarse']
    
    for level in mesh_priority:
        for m in metrics:
            if level in m['name']:
                if 'sigma_z' in m['name']:
                    continue
                if m['analytical'] == m['computed'] and m['relative_error'] == 0:
                    continue
                return abs(m['relative_error'])
    
    # Fallback
    for m in metrics:
        if m['analytical'] == m['computed']:
            continue
        if 'sigma_z' in m['name']:
            continue
        return abs(m['relative_error'])
    return 0


def extract_convergence(metrics):
    """Extract convergence data by mesh level."""
    mesh_sizes = {
        'coarse': 1.0, 'medium': 0.5, 'fine': 0.25, 
        'very_fine': 0.125, 'ultra_fine': 0.0625
    }
    grouped = defaultdict(dict)
    
    for m in metrics:
        name = m['name']
        if 'sigma_z' in name or (m['analytical'] == m['computed'] and m['relative_error'] == 0):
            continue
        
        for level in ['refine_ultra_fine', 'refine_very_fine', 'refine_fine', 'refine_medium', 'refine_coarse',
                      'ultra_fine', 'very_fine', 'fine', 'medium', 'coarse']:
            if name.startswith(level):
                norm_level = level.replace('refine_', '')
                base = name.replace(f'{level}_', '').replace('refine_', '')
                if norm_level in mesh_sizes:
                    grouped[base][norm_level] = abs(m['relative_error'])
                break
    
    return grouped


def plot_summary(results, output_path):
    """Bar chart of errors by benchmark."""
    fig, ax = plt.subplots(figsize=(10, 5))
    
    names = []
    errors = []
    times = []
    
    for r in results:
        name = r['name'].replace('(3D Block)', '').replace('(3D Solids)', '')
        name = name.replace('Hertz Sphere-on-', 'Hertz ').replace(' in Infinite Solid', '')
        name = name.replace(' Half-Space', '').replace('(Circular/Square)', '').strip()
        names.append(name)
        errors.append(get_best_error(r['metrics']) * 100)
        times.append(r['elapsed_ms'])
    
    x = np.arange(len(names))
    colors = ['#2ecc71' if r['passed'] else '#e74c3c' for r in results]
    
    bars = ax.bar(x, errors, color=colors, edgecolor='white')
    
    for bar, err, t in zip(bars, errors, times):
        label = f'{err:.1f}%' if err > 0.1 else f'{err:.2e}%'
        ax.annotate(label, xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                   xytext=(0, 3), textcoords='offset points', ha='center', fontsize=8)
    
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Error (%)')
    ax.set_title('Benchmark Errors (Finest Mesh)')
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 200)
    ax.axhline(1, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(10, color='gray', linestyle=':', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_convergence_grid(results, output_path):
    """Grid of convergence plots."""
    conv_data = [(r, extract_convergence(r['metrics'])) for r in results]
    conv_data = [(r, c) for r, c in conv_data if c]
    
    if not conv_data:
        return
    
    n = len(conv_data)
    cols = 3
    rows = (n + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(12, 3.5 * rows))
    axes = np.atleast_2d(axes)
    
    mesh_sizes = {'coarse': 1.0, 'medium': 0.5, 'fine': 0.25, 'very_fine': 0.125, 'ultra_fine': 0.0625}
    
    idx = 0
    for result, conv in conv_data:
        row, col = idx // cols, idx % cols
        ax = axes[row, col]
        
        has_data = False
        for metric_name, level_errors in conv.items():
            if len(level_errors) < 2:
                continue
            
            has_data = True
            xs, ys = [], []
            for lvl in ['coarse', 'medium', 'fine', 'very_fine', 'ultra_fine']:
                if lvl in level_errors:
                    xs.append(mesh_sizes[lvl])
                    ys.append(level_errors[lvl] * 100)
            
            if len(xs) >= 2:
                pairs = sorted(zip(xs, ys))
                xs, ys = zip(*pairs)
                ax.loglog(xs, ys, 'o-', lw=2, ms=6, label=metric_name.replace('_', ' ').title())
        
        if has_data:
            # Reference lines
            xr = np.array([0.05, 1.2])
            ax.loglog(xr, 10 * xr, '--', color='gray', alpha=0.4, lw=1, label='O(h)')
            ax.loglog(xr, 10 * xr**2, ':', color='gray', alpha=0.4, lw=1, label='O(h²)')
            
            name = result['name'].split('(')[0].strip()
            ax.set_title(name, fontsize=10)
            ax.set_xlabel('Mesh Size (h)')
            ax.set_ylabel('Error (%)')
            ax.legend(fontsize=7, loc='best')
            ax.grid(True, alpha=0.3)
            idx += 1
        else:
            ax.set_visible(False)
    
    for i in range(idx, rows * cols):
        axes[i // cols, i % cols].set_visible(False)
    
    fig.suptitle('Mesh Convergence (log-log)', fontsize=12, y=1.01)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def generate_report(results, output_dir):
    """Generate clean HTML report."""
    output_dir = Path(output_dir)
    
    passed = sum(1 for r in results if r['passed'])
    total = len(results)
    total_time = sum(r['elapsed_ms'] for r in results)
    
    # Find VTK images
    vtk_imgs = {}
    for name in ['torsion_shaft', 'cantilever_beam', 'spherical_cavity', 'hertz_sphere_flat',
                 'torsion_shaft_explicit', 'contact_bar_block']:
        for suffix in ['ultra_fine', 'very_fine', 'fine', 'medium', '']:
            if suffix:
                combined = output_dir / f'{name}_{suffix}_combined.png'
            else:
                combined = output_dir / f'{name}_combined.png'
            if combined.exists():
                vtk_imgs[name] = combined.name
                break
    
    def fmt_error(e):
        return f'{e*100:.2f}%' if e > 0.001 else f'{e:.2e}'
    
    def fmt_time(ms):
        if ms < 1000:
            return f'{ms:.0f} ms'
        return f'{ms/1000:.2f} s'
    
    html = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>RustFEA Benchmark Report</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; 
       max-width: 1200px; margin: 0 auto; padding: 2rem; background: #fff; color: #333; line-height: 1.5; }}
h1 {{ border-bottom: 2px solid #333; padding-bottom: 0.5rem; }}
h2 {{ margin-top: 2rem; color: #555; }}
h3 {{ margin-top: 1.5rem; color: #666; }}
table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.9rem; }}
th, td {{ border: 1px solid #ddd; padding: 0.5rem 0.75rem; text-align: left; }}
th {{ background: #f5f5f5; }}
.pass {{ color: #27ae60; }}
.fail {{ color: #e74c3c; }}
.mono {{ font-family: 'SF Mono', Monaco, monospace; font-size: 0.85rem; }}
img {{ max-width: 100%; height: auto; margin: 1rem 0; border: 1px solid #eee; }}
.summary {{ background: #f9f9f9; padding: 1rem; border-radius: 4px; margin: 1rem 0; }}
.summary span {{ margin-right: 2rem; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 1rem; }}
.note {{ color: #666; font-size: 0.9rem; font-style: italic; }}
</style>
</head>
<body>

<h1>RustFEA Benchmark Report</h1>

<div class="summary">
<span><strong>Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}</span>
<span><strong>Passed:</strong> {passed}/{total}</span>
<span><strong>Total Time:</strong> {fmt_time(total_time)}</span>
</div>

<h2>Summary</h2>

<img src="benchmark_summary.png" alt="Benchmark Summary">

<table>
<tr><th>#</th><th>Benchmark</th><th>Status</th><th>Best Error</th><th>Time</th></tr>
'''
    
    for i, r in enumerate(results, 1):
        err = get_best_error(r['metrics'])
        status = '<span class="pass">PASS</span>' if r['passed'] else '<span class="fail">FAIL</span>'
        html += f'<tr><td>{i}</td><td>{r["name"]}</td><td>{status}</td>'
        html += f'<td class="mono">{fmt_error(err)}</td><td class="mono">{fmt_time(r["elapsed_ms"])}</td></tr>\n'
    
    html += '</table>\n'
    
    # Convergence
    html += '''
<h2>Mesh Convergence</h2>
<img src="convergence_grid.png" alt="Convergence Grid">
'''
    
    # Field Visualizations
    if vtk_imgs:
        html += '<h2>Field Visualizations</h2>\n'
        
        if 'torsion_shaft' in vtk_imgs:
            html += f'''
<h3>Torsion Shaft</h3>
<img src="{vtk_imgs['torsion_shaft']}" alt="Torsion Shaft">
<p class="note">Displacement magnitude (left) and von Mises stress (right).</p>
'''
        
        if 'cantilever_beam' in vtk_imgs:
            html += f'''
<h3>Cantilever Beam</h3>
<img src="{vtk_imgs['cantilever_beam']}" alt="Cantilever Beam">
<p class="note">Displacement magnitude (left) and von Mises stress (right).</p>
'''
        
        if 'spherical_cavity' in vtk_imgs:
            html += f'''
<h3>Spherical Cavity</h3>
<img src="{vtk_imgs['spherical_cavity']}" alt="Spherical Cavity">
<p class="note">Displacement magnitude (left) and von Mises stress (right).</p>
'''
        
        if 'hertz_sphere_flat' in vtk_imgs:
            html += f'''
<h3>Hertz Sphere-on-Flat Contact</h3>
<img src="{vtk_imgs['hertz_sphere_flat']}" alt="Hertz Contact">
<p class="note">Displacement magnitude (left) and von Mises stress (right). Contact region with prescribed Hertzian approach.</p>
'''
        
        if 'torsion_shaft_explicit' in vtk_imgs:
            html += f'''
<h3>Torsion Shaft (Explicit Solver)</h3>
<img src="{vtk_imgs['torsion_shaft_explicit']}" alt="Explicit Torsion">
<p class="note">Explicit time integration result. Displacement magnitude (left) and von Mises stress (right).</p>
'''
        
        if 'contact_bar_block' in vtk_imgs:
            html += f'''
<h3>Explicit Contact (Bar-Block)</h3>
<img src="{vtk_imgs['contact_bar_block']}" alt="Explicit Contact">
<p class="note">Bar pushed into block with explicit solver and penalty contact. Displacement (left) and stress (right).</p>
'''
    
    # Detailed Results
    html += '<h2>Detailed Results</h2>\n'
    
    for r in results:
        err = get_best_error(r['metrics'])
        status_class = 'pass' if r['passed'] else 'fail'
        
        html += f'''
<h3>{r["name"]}</h3>
<p>{r["description"]}</p>
<p><strong>Status:</strong> <span class="{status_class}">{"PASS" if r["passed"] else "FAIL"}</span> | 
<strong>Best Error:</strong> {fmt_error(err)} | 
<strong>Time:</strong> {fmt_time(r["elapsed_ms"])}</p>
'''
        
        if r.get('notes'):
            html += f'<p class="note">{r["notes"]}</p>\n'
        
        # Metrics table (show key ones)
        html += '<table><tr><th>Metric</th><th>Analytical</th><th>Computed</th><th>Error</th></tr>\n'
        
        count = 0
        for m in r['metrics']:
            if m['analytical'] == m['computed'] and m['relative_error'] == 0:
                continue
            if 'sigma_z' in m['name'] and m['computed'] == 0:
                continue
            if count >= 12:  # Limit rows
                html += '<tr><td colspan="4">...</td></tr>\n'
                break
            
            err_pct = abs(m['relative_error']) * 100
            
            def fmt_val(v):
                if abs(v) < 1e-6 or abs(v) > 1e6:
                    return f'{v:.3e}'
                return f'{v:.4f}'
            
            html += f'<tr><td>{m["name"]}</td><td class="mono">{fmt_val(m["analytical"])}</td>'
            html += f'<td class="mono">{fmt_val(m["computed"])}</td><td class="mono">{err_pct:.2f}%</td></tr>\n'
            count += 1
        
        html += '</table>\n'
    
    html += '''
</body>
</html>
'''
    
    out_file = output_dir / 'benchmark_report.html'
    with open(out_file, 'w') as f:
        f.write(html)
    print(f'Saved: {out_file}')


def main():
    script_dir = Path(__file__).parent
    json_path = script_dir / 'output' / 'benchmark_results.json'
    output_dir = script_dir / 'output' / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f'Loading: {json_path}')
    results = load_results(json_path)
    
    print('Generating plots...')
    plot_summary(results, output_dir / 'benchmark_summary.png')
    plot_convergence_grid(results, output_dir / 'convergence_grid.png')
    
    print('Generating report...')
    generate_report(results, output_dir)
    
    print(f'\nDone. Open {output_dir / "benchmark_report.html"}')


if __name__ == '__main__':
    main()
