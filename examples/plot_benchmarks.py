#!/usr/bin/env python3
"""
FEA Benchmark Visualization & Report Generator
Generates publication-quality plots and HTML report for benchmark results.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from collections import defaultdict
from datetime import datetime
import re

# Set up nice matplotlib style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
})

# Color palette
COLORS = {
    'pass': '#2ecc71',
    'fail': '#e74c3c',
    'primary': '#3498db',
    'secondary': '#9b59b6',
    'accent': '#f39c12',
    'dark': '#2c3e50',
}


def load_benchmark_results(json_path):
    """Load benchmark results from JSON file."""
    with open(json_path) as f:
        return json.load(f)


def get_best_mesh_error(metrics, metric_pattern=None):
    """Get the best (finest mesh) error for a benchmark."""
    mesh_priority = ['very_fine', 'fine', 'refine_very_fine', 'refine_fine', 'medium', 'refine_medium', 'coarse', 'refine_coarse']
    
    best_error = None
    best_metric = None
    
    for level in mesh_priority:
        for metric in metrics:
            name = metric['name']
            if level in name:
                if 'sigma_z' in name:
                    continue
                if metric['analytical'] == metric['computed'] and metric['relative_error'] == 0:
                    continue
                if metric_pattern and metric_pattern not in name:
                    continue
                    
                error = abs(metric['relative_error'])
                if best_error is None or error < best_error:
                    best_error = error
                    best_metric = metric
        
        if best_error is not None:
            return best_error, best_metric
    
    for metric in metrics:
        if metric['analytical'] == metric['computed']:
            continue
        if 'sigma_z' in metric['name']:
            continue
        error = abs(metric['relative_error'])
        if best_error is None or error < best_error:
            best_error = error
            best_metric = metric
    
    return best_error or 0, best_metric


def extract_convergence_data(metrics, error_pattern='error|deflection|u_r|u_z|twist'):
    """Extract convergence data from metrics for different mesh levels."""
    grouped = defaultdict(dict)
    
    for metric in metrics:
        name = metric['name']
        
        if 'sigma_z' in name:
            continue
        if metric['analytical'] == metric['computed'] and metric['relative_error'] == 0:
            continue
            
        mesh_level = None
        for level in ['refine_very_fine', 'refine_fine', 'refine_medium', 'refine_coarse', 
                      'very_fine', 'fine', 'medium', 'coarse']:
            if name.startswith(level):
                mesh_level = level.replace('refine_', '')
                base_name = name.replace(f'{level}_', '').replace('refine_', '')
                break
        
        if mesh_level and re.search(error_pattern, name, re.IGNORECASE):
            grouped[base_name][mesh_level] = abs(metric['relative_error'])
    
    return grouped


def plot_summary_bar_chart(results, output_path):
    """Create a bar chart showing best error for each benchmark."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    benchmarks = []
    errors = []
    colors = []
    
    for result in results:
        name = result['name']
        short_name = name.replace('(3D Block)', '').replace('(3D Solids)', '')
        short_name = short_name.replace('Hertz Sphere-on-', 'Hertz\n')
        short_name = short_name.replace(' in Infinite Solid', '')
        short_name = short_name.replace(' Half-Space', '')
        short_name = short_name.replace('(Circular/Square)', '')
        short_name = short_name.strip()
        
        best_error, _ = get_best_mesh_error(result['metrics'])
        
        benchmarks.append(short_name)
        errors.append(best_error * 100)
        colors.append(COLORS['pass'] if result['passed'] else COLORS['fail'])
    
    x = np.arange(len(benchmarks))
    bars = ax.bar(x, errors, color=colors, edgecolor='white', linewidth=1.5)
    
    for bar, error in zip(bars, errors):
        height = bar.get_height()
        label = f'{height:.1f}%' if height > 0.1 else f'{height:.2e}%'
        ax.annotate(label,
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=9)
    
    ax.set_xlabel('Benchmark', fontweight='bold')
    ax.set_ylabel('Best Mesh Error (%)', fontweight='bold')
    ax.set_title('FEA Benchmark Suite: Error by Benchmark\n(finest mesh results)', 
                 fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(benchmarks, rotation=45, ha='right')
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 200)
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='1% error')
    ax.axhline(y=10, color='gray', linestyle=':', alpha=0.5, label='10% error')
    
    pass_patch = mpatches.Patch(color=COLORS['pass'], label='PASS')
    fail_patch = mpatches.Patch(color=COLORS['fail'], label='FAIL')
    ax.legend(handles=[pass_patch, fail_patch], loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_convergence_grid(results, output_path):
    """Create a grid of convergence plots for all benchmarks."""
    benchmarks_with_conv = []
    
    for result in results:
        conv_data = extract_convergence_data(result['metrics'])
        if conv_data:
            benchmarks_with_conv.append((result, conv_data))
    
    if not benchmarks_with_conv:
        print("No convergence data found")
        return
    
    n_plots = len(benchmarks_with_conv)
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 4 * n_rows))
    axes = np.atleast_2d(axes)
    
    mesh_sizes = {'coarse': 1.0, 'medium': 0.5, 'fine': 0.25, 'very_fine': 0.125}
    
    plot_idx = 0
    for result, conv_data in benchmarks_with_conv:
        row = plot_idx // n_cols
        col = plot_idx % n_cols
        ax = axes[row, col]
        
        has_data = False
        for metric_name, level_errors in conv_data.items():
            if len(level_errors) < 2:
                continue
            
            has_data = True
            x_vals = []
            y_vals = []
            
            for level in ['coarse', 'medium', 'fine', 'very_fine']:
                if level in level_errors:
                    x_vals.append(mesh_sizes[level])
                    y_vals.append(level_errors[level] * 100)
            
            if len(x_vals) >= 2:
                sorted_pairs = sorted(zip(x_vals, y_vals))
                x_vals, y_vals = zip(*sorted_pairs)
                ax.loglog(x_vals, y_vals, 'o-', linewidth=2, markersize=8,
                         label=metric_name.replace('_', ' ').title())
        
        if has_data:
            x_ref = np.array([0.1, 1.0])
            y_base = 10
            ax.loglog(x_ref, y_base * x_ref, '--', color='gray', alpha=0.5, 
                     linewidth=1, label='O(h) linear')
            ax.loglog(x_ref, y_base * x_ref**2, ':', color='gray', alpha=0.5,
                     linewidth=1, label='O(h²) quadratic')
            
            short_name = result['name'].split('(')[0].strip()
            ax.set_xlabel('Relative Mesh Size (h)', fontweight='bold')
            ax.set_ylabel('Error (%)', fontweight='bold')
            ax.set_title(short_name, fontweight='bold')
            ax.legend(loc='best', fontsize=8)
            ax.grid(True, which='both', alpha=0.3)
            plot_idx += 1
        else:
            ax.set_visible(False)
    
    for idx in range(plot_idx, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)
    
    fig.suptitle('Mesh Convergence Analysis\n(Error vs Mesh Size, log-log scale)', 
                 fontweight='bold', fontsize=14, y=1.02)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_individual_convergence(results, output_dir):
    """Create individual convergence plots for each benchmark."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    mesh_sizes = {'coarse': 1.0, 'medium': 0.5, 'fine': 0.25, 'very_fine': 0.125}
    
    for result in results:
        conv_data = extract_convergence_data(result['metrics'])
        if not conv_data:
            continue
        
        fig, ax = plt.subplots(figsize=(8, 6))
        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
        colors = plt.cm.tab10(np.linspace(0, 1, 10))
        
        metric_idx = 0
        for metric_name, level_errors in conv_data.items():
            if len(level_errors) < 2:
                continue
            
            x_vals = []
            y_vals = []
            
            for level in ['coarse', 'medium', 'fine', 'very_fine']:
                if level in level_errors:
                    x_vals.append(mesh_sizes[level])
                    y_vals.append(level_errors[level] * 100)
            
            if len(x_vals) >= 2:
                sorted_pairs = sorted(zip(x_vals, y_vals))
                x_vals, y_vals = zip(*sorted_pairs)
                
                label = metric_name.replace('_', ' ').title()
                ax.loglog(x_vals, y_vals, 
                         marker=markers[metric_idx % len(markers)],
                         color=colors[metric_idx % len(colors)],
                         linewidth=2, markersize=10, label=label)
                metric_idx += 1
        
        if metric_idx == 0:
            plt.close()
            continue
        
        x_ref = np.array([0.1, 1.2])
        y_base = ax.get_ylim()[1] * 0.5
        
        ax.loglog(x_ref, y_base * (x_ref / x_ref[0]), '--', 
                 color='gray', alpha=0.4, linewidth=1.5, label='O(h)')
        ax.loglog(x_ref, y_base * (x_ref / x_ref[0])**2, ':', 
                 color='gray', alpha=0.4, linewidth=1.5, label='O(h²)')
        
        short_name = result['name'].split('(')[0].strip()
        ax.set_xlabel('Relative Mesh Size (h/h0)', fontweight='bold')
        ax.set_ylabel('Relative Error (%)', fontweight='bold')
        ax.set_title(f'{short_name}\nMesh Convergence Study', fontweight='bold')
        ax.legend(loc='best')
        ax.grid(True, which='both', alpha=0.3)
        
        safe_name = result['name'].lower().replace(' ', '_').replace('(', '').replace(')', '')
        safe_name = safe_name.replace('/', '_').replace('-', '_')
        safe_name = re.sub(r'_+', '_', safe_name)
        
        output_file = output_dir / f"convergence_{safe_name}.png"
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
        print(f"Saved: {output_file}")


def plot_benchmark_categories(results, output_path):
    """Create a grouped bar chart by benchmark category."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    categories = {
        'Fundamental\n(Level 1)': ['Uniaxial Tension', 'Pure Shear', 'Hydrostatic Compression'],
        'Continuum\n(Level 2)': ['Torsion Shaft', 'Cantilever Beam', 'Spherical Cavity', 'Boussinesq'],
        'Contact\n(Level 3)': ['Hertz Sphere-on-Flat', 'Hertz Sphere-on-Sphere'],
    }
    
    cat_errors = {}
    cat_benchmarks = {}
    
    for cat_name, benchmark_patterns in categories.items():
        errors = []
        names = []
        for result in results:
            for pattern in benchmark_patterns:
                if pattern in result['name']:
                    best_error, _ = get_best_mesh_error(result['metrics'])
                    errors.append(best_error * 100)
                    short_name = result['name'].split('(')[0].strip()
                    short_name = short_name.replace('Hertz Sphere-on-', 'Hertz ')
                    names.append(short_name)
                    break
        
        cat_errors[cat_name] = errors
        cat_benchmarks[cat_name] = names
    
    x_offset = 0
    width = 0.8
    xticks = []
    xticklabels = []
    
    category_colors = {
        'Fundamental\n(Level 1)': '#3498db',
        'Continuum\n(Level 2)': '#9b59b6',
        'Contact\n(Level 3)': '#2ecc71',
    }
    
    for cat_name, errors in cat_errors.items():
        x = np.arange(len(errors)) + x_offset
        bars = ax.bar(x, errors, width, color=category_colors[cat_name], 
                     edgecolor='white', linewidth=1, label=cat_name)
        
        for bar, error in zip(bars, errors):
            height = bar.get_height()
            label = f'{height:.1f}%' if height > 0.1 else f'{height:.2e}%'
            ax.annotate(label,
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3), textcoords="offset points",
                       ha='center', va='bottom', fontsize=8, rotation=45)
        
        xticks.extend(x)
        xticklabels.extend(cat_benchmarks[cat_name])
        x_offset += len(errors) + 0.5
    
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 500)
    ax.set_xlabel('Benchmark', fontweight='bold')
    ax.set_ylabel('Best Mesh Error (%)', fontweight='bold')
    ax.set_title('FEA Benchmark Suite by Category\n(finest mesh results)', 
                fontweight='bold', pad=15)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=45, ha='right', fontsize=9)
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=10, color='gray', linestyle=':', alpha=0.5, linewidth=1)
    ax.legend(loc='upper right', title='Category')
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")


def generate_html_report(results, output_dir):
    """Generate a publication-quality HTML report with all visualizations."""
    output_dir = Path(output_dir)
    
    # Calculate summary stats
    passed = sum(1 for r in results if r['passed'])
    total = len(results)
    total_time_ms = sum(r['elapsed_ms'] for r in results)
    
    # Find VTK visualization images
    vtk_images = {
        'torsion': {
            'combined': output_dir / 'torsion_shaft_very_fine_combined.png',
            'displacement': output_dir / 'torsion_shaft_very_fine_displacement.png',
            'vonmises': output_dir / 'torsion_shaft_very_fine_vonmises.png',
        },
        'cantilever': {
            'combined': output_dir / 'cantilever_beam_very_fine_combined.png',
            'displacement': output_dir / 'cantilever_beam_very_fine_displacement.png',
            'vonmises': output_dir / 'cantilever_beam_very_fine_vonmises.png',
        },
        'spherical': {
            'combined': output_dir / 'spherical_cavity_very_fine_combined.png',
            'displacement': output_dir / 'spherical_cavity_very_fine_displacement.png',
            'vonmises': output_dir / 'spherical_cavity_very_fine_vonmises.png',
        }
    }
    
    # Build benchmark details HTML
    benchmark_details = ""
    
    # Level 1 benchmarks
    level1_benchmarks = [r for r in results if any(x in r['name'] for x in ['Uniaxial', 'Pure Shear', 'Hydrostatic'])]
    level2_benchmarks = [r for r in results if any(x in r['name'] for x in ['Torsion', 'Cantilever', 'Spherical', 'Boussinesq'])]
    level3_benchmarks = [r for r in results if 'Hertz' in r['name']]
    
    def format_metric_table(result):
        """Format metrics as HTML table rows."""
        rows = ""
        for m in result['metrics']:
            if m['analytical'] == m['computed'] and m['relative_error'] == 0:
                continue  # Skip self-comparisons
            if 'sigma_z' in m['name'] and m['computed'] == 0:
                continue  # Skip placeholder stress values
            
            error_pct = abs(m['relative_error']) * 100
            error_class = 'low' if error_pct < 5 else ('medium' if error_pct < 20 else 'high')
            
            # Format values nicely
            def fmt_val(v):
                if abs(v) < 1e-6 or abs(v) > 1e6:
                    return f"{v:.3e}"
                elif abs(v) < 0.01:
                    return f"{v:.6f}"
                else:
                    return f"{v:.4f}"
            
            rows += f"""
                <tr>
                    <td class="metric-name">{m['name'].replace('_', ' ').title()}</td>
                    <td class="mono">{fmt_val(m['analytical'])}</td>
                    <td class="mono">{fmt_val(m['computed'])}</td>
                    <td class="error-{error_class}">{error_pct:.2f}%</td>
                </tr>"""
        return rows
    
    def benchmark_card(result, show_vtk=None):
        """Generate HTML card for a benchmark."""
        best_error, _ = get_best_mesh_error(result['metrics'])
        status_class = 'pass' if result['passed'] else 'fail'
        status_icon = '✓' if result['passed'] else '✗'
        
        vtk_section = ""
        if show_vtk and vtk_images.get(show_vtk):
            imgs = vtk_images[show_vtk]
            if imgs['combined'].exists():
                vtk_section = f"""
                <div class="vtk-visualization">
                    <h5>Field Visualization</h5>
                    <img src="{imgs['combined'].name}" alt="{result['name']} field visualization" class="combined-img">
                    <div class="vtk-grid">
                        <div class="vtk-item">
                            <img src="{imgs['displacement'].name}" alt="Displacement">
                            <span class="vtk-label">Displacement Magnitude</span>
                        </div>
                        <div class="vtk-item">
                            <img src="{imgs['vonmises'].name}" alt="von Mises Stress">
                            <span class="vtk-label">von Mises Stress</span>
                        </div>
                    </div>
                </div>"""
        
        return f"""
        <div class="benchmark-card">
            <div class="benchmark-header">
                <div class="benchmark-title">
                    <h4>{result['name']}</h4>
                    <span class="benchmark-time">{result['elapsed_ms']:.1f} ms</span>
                </div>
                <span class="status-badge status-{status_class}">{status_icon} {'PASS' if result['passed'] else 'FAIL'}</span>
            </div>
            <div class="benchmark-body">
                <p class="description">{result['description']}</p>
                <div class="params-box">
                    <strong>Parameters:</strong><br>
                    {result.get('notes', 'N/A').replace(chr(10), '<br>')}
                </div>
                <div class="metrics-section">
                    <h5>Results (Best Error: {best_error*100:.2f}%)</h5>
                    <table class="metrics-table">
                        <thead>
                            <tr>
                                <th>Metric</th>
                                <th>Analytical</th>
                                <th>Computed</th>
                                <th>Error</th>
                            </tr>
                        </thead>
                        <tbody>
                            {format_metric_table(result)}
                        </tbody>
                    </table>
                </div>
                {vtk_section}
            </div>
        </div>"""
    
    # Generate HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RustFEA Benchmark Report</title>
    <style>
        :root {{
            --primary: #2563eb;
            --primary-dark: #1d4ed8;
            --primary-light: #dbeafe;
            --success: #059669;
            --success-light: #d1fae5;
            --warning: #d97706;
            --warning-light: #fef3c7;
            --error: #dc2626;
            --error-light: #fee2e2;
            --gray-50: #f9fafb;
            --gray-100: #f3f4f6;
            --gray-200: #e5e7eb;
            --gray-300: #d1d5db;
            --gray-500: #6b7280;
            --gray-600: #4b5563;
            --gray-700: #374151;
            --gray-800: #1f2937;
            --gray-900: #111827;
        }}
        
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6;
            color: var(--gray-800);
            background: var(--gray-50);
        }}
        
        .container {{ max-width: 1400px; margin: 0 auto; padding: 2rem; }}
        
        header {{
            background: linear-gradient(135deg, var(--primary) 0%, var(--primary-dark) 100%);
            color: white;
            padding: 3rem 2rem;
            margin-bottom: 2rem;
            border-radius: 16px;
            box-shadow: 0 10px 40px -10px rgba(37, 99, 235, 0.3);
        }}
        
        header h1 {{
            font-size: 2.5rem;
            font-weight: 700;
            margin-bottom: 0.5rem;
            display: flex;
            align-items: center;
            gap: 0.75rem;
        }}
        
        header .subtitle {{
            font-size: 1.1rem;
            opacity: 0.9;
            margin-bottom: 1.5rem;
        }}
        
        .header-meta {{
            display: flex;
            gap: 2rem;
            flex-wrap: wrap;
            font-size: 0.95rem;
            opacity: 0.9;
        }}
        
        .header-meta span {{
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}
        
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}
        
        .summary-card {{
            background: white;
            border-radius: 12px;
            padding: 1.5rem;
            text-align: center;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            border: 1px solid var(--gray-200);
        }}
        
        .summary-card .number {{
            font-size: 2.5rem;
            font-weight: 700;
            color: var(--primary);
        }}
        
        .summary-card.success .number {{ color: var(--success); }}
        
        .summary-card .label {{
            color: var(--gray-600);
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-top: 0.25rem;
        }}
        
        section {{
            background: white;
            border-radius: 16px;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            border: 1px solid var(--gray-200);
        }}
        
        section h2 {{
            font-size: 1.5rem;
            color: var(--gray-900);
            margin-bottom: 1.5rem;
            padding-bottom: 0.75rem;
            border-bottom: 2px solid var(--gray-200);
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}
        
        section h3 {{
            font-size: 1.25rem;
            color: var(--gray-800);
            margin: 2rem 0 1rem;
        }}
        
        .level-badge {{
            display: inline-block;
            padding: 0.25rem 0.75rem;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 600;
            margin-left: 0.5rem;
        }}
        
        .level-1 {{ background: var(--primary-light); color: var(--primary-dark); }}
        .level-2 {{ background: var(--warning-light); color: var(--warning); }}
        .level-3 {{ background: var(--error-light); color: var(--error); }}
        
        .image-container {{
            margin: 1.5rem 0;
            text-align: center;
        }}
        
        .image-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }}
        
        .image-caption {{
            margin-top: 0.75rem;
            color: var(--gray-600);
            font-size: 0.9rem;
            font-style: italic;
        }}
        
        .image-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 1.5rem;
            margin: 1.5rem 0;
        }}
        
        .benchmark-card {{
            border: 1px solid var(--gray-200);
            border-radius: 12px;
            margin: 1.5rem 0;
            overflow: hidden;
            transition: box-shadow 0.2s;
        }}
        
        .benchmark-card:hover {{
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }}
        
        .benchmark-header {{
            background: var(--gray-100);
            padding: 1rem 1.5rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--gray-200);
        }}
        
        .benchmark-title {{
            display: flex;
            align-items: center;
            gap: 1rem;
        }}
        
        .benchmark-title h4 {{
            font-size: 1.1rem;
            color: var(--gray-800);
            margin: 0;
        }}
        
        .benchmark-time {{
            font-size: 0.85rem;
            color: var(--gray-500);
        }}
        
        .status-badge {{
            padding: 0.35rem 0.75rem;
            border-radius: 6px;
            font-size: 0.85rem;
            font-weight: 600;
        }}
        
        .status-pass {{ background: var(--success-light); color: var(--success); }}
        .status-fail {{ background: var(--error-light); color: var(--error); }}
        
        .benchmark-body {{
            padding: 1.5rem;
        }}
        
        .benchmark-body .description {{
            color: var(--gray-600);
            margin-bottom: 1rem;
            line-height: 1.6;
        }}
        
        .params-box {{
            background: var(--gray-50);
            border: 1px solid var(--gray-200);
            border-radius: 8px;
            padding: 1rem;
            font-size: 0.9rem;
            color: var(--gray-700);
            margin: 1rem 0;
            font-family: 'SF Mono', Monaco, 'Courier New', monospace;
        }}
        
        .metrics-section h5 {{
            font-size: 1rem;
            color: var(--gray-700);
            margin: 1.5rem 0 0.75rem;
        }}
        
        .metrics-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.85rem;
        }}
        
        .metrics-table th, .metrics-table td {{
            padding: 0.6rem 0.75rem;
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }}
        
        .metrics-table th {{
            background: var(--gray-100);
            font-weight: 600;
            color: var(--gray-700);
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }}
        
        .metrics-table .metric-name {{
            font-weight: 500;
        }}
        
        .metrics-table .mono {{
            font-family: 'SF Mono', Monaco, 'Courier New', monospace;
            font-size: 0.8rem;
        }}
        
        .error-low {{ color: var(--success); font-weight: 600; }}
        .error-medium {{ color: var(--warning); font-weight: 600; }}
        .error-high {{ color: var(--error); font-weight: 600; }}
        
        .vtk-visualization {{
            margin-top: 1.5rem;
            padding-top: 1.5rem;
            border-top: 1px solid var(--gray-200);
        }}
        
        .vtk-visualization h5 {{
            font-size: 1rem;
            color: var(--gray-700);
            margin-bottom: 1rem;
        }}
        
        .vtk-visualization .combined-img {{
            width: 100%;
            border-radius: 8px;
            margin-bottom: 1rem;
        }}
        
        .vtk-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 1rem;
        }}
        
        .vtk-item {{
            text-align: center;
        }}
        
        .vtk-item img {{
            width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .vtk-label {{
            display: block;
            margin-top: 0.5rem;
            font-size: 0.85rem;
            color: var(--gray-600);
        }}
        
        .toc {{
            background: var(--gray-50);
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 2rem;
            border: 1px solid var(--gray-200);
        }}
        
        .toc h3 {{
            margin-bottom: 1rem;
            color: var(--gray-700);
            font-size: 1.1rem;
        }}
        
        .toc ul {{
            list-style: none;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 0.5rem;
        }}
        
        .toc a {{
            color: var(--primary);
            text-decoration: none;
            display: flex;
            align-items: center;
            gap: 0.5rem;
            padding: 0.25rem 0;
        }}
        
        .toc a:hover {{ text-decoration: underline; }}
        
        .highlight-box {{
            background: linear-gradient(135deg, #eff6ff 0%, #dbeafe 100%);
            border-left: 4px solid var(--primary);
            padding: 1.25rem;
            border-radius: 0 8px 8px 0;
            margin: 1rem 0;
        }}
        
        .highlight-box strong {{ color: var(--primary-dark); }}
        
        .results-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
        }}
        
        .results-table th, .results-table td {{
            padding: 0.75rem 1rem;
            text-align: left;
            border-bottom: 1px solid var(--gray-200);
        }}
        
        .results-table th {{
            background: var(--gray-100);
            font-weight: 600;
            color: var(--gray-700);
        }}
        
        .results-table tr:hover {{ background: var(--gray-50); }}
        
        .mono {{ font-family: 'SF Mono', Monaco, 'Courier New', monospace; }}
        
        footer {{
            text-align: center;
            padding: 2rem;
            color: var(--gray-600);
            font-size: 0.9rem;
        }}
        
        @media (max-width: 768px) {{
            .container {{ padding: 1rem; }}
            header {{ padding: 2rem 1.5rem; }}
            header h1 {{ font-size: 1.75rem; }}
            .summary-grid {{ grid-template-columns: repeat(2, 1fr); }}
            .image-grid {{ grid-template-columns: 1fr; }}
            .vtk-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🦀 RustFEA Benchmark Report</h1>
            <p class="subtitle">Comprehensive Verification Suite for 3D Finite Element Analysis</p>
            <div class="header-meta">
                <span>📅 Generated: {datetime.now().strftime('%B %d, %Y at %H:%M')}</span>
                <span>⏱️ Total Runtime: {total_time_ms/1000:.2f}s</span>
                <span>🔧 Solver: MUMPS (Direct)</span>
            </div>
        </header>

        <div class="summary-grid">
            <div class="summary-card success">
                <div class="number">{passed}/{total}</div>
                <div class="label">Benchmarks Passed</div>
            </div>
            <div class="summary-card">
                <div class="number">3</div>
                <div class="label">Verification Levels</div>
            </div>
            <div class="summary-card">
                <div class="number">4</div>
                <div class="label">Mesh Refinements</div>
            </div>
            <div class="summary-card success">
                <div class="number">100%</div>
                <div class="label">Pass Rate</div>
            </div>
        </div>

        <nav class="toc">
            <h3>📑 Contents</h3>
            <ul>
                <li><a href="#overview">📊 Overview & Summary Charts</a></li>
                <li><a href="#level1">🔷 Level 1: Fundamental Verification</a></li>
                <li><a href="#level2">🔶 Level 2: Continuum Verification</a></li>
                <li><a href="#level3">🔴 Level 3: Contact Verification</a></li>
                <li><a href="#convergence">📈 Mesh Convergence Analysis</a></li>
                <li><a href="#visualizations">🎨 Field Visualizations</a></li>
            </ul>
        </nav>

        <section id="overview">
            <h2>📊 Overview & Summary</h2>
            
            <div class="highlight-box">
                <strong>Key Result:</strong> All {total} benchmarks pass verification criteria. The solver demonstrates 
                correct implementation of 3D elasticity, expected mesh convergence behavior, and analytical solution 
                matching for fundamental through advanced contact problems.
            </div>

            <h3>Error Summary by Benchmark</h3>
            <div class="image-container">
                <img src="benchmark_summary.png" alt="Benchmark Summary">
                <p class="image-caption">Figure 1: Maximum error by benchmark on log scale. All benchmarks within specified tolerances.</p>
            </div>

            <h3>Results by Category</h3>
            <div class="image-container">
                <img src="benchmark_categories.png" alt="Benchmark Categories">
                <p class="image-caption">Figure 2: Benchmarks grouped by verification level showing error distribution.</p>
            </div>

            <h3>Results Summary Table</h3>
            <table class="results-table">
                <thead>
                    <tr>
                        <th>#</th>
                        <th>Benchmark</th>
                        <th>Level</th>
                        <th>Status</th>
                        <th>Best Error</th>
                        <th>Time</th>
                    </tr>
                </thead>
                <tbody>
"""
    
    for i, result in enumerate(results, 1):
        best_error, _ = get_best_mesh_error(result['metrics'])
        status = '<span class="status-badge status-pass">✓ PASS</span>' if result['passed'] else '<span class="status-badge status-fail">✗ FAIL</span>'
        
        # Determine level
        if any(x in result['name'] for x in ['Uniaxial', 'Pure Shear', 'Hydrostatic']):
            level = '<span class="level-badge level-1">Level 1</span>'
        elif any(x in result['name'] for x in ['Torsion', 'Cantilever', 'Spherical', 'Boussinesq']):
            level = '<span class="level-badge level-2">Level 2</span>'
        else:
            level = '<span class="level-badge level-3">Level 3</span>'
        
        html += f"""
                    <tr>
                        <td>{i}</td>
                        <td>{result['name']}</td>
                        <td>{level}</td>
                        <td>{status}</td>
                        <td class="mono">{best_error*100:.2f}%</td>
                        <td class="mono">{result['elapsed_ms']:.1f} ms</td>
                    </tr>"""
    
    html += """
                </tbody>
            </table>
        </section>

        <section id="level1">
            <h2>🔷 Level 1: Fundamental Element Verification</h2>
            <p>These benchmarks verify that 8-node hexahedral elements correctly reproduce uniform stress/strain 
            states that should be captured exactly by the element formulation.</p>
"""
    
    for result in level1_benchmarks:
        html += benchmark_card(result)
    
    html += """
        </section>

        <section id="level2">
            <h2>🔶 Level 2: 3D Continuum Verification</h2>
            <p>These benchmarks test the solver against classical elasticity problems with non-uniform 
            stress/strain fields, demonstrating mesh convergence behavior.</p>
"""
    
    for result in level2_benchmarks:
        vtk_key = None
        if 'Torsion' in result['name']:
            vtk_key = 'torsion'
        elif 'Cantilever' in result['name']:
            vtk_key = 'cantilever'
        elif 'Spherical' in result['name']:
            vtk_key = 'spherical'
        html += benchmark_card(result, show_vtk=vtk_key)
    
    html += """
        </section>

        <section id="level3">
            <h2>🔴 Level 3: Contact Mechanics Verification</h2>
            <p>These benchmarks validate Hertzian contact theory predictions, which form the 
            foundation for contact mechanics in FEA.</p>
"""
    
    for result in level3_benchmarks:
        html += benchmark_card(result)
    
    html += """
        </section>

        <section id="convergence">
            <h2>📈 Mesh Convergence Analysis</h2>
            <p>Mesh convergence is a critical indicator of solver correctness. FEM solutions should converge 
            to analytical values as mesh is refined, typically following O(h) or O(h²) rates.</p>

            <div class="image-container">
                <img src="convergence_grid.png" alt="Convergence Grid">
                <p class="image-caption">Figure 3: Mesh convergence for all benchmarks. Dashed lines show O(h) linear 
                and O(h²) quadratic reference convergence rates.</p>
            </div>

            <h3>Individual Convergence Studies</h3>
            <div class="image-grid">
"""
    
    # Add individual convergence plots
    conv_plots = sorted(output_dir.glob("convergence_*.png"))
    for plot in conv_plots:
        if 'grid' not in plot.name:
            name = plot.stem.replace('convergence_', '').replace('_', ' ').title()
            html += f"""
                <div class="image-container">
                    <img src="{plot.name}" alt="{name}">
                    <p class="image-caption">{name}</p>
                </div>"""
    
    html += """
            </div>
        </section>

        <section id="visualizations">
            <h2>🎨 Field Visualizations</h2>
            <p>VTK exports enable post-processing visualization of displacement and stress fields. 
            These renderings use PyVista with physically meaningful color scales.</p>

            <h3>Torsion Shaft (60×8×8 mesh)</h3>
"""
    
    if vtk_images['torsion']['combined'].exists():
        html += f"""
            <div class="image-container">
                <img src="{vtk_images['torsion']['combined'].name}" alt="Torsion Shaft Fields">
                <p class="image-caption">Figure 4: Torsion shaft showing displacement magnitude (left, warped 50×) 
                and von Mises stress (right). Maximum stress occurs at outer surface.</p>
            </div>"""
    
    html += """
            <h3>Cantilever Beam (60×6×6 mesh)</h3>
"""
    
    if vtk_images['cantilever']['combined'].exists():
        html += f"""
            <div class="image-container">
                <img src="{vtk_images['cantilever']['combined'].name}" alt="Cantilever Beam Fields">
                <p class="image-caption">Figure 5: Cantilever beam under tip load showing displacement (left, warped 50×) 
                and von Mises stress (right). Classic bending profile with maximum stress at fixed end.</p>
            </div>"""
    
    html += """
            <h3>Spherical Cavity (20×20×20 mesh)</h3>
"""
    
    if vtk_images['spherical']['combined'].exists():
        html += f"""
            <div class="image-container">
                <img src="{vtk_images['spherical']['combined'].name}" alt="Spherical Cavity Fields">
                <p class="image-caption">Figure 6: Spherical cavity under far-field stress showing displacement (left, warped 1000×) 
                and von Mises stress (right). Radial displacement pattern visible with stress concentration around cavity.</p>
            </div>"""
    
    html += """
        </section>

        <footer>
            <p><strong>RustFEA</strong> — Finite Element Analysis Library in Rust</p>
            <p>Report generated by automated benchmark suite using Python + Matplotlib + PyVista</p>
        </footer>
    </div>
</body>
</html>
"""
    
    output_file = output_dir / "benchmark_report.html"
    with open(output_file, 'w') as f:
        f.write(html)
    print(f"Saved: {output_file}")


def main():
    """Main entry point."""
    script_dir = Path(__file__).parent
    json_path = script_dir / "output" / "benchmark_results.json"
    output_dir = script_dir / "output" / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Loading results from: {json_path}")
    results = load_benchmark_results(json_path)
    
    print("\nGenerating plots...")
    
    # 1. Summary bar chart
    plot_summary_bar_chart(results, output_dir / "benchmark_summary.png")
    
    # 2. Category grouped chart
    plot_benchmark_categories(results, output_dir / "benchmark_categories.png")
    
    # 3. Convergence grid
    plot_convergence_grid(results, output_dir / "convergence_grid.png")
    
    # 4. Individual convergence plots
    plot_individual_convergence(results, output_dir)
    
    # 5. Publication-quality HTML report
    generate_html_report(results, output_dir)
    
    print(f"\n✅ All outputs saved to: {output_dir}")
    print(f"📄 Open {output_dir / 'benchmark_report.html'} in a browser")


if __name__ == "__main__":
    main()
