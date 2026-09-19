#!/usr/bin/env python3
"""
VTK Visualization Script for FEA Benchmark Results
Renders displacement and von Mises stress fields from VTK files.
"""

import pyvista as pv
from pathlib import Path
import numpy as np

# Configure PyVista for off-screen rendering
pv.OFF_SCREEN = True
pv.global_theme.background = 'white'
pv.global_theme.font.color = 'black'
pv.global_theme.font.size = 12


def render_vtk_fields(vtk_path, output_dir, name):
    """Render displacement and stress fields from a VTK file."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load mesh
    mesh = pv.read(vtk_path)
    print(f"Loaded: {vtk_path}")
    print(f"  Points: {mesh.n_points}, Cells: {mesh.n_cells}")
    print(f"  Arrays: {mesh.array_names}")
    
    # Get displacement magnitude
    if 'displacement_magnitude' in mesh.array_names:
        disp_mag = mesh['displacement_magnitude']
    elif 'displacement' in mesh.array_names:
        disp = mesh['displacement']
        disp_mag = np.linalg.norm(disp, axis=1)
    else:
        print(f"  Warning: No displacement data found")
        disp_mag = None
    
    # Get von Mises stress
    if 'von_mises_stress' in mesh.array_names:
        vm_stress = mesh['von_mises_stress']
    elif 'vm' in mesh.array_names:
        vm_stress = mesh['vm']
    else:
        vm_stress = None
    
    # Scale factor for deformed shape (auto-scale)
    if disp_mag is not None and disp_mag.max() > 0:
        scale = mesh.length / disp_mag.max() * 0.1
    else:
        scale = 1.0
    
    # Create deformed mesh if displacement vectors exist
    if 'displacement' in mesh.array_names:
        warped = mesh.warp_by_vector('displacement', factor=scale)
    else:
        warped = mesh
    
    # 1. Render displacement magnitude
    if disp_mag is not None:
        p = pv.Plotter(off_screen=True, window_size=[1200, 900])
        p.add_mesh(
            warped if 'displacement' in mesh.array_names else mesh,
            scalars=disp_mag,
            cmap='turbo',
            scalar_bar_args={
                'title': 'Displacement (m)',
                'title_font_size': 14,
                'label_font_size': 12,
                'shadow': True,
                'n_labels': 5,
                'fmt': '%.2e',
            },
            show_edges=True,
            edge_color='gray',
            edge_opacity=0.3,
        )
        p.add_text(f"{name}\nDisplacement Field (scale: {scale:.0f}x)", 
                  font_size=12, position='upper_left')
        p.add_axes()
        p.view_isometric()
        p.camera.zoom(1.2)
        
        output_file = output_dir / f"{name}_displacement.png"
        p.screenshot(str(output_file))
        p.close()
        print(f"  Saved: {output_file}")
    
    # 2. Render von Mises stress
    if vm_stress is not None and vm_stress.max() > 0:
        p = pv.Plotter(off_screen=True, window_size=[1200, 900])
        p.add_mesh(
            warped if 'displacement' in mesh.array_names else mesh,
            scalars=vm_stress,
            cmap='jet',
            scalar_bar_args={
                'title': 'von Mises Stress (Pa)',
                'title_font_size': 14,
                'label_font_size': 12,
                'shadow': True,
                'n_labels': 5,
                'fmt': '%.2e',
            },
            show_edges=True,
            edge_color='gray',
            edge_opacity=0.3,
        )
        p.add_text(f"{name}\nvon Mises Stress Field", 
                  font_size=12, position='upper_left')
        p.add_axes()
        p.view_isometric()
        p.camera.zoom(1.2)
        
        output_file = output_dir / f"{name}_vonmises.png"
        p.screenshot(str(output_file))
        p.close()
        print(f"  Saved: {output_file}")
    
    # 3. Multi-view render (for publication)
    if disp_mag is not None:
        p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=[1600, 700])
        
        # Left: Displacement
        p.subplot(0, 0)
        p.add_mesh(
            warped if 'displacement' in mesh.array_names else mesh,
            scalars=disp_mag,
            cmap='turbo',
            scalar_bar_args={'title': 'Displacement (m)', 'fmt': '%.2e'},
            show_edges=True,
            edge_color='gray',
            edge_opacity=0.2,
        )
        p.add_text("Displacement", font_size=11)
        p.view_isometric()
        
        # Right: von Mises stress
        p.subplot(0, 1)
        if vm_stress is not None and vm_stress.max() > 0:
            p.add_mesh(
                warped if 'displacement' in mesh.array_names else mesh,
                scalars=vm_stress,
                cmap='jet',
                scalar_bar_args={'title': 'von Mises (Pa)', 'fmt': '%.2e'},
                show_edges=True,
                edge_color='gray',
                edge_opacity=0.2,
            )
            p.add_text("von Mises Stress", font_size=11)
        else:
            p.add_mesh(mesh, color='lightblue', show_edges=True)
            p.add_text("Stress N/A", font_size=11)
        p.view_isometric()
        
        output_file = output_dir / f"{name}_combined.png"
        p.screenshot(str(output_file))
        p.close()
        print(f"  Saved: {output_file}")
    
    return True


def main():
    """Process all VTK files in the output directory."""
    script_dir = Path(__file__).parent
    vtk_dir = script_dir / "output" / "vtk"
    output_dir = script_dir / "output" / "plots"
    
    if not vtk_dir.exists():
        print(f"VTK directory not found: {vtk_dir}")
        print("Run benchmarks first with: cargo run --bin run_benchmarks")
        return
    
    # Find all VTK files
    vtk_files = list(vtk_dir.glob("*.vtk"))
    
    if not vtk_files:
        print(f"No VTK files found in: {vtk_dir}")
        return
    
    print(f"\nProcessing {len(vtk_files)} VTK files...")
    print("=" * 50)
    
    for vtk_path in sorted(vtk_files):
        name = vtk_path.stem
        print(f"\n{name}:")
        try:
            render_vtk_fields(vtk_path, output_dir, name)
        except Exception as e:
            print(f"  Error: {e}")
    
    print("\n" + "=" * 50)
    print(f"Visualization complete. Images saved to: {output_dir}")


if __name__ == "__main__":
    main()
