# RustFEA GUI

A graphical user interface for RustFEA - Finite Element Analysis in Rust.

**[Try the Web Demo](https://choosedews.github.io/RustFEA/)**

## Features

- **Mesh Import**: Load mesh files in `.inp` (Abaqus), `.bin`, and `.json` formats
- **3D Visualization**: Interactive 3D view with rotation, pan, zoom
- **Section Cuts**: Clipping planes with 2D cross-section view and PNG export
- **Built-in Examples**: Quick-load example models (cantilever beam, torque shaft, contact blocks)
- **Simulation Setup**: 
  - Material definition with presets (Steel, Aluminum, Titanium)
  - Boundary conditions (Fixed, Load, Torque, Contact, Pressure)
  - Solver configuration (Direct/Explicit)
- **Simulation Execution**: Run simulations with progress tracking
- **Results Visualization**: 
  - Displacement/stress/strain color maps
  - Deformation scaling
  - Export to VTK

## Building

### Native (Desktop)

```bash
# From the project root
cargo build -p rust_fea_gui --release

# Run
cargo run -p rust_fea_gui --release
```

### Web (WASM)

Prerequisites:
- Install wasm-bindgen-cli: `cargo install wasm-bindgen-cli`
- Optionally install wasm-opt for smaller binaries: `brew install binaryen`

```bash
cd gui

# Use the build script
./build_wasm.sh

# Serve locally
cd web
python3 -m http.server 8080
```

Then open `http://localhost:8080` in your browser.

**WASM Notes:**
- Direct solver uses nalgebra-sparse Cholesky (pure Rust) instead of faer
- Compressed file formats (.xz, .zst) are not supported in WASM
- Simulations run synchronously (may briefly freeze the UI on large models)

## Usage

### Workflow

1. **Mesh Tab**: Import a mesh file or load a built-in example
2. **Setup Tab**: Configure materials and boundary conditions
3. **Run Tab**: Execute the simulation
4. **Results Tab**: Visualize and export results

### Camera Controls

| Action | Control |
|--------|---------|
| Rotate | Left mouse + drag |
| Pan | Right mouse + drag |
| Zoom | Scroll wheel |

## Architecture

```
gui/
├── src/
│   ├── main.rs          # Entry point (native + web)
│   ├── lib.rs           # Library exports
│   ├── app.rs           # Main application logic
│   ├── state.rs         # Application state management  
│   ├── icons.rs         # RemixIcon Unicode constants
│   ├── renderer.rs      # Software 3D renderer
│   ├── render_cache.rs  # Mesh render caching
│   ├── section_cut.rs   # Section plane visualization
│   ├── examples.rs      # Built-in example models
│   ├── project_io.rs    # Project file handling
│   ├── web_file_io.rs   # WASM file operations
│   └── ui/
│       ├── menu_bar.rs  # Top menu bar
│       ├── side_panel.rs # Left workflow panel
│       ├── status_bar.rs # Bottom status bar
│       ├── viewport.rs  # 3D viewport
│       ├── mesh_panel.rs # Mesh import/info
│       ├── setup_panel.rs # BC/material config
│       ├── run_panel.rs # Simulation execution
│       └── results_panel.rs # Results visualization
├── assets/              # Fonts (Inter, JetBrains Mono, RemixIcon, Noto Symbols)
└── web/
    └── index.html       # Web entry point
```

## Dependencies

- **egui/eframe**: Immediate mode GUI framework
- **rust_fea**: Core FEA library
- **rfd**: Native file dialogs (desktop only)

## License

[Apache License 2.0](../LICENSE)
