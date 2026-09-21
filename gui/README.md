# RustFEA GUI

A graphical user interface for RustFEA - Finite Element Analysis in Rust.

## Features

- **Mesh Import**: Load mesh files in `.inp` (Abaqus), `.bin`, and `.json` formats
- **3D Visualization**: Interactive 3D view with rotation, pan, and zoom
- **Simulation Setup**: 
  - Material definition with presets (Steel, Aluminum, Titanium)
  - Boundary conditions (Fixed, Load, Torque, Contact)
  - Solver configuration (Direct/Explicit)
- **Simulation Execution**: Run simulations with progress tracking
- **Results Visualization**: 
  - Displacement color maps
  - Deformation scaling
  - Export to VTK and CSV

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
- Optionally install wasm-opt for smaller binaries: `brew install binaryen` or download from https://github.com/WebAssembly/binaryen

```bash
cd gui

# Option 1: Use the build script
./build_wasm.sh

# Option 2: Manual build
cargo build --target wasm32-unknown-unknown --features web --no-default-features --release
wasm-bindgen --target web --out-dir web/dist ../target/wasm32-unknown-unknown/release/rust_fea_gui.wasm

# Serve locally
cd web/dist
python3 -m http.server 8080
```

Then open `http://localhost:8080` in your browser.

**WASM Notes:**
- Direct solver uses nalgebra-sparse Cholesky (pure Rust) instead of UMFPACK
- File import is not available in web version - use primitive mesh generation
- Simulations run synchronously (blocking UI briefly)
- Compressed file formats (.xz, .zst) are not supported

## Usage

### Workflow

1. **Mesh Tab**: Import a mesh file (.inp from Gmsh, or RustFEA's binary/JSON formats)
2. **Setup Tab**: 
   - Configure materials (Young's modulus, Poisson's ratio, density)
   - Add boundary conditions to node groups
   - Select solver type
3. **Run Tab**: Execute the simulation
4. **Results Tab**: Visualize and export results

### Camera Controls

- **Left Mouse Button + Drag**: Rotate view
- **Right Mouse Button + Drag**: Pan view  
- **Scroll Wheel**: Zoom in/out

### Keyboard Shortcuts

- Standard egui shortcuts for copy/paste in text fields

## Architecture

```
gui/
├── src/
│   ├── main.rs          # Entry point (native + web)
│   ├── lib.rs           # Library exports
│   ├── app.rs           # Main application logic
│   ├── state.rs         # Application state management
│   ├── renderer.rs      # GPU rendering (placeholder for wgpu)
│   └── ui/
│       ├── mod.rs       # UI module exports
│       ├── menu_bar.rs  # Top menu bar
│       ├── side_panel.rs # Left workflow panel
│       ├── status_bar.rs # Bottom status bar
│       ├── viewport.rs  # 3D viewport
│       ├── mesh_panel.rs # Mesh import/info
│       ├── setup_panel.rs # BC/material config
│       ├── run_panel.rs # Simulation execution
│       └── results_panel.rs # Results visualization
└── web/
    └── index.html       # Web entry point
```

## Dependencies

- **egui/eframe**: Immediate mode GUI framework
- **wgpu**: Cross-platform GPU abstraction (for future GPU rendering)
- **rfd**: Native file dialogs
- **rust_fea**: The core FEA library

## Current Limitations

- 3D rendering uses egui's software renderer (painter's algorithm)
- GPU-accelerated rendering with wgpu is stubbed for future implementation
- Web file import requires additional implementation
- Explicit solver GUI integration is incomplete

## Contributing

The GUI is designed to be modular. Key extension points:

1. **New visualization modes**: Add to `ColorMode` enum in `state.rs`
2. **New boundary conditions**: Add config structs in `state.rs`, UI in `setup_panel.rs`
3. **GPU rendering**: Implement `MeshRenderer` in `renderer.rs` with wgpu pipeline
