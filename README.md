# RustFEA

A Finite Element Analysis (FEA) library and GUI written in Rust. Provides efficient, robust tools for solid mechanics simulations with both direct and explicit solving capabilities.

**[Try the Web Demo](https://choosedews.github.io/RustFEA/)** | [Documentation](https://github.com/ChooseDews/RustFEA)

## Features

- **Pure Rust** - No external C dependencies; builds anywhere Rust builds
- **Desktop & Web** - Native app plus WebAssembly version ([try it](https://choosedews.github.io/RustFEA/))
- **Modern Solvers** - `faer` (Cholesky/LU) for direct solving, explicit time integration for dynamics/contact
- **Multiple Element Types** - C3D8 (8-node brick), C3D20 (20-node quadratic brick), C3D4 (4-node tetrahedra)
- **Contact Mechanics** - Penalty method contact via explicit solver
- **Section Cuts** - Interactive clipping planes with 2D cross-section view and PNG export
- **VTK Export** - Standard VTK and pure-Rust VTK HDF output for ParaView

## Quick Start

```bash
# Build the library
cargo build --release

# Run an example simulation
cargo run --bin read_input -- examples/tube_benchmark.toml -v

# Run the GUI (native desktop)
cargo run -p rust_fea_gui --release
```

See [`gui/README.md`](gui/README.md) for GUI-specific documentation and WASM build instructions.

## Project Structure

```
├── src/              # Core FEA library
│   ├── elements/     # Element formulations (C3D8, C3D20, C3D4)
│   ├── bc/           # Boundary conditions (fixed, load, contact, pressure)
│   ├── simulation/   # Simulation orchestration
│   ├── solver.rs     # Direct sparse solver (faer)
│   └── io/           # File I/O (TOML, VTK, mesh formats)
├── gui/              # egui-based graphical interface
├── examples/         # Example simulations and mesh sources
├── test/             # Integration tests
├── benches/          # Criterion benchmarks
└── reports/          # Validation reports and comparisons
```

## Mesh Support

Import meshes in these formats:
- `.inp` - Abaqus format (recommended via [Gmsh](http://gmsh.info/))
- `.bin.xz`, `.json.xz`, `.bin`, `.json` - RustFEA internal formats

Gmsh geometry files are in `examples/mesh_src/*.geo`.

## Elements

| Element | Nodes | Order | Description |
|---------|-------|-------|-------------|
| C3D8 | 8 | Linear | Hexahedral brick |
| C3D20 | 20 | Quadratic | Hexahedral brick with mid-edge nodes |
| C3D4 | 4 | Linear | Tetrahedral |

## Boundary Conditions

- **Dirichlet** - Fixed displacement (clamped faces/nodes)
- **Neumann** - Applied loads, pressure, torque, body forces
- **Contact** - Node-to-segment penalty contact (explicit solver only)

## Solvers

### Direct Solver
Uses `faer` for sparse Cholesky decomposition (symmetric positive-definite systems) with LU fallback. Pure Rust with SIMD optimization, ~1.5-2x faster than UMFPACK on typical FEA matrices.

### Explicit Solver
Central-difference time integration for dynamics and contact problems. Element-by-element assembly avoids global matrix formation.

## Output Formats

- **VTK** (`.vtk`) - Legacy VTK format for ParaView
- **VTK HDF** (`.hdf`) - Modern HDF5-based format (pure Rust via `hdf5-writer`)

## Benchmarks

Run the validation benchmark suite:

```bash
cargo run --bin run_benchmarks --release
```

Results compare against analytical solutions for cantilever beams, torsion shafts, hollow spheres, and Hertzian contact.

## Building

### Native
```bash
cargo build --release
```

### WebAssembly
```bash
cd gui
./build_wasm.sh
```

See [`gui/README.md`](gui/README.md) for detailed WASM instructions.

## License

[Apache License 2.0](LICENSE)

## Author

A project by [John Dews-Flick](https://johndews.com)
