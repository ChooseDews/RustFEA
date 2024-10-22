# Rust FEA

RustFEA is a Finite Element Analysis (FEA) library written in Rust. This library aims to provide efficient and robust tools for performing finite element analysis, leveraging the safety and performance benefits of Rust.

#### Features
- **Modular Design for Extension**: The library is divided into several modules such as `node`, `elements`, `simulation`, `io`, `utilities`, `mesh`, `solver`, and `bc`.
- **Dependency Integration**: Utilizes popular Rust crates like `nalgebra`, `russell_sparse`, and `serde` for mathematical operations and data serialization.
- **Input Formats**: Supports reading and writing of input files in the .toml format.

#### Limitations
- 3D only
- Solid mechanics focused only


#### Building

To get suitesparse and mumps working on Ubuntu you need to install the following dependencies:

```bash
sudo apt-get install -y --no-install-recommends \
    liblapacke-dev \
    libopenblas-dev \
    libsuitesparse-dev
```

To build the project you need to have `rust` and `cargo` installed then run:

```bash
cargo build
```

#### Usage

Take a look in the `examples/*.toml` folder for some example simulation input files.

Here is a basic example to get you started:

```bash
cargo run --bin read_input -- examples/tube_benchmark.toml -v
```


### Rust FEA Library Overview

Rust FEA is a Finite Element Analysis (FEA) library written in Rust, designed to provide efficient and robust tools for performing finite element analysis.

#### Features
- **Modular Design**: The library includes various modules such as `node`, `elements`, `simulation`, `io`, `utilities`, `mesh`, `solver`, and `bc`. Each module is responsible for different aspects of the FEA process, ensuring a clean and organized codebase.
- **Dependency Integration**: Integration with popular Rust crates like `nalgebra` for mathematical operations, `russell_sparse` for sparse matrix computations, and `serde` for data serialization. This ensures the library is both powerful and flexible.
- **Customization**: The library supports a wide range of configurations and customizations, allowing users to tailor the simulation and solver settings to their specific needs.

#### Mesh Support
This library supports the .inp format for importing meshes. Mesh are represented in a struct called `MeshAssembly` which can also be added together to create a larger mesh or multi body meshes.

- `.inp` format
- Internal Format: `.bin.xz`, `.json.xz`, `.bin`, `.json`

[Gmsh](http://gmsh.info/) is recommended for generating meshes then exporting via `.inp` (selecting export options while saving). You can find many gmsh examples under `examples/mesh_src/*.geo` 

#### Elements
The library includes various finite elements, each represented as a structure with associated methods.
- 8-Node Isoparametric Solid Element (Brick) (3D)
- 4-Node Isoparametric Surface Element (Quad) (2D)

#### Boundary Conditions
Boundary conditions used to apply forces or constraints to the problem.

- Fixed Value (Dirichlet)
    - Prescribed Displacement
- Prescribed Force (Neumann)
    - Load
    - Normal Contact: Penalty method

#### `Simulation`
This struct encapsulates the entire simulation setup, including the mesh, elements, boundary conditions, and solver settings. It provides a high-level interface for configuring and executing the simulation.

#### `Project` File
This struct represents a project. Which can contain mulitple simulations. Any project `.toml` is loaded in as a `Project` struct.

#### Solver

##### Direct Solving
`russell_sparse` and `nalgebra` are used for sparse matrix operations and direct solvers. `russell_sparse` uses UMFPACK and MUMPS to solve large sparse linear systems of equations. 

##### Explicit Solving
`nalgebra` is used for matrix operations. Overall there is a focus on element based assembly. Contact is only supported via explicit solving.

## Todo
- [~] add tests for common elements and assemblies
- [x] add nodal property outputs
- [ ] Use better solver for Ax=b. Upgrade to newest russell library
- [ ] add intergrations with paraview?
- [x] add an explict solving scheme
- [x] contact
    - [x] Multi-body handling -> Multi-body mesh -> some way to merge meshes together
    - [ ] Fix issues with contact on element edges
- [x] Switch config file to use toml format
- [ ] New Boundary Conditions
    - [ ] Pressure / Normal vector based BC
    - [ ] Shear / Tangential vector based BC
- [ ] New Element Types
    - [ ] 3D 4 node tetrahedra
    - [ ] Shell element
- [ ] VTL HDF Export
    - [x] Single time step export
    - [ ] Multi time step export
- [ ] Generalized field handling
- [ ] Body forces
- [ ] Simple material input handling
- [ ] Non-linear material handling
### Validations
- [ ] Beam problems compared to abaqus and exact solutions
### Problems
- [ ] What about nodes that are not part of any element?