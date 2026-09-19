use rust_fea::bc::{FixedCondition, LoadCondition};
use rust_fea::elements::BrickElement;
use rust_fea::io::mesh_reader;
use rust_fea::io::vtk_writer;
use rust_fea::simulation::Simulation;
use rust_fea::node::Node;
use rust_fea::elements::Material;
use std::time::Instant;

const TEST_OUTPUT_DIR: &str = "examples/output";

fn ensure_output_dir() {
    std::fs::create_dir_all(TEST_OUTPUT_DIR).ok();
}

/// Generate a simple bar mesh: n_x elements along X, n_y along Y, n_z along Z
/// Returns nodes and elements for a bar of size (s_x, s_y, s_z)
fn generate_bar_mesh(
    n_x: usize,
    n_y: usize,
    n_z: usize,
    s_x: f64,
    s_y: f64,
    s_z: f64,
) -> (Vec<Node>, Vec<BrickElement>) {
    let mut nodes = Vec::new();
    let mut elements = Vec::new();
    let dx = s_x / n_x as f64;
    let dy = s_y / n_y as f64;
    let dz = s_z / n_z as f64;

    // Generate nodes
    let mut id = 0;
    for k in 0..=n_z {
        for j in 0..=n_y {
            for i in 0..=n_x {
                nodes.push(Node::new(id, i as f64 * dx, j as f64 * dy, k as f64 * dz));
                id += 1;
            }
        }
    }

    // Generate element connectivity (8-node hex elements)
    let mut elem_id = 0;
    for k in 0..n_z {
        for j in 0..n_y {
            for i in 0..n_x {
                let n0 = i + j * (n_x + 1) + k * (n_x + 1) * (n_y + 1);
                let n1 = n0 + 1;
                let n2 = n0 + (n_x + 1) + 1;
                let n3 = n0 + (n_x + 1);
                let n4 = n0 + (n_x + 1) * (n_y + 1);
                let n5 = n4 + 1;
                let n6 = n4 + (n_x + 1) + 1;
                let n7 = n4 + (n_x + 1);

                elements.push(BrickElement::new(
                    elem_id,
                    vec![n0, n1, n2, n3, n4, n5, n6, n7],
                    Material::aluminum(),
                ));
                elem_id += 1;
            }
        }
    }

    (nodes, elements)
}

/// Get node IDs on a specific X-plane (for applying boundary conditions)
fn get_nodes_on_x_plane(nodes: &[Node], x_value: f64, tolerance: f64) -> Vec<usize> {
    nodes
        .iter()
        .filter(|n| (n.position.x - x_value).abs() < tolerance)
        .map(|n| n.id)
        .collect()
}

#[test]
fn write_vtkhdf() {
    let start = Instant::now();
    ensure_output_dir();

    let output_file = format!("{}/test_vtkhdf_out.hdf", TEST_OUTPUT_DIR);

    // Clean up from previous runs
    if std::path::Path::new(&output_file).exists() {
        std::fs::remove_file(&output_file).unwrap();
    }

    // Load the example mesh and create a simulation
    let mesh = mesh_reader::read_file("examples/meshes/example_tube.inp.xz");
    let sim = Simulation::from_mesh(mesh, 3);

    // Write VTKHDF output
    vtk_writer::write_vtkhdf(&output_file, &sim).expect("Failed to write VTKHDF file");
    assert!(
        std::path::Path::new(&output_file).exists(),
        "Expected VTKHDF file to exist at {}",
        output_file
    );

    let duration = start.elapsed();
    println!("write_vtkhdf test took: {:?}", duration);
}

#[test]
fn write_vtk_from_mesh() {
    let start = Instant::now();
    ensure_output_dir();

    let output_file = format!("{}/test_vtk_hdf_out.vtk", TEST_OUTPUT_DIR);

    // Clean up from previous runs
    if std::path::Path::new(&output_file).exists() {
        std::fs::remove_file(&output_file).unwrap();
    }

    // Load the example mesh and create a simulation
    let mesh = mesh_reader::read_file("examples/meshes/example_tube.inp.xz");
    let sim = Simulation::from_mesh(mesh, 3);

    // Write VTK output
    vtk_writer::write_vtk(&output_file, &sim).expect("Failed to write VTK file");
    assert!(
        std::path::Path::new(&output_file).exists(),
        "Expected VTK file to exist at {}",
        output_file
    );

    let duration = start.elapsed();
    println!("write_vtk_from_mesh test took: {:?}", duration);
}

/// Test that runs a simple FEA simulation (direct solver) on a small bar mesh
/// and writes the results to VTK and VTKHDF formats.
#[test]
fn run_simulation_and_write_output() {
    let start = Instant::now();
    ensure_output_dir();

    let vtk_output = format!("{}/test_sim_bar.vtk", TEST_OUTPUT_DIR);
    let hdf_output = format!("{}/test_sim_bar.hdf", TEST_OUTPUT_DIR);

    // Clean up from previous runs
    for path in [&vtk_output, &hdf_output] {
        if std::path::Path::new(path).exists() {
            std::fs::remove_file(path).unwrap();
        }
    }

    // Create a small bar mesh: 5x2x2 elements = 20 elements total
    // Bar dimensions: 10mm x 2mm x 2mm (long thin bar along X)
    let (nodes, elements) = generate_bar_mesh(5, 2, 2, 10.0, 2.0, 2.0);

    println!(
        "Created mesh with {} nodes and {} elements",
        nodes.len(),
        elements.len()
    );

    // Find nodes at x=0 (fixed end) and x=10 (loaded end)
    let fixed_nodes = get_nodes_on_x_plane(&nodes, 0.0, 0.001);
    let loaded_nodes = get_nodes_on_x_plane(&nodes, 10.0, 0.001);

    println!(
        "Fixed nodes (x=0): {}, Loaded nodes (x=10): {}",
        fixed_nodes.len(),
        loaded_nodes.len()
    );

    // Create simulation from nodes and elements
    let boxed_elements: Vec<Box<dyn rust_fea::elements::BaseElement>> =
        elements.into_iter().map(|e| Box::new(e) as Box<dyn rust_fea::elements::BaseElement>).collect();
    let mut sim = Simulation::from_arrays(nodes, boxed_elements, 3);

    // Apply boundary conditions:
    // - Fix all DOFs at x=0
    // - Apply tensile load in +X direction at x=10
    let fixed_bc = FixedCondition::static_3d(fixed_nodes);
    sim.add_boundary_condition(Box::new(fixed_bc));

    // Apply 1000N total force in X direction, distributed across loaded nodes
    let load_bc = LoadCondition::new_from_vec(loaded_nodes, vec![1000.0, 0.0, 0.0]);
    sim.add_boundary_condition(Box::new(load_bc));

    // Run the simulation (direct solver)
    println!("Running simulation...");
    sim.solve();
    println!("Simulation complete!");

    // Write VTK output
    vtk_writer::write_vtk(&vtk_output, &sim).expect("Failed to write VTK file");
    assert!(
        std::path::Path::new(&vtk_output).exists(),
        "Expected VTK file to exist"
    );

    // Write VTKHDF output
    vtk_writer::write_vtkhdf(&hdf_output, &sim).expect("Failed to write VTKHDF file");
    assert!(
        std::path::Path::new(&hdf_output).exists(),
        "Expected VTKHDF file to exist"
    );

    // Basic sanity check: max displacement should be at the loaded end (x=10)
    // and should be positive in X direction for tensile load
    let max_disp_x = sim
        .nodes()
        .iter()
        .map(|n| n.displacement.x)
        .fold(f64::NEG_INFINITY, f64::max);

    println!("Max X displacement: {:.6e} mm", max_disp_x);
    assert!(
        max_disp_x > 0.0,
        "Expected positive X displacement for tensile load"
    );

    let duration = start.elapsed();
    println!("run_simulation_and_write_output test took: {:?}", duration);
}
