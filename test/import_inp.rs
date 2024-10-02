use rust_fea::io::mesh_reader;
use rust_fea::simulation::Simulation;
use rust_fea::io::project::Project;
use rust_fea::io::vtk_writer;
use rust_fea::utilities::Keywords;

#[test]
fn import_inp() {
    let mesh = mesh_reader::read_file("examples/meshes/example_tube.inp.xz"); //Expeected: Mesh: None [Nodes: 7790, Elements: 9656, Element Groups: 15, Node Groups: 4]
    assert_eq!(mesh.nodes.len(), 7790, "Expected 7790 nodes, but got {}", mesh.nodes.len());
    assert_eq!(mesh.elements.len(), 9656, "Expected 9656 elements, but got {}", mesh.elements.len());
    assert_eq!(mesh.node_groups.len(), 4, "Expected 4 node groups, but got {}", mesh.node_groups.len());
    println!("{:}", mesh);
}

#[test]
fn import_compressed_inp_save() {
    let mesh = mesh_reader::read_file("examples/meshes/example_tube.inp.xz");
    let sim = Simulation::from_mesh(mesh, 3);
    let project = Project::new(vec![sim], Keywords::new());
    let output_file = "examples/output/test_tube_out.bin.xz";
    project.save_to_file(output_file);
    assert!(std::path::Path::new(output_file).exists(), "Expected file to exist");
}

#[test]
fn import_save_and_load() {
    let mesh = mesh_reader::read_file("examples/meshes/example_tube.inp.xz");
    let sim = Simulation::from_mesh(mesh, 3);
    let mut project = Project::new(vec![sim], Keywords::new());
    let output_file = "examples/output/test_tube_out_2.bin.xz";
    project.save_to_file(output_file);
    assert!(std::path::Path::new(output_file).exists(), "Expected file to exist");
    let mut project_loaded = Project::load(output_file);
    assert_eq!(
        project.get_simulation().get_node(0).unwrap().position,
        project_loaded.get_simulation().get_node(0).unwrap().position,
        "Expected nodes to be the same"
    );
}

#[test]
fn write_vtk() {
    let mesh = mesh_reader::read_file("examples/meshes/example_tube.inp.xz");
    let sim = Simulation::from_mesh(mesh, 3);
    vtk_writer::write_vtk("examples/output/test_tube_out.vtk", &sim);
    assert!(std::path::Path::new("examples/output/test_tube_out.vtk").exists(), "Expected file to exist");
}

