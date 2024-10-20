use rust_fea::io::{project, vtk_writer};
use rust_fea::io::project::Project;
use std::time::Instant;



#[test]
fn write_vtkhdf() {
    let start = Instant::now(); // Start timing

    //delete output file if it exists
    if std::path::Path::new("temp/out.hdf").exists() {
        std::fs::remove_file("temp/out.hdf").unwrap();
    }

    let mut project = Project::load("/home/john/Projects/RustFEA/temp/hertz_contact.json.xz");
    let sim = project.get_simulation();
    vtk_writer::write_vtkhdf("temp/out.hdf", &sim);
    assert!(std::path::Path::new("temp/out.hdf").exists(), "Expected file to exist");

    let duration = start.elapsed(); // Calculate duration
    println!("write_vtkhdf test took: {:?}", duration); // Report duration
}

#[test]
fn write_vtk() {
    let start = Instant::now(); // Start timing

    //delete output file if it exists
    if std::path::Path::new("temp/out.vtk").exists() {
        std::fs::remove_file("temp/out.vtk").unwrap();
    }

    let mut project = Project::load("/home/john/Projects/RustFEA/temp/hertz_contact.json.xz");
    let sim = project.get_simulation();
    vtk_writer::write_vtk("temp/out.vtk", &sim);
    assert!(std::path::Path::new("temp/out.vtk").exists(), "Expected file to exist");

    let duration = start.elapsed(); // Calculate duration
    println!("write_vtk test took: {:?}", duration); // Report duration
}
