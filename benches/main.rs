use rust_fea::io::project::Project;
use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;


fn project_contact() {
    let mut project = Project::from_input_file("examples/bar_contact_benchmark.toml");
    project.solve();
    project.export_vtk();
    project.save();
    project.print();
}

fn project() {
    let mut project = Project::from_input_file("examples/tube_benchmark.toml");
    project.solve();
    project.export_vtk();
    project.save();
    project.print();
}


fn load_project() {
    let mut project = Project::from_input_file("examples/tube_benchmark.toml");
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("FEA Benchmarks");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(70));
    group.warm_up_time(Duration::from_secs(1));
    group.bench_function("Load Project: Tube", |b: &mut criterion::Bencher<'_>| b.iter(|| load_project()));
    group.bench_function("Load + Direct Solve: Tube", |b| b.iter(|| project()));
    group.bench_function("Explicit Contact Problem", |b| b.iter(|| project_contact()));
    group.finish();
}


criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);