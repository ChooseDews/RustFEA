use rust_fea::io::project::Project;
use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;


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
    let mut group = c.benchmark_group("Direct Solve");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(40));
    group.warm_up_time(Duration::from_secs(1));
    group.bench_function("Load Project: Tube", |b: &mut criterion::Bencher<'_>| b.iter(|| load_project()));
    group.bench_function("Load + Direct Solve: Tube", |b| b.iter(|| project()));
    group.finish();
}


criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);