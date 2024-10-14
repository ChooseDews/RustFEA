use rust_fea::io::project::Project;

use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;


fn project() {
    let mut project = Project::from_input_file("examples/bar_contact_benchmark.toml");
    project.solve();
    project.export_vtk();
    project.save();
    project.print();
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Example Contact Problem");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(120));
    group.warm_up_time(Duration::from_secs(1));
    group.bench_function("Explicit Contact Problem", |b| b.iter(|| project()));
    group.finish();
}


criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);