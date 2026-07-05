use criterion::{criterion_group, criterion_main, Criterion};
use petekstatic::uncertainty::{run, Distribution};

fn bench_mc(c: &mut Criterion) {
    let d = Distribution::lognormal(3.0, 0.4).unwrap();
    c.bench_function("mc run 100k lognormal", |b| {
        b.iter(|| run(100_000, 1, |r| d.sample(r)))
    });
}
criterion_group!(benches, bench_mc);
criterion_main!(benches);
