use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rstats::{Vectors,genvec};

pub fn criterion_benchmark(c: &mut Criterion) {
   let pts = genvec(20,30,3,5);
   c.bench_function("GMedian", |b| b.iter(||
      pts.as_slice().gmedian(black_box(20),black_box(1_e-5))));
   c.bench_function("NMedian", |b| b.iter(||
      pts.as_slice().nmedian(black_box(20),black_box(1_e-5))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
