use rstats::{Vectors,genvec,green};

/*
use criterion::{black_box, criterion_group, criterion_main, Criterion}; //, BenchmarkId}
fn compare_medians(c: &mut Criterion) {
   let mut group = c.benchmark_group("Medians");
   let d = 20;
   let pts = genvec(d,60,5,7);
   for i in [0usize,1usize,2usize].iter() {
      group.bench_with_input(BenchmarkId::new("gmedian", i), i, 
           |b, i| b.iter(|| {
              let p = pts.get(i*d .. (i+1)*d).unwrap();
              p.gmedian(black_box(d),black_box(1_e-7)) })); };
   for i in [0usize,1usize,2usize].iter() {          
      group.bench_with_input(BenchmarkId::new("nmedian", i), i, 
           |b, i| b.iter(|| {
              let p = pts.get(i*d .. (i+1)*d).unwrap();
              p.gmedian(black_box(d),black_box(1_e-7)) }));
   };
   group.finish();
}
criterion_group!(benches, compare_medians);
criterion_main!(benches);

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId}
pub fn criterion_benchmark(c: &mut Criterion) {
   let pts = genvec(2,2000,11,13);
   c.bench_function("GMedian", |b| b.iter(||
      pts.as_slice().gmedian(black_box(2),black_box(1_e-1))));
   c.bench_function("NMedian", |b| b.iter(||
      pts.as_slice().nmedian(black_box(2),black_box(1_e-1))));
}
criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
*/
/*
use devtimer::run_benchmark;
pub fn main() {
   const D:usize = 50;
   const N:usize = 100; 

   let res_gmedian = run_benchmark(10, |i| {
      let pts = genvec(D,N,(i+1) as u32,N as u32);
      let (_dg,_) = pts.as_slice().gmedian(D, 1./((i+1)as f64).powi(i as i32)).unwrap();
   });
   println!("gmedian results");
   res_gmedian.print_stats();

   let res_nmedian = run_benchmark(10, |i| {
      let pts = genvec(D,N,(i+1) as u32,N as u32);
      let (_dg,_) = pts.as_slice().nmedian(D, 1./((i+1)as f64).powi(i as i32)).unwrap();
   });
   println!("nmedian results"); 
   res_nmedian.print_stats()
}
*/
use devtimer::DevTime;

fn main() {

   const ITERATIONS:usize = 100;
   const D:usize = 50; // dimensions
   const N:usize = 500; // number of points
   const EPS:f64 = 1e-3;
   let mut sumg = 0_f64; // sum of error distances

   let mut cmplx = DevTime::new_complex();

   cmplx.create_timer("gmedian").unwrap();
   cmplx.start_timer("gmedian").unwrap();
   for i in 0..ITERATIONS {
      // random generator dimensions and seeds
      let pts = genvec(D,N,(i+1) as u32,N as u32); 
      // test over varying data and inreasing accuracies
      let (dg,_) = pts.as_slice().gmedian(D,EPS).unwrap();
      sumg += dg;  
   };
   cmplx.stop_timer("gmedian").unwrap();
   println!("\nSum of gmedian error distances: {}", green(sumg));
   sumg = 0_f64;
   cmplx.create_timer("nmedian").unwrap();
   cmplx.start_timer("nmedian").unwrap();
   for i in 0..ITERATIONS {
      // random generator dimensions and seeds
      let pts = genvec(D,N,(i+1) as u32,N as u32); 
      // test over varying data and inreasing accuracies
      let (dg,_) = pts.as_slice().nmedian(D,EPS).unwrap();
      sumg += dg;  
   };
   cmplx.stop_timer("nmedian").unwrap();
   println!("\nSum of nmedian error distances: {}", green(sumg));
  
   
/*
  // We can output a benchmark in this way
  println!(" `gmedian` took: {}", cmplx.time_in_micros("gmedian").unwrap());

  // Or we can iterate through all timers
  for (tname, timer) in cmplx.iter() {
    println!("{} - {} ns", tname, timer.time_in_nanos().unwrap());
  }
*/
  // Or we can print results in the default '{timername} - {time} ns' format
  cmplx.print_results();
  cmplx.clear_timers();
}
