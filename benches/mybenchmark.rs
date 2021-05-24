use rstats::{VecVecf64,Vecf64};
use rstats::functions::{GI,genvec};
use devtimer::DevTime;

fn main() {

   const ITERATIONS:usize = 100;
   const D:usize = 50; // dimensions
   const N:usize = 500; // number of points
   const EPS:f64 = 1e-3;
   let mut sumg = 0_f64; // sum of error distances
   let mut sumtime = 0_u128;

   let mut cmplx = DevTime::new_complex();

   cmplx.create_timer("nmedian").unwrap();
   for i in 0..ITERATIONS {
      // random generator dimensions and seeds
      let pts = genvec(D,N,(i+1) as u32,N as u32); 
      // test over varying data and inreasing accuracies
      cmplx.start_timer("nmedian").unwrap();  
      let g = pts.as_slice().nmedian(EPS);
      cmplx.stop_timer("nmedian").unwrap();
      sumtime += cmplx.time_in_nanos("nmedian").unwrap();
      sumg += pts.eccnonmember(&g).vmag();  
   };
   println!("\nSum of nmedian residual errors: {} time: {} ", GI(sumg), GI(sumtime));
  
/*
  // Or we can iterate through all timers
  for (tname, timer) in cmplx.iter() {
    println!("{} - {} ns", tname, timer.time_in_nanos().unwrap());
  }

  // Or we can print results in the default '{timername} - {time} ns' format
  cmplx.print_results();
  cmplx.clear_timers();
}
*/
}