#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use anyhow::{Result};
use rstats::{RStats,MutVectors,Vectors,genvec,green};
use devtimer::DevTime;

#[test]
fn fstats() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   let s1 = v1.as_slice();
   println!("\n{:?}",s1);
   println!("Arithmetic mean:{}",green(s1.amean().unwrap()));
   println!("Geometric mean:\t{}",green(s1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",green(s1.hmean().unwrap()));
   println!("Magnitude:\t{}",green(s1.vmag()));
   println!("Arithmetic\t{}",s1.ameanstd().unwrap());
   println!("Geometric\t{}",s1.gmeanstd().unwrap());
   println!("Autocorrelation:{}",green(s1.autocorr().unwrap()));
   println!("{}",s1.median().unwrap());
   let mut v2 = v1.clone();
   v2.reverse();
   let s2 = v2.as_slice();
   println!("\n{:?}",s2);
   println!("Correlation:{}",green(s1.correlation(s2).unwrap()));   
   println!("Scalar product: {}",green(s1.dotp(s2)));
   println!("Euclidian distance: {}",green(s1.vdist(s2)));
   println!("Magnitude of difference: {}",green(s1.vsub(s2).as_slice().vmag()));   
   println!("Vector difference:\n\x1B[01;92m{:?}\x1B[0m",s1.vsub(s2)); 
   println!("Vector addition:\n\x1B[01;92m{:?}\x1B[0m",s1.vadd(s2));   
   Ok(())
}

#[test]
fn intstats() -> Result<()> { 
   let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   let s1 = v1.as_slice();
   println!("\n{:?}",&s1);
   println!("Arithmetic mean:{}",green(s1.amean().unwrap()));
   println!("Geometric mean:\t{}",green(s1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",green(s1.hmean().unwrap()));
   println!("Arithmetic\t{}",&s1.ameanstd().unwrap());
   println!("Geometric\t{}",s1.gmeanstd().unwrap());
   println!("Autocorrelation:{}",green(s1.autocorr().unwrap()));
   println!("{}",s1.median().unwrap());
   let mut v2 = v1.clone(); v2.reverse();
   println!("Correlation:{}",green(s1.icorrelation(v2.as_slice()).unwrap()));   
   Ok(())
}

#[test]
fn mc() -> Result<()> { 
   let pts = genvec(5,20,3,13);
   let (dist,indx) = pts.as_slice().medoid(5).unwrap();
   println!("Sum of Medoid distances: {} Index: {}",green(dist),indx);
   let centroid = pts.as_slice().arcentroid(5);
   println!("Sum of Centroid distances: {}",green(pts.as_slice().distsum(5,&centroid)));
   let (ds,gm) = pts.as_slice().gmedian(5, 1e-5).unwrap();
   println!("Sum of Gmedian distances: {}",green(ds));
   println!("Gmedian eccentricity (residual error): {}",
      green(pts.as_slice().exteccentr(5,&gm)));
   Ok(())
}
#[test]
fn difficult_data() -> Result<()> {
   let pts:Vec<f64> = vec![0.,0.,0.,0., 1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
      -1.,0.,0.,0., 0.,-1.,0.,0., 0.,0.,-1.,0., 0.,0.,0.,-1.];
   let (ds,_gm) = pts.as_slice().nmedian(4, 1e-5).unwrap();
   println!("\nSum of Nmedian distances:  \x1B[01;92m{}\x1B[0m",ds);
   let (ds,gm) = pts.as_slice().gmedian(4, 1e-5).unwrap();
   println!("Sum of Gmedian distances: {}",green(ds));
   Ok(())
}

#[test]
fn medians() -> Result<()> {
   let mut sumg = 0_f64;
   let mut timer = DevTime::new_simple();

   let pts = genvec(5,20,3,13);
   println!();
   for i in 0 .. 20 {
      let ec = pts.as_slice().eccentr(5,i);
      println!("Eccentricity:{} Index:{}",ec,i) }
 
   timer.start();
   for i in 1..6 {
      let pts = genvec(2,7000,i,2*i);
      let (dg,_) = pts.as_slice().gmedian(2, 1e-2).unwrap();
      sumg += dg;
   }
   timer.stop();
   println!("\nSum of gmedian distances and time in ns: \x1B[01;92m{}\t{}\x1B[0m",
      sumg,timer.time_in_nanos().unwrap());
   sumg = 0_f64;

   timer.start();
   for i in 1..6 {
      let pts = genvec(2,7000,i,2*i);
      let (dg,_) = pts.as_slice().nmedian(2, 1e-2).unwrap();
      sumg += dg;
      }
   timer.stop();
   println!("Sum of nmedian distances and time in ns: \x1B[01;92m{}\t{}\x1B[0m",
      sumg,timer.time_in_nanos().unwrap());  
   Ok(())
}
