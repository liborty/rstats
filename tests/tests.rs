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
   let d = 5_usize;
   let pt = genvec(d,20,3,13);
   let pts = pt.as_slice();
   let (dist,indx) = pts.medoid(d).unwrap();
   println!("Sum of Medoid distances: {} Index: {}",green(dist),indx);
   let centroid = pts.arcentroid(d);
   println!("Sum of Centroid distances: {}",green(pts.distsum(d,&centroid)));
   let (ecc,_) = pts.veccentr(d,&centroid).unwrap();
   println!("Centroid ecentricity: {}",green(ecc));
   let gm = pts.nmedian(d, 1e-5).unwrap();
   println!("Median eccentricity (residual error): {}",green(pts.ecc(d,&gm)));
   Ok(())
}
#[test]
fn difficult_data() -> Result<()> {
   let pts:Vec<f64> = vec![0.,0.,0.,0., 1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
      -1.,0.,0.,0., 0.,-1.,0.,0., 0.,0.,-1.,0., 0.,0.,0.,-1.];
   let gm = pts.as_slice().nmedian(4, 1e-5).unwrap();
   println!("\nMedian residual error: {}",green(pts.as_slice().ecc(4,&gm)));
   Ok(())
}

#[test]
fn medians() -> Result<()> {
   let mut sumg = 0_f64;
   let mut sumtime = 0_u128;
   let mut timer = DevTime::new_simple();

   for i in 1..6 {
      let pts = genvec(2,7000,i,2*i);
      timer.start();
      let gm = pts.as_slice().nmedian(2, 1e-5).unwrap();
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      sumg += pts.as_slice().ecc(2,&gm);
   }
   // timer.stop();
   println!("\nSum of median 7000 residual errors {} time in ns: {}",
      green(sumg),sumtime);     
   Ok(())  
}
