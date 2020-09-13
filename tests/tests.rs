#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use anyhow::{Result};
use rstats::{RStats,GreenIt,GreenVec,MutVectors,Vectors,genvec};
use devtimer::DevTime;



#[test]
fn fstats() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   let s1 = v1.as_slice();
   println!("\n{:?}",s1);
   println!("Arithmetic mean:{}",GreenIt(s1.amean().unwrap()));
   println!("Geometric mean:\t{}",GreenIt(s1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",GreenIt(s1.hmean().unwrap()));
   println!("Magnitude:\t{}",GreenIt(s1.vmag()));
   println!("Arithmetic\t{}",s1.ameanstd().unwrap());
   println!("Geometric\t{}",s1.gmeanstd().unwrap());
   println!("Autocorrelation:{}",GreenIt(s1.autocorr().unwrap()));
   println!("{}",s1.median().unwrap());
   let v2 = vec![1_f64,14.,2.,13.,3.,12.,4.,11.,5.,10.,6.,9.,7.,8.];
   let s2 = v2.as_slice();
   println!("\n{:?}",s2);
   println!("Ranking: {}",GreenVec(s2.ranks().unwrap()));
   println!("Pearson's Correlation:\t{}",GreenIt(s1.correlation(s2).unwrap())); 
   println!("Kendall's Correlation:\t{}",GreenIt(s1.kendalcorr(s2).unwrap()));  
   println!("Spearman's Correlation:\t{}",GreenIt(s1.spearmancorr(s2).unwrap()));     
   println!("Scalar product: {}",GreenIt(s1.dotp(s2)));
   println!("Euclidian distance: {}",GreenIt(s1.vdist(s2)));
   println!("Magnitude of difference: {}",GreenIt(s1.vsub(s2).as_slice().vmag()));   
   println!("Vector difference:\n{}",GreenVec(s1.vsub(s2))); 
   println!("Vector addition:\n{}",GreenVec(s1.vadd(s2)));   
   Ok(())
}

#[test]
fn intstats() -> Result<()> { 
   let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   let s1 = v1.as_slice();
   println!("\n{:?}",&s1);
   println!("Arithmetic mean:{}",GreenIt(s1.amean().unwrap()));
   println!("Geometric mean:\t{}",GreenIt(s1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",GreenIt(s1.hmean().unwrap()));
   println!("Arithmetic\t{}",s1.ameanstd().unwrap());
   println!("Geometric\t{}",s1.gmeanstd().unwrap());
   println!("{}",s1.median().unwrap());   
   Ok(())
}

#[test]
fn multidimensional() -> Result<()> { 
   let d = 5_usize;
   let pt = genvec(d,20,3,13); // random test data 5x20
   let pts = pt.as_slice();
   let (med,medi,outd,outi) = pts.medoid(d);
   let centroid = pts.acentroid(d);
   let median = pts.nmedian(d, 1e-5).unwrap();

   println!("\nSum of Outlier distances:\t{} Index: {}",GreenIt(outd),GreenIt(outi as f64));
   println!("Sum of Medoid distances:\t{} Index: {}",GreenIt(med),GreenIt(medi as f64));
   println!("Sum of Centroid distances:\t{}",GreenIt(pts.distsum(d,&centroid)));
   println!("Sum of Median distances:\t{}\n",GreenIt(pts.distsum(d,&median)));

   println!("Outlier eccentricity:\t{}",GreenIt(pts.eccentr(d,outi)));
   println!("Medoid ecentricity:\t{}",GreenIt(pts.eccentr(d,medi)));
   println!("Centroid ecentricity:\t{}",GreenIt(pts.ecc(d,&centroid)));   
   println!("Median eccentricity:\t{}\n",GreenIt(pts.ecc(d,&median)));
   println!("MOE (Median of eccentricities)\n{}",pts.moe(d));  
   Ok(())
}
#[test]
fn difficult_data() -> Result<()> {
   let pts:Vec<f64> = vec![0.,0.,0.,0., 1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
      -1.,0.,0.,0., 0.,-1.,0.,0., 0.,0.,-1.,0., 0.,0.,0.,-1.];
   let gm = pts.as_slice().nmedian(4, 1e-5).unwrap();
   println!("\nMedian residual error: {}",GreenIt(pts.as_slice().ecc(4,&gm)));
   Ok(())
}

#[test]
fn medians() -> Result<()> {
   const ITERATIONS:u32 = 10;
   let mut sumg = 0_f64;
   let mut sumtime = 0_u128;
   let mut timer = DevTime::new_simple();

   for i in 1..ITERATIONS {
      let pts = genvec(2,7000,i,2*i);
      timer.start();
      let gm = pts.as_slice().nmedian(2, 1e-5).unwrap();
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      sumg += pts.as_slice().ecc(2,&gm);
   }
   println!("\nSum of {} medians residual errors {} in {} ns",
      GreenIt(ITERATIONS),GreenIt(sumg),GreenIt(sumtime));     
   Ok(())  
}
