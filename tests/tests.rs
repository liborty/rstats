#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use anyhow::{Result};
use rstats::{RStats,MutVectors,Vectors};
use rstats::vimpls::{GreenIt,GreenVec,genvec,scalarecc,sortf};
use devtimer::DevTime;

#[test]
fn fstats() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   println!("\n{:?}",v1);
   println!("Arithmetic mean:{}",GreenIt(v1.amean().unwrap()));
   println!("Geometric mean:\t{}",GreenIt(v1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",GreenIt(v1.hmean().unwrap()));
   println!("Magnitude:\t{}",GreenIt(v1.vmag()));
   println!("Arithmetic\t{}",v1.ameanstd().unwrap());
   println!("Geometric\t{}",v1.gmeanstd().unwrap());
   println!("Autocorrelation:{}",GreenIt(v1.autocorr().unwrap()));
   println!("{}",v1.median().unwrap());
   let v2 = vec![1_f64,14.,2.,13.,3.,12.,4.,11.,5.,10.,6.,9.,7.,8.];
   println!("\n{:?}",v2);
   println!("Ranking: {}",GreenVec(v2.ranks().unwrap()));
   println!("Pearson's Correlation:\t{}",GreenIt(v1.correlation(&v2).unwrap())); 
   println!("Kendall's Correlation:\t{}",GreenIt(v1.kendalcorr(&v2).unwrap()));  
   println!("Spearman's Correlation:\t{}",GreenIt(v1.spearmancorr(&v2).unwrap()));     
   println!("Scalar product: {}",GreenIt(v1.dotp(&v2)));
   println!("Euclidian distance: {}",GreenIt(v1.vdist(&v2)));
   println!("Magnitude of difference: {}",GreenIt(v1.vsub(&v2).as_slice().vmag()));   
   println!("Vector difference:\n{}",GreenVec(v1.vsub(&v2))); 
   println!("Vector addition:\n{}",GreenVec(v1.vadd(&v2)));   
   Ok(())
}

#[test]
fn intstats() -> Result<()> { 
   let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   println!("\n{:?}",v1);
   println!("Arithmetic mean:{}",GreenIt(v1.amean().unwrap()));
   println!("Geometric mean:\t{}",GreenIt(v1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",GreenIt(v1.hmean().unwrap()));
   println!("Arithmetic\t{}",v1.ameanstd().unwrap());
   println!("Geometric\t{}",v1.gmeanstd().unwrap());
   println!("{}",v1.median().unwrap());   
   Ok(())
}

#[test]
fn multidimensional() -> Result<()> { 
   let d = 5_usize;
   let pt = genvec(d,20,3,13); // random test data 5x20
   let (med,medi,outd,outi) = pt.medoid(d);
   let (mede,medei,oute,outei) = pt.emedoid(d);
   let centroid = pt.acentroid(d);
   let median = pt.nmedian(d, 1e-5).unwrap();

   println!("\nSum of Outlier distances:\t{} Index: {}",GreenIt(outd),GreenIt(outi));
   println!("Distances from E-Outlier:\t{}",GreenIt(pt.distsuminset(d, outei)));
   println!("Sum of Medoid distances:\t{} Index: {}",GreenIt(med),GreenIt(medi));
   println!("Sum of Centroid distances:\t{}",GreenIt(pt.distsum(d,&centroid)));
   println!("Sum of Median distances:\t{}\n",GreenIt(pt.distsum(d,&median)));

   println!("E-Outlier eccentricity:\t{} Index: {}",GreenIt(oute),GreenIt(outei));
   println!("Outlier eccentricity:\t{}",GreenIt(pt.eccentr(d, outi)));
   println!("Medoid ecentricity:\t{} Index: {}",GreenIt(mede),GreenIt(medei));
   println!("Centroid ecentricity:\t{}",GreenIt(pt.ecc(d,&centroid)));   
   println!("Median eccentricity:\t{}\n",GreenIt(pt.ecc(d,&median)));
   println!("Median of eccentricities (MOE)\n{}",pt.moe(d));  
   Ok(())
}
#[test]
fn difficult_data() -> Result<()> {
   let pts:Vec<f64> = vec![0.,0.,0.,0., 1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
      -1.,0.,0.,0., 0.,-1.,0.,0., 0.,0.,-1.,0., 0.,0.,0.,-1.];
   let gm = pts.nmedian(4, 1e-5).unwrap();
   println!("\nMedian residual error: {}",GreenIt(pts.ecc(4,&gm)));
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
      let gm = pts.nmedian(2, 1e-5).unwrap();
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      sumg += pts.ecc(2,&gm);
   }
   println!("\nSum of {} medians residual errors {} in {} ns",
      GreenIt(ITERATIONS),GreenIt(sumg),GreenIt(sumtime));     
   Ok(())  
}
