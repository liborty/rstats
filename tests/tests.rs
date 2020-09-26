#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use anyhow::{Result};
use rstats::{Stats,MutVectors,Vecf64,VecVec,Indices};
use rstats::functions::{GI,GV,genvec};
use devtimer::DevTime;

#[test]
fn fstats() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.];
   let v1r = v1.ranks().unwrap();
   println!("\n{:?}",v1);
   println!("Lin. trans: {}\n",GV(v1.lintrans()));
   println!("Arithmetic mean:{}",GI(v1.amean().unwrap()));
   println!("Geometric mean:\t{}",GI(v1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",GI(v1.hmean().unwrap()));
   println!("Magnitude:\t{}",GI(v1.vmag()));
   println!("Arithmetic {}",v1.ameanstd().unwrap());
   println!("Geometric {}",v1.gmeanstd().unwrap());
   println!("Autocorrelation:{}",GI(v1.autocorr()));
   println!("{}",v1.median().unwrap());
   let v2 = vec![1_f64,14.,2.,13.,3.,12.,4.,11.,5.,10.,6.,9.,7.,8.,15.];
   let v2r = v2.ranks().unwrap();
   println!("\n{:?}",v2);
   println!("Rank:      {}",GV(v2.ranks().unwrap()));
   println!("Sort index:{}",GV(v2.mergesort(0,v2.len()))); 
   println!("R reversed:{}",GV(v2.mergerank().revindex()));    
   println!("Mergerank: {}",GV(v2.mergerank()));
   println!("Mrg.sorted:{}",GV(v2.mergesort(0,v2.len()).unindex(&v2)));
   println!("Sorted:    {}",GV(v2.sortf()));
   println!("Pearson's Correlation:\t{}",GI(v1.correlation(&v2))); 
   println!("Kendall's Correlation:\t{}",GI(v1.kendalcorr(&v2)));  
   println!("Spearman's Correlation:\t{}",GI(v1.spearmancorr(&v2)));  
   println!("Cosine:\t\t\t{}",GI(v1.cosine(&v2))); 
   println!("Cosine of ranks:\t\t{}",GI(v1r.cosine(&v2r)));        
   println!("Euclidian distance:\t{}",GI(v1.vdist(&v2)));
   println!("Difference magnitude:\t{}",GI(v1.vsub(&v2).as_slice().vmag()));   
   println!("Vector difference: {}",GV(v1.vsub(&v2))); 
   println!("Vector addition:   {}",GV(v1.vadd(&v2)));  
   println!("Scalar product:\t\t{}",GI(v1.dotp(&v2)));
   println!("Parallelogram area:\t{}",GI(v1.varea(&v2))); 
   println!("Arc area:\t\t{}\n",GI(v1.varc(&v2)));
   Ok(())
}

#[test]
fn intstats() -> Result<()> { 
   let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
   println!("\n{:?}",v1);
   println!("Arithmetic mean:{}",GI(v1.amean().unwrap()));
   println!("Geometric mean:\t{}",GI(v1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",GI(v1.hmean().unwrap()));
   println!("Arithmetic\t{}",v1.ameanstd().unwrap());
   println!("Geometric\t{}",v1.gmeanstd().unwrap());
   println!("{}",v1.median().unwrap());   
   Ok(())
}

#[test]
fn vecvec() -> Result<()> { 
   let d = 5_usize;
   let n = 60_usize;
   println!("testing on a random set of {} points in {} dimensional space",GI(n),GI(d));
   let pt = genvec(d,n,7,13); // random test data 
   let (med,medi,outd,outi) = pt.medoid();
   let (mede,medei,oute,outei) = pt.emedoid();
   let centroid = pt.acentroid();
   let median = pt.nmedian(1e-5);
   let outlier = &pt[outi]; 
   let eoutlier = &pt[outei];
   let zmed = pt.translate(&median); // zero median transformed data
  
   println!("\nSum of Outlier distances:\t{} Index: {}",GI(outd),GI(outi));
   println!("Outlier's distance to Median:\t{}",GI(outlier.vdist(&median)));
   println!("Outlier's eccentricity:\t\t{}",GI(pt.eccentrinset(outi)));
   println!("Sum of Medoid's distances:\t{} Index: {}",GI(med),GI(medi));
   println!("Sum of Centroid's distances:\t{}",GI(pt.distsum(&centroid)));
   println!("Sum of Median's distances:\t{}\n",GI(pt.distsum(&median)));

   println!("E-Outlier's eccentricity:\t{} Index: {}",GI(oute),GI(outei));
   println!("E-Outlier's distance to Median:\t{}",GI(eoutlier.vdist(&median)));
   println!("E-Outlier's distances:\t\t{}",GI(pt.distsuminset(outei)));   
   println!("Medoid's eccentricity:\t\t{} Index: {}",GI(mede),GI(medei));
   println!("Centroid's eccentricity:\t{}",GI(pt.ecc(&centroid)));   
   println!("Median's eccentricity:\t\t{}",GI(pt.ecc(&median)));
   println!("Zero med's median magnitude:\t{}",GI(zmed.nmedian(1e-5).vmag()));
   let (mu,med) = pt.moe();
   println!("Eccentricities\t{}",mu);  
   println!("Eccentricities median\n{}",med);  
   Ok(())
}

#[test]
fn trend() -> Result<()> {
   let d = 7_usize;
   let pt1 = genvec(d,28,13,19); // random test data 
   let pt2 = genvec(d,38,23,31);
   println!("\nTrend vector:\n{}",GV(pt1.trend(1_e-5,pt2)));
   Ok(())
}

#[test]
fn medians() -> Result<()> {
   const ITERATIONS:u32 = 10;
   let n = 700_usize;
   let d = 10_usize;
   let mut sumg = 0_f64;
   let mut sumtime = 0_u128;
   let mut timer = DevTime::new_simple();
   println!("timing {} medians of {} points each in {} dimensions",GI(ITERATIONS),GI(n),GI(d)); 
   for i in 1..ITERATIONS {
      let pts = genvec(d,n,i,2*i);
      timer.start();
      let gm = pts.nmedian(1e-5);
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      sumg += pts.ecc(&gm);
   }
   println!("Sum of residual errors: {} in {} ns",GI(sumg),GI(sumtime));     
   Ok(())  
}
