#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use anyhow::{Result};
use rstats::{Stats,MutVectors,Vecf64,VecVec,Vecu8,Indices};
use rstats::functions::{GI,GV,genvec};
use devtimer::DevTime;

pub const EPS:f64 = 1e-7;
#[test]
fn entropy() -> Result<()> {
   let v1 = vec![1_u8,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6]; 
   println!("\n{:?}",v1);
   println!("Entropy: {}",GI(v1.entropy()));
   let v2 = vec![1_u8,2,2,3,3,3,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2]; 
   println!("{:?}",v2);
   println!("Entropy: {}",GI(v2.entropy()));
   println!("Joint E: {}",GI(v1.jointentropy(&v2)));
   println!("Dependence: {}\n",GI(v1.dependence(&v2)));
   Ok(())
}

#[test]
fn fstats() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.];
   println!("\n{:?}",v1);
   println!("Linear transform:\n{}",GV(v1.lintrans()));
   println!("Arithmetic mean:{}",GI(v1.amean().unwrap()));
   println!("Geometric mean:\t{}",GI(v1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",GI(v1.hmean().unwrap()));
   println!("Magnitude:\t{}",GI(v1.vmag()));
   println!("Arithmetic {}",v1.ameanstd().unwrap());
   println!("Geometric  {}",v1.gmeanstd().unwrap());
   println!("Autocorrelation:{}",GI(v1.autocorr()));
   println!("{}\n",v1.median().unwrap());
   Ok(())
}
#[test]
fn intstats() -> Result<()> { 
   let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
   println!("\n{:?}",v1);
   println!("Arithmetic mean:{}",GI(v1.amean().unwrap()));
   println!("Geometric mean:\t{}",GI(v1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",GI(v1.hmean().unwrap()));
   println!("Arithmetic {}",v1.ameanstd().unwrap());
   println!("Geometric {}",v1.gmeanstd().unwrap());
   println!("{}\n",v1.median().unwrap()); 
   Ok(())
}
#[test]
fn vecf64() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.];
   println!("\n{:?}",v1);
   let v2 = vec![1_f64,14.,2.,13.,3.,12.,4.,11.,5.,10.,6.,9.,7.,8.,15.];
   println!("{:?}",v2);
   println!("Rank:      {}",GV(v2.ranks().unwrap()));
   println!("Sort index:{}",GV(v2.mergesort(0,v2.len()))); 
   println!("R reversed:{}",GV(v2.mergerank().revindex()));    
   println!("Mergerank: {}",GV(v2.mergerank()));
   println!("Mrg.sorted:{}",GV(v2.sortm()));
   println!("Sorted:    {}",GV(v2.sortf()));
   println!("Pearson's Correlation:\t{}",GI(v1.correlation(&v2))); 
   println!("Kendall's Correlation:\t{}",GI(v1.kendalcorr(&v2)));  
   println!("Spearman's Correlation:\t{}",GI(v1.spearmancorr(&v2)));  
   println!("Cosine:\t\t\t{}",GI(v1.cosine(&v2))); 
   println!("Cosine of ranks:\t{}",
        GI(v1.ranks().unwrap().cosine(&v2.ranks().unwrap())));        
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
fn vecvec() -> Result<()> { 
   let d = 10_usize;
   let n = 100_usize;
   println!("testing on a random set of {} points in {} dimensional space",GI(n),GI(d));
   let pt = genvec(d,n,5,17); // random test data 
   let (med,medi,outd,outi) = pt.medoid();
   let (mede,medei,oute,outei) = pt.emedoid(EPS);
   let hcentroid = pt.hcentroid();
   let acentroid = pt.acentroid(); 
   let firstp = pt.firstpoint();
   let median = pt.gmedian(EPS);
   let outlier = &pt[outi]; 
   let eoutlier = &pt[outei];
   let zmed = pt.translate(&median); // zero median transformed data
  
   println!("\nSum of Outlier distances:\t{} Index: {}",GI(outd),GI(outi));
   println!("E-Outlier's distances:\t\t{}",GI(pt.distsuminset(outei)));   
   println!("Outlier's distance to Median:\t{}",GI(outlier.vdist(&median)));
   println!("E-Outlier's distance to Median:\t{}",GI(eoutlier.vdist(&median)));      
   println!("Sum of Medoid's distances:\t{} Index: {}",GI(med),GI(medi));
   println!("Sum of HCentroid's distances:\t{}",GI(pt.distsum(&hcentroid)));
   println!("Sum of ACentroid's distances:\t{}",GI(pt.distsum(&acentroid)));  
   println!("Sum of Median's distances:\t{}",GI(pt.distsum(&median)));
   let dists = pt.distsums();
   println!("Distances\t{}",dists.ameanstd().unwrap());
   println!("Distances\t{}\n",dists.median().unwrap());

   println!("Outlier's approx eccentricity:\t{}",GI(pt.eccmember(outi).vmag()));
   println!("E-Outlier's eccentricity:\t{} Index: {}",GI(oute),GI(outei));
   println!("E-Medoid's eccentricity:\t{} Index: {}",GI(mede),GI(medei));
   println!("Centroid's approx eccentricity:\t{}",GI(pt.eccnonmember(&acentroid).vmag()));
   println!("Firstpoint's app eccentricity:\t{}",GI(pt.eccnonmember(&firstp).vmag()));
   println!("Median's ecc (passed epsilon):\t{}",GI(pt.eccnonmember(&median).vmag()));
   println!("Median's error:\t{}",GI(zmed.gmedian(EPS).vmag()));
   let (mu,eccmed) = pt.moe(EPS);
   println!("Eccentricities\t{}",mu);  
   println!("Eccentricities\t{}\n",eccmed);  
   Ok(())
}

#[test]
fn trend() -> Result<()> {
   let d = 7_usize;
   let pt1 = genvec(d,28,13,19); // random test data 
   let pt2 = genvec(d,38,23,31);
   println!("\nTrend vector:\n{}\n",GV(pt1.trend(EPS,pt2)));
   Ok(())
}

#[test]
fn geometric_medians() -> Result<()> {
   const ITERATIONS:u32 = 20;
   let n = 700_usize;
   let d = 10_usize;
   println!("timing {} medians of {} points each in {} dimensions",GI(ITERATIONS),GI(n),GI(d)); 
   
   let mut timer = DevTime::new_simple();
   let mut sumg = 0_f64;
   let mut sumtime = 0_u128; 
   for i in 1..ITERATIONS {
      let pts = genvec(d,n,i,5*i);
      timer.start();
      let gm = pts.gmedian(EPS);
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      sumg += pts.distsum(&gm)    
   }

   println!("Gmedian errors: {} ns:\t{}",GI(sumg),GI(sumtime));   
   sumg = 0_f64;
   sumtime = 0_u128;
   timer = DevTime::new_simple();
 
   for i in 1..ITERATIONS {
      let pts = genvec(d,n,i,5*i);
      timer.start();
      let gm = pts.nmedian(EPS);
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      sumg += pts.distsum(&gm)
   }   
   println!("Nmedian errors: {}  ns:\t{}",GI(sumg),GI(sumtime));  
  
    Ok(())  
 }
