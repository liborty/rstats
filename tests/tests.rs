#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use devtimer::DevTime;
use anyhow::Result;
use indxvec::{Indices,merge::*};
use rstats::{Stats,MutVecg,Vecg,VecVec,VecVecg,Vecu8, Vecf64};
use rstats::{wv,wi,printvv,genvec,genvecu8,i64tof64,tof64};

pub const EPS:f64 = 1e-10;
#[test]
fn u8() -> Result<()> {
   let v1 = vec![1_u8,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6]; 
   println!("\nv1: {}",wv(&v1));
   let v2 = vec![1_u8,2,2,3,3,3,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2]; 
   println!("v2: {}",wv(&v2)); 
   println!("Lexical order v1<v2: {}", wi(&(v1<v2)));
   println!("Entropy 1:\t{}",wi(&v1.entropy()));
   println!("Entropy 2 gen:\t{}",wi(&v2.entropy())); // generic
   println!("Entropy 2 u8:\t{}",wi(&v2.entropyu8())); // u8
   println!("Euclid's dist:\t{}",wi(&v2.vdistu8(&v1)));
   println!("Cityblock dist:\t{}",wi(&v2.cityblockd(&v1)));
   println!("Joint Entropy gen: {}",wi(&v1.jointentropy(&v2)));   
   println!("Joint Entropy u8:  {}",wi(&v1.jointentropyu8(&v2)));
   println!("Generic dep.:\t{}",wi(&v1.dependence(&v2))); // generic
   println!("Dependence u8:\t{}",wi(&v1.dependenceu8(&v2))); // u8
   println!("Cos:\t\t{}",wi(&v1.cosineu8(&v2)));
   println!("Pearson Correlation:  {}",wi(&v1.correlation(&v2)));  
   println!("Spearman Correlation: {}",wi(&v1.spearmancorr(&v2)));
   let d = 5_usize;
   let n = 7_usize;
   println!("Testing on a random set of {} points in {} d space:",wi(&n),wi(&d));
   let pt = genvecu8(d,n,5,7); // random test data 
   let cov = pt.covar(&pt.acentroid());
   let com = pt.covar(&pt.gmedian(EPS));
   println!("Covariances:\n{}",wv(&cov));
   println!("Comediances:\n{}",wv(&com));
   println!("Their Distance: {}\n",wi(&cov.vdist(&com)));
   Ok(())
}
#[test]
fn fstats() -> Result<()> { 
    let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.]; 
    println!("\n{}",wv(&v1)); 
    // println!("Linear transform:\n{}",wv(&v1.lintrans()));
    println!("Arithmetic mean:{}",wi(&v1.amean().unwrap()));
    println!("Geometric mean:\t{}",wi(&v1.gmean().unwrap()));
    println!("Harmonic mean:\t{}",wi(&v1.hmean().unwrap()));
    println!("Magnitude:\t{}",wi(&v1.vmag()));
    println!("Arithmetic {}",wi(&v1.ameanstd().unwrap()));
    println!("Geometric  {}",wi(&v1.gmeanstd().unwrap()));
    println!("Harmonic   {}",wi(&v1.hmeanstd().unwrap()));
    println!("Autocorrelation:{}",wi(&v1.autocorr()));
    println!("{}",v1.median().unwrap());
    println!("Mad:\t\t {}\n",wi(&v1.mad().unwrap()));
    Ok(())
 }
#[test]
fn ustats() -> Result<()> { 
   let v1 = vec![1_u8,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; 
   println!("\n{}",wv(&v1)); 
   // println!("Linear transform:\n{}",wv(&v1.lintrans()));
   println!("Arithmetic mean:{}",wi(&v1.amean().unwrap()));
   println!("Geometric mean:\t{}",wi(&v1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",wi(&v1.hmean().unwrap()));
   println!("Magnitude:\t{}",wi(&v1.vmag()));
   println!("Arithmetic {}",wi(&v1.ameanstd().unwrap()));
   println!("Geometric  {}",wi(&v1.gmeanstd().unwrap()));
   println!("Harmonic   {}",wi(&v1.hmeanstd().unwrap()));
   println!("Autocorrelation:{}",wi(&v1.autocorr()));
   println!("{}\n",v1.median().unwrap());
   Ok(())
}

#[test]
/// &[i64] requires explicit recast
fn intstats() -> Result<()> { 
   let v = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; 
   println!("\n{}",wv(&v));
   let v1 = i64tof64(&v); // conversion here
   // println!("Linear transform:\n{}",wv(&v1.lintrans()));
   println!("Arithmetic mean:{}",wi(&v1.amean().unwrap()));
   println!("Geometric mean:\t{}",wi(&v1.gmean().unwrap()));
   println!("Harmonic mean:\t{}",wi(&v1.hmean().unwrap()));
   // println!("Magnitude:\t{}",wi(&v1.vmag()));
   println!("Arithmetic {}",wi(&v1.ameanstd().unwrap()));
   println!("Geometric  {}",wi(&v1.gmeanstd().unwrap()));
   println!("Autocorrelation:{}",wi(&v1.autocorr()));
   println!("{}\n",v1.median().unwrap());
   Ok(())
}
#[test]
/// Generic implementation
fn genericstats() -> Result<()> { 
   let v = vec![1_i32,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; 
   println!("\n{}",wv(&v)); 
   println!("Arithmetic\t{}",wi(&v.ameanstd().unwrap()));
   println!("Geometric\t{}",wi(&v.gmeanstd().unwrap()));
   println!("Harmonic\t{}",wi(&v.hmeanstd().unwrap()));
   println!("Weighted Arit.\t{}",wi(&v.awmeanstd().unwrap()));
   println!("Weighted Geom.\t{}",wi(&v.gwmeanstd().unwrap()));
   println!("Weighted Harm.\t{}",wi(&v.hwmeanstd().unwrap()));
   println!("Autocorrelation:{}",wi(&v.autocorr()));
   println!("{}\n",&v.median().unwrap()); 
   Ok(())
}
#[test]
fn vecg() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,10.,10.,13.,14.,15.];
   println!("v1: {}",wv(&v1));
   let v2 = vec![1_f64,14.,2.,13.,3.,12.,4.,11.,5.,10.,6.,9.,7.,1.,15.];
   println!("v2: {}",wv(&v2)); 
   println!("Lexical order v1<v2:\t{}", wi(&(v1<v2)));
   println!("Pearson's Correlation:\t{}",wi(&v1.correlation(&v2))); 
   println!("Median Correlation:\t{}",wi(&v1.mediancorr(&v2))); 
   println!("Kendall's Correlation:\t{}",wi(&v1.kendalcorr(&v2)));  
   println!("Spearman's Correlation:\t{}",wi(&v1.spearmancorr(&v2)));  
   println!("Euclidian distance:\t{}",wi(&v1.vdistf64(&v2)));
   println!("Cityblock distance:\t{}",wi(&v1.cityblockd(&v2)));   
   println!("Vector difference: {}",wv(&v1.vsub(&v2))); 
   println!("Vector sum: {}",wv(&v1.vadd(&v2)));  
   println!("Scalar product:\t\t{}",wi(&v1.dotp(&v2)));
   println!("Parallelogram area:\t{}",wi(&v1.varea(&v2)));
   println!("Arc area:\t\t{}",wi(&v1.varc(&v2))); 
   println!("Entropy v1:\t\t{}",wi(&v1.entropy()));
   println!("Entropy v2:\t\t{}",wi(&v2.entropy()));
   println!("Joint Entropy:\t\t{}",wi(&v1.jointentropy(&v2)));
   println!("Dependence:\t\t{}",wi(&v1.dependence(&v2)));
   println!("Cosine:\t\t\t{}",wi(&v1.cosine(&v2))); 
   println!("Cosine of ranks:\t{}",
      wi(&rank(&v1,true).indx_to_f64().cosine(&rank(&v2,true).indx_to_f64())));        
   println!("Cos Similarity [0,1]:\t{}",wi(&v1.vsim(&v2)));
   println!("Cos Dissimilarity:\t{}",wi(&v1.vdisim(&v2))); 
   println!("[1,2,3].kron(&[4,5]):\t{}", wv(&[1,2,3].kron(&[4,5])));
   let outerp = [1,2,3].outer(&[4,5,6,7]);
   println!("[1,2,3].outer(&[4,5,6,7]): "); printvv(outerp);
   // println!("Transposed: "); printvv([1,2,3].outer(&[4,5,6,7]).transpose()); 
   Ok(())
}

#[test]
/// Trend between two data sets in space of the same dimensions but 
/// numbers of points can differ
fn trend() -> Result<()> {
   let d = 7_usize;
   let pts1 = genvec(d,28,13,19); // random test data 
   let pts2 = genvec(d,38,23,31);
   println!("\nTrend vector:\n{}\n",wv(&pts1.trend(EPS,pts2)));
   Ok(())
}

#[test]
fn vecvec() -> Result<()> { 
   let d = 55_usize;
   let n = 13_usize;
   println!("testing on a random set of {} points in {} dimensional space",wi(&n),wi(&d));
   let mut weights = Vec::new();
   for i in 1..n+1 { weights.push(i as f64) }; // create test weights data
   let pt = genvecu8(d,n,5,17); // random u8 test data
   println!("Joint entropy: {}",wi(&pt.jointentropyn()) );  
   println!("Statistical component wise dependence: {}",wi(&pt.dependencen()) ); 
   let outcomes = &genvecu8(d,1,7,19)[0];  
   println!("Dependencies of outcomes: {}",wv(&pt.dependencies(outcomes))); 
   println!("Correlations with outcomes: {}",wv(&pt.correlations(outcomes)));
   let (eccstd,eccmed,eccecc) = pt.eccinfo(EPS); 
   // let me = pt.emedoid(EPS);
   let medoid = &pt[eccecc.minindex];
   let outlier = &pt[eccecc.maxindex];
   let hcentroid = pt.hcentroid();
   let gcentroid = pt.gcentroid();
   let acentroid = pt.acentroid(); 
   let firstp = pt.firstpoint();
   let median = pt.smedian(EPS);
   // let outlier = &pt[md.maxindex]; 
   // let eoutlier = &pt[me.maxindex];
 
   let dists = pt.distsums();
   let md = minmax(&dists);
   println!("\nMedoid and Outlier Total Distances:\n{}",md ); 
   println!("Total Distances {}",dists.ameanstd().unwrap());
   println!("Total Distances {}\n",dists.median().unwrap());
   println!("HCentroid's total distances:\t{}",wi(&pt.distsum(&hcentroid)));
   println!("GCentroid's total distances:\t{}",wi(&pt.distsum(&gcentroid)));
   println!("ACentroid's total distances:\t{}",wi(&pt.distsum(&acentroid)));
   println!("Median's total distances:\t{}",wi(&pt.distsum(&median)));  
   println!("Outlier's distance to Medoid:\t{}",wi(&outlier.vdist(medoid)));      
   println!("Outlier's radius (from Median):\t{}",wi(&outlier.vdist(&median)));  
   println!("Medoid's radius (from Median):\t{}",wi(&medoid.vdist(&median)));

   println!("\nMedoid and outlier radii (eccentricities):\n{}
   \nRadii {}\nRadii {}",eccecc,eccstd,eccmed); 
   println!("HCentroid's radius: {}",wi(&hcentroid.vdist(&median)));
   println!("GCentroid's radius: {}",wi(&gcentroid.vdist(&median)));
   println!("ACentroid's radius:  {}",wi(&acentroid.vdist(&median)));
   println!("Firstpoint's radius: {}",wi(&firstp.vdist(&median)));
   println!("Median's error*{:e}: {}",1_f64/EPS,wi(&(pt.eccnonmember(&median).vmag()/EPS)));
   // let zmed = pt.translate(&median); // zero median transformed data
   // println!("Median's error:\t{}\n",wi(&zmed.gmedian(EPS).vmag()));

   let seccs = pt.sortedeccs(true,&median); 
   // println!("\nSorted eccs: {}\n", wv(&seccs));
   let lqcnt = binsearch(&seccs,eccmed.lquartile);
   println!("Inner quarter of points: {} max radius: {}", wi(&lqcnt), wi(&seccs[lqcnt-1]));
   let medcnt = binsearch(&seccs,eccmed.median);
   println!("Inner half of points:    {} max radius: {}", wi(&medcnt), wi(&seccs[medcnt-1]));   
   let uqcnt = binsearch(&seccs,eccmed.uquartile);
   println!("Inner three quarters:    {} max radius: {}", wi(&uqcnt), wi(&seccs[uqcnt-1]));
   // create pretend median of medians
   // let medmed = vec![0.5_f64;n];
   // let (se, cpdf) = 
   //  pt.wsortedcos(&medmed,&medmed.vunit(), &weights);
   //  println!("Sorted coses:\n{}\ncpdf:\n{}\n",wv(&se),wv(&cpdf));
   Ok(())
}

#[test]

fn geometric_medians() -> Result<()> {
    const ITERATIONS:usize = 10;
    let n = 100_usize;
    let d = 1000_usize;
    println!("timing {} medians of {} points in {} dimensions",wi(&ITERATIONS),wi(&n),wi(&d)); 
  
   let mut sumg = 0_f64;
   let mut sumtime = 0_u128;
   let mut timer = DevTime::new_simple();
 
   for i in 1..ITERATIONS {
      let pts = genvecu8(d,n,i as u32,5*i as u32);
      timer.start();
      let gm = pts.gmedian(EPS);
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      // sumg += pts.distsum(&gm)
      sumg += pts.eccnonmember(&gm).vmag();    
   }
   // sumg /= (ITERATIONS*n*d) as f64;
   println!("Gmedian err/eps: {}\tns: {:>12}",wi(&(sumg/EPS)),&sumtime); 

   sumg = 0_f64;
   sumtime = 0_u128;
   timer = DevTime::new_simple();
 
   for i in 1..ITERATIONS {
      let pts = genvecu8(d,n,i as u32,5*i as u32);
      timer.start();
      let gm = pts.smedian(EPS);
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      // sumg += pts.distsum(&gm)
      sumg += pts.eccnonmember(&gm).vmag();    
   }
   // sumg /= (ITERATIONS*n*d) as f64;
   println!("Smedian err/eps: {}\tns: {:>12}",wi(&(sumg/EPS)),&sumtime); 
   
   sumg = 0_f64;
   sumtime = 0_u128;
   timer = DevTime::new_simple();
 
   for i in 1..ITERATIONS {
      let pts = genvecu8(d,n,i as u32,5*i as u32);
      timer.start();
      let gm = pts.acentroid();
      timer.stop();
      sumtime += timer.time_in_nanos().unwrap();
      // sumg += pts.distsum(&gm)
      sumg += pts.eccnonmember(&gm).vmag(); 
   } 
   // sumg /= (ITERATIONS*n*d) as f64;  
   println!("Acentroid er/eps: {}\tns: {:>12}",wi(&(sumg/EPS)),sumtime); 
    Ok(())  
 }
