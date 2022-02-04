#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use devtimer::DevTime;
use random_number::random_fill;
use anyhow::Result;
use indxvec::{GR,UNGR,here,tof64,Indices,Printing,merge::*};
use rstats::{Stats,MutVecg,Vecg,VecVec,VecVecg,Vecu8,Vecf64};
use rstats::{Med,i64tof64};

pub const EPS:f64 = 1e-10;

pub fn genvec(d: usize, n: usize) -> Vec<Vec<f64>> {
    if n * d < 1 { panic!("{}\n\tnon positive dimensions",here!()) }       
    let mut v: Vec<Vec<f64>> = Vec::with_capacity(n);
    for _i in 0..n {
        let mut pt = vec![0_f64;d];
        random_fill!(&mut pt); 
        v.push(pt)
    } // fills the lot with random numbers
    v
}

pub fn genvecu8(d: usize, n: usize) -> Vec<Vec<u8>> {
    if n * d < 1 { panic!("{}\n\tnon positive dimensions",here!()) }       
    let mut v: Vec<Vec<u8>> = Vec::with_capacity(n);
    for _i in 0..n {
        let mut pt:Vec<u8> = vec![0_u8;d];
        random_fill!(&mut pt); 
        v.push(pt)
    } // fills the lot with random numbers
    v
}

#[test]
fn u8() -> Result<()> {
   let v1 = vec![1_u8,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6]; 
   println!("\nv1: {}",(&v1).gr());
   let v2 = vec![1_u8,2,2,3,3,3,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2]; 
   println!("v2: {}",(&v2).gr()); 
   println!("Lexical order v1<v2: {}", (v1<v2).gr());
   println!("Entropy 1:\t{}",v1.entropy().gr());
   println!("Entropy 2 gen:\t{}",v2.entropy().gr()); // generic
   println!("Entropy 2 u8:\t{}",v2.entropyu8().gr()); // u8
   println!("Euclid's dist:\t{}",v2.vdistu8(&v1).gr());
   println!("Cityblock dist:\t{}",v2.cityblockd(&v1).gr());
   println!("Joint Entropy gen: {}",v1.jointentropy(&v2).gr());   
   println!("Joint Entropy u8:  {}",v1.jointentropyu8(&v2).gr());
   println!("Generic dependence:{}",v1.dependence(&v2).gr()); // generic
   println!("Dependence u8:\t{}",v1.dependenceu8(&v2).gr()); // u8
   println!("Cos:\t\t{}",v1.cosineu8(&v2).gr());
   println!("Median Correlation:  {}",v1.mediancorr(&v2).gr());  
   println!("Pearson Correlation:  {}",v1.correlation(&v2).gr());  
   println!("Spearman Correlation: {}",v1.spearmancorr(&v2).gr());
   let d = 5_usize;
   let n = 7_usize;
   println!("Testing on a random set of {} points in {} d space:",n,d);
   let pt = genvecu8(d,n);
   let cov = pt.covar(&pt.acentroid());
   println!("Covariances:\n{}",cov.gr());
   let com = pt.covar(&pt.gmedian(EPS));
   println!("Comediances:\n{}",com.gr());
   println!("Their Distance: {}",cov.vdist(&com).gr());
   let trpt = pt.transpose();
   println!("Column Dependencies:\n{}",trpt.crossfeatures(|v1,v2| v1.dependence(v2)).gr());
   println!("Column Correlations:\n{}",trpt.crossfeatures(|v1,v2| v1.mediancorr(v2)).gr());
   Ok(())
}
#[test]
fn fstats() -> Result<()> { 
    let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.]; 
    println!("\n{}",(&v1).gr()); 
    let v2 = revs(&v1);
    println!("{}",(&v2).gr()); 
    // println!("Linear transform:\n{}",v1.lintrans()));
    println!("Arithmetic mean:{}",v1.amean().unwrap().gr());
    println!("Geometric mean:\t{}",v1.gmean().unwrap().gr());
    println!("Harmonic mean:\t{}",v1.hmean().unwrap().gr());
    println!("Magnitude:\t{}",v1.vmag().gr());
    println!("Arithmetic {}",v1.ameanstd().unwrap().gr());
    println!("Geometric  {}",v1.gmeanstd().unwrap().gr());
    println!("Harmonic   {}",v1.hmeanstd().unwrap().gr());
    println!("Autocorrelation:{}",v1.autocorr().gr());
    println!("{}",v1.median()?.gr());
    println!("Mad:\t\t {}\n",v1.mad().unwrap().gr());
   println!("Entropy 1:\t{}",v1.entropy().gr());
   println!("Entropy 2:\t{}",v2.entropy().gr()); // generic
   println!("Euclid's dist:\t{}",v2.vdist(&v1).gr());
   println!("Cityblock dist:\t{}",v2.cityblockd(&v1).gr());
   println!("Joint Entropy: {}",v1.jointentropy(&v2).gr());   
   println!("Dependence:\t{}",v1.dependence(&v2).gr()); // generic
   let d = 5_usize;
   let n = 7_usize;
   println!("Testing on a random set of {} points in {} d space:",n,d);
   let pt = genvec(d,n); 
   println!("Covariances:\n{}",pt.covar(&pt.acentroid()).gr());
   println!("Covariances of zero median data:\n{}",pt.covar(&pt.gmedian(EPS)).gr());
   println!("Median covariances of zero median data:\n{}",pt.comed(&pt.gmedian(EPS)).gr());
   let trpt = pt.transpose();
   println!("Column Dependencies:\n{}",trpt.crossfeatures(|v1,v2| v1.dependence(v2)).gr());
   println!("Column Correlations:\n{}",trpt.crossfeatures(|v1,v2| v1.mediancorr(v2)).gr());
    Ok(())
 }

#[test]
fn ustats() -> Result<()> { 
   let v1 = vec![1_u8,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; 
   println!("\n{}",(&v1).gr()); 
   // println!("Linear transform:\n{}",v1.lintrans()));
   println!("Arithmetic mean:{}",v1.amean()?.gr());
   println!("Geometric mean:\t{}",v1.gmean()?.gr());
   println!("Harmonic mean:\t{}",v1.hmean()?.gr());
   println!("Magnitude:\t{}",v1.vmag());
   println!("Arithmetic {}",v1.ameanstd()?.gr());
   println!("Geometric  {}",v1.gmeanstd()?.gr());
   println!("Harmonic   {}",v1.hmeanstd()?.gr());
   println!("Autocorrelation:{}",v1.autocorr());
   println!("{}\n",v1.median()?.gr());
   Ok(())
}

#[test]
/// &[i64] requires explicit recast
fn intstats() -> Result<()> { 
   let v = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; 
   println!("\n{}",(&v).gr());
   let v1 = i64tof64(&v); // conversion here
   // println!("Linear transform:\n{}",v1.lintrans()));
   println!("Arithmetic mean:{}",v1.amean()?.gr());
   println!("Geometric mean:\t{}",v1.gmean()?.gr());
   println!("Harmonic mean:\t{}",v1.hmean()?.gr());
   // println!("Magnitude:\t{}",v1.vmag()));
   println!("Arithmetic {}",v1.ameanstd()?.gr());
   println!("Geometric  {}",v1.gmeanstd()?.gr());
   println!("Autocorrelation:{}",v1.autocorr().gr());
   println!("{}\n",v1.median()?.gr());
   Ok(())
}
#[test]
/// Generic implementation
fn genericstats() -> Result<()> { 
   let v = vec![1_i32,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; 
   println!("\n{}",(&v).gr()); 
   println!("Arithmetic\t{}",v.ameanstd()?.gr());
   println!("Geometric\t{}",v.gmeanstd()?.gr());
   println!("Harmonic\t{}",v.hmeanstd()?.gr());
   println!("Weighted Arit.\t{}",v.awmeanstd()?.gr());
   println!("Weighted Geom.\t{}",v.gwmeanstd()?.gr());
   println!("Weighted Harm.\t{}",v.hwmeanstd()?.gr());
   println!("Autocorrelation:{}",v.autocorr().gr());
   println!("{}\n",&v.median()?.gr()); 
   Ok(())
}
#[test]
fn vecg() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,10.,10.,13.,14.,15.];
   println!("v1: {}",(&v1).gr());
   let v2 = vec![1_f64,14.,2.,13.,3.,12.,4.,11.,5.,10.,6.,6.,7.,1.,15.];
   println!("v2: {}",(&v2).gr()); 
   println!("Lexical order v1<v2:\t{}", (v1<v2).gr());
   println!("Pearson's Correlation:\t{}",v1.correlation(&v2).gr()); 
   println!("Median Correlation:\t{}",v1.mediancorr(&v2).gr()); 
   println!("Kendall's Correlation:\t{}",v1.kendalcorr(&v2).gr());  
   println!("Spearman's Correlation:\t{}",v1.spearmancorr(&v2).gr());  
   println!("Euclidian distance:\t{}",v1.vdistf64(&v2).gr());
   println!("Cityblock distance:\t{}",v1.cityblockd(&v2).gr());   
   println!("Vector difference: {}",v1.vsub(&v2).gr()); 
   println!("Vector sum: {}",v1.vadd(&v2).gr());  
   println!("Scalar product:\t\t{}",v1.dotp(&v2).gr());
   println!("Parallelogram area:\t{}",v1.varea(&v2).gr());
   println!("Arc area:\t\t{}",v1.varc(&v2).gr()); 
   println!("Entropy v1:\t\t{}",v1.entropy().gr());
   println!("Entropy v2:\t\t{}",v2.entropy().gr());
   println!("Joint Entropy:\t\t{}",v1.jointentropy(&v2).gr());
   println!("Dependence:\t\t{}",v1.dependence(&v2).gr());
   println!("Cosine:\t\t\t{}",v1.cosine(&v2).gr()); 
   println!("Cosine of ranks:\t{}",
      rank(&v1,true).indx_to_f64().cosine(&rank(&v2,true).indx_to_f64()).gr());        
   println!("Cos Similarity [0,1]:\t{}",v1.vsim(&v2).gr());
   println!("Cos Dissimilarity:\t{}",v1.vdisim(&v2).gr()); 
   println!("[1,2,3].kron(&[4,5]):\t{}", [1,2,3].kron(&[4,5]).gr());
   let outerp = [1,2,3].outer(&[4,5,6,7]);
   println!("[1,2,3].outer(&[4,5,6,7]):\n{}", outerp.gr());
   // println!("Transposed: "); printvv([1,2,3].outer(&[4,5,6,7]).transpose()); 
   Ok(())
}

#[test]
/// Trend between two data sets in space of the same dimensions but 
/// numbers of points can differ
fn trend() -> Result<()> {
   let d = 7_usize;
   let pts1 = genvecu8(d,28);
   let pts2 = genvecu8(d,33); 
   println!("\nTrend vector:\n{}\n",pts1.trend(EPS,pts2).gr());
   Ok(())
}

#[test]
fn vecvec() -> Result<()> { 
   let d = 10_usize;
   let n = 90_usize;
   println!("testing on a random set of {} points in {} dimensional space",n,d);
   let mut weights = Vec::new();
   for i in 1..n+1 { weights.push(i as f64) }; // create test weights data
   let pt = genvec(d,n); 
   // println!("{}",pt.gr());
   println!("Set joint entropy: {}",pt.jointentropyn().gr());  
   println!("Set dependence:    {}",pt.dependencen().gr()); 
   let mut outcomes:Vec<u8> = vec![0_u8;n];
   random_fill!(&mut outcomes); // column vector 
   let transppt = pt.transpose(); 
   println!("Dependencies of outcomes: {}",transppt.dependencies(&outcomes).gr()); 
   println!("Correlations with outcomes: {}",transppt.correlations(&outcomes).gr());
   let (eccstd,eccmed,eccecc) = pt.eccinfo(EPS); 
   // let me = pt.emedoid(EPS);
   let medoid = &pt[eccecc.minindex];
   let outlier = &pt[eccecc.maxindex];
   let hcentroid = pt.hcentroid();
   let gcentroid = pt.gcentroid();
   let acentroid = pt.acentroid(); 
   let firstp = pt.firstpoint();
   let median = pt.gmedian(EPS);
   // let outlier = &pt[md.maxindex]; 
   // let eoutlier = &pt[me.maxindex];
 
   let dists = pt.distsums();
   let md = minmax(&dists);
   println!("\nMedoid and Outlier Total Distances:\n{}",md ); 
   println!("Total Distances {}",dists.ameanstd()?.gr());
   println!("Total Distances {}\n",dists.median()?);
   println!("HCentroid's total distances:\t{}",pt.distsum(&hcentroid).gr());
   println!("GCentroid's total distances:\t{}",pt.distsum(&gcentroid).gr());
   println!("ACentroid's total distances:\t{}",pt.distsum(&acentroid).gr());
   println!("Median's total distances:\t{}",pt.distsum(&median).gr());  
   println!("Outlier's distance to Medoid:\t{}",outlier.vdist(medoid).gr());      
   println!("Outlier's radius (from Median):\t{}",outlier.vdist(&median).gr());  
   println!("Medoid's radius (from Median):\t{}",medoid.vdist(&median).gr());

   println!("\nMedoid and outlier radii (eccentricities):\n{}Radii {}\nRadii {}",eccecc,eccstd,eccmed); 
   println!("HCentroid's radius: {}",hcentroid.vdist(&median).gr());
   println!("GCentroid's radius: {}",gcentroid.vdist(&median).gr());
   println!("ACentroid's radius:  {}",acentroid.vdist(&median).gr());
   println!("Firstpoint's radius: {}",firstp.vdist(&median).gr());
   println!("Median's error*{:e}: {}",1_f64/EPS,(pt.eccnonmember(&median).vmag()/EPS).gr());
   // let zmed = pt.translate(&median); // zero median transformed data
   // println!("Median's error:\t{}\n",zmed.gmedian(EPS).vmag()));

   let seccs = pt.sortedeccs(true,&median); 
   // println!("\nSorted eccs: {}\n", seccs));
   let lqcnt = binsearch(&seccs,eccmed.lquartile);
   println!("Inner quarter of points: {} within radius: {}", lqcnt.gr(), seccs[lqcnt-1].gr());
   let medcnt = binsearch(&seccs,eccmed.median);
   println!("Inner half of points:    {} within radius: {}", medcnt.gr(), seccs[medcnt-1].gr());   
   let uqcnt = binsearch(&seccs,eccmed.uquartile);
   println!("Inner three quarters:    {} within radius: {}", uqcnt.gr(), seccs[uqcnt-1].gr());
   // create pretend median of medians
   // let medmed = vec![0.5_f64;n];
   // let (se, cpdf) = 
   //  pt.wsortedcos(&medmed,&medmed.vunit(), &weights);
   //  println!("Sorted coses:\n{}\ncpdf:\n{}\n",se),cpdf));
   Ok(())
}

#[test]

fn geometric_medians() -> Result<()> {
    const ITERATIONS:usize = 10;
    let n = 100_usize;
    let d = 1000_usize;
    println!("timing {} medians of {} points in {} dimensions",ITERATIONS,n,d);     
  
   let mut sumg = 0_f64;
   let mut timerg = DevTime::new_simple();
   let mut sumq = 0_f64;
   let mut timerq = DevTime::new_simple();
   let mut summ = 0_f64;
   let mut timerm = DevTime::new_simple(); 
   let mut gm:Vec<f64>;
   for _i in 1..ITERATIONS {   
      let pts = genvec(d,n);
      let trpts = pts.transpose();
      timerg.start();
      gm = pts.gmedian(EPS);
      timerg.stop(); 
      sumg += pts.eccnonmember(&gm).vmag();
      timerq.start();
      gm = trpts.iter().map(|p| { let Med{median,..} = p.median().unwrap(); median }).collect();
      timerq.stop();   
      sumq += pts.eccnonmember(&gm).vmag(); 
      timerm.start();
      gm = pts.acentroid();
      timerm.stop(); 
      summ += pts.eccnonmember(&gm).vmag();       
   }  
   println!("Geometric md err/eps: {GR}{:17.5}\ts: {:9}{UNGR}",sumg/EPS,timerg.time_in_nanos().unwrap()as f64/1e9);
   println!("Arithm. mean err/eps: {GR}{:17.5}\ts: {:9}{UNGR}",summ/EPS,timerm.time_in_nanos().unwrap()as f64/1e9);
   println!("Quasi median err/eps: {GR}{:17.5}\ts: {:9}{UNGR}",sumq/EPS,timerq.time_in_nanos().unwrap()as f64/1e9);
    Ok(())  
 }
