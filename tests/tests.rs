use indxvec::{printing::*, Indices, Vecops,Printing};
use rstats::{error::RError,i64tof64, Stats, VecVec, VecVecg, Vecg, Vecu8};
use ran::{*,set_seeds};
use medians::{Median};
use times::{benchvvf64};

pub const EPS: f64 = 1e-3;

#[cfg(test)]

#[test]
fn u8() {
    let v1 = vec![
        1_u8, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
    ];
    println!("\nv1: {}", (&v1).gr());
    let v2 = vec![
        1_u8, 2, 2, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2,
    ];
    println!("v2: {}", (&v2).gr());
    println!("v1*v2:\t{}",v1.dotp(&v2).gr());
    println!("v1+v2: {}",v1.vadd(&v2).gr());   
    println!("v1-v2: {}",v1.vsub(&v2).gr());    
    println!("Lexical order v1<v2: {}", (v1 < v2).gr());
    println!("Entropy v1:\t{}", v1.entropy().gr());
    println!("Entropyu8 v1:\t{}", v1.entropyu8().gr());
    println!("Entropy v2:\t{}", v2.entropy().gr()); // generic
    println!("Entropyu8 v2:\t{}", v2.entropyu8().gr()); // u8 
    println!("Joint Entropy:  {}", v1.jointentropy(&v2).gr());
    println!("Joint Entropyu8:{}", v1.jointentropyu8(&v2).gr());
    println!("Dependence:   {}", v1.dependence(&v2).gr()); // generic
    println!("Dependenceu8: {}", v1.dependenceu8(&v2).gr()); // u8 
    let med =  v1.as_slice().median();
    println!("Median v1:    {} +-{}",med.gr(),v1.mad(med).gr());
    println!("{}",v1.medinfo());   
    let d = 5_usize;
    let n = 7_usize;
    println!("Testing on a random set of {} points in {} d space:", n, d);
    set_seeds(77777);
    let pt = Rnum::newu8().ranvv(d, n).getvvu8();
    let cov = pt.covar(&pt.acentroid());
    println!("Covariances:\n{}", cov.gr());
    let com = pt.covar(&pt.gmedian(EPS));
    println!("Comediances:\n{}", com.gr());
    println!("Their Distance: {}", cov.vdist(&com).gr());
    let trpt = pt.transpose();
    println!(
        "Column Dependencies:\n{}",
        trpt.crossfeatures(|v1, v2| v1.dependence(v2)).gr()
    );
    println!(
        "Column Correlations:\n{}",
        trpt.crossfeatures(|v1, v2| v1.mediancorr(v2)).gr()
    ); 
}

#[test]
fn fstats() {
    let v1 = vec![
        1_f64, 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.,
    ];
    println!("\n{}", (&v1).gr());
    let v2 = v1.revs();
    println!("{}", (&v2).gr());
    // println!("Linear transform:\n{}",v1.lintrans()));
    println!("Arithmetic mean:{}", v1.amean().unwrap().gr());
    println!("Geometric mean:\t{}", v1.gmean().unwrap().gr());
    println!("Harmonic mean:\t{}", v1.hmean().unwrap().gr());
    println!("Magnitude:\t{}", v1.vmag().gr());
    println!("Arithmetic {}", v1.ameanstd().unwrap().gr());
    println!("Geometric  {}", v1.gmeanstd().unwrap().gr());
    println!("Harmonic   {}", v1.hmeanstd().unwrap().gr());
    println!("Autocorrelation:{}", v1.autocorr().gr());
    let med =  v1.as_slice().median();
    println!("Median:\t\t{} +- {}",med.gr(),v1.mad(med).gr());
    println!("Entropy 1:\t{}", v1.entropy().gr());
    println!("Entropy 2:\t{}", v2.entropy().gr()); // generic
    println!("Euclid's dist:\t{}", v2.vdist(&v1).gr());
    println!("Cityblock dist:\t{}", v2.cityblockd(&v1).gr());
    println!("Joint Entropy:  {}", v1.jointentropy(&v2).gr());
    println!("Dependence:\t{}", v1.dependence(&v2).gr()); // generic
    let d = 5_usize;
    let n = 7_usize;
    println!("Testing on a random set of {} points in {} d space:", n, d);
    let pt = Rnum::newf64().ranvv(d, n).getvvf64();
    println!("Classical Covariances:\n{}", pt.covar(&pt.acentroid()).gr());
    println!(
        "Covariances of zero median data:\n{}",
        pt.covar(&pt.gmedian(EPS)).gr()
    );
    println!(
        "Comediances:\n{}",
        pt.comed(&pt.gmedian(EPS)).gr()
    );
    let trpt = pt.transpose();
    println!(
        "Column Dependencies:\n{}",
        trpt.crossfeatures(|v1, v2| v1.dependence(v2)).gr()
    );
    println!(
        "Column Correlations:\n{}",
        trpt.crossfeatures(|v1, v2| v1.mediancorr(v2)).gr()
    );
}

#[test]
fn ustats() -> Result<(),RError> { 
    set_seeds(1234567);
    let v1 = Rnum::newu8().ranv(20).getvu8(); 
    println!("\n{}", (&v1).gr());
    // println!("Linear transform:\n{}",v1.lintrans()));
    println!("Arithmetic mean:{}", v1.amean()?.gr());
    println!("Geometric mean:\t{}", v1.gmean()?.gr());
    println!("Harmonic mean:\t{}", v1.hmean()?.gr());
    println!("Magnitude:\t{}", v1.vmag());
    println!("Arithmetic {}", v1.ameanstd()?.gr());
    println!("Geometric  {}", v1.gmeanstd()?.gr());
    println!("Harmonic   {}", v1.hmeanstd()?.gr());
    println!("Autocorrelation:{}", v1.autocorr().gr());
    println!("{}\n", v1.as_slice().medinfo());
    Ok(())
}

#[test]
/// &[i64] requires explicit recast
fn intstats() -> Result<(),RError> {
    let v = vec![1_i64, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    println!("\n{}", (&v).gr());
    let v1 = i64tof64(&v); // downcast to f64 here
    println!("Linear transform:\n{}",v1.lintrans().gr());
    println!("Arithmetic mean:{}", v1.amean()?.gr());
    println!("Geometric mean:\t{}", v1.gmean()?.gr());
    println!("Harmonic mean:\t{}", v1.hmean()?.gr());
    // println!("Magnitude:\t{}",v1.vmag()));
    println!("Arithmetic {}", v1.ameanstd()?.gr());
    println!("Geometric  {}", v1.gmeanstd()?.gr());
    println!("Harmonic   {}", v1.hmeanstd()?.gr());
    println!("Autocorrelation:{}", v1.autocorr().gr());
    println!("{}\n", v1.as_slice().medinfo());
    Ok(())
}

#[test]
/// Generic implementation
fn genericstats() -> Result<(),RError> {
    let v = vec![1_i32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    println!("\n{}", (&v).gr());
    println!("Arithmetic\t{}", v.ameanstd()?.gr());
    println!("Geometric\t{}", v.gmeanstd()?.gr());
    println!("Harmonic\t{}", v.hmeanstd()?.gr());
    println!("Weighted Arit.\t{}", v.awmeanstd()?.gr());
    println!("Weighted Geom.\t{}", v.gwmeanstd()?.gr());
    println!("Weighted Harm.\t{}", v.hwmeanstd()?.gr());
    println!("Autocorrelation:{}", v.autocorr().gr());
    println!("{}\n", &v.as_slice().medinfo());
    Ok(())
}

#[test]
fn vecg() {
    let v1 = vec![
        1_f64, 2., 3., 4., 5., 6., 7., 8., 9., 10., 10., 10., 13., 14., 15.,
    ];
    println!("v1: {}", (&v1).gr());
    let v2 = vec![
        1_f64, 14., 2., 13., 3., 12., 4., 11., 5., 10., 6., 6., 7., 1., 15.,
    ];
    println!("v2: {}", (&v2).gr());
    println!("Lexical order v1<v2:\t{}", (v1 < v2).gr());
    println!("Pearson's Correlation:\t{}", v1.correlation(&v2).gr());
    println!("Median Correlation:\t{}", v1.mediancorr(&v2).gr());
    println!("Kendall's Correlation:\t{}", v1.kendalcorr(&v2).gr());
    println!("Spearman's Correlation:\t{}", v1.spearmancorr(&v2).gr());
    println!("Euclidian distance:\t{}", v1.vdist::<f64>(&v2).gr());
    println!("Cityblock distance:\t{}", v1.cityblockd(&v2).gr());
    println!("Vector difference: {}", v1.vsub(&v2).gr());
    println!("Vector sum: {}", v1.vadd(&v2).gr());
    println!("Scalar product:\t\t{}", v1.dotp(&v2).gr());
    println!("Parallelogram area:\t{}", v1.varea(&v2).gr());
    println!("Arc area:\t\t{}", v1.varc(&v2).gr());
    println!("Entropy v1:\t\t{}", v1.entropy().gr());
    println!("Entropy v2:\t\t{}", v2.entropy().gr());
    println!("Joint Entropy:\t\t{}", v1.jointentropy(&v2).gr());
    println!("Dependence:\t\t{}", v1.dependence(&v2).gr());
    println!("Independence:\t\t{}", v1.independence(&v2).gr());
    println!("Cosine:\t\t\t{}", v1.cosine(&v2).gr());
    println!("Cosine of ranks:\t{}",
        v1.rank(true)
            .indx_to_f64()
            .cosine(&v2.rank(true).indx_to_f64())
            .gr() );
    println!("Cos Similarity [0,1]:\t{}", v1.vsim(&v2).gr());
    println!("Cos Dissimilarity:\t{}", v1.vdisim(&v2).gr());
    println!("[1,2,3].kron(&[4,5]):\t{}", [1, 2, 3].kron(&[4, 5]).gr());
    let outerp = [1, 2, 3].outer(&[4, 5]);
    println!("[1,2,3].outer(&[4,5]):\n{}", outerp.gr());
    // println!("Transposed: "); printvv([1,2,3].outer(&[4,5,6,7]).transpose()); 
}

#[test]
/// Trend between two data sets in space of the same dimensions but
/// numbers of points can differ
fn trend() {
    let d = 7_usize;
    set_seeds(777);
    let rf64 = Rnum::newf64(); 
    let pts1 = rf64.ranvv(d, 37).getvvf64();
    let pts2 = rf64.ranvv(d, 33).getvvf64();
    println!("\nTrend vector:\n{}\n", pts1.trend(EPS, pts2).gr());
}

#[test]
fn vecvec() -> Result<(),RError> {
    let d = 10_usize;
    let n = 90_usize;
    println!("Testing on a random set of {} points in {} dimensional space",n,d);
    set_seeds(113);
    let ru = Rnum::newu8();
    let pts = ru.ranvv(d,n).getvvu8(); 
    // println!("\nTest data:\n{}",pts.gr());
    println!("Set joint entropy: {}", pts.jointentropyn().gr());
    println!("Set dependence:    {}", pts.dependencen().gr());
    // println!("\nTest outcomes:\n{}",pts.gr());  
    let outcomes = ru.ranv(n).getvu8();
    let transppt = pts.transpose();
    println!(
        "\nDependencies of outcomes:\n{}",
        transppt.dependencies(&outcomes).gr()
    );
    println!(
        "Correlations with outcomes:\n{}",
        transppt.correlations(&outcomes).gr()
    );
    
    let (gm,_vsum,recips) = pts.gmparts(EPS);
    let (eccstd, eccmed, eccecc) = pts.eccinfo(&gm); 
    let medoid = &pts[eccecc.minindex];
    let outlier = &pts[eccecc.maxindex];
    let hcentroid = pts.hcentroid();
    let gcentroid = pts.gcentroid();
    let acentroid = pts.acentroid();
    let firstp = pts.firstpoint();
    let median = pts.gmedian(EPS);  

    println!("\nMean reciprocal to gm: {}",(recips/d as f64).gr() );
    println!("Tukeyvec for outlier:\n{}",pts.tukeyvec(&outlier.tof64()).gr());    
    println!("Magnitude of Tukey vec for gm: {}",pts.tukeyvec(&median).vmag().gr());
    println!("Mag of Tukeyvec for acentroid: {}",pts.tukeyvec(&acentroid).vmag().gr());
    println!("Mag of Tukeyvec for outlier:   {}",pts.tukeyvec(&outlier.tof64()).vmag().gr());
    // let testvec = ru.ranv(d).getvu8();
    let dists = pts.distsums();
    let md = dists.minmax(); 
    println!("\nMedoid and Outlier Total Distances:\n{}", md);
    println!("Total Distances {}", dists.ameanstd()?.gr());
    println!("Total distances {}", dists.as_slice().medinfo());
    println!("GM's total distances:\t{}", pts.distsum(&median).gr()); 
    println!("ACentroid's total distances:\t{}",pts.distsum(&acentroid).gr());
    println!("HCentroid's total distances:\t{}",pts.distsum(&hcentroid).gr());
    println!("GCentroid's total distances:\t{}",pts.distsum(&gcentroid).gr());
 
    println!(
        "\nMedoid, outlier and radii summary:\n{}\nRadii {}\nRadii {}",
        eccecc, eccstd, eccmed
    );
    let (ngm,_,_) = pts.nxnonmember(&median);
    println!("MADGM: {}", pts.madgm(&median).gr());
    println!("Median's error:      {GR}{:e}{UN}",ngm.vdist::<f64>(&median));
    println!("ACentroid's radius:  {}", acentroid.vdist(&median).gr());
    println!("Firstpoint's radius: {}", firstp.vdist(&median).gr());
    println!("Medoid's radius:     {}",medoid.vdist(&median).gr());
    println!("HCentroid's radius:  {}", hcentroid.vdist(&median).gr());
    println!("GCentroid's radius:  {}", gcentroid.vdist(&median).gr());
    println!("Outlier's radius:    {}",outlier.vdist(&median).gr()); 
    println!("Outlier from Medoid: {}",outlier.vdist(medoid).gr());

    let seccs = pts.sortedeccs(true, &median);
    // println!("\nSorted eccs: {}\n", seccs));
    let lqcnt = seccs.binsearch(eccmed.lq);
    println!(
        "Inner quarter of points: {} within radius: {}",
        lqcnt.gr(),
        seccs[lqcnt - 1].gr()
    );
    let medcnt = pts.len()/2; 
    // seccs.binsearch(eccmed.median);
    println!(
        "Inner half of points:    {} within radius: {}",
        medcnt.gr(),
        seccs[medcnt - 1].gr()
    );
    let uqcnt = seccs.binsearch(eccmed.uq);
    println!(
        "Inner three quarters:    {} within radius: {}",
        uqcnt.gr(),
        seccs[uqcnt - 1].gr()
    );

    println!("\nContribution of adding acentroid:   {}",acentroid.contrib_newpt(&gm,recips).gr() );
    println!("Contribution of adding gcentroid:   {}",gcentroid.contrib_newpt(&gm,recips).gr() );
    println!("Contribution of removing outlier:  {}",outlier.contrib_oldpt(&gm,recips).gr() );
    let contribs = pts.iter().map(|p| p.contrib_oldpt(&gm, recips)).collect::<Vec<f64>>();
    println!("\nContributions of Data Points, Summary:\n{}\n{}\n{}",contribs.minmax(),contribs.ameanstd().unwrap(),contribs.medinfo());
    // create pretend median of medians
    // let medmed = vec![0.5_f64;n];
    // let (se, cpdf) =
    //  pt.wsortedcos(&medmed,&medmed.vunit(), &weights);
    //  println!("Sorted coses:\n{}\ncpdf:\n{}\n",se),cpdf));
    Ok(())
}

const NAMES:[&str;4] = [ "gmedian","pmedian","quasimedian","acentroid" ];

const CLOSURESU8:[fn(&[Vec<f64>]);4] = [
    |v:&[_]| { v.gmedian(EPS); }, 
    |v:&[_]| { v.pmedian(EPS); }, 
    |v:&[_]| { v.quasimedian(); },
    |v:&[_]| { v.acentroid(); } ];

#[test]
fn timing_gm() {
    set_seeds(7777777777_u64);   // intialise random numbers generator
    // Rnum encapsulates the type of the data items
   benchvvf64(Rnum::newf64(),100,5,10,&NAMES,&CLOSURESU8); 
}

#[test]
fn geometric_medians() {
    const ITERATIONS: usize = 10;
    let n = 100_usize;
    let d = 1000_usize;
    set_seeds(7777777);
    println!(
        "{} medians of {} points in {} dimensions", ITERATIONS, n, d
    );
    let mut sumg = 0_f64; 
    let mut sump = 0_f64; 
    let mut sumq = 0_f64; 
    let mut summ = 0_f64; 
    let mut gm: Vec<f64>;
     for _i in 1..ITERATIONS { 
        let pts = Rnum::newf64().ranvv(d, n).getvvf64();
        gm = pts.gmedian(EPS);
        sumg += pts.gmerror(&gm);
        gm = pts.pmedian(EPS);
        sump += pts.gmerror(&gm);    
        gm = pts.quasimedian(); 
        sumq += pts.gmerror(&gm);
        gm = pts.acentroid();
        summ += pts.gmerror(&gm);
    }
    println!( "Pmedian      err: {GR}{:.5e}{UN}", sump);
    println!( "Gmedian      err: {GR}{:.5e}{UN}", sumg);
    println!( "Acentroid    err: {GR}{:.5e}{UN}", summ);
    println!( "Quasi median err: {GR}{:.5e}{UN}", sumq);
}
