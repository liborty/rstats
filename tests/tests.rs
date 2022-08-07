use indxvec::{printing::*, Indices, Vecops,Printing};
use rstats::{RE, i64tof64, identity_lmatrix, Stats, VecVec, VecVecg, Vecg, Vecu8};
use ran::{*,set_seeds};
use medians::{Median};
use times::{benchvvf64};

pub const EPS: f64 = 1e-3;

#[cfg(test)]

#[test]
fn u8() -> Result<(),RE> {
    let v1 = vec![
        1_u8,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,
    ];
    println!("\nv1: {}", (&v1).gr());
    let v2 = vec![
        1_u8,2,2,3,3,3,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,
    ];
    println!("v2: {}", (&v2).gr());
    println!("v1+v2: {}",v1.vadd(&v2).gr());   
    println!("v1-v2: {}",v1.vsub(&v2).gr());    
    println!("Lexical order v1<v2: {}", (v1 < v2).gr());
    println!("v1*v2:\t{}",v1.dotp(&v2).gr());
    println!("Entropy v1:\t{}", v1.entropy().gr());
    println!("Entropyu8 v1:\t{}", v1.entropyu8().gr());
    println!("Entropy v2:\t{}", v2.entropy().gr()); // generic
    println!("Entropyu8 v2:\t{}", v2.entropyu8().gr()); // u8 
    println!("Joint Entropy:  {}", v1.jointentropy(&v2)?.gr());
    println!("Joint Entropyu8:{}", v1.jointentropyu8(&v2)?.gr());
    println!("Dependence:   {}", v1.dependence(&v2)?.gr()); // generic
    println!("Dependenceu8: {}", v1.dependenceu8(&v2)?.gr()); // u8 
    println!("Median v1: {}",v1.medstats());
    println!("Median v2: {}",v2.medstats());
    println!("v1 {}",v1.medinfo());   
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
        trpt.crossfeatures(|v1, v2| v1.dependence(v2))?.gr()
    );
    println!(
        "Column Correlations:\n{}",
        trpt.crossfeatures(|v1, v2| Ok(v1.mediancorr(v2)))?.gr()
    ); 
    Ok(())
}

#[test]
fn fstats() -> Result<(),RE> {
    let v1 = vec![
        1_f64, 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.,
    ];
    println!("\n{}", (&v1).gr());
    let v2 = v1.revs();
    println!("{}", (&v2).gr());
    // println!("Linear transform:\n{}",v1.lintrans()));
    println!("Geometric mean:\t{}", v1.gmean()?.gr());
    println!("Harmonic mean:\t{}", v1.hmean()?.gr());
    println!("Magnitude:\t{}", v1.vmag().gr());
    println!("Arithmetic {}", v1.ameanstd()?);
    println!("Geometric  {}", v1.gmeanstd()?);
    println!("Harmonic   {}", v1.hmeanstd()?);
    println!("Median     {}", v1.medstats());
    println!("Autocorrelation:{}", v1.autocorr()?.gr());
    println!("Entropy 1:\t{}", v1.entropy().gr());
    println!("Entropy 2:\t{}", v2.entropy().gr()); // generic 
    println!("Joint Entropy:  {}", v1.jointentropy(&v2)?.gr());
    println!("Dependence:\t{}", v1.dependence(&v2)?.gr()); // generic
    println!("Euclidean dist:\t{}", v2.vdist(&v1).gr());
    println!("Cityblock dist:\t{}", v2.cityblockd(&v1).gr());
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
        trpt.crossfeatures(|v1, v2| v1.dependence(v2))?.gr()
    );
    println!(
        "Column Correlations:\n{}",
        trpt.crossfeatures(|v1, v2| Ok(v1.mediancorr(v2)))?.gr()
    );
    Ok(())
}

#[test]
fn ustats() -> Result<(),RE> { 
    set_seeds(1234567);
    let v1 = Rnum::newu8().ranv(20).getvu8(); 
    println!("\n{}", (&v1).gr());
    // println!("Linear transform:\n{}",v1.lintrans()));
    println!("Arithmetic mean: {GR}{:>14.10}{UN}", v1.amean()?);
    println!("Median:          {GR}{:>14.10}{UN}", v1.median());   
    println!("Geometric mean:  {GR}{:>14.10}{UN}", v1.gmean()?);
    println!("Harmonic mean:   {GR}{:>14.10}{UN}", v1.hmean()?);
    println!("Magnitude:       {GR}{:>14.10}{UN}", v1.vmag());
    println!("Arithmetic {}", v1.ameanstd()?);
    println!("Median     {}", v1.medstats());
    println!("Geometric  {}", v1.gmeanstd()?);
    println!("Harmonic   {}", v1.hmeanstd()?);
    println!("Autocorrelation:{}", v1.autocorr()?.gr());
    println!("{}\n", v1.medinfo());
    Ok(())
}

#[test]
/// &[i64] requires explicit recast
fn intstats() -> Result<(),RE> {
    let v = vec![1_i64, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    println!("\n{}", (&v).gr());
    let v1 = i64tof64(&v); // downcast to f64 here
    println!("Linear transform:\n{}",v1.lintrans()?.gr());
    println!("Arithmetic mean:{}",v1.amean()?.gr());
    println!("Median:       {GR}{:>14.10}{UN}",v1.median());   
    println!("Geometric mean:\t{}", v1.gmean()?.gr());
    println!("Harmonic mean:\t{}", v1.hmean()?.gr());
    // println!("Magnitude:\t{}",v1.vmag()));
    println!("Arithmetic {}", v1.ameanstd()?);
    println!("Median     {}", v1.medstats());
    println!("Geometric  {}", v1.gmeanstd()?);
    println!("Harmonic   {}", v1.hmeanstd()?);
    println!("Autocorrelation:{}", v1.autocorr()?.gr());
    println!("{}\n", v1.medinfo());
    Ok(())
}

#[test]
/// Generic implementation
fn genericstats() -> Result<(),RE> {
    let v = vec![1_i32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    println!("\n{}", (&v).gr());
    println!("Arithmetic\t{}", v.ameanstd()?.gr());
    println!("Median\t\t{GR}{:>14.10}{UN}",v.medstats());
    println!("Geometric\t{}", v.gmeanstd()?.gr());
    println!("Harmonic\t{}", v.hmeanstd()?.gr());
    println!("Weighted Arit.\t{}", v.awmeanstd()?.gr());
    println!("Weighted Geom.\t{}", v.gwmeanstd()?.gr());
    println!("Weighted Harm.\t{}", v.hwmeanstd()?.gr());
    println!("Autocorrelation:{}", v.autocorr()?.gr());
    println!("{}\n", &v.medinfo());
    Ok(())
}

#[test]
fn vecg() -> Result<(),RE> {
    let v1 = vec![
        1_f64, 2., 3., 4., 5., 6., 7., 8., 9., 10., 10., 10., 13., 14., 15.,
    ];
    println!("v1: {}", (&v1).gr());
    let v2 = vec![
        1_f64, 14., 2., 13., 3., 12., 4., 11., 5., 10., 6., 6., 7., 1., 15.,
    ];
    println!("v2: {}", (&v2).gr());
    println!("Lexical order v1<v2:\t{}", (v1 < v2).gr());
    println!("Median Correlation:\t{}", v1.mediancorr(&v2).gr());
    println!("Pearson's Correlation:\t{}", v1.correlation(&v2).gr());
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
    println!("Joint Entropy:\t\t{}", v1.jointentropy(&v2)?.gr());
    println!("Dependence:\t\t{}", v1.dependence(&v2)?.gr());
    println!("Independence:\t\t{}", v1.independence(&v2)?.gr());
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
    Ok(()) 
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
fn matrices() -> Result<(),RE> {
    println!("5x5 identity lower triangular matrix in scan order:\n{}",identity_lmatrix(5).gr());
    let d = 10_usize;
    let n = 90_usize;
    println!("Testing on a random set of {} points in {} dimensional space",n,d);
    set_seeds(1133);
    let ru = Rnum::newf64();
    let pts = ru.ranvv_in(d,n,0.0,4.0).getvvf64(); 
    // println!("\nTest data:\n{}",pts.gr());
    // let transppt = pts.transpose();  
    let cov = pts.covar(&pts.pmedian(EPS));
    println!("Triangular comediance matrix in scan order, size {}:\n{}",cov.len(),cov.gr());
    let cholmat = cov.cholesky()?; 
    println!("\nCholesky L matrix, size {}:\n{}",cholmat.len(),cholmat.gr());
    // set_seeds(77777);
    let pta = ru.ranv(d).getvf64();
    let ptb = ru.ranv(d).getvf64();
    let d = pta.vsub(&ptb);
    let dmag = d.vmag();
    let mahamag =cholmat.mahalanobis(&d)?;
    println!("Test vector d = a-b:\n{}",d.gr());
    println!("Euclidian magnitude   {GR}{:>8.4}{UN}\
        \nMahalanobis magnitude {GR}{:>8.4}{UN}\
        \nScale factor: {GR}{:>0.8}{UN}",
        dmag,mahamag,mahamag/dmag);
    Ok(())
}

#[test]
fn vecvec() -> Result<(),RE> {
    let d = 10_usize;
    let n = 90_usize;
    println!("Testing on a random set of {} points in {} dimensional space",n,d);
    set_seeds(113);
    let ru = Rnum::newu8();
    let pts = ru.ranvv_in(d,n,0.,4.).getvvu8(); 
    // println!("\nTest data:\n{}",pts.gr());
    println!("Set joint entropy: {}", pts.jointentropyn()?.gr());
    println!("Set dependence:    {}", pts.dependencen()?.gr());
    // println!("\nTest outcomes:\n{}",pts.gr());  
    let outcomes = ru.ranv(n).getvu8();
    let transppt = pts.transpose();
    println!(
        "\nDependencies of outcomes:\n{}",
        transppt.dependencies(&outcomes)?.gr()
    );
    println!(
        "Correlations with outcomes:\n{}",
        transppt.correlations(&outcomes).gr()
    );    
    let (median,_vsum,recips) = pts.gmparts(EPS);
    let (eccstd, eccmed, eccecc) = pts.eccinfo(&median); 
    let medoid = &pts[eccecc.minindex];
    let outlier = &pts[eccecc.maxindex];
    let hcentroid = pts.hcentroid();
    let gcentroid = pts.gcentroid();
    let acentroid = pts.acentroid();
    let firstp = pts.firstpoint(); 
  
    println!("\nMean reciprocal of radius: {}",(recips/d as f64).gr() );
    println!("Tukeyvec for outlier:\n{}",pts.tukeyvec(&outlier.tof64()).gr());    
    println!("Magnitude of Tukey vec for gm: {}",pts.tukeyvec(&median).vmag().gr());
    println!("Mag of Tukeyvec for acentroid: {}",pts.tukeyvec(&acentroid).vmag().gr());
    println!("Mag of Tukeyvec for outlier:   {}",pts.tukeyvec(&outlier.tof64()).vmag().gr());
    // let testvec = ru.ranv(d).getvu8();
    let dists = pts.distsums();
    let md = dists.minmax(); 
    println!("\nMedoid and Outlier Total Distances:\n{}", md);
    println!("Total Distances {}", dists.ameanstd()?);
    println!("Total distances {}", dists.medinfo());
    println!("GM's total distances:        {}", pts.distsum(&median).gr()); 
    println!("ACentroid's total distances: {}",pts.distsum(&acentroid).gr());
    println!("HCentroid's total distances: {}",pts.distsum(&hcentroid).gr());
    println!("GCentroid's total distances: {}",pts.distsum(&gcentroid).gr());
 
    println!(
        "\nMedoid, outlier and radii summary:\n{}\nRadii {}\nRadii {}",
        eccecc, eccstd, eccmed
    );
    println!("Madgm:               {}", pts.madgm(&median).gr());
    println!("Median's error:      {GR}{:e}{UN}",pts.gmerror(&median));
    println!("ACentroid's radius:  {}", acentroid.vdist(&median).gr());
    println!("Firstpoint's radius: {}", firstp.vdist(&median).gr());
    println!("Medoid's radius:     {}",medoid.vdist(&median).gr());
    println!("HCentroid's radius:  {}", hcentroid.vdist(&median).gr());
    println!("GCentroid's radius:  {}", gcentroid.vdist(&median).gr());
    println!("Outlier's radius:    {}",outlier.vdist(&median).gr()); 
    println!("Outlier to Medoid:   {}",outlier.vdist(medoid).gr());

    let seccs = pts.radii(&median).sorth(true);
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

    println!("\nContribution of adding acentroid:   {}",acentroid.contrib_newpt(&median,recips).gr() );
    println!("Contribution of adding gcentroid:   {}",gcentroid.contrib_newpt(&median,recips).gr() );
    println!("Contribution of removing outlier:  {}",outlier.contrib_oldpt(&median,recips).gr() );
    let contribs = pts.iter().map(|p| p.contrib_oldpt(&median,recips)).collect::<Vec<f64>>();
    println!("\nContributions of Data Points, Summary:\n{}\n{}\n{}",contribs.minmax(),contribs.ameanstd().unwrap(),contribs.medinfo());
    
    let weights:Vec<usize> = Vec::from_iter(0..n);
    let fweights = weights.indx_to_f64(); 
    let wmedian = pts.wgmedian(&fweights,EPS);
    println!("Sorted weighted rads:\n{}",pts.wsortedrads(&fweights,&wmedian).gr());
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
   benchvvf64(Rnum::newf64(),100,4,10,&NAMES,&CLOSURESU8); 
}

#[test]
fn geometric_medians() {
    const ITERATIONS: usize = 10;
    let n = 100_usize;
    let d = 1000_usize;
    set_seeds(7777777);
    println!(
        "{} repeats of {} points in {} dimensions", ITERATIONS, n, d
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
    println!( "Pmedian      total error: {GR}{:.5e}{UN}", sump);
    println!( "Gmedian      total error: {GR}{:.5e}{UN}", sumg);
    println!( "Acentroid    total error: {GR}{:.5e}{UN}", summ);
    println!( "Quasi median total errot: {GR}{:.5e}{UN}", sumq);
}
