use indxvec::{printing::*, Indices, Printing, Vecops};
use medians::{Median,Medianf64};
use ran::{set_seeds, Rnum };
use rstats::{st_error, noop, fromop, unit_matrix, Stats, TriangMat, VecVec, VecVecg, Vecg, Vecu8, RE};
use times::benchvvf64;

pub const EPS: f64 = 1e-3;

#[cfg(test)]
#[test]
fn u8() -> Result<(), RE> {
    let v1 = vec![
        1_u8, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
    ];
    println!("\nv1: {}", (&v1).gr());
    let v2 = vec![
        1_u8, 2, 2, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2,
    ];
    println!("v2: {}", (&v2).gr());
    println!("v1+v2: {}", v1.vadd(&v2).gr());
    println!("v1-v2: {}", v1.vsub(&v2).gr());
    println!("Lexical order v1<v2: {}", (v1 < v2).gr());
    println!("v1*v2:\t{}", v1.dotp(&v2).gr());
    println!("Entropy v1:\t{}", v1.entropy().gr());
    println!("Entropyu8 v1:\t{}", v1.entropyu8().gr());
    println!("Entropy v2:\t{}", v2.entropy().gr()); // generic
    println!("Entropyu8 v2:\t{}", v2.entropyu8().gr()); // u8
    println!("Joint Entropy:  {}", v1.jointentropy(&v2)?.gr());
    println!("Joint Entropyu8:{}", v1.jointentropyu8(&v2)?.gr());
    println!("Dependence:   {}", v1.dependence(&v2)?.gr()); // generic
    println!("Dependenceu8: {}", v1.dependenceu8(&v2)?.gr()); // u8
    println!("Median v1: {}", v1.medstats(&mut |&x| x.into())?);
    println!("Median v2: {}", v2.medstats(&mut fromop)?);
    let d = 5_usize;
    let n = 7_usize;
    println!("Testing on a random set of {}points in {}d space:", n.yl(), d.yl());
    set_seeds(77777);
    let pt = Rnum::newu8().ranvv(d, n)?.getvvu8()?;
    let cov = pt.covar(&pt.acentroid())?;
    println!("Covariances:\n{cov}");
    let com = pt.covar(&pt.gmedian(EPS))?;
    println!("Their Distance: {}", cov.data.vdist(&com.data)); 
    println!(
        "Column Median Correlations:\n{}",
        pt.transpose().crossfeatures( |v1, v2| 
            v1.mediancorrf64(v2).expect("median corr: crossfeatures u8\n"))?
    );
    Ok(())
}

#[test]
fn fstats() -> Result<(), RE> {
    let v1 = vec![
        1_f64, 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.,
    ];
    println!("\n{}", (&v1).gr());
    let v2 = v1.revs();
    println!("{}", (&v2).gr());
    // println!("Linear transform:\n{}",v1.lintrans()));
    println!("Standard arithmetic error of {}: {}",5.yl(),st_error(5.,v1.ameanstd()?).gr());
    println!("Standard harmonic error of {}:\t {}",5.yl(),st_error(5.,v1.hmeanstd()?).gr());
    println!("Standard median error of {}:\t{}",5.yl(),st_error(5.,v1.medstatsf64()?).gr());
    println!("Geometric mean:\t{}", v1.gmean()?.gr());
    println!("Harmonic mean:\t{}", v1.hmean()?.gr());
    println!("Magnitude:\t{}", v1.vmag().gr());
    println!("Arithmetic {}", v1.ameanstd()?);
    println!("Geometric  {}", v1.gmeanstd()?);
    println!("Harmonic   {}", v1.hmeanstd()?);
    println!("Median     {}", v1.medstatsf64()?);
    println!("Autocorrelation:{}", v1.autocorr()?.gr());
    println!("Entropy 1:\t{}", v1.entropy().gr());
    println!("Entropy 2:\t{}", v2.entropy().gr()); // generic
    println!("Joint Entropy:  {}", v1.jointentropy(&v2)?.gr());
    println!("Dependence:\t{}", v1.dependence(&v2)?.gr()); // generic
    println!("Euclidean dist:\t{}", v2.vdist(&v1).gr());
    println!("Cityblock dist:\t{}", v2.cityblockd(&v1).gr());
    let d = 5_usize;
    let n = 7_usize;
    println!("Testing on a random set of {n} points in {d} d space:");
    let pt = Rnum::newf64().ranvv(d, n)?.getvvf64()?;
    println!(
        "Classical Covariances:\n{}",
        pt.covar(&pt.acentroid())?.gr()
    );
    println!(
        "Comediances (covariances of zero median data):\n{}",
        pt.covar(&pt.gmedian(EPS))?.gr()
    ); 
    println!(
        "Column Median Correlations:\n{}",
        pt.transpose().crossfeatures(|v1,v2| v1.mediancorr(v2,&mut noop).expect("median corr: crossfeatures f64\n"))?
    );
    Ok(())
}

#[test]
fn ustats() -> Result<(), RE> {
    set_seeds(1234567);
    let v1 = Rnum::newu8().ranv(20)?.getvu8()?;
    println!("\n{}", (&v1).gr());
    // println!("Linear transform:\n{}",v1.lintrans()));
    println!("Arithmetic mean: {GR}{:>14.10}{UN}", v1.amean()?);
    println!("Median:          {GR}{:>14.10}{UN}", v1.median(&mut fromop)?);
    println!("Geometric mean:  {GR}{:>14.10}{UN}", v1.gmean()?);
    println!("Harmonic mean:   {GR}{:>14.10}{UN}", v1.hmean()?);
    println!("Magnitude:       {GR}{:>14.10}{UN}", v1.vmag());
    println!("Arithmetic {}", v1.ameanstd()?);
    println!("Median     {}", v1.medstats(&mut fromop)?);
    println!("Geometric  {}", v1.gmeanstd()?);
    println!("Harmonic   {}", v1.hmeanstd()?);
    println!("Autocorrelation:{}", v1.autocorr()?.gr());
    Ok(())
}

#[test]
/// &[i64] requires explicit recast
fn intstats() -> Result<(), RE> {
    let v = vec![1_i64, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    println!("\n{}", (&v).gr());
    let v1:Vec<f64> = v.iter().map(|i| *i as f64).collect(); // downcast to f64 here
    println!("Linear transform:\n{}", v1.lintrans()?.gr());
    println!("Arithmetic mean:{}", v1.amean()?.gr());
    println!("Median:       {GR}{:>14.10}{UN}", v1.medianf64()?);
    println!("Geometric mean:\t{}", v1.gmean()?.gr());
    println!("Harmonic mean:\t{}", v1.hmean()?.gr());
    // println!("Magnitude:\t{}",v1.vmag()));
    println!("Arithmetic {}", v1.ameanstd()?);
    println!("Median     {}", v1.medstats(&mut noop)?);
    println!("Geometric  {}", v1.gmeanstd()?);
    println!("Harmonic   {}", v1.hmeanstd()?);
    println!("Autocorrelation:{}", v1.autocorr()?.gr()); 
    Ok(())
}

#[test]
/// Generic implementation
fn genericstats() -> Result<(), RE> {
    let mut v = vec![1_i32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    println!("\n{}", (&v).gr());
    println!("Arithmetic\t{}", v.ameanstd()?);
    println!("Median\t\t{}", v.medstats(&mut fromop)?);
    println!("Geometric\t{}", v.gmeanstd()?);
    println!("Harmonic\t{}", v.hmeanstd()?);
    println!("Weighted Arit.\t{}", v.awmeanstd()?);
    println!("Weighted Geom.\t{}", v.gwmeanstd()?);
    println!("Weighted Harm.\t{}", v.hwmeanstd()?);
    println!("Autocorrelation: {}", v.autocorr()?.gr());
    println!("dfdt:\t\t {}", v.dfdt()?.gr());
    v.reverse();
    println!("rev dfdt:\t{}", v.dfdt()?.gr());
    Ok(())
}

#[test]
fn vecg() -> Result<(), RE> {
    let v1 = vec![
        1_f64, 2., 3., 4., 5., 6., 7., 8., 9., 10., 10., 10., 13., 14., 15.,
    ];
    println!("v1: {}", (&v1).gr());
    let v2 = vec![
        1_f64, 14., 2., 13., 3., 12., 4., 11., 5., 10., 6., 6., 7., 1., 15.,
    ];
    println!("v2: {}", (&v2).gr());
    println!("Lexical order v1<v2:\t{}", (v1 < v2).gr());
    println!("Median Correlation:\t{}", v1.mediancorr(&v2,&mut noop)?.gr());
    println!("Pearson's Correlation:\t{}", v1.correlation(&v2).gr());
    println!("Kendall's Correlation:\t{}", v1.kendalcorr(&v2).gr());
    println!("Spearman's Correlation:\t{}", v1.spearmancorr(&v2).gr());
    println!("Euclidian distance:\t{}", v1.vdist(&v2).gr());
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
    println!(
        "Cosine of ranks:\t{}",
        v1.rank(true)
            .indx_to_f64()
            .cosine(&v2.rank(true).indx_to_f64())
            .gr()
    );
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
fn trend() -> Result<(), RE> {
    let d = 7_usize;
    set_seeds(777);
    let rf64 = Rnum::newf64();
    let pts1 = rf64.ranvv(d, 37)?.getvvf64()?;
    let pts2 = rf64.ranvv(d, 50)?.getvvf64()?;
    println!("\nTrend vector:\n{}\n", pts1.trend(EPS, pts2)?.gr());
    Ok(())
}

#[test]
fn triangmat() -> Result<(), RE> {
    println!("\n{}", TriangMat::unit(5).gr());
    println!("{}", TriangMat::unit(7).gr());
    let d = 10_usize;
    let n = 90_usize;
    println!("Testing on a random set of {n} points in {d} dimensional space");
    set_seeds(1133);
    let ru = Rnum::newf64();
    let pts = ru.ranvv_in(d, n, 0.0, 4.0)?.getvvf64()?;
    // println!("\nTest data:\n{}",pts.gr());
    // let transppt = pts.transpose();
    let cov = pts.covar(&pts.par_gmedian(EPS))?;
    println!("Comediance matrix:\n{cov}");
    let chol = cov.cholesky()?;
    println!("Cholesky L matrix:\n{chol}"); 
    // chol.transpose();   
    // println!("L matrix trivially transposed:\n{}",chol.to_full().gr()); 
    let pta = ru.ranv(d)?.getvf64()?;
    let ptb = ru.ranv(d)?.getvf64()?;
    let d = pta.vsub(&ptb);
    let dmag = d.vmag();
    let mahamag = chol.mahalanobis(&d)?;
    println!("Test vector d = a-b:\n{}", d.gr());
    println!(
        "Euclidian magnitude   {GR}{:>8.4}{UN}\
        \nMahalanobis magnitude {GR}{:>8.4}{UN}\
        \nScale factor: {GR}{:>0.8}{UN}",
        dmag,
        mahamag,
        mahamag / dmag
    );
    Ok(())
}

#[test]
fn mat() -> Result<(), RE> { 
    let d = 10_usize;
    let n = 12_usize;
    println!("Testing on a random set of {n} points in {d} dimensional space");
    set_seeds(1133);
    let ru = Rnum::newf64();
    let m = ru.ranvv(d, n)?.getvvf64()?;
    println!("\nTest matrix M:\n{}",m.gr());
    let t = m.transpose();
    println!("\nTransposed matrix T:\n{}",t.gr());
    let v = ru.ranv(d)?.getvf64()?;
    println!("\nVector V:\n{}",v.gr()); 
    println!("\nMV:\n{}",m.leftmultv(&v)?.gr());  
    println!("\nVT:\n{}",t.rightmultv(&v)?.gr()); 
    println!("\nMT:\n{}",t.matmult(&m)?.gr());
    println!("\nTM:\n{}",m.matmult(&t)?.gr());                
    Ok(())
}


#[test]
fn vecvec() -> Result<(), RE> {
    let d = 10_usize;
    let n = 90_usize;
    println!("Testing on a random set of {n} points in {d} dimensional space");
    set_seeds(113);
    let ru = Rnum::newu8();
    let pts = ru.ranvv_in(d, n, 0., 4.)?.getvvu8()?;
    // println!("\nTest data:\n{}",pts.gr());
    println!("Set joint entropy: {}", pts.jointentropyn()?.gr());
    println!("Set dependence:    {}", pts.dependencen()?.gr());
    // println!("\nTest outcomes:\n{}",pts.gr());
    let outcomes = ru.ranv(n)?.getvu8()?;
    let transppt = pts.transpose();
    println!(
        "\nDependencies of outcomes:\n{}",
        transppt.dependencies(&outcomes)?.gr()
    );
    println!(
        "Correlations with outcomes:\n{}",
        transppt.correlations(&outcomes)?.gr()
    );
    let (median, _vsum, recips) = pts.gmparts(EPS);    
    let (eccstd, eccmed, eccecc) = pts.eccinfo(&median[..])?;
    let medoid = &pts[eccecc.minindex];
    println!("Medoid: {}",medoid.gr());
    let outlier = &pts[eccecc.maxindex];
    let hcentroid = pts.hcentroid();
    let gcentroid = pts.gcentroid();
    let acentroid = pts.acentroid();
    let firstp = pts.firstpoint();
    let idx = Vec::from_iter(0..n);

    println!("\nMean reciprocal of radius: {}", (recips / d as f64).gr());

    println!(
        "Magnitude of Tukey vec for gm: {}",
        pts.translate(&median)?.tukeyvec(&idx)?.vmag().gr()
    );
    println!(
        "Mag of Tukeyvec for acentroid: {}",
        pts.translate(&acentroid)?.tukeyvec(&idx)?.vmag().gr()
    );
    println!(
        "Mag of Tukeyvec for outlier:   {}",
        pts.translate(outlier)?.tukeyvec(&idx)?.vmag().gr()
    );
    // let testvec = ru.ranv(d).getvu8();
    let dists = pts.distsums();
    let md = dists.minmax();
    println!("\nMedoid and Outlier Total Distances:\n{md}");
    println!("Total Distances {}", dists.ameanstd()?);
    println!("Total distances {}", dists.medstatsf64()?);
    println!(
        "GM's total distances:        {}",
        pts.distsum(&median)?.gr()
    );
    println!(
        "ACentroid's total distances: {}",
        pts.distsum(&acentroid)?.gr()
    );
    println!(
        "HCentroid's total distances: {}",
        pts.distsum(&hcentroid)?.gr()
    );
    println!(
        "GCentroid's total distances: {}",
        pts.distsum(&gcentroid)?.gr()
    );

    println!(
        "\nMedoid, outlier and radii summary:\n{eccecc}\nRadii {eccstd}\nRadii {eccmed}");
    let radsindex = pts.radii(&median).hashsort_indexed(&mut |x| *x);
    println!(
        "Radii ratio:\t {GR}{}{UN}",
        pts.radius(radsindex[0], &median)? / pts.radius(radsindex[radsindex.len() - 1], &median)?
    );
    println!("Madgm:               {}", pts.madgm(&median)?.gr());
    println!("Median's error:      {GR}{:e}{UN}", pts.gmerror(&median));
    println!("ACentroid's radius:  {}", acentroid.vdist(&median).gr());
    println!("Firstpoint's radius: {}", firstp.vdist(&median).gr());
    println!("Medoid's radius:     {}", medoid.vdist(&median).gr());
    println!("HCentroid's radius:  {}", hcentroid.vdist(&median).gr());
    println!("GCentroid's radius:  {}", gcentroid.vdist(&median).gr());
    println!("Outlier's radius:    {}", outlier.vdist(&median).gr());
    println!("Outlier to Medoid:   {}", outlier.vdist(medoid).gr());

    let seccs = pts.radii(&median).sorth(&mut noop,true);
    // println!("\nSorted eccs: {}\n", seccs));
    let lqcnt = seccs.binsearch(&(eccmed.centre-eccmed.dispersion));
    println!(
        "Inner quarter of points: {} within radius: {}",
        lqcnt.start.gr(),
        seccs[lqcnt.start - 1].gr()
    );
    let medcnt = pts.len() / 2;
    // seccs.binsearch(eccmed.median);
    println!(
        "Inner half of points:    {} within radius: {}",
        medcnt.gr(),
        seccs[medcnt - 1].gr()
    );
    let uqcnt = seccs.binsearch(&(eccmed.centre+eccmed.dispersion));
    println!(
        "Inner three quarters:    {} within radius: {}",
        uqcnt.start.gr(),
        seccs[uqcnt.start - 1].gr()
    );

    println!(
        "\nContribution of adding acentroid:   {}",
        acentroid.contrib_newpt(&median, recips).gr()
    );
    println!(
        "Contribution of adding gcentroid:   {}",
        gcentroid.contrib_newpt(&median, recips).gr()
    );
    println!(
        "Contribution of removing outlier:  {}",
        outlier.contrib_oldpt(&median, recips).gr()
    );
    let contribs = pts
        .iter()
        .map(|p| p.contrib_oldpt(&median, recips))
        .collect::<Vec<f64>>();
    println!(
        "\nContributions of Data Points, Summary:\n{}\n{}\n{}",
        contribs.minmax(),
        contribs.ameanstd()?,
        contribs.medstatsf64()?
    );
    Ok(())
}

#[test]
fn hulls() -> Result<(), RE> {
    let d = 5_usize;
    let n = 100_usize;
    println!(
        "Testing on a random set of {n} points in {d} dimensional space");
    // set_seeds(113);
    let rf = Rnum::newf64();
    let pts = rf.ranvv(d, n)?.getvvf64()?;
    // let wts = rf.ranv_in(n, 0., 100.).getvf64()?;
    let median = pts.gmedian(EPS);
    let zeropts = pts.translate(&median)?;
    let (innerhull,outerhull) = zeropts.hulls();
    let mad = pts.madgm(&median)?;
    println!("\nMADGM: {}",mad.gr());

    println!(
        "\nInner hull has {}/{} points:\n{}",
        innerhull.len().gr(),
        pts.len().gr(),
        innerhull.yl()
    );
    println!(
        "Inner hull min max radii: {} {}\nStandard errors:\t  {} {}",
        zeropts[*innerhull.first().expect("Empty hullidx")].vmag().gr(), 
        zeropts[*innerhull.last().expect("Empty hullidx")].vmag().gr(),
        (pts[*innerhull.first().unwrap()].vdist(&median)/mad).gr(),
        pts[*innerhull.last().unwrap()].st_error(&median,mad)?.gr(),
    );

    println!(
        "\nOuter hull has {}/{} points:\n{}",
        outerhull.len().gr(),
        pts.len().gr(),
        outerhull.yl()
    );
    println!(
        "Outer hull min max radii: {} {}\nStandard errors:\t  {} {}",
        zeropts[*outerhull.first().expect("Empty hullidx")].vmag().gr(),
        zeropts[*outerhull.last().expect("Empty hullidx")].vmag().gr(),
        pts[*outerhull.first().unwrap()].st_error(&median,mad)?.gr(),
        pts[*outerhull.last().unwrap()].st_error(&median,mad)?.gr(),
    );
    let tukeyvec = zeropts
        .tukeyvec(&innerhull)?;
    println!("\nInner hull tukeyvec: {}", tukeyvec.gr());
    println!(
        "Dottukey mapped: {}",
        innerhull
            .iter()
            .map(|&hi| pts[hi].vsub(&median).dottukey(&tukeyvec))
            .collect::<Result<Vec<f64>,RE>>()?
            .gr()
    );
    let allptstukv = zeropts.tukeyvec(&Vec::from_iter(0..zeropts.len()))?;
    println!("All points tukeyvec: {}", allptstukv.gr());
    println!(
        "Dottukey mapped: {}",
        zeropts
            .iter()
            .map(|pt| pt.dottukey(&allptstukv))
            .collect::<Result<Vec<f64>,RE>>()?
            .gr()
    );

    Ok(())
}

#[test]
fn householder() -> Result<(), RE> {
    let a = &[
        vec![35., 1., 6., 26., 19., 24.],
        vec![3., 32., 7., 21., 23., 25.],
        vec![31., 9., 2., 22., 27., 20.],
        vec![8., 28., 33., 17., 10., 15.],
        vec![30., 5., 34., 12., 14., 16.],
        vec![4., 36., 29., 13., 18., 11.],
    ];
    let atimesunit = a.matmult(&unit_matrix(a.len()))?;
    println!("Matrix a:\n{}", atimesunit.gr());
    let (u, r) = a.house_ur(); 
    println!("house_ur u' {u}");
    println!("house_ur r'  {r}");
    let q = u.house_uapply(&unit_matrix(a.len().min(a[0].len())));
    println!(
        "Q matrix\n{}\nOthogonality of Q check (Q'*Q = I):\n{}",
        q.gr(),
        q.transpose().matmult(&q)?.gr()
    );
    println!("Matrix a = QR recreated:\n{}", q.matmult(&r.to_full())?.gr());
    Ok(())
}

#[test]
fn geometric_medians() -> Result<(),RE> {
    const NAMES: [&str; 4] = ["par_gmedian", "gmedian", "quasimedian", "acentroid"];
    const CLOSURESU8: [fn(&[Vec<f64>]); 4] = [
        |v: &[_]| {
            v.par_gmedian(EPS);
        },
        |v: &[_]| {
            v.gmedian(EPS);
        },
        |v: &[_]| {
            v.quasimedian().expect("quasimedian failed");
        },
        |v: &[_]| {
            v.acentroid();
        },
    ];
    set_seeds(7777777777_u64); // intialise random numbers generator
                               // Rnum specifies the type of the random numbers required
    println!("\n{YL}Timing Comparisons{UN}");
    benchvvf64(
        Rnum::newf64(),
        100,
        1000..1500,
        200,
        10,
        &NAMES,
        &CLOSURESU8,
    ); 
    const ITERATIONS: usize = 10;
    let n = 100_usize;
    let d = 1000_usize;
    set_seeds(7777777);  
    println!("\n{YL}Total errors for {ITERATIONS} repeats of {n} points in {d} dimensions{UN}\n");
    let mut sumg = 0_f64; 
    let mut sumr = 0_f64;
    let mut sumq = 0_f64;
    let mut summ = 0_f64;
    let mut gm: Vec<f64>;
    for _i in 1..ITERATIONS {
        let pts = Rnum::newf64().ranvv(d, n)?.getvvf64()?;
        gm = pts.gmedian(EPS);
        sumg += pts.gmerror(&gm); 
        gm = pts.par_gmedian(EPS);
        sumr += pts.gmerror(&gm);
        gm = pts.quasimedian()?;
        sumq += pts.gmerror(&gm);
        gm = pts.acentroid();
        summ += pts.gmerror(&gm);
    }
    println!("{MG}par_gmedian  {GR}{sumr:.10}{UN}");
    println!("{MG}gmedian      {GR}{sumg:.10}{UN}");
    println!("{MG}acentroid    {GR}{summ:.10}{UN}");
    println!("{MG}quasimedian  {GR}{sumq:.10}{UN}\n");
    Ok(())
}
