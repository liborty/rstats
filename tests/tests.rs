use indxvec::{printing::*, Indices, Printing, Vecops};
use medians::{Median, Medianf64};
use ran::*;
use rstats::*;
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
    println!("Median v1: {}", v1.qmedian_by(&mut |a:&u8,b| a.cmp(b), fromop)?);
    println!("Median v2: {}", v2.qmedian_by(&mut |a:&u8,b| a.cmp(b), fromop)?);
    let d = 5_usize;
    let n = 7_usize;
    println!(
        "{YL}Testing on a random set of {}points in {}d space:{UN}",
        n.yl(),
        d.yl()
    );
    set_seeds(77777);
    let pt = ranvv_u8(n,d)?;
    println!("Acentroid:\n{}", pt.acentroid().gr());
    println!("Geometric median :\n{}", pt.gmedian(EPS).gr());
    let cov = pt.covar(&pt.acentroid())?;
    println!("Covariances:\n{cov}");
    let com = pt.covar(&pt.gmedian(EPS))?;
    println!("Comediances:\n{com}");
    println!("Their Distance: {}", cov.data.vdist(&com.data));
    println!(
        "Median correlations of data columns:\n{}",
        pt.transpose()
            .crossfeatures(|v1, v2| v1.medf_correlation(v2).expect("median corr: crossfeatures u8\n"))?
    );
    Ok(())
}

#[test]
fn fstats() -> Result<(), RE> {
    let v1 = vec![
        1_f64, 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 50., 52.,
    ];
    println!("\n{}", (&v1).gr());
    let v2 = v1.revs();
    println!("{}", (&v2).gr());
    println!("Reciprocals of v1:\n{}", v1.vreciprocal()?.gr());
    println!("Unit v1:\n{}", v1.vunit()?.gr());
    println!("Inverse magnitude v1:\n{}", v1.vinverse()?.gr());
    println!("Linear transform of v1:\n{}\n", v1.lintrans()?.gr());
    println!("Magnitudes: {} {}", v1.vmag().gr(), v2.vmag().gr());
    println!("Harmonic spread  {}", v1.hmad()?.gr());
    println!("Arithmetic Mean  {}", v1.ameanstd()?);
    println!("Median & Mad     {}", v1.medmad()?);
    println!("Geometric  Mean  {}", v1.gmeanstd()?);
    println!("Harmonic   Mean  {}", v1.hmeanstd()?);
    println!(
        "tm_stat of 5 against median {}",
        tm_stat(5., v1.medmad()?).gr()
    );
    println!(
        "tm_stat of 5 against amean  {}",
        tm_stat(5., v1.ameanstd()?).gr()
    );
    println!(
        "tm_stat of 5 against gmean  {}",
        tm_stat(5., v1.gmeanstd()?).gr()
    );
    println!(
        "tm_stat of 5 against hmean   {}",
        tm_stat(5., v1.hmeanstd()?).gr()
    );
    println!("Autocorr1:\t{}", v1.autocorr()?.gr());
    println!("Autocorr2:\t{}", v2.autocorr()?.gr());
    println!("Entropy 1:\t{}", v1.entropy().gr());
    println!("Entropy 2:\t{}", v2.entropy().gr()); // generic
    println!("Joint Entropy:  {}", v1.jointentropy(&v2)?.gr());
    println!("Dependence:\t{}", v1.dependence(&v2)?.gr()); // generic
    println!("Euclidean dist:\t{}", v2.vdist(&v1).gr());
    println!("Cityblock dist:\t{}", v2.cityblockd(&v1).gr());
    let d = 5_usize;
    let n = 9_usize;
    println!("{YL}Testing on a random set of {n} points in {d} d space:{UN}");
    let pt = ranvv_f64(n,d)?;
    println!(
        "Classical Covariances (multithreading implementation):\n{}",
        pt.covar(&pt.acentroid())?.gr()
    );
    println!(
        "Classical Covariances (serial implementation):\n{}",
        pt.serial_covar(&pt.acentroid())?.gr()
    );
    println!(
        "Comediances (covariances of zero median data):\n{}",
        pt.covar(&pt.gmedian(EPS))?.gr()
    );
    println!(
        "Median Correlations of data columns:\n{}",
        pt.transpose()
            .crossfeatures(|v1, v2| v1.medf_correlation(v2).expect("median corr: crossfeatures f64\n"))?
    );
    Ok(())
}

#[test]
fn ustats() -> Result<(), RE> {
    set_seeds(1234567);
    let v1 = ranv_u8(20)?;
    println!("\n{}", (&v1).gr());
    println!("Arithmetic mean: {GR}{:>14.10}{UN}", v1.amean()?);
    println!(
        "Median:          {GR}{:>14.10}{UN}",
        v1.qmedian_by(&mut |a:&u8,b| a.cmp(b), fromop)?
    );
    println!("Geometric mean:  {GR}{:>14.10}{UN}", v1.gmean()?);
    println!("Harmonic mean:   {GR}{:>14.10}{UN}", v1.hmean()?);
    println!("Magnitude:       {GR}{:>14.10}{UN}", v1.vmag());
    println!("Arithmetic {}", v1.ameanstd()?);
    println!("Median     {}", v1.medmad()?);
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
    let v1: Vec<f64> = v.iter().map(|&f| f as f64).collect(); // downcast to f64 here
    println!("Linear transform:\n{}", v1.lintrans()?.gr());
    println!("Arithmetic\t{}", v1.ameanstd()?);
    println!("Median\t\t{}",v1.medmad()?);
    println!("Geometric:\t{}", v1.gmeanstd()?);
    println!("Harmonic:\t{}", v1.hmeanstd()?);
    println!("Autocorrelation:{}", v1.autocorr()?.gr());
    Ok(())
}

#[test]
/// Generic implementation
fn genericstats() -> Result<(), RE> {
    let mut v = vec![1_i32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    println!("\n{}", (&v).gr());
    println!("Arithmetic\t{}", v.ameanstd()?);
    println!("Median\t\t{}",v.medmad()?);
    println!("Geometric\t{}", v.gmeanstd()?);
    println!("Harmonic\t{}", v.hmeanstd()?);
    println!("Weighted Arit.\t{}", v.awmeanstd()?);
    println!("Weighted Geom.\t{}", v.gwmeanstd()?);
    println!("Weighted Harm.\t{}", v.hwmeanstd()?);
    println!("Autocorrelation: {}", v.autocorr()?.gr());
    let median = v.qmedian_by(&mut |a,b| a.cmp(b), fromop)?;
    println!("dfdt:\t\t {}", v.dfdt(median)?.gr());
    v.reverse();
    println!("dfdt(reversed):\t{}", v.dfdt(median)?.gr());
    Ok(())
}

#[test]
fn vecg() -> Result<(), RE> {
    let v1 = vec![
        1_f64, 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.,
    ];
    println!("v1: {}", (&v1).gr());
    let v2 = vec![ 1_f64, -2., 3., - 4., 5., -6., 7., -8., 9., -10., 11., -12., 13., -14., 15.,
    ];
    println!("v2: {}", (&v2).gr());
    println!("Lexical order v1<v2:\t{}", (v1 < v2).gr());
    println!("Median Correlation:\t{}", v1.medf_correlation(&v2)?.gr());
    println!("Pearson's Correlation:\t{}", v1.correlation(&v2).gr());
    println!("Kendall's Correlation:\t{}", v1.kendalcorr(&v2).gr());
    println!("Spearman's Correlation:\t{}", v1.spearmancorr(&v2).gr());
    println!("Euclidian distance:\t{}", v1.vdist(&v2).gr());
    println!("Cityblock distance:\t{}", v1.cityblockd(&v2).gr());
    println!("Vector difference: {}", v1.vsub(&v2).gr());
    println!("Vector sum:        {}", v1.vadd(&v2).gr());
    println!("Scalar product:\t\t{}", v1.dotp(&v2).gr());
    println!("Parallelogram area:\t{}", v1.varea(&v2).gr());
    println!("Arc area:\t\t{}", v1.varc(&v2).gr());
    println!("Entropy v1:\t\t{}", v1.entropy().gr());
    println!("Entropy v2:\t\t{}", v2.entropy().gr());
    println!("Joint Entropy:\t\t{}", v1.jointentropy(&v2)?.gr());
    println!("Dependence:\t\t{}", v1.dependence(&v2)?.gr());
    println!("Independence:\t\t{}", v1.independence(&v2)?.gr());
    println!("\nWedge product:\n{}",v1.wedge(&v2).gr());
    println!("Geometric product:\n{}",v1.geometric(&v2).gr());
    println!("Sine v1v2: {}  v2v1: {} check: {}",v1.sine(&v2).gr(),v2.sine(&v1).gr(),(v1.varea(&v2)/v1.vmag()/v2.vmag()).gr());
    println!("Cosine:\t\t\t{}", v1.cosine(&v2).gr());
    println!("cos^2+sin^2 check:\t{}", (v1.cosine(&v2).powi(2)+v1.sine(&v2).powi(2)).gr());
    println!(
        "Cosine of ranks:\t{}",
        v1.rank(true)
            .indx_to_f64()
            .cosine(&v2.rank(true).indx_to_f64())
            .gr()
    );
    println!("Cos Similarity [0,1]:\t{}", v1.vsim(&v2).gr());
    println!("Cor Similarity [0,1]:\t{}", v1.vcorrsim(&v2)?.gr());
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
    // set_seeds(777);
    let pts1 = ranvv_f64(37,d)?;
    let pts2 = ranvv_f64(50,d)?;
    println!("\nTrend vector (of new random data):\n{}\n", pts1.trend(EPS, pts2)?.gr());
    Ok(())
}

#[test]
fn triangmat() -> Result<(), RE> {
    println!("\n{}", TriangMat::unit(7).gr());
    println!("\n{}", TriangMat::unit(7).to_full().gr());   
    println!("Diagonal: {}",TriangMat::unit(7).diagonal().gr());
    let d = 10_usize;
    let n = 90_usize;
    println!("Testing on a random set of {n} points in {d} dimensional space"); 
    let pts = ranvv_f64_range(n,d, 0.0..=4.0)?;
    // println!("\nTest data:\n{}",pts.gr());
    // let transppt = pts.transpose();
    let cov = pts.covar(&pts.par_gmedian(EPS))?;
    println!("Comediance matrix:\n{cov}");
    println!("Projected to subspace given by [0,2,4,6,9]:\n{}",cov.project(&[0,2,4,6,9]));
    let chol = cov.cholesky()?;
    println!("Cholesky L matrix:\n{chol}");
    println!("Eigenvalues by Cholesky decomposition:\n{}",
        chol.eigenvalues().gr());
    println!("Determinant (their product): {}",chol.determinant().gr());  
    let d = ranv_f64(d)?;
    let dmag = d.vmag();
    let mahamag = chol.mahalanobis(&d)?;
    println!("Random test vector:\n{}", d.gr());
    println!(
        "Its Euclidian magnitude   {GR}{:>8.8}{UN}\
        \nIts Mahalanobis magnitude {GR}{:>8.8}{UN}\
        \nScale factor: {GR}{:>0.8}{UN}",
        dmag,
        mahamag,
        mahamag / dmag
    );
    let (evecs,index) = chol.eigenvectors()?;
    println!("Eigenvectors:\n{}Their sort index by eigenvalues:\n{}",
        evecs.gr(),index.gr());
    println!("Original data PCA reduced to 3 dimensions:\n{}",
        chol.pca_reduction(&pts,3)?.gr());
    Ok(())
}

#[test]
fn mat() -> Result<(), RE> {
    let d = 10_usize;
    let n = 12_usize;
    println!("Testing on a random set of {n} points in {d} dimensional space");
    // set_seeds(1133); 
    let m = ranvv_f64(n,d)?;
    println!("\nTest matrix M:\n{}", m.gr());
    let t = m.transpose();
    println!("\nTransposed matrix M':\n{}", t.gr());
    let v = ranv_f64(d)?;
    println!("\nVector V:\n{}", v.gr());
    println!("\nMV:\n{}", m.leftmultv(&v)?.gr());
    println!("\nVM':\n{}", t.rightmultv(&v)?.gr());
    println!("\nMM':\n{}", t.matmult(&m)?.gr());
    println!("\nM'M:\n{}", m.matmult(&t)?.gr());
    Ok(())
}

#[test]
fn vecvec() -> Result<(), RE> {
    let d = 10_usize;
    let n = 120_usize;
    println!("Testing on a random set of {n} points in {d} dimensional space");
    // set_seeds(113);
    let pts = ranvv_u8(n,d)?;
    println!("First data vector:\n{}",pts[0].gr());
    println!("Joint entropy: {}", pts.jointentropyn()?.gr());
    println!("Dependence:    {}", pts.dependencen()?.gr());
    let (median,recipsum) = pts.gmparts(EPS);
    println!("Approximate dv/dt:\n{}", pts.dvdt(&median)?.gr());
    let outcomes = ranv_u8(n)?;
    println!("\nRandom testing outcomes:\n{}",outcomes.gr());
    println!("wdvdt using outcomes as weigths:\n{}", pts.wdvdt(&outcomes,&median)?.gr());
    println!("wdvdt as wgmedian-gmedian:\n{}", pts.wgmedian(&outcomes,EPS)?.vsub(&median).gr());
    println!("wdvdt as wacentroid-acentroid:\n{}", pts.wacentroid(&outcomes).vsub(&pts.acentroid()).gr());
    let transppt = pts.transpose();
    println!(
        "\nDependencies of columns with test outcomes:\n{}",
        transppt.dependencies(&outcomes)?.gr()
    );
    println!(
        "Correlations with outcomes:\n{}",
        transppt.scalar_fn(|column| Ok(column.correlation(&outcomes)))?.gr());
    let (eccstd, eccmed, eccecc) = pts.eccinfo(&median[..])?;
    let medoid = &pts[eccecc.minindex];
    println!("Medoid: {}", medoid.gr());
    let outlier = &pts[eccecc.maxindex];
    let hcentroid = pts.hcentroid()?;
    let gcentroid = pts.gcentroid()?;
    let acentroid = pts.acentroid();
    let quasimed = pts.quasimedian()?;
    let dists = pts.distsums();
    let md = dists.minmax();
    println!("Medoid and Outlier Total Distances:\n{md}");
    println!("Centroid of total Distances {}", dists.ameanstd()?);
    println!("Median of total distances   {}", dists.medf_unchecked());
    println!(
        "GM's total distances:        {}",
        pts.distsum(&median)?.gr()
    );
    println!(
        "ACentroid's total distances: {}",
        pts.distsum(&acentroid)?.gr()
    );
    println!(
        "GCentroid's total distances: {}",
        pts.distsum(&gcentroid)?.gr()
    );
    println!(
        "HCentroid's total distances: {}",
        pts.distsum(&hcentroid)?.gr()
    );

    println!(
        "\nMedoid, outlier and radii summary:\n{eccecc}\nRadii centroid {eccstd}\nRadii median   {eccmed}");
    let radsindex = pts.radii(&median)?.hashsort_indexed(|&x| x);
    println!(
        "Radii ratio:         {GR}{}{UN}",
        pts.radius(radsindex[0], &median)? / pts.radius(radsindex[radsindex.len() - 1], &median)?
    );
    println!("Madgm:               {}", pts.madgm(&median)?.gr());
    println!("Median's error:      {}", pts.gmerror(&median)?.gr());
    println!("Stdgm:               {}", pts.stdgm(&median)?.gr());
    println!("ACentroid's radius:  {}", acentroid.vdist(&median).gr());
    println!("Quasimed's radius:   {}", quasimed.vdist(&median).gr());
    println!("GCentroid's radius:  {}", gcentroid.vdist(&median).gr());
    println!("HCentroid's radius:  {}", hcentroid.vdist(&median).gr());
    println!("Medoid's radius:     {}", medoid.vdist(&median).gr());
    println!("Outlier's radius:    {}", outlier.vdist(&median).gr());
    println!("Outlier to Medoid:   {}", outlier.vdist(medoid).gr());

    let seccs = pts.radii(&median)?.sorth(|&f| f, true);
    // println!("\nSorted eccs: {}\n", seccs));
    let lqcnt = seccs.binsearch(&(eccmed.centre - eccmed.spread));
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
    let uqcnt = seccs.binsearch(&(eccmed.centre + eccmed.spread));
    println!(
        "Inner three quarters:    {} within radius: {}",
        uqcnt.start.gr(),
        seccs[uqcnt.start - 1].gr()
    );

    let nf = n as f64;

    println!(
        "\nContribution of adding acentroid:    {}",
        acentroid.contrib_newpt(&median, recipsum, nf)?.gr()
    );
    println!(
        "Contribution of adding gcentroid:    {}",
        gcentroid.contrib_newpt(&median, recipsum, nf)?.gr()
    );
    println!(
        "Contribution of removing gcentroid: {}",
        gcentroid
            .contrib_oldpt(&median, recipsum + 1.0 / median.vdist(&gcentroid), nf)? 
            .gr()
    );
    let contribs = pts
        .iter()
        .map(|p|-> Result<f64,RE> { p.contrib_oldpt(&median, recipsum, nf)})
        .collect::<Result<Vec<f64>,RE>>()?;
    println!(
        "\nContributions of removing data points, summary:\n{}\nCentroid: {}\nMedian: {}",
        contribs.minmax(),
        contribs.ameanstd()?,
        contribs.medf_unchecked()
    );
    println!("\nWeighted madgm: {}",pts.wmadgm(&outcomes,&median)?.gr());
    println!("Weighted divs median: {}",pts.wdivsmed(&outcomes,&median)?.gr()); 
    let (divs,wsum) = pts.wdivs(&outcomes,&median)?;
    println!("Weighted divs mean:   {}",(divs.iter().sum::<f64>()/wsum).gr());    
    Ok(())
}

#[test]
fn hulls() -> Result<(), RE> {
    let d = 3_usize;
    let n = 777_usize;
    println!("Testing on a random set of {n} points in {d} dimensional space");
    // set_seeds(77777);
    let pts = ranvv_f64(n,d)?;
    // let wts = rf.ranv_in(n, 0., 100.).getvf64()?;
    let median = pts.gmedian(EPS);
    let zeropts = pts.translate(&median)?;
    let (innerhull, outerhull) = zeropts.hulls();
    if innerhull.is_empty() || outerhull.is_empty() {
        return arith_error("no hull points found"); };
    let mad = zeropts.madgm(&median)?;
    println!("Madgm of zeropts: {}", mad.gr());
    println!(
        "\nInner hull has {}/{} points:\n{}",
        innerhull.len().gr(),
        pts.len().gr(),
        innerhull.yl()
    ); 
    println!(
        "Inner hull min max radii: {} {}\nTheir tm_statistics:\t  {} {}",
        zeropts[*innerhull.first().expect("Empty hullidx")]
            .vmag()
            .gr(),
        zeropts[*innerhull.last().expect("Empty hullidx")]
            .vmag()
            .gr(),
        pts[*innerhull.first().unwrap()]
            .tm_statistic(&median, mad)?
            .gr(),
        pts[*innerhull.last().unwrap()]
            .tm_statistic(&median, mad)?
            .gr()
    );
    let sqradii = zeropts.scalar_fn(|p|Ok(p.vmagsq()))?; 
    let mut radindex = sqradii.mergesort_indexed();
    radindex.reverse();
    println!("Depths of innerhull points:\n{}",
    innerhull
        .iter()
        .map(|&p| zeropts.depth(&radindex,&zeropts[p]))
        .collect::<Result<Vec<f64>,RE>>()?
        .gr()
    );
    println!("Depths ratios of innerhull points:\n{}",
    innerhull
        .iter()
        .map(|&p| zeropts.depth_ratio(&radindex,&zeropts[p]))
        .collect::<Vec<f64>>()
        .gr()
    );


    let sigvec = zeropts.sigvec(&innerhull)?;
    println!(
        "Inner hull sigvec: {}",
        sigvec.gr()
    );

    println!(
        "\nOuter hull has {}/{} points:\n{}",
        outerhull.len().gr(),
        pts.len().gr(),
        outerhull.yl()
    );
    println!(
        "Outer hull min max radii: {} {}\nTheir tm_statistics:\t  {} {}",
        zeropts[*outerhull.last().expect("Empty hullidx")]
            .vmag()
            .gr(),
        zeropts[*outerhull.first().expect("Empty hullidx")]
            .vmag()
            .gr(),
        pts[*outerhull.last().unwrap()]
            .tm_statistic(&median, mad)?
            .gr(),
        pts[*outerhull.first().unwrap()]
            .tm_statistic(&median, mad)?
            .gr()
    );
    println!("Depths of outerhull points: {}",
    outerhull
        .iter()
        .map(|&p| zeropts.depth(&radindex,&zeropts[p]))
        .collect::<Result<Vec<f64>,RE>>()?
        .gr()
    );
    let sigvec = zeropts.sigvec(&outerhull)?;
    println!(
        "Outer hull sigvec: {}",
        sigvec.gr()
    );

    let allptsig = zeropts.sigvec(&Vec::from_iter(0..zeropts.len()))?;
    println!(
        "\nSigvec for all points: {} mod: {}",
        allptsig.gr(), allptsig.vmag().gr()
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
    let (u, r) = a.house_ur()?;
    println!("house_ur u' {u}");
    println!("house_ur r'  {r}");
    let q = u.house_uapply(&unit_matrix(a.len().min(a[0].len())));
    println!(
        "Q matrix\n{}\nOthogonality of Q check (Q'*Q = I):\n{}",
        q.gr(),
        q.transpose().matmult(&q)?.gr()
    );
    println!(
        "Matrix a = QR recreated:\n{}",
        q.matmult(&r.to_full())?.gr()
    );
    Ok(())
}

#[test]
fn geometric_medians() -> Result<(), RE> {
    const NAMES: [&str; 5] = [
        "par_gmedian",
        "gmedian",
        "quasimedian",
        "acentroid",
        "par_acentroid",
    ];
    const CLOSURESU8: [fn(&[Vec<f64>]); 5] = [
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
        |v: &[_]| {
            v.par_acentroid();
        },
    ];
    set_seeds(7777777777_u64); // intialise random numbers generator
                               // Rnum specifies the type of the random numbers required
    println!("\n{YL}Timing Comparisons (in nanoseconds):   {UN}");
    benchvvf64(
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
    println!("\n{RD}Total errors for {ITERATIONS} repeats of {n} points in {d} dimensions:{UN}\n");
    let mut sumg = 0_f64;
    let mut sumr = 0_f64;
    let mut sumq = 0_f64;
    let mut summ = 0_f64;
    let mut sump = 0_f64;
    let mut gm: Vec<f64>;
    for _i in 1..ITERATIONS {
        let pts = ranvv_f64(n,d)?;
        gm = pts.gmedian(EPS);
        sumg += pts.gmerror(&gm)?;
        gm = pts.par_gmedian(EPS);
        sumr += pts.gmerror(&gm)?;
        gm = pts.quasimedian()?;
        sumq += pts.gmerror(&gm)?;
        gm = pts.acentroid();
        summ += pts.gmerror(&gm)?;
        gm = pts.par_acentroid();
        sump += pts.gmerror(&gm)?;
    }
    println!("{MG}par_gmedian   {GR}{sumr:.10}{UN}");
    println!("{MG}gmedian       {GR}{sumg:.10}{UN}");
    println!("{MG}acentroid     {GR}{summ:.10}{UN}");
    println!("{MG}par_acentroid {GR}{sump:.10}{UN}");
    println!("{MG}quasimedian   {GR}{sumq:.10}{UN}\n");
    Ok(())
}
