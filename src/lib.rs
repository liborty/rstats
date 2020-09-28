mod statsf64;
mod statsi64;
mod vecf64impls;
mod vecu8impls;
mod mutvecimpls;
mod vecvecimpls;
pub mod functions;

use crate::functions::GI;
use anyhow::{Result,bail};

/// Median and quartiles
#[derive(Default)]
pub struct Med {
    pub lquartile: f64,
    pub median: f64,
    pub uquartile: f64,
}
impl std::fmt::Display for Med {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "(Lower Q: {}, Median: {}, Upper Q: {})",
            GI(self.lquartile),
            GI(self.median),
            GI(self.uquartile)
        )
    }
}

/// Mean and standard deviation (or std ratio for geometric mean).
#[derive(Default)]
pub struct MStats {
    pub mean: f64,
    pub std: f64,
}
impl std::fmt::Display for MStats {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "mean±std: {}±{}", GI(self.mean), GI(self.std))
    }
}

/// Basic one dimensional (1-d) statistical measures.
/// These methods operate on just one vector (of data) and take no arguments.
pub trait Stats {

    /// Arithmetic mean
    fn amean(self) -> Result<f64> 
        where Self: std::marker::Sized { bail!("amean not implemented for this type")}
    /// Arithmetic mean and standard deviation
    fn ameanstd(self) -> Result<MStats> 
        where Self: std::marker::Sized { bail!("ameanstd not implemented for this type")}
    /// Weighted arithmetic mean
    fn awmean(self) -> Result<f64> 
        where Self: std::marker::Sized { bail!("awmean not implemented for this type")}
    /// Weighted arithmetic men and standard deviation
    fn awmeanstd(self) -> Result<MStats>
        where Self: std::marker::Sized { bail!("awmeanstd not implemented for this type")}
    /// Harmonic mean
    fn hmean(self) -> Result<f64>
        where Self: std::marker::Sized { bail!("hmean not implemented for this type")}
    /// Weighted harmonic mean
    fn hwmean(self) -> Result<f64> 
        where Self: std::marker::Sized { bail!("hwmean not implemented for this type")}
    /// Geometric mean
    fn gmean(self) -> Result<f64>
        where Self: std::marker::Sized { bail!("gmean not implemented for this type")}
    /// Geometric mean and standard deviation ratio
    fn gmeanstd(self) -> Result<MStats>
        where Self: std::marker::Sized { bail!("gmeanstd not implemented for this type")}
    /// Weighed geometric mean
    fn gwmean(self) -> Result<f64> 
        where Self: std::marker::Sized { bail!("gwmean not implemented for this type")}
    /// Weighted geometric mean and standard deviation ratio
    fn gwmeanstd(self) -> Result<MStats>
        where Self: std::marker::Sized { bail!("gwmeanstd not implemented for this type")}
    /// Median and quartiles
    fn median(self) -> Result<Med>
        where Self: std::marker::Sized { bail!("median not implemented for this type")}
    /// Creates vector of ranks for values in self
    fn ranks(self) -> Result<Vec<f64>>
        where Self: std::marker::Sized { bail!("ranks not implemented for this type")}
    /// Creates vector of integer ranks for values in self
    fn iranks(self) -> Result<Vec<i64>>
        where Self: std::marker::Sized { bail!("iranks not implemented for this type")}    
}

/// Vector algebra on one or two vectors.
pub trait Vecf64 {

   /// Scalar multiplication with a vector
    fn smult(self, s: f64) -> Vec<f64>;
    /// Scalar addition to vector
    fn sadd(self, s: f64) -> Vec<f64>; 
    /// Scalar product of two vectors
    fn dotp(self, v: &[f64]) -> f64;
    /// Cosine = a.dotp(b)/(a.vmag*b.vmag)
    fn cosine(self, _v: &[f64]) -> f64; 
    /// Vector subtraction
    fn vsub(self, v: &[f64]) -> Vec<f64>;
    /// Vector addition
    fn vadd(self, v: &[f64]) -> Vec<f64>;
    /// Vector magnitude
    fn vmag(self) -> f64;
    /// Vector magnitude squared
    fn vmagsq(self) -> f64;
    /// Euclidian distance between two points
    fn vdist(self, v: &[f64]) -> f64;
    
    /// Unit vector
    fn vunit(self) -> Vec<f64>;
    /// Area of parallelogram between two vectors (magnitude of cross product)
    fn varea(self, v:&[f64]) -> f64;
    /// Area proportional to the swept arc
    fn varc(self, v:&[f64]) -> f64; 
 
    /// Correlation
    fn correlation(self, _v: &[f64]) -> f64; 
    /// Kendall's tau-b (rank order) correlation
    fn kendalcorr(self, _v: &[f64]) -> f64;
    /// Spearman's rho (rank differences) correlation
    fn spearmancorr(self, _v: &[f64]) -> f64;
    /// Kazutsugi Spearman's corelation against just five distances (to outcomes classes)
    fn kazutsugi(self) -> f64;
    /// Autocorrelation
    fn autocorr(self) -> f64;
 
    /// Minimum, minimum's index, maximum, maximum's index.
    fn minmax(self) -> (f64, usize, f64, usize); 
    /// Linear transformation to [0,1]
    fn lintrans(self) -> Vec<f64>;
    /// Sorted vector
    fn sortf(self) -> Vec<f64>;
    /// Reursive merge sort building ranks in n*log(n)
    fn mergerank(self) -> Vec<usize>;
    /// Immutable merge sort, makes an index
    fn mergesort(self, i:usize, n:usize) -> Vec<usize>;
}

/// Minimal support also for Vec[u8]
pub trait Vecu8 {

    /// Scalar multiplication with a vector
    fn smult(self, s: f64) -> Vec<f64>;
    /// Scalar addition to vector
    fn sadd(self, s: f64) -> Vec<f64>;
    /// Scalar product of two vectors
    fn dotp(self, v: &[f64]) -> f64;
    /// Vector magnitude squared (sum of squares)
    fn vmagsq(self) -> f64;
    /// Area proportional to the swept arc
    fn varc(self, v:&[f64]) -> f64;
    /// Euclidian distance 
    fn vdist(self, v: &[f64]) -> f64;  
    
}

/// Mutable primitive vector operations.
/// Some of the Vectors trait methods reimplemented for efficiency to mutate in-place.
pub trait MutVectors {
    /// mutable multiplication by a scalar
    fn mutsmult(self, _s: f64) where Self: std::marker::Sized {}  
    /// mutable vector subtraction
    fn mutvsub(self, _v: &[f64]) where Self: std::marker::Sized {}
    fn mutvsubu8(self, _v: &[u8]) where Self: std::marker::Sized {} 
    /// mutable vector addition
    fn mutvadd(self, _v: &[f64]) where Self: std::marker::Sized {}
    fn mutvaddu8(self, _v: &[u8]) where Self: std::marker::Sized {}
     /// mutably makes into a unit vector
    fn mutvunit(self) where Self: std::marker::Sized {}
    /// sort in place
    fn mutsortf(self) where Self: std::marker::Sized {} 
}

/// Minimal support also for sets of bytes of type Vec<Vec<u8>>
pub trait VecVecu8 {
    /// Centroid = euclidian mean of a set of points  
    fn acentroid(self) -> Vec<f64>; 
    fn nmedian(self, eps:f64) -> Vec<f64>;
    fn betterpoint(self, v: &[f64]) -> (f64, Vec<f64>);
}

/// Methods applicable to sets of vectors.
pub trait VecVec {
    /// Centroid = euclidian mean of a set of points
    fn acentroid(self) -> Vec<f64>;
    /// Sums of distances from each point to all other points
    fn distsums(self) -> Vec<f64>;
    /// Sum of distances from one point given by indx
    fn distsuminset(self, indx: usize) -> f64;
    /// Sum of distances from arbitrary point (v) to all the points in self   
    fn distsum(self, v: &[f64]) -> f64;
    /// Individual distances from any point v (typically not in self) to all the points in self.    
    fn dists(self, v: &[f64]) -> Vec<f64>;
    /// Medoid and Outlier (by distance) of a set of points
    fn medoid(self) -> (f64, usize, f64, usize);

    /// Eccentricity vectors from each point
    fn eccentricities(self) -> Vec<Vec<f64>>;
    /// Ecentricity scalar measure of an internal point given by indx
    fn eccentrinset(self, indx: usize) -> f64;
    /// Eccentricity scalar measure and vector of any point     
    fn veccentr(self, thisp: &[f64]) -> (f64, Vec<f64>);
    /// Eccentricity scalar measure only, of any point
    fn ecc(self, v: &[f64]) -> f64;
    /// magnitudes of a set of vectors
    fn mags(self) -> Vec<f64>;
    /// scaled magnitudes (typically of eccentricities measures)
    fn scalarecc(self) -> Vec<f64>;
    /// Median and quartiles of eccentricities (new robust measure of spread of a multivariate sample)
    fn moe(self) -> (MStats,Med);
    /// Medoid and Outlier as defined by eccentricities.
    fn emedoid(self) -> (f64, usize, f64, usize);

    /// Geometric median of a set
    fn nmedian(self, eps: f64) -> Vec<f64>;
    /// Betterpoint gives new approximation to nmedian
    fn betterpoint(self, v: &[f64]) -> (f64, Vec<f64>);
    /// Trend between two sets
    fn trend(self, eps: f64, v: Vec<Vec<f64>>) -> Vec<f64>;
    /// Subtract m from all points - e.g. transform to zero median form
    fn translate(self, m: &[f64]) -> Vec<Vec<f64>>;
}

/// Methods to manipulate indices
pub trait Indices {

    /// Reverse index
    fn revindex(self) -> Vec<usize>;
    /// Collects values from `v` as per indices in self.
    fn unindex(self, v:&[f64]) -> Vec<f64>;
    /// Pearson's correlation coefficient of a two $[usize] slices,
    /// typically the ranked data.  
    fn ucorrelation(self, v: &[usize]) -> f64;
  
}

