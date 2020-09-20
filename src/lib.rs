mod f64impls;
pub mod functions;
mod i64impls;
mod vecimpls;
mod mutvecimpls;
mod vecvecimpls;

use crate::functions::{GreenIt};
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
            GreenIt(self.lquartile),
            GreenIt(self.median),
            GreenIt(self.uquartile)
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
        write!(f, "mean±std: {}±{}", GreenIt(self.mean), GreenIt(self.std))
    }
}

/// Basic one dimensional (1-d) statistical measures.
/// All these methods operate on only one vector (of data). They take no arguments.
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
    fn iranks(self) -> Result<Vec<i64>>
        where Self: std::marker::Sized { bail!("iranks not implemented for this type")}    
}

/// Vector algebra on one or two vectors.
pub trait Vectors {

    /// Scalar product of two vectors
    fn dotp(self, v: &[f64]) -> f64;
    /// Cosine = a.dotp(b)/(a.vmag*b.vmag)
    fn  cosine(self, _v: &[f64]) -> f64; 
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
    /// Scalar multiplication
    fn smult(self, s: f64) -> Vec<f64>;
    /// Scalar addition to vector
    fn sadd(self, s: f64) -> Vec<f64>;
    /// Unit vector
    fn vunit(self) -> Vec<f64>;
    /// Area of parallelogram between two vectors (magnitude of cross product)
    fn varea(self, v:&[f64]) -> f64;
    /// Area proportional to the swept arc
    fn dv(self, v:&[f64]) -> f64; 

    /// Correlation
    fn correlation(self, _v: &[f64]) -> Result<f64>
    where Self: std::marker::Sized { bail!("correlation not implemented for this type") }
 
    /// Kendall's tau-b (rank order) correlation
    fn kendalcorr(self, _v: &[f64]) -> Result<f64>
    where Self: std::marker::Sized { bail!("kendalcorr not implemented for this type") }
    /// Spearman's rho (rank differences) correlation
    fn spearmancorr(self, _v: &[f64]) -> Result<f64>
    where Self: std::marker::Sized { bail!("spearmancorr not implemented for this type") }
    /// Autocorrelation
    fn autocorr(self) -> Result<f64>
    where Self: std::marker::Sized { bail!("autocorr not implemented for this type") }

    /// Minimum, minimum's index, maximum, maximum's index.
    fn minmax(self) -> (f64, usize, f64, usize); 
    /// Linear transformation to [0,1]
    fn lintrans(self) -> Vec<f64>;
    /// Sorted vector
    fn sortf(self) -> Vec<f64>;
}

/// Mutable primitive vector operations.  
/// Some of the Vec trait methods reimplemented to mutate in-place (for efficiency).
pub trait MutVectors {
    /// mutable multiplication by a scalar
    fn mutsmult(self, _s: f64) where Self: std::marker::Sized {}  
    /// mutable vector subtraction
    fn mutvsub(self, _v: &[f64]) where Self: std::marker::Sized {}
    /// mutable vector addition
    fn mutvadd(self, _v: &[f64]) where Self: std::marker::Sized {}
    /// mutably makes into a unit vector
    fn mutvunit(self) where Self: std::marker::Sized {}
    /// sort in place
    fn mutsortf(self) where Self: std::marker::Sized {}
}

/// Methods applicable to sets of vectors
pub trait VecVec {
    /// Centroid = euclidian mean of a set of points
    fn acentroid(self) -> Vec<f64>;
    /// Sums of distances from each point to all other points
    fn distances(self) -> Vec<f64>;
    /// Sum of distances from one point given by indx
    fn distsuminset(self, indx: usize) -> f64;
    /// Sum of distances from arbitrary point (v) to all the points in self   
    fn distsum(self, v: &[f64]) -> f64;
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
