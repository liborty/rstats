#![warn(missing_docs)]
//! Statistics, Vector Algebra, 
//! Characterising Multidimensional Data, Machine Learning,
//! Data Analysis

/// Custom Error RError
pub mod error;
/// Basic statistics on a single generic vector
pub mod statsg;
/// Vector Algebra on two generic vectors
pub mod vecg;
/// Stats and Vector Algebra on one or two u8 vectors 
pub mod vecu8;
/// Vector Algerba mutating an f64 vector
pub mod mutvec;
/// Multidimensional operations on sets of vectors
pub mod vecvec;
/// Multidimensional operations on sets of vectors, with additional inputs
pub mod vecvecg;

// reexporting useful related methods
pub use indxvec::{MinMax,F64,Printing,here};
pub use medians::{Med,Median};
/// simple error handling
// use anyhow::{Result,bail}; 
use crate::error::RError;
use core::iter::FromIterator;

/// Mean and standard deviation (or std ratio for geometric mean)
#[derive(Default)]
pub struct MStats {
    /// the mean of some kind (geometric, arithmetic, harmonic...)
    pub mean: f64,
    /// standard deviation (measure of data spread)
    pub std: f64
}
impl std::fmt::Display for MStats {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "mean±std: {}±{}", self.mean.gr(), self.std.gr())
    }
}

// Auxiliary Functions //

/// Downcast of a whole i64 slice to f64  
pub fn i64tof64(s: &[i64]) -> Vec<f64> {
    s.iter().map(| &x | x as f64).collect()
}

/// Sum of natural numbers 1..n.
/// Also the size of an upper or lower triangle 
/// of a square array (including the diagonal)
/// to exclude the diagonal, use `sumn(n-1)`
pub fn sumn(n: usize) -> f64 {
    (n * (n + 1)) as f64 / 2.
}


// Traits

/// Statistical measures of a single variable (one generic vector of data) and 
/// vector algebra applicable to a single (generic) vector. 
/// Thus these methods take no arguments.
/// There is just one limitation: data of end type `i64` has to be explicitly converted to `f64`.
/// That is to raise awareness that, in this particular case, some precision may be lost. 
/// Function `statsg::i64tof64(&s)` will convert the whole slice.
pub trait Stats { 

    /// Vector magnitude
    fn vmag(self) -> f64;
    /// Vector magnitude squared (sum of squares)
    fn vmagsq(self) -> f64;
    /// vector with reciprocal components
    fn vreciprocal(self) -> Result<Vec<f64>,RError>;
    /// vector with inverse magnitude 
    fn vinverse(self) -> Result<Vec<f64>,RError>;
    /// negated vector (all components swap sign)
    fn negv(self) -> Vec<f64>;
    /// Unit vector
    fn vunit(self) -> Vec<f64>;    
    /// Arithmetic mean
    fn amean(self) -> Result<f64,RError>; 
       // where Self: std::marker::Sized { bail!("amean not implemented for this type")}
    /// Arithmetic mean and standard deviation
    fn ameanstd(self) -> Result<MStats,RError>; 
       // where Self: std::marker::Sized { bail!("ameanstd not implemented for this type")}
    /// Weighted arithmetic mean
    fn awmean(self) -> Result<f64,RError>;  
       // where Self: std::marker::Sized { bail!("awmean not implemented for this type")}
    /// Weighted arithmetic men and standard deviation
    fn awmeanstd(self) -> Result<MStats,RError>;
        // where Self: std::marker::Sized { bail!("awmeanstd not implemented for this type")}
    /// Harmonic mean
    fn hmean(self) -> Result<f64,RError>;
        // where Self: std::marker::Sized { bail!("hmean not implemented for this type")}
    /// Harmonic mean and experimental standard deviation 
    fn hmeanstd(self) -> Result<MStats,RError>;
        // where Self: std::marker::Sized { bail!("hmeanstd not implemented for this type")}
    /// Weighted harmonic mean
    fn hwmean(self) -> Result<f64,RError>; 
        // where Self: std::marker::Sized { bail!("hwmean not implemented for this type")}
    /// Weighted harmonic mean and standard deviation 
    fn hwmeanstd(self) -> Result<MStats,RError>;
       // where Self: std::marker::Sized { bail!("hwgmeanstd not implemented for this type")}
    /// Geometric mean
    fn gmean(self) -> Result<f64,RError>;
        // where Self: std::marker::Sized { bail!("gmean not implemented for this type")}
    /// Geometric mean and standard deviation ratio
    fn gmeanstd(self) -> Result<MStats,RError>;
        // where Self: std::marker::Sized { bail!("gmeanstd not implemented for this type")}
    /// Weighed geometric mean
    fn gwmean(self) -> Result<f64,RError>; 
        // where Self: std::marker::Sized { bail!("gwmean not implemented for this type")}
    /// Weighted geometric mean and standard deviation ratio
    fn gwmeanstd(self) -> Result<MStats,RError>;
        // where Self: std::marker::Sized { bail!("gwmeanstd not implemented for this type")}
    /// Zero median data, obtained by subtracting the median
    fn zeromedian(self) -> Result<Vec<f64>,RError>;
       // where Self: std::marker::Sized { bail!("gwmeanstd not implemented for this type")}
    /// Probability density function of a sorted slice
    fn pdf(self) -> Vec<f64>;
    /// Information (entropy) in nats
    fn entropy(self) -> f64;
    /// (Auto)correlation coefficient of pairs of successive values of (time series) variable.
    fn autocorr(self) -> f64; 
    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Vec<f64>;
    /// Reconstructs the full symmetric matrix from its lower diagonal compact form
    fn symmatrix(self) -> Vec<Vec<f64>>;
   }

/// Vector Algebra on two vectors (represented here as generic slices).
/// Also included are scalar operations on the `self` vector.
pub trait Vecg { 
    /// Scalar addition to vector
    fn sadd<U>(self, s:U) -> Vec<f64> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Scalar multiplication of vector, creates new vec
     fn smult<U>(self, s:U) -> Vec<f64> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Scalar product 
    fn dotp<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Cosine of angle between two slices
    fn cosine<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Vectors' subtraction (difference)
    fn vsub<U>(self, v:&[U]) -> Vec<f64> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Vectors difference as unit vector (done together for efficiency)
    fn vsubunit<U>(self, v: &[U]) -> Vec<f64> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>; 
    /// Vector addition
    fn vadd<U>(self, v:&[U]) -> Vec<f64> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>; 
    /// Euclidian distance 
    fn vdist<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Weighted distance of self:&[T],weighted by ws:&[U],to v:&[V] 
    fn wvdist<U,V>(self, ws:&[U],v:&[V]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>, V: Copy, f64:From<V>;
    /// Euclidian distance squared
    fn vdistsq<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>; 
     /// cityblock distance
    fn cityblockd<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Magnitude of the cross product |a x b| = |a||b|sin(theta).
    /// Attains maximum `|a|.|b|` when the vectors are othogonal.
    fn varea<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Area proportional to the swept arc
     fn varc<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>; 
    /// Vector similarity S in the interval [0,1]: S = (1+cos(theta))/2
    fn vsim<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Positive dotp [0,2|a||b|] 
    fn pdotp<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// We define vector dissimilarity D in the interval [0,1]: D = 1-S = (1-cos(theta))/2
    fn vdisim<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Lower triangular part of a covariance matrix for a single f64 vector.
    fn covone<U>(self, m:&[U]) -> Vec<f64> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Kronecker product of two vectors 
    fn kron<U>(self, m:&[U]) -> Vec<f64> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>; 
    /// Outer product of two vectors 
    fn outer<U>(self, m:&[U]) -> Vec<Vec<f64>> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>; 
    /// Joint probability density function 
    fn jointpdf<U>(self,v:&[U]) -> Vec<f64> where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;     
    /// Joint entropy of &[T],&[U] in nats (units of e)
    fn jointentropy<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Statistical pairwise dependence of &[T] &[U] variables in the range [0,1] 
    fn dependence<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;   
    /// Statistical pairwise independence in the range [0,1] based on joint entropy
    fn independence<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>; 
    /// Cholesky decomposition of positive definite matrix into LLt
    // fn cholesky(self) -> Vec<f64>;
    /// Pearson's correlation.  
    fn correlation<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Median based correlation
    fn mediancorr<U>(self, v: &[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Kendall Tau-B correlation. 
    fn kendalcorr<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>;
    /// Spearman rho correlation.
    fn spearmancorr<U>(self, v:&[U]) -> f64 where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>; 
    /// Change to gm that adding point p will cause
    fn contribvec_newpt(self,gm:&[f64],recips:f64) -> Vec<f64>;
    /// Magnitude of change to gm that adding point p will cause
    fn contrib_newpt(self,gm:&[f64],recips:f64) -> f64;
    /// Contribution an existing set point p has made to the gm
    fn contribvec_oldpt(self,gm:&[f64],recips:f64) -> Vec<f64>;
    /// Contribution removing an existing p will make (as a negative number)
    fn contrib_oldpt(self,gm:&[f64],recips:f64) -> f64;  
}

/// Mutable vector operations that take one generic argument. 
/// A few of the essential `Vecg` methods are reimplemented here 
/// to mutate `self` in-place (only for f64). 
/// This is for efficiency and convenience, for example, in
/// vector iterative methods.
pub trait MutVecg {
    /// mutable multiplication by a scalar
    fn mutsmult<U>(self, _s:U) where U: Copy+PartialOrd, f64: From<U>;  
    /// mutable vector subtraction
    fn mutvsub<U>(self, _v: &[U]) where U: Copy+PartialOrd, f64: From<U>; 
    /// mutable vector addition
    fn mutvadd<U>(self, _v: &[U]) where U: Copy+PartialOrd, f64: From<U>;  

    /// Invert the magnitude
    fn minvert(self); 
    /// Negate the vector (all components swap sign)
    fn mneg(self);
    /// Make into a unit vector
    fn munit(self);
    /// Linearly transform to interval [0,1]
    fn mlintrans(self);
}

/// Methods specialised to, or more efficient for `&[u8]`
pub trait Vecu8 {    
    /// Probability density function of bytes data
    fn pdfu8(self) -> Vec<f64>;
    /// Information (entropy) of &[u8] (in nats)
    fn entropyu8(self) -> f64;
    /// Joint probability density function (here just co-occurence counts) 
    /// of paired values in two vectors of bytes of the same length.
    /// Needs n^2 x 32bits of memory. 
    fn jointpdfu8(self, v:&[u8]) -> Vec<Vec<u32>>;
    /// Joint entropy of &[u8],&[u8] (in nats)
    fn jointentropyu8(self, v:&[u8]) -> f64;
    /// Statistical pairwise dependence of two &[u8] variables in the range [0,1] 
    fn dependenceu8(self, v:&[u8]) -> f64;   
    /// Independence in the range [1,2] of two &[u8] variables 
    fn independenceu8(self, v:&[u8]) -> f64;
}

/// Methods applicable to a single argument: a vector of vectors of generic end type.
/// Operations on a set of multidimensional vectors.
pub trait VecVec<T> {
    /// Transpose vec of vecs like a classical array
    fn transpose(self) -> Vec<Vec<T>>;
    /// Joint probability density function of n matched slices of the same length
    fn jointpdfn(self) -> Vec<f64>;
    /// Joint entropy between a set of vectors of the same length
    fn jointentropyn(self) -> f64;
    /// Independence (component wise) of a set of vectors.
    fn dependencen(self) -> f64; 
    /// Flattened lower triangular relations matrix between columns of self 
    fn crossfeatures(self,f:fn(&[T],&[T])->f64) -> Vec<f64>;
    /// Sum of nd points (or vectors)
    fn sumv(self) -> Vec<f64>; 
    /// Arithmetic Centre = euclidian mean of a set of points
    fn acentroid(self) -> Vec<f64>;
    /// Geometric Centroid
    fn gcentroid(self) -> Vec<f64>;
    /// Harmonic Centroid = harmonic mean of a set of points
    fn hcentroid(self) -> Vec<f64>; 
    /// Possible first iteration point for geometric medians
    fn firstpoint(self) -> Vec<f64>;
    /// Sums of distances from each point to all other points
    fn distsums(self) -> Vec<f64>;    
    /// Medoid distance, its index, outlier distance, its index
    fn medout(self,gm:&[f64]) -> MinMax<f64>; 
    /// Fast sums of distances from each point to all other points 
    fn distsuminset(self, indx: usize) -> f64;
    /// Next approx gm computed at a member point given by its indx
    fn nxmember(self, indx: usize) -> Vec<f64>;
    /// Like gmparts, except only does one iteration from any non-member point g
    fn nxnonmember(self, g:&[f64]) -> (Vec<f64>,Vec<f64>,f64);
    /// Estimated eccentricity vectors from each member point
    fn eccentricities(self) -> Vec<Vec<f64>>;
    /// Exact eccentricity vectors to each member point by using the gm 
    fn exacteccs(self, gm:&[f64]) -> Vec<Vec<f64>>; 
    /// Median and quartiles of radii (new robust measure of spread of a multivariate sample)
    fn eccinfo(self, gm:&[f64]) -> (MStats,Med,MinMax<f64>) where Vec<f64>:FromIterator<f64>;
    /// Quasi median, recommended only for comparison purposes
    fn quasimedian(self) -> Vec<f64>;
    /// Geometric median estimate's error
    fn gmerror(self,gm:&[f64]) -> f64;
    /// MADGM, absolute deviations from geometric median: stable nd data spread estimator
    fn madgm(self, gm: &[f64]) -> f64;
    /// Proportions of points found along each axis
    fn tukeyvec(self, gm: &[f64]) -> Vec<f64>;
    /// Returns sorted eccentricities magnitudes
    fn sortedeccs(self, ascending:bool, gm:&[f64]) -> Vec<f64> where F64:From<T>;
    /// New algorithm for geometric median, to accuracy eps    
    fn gmedian(self, eps: f64) -> Vec<f64>;
    /// Point-by-point geometric median
    fn pmedian(self, eps: f64) -> Vec<f64>;
    /// Like `gmedian` but returns the sum of unit vecs and the sum of reciprocals of distances.
    fn gmparts(self, eps: f64) -> (Vec<f64>,Vec<f64>,f64);
}

/// Methods applicable to vector of vectors of generic end type and one argument
/// of a similar kind.
pub trait VecVecg<T,U> { 

    /// Weighted sum of nd points (or vectors)
    fn wsumv(self,ws: &[U]) -> Vec<f64>;
    /// Weighted Arithmetic Centre = weighted euclidian mean of a set of points
    fn wacentroid(self,ws: &[U]) -> Vec<f64>;
    /// Trend between two sets
    fn trend(self, eps: f64, v: Vec<Vec<U>>) -> Vec<f64>;
    /// Subtract m from all points - e.g. transform to zero median form
    fn translate(self, m: &[U]) -> Vec<Vec<f64>>;
    /// Transform nd data to zeromedian for 
    fn zerogm(self, gm: &[f64]) -> Vec<Vec<f64>>; 
    /// Dependencies of vector m on each vector in self
    fn dependencies(self, m: &[U]) -> Vec<f64>;
    /// (Median) correlations of m with each vector in self
    fn correlations(self, m: &[U]) -> Vec<f64>;
    /// Sum of distances from arbitrary point (v) to all the points in self      
    fn distsum(self, v: &[U]) -> f64;
    /// Individual distances from any point v (typically not in self) to all the points in self.    
    fn dists(self, v: &[U]) -> Vec<f64>;
    /// ( wgm, sorted eccentricities magnitudes, associated cpdf )
    fn wsortedeccs(self,ws:&[U],gm:&[f64]) -> (Vec<f64>,Vec<f64>) where F64:From<T>; 
    /// Like wgmparts, except only does one iteration from any non-member point g
    fn wnxnonmember(self, ws:&[U], g:&[f64]) -> (Vec<f64>,Vec<f64>,f64); 
    /// The weighted geometric median to accuracy eps 
    fn wgmedian(self, ws: &[U], eps: f64) -> Vec<f64>;
    /// Like `wgmedian` but returns also the sum of unit vecs and the sum of reciprocals. 
    fn wgmparts(self, ws:&[U],eps: f64) -> (Vec<f64>,Vec<f64>,f64);
    /// wmadgm median of weighted absolute deviations from weighted gm: stable nd data spread estimator
    fn wmadgm(self, ws: &[U], wgm: &[f64]) -> f64;     
    /// Flattened lower triangular part of a covariance matrix of a Vec of f64 vectors.
    fn covar(self, med:&[U]) -> Vec<f64>;  
    /// Flattened lower triangular part of a covariance matrix for weighted f64 vectors.
    fn wcovar(self, ws:&[U], m:&[f64]) -> Vec<f64>;
    /// Flattened comediance matrix for f64 vectors in self.
    /// Similar to `covar` above but medians instead of means are returned.
    fn comed(self, m:&[U]) -> Vec<f64>;
    /// Flatteened comediance matrix for weighted f64 vectors.
    /// Similar to `wcovar` above but medians instead of means are returned.
    fn wcomed(self, ws:&[U], m:&[f64]) -> Vec<f64>;
}
