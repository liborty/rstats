#![warn(missing_docs)]
//! Statistics, Linear Algebra, Information Measures, Cholesky Matrix Decomposition,
//! Mahalanobis Distance, Multidimensional Data Analysis, Machine Learning and more ...

/// Custom error RError
pub mod error;
/// Vector algebra mutating an f64 vector
pub mod mutvec;
/// Basic statistics on a single generic vector
pub mod stats;
/// Associated functions implemented for struct TriangMat
pub mod triangmat;
/// Vector algebra on two generic vectors
pub mod vecg;
/// Stats and vector algebra on one or two u8 vectors
pub mod vecu8;
/// Multidimensional operations on sets of vectors
pub mod vecvec;
/// Multidimensional operations on sets of vectors, with additional inputs
pub mod vecvecg;

use crate::error::RError;
// reexporting useful related methods
pub use indxvec::{printing::*, MinMax, Printing};
pub use medians::{MStats, MedError, Median, Medianf64};

/// Shorthand type for returned errors with message payload
pub type RE = RError<String>;

// Auxiliary Functions

/// Convenience dummy function for quantify closure
pub fn noop(f: &f64) -> f64 {
    *f
}

/// Convenience From quantification invocation
pub fn fromop<T>(f: &T) -> f64
where
    T: Clone,
    f64: From<T>,
{
    f64::from(f.clone())
}

/// t_statistic in 1d: (value-centre)/dispersion
/// generalized to any measure of central tendency and dispersion
pub fn t_stat(val: f64, mstats: MStats) -> f64 {
    (val - mstats.centre) / mstats.dispersion
}

/// Sum of natural numbers 1..n.
/// Also the size of an upper or lower triangle
/// of a square array (including the diagonal)
/// to exclude the diagonal, use `sumn(n-1)`
pub fn sumn(n: usize) -> usize {
    n * (n + 1) / 2
}

/// Generates full nxn unit (identity) matrix
pub fn unit_matrix(n: usize) -> Vec<Vec<f64>> {
    let mut res: Vec<Vec<f64>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = vec![0.; n];
        row[i] = 1.0;
        res.push(row);
    }
    res
}

/// Compact Triangular Matrix.
/// TriangMat is typically result of some matrix calculations,
/// so concrete end-type f64 is used for simplicity and accuracy.
/// Data is of length `n*(n+1)/2` instead of `n*n`, saving memory.  
/// `.kind == 0` is plain lower triangular matrix.  
/// `.kind == 1` is antisymmetric square matrix.  
/// `.kind == 2` is symmetric square matrix.  
/// `.kind == 3` is upper triangular matrix (transposed lower).  
/// `.kind == 4` is upper (transposed lower), antisymmetric.  
/// `.kind == 5` is unnecessary, as transposed symmetric matrix is unchanged.  
/// Simply adding (or subtracting) 3 to .kind implicitly transposes the matrix.
/// `.kind > 2 are all transposed, individual variants are determined by kind % 3.
/// The size of the implied full square matrix, nxn, is not explicitly stored.
/// It is obtained by solving the quadratic equation:
/// `((((8 * s + 1) as f64).sqrt() - 1.) / 2.) as usize;`
/// where `s = triangmat.len()` = `n*(n+1)/2`
#[derive(Default, Clone)]
pub struct TriangMat {
    /// Matrix kind encoding: 0=plain, 1=antisymmetric, 2=symmetric, +3=transposed
    pub kind: usize,
    /// Packed 1d vector of triangular matrix data of size `(n+1)*n/2`
    pub data: Vec<f64>,
}

// Traits

/// Statistical measures of a single variable (one generic vector of data) and
/// vector algebra applicable to a single (generic) vector.
/// Thus these methods take no arguments.
pub trait Stats {
    /// Vector magnitude
    fn vmag(self) -> f64;
    /// Vector magnitude squared (sum of squares)
    fn vmagsq(self) -> f64;
    /// vector with reciprocal components
    fn vreciprocal(self) -> Result<Vec<f64>, RE>;
    /// vector with inverse magnitude
    fn vinverse(self) -> Result<Vec<f64>, RE>;
    /// negated vector (all components swap sign)
    fn negv(self) -> Result<Vec<f64>, RE>;
    /// Unit vector
    fn vunit(self) -> Result<Vec<f64>, RE>;
    /// Harmonic spread from median
    fn hmad(self) -> Result<f64, RE>; 
    /// Arithmetic mean
    fn amean(self) -> Result<f64, RE>;
    /// Arithmetic mean and standard deviation
    fn ameanstd(self) -> Result<MStats, RE>;
    /// Weighted arithmetic mean
    fn awmean(self) -> Result<f64, RE>;
    /// Weighted arithmetic men and standard deviation
    fn awmeanstd(self) -> Result<MStats, RE>;
    /// Harmonic mean
    fn hmean(self) -> Result<f64, RE>; 
    /// Harmonic mean and experimental standard deviation
    fn hmeanstd(self) -> Result<MStats, RE>;
    /// Weighted harmonic mean
    fn hwmean(self) -> Result<f64, RE>;
    /// Weighted harmonic mean and standard deviation
    fn hwmeanstd(self) -> Result<MStats, RE>;
    /// Geometric mean
    fn gmean(self) -> Result<f64, RE>;
    /// Geometric mean and standard deviation ratio
    fn gmeanstd(self) -> Result<MStats, RE>;
    /// Weighed geometric mean
    fn gwmean(self) -> Result<f64, RE>;
    /// Weighted geometric mean and standard deviation ratio
    fn gwmeanstd(self) -> Result<MStats, RE>;
    /// Probability density function of a sorted slice
    fn pdf(self) -> Vec<f64>;
    /// Information (entropy) in nats
    fn entropy(self) -> f64;
    /// (Auto)correlation coefficient of pairs of successive values of (time series) variable.
    fn autocorr(self) -> Result<f64, RE>;
    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Result<Vec<f64>, RE>;
    /// Linearly weighted approximate time series derivative at the last (present time) point.
    fn dfdt(self) -> Result<f64, RE>;
    /// Householder reflection
    fn house_reflector(self) -> Vec<f64>;
}

/// Vector Algebra on two vectors (represented here as generic slices).
/// Also included are scalar operations on the `self` vector.
pub trait Vecg {
    /// nd t_statistic of self against geometric median and madgm spread
    fn t_statistic(self, gm: &[f64], madgm: f64) -> Result<f64, RE>;
    /// Dot product of vector self with column c of matrix v
    fn columnp<U>(self, c: usize, v: &[Vec<U>]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Scalar addition to vector
    fn sadd<U>(self, s: U) -> Vec<f64>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Scalar multiplication of vector, creates new vec
    fn smult<U>(self, s: U) -> Vec<f64>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Scalar product
    fn dotp<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Product with Tukeyvec of hemispheric counts.
    fn dottukey(self, tukey: &[f64]) -> Result<f64, RE>;
    /// Cosine of angle between two slices
    fn cosine<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Vectors' subtraction (difference)
    fn vsub<U>(self, v: &[U]) -> Vec<f64>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Vectors difference as unit vector (done together for efficiency)
    fn vsubunit<U>(self, v: &[U]) -> Vec<f64>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Vector addition
    fn vadd<U>(self, v: &[U]) -> Vec<f64>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Euclidian distance
    fn vdist<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Weighted arithmetic mean of `self:&[T]`, scaled by `ws:&[U]`
    fn wvmean<U>(self, ws: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Weighted distance of self:&[T],weighted by ws:&[U],to v:&[V]
    fn wvdist<U, V>(self, ws: &[U], v: &[V]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>,
        V: Copy,
        f64: From<V>;
    /// Euclidian distance squared
    fn vdistsq<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// cityblock distance
    fn cityblockd<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Magnitude of the cross product |a x b| = |a||b|sin(theta).
    /// Attains maximum `|a|.|b|` when the vectors are othogonal.
    fn varea<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Area proportional to the swept arc
    fn varc<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Vector similarity S in the interval [0,1]: S = (1+cos(theta))/2
    fn vsim<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Positive dotp [0,2|a||b|]
    fn pdotp<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// We define vector dissimilarity D in the interval [0,1]: D = 1-S = (1-cos(theta))/2
    fn vdisim<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Lower triangular part of a covariance matrix for a single f64 vector.
    fn covone<U>(self, m: &[U]) -> TriangMat
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Kronecker product of two vectors
    fn kron<U>(self, m: &[U]) -> Vec<f64>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Outer product of two vectors
    fn outer<U>(self, m: &[U]) -> Vec<Vec<f64>>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Joint probability density function
    fn jointpdf<U>(self, v: &[U]) -> Result<Vec<f64>, RE>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Joint entropy of &[T],&[U] in nats (units of e)
    fn jointentropy<U>(self, v: &[U]) -> Result<f64, RE>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Statistical pairwise dependence of &[T] &[U] variables in the range [0,1]
    fn dependence<U>(self, v: &[U]) -> Result<f64, RE>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Statistical pairwise independence in the range [0,1] based on joint entropy
    fn independence<U>(self, v: &[U]) -> Result<f64, RE>
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Pearson's correlation.  
    fn correlation<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Kendall Tau-B correlation.
    fn kendalcorr<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Spearman rho correlation.
    fn spearmancorr<U>(self, v: &[U]) -> f64
    where
        U: Copy + PartialOrd + Into<U>,
        f64: From<U>;
    /// Change to gm that adding point self will cause
    fn contribvec_newpt(self, gm: &[f64], recips: f64) -> Vec<f64>;
    /// Normalized magnitude of change to gm that adding point self will cause
    fn contrib_newpt(self, gm: &[f64], recips: f64, nf:f64) -> f64;
    /// Contribution of removing point self
    fn contribvec_oldpt(self, gm: &[f64], recips: f64) -> Vec<f64>;
    /// Normalized contribution of removing point self (as negative scalar)
    fn contrib_oldpt(self, gm: &[f64], recips: f64, nf: f64) -> f64;
    /// Householder reflect
    fn house_reflect<U>(self, x: &[U]) -> Vec<f64>
    where
        U: Copy + PartialOrd + std::fmt::Display,
        f64: From<U>;
}

/// Mutable operations on one generic slice.
/// A few of the essential `Vecg` methods are reimplemented here
/// to mutate `self`. This is for efficiency and convenience.
/// For example, in vector iterative methods.
pub trait MutVecg {
    /// mutable multiplication by a scalar
    fn mutsmult<U>(self, _s: U)
    where
        U: Copy + PartialOrd,
        f64: From<U>;
    /// mutable vector subtraction
    fn mutvsub<U>(self, _v: &[U])
    where
        U: Copy + PartialOrd,
        f64: From<U>;
    /// mutable vector addition
    fn mutvadd<U>(self, _v: &[U])
    where
        U: Copy + PartialOrd,
        f64: From<U>;

    /// Invert the magnitude
    fn minvert(self);
    /// Negate the vector (all components swap sign)
    fn mneg(self);
    /// Make into a unit vector
    fn munit(self);
    /// Linearly transform to interval [0,1]
    fn mlintrans(self);
}

/// Methods specialised to and more efficient, for `&[u8]`
pub trait Vecu8 {
    /// Probability density function of bytes data
    fn pdfu8(self) -> Vec<f64>;
    /// Information (entropy) of &[u8] (in nats)
    fn entropyu8(self) -> f64;
    /// Joint probability density function (here just co-occurence counts)
    /// of paired values in two vectors of bytes of the same length.
    /// Needs n^2 x 32bits of memory.
    fn jointpdfu8(self, v: &[u8]) -> Result<Vec<Vec<u32>>, RE>;
    /// Joint entropy of &[u8],&[u8] (in nats)
    fn jointentropyu8(self, v: &[u8]) -> Result<f64, RE>;
    /// Statistical pairwise dependence of two &[u8] variables in the range [0,1]
    fn dependenceu8(self, v: &[u8]) -> Result<f64, RE>;
    /// Independence in the range [1,2] of two &[u8] variables
    fn independenceu8(self, v: &[u8]) -> Result<f64, RE>;
}

/// Methods applicable to a slice of vectors of generic end type.
/// Operations on a whole set of multidimensional vectors.
pub trait VecVec<T> {
    /// Selects a column by number
    fn column(self, cnum: usize) -> Vec<f64>;
    /// Transpose slice of vecs matrix
    fn transpose(self) -> Vec<Vec<f64>>;
    /// Normalize columns, so that they are all unit vectors
    fn normalize(self) -> Result<Vec<Vec<f64>>,RE>;
    /// Householder's method returning matrices (U,R)
    fn house_ur(self) -> Result<(TriangMat, TriangMat),RE>; 
    /// Joint probability density function of n matched slices of the same length
    fn jointpdfn(self) -> Result<Vec<f64>, RE>;
    /// Joint entropy between a set of vectors of the same length
    fn jointentropyn(self) -> Result<f64, RE>;
    /// Independence (component wise) of a set of vectors.
    fn dependencen(self) -> Result<f64, RE>;
    /// Flattened lower triangular relations matrix between columns of self
    fn crossfeatures(self, f: fn(&[T], &[T]) -> f64) -> Result<TriangMat, RE>;
    /// Sum of nd points (or vectors)
    fn sumv(self) -> Vec<f64>;
    /// Arithmetic Centre = euclidian mean of a set of points
    fn acentroid(self) -> Vec<f64>;
    /// Multithreaded Arithmetic Centre = euclidian mean of a set of points
    fn par_acentroid(self) -> Vec<f64>;
    /// Geometric Centroid
    fn gcentroid(self) -> Result<Vec<f64>,RE>;
    /// Harmonic Centroid = harmonic mean of a set of points
    fn hcentroid(self) -> Result<Vec<f64>,RE>;
    /// Possible first iteration point for geometric medians
    fn firstpoint(self) -> Vec<f64>;
    /// Sums of distances from each point to all other points
    fn distsums(self) -> Vec<f64>;
    /// Medoid distance, its index, outlier distance, its index
    fn medout(self, gm: &[f64]) -> MinMax<f64>;
    /// Like gmparts, except only does one iteration from any non-member point g
    fn nxnonmember(self, g: &[f64]) -> (Vec<f64>, Vec<f64>, f64);
    /// Radius of a point specified by its subscript.    
    fn radius(self, i: usize, gm: &[f64]) -> Result<f64, RE>;
    /// Exact radii vectors to each member point by using the gm
    fn radii(self, gm: &[f64]) -> Vec<f64>;
    /// Arith mean and std (in MStats struct), Median and mad, Medoid and Outlier (in MinMax struct)
    fn eccinfo(self, gm: &[f64]) -> Result<(MStats, MStats, MinMax<f64>), RE>
    where
        Vec<f64>: FromIterator<f64>;
    /// Quasi median, recommended only for comparison purposes
    fn quasimedian(self) -> Result<Vec<f64>, RE>;
    /// Geometric median estimate's error
    fn gmerror(self, gm: &[f64]) -> f64;
    /// Proportions of points in idx along each +/-axis (hemisphere)
    /// Self will normally be zero median vectors
    fn tukeyvec(self, idx: &[usize]) -> Result<Vec<f64>, RE>;
    /// MADGM, absolute deviations from geometric median: stable nd data spread estimator
    fn madgm(self, gm: &[f64]) -> Result<f64, RE>;
    /// Collects indices of outer and inner hull points, from zero median data
    fn hulls(self) -> (Vec<usize>, Vec<usize>);
    /// New algorithm for geometric median, to accuracy eps    
    fn gmedian(self, eps: f64) -> Vec<f64>;
    /// Parallel (multithreaded) implementation of Geometric Median. Possibly the fastest you will find.
    fn par_gmedian(self, eps: f64) -> Vec<f64>;
    /// Like `gmedian` but returns the sum of unit vecs and the sum of reciprocals of distances.
    fn gmparts(self, eps: f64) -> (Vec<f64>, Vec<f64>, f64);
}

/// Methods applicable to slice of vectors of generic end type, plus one other argument
/// of a similar kind
pub trait VecVecg<T, U> {
    /// Leftmultiply (column) vector v by (rows of) matrix self
    fn leftmultv(self, v: &[U]) -> Result<Vec<f64>, RE>;
    /// Rightmultiply (row) vector v by (columns of) matrix self
    fn rightmultv(self, v: &[U]) -> Result<Vec<f64>, RE>;
    /// Matrix multiplication self * m
    fn matmult(self, m: &[Vec<U>]) -> Result<Vec<Vec<f64>>, RE>;
    /// Weighted sum of nd points (or vectors)
    fn wsumv(self, ws: &[U]) -> Vec<f64>;
    /// Weighted Arithmetic Centre = weighted euclidian mean of a set of points
    fn wacentroid(self, ws: &[U]) -> Vec<f64>;
    /// Trend between two sets
    fn trend(self, eps: f64, v: Vec<Vec<U>>) -> Result<Vec<f64>, RE>;
    /// Subtract m from all points - e.g. transform to zero median form
    fn translate(self, m: &[U]) -> Result<Vec<Vec<f64>>, RE>;
    /// Proportions of points along each +/-axis (hemisphere)
    fn wtukeyvec(self, idx: &[usize], ws: &[U]) -> Result<Vec<f64>, RE>;
    /// Dependencies of vector m on each vector in self
    fn dependencies(self, m: &[U]) -> Result<Vec<f64>, RE>;
    /// (Median) correlations of m with each vector in self
    fn correlations(self, m: &[U]) -> Result<Vec<f64>, RE>;
    /// Sum of distances from arbitrary point (v) to all the points in self      
    fn distsum(self, v: &[U]) -> Result<f64, RE>;
    /// Individual distances from any point v (typically not in self) to all the points in self.    
    fn dists(self, v: &[U]) -> Result<Vec<f64>, RE>;
    /// Weighted sorted weighted radii magnitudes, normalised
    fn wsortedrads(self, ws: &[U], gm: &[f64]) -> Result<Vec<f64>, RE>;
    /// Like wgmparts, except only does one iteration from any non-member point g
    fn wnxnonmember(self, ws: &[U], g: &[f64]) -> Result<(Vec<f64>, Vec<f64>, f64), RE>;
    /// The weighted geometric median to accuracy eps
    fn wgmedian(self, ws: &[U], eps: f64) -> Result<Vec<f64>, RE>;
    /// Parallel (multithreaded) implementation of the weighted Geometric Median.  
    fn par_wgmedian(self, ws: &[U], eps: f64) -> Result<Vec<f64>, RE>;
    /// Like `wgmedian` but returns also the sum of unit vecs and the sum of reciprocals.
    fn wgmparts(self, ws: &[U], eps: f64) -> Result<(Vec<f64>, Vec<f64>, f64), RE>;
    /// wmadgm median of weighted absolute deviations from (weighted) gm: stable nd data spread estimator
    fn wmadgm(self, ws: &[U], wgm: &[f64]) -> Result<f64, RE>;
    /// Flattened lower triangular part of a covariance matrix of a Vec of f64 vectors.
    fn covar(self, med: &[U]) -> Result<TriangMat, RE>;
    /// Flattened lower triangular part of a covariance matrix for weighted f64 vectors.
    fn wcovar(self, ws: &[U], m: &[f64]) -> Result<TriangMat, RE>;
}
