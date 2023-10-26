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

pub use crate::error::{re_error, RError, RE};
// reexporting useful related methods
pub use indxvec::{here, printing::*, MinMax, Printing};
pub use medians::{error::MedError, Median, Medianf64};

// Auxiliary Functions

/// Convenience dummy function for quantify closure
pub fn noop(f: &f64) -> f64 { *f }

/// Convenience From quantification invocation
pub fn fromop<T: Clone + Into<f64>>(f: &T) -> f64 {
    (*f).clone().into()
}

/// tm_statistic in 1d: (value-centre)/dispersion
/// generalized to any measure of central tendency and dispersion
pub fn tm_stat(val: f64, params: Params) -> f64 {
    (val - params.centre) / params.spread
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

/// Holds measures of central tendency and spread.
/// Usually some kind of mean and its associated standard deviation, or median and its MAD
#[derive(Default)]
pub struct Params {
    /// central tendency - (geometric, arithmetic, harmonic means or median)
    pub centre: f64,
    /// measure of data spread, typically standard deviation or MAD
    pub spread: f64,
}
impl std::fmt::Display for Params {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{YL}centre: {GR}{:.5}{YL} ± spread: {GR}{:.5}{UN}",
            self.centre, self.spread
        )
    }
}

/// Compact Triangular Matrix.  
/// TriangMat stores and manipulates triangular matrices.
/// It also compactly represents square symmetric and antisymmetric matrices.
/// TriangMat is typically produced by some matrix calculations,
/// so the end-type for its data is f64.  
/// The actual length of its data is `triangmat.len() = n*(n+1)/2`.  
/// The dimension of the implied nxn matrix is `n = triangmat.dim()`.  
/// The kind (field) of the TriangMat is encoded as 0..5. This enables
/// trivial transpositions by: `(kind+3) % 6`.
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
    /// Median and Mad packed into Params struct
    fn medmad(self) -> Result<Params, RE>;
    /// Arithmetic mean and standard deviation
    fn ameanstd(self) -> Result<Params, RE>;
    /// Weighted arithmetic mean
    fn awmean(self) -> Result<f64, RE>;
    /// Weighted arithmetic men and standard deviation
    fn awmeanstd(self) -> Result<Params, RE>;
    /// Harmonic mean
    fn hmean(self) -> Result<f64, RE>;
    /// Harmonic mean and experimental standard deviation
    fn hmeanstd(self) -> Result<Params, RE>;
    /// Weighted harmonic mean
    fn hwmean(self) -> Result<f64, RE>;
    /// Weighted harmonic mean and standard deviation
    fn hwmeanstd(self) -> Result<Params, RE>;
    /// Geometric mean
    fn gmean(self) -> Result<f64, RE>;
    /// Geometric mean and standard deviation ratio
    fn gmeanstd(self) -> Result<Params, RE>;
    /// Weighed geometric mean
    fn gwmean(self) -> Result<f64, RE>;
    /// Weighted geometric mean and standard deviation ratio
    fn gwmeanstd(self) -> Result<Params, RE>;
    /// Probability density function of a sorted slice
    fn pdf(self) -> Vec<f64>;
    /// Information (entropy) in nats
    fn entropy(self) -> f64;
    /// (Auto)correlation coefficient of pairs of successive values of (time series) variable.
    fn autocorr(self) -> Result<f64, RE>;
    /// Linear transform to interval `[0,1]`
    fn lintrans(self) -> Result<Vec<f64>, RE>;
    /// Linearly weighted approximate time series derivative at the last (present time) point.
    fn dfdt(self, centre: f64) -> Result<f64, RE>;
    /// Householder reflection
    fn house_reflector(self) -> Vec<f64>;
}

/// Vector Algebra on two vectors (represented here as generic slices).
/// Also included are scalar operations on the `self` vector.
pub trait Vecg {
    /// nd tm_statistic of self against centre (geometric median) spread (madgm)
    fn tm_statistic(self, centre: &[f64], spread: f64) -> Result<f64, RE>;
    /// Dot product of vector self with column c of matrix v
    fn columnp<U: Clone + Into<f64>>(self, c: usize, v: &[Vec<U>]) -> f64;
    /// Scalar addition to vector
    fn sadd<U: Into<f64>>(self, s: U) -> Vec<f64>;
    /// Scalar multiplication of vector, creates new vec
    fn smult<U: Into<f64>>(self, s: U) -> Vec<f64>;
    /// Scalar product
    fn dotp<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Sigvec product with (zero median) vector self. Cloud density d in its direction:`0 <=  d <= |self|`
    fn dotsig(self, sig: &[f64]) -> Result<f64, RE>;
    /// Cosine of angle between two slices
    fn cosine<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Sine of an angle with correct sign in any number of dimensions
    fn sine<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Vectors' subtraction (difference)
    fn vsub<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<f64>;
    /// Vectors difference as unit vector (done together for efficiency)
    fn vsubunit<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<f64>;
    /// Vector addition
    fn vadd<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<f64>;
    /// Euclidian distance
    fn vdist<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Weighted arithmetic mean of `self:&[T]`, scaled by `ws:&[U]`
    fn wvmean<U: Clone + Into<f64>>(self, ws: &[U]) -> f64;
    /// Weighted distance of `self:&[T]` to `v:&[V]`, weighted by `ws:&[U]`
    fn wvdist<U: Clone + Into<f64>, V: Clone + Into<f64>>(self, ws: &[U], v: &[V]) -> f64;
    /// Euclidian distance squared
    fn vdistsq<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// cityblock distance
    fn cityblockd<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Area spanned by two vectors always over their concave angle
    fn varea<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> f64;
    /// Area proportional to the swept arc
    fn varc<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> f64;
    /// Vector similarity S in the interval `[0,1]: S = (1+cos(theta))/2`
    fn vsim<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Median correlation similarity of vectors, range [0,1]
    fn vcorrsim(self, v:Self) -> Result<f64, RE>;
    /// Positive dotp [0,2|a||b|]
    fn pdotp<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> f64;
    /// We define vector dissimilarity D in the interval `[0,1]: D = 1-S = (1-cos(theta))/2`
    fn vdisim<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Lower triangular part of a covariance matrix for a single f64 vector.
    fn covone<U: Clone + Into<f64>>(self, m: &[U]) -> TriangMat;
    /// Kronecker product of two vectors
    fn kron<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<f64>;
    /// Outer product: matrix multiplication of column self with row v.
    fn outer<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<Vec<f64>>;
    /// Exterior (Grassman) algebra product: produces 2-blade components
    fn wedge<U: Clone + Into<f64>>(self, v: &[U]) -> TriangMat;
    /// Geometric (Clifford) algebra product in matrix form: **a*b** + **a∧b**
    fn geometric<U: Clone + Into<f64>>(self, b: &[U]) -> TriangMat;
    /// Joint probability density function
    fn jointpdf<U: Clone + Into<f64>>(self, v: &[U]) -> Result<Vec<f64>, RE>;
    /// Joint entropy of `&[T],&[U]` in nats (units of e)
    fn jointentropy<U: Clone + Into<f64>>(self, v: &[U]) -> Result<f64, RE>;
    /// Statistical pairwise dependence of `&[T] &[U]` variables in the range `[0,1]`
    fn dependence<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> Result<f64, RE>;
    /// Statistical pairwise independence in the range `[0,1]` based on joint entropy
    fn independence<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> Result<f64, RE>;
    /// Pearson's correlation.  
    fn correlation<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Kendall Tau-B correlation.
    fn kendalcorr<U: Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Spearman rho correlation.
    fn spearmancorr<U: PartialOrd + Clone + Into<f64>>(self, v: &[U]) -> f64;
    /// Change to gm that adding point self will cause
    fn contribvec_newpt(self, gm: &[f64], recips: f64) -> Result<Vec<f64>,RE>;
    /// Normalized magnitude of change to gm that adding point self will cause
    fn contrib_newpt(self, gm: &[f64], recips: f64, nf: f64) -> Result<f64,RE>;
    /// Contribution of removing point self
    fn contribvec_oldpt(self, gm: &[f64], recips: f64) -> Result<Vec<f64>,RE>;
    /// Normalized contribution of removing point self (as negative scalar)
    fn contrib_oldpt(self, gm: &[f64], recips: f64, nf: f64) -> Result<f64,RE>;
    /// Householder reflect
    fn house_reflect<U: Clone + PartialOrd + Into<f64>>(self, x: &[U]) -> Vec<f64>;
}

/// Mutable operations on one generic slice.
/// A few of the essential `Vecg` methods are reimplemented here
/// to mutate `self`. This is for efficiency and convenience.
/// For example, in vector iterative methods.
pub trait MutVecg {
    /// mutable multiplication by a scalar
    fn mutsmult<U: PartialOrd + Into<f64>>(self, _s: U);
    /// mutable vector subtraction
    fn mutvsub<U: Clone + PartialOrd + Into<f64>>(self, _v: &[U]);
    /// mutable vector addition
    fn mutvadd<U: Clone + PartialOrd + Into<f64>>(self, _v: &[U]);
    /// Invert the magnitude
    fn minvert(self);
    /// Negate the vector (all components swap sign)
    fn mneg(self);
    /// Make into a unit vector
    fn munit(self);
    /// Linearly transform to interval `[0,1]`
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
    /// Statistical pairwise dependence of two `&[u8]` variables in the range `[0,1]`
    fn dependenceu8(self, v: &[u8]) -> Result<f64, RE>;
    /// Independence in the range `[1,2]` of two `&[u8]` variables
    fn independenceu8(self, v: &[u8]) -> Result<f64, RE>;
}

/// Methods applicable to a slice of vectors of generic end type.
/// Operations on a whole set of multidimensional vectors.
pub trait VecVec<T> {
    /// Linearly weighted approximate time series derivative at the last point (present time).
    fn dvdt(self, centre: &[f64]) -> Result<Vec<f64>, RE>; 
    /// Maps a scalar valued closure onto all vectors in self
    fn scalar_fn(self, f: impl Fn(&[T]) -> Result<f64, RE>) -> Result<Vec<f64>, RE>;
    /// Maps vector valued closure onto all vectors in self and collects
    fn vector_fn(self, f: impl Fn(&[T]) -> Result<Vec<f64>, RE>) -> Result<Vec<Vec<f64>>, RE>;
    /// Exact radii magnitudes to all member points from the Geometric Median.
    fn radii(self, gm: &[f64]) -> Result<Vec<f64>, RE>;
    /// Selects a column by number
    fn column(self, cnum: usize) -> Vec<f64>;
    /// Transpose slice of vecs matrix
    fn transpose(self) -> Vec<Vec<f64>>;
    /// Normalize columns, so that they are all unit vectors
    fn normalize(self) -> Result<Vec<Vec<f64>>, RE>;
    /// Householder's method returning matrices (U,R)
    fn house_ur(self) -> Result<(TriangMat, TriangMat), RE>;
    /// Joint probability density function of n matched slices of the same length
    fn jointpdfn(self) -> Result<Vec<f64>, RE>;
    /// Joint entropy between a set of vectors of the same length
    fn jointentropyn(self) -> Result<f64, RE>;
    /// Dependence (component wise) of a set of vectors.
    fn dependencen(self) -> Result<f64, RE>;
    /// Binary relations between columns of self
    fn crossfeatures(self, f: fn(&[T], &[T]) -> f64) -> Result<TriangMat, RE>;
    /// Sum of nd points (or vectors)
    fn sumv(self) -> Vec<f64>;
    /// Arithmetic Centre = euclidian mean of a set of points
    fn acentroid(self) -> Vec<f64>;
    /// Multithreaded Arithmetic Centre = euclidian mean of a set of points
    fn par_acentroid(self) -> Vec<f64>;
    /// Geometric Centroid
    fn gcentroid(self) -> Result<Vec<f64>, RE>;
    /// Harmonic Centroid = harmonic mean of a set of points
    fn hcentroid(self) -> Result<Vec<f64>, RE>;
    /// Possible first iteration point for geometric medians
    fn firstpoint(self) -> Vec<f64>;
    /// Sums of distances from each point to all other points
    fn distsums(self) -> Vec<f64>;
    /// Medoid distance, its index, outlier distance, its index
    fn medout(self, gm: &[f64]) -> Result<MinMax<f64>, RE>;
    /// Like gmparts, except only does one iteration from any non-member point g
    fn nxnonmember(self, g: &[f64]) -> (Vec<f64>, Vec<f64>, f64);
    /// Radius of a point specified by its subscript.    
    fn radius(self, i: usize, gm: &[f64]) -> Result<f64, RE>;
    /// Arith mean and std (in Params struct), Median and mad, Medoid and Outlier (in MinMax struct)
    fn eccinfo(self, gm: &[f64]) -> Result<(Params, Params, MinMax<f64>), RE>
    where
        Vec<f64>: FromIterator<f64>;
    /// Quasi median, recommended only for comparison purposes
    fn quasimedian(self) -> Result<Vec<f64>, RE>;
    /// Geometric median estimate's error
    fn gmerror(self, gm: &[f64]) -> f64;
    /// Proportional projections on each +/- axis (by hemispheres)
    fn sigvec(self, idx: &[usize]) -> Result<Vec<f64>, RE>;
    /// madgm, median of radii from geometric median: stable nd data spread estimator
    fn madgm(self, gm: &[f64]) -> Result<f64, RE>;
    /// stdgm mean of radii from gm: nd data spread estimator
    fn stdgm(self, gm: &[f64]) -> Result<f64, RE>;
    /// Outer hull subscripts from their square radii and their sort index. 
    fn outer_hull(self, sqrads: &[f64], sindex: &[usize]) -> Vec<usize>; 
    /// Inner hull subscripts from their square radii and their sort index.  
    fn inner_hull(self, sqrads: &[f64], sindex: &[usize]) -> Vec<usize>; 
    /// Measure of likelihood of zero median point **p** belonging to zero median data cloud `self`.
    fn depth(self, descending_index: &[usize], p: &[f64]) -> Result<f64,RE>; 
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
    /// Applies scalar valued closure to all vectors in self and multiplies by their weights
    fn scalar_wfn(
        self,
        ws: &[U],
        f: impl Fn(&[T]) -> Result<f64, RE>,
    ) -> Result<(Vec<f64>, f64), RE>;
    /// Applies vector valued closure to all vectors in self and multiplies by their weights
    fn vector_wfn(
        self,
        v: &[U],
        f: impl Fn(&[T]) -> Result<Vec<f64>, RE>,
    ) -> Result<(Vec<Vec<f64>>, f64), RE>;
    /// Individually weighted time series derivative of vectors
    fn wdvdt(self, ws: &[U], centre: &[f64]) -> Result<Vec<f64>, RE>;
    /// 1.0-dotproducts with **v**, in range [0,2]
    fn divs(self, v: &[U]) -> Result<Vec<f64>, RE>;
    /// weighted 1.0-dotproduct of **v**, with all in self
    fn wdivs(self, ws: &[U], v: &[f64]) -> Result<(Vec<f64>, f64), RE>;
    /// median of weighted cos deviations from **v**
    fn wdivsmed(self, ws: &[U], v: &[f64]) -> Result<f64,RE>;
    /// weighted radii to all points in self
    fn wradii(self, ws:&[U], gm: &[f64]) -> Result<(Vec<f64>,f64),RE>;
    /// wmadgm median of weighted radii from (weighted) gm: stable nd data spread estimator
    fn wmadgm(self, ws: &[U], wgm: &[f64]) -> Result<f64, RE>;
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
    /// Weighted sums of points in each hemisphere
    fn wsigvec(self, idx: &[usize], ws: &[U]) -> Result<Vec<f64>, RE>;
    /// Weighted likelihood of zero median point **p** belonging to zero median data cloud `self`.
    fn wdepth(self, descending_index: &[usize], ws:&[U], p: &[f64]) -> Result<f64,RE>;
    /// Dependencies of vector m on each vector in self
    fn dependencies(self, m: &[U]) -> Result<Vec<f64>, RE>;
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
    /// Flattened lower triangular part of a covariance matrix of a Vec of f64 vectors.
    fn covar(self, med: &[U]) -> Result<TriangMat, RE>;
    /// Flattened lower triangular part of a covariance matrix for weighted f64 vectors.
    fn wcovar(self, ws: &[U], m: &[f64]) -> Result<TriangMat, RE>;
}
