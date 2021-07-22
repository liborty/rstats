pub mod statsg;
pub mod mutstats;
pub mod vecg;
pub mod vecu8;
pub mod vecvecu8;
pub mod mutvec;
pub mod vecvec;
pub mod functions;

// reexporting to avoid duplication and for backwards compatibility
pub use indxvec::{here,wi,wv}; 
/// simple error handling
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
            "median:\n\tLower Q: {}\n\tMedian:  {}\n\tUpper Q: {}",
            wi(&self.lquartile),
            wi(&self.median),
            wi(&self.uquartile)
        )
    }
}

/// Mean and standard deviation (or std ratio for geometric mean)
#[derive(Default)]
pub struct MStats {
    pub mean: f64,
    pub std: f64,
}
impl std::fmt::Display for MStats {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "mean±std: {}±{}", wi(&self.mean), wi(&self.std))
    }
}

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
    /// Vector with inverse magnitude
    fn vinverse(self) -> Vec<f64>;
    // negated vector (all components swap sign)
    fn negv(self) -> Vec<f64>;
    /// Unit vector
    fn vunit(self) -> Vec<f64>;
    
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
    /// Probability density function of a sorted slice
    fn pdf(self) -> Vec<f64>;
    /// Information (entropy) in nats
    fn entropy(self) -> f64;
    /// (Auto)correlation coefficient of pairs of successive values of (time series) variable.
    fn autocorr(self) -> f64; 
    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Vec<f64>;
    /// Reconstructs the full symmetric square matrix from its lower diagonal compact form,
    /// produced by covar, covone, wcovar
    fn symmatrix(self) -> Vec<Vec<f64>>;
   }

/// A few of the `Stats` methods are reimplemented here
/// (only for f64), so that they mutate `self` in-place.
/// This is more efficient and convenient in some circumstances.
pub trait MutStats {
    /// Invert the magnitude
    fn minvert(self); 
    // negate vector (all components swap sign)
    fn mneg(self);
    /// Make into a unit vector
    fn munit(self);
    /// Linearly transform to interval [0,1]
    fn mlintrans(self);
    /// Sort in place  
    fn msortf(self); 
}

/// Vector Algebra on two vectors (represented here as generic slices).
pub trait Vecg<_T,U> {
    /// Scalar multiplication of a vector
    fn smult(self, s:U) -> Vec<f64>;
    /// Scalar multiplication by f64
    fn smultf64(self, s:f64) -> Vec<f64>;
    /// Scalar addition to vector
    fn sadd(self, s:U) -> Vec<f64>;
    /// Scalar product 
    fn dotp(self, v:&[U]) -> f64;
    /// Cosine of angle between two slices
    fn cosine(self, v:&[U]) -> f64;
    /// Vectors' subtraction (difference)
    fn vsub(self, v:&[U]) -> Vec<f64>;
    /// Vector subtraction of `&[f64]`
    fn vsubf64(self, v:&[f64]) -> Vec<f64>;
    /// Vectors difference as unit vector (done together for efficiency)
    fn vsubunit(self, v: &[U]) -> Vec<f64>; 
    /// Vector addition
    fn vadd(self, v:&[U]) -> Vec<f64>; 
    /// Adding `&[f64]` 
    fn vaddf64(self, v:&[f64]) -> Vec<f64>; 
    /// Euclidian distance 
    fn vdist(self, v:&[U]) -> f64;
    /// Euclidian distance to `&[f64]`
    fn vdistf64(self, v:&[f64]) -> f64;
    /// Euclidian distance squared
    fn vdistsq(self, v:&[U]) -> f64; 
    /// cityblock distance
    fn cityblockd(self, v:&[U]) -> f64;
    /// Magnitude of the cross product |a x b| = |a||b|sin(theta).
    /// Attains maximum `|a|.|b|` when the vectors are othogonal.
    fn varea(self, v:&[U]) -> f64;
    /// Area proportional to the swept arc
     fn varc(self, v:&[U]) -> f64; 
    /// Vector similarity S in the interval [0,1]: S = (1+cos(theta))/2
    fn vsim(self, v:&[U]) -> f64;
    /// We define vector dissimilarity D in the interval [0,1]: D = 1-S = (1-cos(theta))/2
    fn vdisim(self, v:&[U]) -> f64;
    /// Lower triangular part of a covariance matrix for a single f64 vector.
    fn covone(self, m:&[U]) -> Vec<f64>;

    /// Joint entropy of &[u8],&[u8] in nats 
    fn jointentropy(self, v:&[U]) -> f64;
    /// Statistical independence measure based on joint entropy
    fn dependence(self, v:&[U]) -> f64; 

    /// Pearson's correlation.  
    fn correlation(self, v:&[U]) -> f64; 
    /// Kendall Tau-B correlation. 
    fn kendalcorr(self, v:&[U]) -> f64;
    /// Spearman rho correlation.
    fn spearmancorr(self, v:&[U]) -> f64;   
}

/// Mutable vector operations that take one generic argument. 
/// A few of the essential `Vecg` methods are reimplemented here 
/// to mutate `self` in-place (only for f64). 
/// This is for efficiency and convenience, for example, in
/// vector iterative methods.
pub trait MutVecg<U> {
    /// mutable multiplication by a scalar
    fn msmult(self, _s:U);  
    /// mutable vector subtraction
    fn mvsub(self, _v: &[U]); 
    /// mutable vector addition
    fn mvadd(self, _v: &[U]);  
}
/// Vector mutating operations that take one argument of f64 end type.
pub trait MutVecf64 {
    /// mutable multiplication by a scalar
    fn mutsmult(self, _s:f64); 
    /// mutable vector subtraction 
    fn mutvsub(self, _v: &[f64]);
    /// mutable vector addition
    fn mutvadd(self, _v: &[f64]); 
}

/// Methods specialised to, or more efficient for `&[u8]`
pub trait Vecu8 {
    /// Scalar product of two (positive) u8 slices.   
    /// Must be of the same length - no error checking (for speed)
    fn dotpu8(self, v: &[u8]) -> u64; 
    /// Cosine between two (positive) u8 slices.
    fn cosineu8(self, v: &[u8]) -> f64; 
    /// Vector subtraction (converts results to f64 as they can be negative)
    fn vsubu8(self, v: &[u8]) -> Vec<f64>;
    /// Vector addition ( converts results to f64, as they can exceed 255 )
    fn vaddu8(self, v: &[u8]) -> Vec<f64>;
    /// Euclidian distance between self &[u8] and v:&[u8].  
    /// Faster than vsub followed by vmag, as both are done in one loop
    fn vdistu8(self, v: &[u8]) -> f64; 
    /// cityblock distance
    fn cityblockdu8(self, v:&[u8]) -> f64;
    ///Euclidian distance squared, the arguments are both of &[u8] type  
    fn vdistsqu8(self, v: &[u8]) -> u64; 
    /// Probability density function of bytes data
    fn pdfu8(self) -> Vec<f64>;
    /// Information (entropy) of &[u8] (in nats)
    fn entropyu8(self) -> f64;
    /// Joint probability density function (here just co-occurence counts) 
    /// of paired values in two vectors of bytes of the same length.
    /// Needs n^2 x 32bits of memory. Do not use for very long vectors, 
    /// those need hashing implementation.
    fn jointpdfu8(self, v:&[u8]) -> Vec<Vec<u32>>;
    /// Joint entropy of &[u8],&[u8] (in nats)
    fn jointentropyu8(self, v:&[u8]) -> f64;
    /// Dependence of two &[u8] variables, the range is [0,1],
    /// i.e. it returns 0 iff they are statistically independent
    /// and 1 when they are identical
    fn dependenceu8(self, v:&[u8]) -> f64;
}

/// A few specialised methods applicable to `Vec<Vec<u8>>` (vector of vectors of bytes).
pub trait VecVecu8 { 
    /// Centroid = euclidian mean of a set of points  
    fn acentroid(self) -> Vec<f64>;
    /// Weighted Centre
    fn wacentroid(self,ws: &[u8]) -> Vec<f64>;
    /// Eccentricity vector added to a non member point,
    fn nxnonmember(self, p:&[f64]) -> Vec<f64>;
    /// Weighted eccentricity vector for a non member point
    fn wnxnonmember(self, ws:&[u8], p:&[f64]) -> Vec<f64>; 
    /// Weighted geometric median, sorted eccentricities magnitudes, cpdf of the weights
    fn gmedian(self, eps:f64) -> Vec<f64>; 
    /// The weighted geometric median
    fn wgmedian(self, ws:&[u8], eps: f64) -> Vec<f64>;
    /// Lower triangular part of a covariance matrix of a Vec of u8 vectors.
    fn covar(self, med:&[f64]) -> Vec<f64>; 
    fn wcovar(self, ws:&[u8], m:&[f64]) -> Vec<f64>;
}

/// Methods applicable to vector of vectors of generic end type
pub trait VecVec<T> {

    /// Arithmetic Centre = euclidian mean of a set of points
    fn acentroid(self) -> Vec<f64>;
    /// Weighted Arithmetic Centre = weighted euclidian mean of a set of points
    fn wacentroid(self,ws: &[T]) -> Vec<f64>;
    /// Geometric Centroid
    fn gcentroid(self) -> Vec<f64>;
    /// Harmonic Centroid = harmonic mean of a set of points
    fn hcentroid(self) -> Vec<f64>; 
    /// Trend between two sets
    fn trend(self, eps: f64, v: Vec<Vec<T>>) -> Vec<f64>;
    /// Subtract m from all points - e.g. transform to zero median form
    fn translate(self, m: &[T]) -> Vec<Vec<f64>>;

    /// Sums of distances from each point to all other points.
     fn distsums(self) -> Vec<f64>;
    /// Fast sums of distances from each point to all other points 
    fn distsuminset(self, indx: usize) -> f64;
    /// Sum of distances from arbitrary point (v) to all the points in self   
    fn distsum(self, v: &[T]) -> f64;
    /// Individual distances from any point v (typically not in self) to all the points in self.    
    fn dists(self, v: &[T]) -> Vec<f64>;
    /// Medoid and Outlier (by distance) of a set of points
    fn medoid(self) -> (f64, usize, f64, usize); 
    /// Eccentricity vectors from each point
    fn eccentricities(self) -> Vec<Vec<f64>>;
    /// Exact eccentricity vectors from all member points by first finding the Geometric Median.
    /// As well as being more accurate, it is usually faster than `eccentricities` above, 
    /// especially for large numbers of points.
    fn exacteccs(self, eps: f64) -> Vec<Vec<f64>>;
    /// Returns ( gm, sorted eccentricities magnitudes )
    fn sortedeccs(self, ascending:bool, eps:f64) -> ( Vec<f64>,Vec<f64> );
    /// ( wgm, sorted eccentricities magnitudes, associated cpdf )
    fn wsortedeccs(self, ws: &[T], eps:f64) -> ( Vec<f64>,Vec<f64>,Vec<f64> ); 
    /// Sorted cosines magnitudes and cpdf, needs central median
    fn wsortedcos(self, medmed: &[T], med: &[T], ws: &[T]) -> ( Vec<f64>,Vec<f64> ); 
    /// Next approx median point from this member point given by its indx
    fn nxmember(self, indx: usize) -> Vec<f64>;
    /// Ecentricity of a member point given by its indx
    fn eccmember(self, indx: usize) -> Vec<f64>;
    /// Next approx median point from this nonmember point
    fn nxnonmember(self, p:&[f64]) -> Vec<f64>;
    /// Error vector, i.e. unscaled eccentricity vector
    fn errorv(self, p:&[T]) -> Vec<f64>;
    /// Eccentricity vector for a non member point
    fn eccnonmember(self, p:&[f64]) -> Vec<f64>; 
    /// Weighted eccentricity vector for a non member point
    fn wnxnonmember(self, ws:&[T], p:&[f64]) -> Vec<f64>; 
    /// magnitudes of a set of vectors
    fn mags(self) -> Vec<f64>; 
    /// Median and quartiles of eccentricities (new robust measure of spread of a multivariate sample)
    fn moe(self, eps: f64) -> (MStats,Med);
    /// Medoid and Outlier as defined by eccentricities.
    fn emedoid(self, eps: f64) -> (f64, usize, f64, usize);

    /// Geometric medians of a set
    /// First iteration point for geometric medians
    fn firstpoint(self) -> Vec<f64>;
    /// Improved Weizsfeld's Algorithm for geometric median
    fn nmedian(self, eps: f64) -> Vec<f64>;
    /// New secant algorithm for geometric median
    fn gmedian(self, eps: f64) -> Vec<f64>; 
    /// The weighted geometric median
    fn wgmedian(self, ws: &[T],eps: f64) -> Vec<f64>; 
    /// Flattened lower triangular part of a covariance matrix of a Vec of f64 vectors.
    fn covar(self, med:&[f64]) -> Vec<f64>; 
    /// Flattened lower triangular part of a comediance matrix of a Vec of f64 vectors.
    fn comed(self, m:&[f64], eps:f64) -> Vec<f64>;
    /// Flattened lower triangular part of a covariance matrix for weighted f64 vectors.
    fn wcovar(self, ws:&[f64], m:&[f64]) -> Vec<f64>;
    /// Flattened lower triangular part of a comediance matrix for weighted f64 vectors.
    fn wcomed(self, ws:&[f64], m:&[f64], eps:f64) -> Vec<f64>;
}
