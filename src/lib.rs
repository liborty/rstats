pub mod i64impls;
pub mod f64impls;
pub mod vimpls;
pub mod functions;

use anyhow::Result;
use crate::functions::GreenIt;


/// Median and quartiles
#[derive(Default)]
pub struct Med {
    pub lquartile: f64,
    pub median: f64,
    pub uquartile: f64
}
impl std::fmt::Display for Med {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f,"(Lower Q: {}, Median: {}, Upper Q: {})", 
        GreenIt(self.lquartile), GreenIt(self.median), GreenIt(self.uquartile))
    }
}

/// Mean and standard deviation (or std ratio for geometric mean).
#[derive(Default)]
pub struct MStats {
    pub mean: f64,
    pub std: f64
}
impl std::fmt::Display for MStats {
   fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
      write!(f,"mean±std: {}±{}",GreenIt(self.mean),GreenIt(self.std))
   }
}

/// Implementing basic statistical measures.
/// All these methods operate on only one vector (of data),
/// so they take no arguments.
pub trait RStats { 
   /// Arithmetic mean
   fn amean(self) -> Result<f64>;
   /// Arithmetic mean and standard deviation
   fn ameanstd(self) -> Result<MStats>;
   /// Weighted arithmetic mean
   fn awmean(self) -> Result<f64>;
   /// Weighted arithmetic men and standard deviation
   fn awmeanstd(self) -> Result<MStats>;
   /// Harmonic mean
   fn hmean(self) -> Result<f64>;
   /// Weighted harmonic mean
   fn hwmean(self) -> Result<f64>;
   /// Geometric mean
   fn gmean(self) -> Result<f64>;
   /// Geometric mean and stndard deviation ratio
   fn gmeanstd(self) -> Result<MStats>;
   /// Weighed geometric mean
   fn gwmean(self) -> Result<f64>;
   /// Weighted geometric mean and standard deviation ratio
   fn gwmeanstd(self) -> Result<MStats>;
   /// Median and quartiles
   fn median(self) -> Result<Med>; 
   /// Creates vector of ranks for values in self 
   fn ranks(self) -> Result<Vec<f64>>;
}

/// Mutable primitive vector operations (for efficiency).
pub trait MutVectors {
   /// mutable multiplication by a scalar
   fn mutsmult(self, s:f64);
   /// mutable vector subtraction
   fn mutvsub(self, v:&[f64]);
   /// mutable vector addition
   fn mutvadd(self, v:&[f64]);
   /// mutably makes into a unit vector
   fn mutvunit(self);
   /// transforms into zero median form  
   fn mutzeromd(self, d:usize, eps:f64);
}

/// Implementing basic vector algebra and safe geometric median.
pub trait Vectors {

   /// Utility method to retrieve a sub-slice from multidimensional flat slice.
   fn point(&self,d:usize,i:usize) -> &[f64];

   /// Scalar product of two vectors
   fn dotp(self, v:&[f64]) -> f64;
   /// Vector subtraction
   fn vsub(self, v:&[f64]) -> Vec<f64>;
   /// Vector addition
   fn vadd(self, v:&[f64]) -> Vec<f64>;
   /// Vector magnitude
   fn vmag(self) -> f64;
   /// Euclidian distance between two points
   fn vdist(self, v:&[f64]) -> f64;
   /// Scalar multiplication
   fn smult(self, s:f64) -> Vec<f64>;
   /// Unit vector
   fn vunit(self) -> Vec<f64>;
  
   /// Correlation
   fn correlation(self, v:&[f64]) -> Result<f64>;
   /// Kendall's tau-b (rank order) correlation
   fn kendalcorr(self, v:&[f64]) -> Result<f64>;
   /// Spearman's rho (rank differences) correlation
   fn spearmancorr(self,v:&[f64]) -> Result<f64>;
   /// Autocorrelation
   fn autocorr(self) -> Result<f64>;
   /// Minimum, minimum's index, maximum, maximum's index.
   fn minmax(self) -> (f64,usize,f64,usize);

   /// Centroid = euclidian mean of a set of points
   fn acentroid(self, d:usize) -> Vec<f64>;

   /// Sums of distances from each point to all other points
   fn distances(self, d:usize) -> Result<Vec <f64>>;
   /// Sum of distances from one point given by indx
   fn distsuminset(self, d:usize, indx:usize) -> f64;
   /// Sum of distances from arbitrary point (v) to all the points in self   
   fn distsum(self, d:usize, v:&[f64] ) -> f64;
   /// Medoid and Outlier (by distance) of a set of points
   fn medoid(self, d:usize) -> (f64,usize,f64,usize);

   /// Eccentricity vectors from each point
   fn eccentricities(self, d:usize) -> Result<Vec<Vec<f64>>>;
   /// Ecentricity scalar measure of an internal point given by indx
   fn eccentrinset(self, d:usize, indx:usize) -> f64;
   /// Eccentricity scalar measure and vector of any point     
   fn veccentr(self, d:usize, thisp:&[f64]) -> Result<(f64,Vec<f64>)>;
   /// Eccentricity scalar measure only, of any point 
   fn ecc(self, d:usize, v:&[f64]) -> f64;
   /// Median and quartiles of eccentricities (new robust measure of spread of a multivariate sample)
   fn moe(self, d:usize) -> Med;
   /// Medoid and Outlier as defined by eccentricities.
   fn emedoid(self, d:usize) -> (f64,usize,f64,usize);
     
   /// Geometric median of the set
   fn nmedian(self, d:usize, eps:f64) -> Result<Vec<f64>>; 
   /// Trend between two sets
   fn trend(self, d:usize, eps:f64, v:&[f64]) -> Vec<f64>;
   /// Transform to zero median form. or subtract any other vector `m`
   fn setsub(self, d:usize, m:&[f64]) -> Vec<f64>;
}