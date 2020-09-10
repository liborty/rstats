pub mod i64impls;
pub mod f64impls;
pub mod vimpls;

use std::cmp::Ordering::Equal;
use anyhow::Result;

/// Median and quartiles
#[derive(Default)]
pub struct Med {
    pub lquartile: f64,
    pub median: f64,
    pub uquartile: f64
}
impl std::fmt::Display for Med {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f,"Lower Quartile: {}, Median: {}, Upper Qartile: {}", 
        green(self.lquartile), green(self.median), green(self.uquartile))
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
      write!(f,"mean±std: {}±{}",green(self.mean),green(self.std))
   }
}
/// Implementing basic statistical measures.
pub trait RStats { 
   /// Arithmetic mean
   fn amean(&self) -> Result<f64>;
   /// Arithmetic mean and standard deviation
   fn ameanstd(&self) -> Result<MStats>;
   /// Weighted arithmetic mean
   fn awmean(&self) -> Result<f64>;
   /// Weighted arithmetic men and standard deviation
   fn awmeanstd(&self) -> Result<MStats>;
   /// Harmonic mean
   fn hmean(&self) -> Result<f64>;
   /// Weighted harmonic mean
   fn hwmean(&self) -> Result<f64>;
   /// Geometric mean
   fn gmean(&self) -> Result<f64>;
   /// Geometric mean and stndard deviation ratio
   fn gmeanstd(&self) -> Result<MStats>;
   /// Weighed geometric mean
   fn gwmean(&self) -> Result<f64>;
   /// Weighted geometric mean and standard deviation ratio
   fn gwmeanstd(&self) -> Result<MStats>;
   /// Median and quartiles
   fn median(&self) -> Result<Med>;
   /// Correlation of i64 variables
   fn icorrelation(&self, v:&[i64]) -> Result<f64>;
   /// Correlation of f64 variables
   fn correlation(&self, v:&[f64]) -> Result<f64>;
   /// Autocorrelation
   fn autocorr(&self) -> Result<f64>;
  
}

/// Mutable primitive vector operations (for efficiency).
pub trait MutVectors {
   /// mutable multiplication by a scalar
   fn mutsmult(&mut self, s:f64);
   /// mutable vector subtraction
   fn mutvsub(&mut self, v:&[f64]);
   /// mutable vector addition
   fn mutvadd(&mut self, v:&[f64]);
   /// mutably makes into a unit vector
   fn mutvunit(&mut self);
   /// magnitude of a mutable vector (vector unchanged)
   fn mutvmag(&mut self) -> f64;

   }

/// Implementing basic vector algebra and safe geometric median.
pub trait Vectors {
   /// Scalar product of two vectors
   fn dotp(&self, v:&[f64]) -> f64;
   /// Vector subtraction
   fn vsub(&self, v:&[f64]) -> Vec<f64>;
   /// Vector addition
   fn vadd(&self, v:&[f64]) -> Vec<f64>;
   /// Vector magnitude
   fn vmag(&self) -> f64;
   /// Euclidian distance between two points
   fn vdist(&self, v:&[f64]) -> f64;
   /// Scalar multiplication
   fn smult(&self, s:f64) -> Vec<f64>;
   /// Unit vector
   fn vunit(&self) -> Vec<f64>;

   /// Centroid = euclidian mean of a set of points
   fn acentroid(&self, d:usize) -> Vec<f64>;
   /// Medoid of a set of points (most central of the points)
   fn medoid(&self, d:usize) -> Result<(f64,usize,f64,usize)>;
   /// Sum of distances from all the points in a set to v
   fn distsum(&self, d:usize, v:&[f64] ) -> f64;
   /// Ecentricity measure (0,1) of an internal point given by indx, w.r.t. the set
   fn eccentr(&self, d:usize, indx:usize) -> f64;
   /// Eccentricity measure (0,1) and vector of any point, w.r.t. the set   
   fn veccentr(&self, d:usize, thisp:&[f64]) -> Result<(f64,Vec<f64>)>;
   /// Eccentricity vecor of any point w.r.t. the set
   fn ecc(&self, d:usize, v:&[f64]) -> f64;     
   /// Geometric median of the set
   fn nmedian(&self, d:usize, eps:f64) -> Result<Vec<f64>>;
 
}

/// private helper function for formatting error messages
fn emsg(file:&'static str, line:u32, msg:&'static str)-> String {
   format!("{}:{} rstats {}",file,line,msg)
}

/// Private sum of linear weights 
fn wsum(n: usize) -> f64 { (n*(n+1)) as f64/2. }


/// green terminal text output to make test results stand out.
pub fn green(n:f64) -> String { format!("\x1B[01;92m{}\x1B[0m",n) }


/// Sorts a mutable `Vec<f64>` in place.  
/// It is the responsibility of the user to ensure that there are no NaNs etc.
pub fn sortf(v: &mut [f64]) { 
   v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
}

/// Generates a random f64 vector of size d x n suitable for testing. It needs two seeds.  
/// Uses local closure `rand` to generate random numbers (avoids dependencies).  
/// Random numbers are in the open interval 0..1 with uniform distribution.  
pub fn genvec(d:usize, n:usize, s1:u32, s2:u32 ) -> Vec<f64> {
   let size = d*n;
   // change the seeds as desired
   let mut m_z = s1 as u32;
   let mut m_w = s2 as u32;
   let mut rand = || {
      m_z = 36969 * (m_z & 65535) + (m_z >> 16);
      m_w = 18000 * (m_w & 65535) + (m_w >> 16);
      (((m_z << 16) & m_w) as f64 + 1.0)*2.328306435454494e-10
   };
   let mut v = Vec::with_capacity(size); 
   for _i in 0..size { v.push(rand()) }; // fills the lot with random numbers
   return v
}
