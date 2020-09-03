// pub mod tests;
pub mod i64impls;
pub mod f64impls;
pub mod vimpls;
pub mod tests;

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
        write!(f, "(LQ: {}, M: {}, UQ: {})", self.lquartile, self.median, self.uquartile)
    }
}

/// Mean and standard deviation (or std ratio for geometric mean)
#[derive(Default)]
pub struct MStats {
    pub mean: f64,
    pub std: f64
}
impl std::fmt::Display for MStats {
   fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
      write!(f, "Mean:\t{}\nStd:\t{}", self.mean, self.std)
   }
}
/// Implementing basic statistical measures.
pub trait RStats { 

   fn amean(&self) -> Result<f64>;
   fn ameanstd(&self) -> Result<MStats>;
   fn awmean(&self) -> Result<f64>;
   fn awmeanstd(&self) -> Result<MStats>;
   fn hmean(&self) -> Result<f64>;
   fn hwmean(&self) -> Result<f64>;
   fn gmean(&self) -> Result<f64>;
   fn gwmean(&self) -> Result<f64>;
   fn gmeanstd(&self) -> Result<MStats>;
   fn gwmeanstd(&self) -> Result<MStats>;
   fn median(&self) -> Result<Med>;
   fn icorrelation(&self, other:&[i64]) -> Result<f64>;
   fn correlation(&self, other:&[f64]) -> Result<f64>;
   fn autocorr(&self) -> Result<f64>;
  
}

/// Implementing basic vector algebra.
pub trait Vectors {

   fn dotp(&self, other:&[f64]) -> f64;
   fn vsub(&self, other:&[f64]) -> Vec<f64>;
   fn vadd(&self, other:&[f64]) -> Vec<f64>;
   fn vmag(&self) -> f64;
   fn vdist(&self, other:&[f64]) -> f64;
   fn smult(&self, s:f64) -> Vec<f64>;
   fn vunit(&self) -> Vec<f64>;
   fn medoid(&self, other:usize) -> Result<(usize,f64)>;
   fn distances(&self, other: usize, other: &[f64] ) -> f64;   

}

/// Private helper function for formatting error messages
fn emsg(file:&'static str, line:u32, msg:&'static str)-> String {
   format!("{}:{} rstats {}",file,line,msg)
}

/// Private sum of linear weights 
fn wsum(n: usize) -> f64 { (n*(n+1)) as f64/2. }

/// Sorts a mutable `Vec<f64>` in place.  
/// It is the responsibility of the user to ensure that there are no NaNs etc.
pub fn sortf(v: &mut [f64]) { 
   v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
}

/// Generates a random f64 vector of `size` suitable for testing.  
/// Uses closure `rand` to generate random numbers to avoid dependencies.  
/// Random numbers are in the open interval 0..1 with uniform distribution.  
/// Requires two seeds.
pub fn genvec(size: usize , s1: usize, s2: usize) -> Vec<f64> {
   let mut m_z = s1 as u32;
   let mut m_w = s2 as u32;
   // returns f64 in the open interval 0..1 with uniform random distribution
   let mut rand = || {
      m_z = 36969 * (m_z & 65535) + (m_z >> 16);
      m_w = 18000 * (m_w & 65535) + (m_w >> 16);
      (((m_z << 16) & m_w) as f64 + 1.0)*2.328306435454494e-10
   };
   let mut v = Vec::with_capacity(size); 
   for _i in 0..size { v.push(rand()) }; // fills the lot with random numbers
   return v
}
