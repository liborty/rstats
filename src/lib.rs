// pub mod tests;
pub mod i64impls;
pub mod f64impls;
pub mod vimpls;
pub mod tests;


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

pub trait Vectors {

   fn sortf(&mut self);
   fn dotp(&self, other:&[f64]) -> f64;
   fn vsub(&self, other:&[f64]) -> Vec<f64>;
   fn vmag(&self) -> f64;
}

/// Private helper function for formatting error messages
fn emsg(file:&'static str, line:u32, msg:&'static str)-> String {
   format!("{}:{} rstats {}",file,line,msg)
}

/// Private sum of linear weights 
fn wsum(n: usize) -> f64 { (n*(n+1)) as f64/2. }
