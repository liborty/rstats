pub mod tests;
pub mod i64impls;
pub mod f64impls;

use anyhow::{Result,ensure};
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
//   fn correlation(&self, &other) -> Result<f64>;
//  fn autocorr(&self) -> Result<f64>;

}

/// Private helper function for formatting error messages
fn emsg(file:&'static str, line:u32, msg:&'static str)-> String {
   format!("{}:{} rstats {}",file,line,msg)
}

/// Private sum of linear weights 
fn wsum(n: usize) -> f64 { (n*(n+1)) as f64/2. }


/// Correlation coefficient of a sample of two integer variables.
/// # Example
/// ```
/// use rstats::correlation;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// const VEC2:[i64;14] = [14,13,12,11,10,9,8,7,6,5,4,3,2,1];
/// assert_eq!(correlation(&VEC1,&VEC2).unwrap(),-1_f64);
/// ```
pub fn correlation(v1:&[i64],v2:&[i64]) -> Result<f64> {
   let n = v1.len();
   ensure!(n>0,emsg(file!(),line!(),"correlation - first sample is empty"));
   ensure!(n==v2.len(),emsg(file!(),line!(),"correlation - samples are not of the same size"));
   let (mut sy,mut sxy,mut sx2,mut sy2) = (0,0,0,0);
   let sx:i64 = v1.iter().enumerate().map(|(i,&x)| {
      let y = v2[i]; sy += y; sxy += x*y; sx2 += x*x; sy2 += y*y; x    
      }).sum();
   let (sxf,syf,sxyf,sx2f,sy2f,nf) = 
       (sx as f64,sy as f64,sxy as f64,sx2 as f64,sy2 as f64,n as f64);
   Ok( (sxyf-sxf/nf*syf)/(((sx2f-sxf/nf*sxf)*(sy2f-syf/nf*syf)).sqrt()) )
}

/// (Auto)correlation coefficient of pairs of successive values of (time series) integer variable.
/// # Example
/// ```
/// use rstats::autocorr;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// assert_eq!(autocorr(&VEC1).unwrap(),0.9984603532054123_f64);
/// ```
pub fn autocorr(v1:&[i64]) -> Result<f64> {
   let n = v1.len();
   ensure!(n>=2,emsg(file!(),line!(),"autocorr - sample is too small"));
   let (mut sx,mut sy,mut sxy,mut sx2,mut sy2) = (0,0,0,0,0);
   for i in 0..n-1 {
      let x = v1[i]; let y = v1[i+1]; 
      sx += x; sy += y; sxy += x*y; sx2 += x*x; sy2 += y*y }    
   let (sxf,syf,sxyf,sx2f,sy2f,nf) = 
       (sx as f64,sy as f64,sxy as f64,sx2 as f64,sy2 as f64,n as f64);
   Ok( (sxyf-sxf/nf*syf)/(((sx2f-sxf/nf*sxf)*(sy2f-syf/nf*syf)).sqrt()) )
}
