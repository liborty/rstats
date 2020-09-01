pub mod tests;
pub mod i64impls;
pub mod f64impls;

use anyhow::{Result,Context,ensure};
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
//  fn gwmeanstd(&self) -> Result<MStats>;
// fn median(&self) -> Result<Med>;
//   fn correlation(&self, &other) -> Result<f64>;
//  fn autocorr(&self) -> Result<f64>;

}

/// Private helper function for formatting error messages
fn emsg(file:&'static str, line:u32, msg:&'static str)-> String {
   format!("{}:{} rstats {}",file,line,msg)
}

/// Private sum of linear weights 
fn wsum(n: usize) -> f64 { (n*(n+1)) as f64/2. }




/// Linearly weighted version of gmeanstd.
/// # Example
/// ```
/// use rstats::gwmeanstd;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// let res = gwmeanstd(&VEC1).unwrap();
/// assert_eq!(res.mean,4.144953510241978_f64);
/// assert_eq!(res.std,2.1572089236412597_f64);
/// ```
pub fn gwmeanstd(dvec: &[i64]) -> Result<MStats> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} gwmeanstd - supplied sample is empty!",file!(),line!());
   let mut iw = n as i64; // descending weights
   let mut sum = 0f64;
   let mut sx2 = 0f64;
   for &x in dvec { 
      ensure!(x!=0i64,
         "{}:{} gwmeanstd does not accept zero valued data!",file!(),line!());  
      let lx = (x as f64).ln();
      sum += (iw as f64)*lx;
      sx2 += (iw as f64)*lx*lx;
      iw -= 1;
   }
   sum /= wsum(n);
   Ok( MStats { 
      mean : sum.exp(),
      std : (sx2 as f64/wsum(n) - sum.powi(2)).sqrt().exp() }
    )
}	

/// Fast median (avoids sorting).  
/// The data values must be within a moderate range not exceeding u16size (65535).
/// # Example
/// ```
/// use rstats::median;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// let res = median(&VEC1).unwrap();
/// assert_eq!(res.median,7.5_f64);
/// assert_eq!(res.lquartile,4_f64);
/// assert_eq!(res.uquartile,11_f64);
/// ```
pub fn median(data: &[i64]) -> Result<Med> {
   let max = *data.iter().max().with_context(||emsg(file!(),line!(),"median failed to find maximum"))?;
   let min = *data.iter().min().with_context(||emsg(file!(),line!(),"median failed to find minimum"))?;
   let range =  (max-min+1) as usize;
   ensure!(range <= u16::max_value() as usize, // range too big to use as subscripts
      "{}:{} median range {} of values exceeds u16",file!(),line!(),range);
	let mut acc = vec![0_u16; range]; // min max values inclusive
   for &item in data { acc[(item-min) as usize] += 1_u16 } // frequency distribution
   let mut result: Med = Default::default();
   let rowlength = data.len();
   let mut cumm = 0_usize;
   let mut i2;

   for i in 0..range {
      // find the lower quartile
      cumm += (acc[i]) as usize; // accummulate frequencies
      if 4 * cumm >= rowlength {
         result.lquartile = (i as i64 + min) as f64; // restore min value
         break;
      }
   }
   cumm = 0usize;
   for i in (0..range).rev() {
      // find the upper quartile
      cumm += (acc[i]) as usize; // accummulate frequencies
      if 4 * cumm >= rowlength {
         result.uquartile = (i as i64 + min) as f64;
         break;
      }
   }
   cumm = 0usize;
   for i in 0..range {
   // find the midpoint of the frequency distribution
      cumm += (acc[i]) as usize; // accummulate frequencies
      if 2 * cumm == rowlength {
         // even, the other half must have the same value
         i2 = i + 1;
         while acc[i2] == 0 { i2 += 1 }
         // first next non-zero acc[i2] must represent the other half
         result.median = ((i + i2) as i64 + 2*min) as f64 / 2.;
         break;
      }
      if 2 * cumm > rowlength {
         result.median = (i as i64 + min) as f64;
         break;
      }
      // first over the half, this must be the odd midpoint
   }
   Ok(result)
}

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
