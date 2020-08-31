mod tests;
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
}

impl RStats for Vec<i64> {
   
   /// Arithmetic mean of an i64 slice
   /// # Example
   /// ```
   /// use crate::rstats::RStats;
   /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   /// assert_eq!(v1.amean().unwrap(),7.5_f64);
   /// ```
   fn amean(&self) -> Result<f64> { 
      let n = self.len();
      ensure!(n > 0, "{}:{} amean - supplied sample is empty!",file!(),line!() );
      Ok( self.iter().map(|&x| x as f64).sum::<f64>() / (n as f64) )
   } 

   /// Arithmetic mean and standard deviation of an i64 slice
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   /// let res = v1.ameanstd().unwrap();
   /// assert_eq!(res.mean,7.5_f64);
   /// assert_eq!(res.std,4.031128874149275_f64);
   /// ```
   fn ameanstd(&self) -> Result<MStats> {
   let n = self.len();
   ensure!(n > 0,"{}:{} ameanstd - supplied sample is empty!",file!(),line!());
   let mut sx2 = 0_f64;
   let mean = self.iter().map(|&x|{ let lx = x as f64;sx2+=lx*lx; lx}).sum::<f64>() / (n as f64);
   Ok( MStats { 
      mean : mean, 
      std : (sx2 /(n as f64) - mean.powi(2)).sqrt() } )
}

}

impl RStats for Vec<f64> { 

   /// Arithmetic mean of an f64 slice
   /// # Example
   /// ```
   /// use crate::rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// assert_eq!(v1.amean().unwrap(),7.5_f64);
   /// ```
   fn amean(&self) -> Result<f64> { 
      let n = self.len();
      ensure!(n > 0, "{}:{} amean - supplied sample is empty!",file!(),line!() );
      Ok( self.iter().sum::<f64>() / (n as f64) )
   }  
   
     /// Arithmetic mean and standard deviation of an f64 slice
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// let res = v1.ameanstd().unwrap();
   /// assert_eq!(res.mean,7.5_f64);
   /// assert_eq!(res.std,4.031128874149275_f64);
   /// ```
   fn ameanstd(&self) -> Result<MStats> {
      let n = self.len();
      ensure!(n > 0,"{}:{} ameanstd - supplied sample is empty!",file!(),line!());
      let mut sx2 = 0_f64;
      let mean = self.iter().map(|&x|{ sx2+=x*x; x}).sum::<f64>() / (n as f64);
      Ok( MStats { 
         mean : mean, 
         std : (sx2 /(n as f64) - mean.powi(2)).sqrt() } )
   }
}

/// Private helper function for formatting error messages
fn cmsg(file:&'static str, line:u32, msg:&'static str)-> String {
   format!("{}:{} stats {}",file,line,msg)
}

/// Private sum of linear weights 
fn wsum(n: usize) -> f64 { (n*(n+1)) as f64/2. }

/// Arithmetic mean of an i64 slice
/// # Example
/// ```
/// use rstats::amean;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// assert_eq!(amean(&VEC1).unwrap(),7.5_f64);
/// ```
pub fn amean(dvec: &[i64]) -> Result<f64> { 
   let n = dvec.len();
   ensure!(n > 0, "{}:{} amean - supplied sample is empty!",file!(),line!() );
   Ok( dvec.iter().map(|&x| x as f64).sum::<f64>() / (n as f64) )
}

/// Arithmetic mean and standard deviation of an i64 slice
/// # Example
/// ```
/// use rstats::ameanstd;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// let res = ameanstd(&VEC1).unwrap();
/// assert_eq!(res.mean,7.5_f64);
/// assert_eq!(res.std,4.031128874149275_f64);
/// ```
pub fn ameanstd(dvec: &[i64]) -> Result<MStats> {
   let n = dvec.len();
   ensure!(n > 0,"{}:{} ameanstd - supplied sample is empty!",file!(),line!());
   let mut sx2:i64 = 0;
   let mean = dvec.iter().map(|&x|{ sx2+=x*x; x}).sum::<i64>() as f64 / (n as f64);
   Ok( MStats { 
      mean : mean, 
      std : (sx2 as f64/(n as f64) - mean.powi(2)).sqrt() } )
}

/// Linearly weighted arithmetic mean of an i64 slice.     
/// Linearly descending weights from n down to one.    
/// Time dependent data should be in the stack order - the last being the oldest.
/// # Example
/// ```
/// use rstats::awmean;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// assert_eq!(awmean(&VEC1).unwrap(),5.333333333333333_f64);
/// ```
pub fn awmean(dvec: &[i64]) -> Result<f64> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} awmean - supplied sample is empty!",file!(),line!());
	let mut iw = dvec.len() as i64 + 1; // descending linear weights
	Ok( dvec.iter().map(|&x| { iw -= 1; iw*x }).sum::<i64>() as f64 / wsum(n))
}

/// Liearly weighted arithmetic mean and standard deviation of an i64 slice.    
/// Linearly descending weights from n down to one.    
/// Time dependent data should be in the stack order - the last being the oldest.
/// # Example
/// ```
/// use rstats::awmeanstd;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// let res = awmeanstd(&VEC1).unwrap();
/// assert_eq!(res.mean,5.333333333333333_f64);
/// assert_eq!(res.std,3.39934634239519_f64);
/// ```
pub fn awmeanstd(dvec: &[i64]) -> Result<MStats> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} awmeanstd - supplied sample is empty!",file!(),line!());
   let mut sx2 = 0f64;
   let mut iw = n as f64; // descending linear weights
   let mean = dvec.iter().map( |&x| { 
      let wx = iw*x as f64;
      sx2 += wx*x as f64;
      iw -= 1.; 
      wx } ).sum::<f64>() as f64 / wsum(n);
   Ok( MStats { 
      mean : mean, 
      std : (sx2 as f64/wsum(n) - mean.powi(2)).sqrt() } )  
}

/// Harmonic mean of an i64 slice.
/// # Example
/// ```
/// use rstats::hmean;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// assert_eq!(hmean(&VEC1).unwrap(),4.305622526633627_f64);
/// ```
pub fn hmean(dvec: &[i64]) -> Result<f64> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} hmean - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   for &x in dvec {
      ensure!(x != 0i64,"{}:{} hmean does not accept zero valued data!",file!(),line!());  
      sum += 1.0/(x as f64) 
   }
   Ok ( n as f64 / sum )
}

/// Linearly weighted harmonic mean of an i64 slice.    
/// Linearly descending weights from n down to one.    
/// Time dependent data should be in the stack order - the last being the oldest.
/// # Example
/// ```
/// use rstats::hwmean;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// assert_eq!(hwmean(&VEC1).unwrap(),3.019546395306663_f64);
/// ```
pub fn hwmean(dvec: &[i64]) -> Result<f64> {
   let mut n = dvec.len();
   ensure!(n>0,"{}:{} hwmean - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   for &x in dvec {
      ensure!(x!=0i64,
         "{}:{} hwmean does not accept zero valued data!",file!(),line!());  
      sum += n as f64/x as f64;
      n -= 1; 
   }
   Ok( wsum(dvec.len()) / sum )
}

/// Geometric mean of an i64 slice.  
/// The geometric mean is just an exponential of an arithmetic mean
/// of log data (natural logarithms of the data items).  
/// The geometric mean is less sensitive to outliers near maximal value.  
/// Zero valued data is not allowed.
/// # Example
/// ```
/// use rstats::gmean;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// assert_eq!(gmean(&VEC1).unwrap(),6.045855171418503_f64);
/// ```
pub fn gmean(dvec: &[i64]) -> Result<f64> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} gmean - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   for &x in dvec {   
      ensure!(x!=0i64,
         "{}:{} gmean does not accept zero valued data!",file!(),line!()); 
      sum += (x as f64).ln()
   }
   Ok( (sum/(n as f64)).exp() )
}

/// Time linearly weighted geometric mean of an i64 slice.  
/// Linearly descending weights from n down to one.  
/// Time dependent data should be in the stack order - the last being the oldest.  
/// The geometric mean is just an exponential of an arithmetic mean
/// of log data (natural logarithms of the data items).  
/// The geometric mean is less sensitive to outliers near maximal value.  
/// Zero data is not allowed - would at best only produce zero result.
/// # Example
/// ```
/// use rstats::gwmean;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// assert_eq!(gwmean(&VEC1).unwrap(),4.144953510241978_f64);
/// ```
pub fn gwmean(dvec: &[i64]) -> Result<f64> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} gwmean - supplied sample is empty!",file!(),line!());
   let mut iw = n as i64; // descending weights
   let mut sum = 0f64;
   for &x in dvec {  
      ensure!(x!=0i64,
         "{}:{} gwmean does not accept zero valued data!",file!(),line!()); 
      sum += (iw as f64)*(x as f64).ln();
      iw -= 1;
   }
   Ok( (sum/wsum(n)).exp() )
}	
/// Geometric mean and std ratio of an i64 slice.  
/// Zero valued data is not allowed.  
/// Std of ln data becomes a ratio after conversion back.
/// # Example
/// ```
/// use rstats::gmeanstd;
/// const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// let res = gmeanstd(&VEC1).unwrap();
/// assert_eq!(res.mean,6.045855171418503_f64);
/// assert_eq!(res.std,2.1084348239406303_f64);
/// ```
pub fn gmeanstd(dvec: &[i64]) -> Result<MStats> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} gmeanstd - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   let mut sx2 = 0f64;
   for &x in dvec { 
      ensure!(x!=0i64,
         "{}:{} gmeanstd does not accept zero valued data!",file!(),line!());   
      let lx = (x as f64).ln();
      sum += lx;
      sx2 += lx*lx    
   }
   sum /= n as f64;
   Ok( MStats { 
      mean: sum.exp(), 
      std: (sx2/(n as f64) - sum.powi(2)).sqrt().exp() }
    )
}

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
   let max = *data.iter().max().with_context(||cmsg(file!(),line!(),"median failed to find maximum"))?;
   let min = *data.iter().min().with_context(||cmsg(file!(),line!(),"median failed to find minimum"))?;
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
   ensure!(n>0,cmsg(file!(),line!(),"correlation - first sample is empty"));
   ensure!(n==v2.len(),cmsg(file!(),line!(),"correlation - samples are not of the same size"));
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
   ensure!(n>=2,cmsg(file!(),line!(),"autocorr - sample is too small"));
   let (mut sx,mut sy,mut sxy,mut sx2,mut sy2) = (0,0,0,0,0);
   for i in 0..n-1 {
      let x = v1[i]; let y = v1[i+1]; 
      sx += x; sy += y; sxy += x*y; sx2 += x*x; sy2 += y*y }    
   let (sxf,syf,sxyf,sx2f,sy2f,nf) = 
       (sx as f64,sy as f64,sxy as f64,sx2 as f64,sy2 as f64,n as f64);
   Ok( (sxyf-sxf/nf*syf)/(((sx2f-sxf/nf*sxf)*(sy2f-syf/nf*syf)).sqrt()) )
}
