use anyhow::{Result,ensure};
use crate::{MStats,RStats,wsum};

impl RStats for Vec<i64> {
   
   /// Arithmetic mean of an i64 slice
   /// # Example
   /// ```
   /// use rstats::RStats;
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

   /// Linearly weighted arithmetic mean of an i64 slice.     
   /// Linearly descending weights from n down to one.    
   /// Time dependent data should be in the stack order - the last being the oldest.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   /// assert_eq!(v1.awmean().unwrap(),5.333333333333333_f64);
   /// ```
   fn awmean(&self) -> Result<f64> {
   let n = self.len();
   ensure!(n>0,"{}:{} awmean - supplied sample is empty!",file!(),line!());
	let mut w = (n+1)as f64; // descending linear weights
	Ok( self.iter().map(|&x| { w -= 1.; w*x as f64 }).sum::<f64>() / wsum(n))
   }

   /// Liearly weighted arithmetic mean and standard deviation of an i64 slice.    
   /// Linearly descending weights from n down to one.    
   /// Time dependent data should be in the stack order - the last being the oldest.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   /// let res = v1.awmeanstd().unwrap();
   /// assert_eq!(res.mean,5.333333333333333_f64);
   /// assert_eq!(res.std,3.39934634239519_f64);
   /// ```
   fn awmeanstd(&self) -> Result<MStats> {
      let n = self.len();
      ensure!(n>0,"{}:{} awmeanstd - supplied sample is empty!",file!(),line!());
      let mut sx2 = 0f64;
      let mut w = n as f64; // descending linear weights
      let mean = self.iter().map( |&x| {
         let lx = x as f64;
         let wx = w*lx;
         sx2 += wx*lx;
         w -= 1.; 
         wx } ).sum::<f64>() as f64 / wsum(n);
   Ok( MStats { 
      mean : mean, 
      std : (sx2/wsum(n) - mean.powi(2)).sqrt() } )  
   }
   /// Harmonic mean of an i64 slice.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   /// assert_eq!(v1.hmean().unwrap(),4.305622526633627_f64);
   /// ```
   fn hmean(&self) -> Result<f64> {
   let n = self.len();
   ensure!(n>0,"{}:{} hmean - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   for &x in self {
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
/// use rstats::RStats;
/// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
/// assert_eq!(v1.hwmean().unwrap(),3.019546395306663_f64);
/// ```
fn hwmean(&self) -> Result<f64> {
   let mut n = self.len();
   ensure!(n>0,"{}:{} hwmean - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   for &x in self {
      ensure!(x!=0i64,
         "{}:{} hwmean does not accept zero valued data!",file!(),line!());  
      sum += n as f64/x as f64;
      n -= 1; 
   }
   Ok( wsum(self.len()) / sum )
}

}

impl RStats for Vec<f64> { 

   /// Arithmetic mean of an f64 slice
   /// # Example
   /// ```
   /// use rstats::RStats;
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
   /// Linearly weighted arithmetic mean of an f64 slice.     
   /// Linearly descending weights from n down to one.    
   /// Time dependent data should be in the stack order - the last being the oldest.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// assert_eq!(v1.awmean().unwrap(),5.333333333333333_f64);
   /// ```
   fn awmean(&self) -> Result<f64> {
      let n = self.len();
      ensure!(n>0,"{}:{} awmean - supplied sample is empty!",file!(),line!());
	   let mut iw = (n+1) as f64; // descending linear weights
	   Ok( self.iter().map(|&x| { iw -= 1.; iw*x }).sum::<f64>()/wsum(n))
   }

   /// Liearly weighted arithmetic mean and standard deviation of an f64 slice.    
   /// Linearly descending weights from n down to one.    
   /// Time dependent data should be in the stack order - the last being the oldest.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// let res = v1.awmeanstd().unwrap();
   /// assert_eq!(res.mean,5.333333333333333_f64);
   /// assert_eq!(res.std,3.39934634239519_f64);
   /// ```
   fn awmeanstd(&self) -> Result<MStats> {
      let n = self.len();
      ensure!(n>0,"{}:{} awmeanstd - supplied sample is empty!",file!(),line!());
      let mut sx2 = 0f64;
      let mut w = n as f64; // descending linear weights
      let mean = self.iter().map( |&x| {
         let wx = w*x;
         sx2 += wx*x;
         w -= 1.; 
         wx } ).sum::<f64>() / wsum(n);
   Ok( MStats { 
      mean : mean, 
      std : (sx2/wsum(n) - mean.powi(2)).sqrt() } )  
   }
   /// Harmonic mean of an f64 slice.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// assert_eq!(v1.hmean().unwrap(),4.305622526633627_f64);
   /// ```
   fn hmean(&self) -> Result<f64> {
      let n = self.len();
      ensure!(n>0,"{}:{} hmean - supplied sample is empty!",file!(),line!());
      let mut sum = 0f64;
      for &x in self {
         ensure!(x != 0f64,"{}:{} hmean does not accept zero valued data!",file!(),line!());  
         sum += 1.0/x
      }
      Ok ( n as f64 / sum )
   }
   /// Linearly weighted harmonic mean of an f64 slice.    
   /// Linearly descending weights from n down to one.    
   /// Time dependent data should be in the stack order - the last being the oldest.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
/// assert_eq!(v1.hwmean().unwrap(),3.019546395306663_f64);
/// ```
   fn hwmean(&self) -> Result<f64> {
      let mut n = self.len();
      ensure!(n>0,"{}:{} hwmean - supplied sample is empty!",file!(),line!());
      let mut sum = 0f64;
      for &x in self {
         ensure!(x!=0f64,
         "{}:{} hwmean does not accept zero valued data!",file!(),line!());  
         sum += n as f64/x;
         n -= 1; 
      }
   Ok( wsum(self.len()) / sum )
   }
}
