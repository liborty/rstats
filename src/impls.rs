use anyhow::{Result,ensure};
use crate::{Med,MStats,RStats,wsum};

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
   /// Linearly weighted arithmetic mean of an i64 slice.     
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
}
