use anyhow::{Result,ensure,bail};
use crate::{RStats,MStats,Med,wsum,emsg};

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
      ensure!(n>0,emsg(file!(),line!(),"amean - sample is empty!"));
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
      ensure!(n>0,emsg(file!(),line!(),"ameanstd - sample is empty!"));
      let mut sx2 = 0_f64;
      let mean = self.iter().map(|&x|{ sx2+=x*x; x}).sum::<f64>() / (n as f64);
      Ok( MStats { mean : mean, std : (sx2 /(n as f64) - mean.powi(2)).sqrt() } )
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
      ensure!(n>0,emsg(file!(),line!(),"awmean - sample is empty!"));
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
      ensure!(n>0,emsg(file!(),line!(),"awmeanstd - sample is empty!"));
      let mut sx2 = 0_f64;
      let mut w = n as f64; // descending linear weights
      let mean = self.iter().map( |&x| {
         let wx = w*x;
         sx2 += wx*x; w -= 1_f64; wx } ).sum::<f64>() / wsum(n);
   Ok( MStats { mean : mean, std : (sx2/wsum(n) - mean.powi(2)).sqrt() } )  
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
      ensure!(n>0,emsg(file!(),line!(),"hmean - sample is empty!"));
      let mut sum = 0_f64;
      for &x in self {
         ensure!(x.is_normal(),emsg(file!(),line!(),"hmean does not accept zero valued data!"));  
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
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"hwmean - sample is empty!"));
      let mut sum = 0_f64; let mut w = n as f64;
      for &x in self {
         ensure!(x.is_normal(),emsg(file!(),line!(),"hwmean does not accept zero valued data!")); 
         sum += w/x; w -= 1_f64; 
      }
   Ok( wsum(n)/sum )
   }

   /// Geometric mean of an i64 slice.  
   /// The geometric mean is just an exponential of an arithmetic mean
   /// of log data (natural logarithms of the data items).  
   /// The geometric mean is less sensitive to outliers near maximal value.  
   /// Zero valued data is not allowed.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// assert_eq!(v1.gmean().unwrap(),6.045855171418503_f64);
   /// ```
   fn gmean(&self) -> Result<f64> {
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"gmean - sample is empty!"));
      let mut sum = 0_f64;
      for &x in self {   
         ensure!(x.is_normal(),emsg(file!(),line!(),"gmean does not accept zero valued data!"));
         sum += x.ln()
      }
   Ok( (sum/(n as f64)).exp() )
   }

   /// Linearly weighted geometric mean of an i64 slice.  
   /// Descending weights from n down to one.    
   /// Time dependent data should be in the stack order - the last being the oldest.  
   /// The geometric mean is just an exponential of an arithmetic mean
   /// of log data (natural logarithms of the data items).  
   /// The geometric mean is less sensitive to outliers near maximal value.  
   /// Zero data is not allowed - would at best only produce zero result.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// assert_eq!(v1.gwmean().unwrap(),4.144953510241978_f64);
   /// ```
   fn gwmean(&self) -> Result<f64> {
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"gwmean - sample is empty!"));
      let mut w = n as f64; // descending weights
      let mut sum = 0_f64;
      for &x in self {  
         ensure!(x.is_normal(),emsg(file!(),line!(),"gwmean does not accept zero valued data!"));
         sum += w*x.ln();
         w -= 1_f64;
      }
   Ok( (sum/wsum(n)).exp() )
   }

   /// Geometric mean and std ratio of an f64 slice.  
   /// Zero valued data is not allowed.  
   /// Std of ln data becomes a ratio after conversion back.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// let res = v1.gmeanstd().unwrap();
   /// assert_eq!(res.mean,6.045855171418503_f64);
   /// assert_eq!(res.std,2.1084348239406303_f64);
   /// ```
   fn gmeanstd(&self) -> Result<MStats> {
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"gmeanstd - sample is empty!"));
      let mut sum = 0_f64; let mut sx2 = 0_f64;
      for &x in self { 
         ensure!(x.is_normal(),emsg(file!(),line!(),"gmeanstd does not accept zero valued data!"));
         let lx = x.ln();
         sum += lx; sx2 += lx*lx    
      }
      sum /= n as f64;
      Ok( MStats { mean: sum.exp(), std: (sx2/(n as f64) - sum.powi(2)).sqrt().exp() } )
   }

   /// Linearly weighted version of gmeanstd.
   /// # Example
   /// ```
   /// use rstats::RStats;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// let res = v1.gwmeanstd().unwrap();
   /// assert_eq!(res.mean,4.144953510241978_f64);
   /// assert_eq!(res.std,2.1572089236412597_f64);
   /// ```
   fn gwmeanstd(&self) -> Result<MStats> {
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"gwmeanstd - sample is empty!"));
      let mut w = n as f64; // descending weights
      let mut sum = 0_f64; let mut sx2 = 0_f64;
      for &x in self { 
         ensure!(x.is_normal(),emsg(file!(),line!(),"gwmeanstd does not accept zero valued data!"));
         let lnx = x.ln();
         sum += w*lnx; sx2 += w*lnx*lnx;
         w -= 1_f64;
      }
   sum /= wsum(n);
   Ok( MStats { mean : sum.exp(), std : (sx2 as f64/wsum(n) - sum.powi(2)).sqrt().exp() } )
   }

   fn median(&self) -> Result<Med> {
      bail!(emsg(file!(),line!(),"gwmeanstd - sample is empty!"));
   }

}
