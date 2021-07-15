use crate::{MStats, Med, Stats, functions::wsum, here};
use anyhow::{ensure, Result};
use std::ops::*;

pub use indxvec::merge::sortm;

pub struct GSlice<'a,T> ( &'a[T] ) where T:Copy;

impl<'a,T> Deref for GSlice<'a,T> where T:Copy {
    type Target = &'a[T];
    #[inline] 
    fn deref(&self) -> &Self::Target { &self.0 }
}
impl<'a,T> DerefMut for GSlice<'a,T> where T:Copy{ 
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<'a,T> Stats for GSlice<'a,T> 
    where T: Copy+PartialOrd
        +Add::<Output = T>
        +Mul::<Output = T>,
        f64: From<T> {  
    /// Arithmetic mean of an f64 slice
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().amean().unwrap(),7.5_f64);
    /// ```
    fn amean(self) -> Result<f64> {
        let n = self.0.len();
        ensure!(n > 0, "{} sample is empty!",here!());
        Ok(self.0.iter().map(|&x| f64::from(x)).sum::<f64>() / (n as f64))
    }

    /// Arithmetic mean and (population) standard deviation of an f64 slice
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.as_slice().ameanstd().unwrap();
    /// assert_eq!(res.mean,7.5_f64);
    /// assert_eq!(res.std,4.031128874149275_f64);
    /// ```
    fn ameanstd(self) -> Result<MStats> {
        let n = self.0.len();
        ensure!(n > 0, "{} sample is empty!",here!());
        let mut sx2 = 0_f64;
        let mean = self.0
            .iter()
            .map(|&x| {
                sx2 += f64::from(x) * f64::from(x);
                f64::from(x)
            })
            .sum::<f64>() / (n as f64);
        Ok(MStats {
            mean,
            std: (sx2 / (n as f64) - mean.powi(2)).sqrt(),
        })
    }

    /// Linearly weighted arithmetic mean of an f64 slice.     
    /// Linearly descending weights from n down to one.    
    /// Time dependent data should be in the stack order - the last being the oldest.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().awmean().unwrap(),5.333333333333333_f64);
    /// ```
    fn awmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut iw = (n + 1) as f64; // descending linear weights
        Ok(self
            .iter()
            .map(|&x| {
                iw -= 1.;
                iw * f64::from(x)
            })
            .sum::<f64>()
            / wsum(n))
    }

    /// Liearly weighted arithmetic mean and standard deviation of an f64 slice.    
    /// Linearly descending weights from n down to one.    
    /// Time dependent data should be in the stack order - the last being the oldest.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.as_slice().awmeanstd().unwrap();
    /// assert_eq!(res.mean,5.333333333333333_f64);
    /// assert_eq!(res.std,3.39934634239519_f64);
    /// ```
    fn awmeanstd(self) -> Result<MStats> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());        
        let mut sx2 = 0_f64;
        let mut w = n as f64; // descending linear weights
        let mean = self
            .iter()
            .map(|&x| {
                let wx = w * f64::from(x);
                sx2 += wx * f64::from(x);
                w -= 1_f64;
                wx
            })
            .sum::<f64>()
            / wsum(n);
        Ok(MStats {
            mean: mean,
            std: (sx2 / wsum(n) - mean.powi(2)).sqrt(),
        })
    }

    /// Harmonic mean of an f64 slice.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().hmean().unwrap(),4.305622526633627_f64);
    /// ```
    fn hmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());        
        let mut sum = 0_f64;
        for &x in *self {
            let fx = f64::from(x);
            ensure!( fx.is_normal(),"{} does not accept zero valued data!",here!());         
            sum += 1.0 / fx
        }
        Ok(n as f64 / sum)
    }
    /// Linearly weighted harmonic mean of an f64 slice.    
    /// Linearly descending weights from n down to one.    
    /// Time dependent data should be in the stack order - the last being the oldest.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().hwmean().unwrap(),3.019546395306663_f64);
    /// ```
    fn hwmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut sum = 0_f64;
        let mut w = n as f64;
        for &x in *self {
            let fx = f64::from(x);
            ensure!(fx.is_normal(),"{} does not accept zero valued data!",here!());
            sum += w / fx;
            w -= 1_f64;
        }
        Ok(wsum(n) / sum)
    }

    /// Geometric mean of an i64 slice.  
    /// The geometric mean is just an exponential of an arithmetic mean
    /// of log data (natural logarithms of the data items).  
    /// The geometric mean is less sensitive to outliers near maximal value.  
    /// Zero valued data is not allowed.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().gmean().unwrap(),6.045855171418503_f64);
    /// ```
    fn gmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut sum = 0_f64;
        for &x in *self {
            let fx = f64::from(x);
            ensure!( fx.is_normal(),"{} does not accept zero valued data!",here!());        
            sum += fx.ln()
        }
        Ok((sum / (n as f64)).exp())
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
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().gwmean().unwrap(),4.144953510241978_f64);
    /// ```
    fn gwmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut w = n as f64; // descending weights
        let mut sum = 0_f64;
        for &x in *self {
            let fx = f64::from(x);
            ensure!(fx.is_normal(),"{} does not accept zero valued data!",here!());
            sum += w * fx.ln();
            w -= 1_f64;
        }
        Ok((sum / wsum(n)).exp())
    }

    /// Geometric mean and std ratio of an f64 slice.  
    /// Zero valued data is not allowed.  
    /// Std of ln data becomes a ratio after conversion back.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.as_slice().gmeanstd().unwrap();
    /// assert_eq!(res.mean,6.045855171418503_f64);
    /// assert_eq!(res.std,2.1084348239406303_f64);
    /// ```
    fn gmeanstd(self) -> Result<MStats> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut sum = 0_f64;
        let mut sx2 = 0_f64;
        for &x in *self {
            let fx = f64::from(x);
            ensure!(fx.is_normal(),"{} does not accept zero valued data!",here!());
            let lx = fx.ln();
            sum += lx;
            sx2 += lx * lx
        }
        sum /= n as f64;
        Ok(MStats {
            mean: sum.exp(),
            std: (sx2 / (n as f64) - sum.powi(2)).sqrt().exp(),
        })
    }

    /// Linearly weighted version of gmeanstd.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.as_slice().gwmeanstd().unwrap();
    /// assert_eq!(res.mean,4.144953510241978_f64);
    /// assert_eq!(res.std,2.1572089236412597_f64);
    /// ```
    fn gwmeanstd(self) -> Result<MStats> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!()); 
        let mut w = n as f64; // descending weights
        let mut sum = 0_f64;
        let mut sx2 = 0_f64;
        for &x in *self {
            let fx = f64::from(x);
            ensure!(fx.is_normal(),"{} does not accept zero valued data!",here!());
            let lnx = fx.ln();
            sum += w * lnx;
            sx2 += w * lnx * lnx;
            w -= 1_f64;
        }
        sum /= wsum(n);
        Ok(MStats {
            mean: sum.exp(),
            std: (sx2 as f64 / wsum(n) - sum.powi(2)).sqrt().exp(),
        })
    }

    /// Median of a &[T] slice
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_u8,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// let res = v1.as_slice().median().unwrap();
    /// assert_eq!(res.median,7.5_f64);
    /// assert_eq!(res.lquartile,4.25_f64);
    /// assert_eq!(res.uquartile,10.75_f64);
    /// ```
    fn median(self) -> Result<Med> {
        let gaps = self.0.len()-1;
        let mid = gaps / 2;
        let quarter = gaps / 4;
        let threeq = 3 * gaps / 4;
        let v = sortm(self.0,true);     
        let mut result: Med = Default::default();
        result.median = if 2*mid < gaps { (f64::from(v[mid]) + f64::from(v[mid + 1])) / 2.0 }
            else { f64::from(v[mid]) };
        match gaps % 4 {
        0 => {
            result.lquartile = f64::from(v[quarter]);
            result.uquartile = f64::from(v[threeq]);
            return Ok(result) },
        1 => {
            result.lquartile = (3.*f64::from(v[quarter]) + f64::from(v[quarter+1])) / 4_f64;
            result.uquartile = (f64::from(v[threeq]) + 3.*f64::from(v[threeq+1])) / 4_f64;
            return Ok(result) },
        2 => {
            result.lquartile = (f64::from(v[quarter]) + f64::from(v[quarter+1])) / 2.;
            result.uquartile = (f64::from(v[threeq]) + f64::from(v[threeq+1])) / 2.;
            return Ok(result) },
        3 => {
            result.lquartile = (f64::from(v[quarter]) + 3.*f64::from(v[quarter+1])) / 4.;
            result.uquartile = (3.*f64::from(v[threeq]) + f64::from(v[threeq+1])) / 4.
            },
        _ => { }  
        }
        Ok(result)       
    }   
}
