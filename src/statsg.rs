use crate::{ here,functions::wsum, MStats, Med, Stats };
use anyhow::{ensure, Result};
// use std::ops::Sub;

pub use indxvec::merge::{sortm,minmax};          

impl<T> Stats for &[T] 
    where T: Copy+PartialOrd, // +Sub::<Output = T>,
        f64: From<T> {  

    /// Vector magnitude
    fn vmag(self) -> f64 {
        self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>().sqrt()
    }

    /// Vector magnitude squared (sum of squares)
    fn vmagsq(self) -> f64  {
        self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>()
    }

    /// Vector with inverse magnitude
    fn vinverse(self) -> Vec<f64> {
        let sf = 1.0/self.vmagsq();
        self.iter().map(|&x| sf*(f64::from(x))).collect()     
    }

    // negated vector (all components swap sign)
    fn negv(self) -> Vec<f64> { 
        self.iter().map(|&x| (-f64::from(x))).collect()
    }

    /// Unit vector
    fn vunit(self) -> Vec<f64> {
        let m = 1.0 / self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>().sqrt();
        self.iter().map(|&x| m*(f64::from(x))).collect() 
    }

    /// Arithmetic mean of an f64 slice
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().amean().unwrap(),7.5_f64);
    /// ```
    fn amean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!",here!());
        Ok(self.iter().map(|&x| f64::from(x)).sum::<f64>() / (n as f64))
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
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!",here!());
        let mut sx2 = 0_f64;
        let mean = self
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

    /// Linearly weighted arithmetic mean and standard deviation of an f64 slice.    
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
        for &x in self {
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
        for &x in self {
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
        for &x in self {
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
        for &x in self {
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
        for &x in self {
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
        for &x in self {
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
    /// use rstats::{Stats};
    /// let v1 = vec![1_u8,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// let res = &v1.median().unwrap();
    /// assert_eq!(res.median,7.5_f64);
    /// assert_eq!(res.lquartile,4.25_f64);
    /// assert_eq!(res.uquartile,10.75_f64);
    /// ```
    fn median(self) -> Result<Med> {
        let gaps = self.len()-1;
        let mid = gaps / 2;
        let quarter = gaps / 4;
        let threeq = 3 * gaps / 4;
        let v = sortm(self,true);     
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

    /// Probability density function of a sorted slice with repeats
    fn pdf(self) -> Vec<f64> {     
        let n = self.len();   
        let mut res:Vec<f64> = vec![1.0;n];  
        let mut rcount = 1_f64; // running count
        let mut lastval = self[0]; // first item

        for i in 1..n { // first pass to count same values 
            if self[i] > lastval { lastval = self[i]; rcount = 1.0 } // new value                 
            else { rcount += 1.0 } // same value
            res[i] = rcount            
        }
        let nf = n as f64;
        for i in (0..n).rev() {  // second reverse pass to propagate maxima
            if self[i] < lastval { lastval = self[i]; rcount = res[i] }; // new value                
            res[i] = rcount/nf // propagate maximum count from previous item           
        } 
        // let ressum = res.iter().sum::<f64>();
        // eprintln!("Sum of probs: {}",ressum); 
        res
    } 

    /// Information (entropy) (in nats)
    fn entropy(self) -> f64 {
        let pdfv = sortm(self,true).pdf();      
        pdfv.iter().map( |&x| -x*(x.ln()) ).sum()                 
    }

    /// (Auto)correlation coefficient of pairs of successive values of (time series) f64 variable.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.autocorr(),0.9984603532054123_f64);
    /// ```
    fn autocorr(self) -> f64 {
        let (mut sx, mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64, 0_f64);
        let n = self.len();
        if n < 2 { panic!("{} vector is too short",here!()) }
        let mut x = f64::from(self[0]);    
        for i in 1..n {
            let y = f64::from(self[i]);
            sx += x;
            sy += y;
            sxy += x * y;
            sx2 += x * x;
            sy2 += y * y;
            x = y
        }        
        let nf = n as f64;
        (sxy - sx / nf * sy) / ((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()
    }
    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Vec<f64> {
        let (min,_,max,_) = minmax(self);
        let range = f64::from(max)-f64::from(min);
        self.iter().map(|&x|(f64::from(x)-f64::from(min))/range).collect()        
    }
    /// Reconstructs the full symmetric square matrix from its lower diagonal compact form,
    /// as produced by covar, covone, wcovar
    fn symmatrix(self) -> Vec<Vec<f64>> {
        // solve quadratic equation to find the dimension of the square matrix
        let n = (((8*self.len()+1) as f64).sqrt() as usize - 1)/2;
        let mut mat = vec![vec![0_f64;n];n]; // create the square matrix 
        let mut selfindex = 0;
        for row in 0..n {     
            for column in 0..row { // excludes the diagonal 
                let this = f64::from(self[selfindex]); 
                mat[row][column] = this; // just copy the value into the lower triangle
                mat[column][row] = this; // and into transposed upper position 
                selfindex += 1 // move to the next input value
            } // this row of lower triangle finished
            mat[row][row] = f64::from(self[selfindex]);  // now the diagonal element, no transpose
            selfindex += 1 // move to the next input value
        }
        mat
    }    
 
}
