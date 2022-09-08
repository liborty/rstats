use crate::{error::RError,RE,sumn,Vecg,MutVecg,MStats,Stats};
use indxvec::Vecops;
use medians::Median; 

impl<T> Stats for &[T] 
    where T: Copy+PartialOrd+std::fmt::Display,f64:From<T> {  

    /// Vector magnitude
    fn vmag(self) -> f64 {
        match self.len() {
            0 => 0_f64,
            1 => f64::from(self[0]),
            _ =>  self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>().sqrt()
        }
    }

    /// Vector magnitude squared (sum of squares)
    fn vmagsq(self) -> f64 {
        match self.len() {
            0 => 0_f64,
            1 => f64::from(self[0]).powi(2),
            _=> self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>()
        }
    }

    /// Vector with reciprocal components
    fn vreciprocal(self) -> Result<Vec<f64>,RE> {
        if self.is_empty() { return  Err(RError::NoDataError("empty self vec"))  };
        for &component in self {
           let c = f64::from(component); 
           if !c.is_normal() { return Err(RError::ArithError("no reciprocal with zero component")); }; 
        }     
        Ok( self.iter().map(|&x| 1.0/(f64::from(x))).collect() )     
    }

    /// Vector with inverse magnitude
    fn vinverse(self) -> Result<Vec<f64>,RE> {
        if self.is_empty() { return  Err(RError::NoDataError("empty self vec"))  };
        let mag = self.vmagsq();
        if mag > 0.0 {  
            Ok( self.iter().map(|&x| f64::from(x)/mag).collect() ) }
        else { Err(RError::DataError("no inverse of zero vector magnitude")) }    
    }

    // negated vector (all components swap sign)
    fn negv(self) -> Vec<f64> { 
        if self.is_empty() { return Vec::new();  };
        self.iter().map(|&x| (-f64::from(x))).collect()
    }

    /// Unit vector
    fn vunit(self) -> Vec<f64> {
        if self.is_empty() { return Vec::new();  };
        let m = 1.0 / self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>().sqrt();
        self.iter().map(|&x| m*(f64::from(x))).collect() 
    }

    /// Arithmetic mean 
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.amean().unwrap(),7.5_f64);
    /// ```
    fn amean(self) -> Result<f64,RE> {
        let n = self.len();
        if n > 0 {
            Ok(self.iter().map(|&x| f64::from(x)).sum::<f64>() / (n as f64)) }
        else { Err(RError::NoDataError("empty self vec")) }    
    }

    /// Arithmetic mean and (population) standard deviation 
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.ameanstd().unwrap();
    /// assert_eq!(res.centre,7.5_f64);
    /// assert_eq!(res.dispersion,4.031128874149275_f64);
    /// ```
    fn ameanstd(self) -> Result<MStats,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); };
        let nf = n as f64;
        let mut sx2 = 0_f64;
        let mean = self
            .iter()
            .map(|&x| {
                let fx = f64::from(x);
                sx2 += fx * fx;
                fx
            })
            .sum::<f64>()/nf; 
        Ok(MStats {
            centre: mean,
            dispersion: (sx2/nf - mean.powi(2)).sqrt()
        })
    }

    /// Linearly weighted arithmetic mean of an f64 slice.     
    /// Linearly ascending weights from 1 to n.    
    /// Time dependent data should be in the order of time increasing.
    /// Then the most recent gets the most weight.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.awmean().unwrap(),9.666666666666666_f64);
    /// ```
    fn awmean(self) -> Result<f64,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); };
        let mut iw = 0_f64; // descending linear weights
        Ok(self.iter()
            .map(|&x| {
                iw += 1_f64;
                iw * f64::from(x)
            }).sum::<f64>() / (sumn(n) as f64))
    }

    /// Linearly weighted arithmetic mean and standard deviation of an f64 slice.    
    /// Linearly ascending weights from 1 to n.    
    /// Time dependent data should be in the order of time increasing.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.awmeanstd().unwrap();
    /// assert_eq!(res.centre,9.666666666666666_f64);
    /// assert_eq!(res.dispersion,3.399346342395192_f64);
    /// ```
    fn awmeanstd(self) -> Result<MStats,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); }; 
        let mut sx2 = 0_f64;
        let mut w = 0_f64; // descending linear weights
        let nf = sumn(n) as f64;
        let centre = self
            .iter()
            .map(|&x| {
                let fx = f64::from(x);
                w += 1_f64;
                let wx = w * fx;
                sx2 += wx * fx;            
                wx
            })
            .sum::<f64>()/ nf;
        Ok(MStats { centre, dispersion:(sx2 / nf - centre.powi(2)).sqrt()})
    }

    /// Harmonic mean of an f64 slice.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.hmean().unwrap(),4.305622526633627_f64);
    /// ```
    fn hmean(self) -> Result<f64,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); };    
        let mut sum = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError("attempt to divide by zero")); };      
            sum += 1.0 / fx
        }
        Ok(n as f64 / sum)
    }

    /// Harmonic mean and standard deviation 
    /// std is based on reciprocal moments
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.hmeanstd().unwrap();
    /// assert_eq!(res.centre,4.305622526633627_f64);
    /// assert_eq!(res.dispersion,1.1996764516690959_f64);
    /// ```
    fn hmeanstd(self) -> Result<MStats,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); };
        let nf = n as f64;
        let mut sx2 = 0_f64;        
        let mut sx = 0_f64;
        for &x in self { 
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError("attempt to divide by zero")); };   
            let rx = 1_f64/fx;  // work with reciprocals
            sx2 += rx*rx;
            sx += rx;   
            }; 
        let recipmean = sx/nf;
        Ok(MStats {
            centre: 1.0/recipmean,
            dispersion: ((sx2/nf-recipmean.powi(2))/(recipmean.powi(4))/nf).sqrt()
        })
    }
    /// Linearly weighted harmonic mean of an f64 slice.    
    /// Linearly ascending weights from 1 to n.    
    /// Time dependent data should be ordered by increasing time.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.hwmean().unwrap(),7.5_f64);
    /// ```
    fn hwmean(self) -> Result<f64,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); };
        let mut sum = 0_f64;
        let mut w = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError("attempt to divide by zero")); };
            w += 1_f64;
            sum += w / fx;
        }
        Ok(sumn(n) as f64 / sum) // reciprocal of the mean of reciprocals
    }

    /// Weighted harmonic mean and standard deviation 
    /// std is based on reciprocal moments
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.hmeanstd().unwrap();
    /// assert_eq!(res.centre,4.305622526633627_f64);
    /// assert_eq!(res.dispersion,1.1996764516690959_f64);
    /// ```
    fn hwmeanstd(self) -> Result<MStats,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); };
        let nf = sumn(n) as f64;
        let mut sx2 = 0_f64;
        let mut sx = 0_f64;
        let mut w = 0_f64;        
        for &x in self {
                w += 1_f64;
                let fx = f64::from(x);
                if !fx.is_normal() { return Err(RError::ArithError("attempt to divide by zero")); };
                sx += w/fx;  // work with reciprocals
                sx2 += w/(fx*fx);
            }; 
        let recipmean = sx/nf; 
        Ok(MStats {
            centre: 1.0/recipmean,
            dispersion: ((sx2/nf-recipmean.powi(2))/(recipmean.powi(4))/nf).sqrt()
        })
    }

    /// Geometric mean of an i64 slice.  
    /// The geometric mean is just an exponential of an arithmetic mean
    /// of log data (natural logarithms of the data items).  
    /// The geometric mean is less sensitive to outliers near maximal value.  
    /// Zero valued data is not allowed!
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.gmean().unwrap(),6.045855171418503_f64);
    /// ```
    fn gmean(self) -> Result<f64,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); }; 
        let mut sum = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError("gmean attempt to take ln of zero")); };        
            sum += fx.ln()
        }
        Ok((sum / (n as f64)).exp())
    }

    /// Geometric mean and std ratio of an f64 slice.  
    /// Zero valued data is not allowed.  
    /// Std of ln data becomes a ratio after conversion back.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.gmeanstd().unwrap();
    /// assert_eq!(res.centre,6.045855171418503_f64);
    /// assert_eq!(res.dispersion,2.1084348239406303_f64);
    /// ```
    fn gmeanstd(self) -> Result<MStats,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); }; 
        let mut sum = 0_f64;
        let mut sx2 = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError("gmeanstd attempt to take ln of zero")); };
            let lx = fx.ln();
            sum += lx;
            sx2 += lx * lx
        }
        sum /= n as f64;
        Ok(MStats {
            centre: sum.exp(),
            dispersion: (sx2 / (n as f64) - sum.powi(2)).sqrt().exp(),
        })
    }

    /// Linearly weighted geometric mean of an i64 slice.  
    /// Ascending weights from 1 down to n.    
    /// Time dependent data should be in time increasing order.  
    /// The geometric mean is an exponential of an arithmetic mean
    /// of log data (natural logarithms of the data items).  
    /// The geometric mean is less sensitive to outliers near maximal value.  
    /// Zero valued data is not allowed!
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.gwmean().unwrap(),8.8185222496341_f64);
    /// ```
    fn gwmean(self) -> Result<f64,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); }; 
        let mut w = 0_f64; // ascending weights
        let mut sum = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError("gwmean attempt to take ln of zero")); }; 
            w += 1_f64;
            sum += w * fx.ln();

        }
        Ok((sum/sumn(n) as f64).exp())
    }

    /// Linearly weighted version of gmeanstd.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.gwmeanstd().unwrap();
    /// assert_eq!(res.centre,8.8185222496341_f64);
    /// assert_eq!(res.dispersion,1.626825493266009_f64);
    /// ```
    fn gwmeanstd(self) -> Result<MStats,RE> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError("empty self vec")); }; 
        let mut w = 0_f64; // ascending weights
        let mut sum = 0_f64;
        let mut sx2 = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError("gwmeanstd attempt to take ln of zero")); }; 
            let lnx = fx.ln();
            w += 1_f64;
            sum += w * lnx;
            sx2 += w * lnx * lnx; 
        }
        let nf = sumn(n) as f64;
        sum /= nf;
        Ok(MStats {
            centre: sum.exp(),
            dispersion: (sx2 as f64 / nf - sum.powi(2)).sqrt().exp(),
        })
    }

    /// Zero median data produced by subtracting the median.
    /// Analogous to zero mean data when subtracting the mean.
    fn zeromedian(self) -> Vec<f64> {
        let median = self.median(); 
        self.iter().map(|&s| f64::from(s)-median).collect()
    }

    /// Probability density function of a sorted slice with repeats. 
    /// Repeats are counted and removed
    fn pdf(self) -> Vec<f64> {     
        let nf = self.len() as f64;   
        let mut res:Vec<f64> = Vec::new();  
        let mut count = 1_usize; // running count
        let mut lastval = self[0];
        self.iter().skip(1).for_each(|&s| (
            if s > lastval { // new value encountered
                res.push((count as f64)/nf); // save previous probability
                lastval = s; // new value
                count = 1_usize; // reset counter
            } else { count += 1; }));
        res.push((count as f64)/nf);  // flush the rest!
        res
    } 

    /// Information (entropy) (in nats)
    fn entropy(self) -> f64 {
        let pdfv = self.sortm(true).pdf();      
        pdfv.iter().map(|&x| -x*(x.ln()) ).sum()                 
    }

    /// (Auto)correlation coefficient of pairs of successive values of (time series) f64 variable.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.autocorr().unwrap(),0.9984603532054123_f64);
    /// ```
    fn autocorr(self) -> Result<f64,RE> {
        let (mut sx, mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64, 0_f64);
        let n = self.len();
        if n < 2 { return Err(RError::NoDataError("autocorr needs a Vec of at least two items")); }; 
        let mut x = f64::from(self[0]);    
        self.iter().skip(1).for_each(|&si| {
            let y = f64::from(si);
            sx += x;
            sy += y;
            sxy += x * y;
            sx2 += x * x;
            sy2 += y * y;
            x = y
        });        
        let nf = n as f64;
        Ok((sxy - sx / nf * sy) / ((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt())
    }

    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Result<Vec<f64>,RE> {
        let mm = self.minmax();
        let range = f64::from(mm.max)-f64::from(mm.min);
        if range == 0_f64 { return  Err(RError::ArithError("lintrans self has zero range")); };
        Ok(self.iter().map(|&x|(f64::from(x)-f64::from(mm.min))/range).collect())        
    }
 
    /// Householder reflector
    fn house_reflector(self) -> Vec<f64> {
        let norm = self.vmag();
        if norm.is_normal() {
            let mut u = self.smult::<f64>(1./norm); 
            if u[0] < 0. { u[0] -= 1.; } else { u[0] += 1.; };
            let uzero = 1.0/(u[0].abs().sqrt());
            u.mutsmult::<f64>(uzero);
            return u;
        };
        let mut u = vec![0.;self.len()];
        u[0] = std::f64::consts::SQRT_2;
        u
    }
}
