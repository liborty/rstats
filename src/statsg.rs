use crate::{ sumn, MStats, Stats, error::RError };
// use anyhow::{ensure, Result};

use indxvec::{Vecops};
use medians::{Median};    

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
    fn vreciprocal(self) -> Result<Vec<f64>,RError> {
        if self.is_empty() { return  Err(RError::NoDataError)  };
        for &component in self {
           let c = f64::from(component); 
           if !c.is_normal() { return Err(RError::ArithError); }; 
        }     
        Ok( self.iter().map(|&x| 1.0/(f64::from(x))).collect() )     
    }

    /// Vector with inverse magnitude
    fn vinverse(self) -> Result<Vec<f64>,RError> {
        if self.is_empty() { return  Err(RError::NoDataError)  };
        let mag = self.vmagsq();
        if mag > 0.0 {  
            Ok( self.iter().map(|&x| f64::from(x)/mag).collect() ) }
        else { Err(RError::DataError) }    
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
    /// assert_eq!(v1.as_slice().amean().unwrap(),7.5_f64);
    /// ```
    fn amean(self) -> Result<f64,RError> {
        let n = self.len();
        if n > 0 {
            Ok(self.iter().map(|&x| f64::from(x)).sum::<f64>() / (n as f64)) }
        else { Err(RError::NoDataError) }    
    }

    /// Arithmetic mean and (population) standard deviation 
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.as_slice().ameanstd().unwrap();
    /// assert_eq!(res.mean,7.5_f64);
    /// assert_eq!(res.std,4.031128874149275_f64);
    /// ```
    fn ameanstd(self) -> Result<MStats,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); };
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
            mean,
            std: (sx2/nf - mean.powi(2)).sqrt()
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
    /// assert_eq!(v1.as_slice().awmean().unwrap(),9.666666666666666_f64);
    /// ```
    fn awmean(self) -> Result<f64,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); };
        let mut iw = 0_f64; // descending linear weights
        Ok(self.iter()
            .map(|&x| {
                iw += 1_f64;
                iw * f64::from(x)
            }).sum::<f64>() / sumn(n))
    }

    /// Linearly weighted arithmetic mean and standard deviation of an f64 slice.    
    /// Linearly ascending weights from 1 to n.    
    /// Time dependent data should be in the order of time increasing.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.as_slice().awmeanstd().unwrap();
    /// assert_eq!(res.mean,9.666666666666666_f64);
    /// assert_eq!(res.std,3.399346342395192_f64);
    /// ```
    fn awmeanstd(self) -> Result<MStats,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); }; 
        let mut sx2 = 0_f64;
        let mut w = 0_f64; // descending linear weights
        let mean = self
            .iter()
            .map(|&x| {
                let fx = f64::from(x);
                w += 1_f64;
                let wx = w * fx;
                sx2 += wx * fx;            
                wx
            })
            .sum::<f64>()
            / sumn(n);
        Ok(MStats { mean,std:(sx2 / sumn(n) - mean.powi(2)).sqrt()})
    }

    /// Harmonic mean of an f64 slice.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().hmean().unwrap(),4.305622526633627_f64);
    /// ```
    fn hmean(self) -> Result<f64,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); };    
        let mut sum = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError); };      
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
    /// let res = v1.as_slice().hmeanstd().unwrap();
    /// assert_eq!(res.mean,4.305622526633627_f64);
    /// assert_eq!(res.std,1.1996764516690959_f64);
    /// ```
    fn hmeanstd(self) -> Result<MStats,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); };
        let nf = n as f64;
        let mut sx2 = 0_f64;        
        let mut sx = 0_f64;
        for &x in self { 
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError); };   
            let rx = 1_f64/fx;  // work with reciprocals
            sx2 += rx*rx;
            sx += rx;   
            }; 
        let recipmean = sx/nf;
        Ok(MStats {
            mean: 1.0/recipmean,
            std: ((sx2/nf-recipmean.powi(2))/(recipmean.powi(4))/nf).sqrt()
        })
    }
    /// Linearly weighted harmonic mean of an f64 slice.    
    /// Linearly ascending weights from 1 to n.    
    /// Time dependent data should be ordered by increasing time.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.as_slice().hwmean().unwrap(),7.5_f64);
    /// ```
    fn hwmean(self) -> Result<f64,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); };
        let mut sum = 0_f64;
        let mut w = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError); };
            w += 1_f64;
            sum += w / fx;
        }
        Ok(sumn(n) / sum)
    }

    /// Weighted harmonic mean and standard deviation 
    /// std is based on reciprocal moments
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.as_slice().hmeanstd().unwrap();
    /// assert_eq!(res.mean,4.305622526633627_f64);
    /// assert_eq!(res.std,1.1996764516690959_f64);
    /// ```
    fn hwmeanstd(self) -> Result<MStats,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); };
        let nf = sumn(n);
        let mut sx2 = 0_f64;
        let mut sx = 0_f64;
        let mut w = 0_f64;        
        for &x in self {
                w += 1_f64;
                let fx = f64::from(x);
                if !fx.is_normal() { return Err(RError::ArithError); };
                sx += w/fx;  // work with reciprocals
                sx2 += w/(fx*fx);
            }; 
        let recipmean = sx/nf; 
        Ok(MStats {
            mean: 1.0/recipmean,
            std: ((sx2/nf-recipmean.powi(2))/(recipmean.powi(4))/nf).sqrt()
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
    /// assert_eq!(v1.as_slice().gmean().unwrap(),6.045855171418503_f64);
    /// ```
    fn gmean(self) -> Result<f64,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); }; 
        let mut sum = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError); };        
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
    /// let res = v1.as_slice().gmeanstd().unwrap();
    /// assert_eq!(res.mean,6.045855171418503_f64);
    /// assert_eq!(res.std,2.1084348239406303_f64);
    /// ```
    fn gmeanstd(self) -> Result<MStats,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); }; 
        let mut sum = 0_f64;
        let mut sx2 = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError); };
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
    /// assert_eq!(v1.as_slice().gwmean().unwrap(),8.8185222496341_f64);
    /// ```
    fn gwmean(self) -> Result<f64,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); }; 
        let mut w = 0_f64; // ascending weights
        let mut sum = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError); }; 
            w += 1_f64;
            sum += w * fx.ln();

        }
        Ok((sum/sumn(n)).exp())
    }

    /// Linearly weighted version of gmeanstd.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let res = v1.as_slice().gwmeanstd().unwrap();
    /// assert_eq!(res.mean,8.8185222496341_f64);
    /// assert_eq!(res.std,1.626825493266009_f64);
    /// ```
    fn gwmeanstd(self) -> Result<MStats,RError> {
        let n = self.len();
        if n == 0 { return Err(RError::NoDataError); }; 
        let mut w = 0_f64; // ascending weights
        let mut sum = 0_f64;
        let mut sx2 = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            if !fx.is_normal() { return Err(RError::ArithError); }; 
            let lnx = fx.ln();
            w += 1_f64;
            sum += w * lnx;
            sx2 += w * lnx * lnx; 
        }
        sum /= sumn(n);
        Ok(MStats {
            mean: sum.exp(),
            std: (sx2 as f64 / sumn(n) - sum.powi(2)).sqrt().exp(),
        })
    }

    /// Zero median data produced by subtracting the median.
    /// Analogous to zero mean data when subtracting the mean.
    fn zeromedian(self) -> Result<Vec<f64>,RError> {
        let median = self.median(); 
        Ok(self.iter().map(|&s| f64::from(s)-median).collect())
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
    fn autocorr(self) -> Result<f64,RError> {
        let (mut sx, mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64, 0_f64);
        let n = self.len();
        if n < 2 { return Err(RError::NoDataError); }; 
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
    fn lintrans(self) -> Result<Vec<f64>,RError> {
        let mm = self.minmax();
        let range = f64::from(mm.max)-f64::from(mm.min);
        if range == 0_f64 { return  Err(RError::ArithError); };
        Ok(self.iter().map(|&x|(f64::from(x)-f64::from(mm.min))/range).collect())        
    }

    /// Reconstructs the full symmetric square matrix from its lower diagonal compact form,
    /// as produced by covar, covone, wcovar
    fn symmatrix(self) -> Vec<Vec<f64>> {

         fn trseqtosubs(s:usize) -> (usize,usize) { 
            // solution of quadratic equation to find the dimension of the full square matrix
            let row = ((((8*s+1) as f64).sqrt() - 1.)/2.) as usize; // cast truncates, like .floor()
            let column = s - row*(row+1)/2; // subtracting previous triangular number
            (row,column)
        }

        let (n,_) = trseqtosubs(self.len());
        let mut mat = vec![vec![0_f64;n];n]; // create the square matrix 
        self.iter().enumerate().for_each(|(i,&s)| {
            let (row,column) = trseqtosubs(i);
            if row < column { mat[column][row] = f64::from(s) as f64; }; // symmetrical reflection
            // also set values in lower triangular region, including the diagonal
            mat[row][column] = f64::from(s); } ); 
        mat
    }    
 
}
