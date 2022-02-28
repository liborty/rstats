use crate::{ wsum,MStats,Med,Stats};
use anyhow::{ensure, Result};
// use std::ops::Sub;
pub use indxvec::{here,Printing,merge::{sortm,minmax}};          

impl<T> Stats for &[T] 
    where T: Copy+PartialOrd+std::fmt::Display, // +Sub::<Output = T>,
        f64: From<T> {  

    /// Vector magnitude
    fn vmag(self) -> f64 {
        self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>().sqrt()
    }

    /// Vector magnitude squared (sum of squares)
    fn vmagsq(self) -> f64  {
        self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>()
    }

    /// Vector with reciprocal components
    fn vreciprocal(self) -> Result<Vec<f64>> {
        for &component in self {
           let c = f64::from(component); 
           ensure!(c.is_normal(),
            "{} reciprocal not allowed for zero components!\n{}\n",here!(),self.gr()); 
        }     
        Ok( self.iter().map(|&x| 1.0/(f64::from(x))).collect() )     
    }

    /// Vector with inverse magnitude
    fn vinverse(self) -> Result<Vec<f64>> {
        let mag = self.vmagsq();
        ensure!(mag > 0.0,"{} zero vector can not be inverted!",here!());     
        Ok( self.iter().map(|&x| f64::from(x)/mag).collect() )     
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

    /// Arithmetic mean 
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

    /// Arithmetic mean and (population) standard deviation 
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
    fn awmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut iw = 0_f64; // descending linear weights
        Ok(self.iter()
            .map(|&x| {
                iw += 1_f64;
                iw * f64::from(x)
            })
            .sum::<f64>()
            / wsum(n))
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
    fn awmeanstd(self) -> Result<MStats> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());        
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
            / wsum(n);
        Ok(MStats { mean,std:(sx2 / wsum(n) - mean.powi(2)).sqrt()})
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
    fn hmeanstd(self) -> Result<MStats> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!",here!());
        let nf = n as f64;
        let mut sx2 = 0_f64;        
        let rmean = self
            .iter()
            .map(|&x| {
                let fx = f64::from(x);
                assert!(fx.is_normal(),"{} does not accept zero valued data!",here!());     
                let rx = 1_f64/fx;  // work with reciprocals
                sx2 += rx * rx;
                rx   
            }).sum::<f64>()/nf;    
        Ok(MStats {
            mean: 1.0/rmean,
            std: ((sx2/nf-rmean.powi(2))/(nf*rmean.powi(4))).sqrt()
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
    fn hwmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut sum = 0_f64;
        let mut w = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            ensure!(fx.is_normal(),"{} does not accept zero valued data!",here!());
            w += 1_f64;
            sum += w / fx;
        }
        Ok(wsum(n) / sum)
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
    fn hwmeanstd(self) -> Result<MStats> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!",here!());
        let nf = wsum(n);
        let mut sx2 = 0_f64;
        let mut w = 0_f64;        
        let sx = self
            .iter()
            .map(|&x| {
                w += 1_f64;
                let fx = f64::from(x);
                if !fx.is_normal() { panic!("{} does not accept zero valued data!",here!()) };     
                let rx = w/fx;  // work with reciprocals
                sx2 += w/(fx*fx); 
                rx   
            }).sum::<f64>()/nf;  
        Ok(MStats {
            mean: 1.0/sx,
             std: ((sx2/nf-sx.powi(2))/(nf*sx.powi(4))).sqrt() 
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
    fn gwmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut w = 0_f64; // ascending weights
        let mut sum = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            ensure!(fx.is_normal(),"{} does not accept zero valued data!",here!());
            w += 1_f64;
            sum += w * fx.ln();

        }
        Ok((sum/wsum(n)).exp())
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
    fn gwmeanstd(self) -> Result<MStats> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!()); 
        let mut w = 0_f64; // ascending weights
        let mut sum = 0_f64;
        let mut sx2 = 0_f64;
        for &x in self {
            let fx = f64::from(x);
            ensure!(fx.is_normal(),"{} does not accept zero valued data!",here!());
            let lnx = fx.ln();
            w += 1_f64;
            sum += w * lnx;
            sx2 += w * lnx * lnx; 
        }
        sum /= wsum(n);
        Ok(MStats {
            mean: sum.exp(),
            std: (sx2 as f64 / wsum(n) - sum.powi(2)).sqrt().exp(),
        })
    }

    /// Median of a &[T] slice by sorting
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
        let mut result = Med { median: 
            if 2*mid < gaps { (f64::from(v[mid]) + f64::from(v[mid + 1])) / 2.0 }
            else { f64::from(v[mid]) }, ..Default::default() };
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

    /// MAD median absolute deviation: data spread estimator that is more stable than variance
    /// and more precise than quartiles
    fn mad(self) -> Result<f64> {
        let Med{median,..} = self.median() // ignore quartile fields
            .unwrap_or_else(|_|panic!("{} failed to obtain median",here!()));
        let diffs:Vec<f64> = self.iter().map(|&s| (f64::from(s)-median).abs()).collect();
        let Med{median:res,..} = diffs.median()
            .unwrap_or_else(|_|panic!("{} failed to obtain median",here!()));
        Ok(res)
    }

    /// Zero median data produced by subtracting the median.
    /// Analogous to zero mean data when subtracting the mean.
    fn zeromedian(self) -> Result<Vec<f64>> {
        let Med{median,..} = self.median() // ignore quartile fields
            .unwrap_or_else(|_|panic!("{} failed to obtain median",here!()));
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
        let pdfv = sortm(self,true).pdf();      
        pdfv.iter().map(|&x| -x*(x.ln()) ).sum()                 
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
        (sxy - sx / nf * sy) / ((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()
    }
    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Vec<f64> {
        let mm = minmax(self);
        let range = f64::from(mm.max)-f64::from(mm.min);
        self.iter().map(|&x|(f64::from(x)-f64::from(mm.min))/range).collect()        
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
