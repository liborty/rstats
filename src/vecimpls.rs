use crate::{MutVectors, Stats, Vectors};
use crate::functions::emsg;
use anyhow::{ensure, Result};


impl Vectors for &[f64] {
 
    /// Scalar multiplication of a vector, creates new vec
    fn smult(self, s: f64) -> Vec<f64> {
        self.iter().map(|&x| s*x).collect()
    }

    /// Scalar addition to a vector, creates new vec
    fn sadd(self, s:f64) -> Vec<f64> {
        self.iter().map(|&x| s+x).collect()
    }

    /// Scalar product of two f64 slices.   
    /// Must be of the same length - no error checking (for speed)
    fn dotp(self, v: &[f64]) -> f64 {
        self.iter().zip(v).map(|(&xi, &vi)| xi * vi).sum::<f64>()
    }

    /// Vector subtraction, creates a new Vec result
    fn vsub(self, v: &[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| xi - vi).collect()
    }

    /// Vector addition, creates a new Vec result
    fn vadd(self, v: &[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| xi + vi).collect()
    }

    /// Euclidian distance between two n dimensional points (vectors).  
    /// Slightly faster than vsub followed by vmag, as both are done in one loop
    fn vdist(self, v: &[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (xi - vi).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Vector magnitude
    fn vmag(self) -> f64 {
        self.iter().map(|&x| x.powi(2)).sum::<f64>().sqrt()
    }

    /// Unit vector - creates a new one
    fn vunit(self) -> Vec<f64> {
        self.smult(1. / self.iter().map(|x| x.powi(2)).sum::<f64>().sqrt())
    }

    /// Correlation coefficient of a sample of two f64 variables.
    /// # Example
    /// ```
    /// use rstats::Vectors;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.];
    /// assert_eq!(v1.correlation(&v2).unwrap(),-1_f64);
    /// ```
    fn correlation(self, v: &[f64]) -> Result<f64> {
        let n = self.len();
        ensure!(
            n > 0,
            emsg(file!(), line!(), "correlation - first sample is empty")
        );
        ensure!(
            n == v.len(),
            emsg(
                file!(),
                line!(),
                "correlation - samples are not of the same size"
            )
        );
        let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let sx: f64 = self
            .iter()
            .zip(v)
            .map(|(&x, &y)| {
                sy += y;
                sxy += x * y;
                sx2 += x * x;
                sy2 += y * y;
                x
            })
            .sum();
        let nf = n as f64;
        Ok((sxy - sx / nf * sy) / (((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()))
    }

    /// Kendall Tau-B correlation coefficient of a sample of two f64 variables.
    /// Defined by: tau = (conc - disc) / sqrt((conc + disc + tiesx) * (conc + disc + tiesy))
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vectors;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.];
    /// assert_eq!(v1.kendalcorr(&v2).unwrap(),-1_f64);
    /// ```
    fn kendalcorr(self, v: &[f64]) -> Result<f64> {
        let n = self.len();
        ensure!(
            n > 0,
            emsg(file!(), line!(), "kendalcorr - first sample is empty")
        );
        ensure!(
            n == v.len(),
            emsg(
                file!(),
                line!(),
                "kendalcorr - samples are not of the same size"
            )
        );
        let (mut conc, mut disc, mut tiesx, mut tiesy) = (0_i64, 0_i64, 0_i64, 0_i64);
        for i in 1..n {
            let x = self[i];
            let y = v[i];
            for j in 0..i {
                let xd = x - self[j];
                let yd = y - v[j];
                if !xd.is_normal() {
                    if !yd.is_normal() {
                        continue;
                    } else {
                        tiesx += 1;
                        continue;
                    }
                };
                if !yd.is_normal() {
                    tiesy += 1;
                    continue;
                };
                if (xd * yd).signum() > 0_f64 {
                    conc += 1
                } else {
                    disc += 1
                }
            }
        }
        Ok((conc - disc) as f64 / (((conc + disc + tiesx) * (conc + disc + tiesy)) as f64).sqrt())
    }
    /// Spearman rho correlation coefficient of a sample of two f64 variables.
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vectors;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.];
    /// assert_eq!(v1.spearmancorr(&v2).unwrap(),-1_f64);
    /// ```
    fn spearmancorr(self, v: &[f64]) -> Result<f64> {
        let n = self.len();
        ensure!(
            n > 0,
            emsg(file!(), line!(), "spearmancorr - first sample is empty")
        );
        ensure!(
            n == v.len(),
            emsg(
                file!(),
                line!(),
                "spearmancorr - samples are not of the same size"
            )
        );
        let xvec = self.ranks().unwrap();
        let yvec = v.ranks().unwrap();
        let mx = xvec.ameanstd().unwrap();
        let my = yvec.ameanstd().unwrap();
        let mut covar = 0_f64;
        for i in 0..n {
            covar += (xvec[i] - mx.mean) * (yvec[i] - my.mean);
        }
        covar /= mx.std * my.std * (n as f64);
        // remove small truncation errors
        if covar > 1.0 {
            covar = 1_f64
        } else if covar < -1_f64 {
            covar = -1.0
        };
        Ok(covar)
    }

    /// (Auto)correlation coefficient of pairs of successive values of (time series) f64 variable.
    /// # Example
    /// ```
    /// use rstats::Vectors;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.autocorr().unwrap(),0.9984603532054123_f64);
    /// ```
    fn autocorr(self) -> Result<f64> {
        let n = self.len();
        ensure!(
            n >= 2,
            emsg(file!(), line!(), "autocorr - sample is too small")
        );
        let (mut sx, mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64, 0_f64);
        for i in 0..n - 1 {
            let x = self[i];
            let y = self[i + 1];
            sx += x;
            sy += y;
            sxy += x * y;
            sx2 += x * x;
            sy2 += y * y
        }
        let nf = n as f64;
        Ok((sxy - sx / nf * sy) / (((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()))
    }

    /// Finds minimum, minimum's index, maximum, maximum's index of &[f64]
    /// Here self is usually some data, rather than a vector
    fn minmax(self) -> (f64, usize, f64, usize) {
        let mut min = self[0]; // initialise to the first value
        let mut mini = 0;
        let mut max = self[0]; // initialised as min, allowing 'else' below
        let mut maxi = 0;
        for i in 1..self.len() {
            let x = self[i];
            if x < min {
                min = x;
                mini = i
            } else if x > max {
                max = x;
                maxi = i
            }
        }
        (min, mini, max, maxi)
    }
    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Vec<f64> {
        let (min,_,max,_) = self.minmax();
        let range = max-min;
        self.iter().map(|&x|(x-min)/range).collect()        
    }
    /// Sorted vector
    fn sortf(self) -> Vec<f64> {
        let mut sorted:Vec<f64> = self.to_vec();
        sorted.mutsortf();
        sorted      
    }
}
