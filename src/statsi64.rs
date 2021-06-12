use crate::{MStats, Med, Stats, functions::wsum, here};
use anyhow::{ensure, Result};

impl Stats for &[i64] {
    /// Arithmetic mean of an i64 slice
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1:Vec<i64> = vec![1,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// assert_eq!(v1.as_slice().amean().unwrap(),7.5_f64);
    /// ```
    fn amean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        Ok(self.iter().map(|&x| x as f64).sum::<f64>() / (n as f64))
    }

    /// Arithmetic mean and standard deviation of an i64 slice
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// let res = v1.as_slice().ameanstd().unwrap();
    /// assert_eq!(res.mean,7.5_f64);
    /// assert_eq!(res.std,4.031128874149275_f64);
    /// ```
    fn ameanstd(self) -> Result<MStats> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut sx2 = 0_f64;
        let mean = self
            .iter()
            .map(|&x| {
                let lx = x as f64;
                sx2 += lx * lx;
                lx
            })
            .sum::<f64>()
            / (n as f64);
        Ok(MStats {
            mean: mean,
            std: (sx2 / (n as f64) - mean.powi(2)).sqrt(),
        })
    }

    /// Linearly weighted arithmetic mean of an i64 slice.     
    /// Linearly descending weights from n down to one.    
    /// Time dependent data should be in the stack order - the last being the oldest.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// assert_eq!(v1.as_slice().awmean().unwrap(),5.333333333333333_f64);
    /// ```
    fn awmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut w = (n + 1) as f64; // descending linear weights
        Ok(self
            .iter()
            .map(|&x| {
                w -= 1.;
                w * x as f64
            })
            .sum::<f64>()
            / wsum(n))
    }

    /// Liearly weighted arithmetic mean and standard deviation of an i64 slice.    
    /// Linearly descending weights from n down to one.    
    /// Time dependent data should be in the stack order - the last being the oldest.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
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
                let lx = x as f64;
                let wx = w * lx;
                sx2 += wx * lx;
                w -= 1.;
                wx
            })
            .sum::<f64>() as f64
            / wsum(n);
        Ok(MStats {
            mean: mean,
            std: (sx2 / wsum(n) - mean.powi(2)).sqrt(),
        })
    }

    /// Harmonic mean of an i64 slice.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// assert_eq!(v1.as_slice().hmean().unwrap(),4.305622526633627_f64);
    /// ```
    fn hmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut sum = 0_f64;
        for &x in self {
            ensure!(x!=0_i64,"{} does not accept zero valued data!",here!());
            sum += 1.0 / (x as f64)
        }
        Ok(n as f64 / sum)
    }

    /// Linearly weighted harmonic mean of an i64 slice.    
    /// Linearly descending weights from n down to one.    
    /// Time dependent data should be in the stack order - the last being the oldest.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// assert_eq!(v1.as_slice().hwmean().unwrap(),3.019546395306663_f64);
    /// ```
    fn hwmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut sum = 0_f64;
        let mut w = n as f64;
        for &x in self {
            ensure!(x!=0_i64,"{} does not accept zero valued data!",here!());
            sum += w / x as f64;
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
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// assert_eq!(v1.as_slice().gmean().unwrap(),6.045855171418503_f64);
    /// ```
    fn gmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut sum = 0_f64;
        for &x in self {
            ensure!(x!=0_i64,"{} does not accept zero valued data!",here!());
            sum += (x as f64).ln()
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
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
    /// assert_eq!(v1.as_slice().gwmean().unwrap(),4.144953510241978_f64);
    /// ```
    fn gwmean(self) -> Result<f64> {
        let n = self.len();
        ensure!(n > 0, "{} sample is empty!", here!());
        let mut w = n as f64; // descending weights
        let mut sum = 0_f64;
        for &x in self {
            ensure!(x!=0_i64,"{} does not accept zero valued data!",here!());
            sum += w * (x as f64).ln();
            w -= 1_f64;
        }
        Ok((sum / wsum(n)).exp())
    }

    /// Geometric mean and std ratio of an i64 slice.  
    /// Zero valued data is not allowed.  
    /// Std of ln data becomes a ratio after conversion back.
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
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
            ensure!(x!=0_i64,"{} does not accept zero valued data!",here!());
            let lx = (x as f64).ln();
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
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
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
            ensure!(x!=0_i64,"{} does not accept zero valued data!",here!());
            let lnx = (x as f64).ln();
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

    /// Median. 
    /// # Example
    /// ```
    /// use rstats::Stats;
    /// let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
    /// let res = v1.as_slice().median().unwrap();
    /// assert_eq!(res.median,8.);
    /// assert_eq!(res.lquartile,4.5);
    /// assert_eq!(res.uquartile,11.5);
    /// ```
    fn median(self) -> Result<Med> {
        let gaps = self.len()-1;
        let mid = gaps / 2;
        let quarter = gaps / 4;
        let threeq = 3 * gaps / 4;
        let qrem = gaps % 4;
        let mut v:Vec<i64> = self.to_vec();
        v.sort_unstable();     
        let mut result: Med = Default::default();
        result.median = if 2*mid < gaps { (v[mid] + v[mid + 1]) as f64 / 2.0 }
            else { v[mid] as f64 };
        if qrem == 0 {
            result.lquartile = v[quarter] as f64;
            result.uquartile = v[threeq] as f64;
            return Ok(result) };
        if qrem == 1 {
            result.lquartile = (3*v[quarter] + v[quarter+1]) as f64 / 4.;
            result.uquartile = (v[threeq] + 3*v[threeq+1]) as f64 / 4.;
            return Ok(result) };
        if qrem == 2 {
            result.lquartile = (v[quarter]+v[quarter+1]) as f64 / 2.;
            result.uquartile = (v[threeq] + v[threeq+1]) as f64 / 2.;
            return Ok(result) };
        if qrem == 3 {
            result.lquartile = (v[quarter] + 3*v[quarter+1]) as f64 / 4.;
            result.uquartile = (3*v[threeq] + v[threeq+1]) as f64 / 4.
        };
        Ok(result)       
    }    

    /// Returns vector of ranks, 
    /// ranked from the smallest number in self (rank 0) to the biggest (rank n-1).
    /// Equalities lead to fractional ranks (hence Vec<f64> output) and the range of rank values is reduced.
    fn ranks(self) -> Result<Vec<f64>> {
        let n = self.len();
        let mut rank = vec![0_f64; n];
        // make each of n*(n-1)/2 comparisons just once
        for i in 1..n {
            let x = self[i];
            for j in 0..i {
                if x > self[j] {
                    rank[i] += 1_f64; // demoting i
                    continue;
                };
                if x < self[j] {
                    rank[j] += 1_f64; // demoting j
                    continue;
                };
                rank[i] += 0.5;
                rank[j] += 0.5;
            }
        }
        Ok(rank)
    }
 
    /// Returns vector of ranks, 
    /// ranked from the smallest number in self (rank 0) to the biggest (rank n-1).
    fn iranks(self) -> Result<Vec<i64>> {
        let n = self.len();
        let mut rank = vec![0_i64; n];
        // make each of n*(n-1)/2 comparisons just once
        for i in 1..n {
            let x = self[i];
            for j in 0..i {
                if x > self[j] {
                    rank[i] += 1_i64; // demoting i
                } else if x < self[j] {
                    rank[j] += 1_i64; // demoting j
                };
                // else items are equal, not demoting any
            }
        }
        Ok(rank)
    } 
}
