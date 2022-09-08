// use core::slice::SlicePattern;

use crate::{error::RError,RE,Stats,Vecg, TriangMat};
use indxvec::{Indices,Vecops};

impl<T> Vecg for &[T] 
    where 
        T: Copy+PartialOrd+Into<T>+std::fmt::Display, f64:From<T> {

     /// Scalar addition to a vector, creates new vec
     fn sadd<U>(self, s:U) -> Vec<f64> 
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> { 
        let sf = f64::from(s);
        self.iter().map(|&x| sf+(f64::from(x))).collect()
     }

    /// Scalar addition to a vector, creates new vec
     fn smult<U>(self, s:U) -> Vec<f64>
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> { 
        let sf = f64::from(s);
        self.iter().map(|&x| sf*(f64::from(x))).collect()
     }

    /// Scalar product.   
    /// Must be of the same length - no error checking (for speed)
    fn dotp<U>(self, v: &[U]) -> f64 
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)*f64::from(vi)).sum::<f64>()
    }

    /// Product with Tukeyvec of hemispheric counts. Self is unitised, only using its direction. 
    /// It should have had subtracted the same reference point as was used in the construction of tukeyvec,
    /// typically the geometric median.
    /// Similar result could be obtained by projecting onto it all points but that is much slower.
    fn dottukey(self, tukey:&[f64]) -> Result<f64,RE> {
        let dims = self.len();
        if 2*dims != tukey.len() { return Err(RError::DataError("{dottukey: tukey vec must have double the dimensions!")); }
        let mut ressum = 0_f64;
        for (i,&scomp) in self.vunit().iter().enumerate() {
            if scomp > 0_f64 { ressum += scomp*tukey[i]; continue };
            if scomp < 0_f64 { ressum -= scomp*tukey[dims+i]; };
        }
        Ok(ressum)
    }   
    
    /// Cosine of angle between the two slices. 
    /// Done in one iteration for efficiency.
    fn cosine<U>(self, v:&[U]) -> f64 
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        let (mut sxy, mut sy2) = (0_f64, 0_f64);
        let sx2: f64 = self
            .iter()
            .zip(v)
            .map(|(&tx, &uy)| {
                let x  = f64::from(tx);
                let y =  f64::from(uy);           
                sxy += x * y;
                sy2 += y * y;
                x*x
            })
            .sum();
        sxy / (sx2*sy2).sqrt()
    }

    /// Vector subtraction 
    fn vsub<U>(self, v:&[U]) -> Vec<f64> 
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)-f64::from(vi)).collect()
    }

    /// Vectors difference unitised (done together for efficiency)
    fn vsubunit<U>(self, v: &[U]) -> Vec<f64>
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        let mut sumsq = 0_f64;
        let dif = self.iter().zip(v).map(
            |(&xi, &vi)| { let d = f64::from(xi) - f64::from(vi); sumsq += d*d; d }
        ).collect::<Vec<f64>>();
        dif.smult(1_f64/sumsq.sqrt())
    } 

    /// Vector addition
    fn vadd<U>(self, v:&[U]) -> Vec<f64>
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)+f64::from(vi)).collect()
    }

    /// Euclidian distance   
    fn vdist<U>(self, v:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (f64::from(xi)-f64::from(vi)).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Weighted arithmetic mean of `self:&[T]`, scaled by `ws:&[U]`
    fn wvmean<U>(self,ws:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        let mut wsum:f64 = 0.;
        let mut sum:f64 = 0.;
        for (&s,&w) in self.iter().zip(ws) { 
            let fw = f64::from(w);
            sum += fw*(f64::from(s));
            wsum += fw;
        };
        sum/wsum
    } 

    /// Weighted distance of `self:&[T]` to `v:&[V]`, scaled by `ws:&[U]`
    /// allows all three to be of different types 
    fn wvdist<U,V>(self,ws:&[U],v:&[V]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U>, 
              V: Copy, f64:From<V> {
        self.iter().enumerate() 
            .map(|(i, &xi)| (f64::from(ws[i])*(f64::from(xi)-f64::from(v[i])).powi(2)))
            .sum::<f64>()
            .sqrt()
    } 

    /// Euclidian distance squared  
    fn vdistsq<U>(self, v:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (f64::from(xi)-f64::from(vi)).powi(2))
            .sum::<f64>()
    }

    /// cityblock distance
    fn cityblockd<U>(self, v:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)|  (f64::from(xi)-f64::from(vi)).abs()) 
            .sum::<f64>()      
    }

    /// Magnitude of the cross product |a x b| = |a||b|sin(theta).
    /// Attains maximum `|a|.|b|` when the vectors are orthogonal.
    fn varea<U>(self, v:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        (self.vmagsq()*v.vmagsq() - self.dotp(v).powi(2)).sqrt()
    }

    /// Area of swept arc 
    /// = |a||b|(1-cos(theta)) = 2|a||b|D
    fn varc<U>(self, v:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {        
        ( v.vmagsq() * self.vmagsq() ).sqrt() - self.dotp(v)
    } 

    /// Positive dotp [0,2|a||b|]
    /// = |a||b|(1+cos(theta)) = 2|a||b|S
    fn pdotp<U>(self, v:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> { 
        ( v.vmagsq() * self.vmagsq() ).sqrt() + self.dotp(v) 
    }

    /// We define vector similarity S in the interval [0,1] as
    /// S = (1+cos(theta))/2
    fn vsim<U>(self, v:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> { 
        (1.0+self.cosine(v))/2.0 }

    /// We define vector dissimilarity D in the interval [0,1] as
    /// D = 1-S = (1-cos(theta))/2
    fn vdisim<U>(self, v:&[U]) -> f64
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        (1.0-self.cosine(v))/2.0 }

    /// Lower triangular covariance matrix for a single vector. 
    /// Where m is either mean or median vector (to be subtracted).
    /// Covariance matrix is symmetric (positive semi definite). 
    fn covone<U>(self, m:&[U]) -> TriangMat 
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> { 
        let mut cov:Vec<f64> = Vec::new(); // flat lower triangular result array
        let vm = self.vsub(m); // zero mean/median vector
        vm.iter().enumerate().for_each(|(i,&thisc)|
            // generate its products up to and including the diagonal (itself)
            vm.iter().take(i+1).for_each(|&component| cov.push(thisc*component)) );
        TriangMat{ trans:false,symmetric:true, data:cov }
    }

    /// Kronecker product of two vectors.   
    /// The indexing is always assumed to be in this order: row,column. 
    fn kron<U>(self, m:&[U]) -> Vec<f64>
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> { 
        let mut krn:Vec<f64> = Vec::new(); // result vector 
        for &a in self {
            for &b in m { krn.push(f64::from(a)*f64::from(b)) }
        }
        krn
    }

    /// Outer product of two vectors.   
    /// The indexing is always assumed to be in this order: row,column. 
    fn outer<U>(self, m:&[U]) -> Vec<Vec<f64>>
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> { 
        let mut out:Vec<Vec<f64>> = Vec::new(); // result vector 
        for &s in self { out.push(m.smult(s)) }      
        out
    }

    /// Joint probability density function of two pairwise matched slices 
    fn jointpdf<U>(self,v:&[U]) -> Result<Vec<f64>,RE> 
        where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {     
        let n = self.len();
        if v.len() != n { 
            return Err(RError::DataError("{jointpdf argument vectors must be of equal length!"));
            }; 
        let nf = n as f64;              
        let mut res:Vec<f64> = Vec::new();
        // collect successive pairs, upgrading all end types to common f64
        let mut spairs:Vec<Vec<f64>> = self.iter().zip(v).map(|(&si,&vi)|
            vec![f64::from(si),f64::from(vi)]).collect(); 
        // sort them to group all same pairs together for counting    
        spairs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap()); 
        let mut count = 1_usize; // running count
        let mut lastindex = 0;
        spairs.iter().enumerate().skip(1).for_each(|(i,si)| 
            if si > &spairs[lastindex] { // new pair encountered
                res.push((count as f64)/nf); // save previous probability
                lastindex = i; // current index becomes the new one
                count = 1_usize; // reset counter
            } else { count += 1; }); 
        res.push((count as f64)/nf);  // flush the rest!
        Ok(res)
    } 

    /// Joint entropy of two sets of the same length
    fn jointentropy<U>(self, v:&[U]) -> Result<f64,RE>
    where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        let jpdf = self.jointpdf(v)?;
        Ok(jpdf.iter().map(|&x|-x*(x.ln())).sum()) 
    }

    /// Dependence of &[T] &[U] variables in the range [0,1]
    /// returns 0 iff they are statistically component wise independent
    /// returns 1 when they are identical or all their values are unique
    fn dependence<U>(self, v:&[U]) -> Result<f64,RE>
    where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {   
        Ok((self.entropy() + v.entropy())/self.jointentropy(v)?-1.0) 
    }

    /// Independence of &[T] &[U] variables in the range [0,1]
    /// returns 1 iff they are statistically component wise independent
    fn independence<U>(self, v:&[U]) -> Result<f64,RE>
    where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {   
        Ok(2.0 * self.jointentropy(v)? / (self.entropy() + v.entropy())-1.0)
    }

    /// We define median based correlation as cosine of an angle between two
    /// zero median vectors (analogously to Pearson's zero mean vectors) 
    /// # Example
    /// ```
    /// use rstats::Vecg;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.correlation(&v2),-0.1076923076923077);
    /// ```
    fn mediancorr<U>(self, v: &[U]) -> f64
    where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        // let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let zeroself = self.zeromedian();
        let zerov = v.zeromedian();
        zeroself.cosine(&zerov)
    }        

    /// Pearson's (most common) correlation. 
    /// # Example
    /// ```
    /// use rstats::Vecg;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.correlation(&v2),-0.1076923076923077);
    /// ```
    fn correlation<U>(self, v: &[U]) -> f64
    where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let sx: f64 = self
            .iter()
            .zip(v)
            .map(|(&xt, &yu)| {
                let x = f64::from(xt);
                let y = f64::from(yu);
                sy += y;
                sxy += x * y;
                sx2 += x * x;
                sy2 += y * y;
                x
            })
            .sum();
        let nf = self.len() as f64;
        (sxy - sy*sx / nf) / ((sx2 - sx*sx / nf) * (sy2 - sy*sy / nf)).sqrt()
    }
    /// Kendall Tau-B correlation.
    /// Defined by: tau = (conc - disc) / sqrt((conc + disc + tiesx) * (conc + disc + tiesy))
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vecg;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.kendalcorr(&v2),-0.07692307692307693);
    /// ```
    fn kendalcorr<U>(self, v:&[U]) -> f64
    where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        let (mut conc, mut disc, mut tiesx, mut tiesy) = (0_i64, 0_i64, 0_i64, 0_i64); 
        for i in 1..self.len() {
            let x = f64::from(self[i]);
            let y = f64::from(v[i]);
            for j in 0..i {
                let xd = x - f64::from(self[j]);
                let yd = y - f64::from(v[j]);
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
        (conc - disc) as f64 / (((conc + disc + tiesx) * (conc + disc + tiesy)) as f64).sqrt()
    }
    /// Spearman rho correlation.
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vecg;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.spearmancorr(&v2),-0.1076923076923077);
    /// ```
    fn spearmancorr<U>(self, v: &[U]) -> f64
    where U: Copy+PartialOrd+Into<U>+std::fmt::Display, f64:From<U> {
        let xvec = self.rank(true);
        let yvec = v.rank(true); // rank from crate idxvec::merge
        // It is just Pearson's correlation of usize ranks
        xvec.ucorrelation(&yvec) // using Indices trait from idxvec
    }

    /// Change to gm that adding point p will cause
    fn contribvec_newpt(self,gm:&[f64],recips:f64) -> Vec<f64> {
        let dv = self.vsub::<f64>(gm);
        let mag = dv.vmag();
        if !mag.is_normal() { return dv; }; 
        let recip = 1f64/mag; // first had to test for division by zero
        // adding new unit vector (to approximate zero vector)
        dv.smult::<f64>(recip/(recips+recip)) // to unit v. and scaling by new sum of reciprocals 
    }
    
    /// Magnitude of change to gm that adding point p will cause
    fn contrib_newpt(self,gm:&[f64],recips:f64) -> f64 {
        let mag = self.vdist::<f64>(gm);
        if !mag.is_normal() { return 0_f64; }; 
        let recip = 1f64/mag; // first had to test for division by zero
        1_f64 / (recips + recip)
    }    
    
    /// Contribution an existing set point p has made to the gm
    fn contribvec_oldpt(self,gm:&[f64],recips:f64) -> Vec<f64> {
        let dv = self.vsub::<f64>(gm);
        let mag = dv.vmag();
        if !mag.is_normal() { return dv; };
        let recip = 1f64/mag; // first had to test for division by zero 
        dv.smult::<f64>(recip/(recip - recips)) // scaling
    }
        
    /// Contribution removing an existing set point p will make
    /// Is a negative number
    fn contrib_oldpt(self,gm:&[f64],recips:f64) -> f64 {
        let mag = self.vdist::<f64>(gm);
        if !mag.is_normal() { return 0_f64; }; 
        let recip = 1f64/mag; // first had to test for division by zero
        1_f64 / (recip - recips) 
        // self.contribvec_oldpt(gm,recips,p).vmag()
    } 

    /// Householder reflect
    fn house_reflect<U>(self,x:&[U]) -> Vec<f64>
    where U: Copy+PartialOrd+std::fmt::Display, f64:From<U> {
        x.vsub(&self.smult(x.dotp(self)))
    }

}
