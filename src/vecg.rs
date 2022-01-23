use crate::{here,Stats,Vecg,Vecf64};
pub use indxvec::{Indices,merge::{sortm,rank}};

impl<T,U> Vecg<T,U> for &[T] 
    where T: Copy+PartialOrd+std::fmt::Display,
        U: Copy+PartialOrd+std::fmt::Display,
        f64: From<T>, f64: From<U> {

     /// Scalar addition to a vector, creates new vec
     fn sadd(self, s:U) -> Vec<f64> {
        let sf = f64::from(s);
        self.iter().map(|&x| sf+(f64::from(x))).collect()
     }

    /// Scalar addition to a vector, creates new vec
     fn smult(self, s:U) -> Vec<f64> {
        let sf = f64::from(s);
        self.iter().map(|&x| sf*(f64::from(x))).collect()
     }

    /// Scalar product.   
    /// Must be of the same length - no error checking (for speed)
    fn dotp(self, v: &[U]) -> f64 {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)*f64::from(vi)).sum::<f64>()
    }
    
    /// Cosine of angle between the two slices. 
    /// Done in one iteration for efficiency.
    fn cosine(self, v:&[U]) -> f64 {
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
    fn vsub(self, v:&[U]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)-f64::from(vi)).collect()
    }

    /// Vectors difference unitised (done together for efficiency)
    fn vsubunit(self, v: &[U]) -> Vec<f64> {
        let mut sumsq = 0_f64;
        let dif = self.iter().zip(v).map(
            |(&xi, &vi)| { let d = f64::from(xi) - f64::from(vi); sumsq += d*d; d }
        ).collect::<Vec<f64>>();
        dif.smult(1_f64/sumsq.sqrt())
    } 

    /// Vector addition
    fn vadd(self, v:&[U]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)+f64::from(vi)).collect()
    }

    /// Euclidian distance   
    fn vdist(self, v:&[U]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (f64::from(xi)-f64::from(vi)).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Euclidian distance squared  
    fn vdistsq(self, v:&[U]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (f64::from(xi)-f64::from(vi)).powi(2))
            .sum::<f64>()
    }

    /// cityblock distance
    fn cityblockd(self, v:&[U]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)|  (f64::from(xi)-f64::from(vi)).abs()) 
            .sum::<f64>()      
    }

    /// Magnitude of the cross product |a x b| = |a||b|sin(theta).
    /// Attains maximum `|a|.|b|` when the vectors are orthogonal.
    fn varea(self, v:&[U]) -> f64 {
        (self.vmagsq()*v.vmagsq() - self.dotp(v).powi(2)).sqrt()
    }

    /// Area of swept arc 
    /// = |a||b|(1-cos(theta)) = 2|a||b|D
    fn varc(self, v:&[U]) -> f64 {        
        ( v.vmagsq() * self.vmagsq() ).sqrt() - self.dotp(v)
    } 

    /// We define vector similarity S in the interval [0,1] as
    /// S = (1+cos(theta))/2
    fn vsim(self, v:&[U]) -> f64 { (1.0+self.cosine(v))/2.0 }

    /// We define vector dissimilarity D in the interval [0,1] as
    /// D = 1-S = (1-cos(theta))/2
    fn vdisim(self, v:&[U]) -> f64 { (1.0-self.cosine(v))/2.0 }

    /// Flattened lower triangular part of a covariance matrix. 
    /// m can be either mean or median vector. 
    /// Since covariance matrix is symmetric (positive semi definite), 
    /// the upper triangular part can be trivially added for all j>i by: c(j,i) = c(i,j).
    /// N.b. the indexing is always assumed to be in this order: row,column.
    /// The items of the resulting lower triangular array c[i][j] are pushed flat
    /// into a single vector in this double loop order: left to right, top to bottom 
    fn covone(self, m:&[U]) -> Vec<f64> { 
        let mut cov:Vec<f64> = Vec::new(); // flat lower triangular result array
        let vm = self.vsub(m); // zero mean vector
        vm.iter().enumerate().for_each(|(i,&thisc)|
            // generate its products up to and including the diagonal (itself)
            vm.iter().take(i+1).for_each(|&component| cov.push(thisc*component)) );
        cov
    }

    /// Kronecker product of two vectors.   
    /// The indexing is always assumed to be in this order: row,column. 
    fn kron(self, m:&[U]) -> Vec<f64> { 
        let mut krn:Vec<f64> = Vec::new(); // result vector 
        for &a in self {
            for &b in m { krn.push(f64::from(a)*f64::from(b)) }
        }
        krn
    }

    /// Outer product of two vectors.   
    /// The indexing is always assumed to be in this order: row,column. 
    fn outer(self, m:&[U]) -> Vec<Vec<f64>> { 
        let mut out:Vec<Vec<f64>> = Vec::new(); // result vector 
        for &s in self { out.push(m.smult(s)) }      
        out
    }
/*
    /// Cholesky decomposition of positive definite matrix into LL^T
    fn cholesky(self) -> Vec<f64> {
    let n = self.len();
    let mut res = vec![0.0; n];
    for i in 0..n {
        for j in 0..(i+1){
            let mut s = 0.0;
            for k in 0..j {
                s += res[i * n + k] * res[j * n + k];
            }
            res[i * n + j] = if i == j { (f64::from(self[i * n + i]) - s).sqrt() } 
            else { 1.0 / res[j * n + j] * (f64::from(self[i * n + j]) - s) };
        }
    }
    res
}
*/
    /// Joint probability density function of two pairwise matched slices 
    fn jointpdf(self,v:&[U]) -> Vec<f64> {     
        let n = self.len();
        if v.len() != n { panic!("{} argument vectors must be of equal length!",here!()) }; 
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
        res
    } 

    /// Joint entropy of two sets of the same length
    fn jointentropy(self, v:&[U]) -> f64 {
        let jpdf = self.jointpdf(v);
        jpdf.iter().map(|&x| -x*(x.ln()) ).sum() 
    }

    /// Independence of &[T] &[U] variables in the range [1,2]
    /// returns 2 iff they are statistically component wise independent
    fn independence(self, v:&[U]) -> f64 {   
        2.0 * self.jointentropy(v) / (self.entropy() + v.entropy()) 
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
    fn mediancorr(self, v: &[U]) -> f64 {
        // let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let zeroself = self.zeromedian()
            .unwrap_or_else(|_| panic!("{} failed to find the median",here!()));
        let zerov = v.zeromedian()
            .unwrap_or_else(|_| panic!("{} failed to find the median",here!()));
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
    fn correlation(self, v: &[U]) -> f64 {
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
    fn kendalcorr(self, v:&[U]) -> f64 {
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
    fn spearmancorr(self, v: &[U]) -> f64 {
        let xvec = rank(self,true);
        let yvec = rank(v,true); // rank from crate idxvec::merge
        // It is just Pearson's correlation of usize ranks
        xvec.ucorrelation(&yvec) // using Indices trait from idxvec
    }
}

impl<T> Vecf64<T> for &[T] 
    where T: Copy+PartialOrd+std::fmt::Display,         
        f64: From<T>, f64: From<T> {

    /// scalar multiplication by f64
    fn smultf64(self, s:f64) -> Vec<f64> { 
        self.iter().map(|&x| s*(f64::from(x))).collect()
     }

    /// Vector subtraction of `&[f64]`
    fn vsubf64(self, v:&[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)-vi).collect()
    }

    /// Vectors difference unitised (done together for efficiency)
    fn vsubunitf64(self, v:&[f64]) -> Vec<f64> {
        let mut sumsq = 0_f64;
        let dif = self.iter().zip(v).map(
            |(&xi, &vi)| { let d = f64::from(xi) - vi; sumsq += d*d; d }
        ).collect::<Vec<f64>>();
        dif.smultf64(1_f64/sumsq.sqrt())
    }

    /// Addition of `&[f64]` slice
    fn vaddf64(self, v:&[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)+vi).collect()
    }
  
    /// Euclidian distance to `&[f64]`  
    fn vdistf64(self, v:&[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (f64::from(xi)-vi).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Euclidian distance to `&[f64]` squared  
    fn vdistsqf64(self, v:&[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (f64::from(xi)-vi).powi(2))
            .sum::<f64>()
    }
}
