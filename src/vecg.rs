use crate::{Stats,Vecg};
pub use indxvec::{Indices,merge::{sortm,rank}};

impl<T,U> Vecg<T,U> for &[T] 
    where T: Copy+PartialOrd,
          U: Copy+PartialOrd,
          f64: From<T>, f64: From<U> {

    /// Scalar multiplication of a vector, creates new vec
    fn smult(self, s:U) -> Vec<f64> {
        let sf = f64::from(s);
        self.iter().map(|&x| sf*(f64::from(x))).collect()
     }
    /// scalar multiplication by f64
    fn smultf64(self, s:f64) -> Vec<f64> { 
        self.iter().map(|&x| s*(f64::from(x))).collect()
     }
     /// Scalar addition to a vector, creates new vec
     fn sadd(self, s:U) -> Vec<f64> {
         let sf = f64::from(s);
        self.iter().map(|&x| sf+(f64::from(x))).collect()
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
    /// Vector subtraction of `&[f64]`
    fn vsubf64(self, v:&[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)-vi).collect()
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
    /// Addition of `&[f64]` slice
    fn vaddf64(self, v:&[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| f64::from(xi)+vi).collect()
    }
    /// Euclidian distance   
    fn vdist(self, v:&[U]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (f64::from(xi)-f64::from(vi)).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    /// Euclidian distance to `&[f64]`  
    fn vdistf64(self, v:&[f64]) -> f64 {
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
        let n = self.len(); // dimension of the vector
        let mut cov:Vec<f64> = Vec::new(); // flat lower triangular result array
        let vm = self.vsub(&m); // zero mean vector
        for i in 0..n {
            let thisc = vm[i]; // take this component
            // generate its products up to and including the diagonal (itself)
            for j in 0..i+1 { cov.push(thisc*vm[j]) }
        }
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
 
    /// Joint entropy 
    fn jointentropy(self, v:&[U]) -> f64 {      
        let pdfs = sortm(self,true).pdf();
        let pdfv = sortm(v,true).pdf();
        let mut entr = 0_f64;
        for i in 0..pdfs.len() {
            for j in 0..pdfv.len() {
                let p = pdfs[i]*pdfv[j]; 
                entr -= p*(p.ln()) 
            }
        }                
        entr           
    }
    /// Dependence of &[T] &[U] variables.
    /// i.e. `dependence` returns 0 iff they are statistically independent
    /// and 1 when they are identical
    fn dependence(self, v:&[U]) -> f64 {   
        2.0 - (self.entropy() + v.entropy())/self.jointentropy(v) 
    }         

    /// Pearson's correlation. 
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
        let xvec = rank(&self,true);
        let yvec = rank(&v,true); // rank from crate idxvec::merge
        // It is just Pearson's correlation of ranks
        xvec.ucorrelation(&yvec) // using Indices trait from idxvec
    }
}
