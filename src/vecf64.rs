use crate::{Vecf64};
use indxvec::{merge::*,Indices};

impl Vecf64 for &[f64] {

    /// Scalar multiplication of a vector, creates new vec
    fn smult(self, s:f64) -> Vec<f64> {
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

    fn vinverse(self) -> Vec<f64> {
        self.smult(1.0/self.vmagsq())
    }

    // negated vector (components with opposite sign)
    fn negv(self) -> Vec<f64> { self.smult(-1.0) }

    /// Cosine of an angle between two vectors.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.cosine(&v2),0.7517241379310344);
    /// ```
    fn cosine(self, v: &[f64]) -> f64 {
        let (mut sxy, mut sy2) = (0_f64, 0_f64);
        let sx2: f64 = self
            .iter()
            .zip(v)
            .map(|(&x, &y)| {
                sxy += x * y;
                sy2 += y * y;
                x*x
            })
            .sum();
        sxy / (sx2*sy2).sqrt()
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
    /// Euclidian distance squared between two n dimensional points (vectors).  
    /// Slightly faster than vsub followed by vmasq, as both are done in one loop
    /// Same as vdist without taking the square root
    fn vdistsq(self, v: &[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (xi - vi).powi(2))
            .sum::<f64>() 
    } 

    /// cityblock distance
    fn cityblockd(self, v:&[f64]) -> f64 {
        self.iter()
        .zip(v)
        .map(|(&xi, &vi)| (xi-vi).abs()) 
        .sum::<f64>()      
    }

    /// Vector magnitude
    fn vmag(self) -> f64 {
        self.iter().map(|&x| x.powi(2)).sum::<f64>().sqrt()
    }

    /// Vector magnitude squared
    fn vmagsq(self) -> f64 {
        self.iter().map(|&x| x.powi(2)).sum::<f64>()
    }

    /// Unit vector
    fn vunit(self) -> Vec<f64> {
        self.smult(1. / self.iter().map(|x| x.powi(2)).sum::<f64>().sqrt())
    }

    /// Unit vectors difference (done together for efficiency)
    fn vsubunit(self, v: &[f64]) -> Vec<f64> {
        let mut sumsq = 0_f64;
        let dif = self.iter().zip(v).map(
            |(&xi, &vi)| { let d = xi - vi; sumsq += d*d; d }
        ).collect::<Vec<f64>>();
        dif.smult(1_f64/sumsq.sqrt())
    }    

    /// Area of a parallelogram between two vectors.
    /// Same as the magnitude of their cross product |a ^ b| = |a||b|sin(theta).
    /// Attains maximum `|a|.|b|` when the vectors are othogonal.
    fn varea(self, v:&[f64]) -> f64 {
        (self.vmagsq()*v.vmagsq() - self.dotp(v).powi(2)).sqrt()
    }

    /// We define vector similarity S in the interval [0,1] as
    /// S = (1+cos(theta))/2
    fn vsim(self, v:&[f64]) -> f64 { (1.0+self.cosine(&v))/2.0 }

    /// We define vector dissimilarity D in the interval [0,1] as
    /// D = 1-S = (1-cos(theta))/2
    fn vdisim(self, v:&[f64]) -> f64 { (1.0-self.cosine(&v))/2.0 }

    /// Area proportional to the swept arc up to angle theta. 
    /// Attains maximum of `2|a||b|` when the vectors have opposite orientations.
    /// This is really |a||b|(1-cos(theta)) = 2|a||b|D
    fn varc(self, v:&[f64]) -> f64 { 
        (self.vmagsq()*v.vmagsq()).sqrt() - self.dotp(v)
   }
 
    /// Pearson's correlation coefficient of a sample of two f64 variables.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.correlation(&v2),-0.1076923076923077);
    /// ```
    fn correlation(self, v: &[f64]) -> f64 {
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
        let nf = self.len() as f64;
        (sxy - sx / nf * sy) / ((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()
    }

    /// Kendall Tau-B correlation coefficient of a sample of two f64 variables.
    /// Defined by: tau = (conc - disc) / sqrt((conc + disc + tiesx) * (conc + disc + tiesy))
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.kendalcorr(&v2),-0.07692307692307693);
    /// ```
    fn kendalcorr(self, v: &[f64]) -> f64 {
        let (mut conc, mut disc, mut tiesx, mut tiesy) = (0_i64, 0_i64, 0_i64, 0_i64);
        for i in 1..self.len() {
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
        (conc - disc) as f64 / (((conc + disc + tiesx) * (conc + disc + tiesy)) as f64).sqrt()
    }
    /// Spearman rho correlation coefficient of two f64 variables.
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.spearmancorr(&v2),-0.1076923076923077);
    /// ```
    fn spearmancorr(self, v: &[f64]) -> f64 {
        let xvec = rank(&self,true);
        let yvec = rank(&v,true);
        // It is just Pearson's correlation of ranks
        xvec.ucorrelation(&yvec)
    }

    /// (Auto)correlation coefficient of pairs of successive values of (time series) f64 variable.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.autocorr(),0.9984603532054123_f64);
    /// ```
    fn autocorr(self) -> f64 {
        let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let sx:f64 = self.windows(2).map(|w| {
            let x = w[0];
            let y = w[1];
            sy += y;
            sxy += x * y;
            sx2 += x * x;
            sy2 += y * y;
            x }).sum();
        let nf = self.len() as f64;
        (sxy - sx / nf * sy) / ((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()
    }
 
    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Vec<f64> {
        let (min,_,max,_) = minmax(self);
        let range = max-min;
        self.iter().map(|&x|(x-min)/range).collect()        
    }
/*
    /// New sorted vector. Immutable sort.
    /// Copies self and then sorts it in place, leaving self unchanged.
    /// Calls mutsortf and that calls the standard self.sort_unstable_by.
    /// Consider using our `sortm` instead.
    fn sortf(self) -> Vec<f64> {
        let mut sorted:Vec<f64> = self.to_vec();
        sorted.mutsortf();
        sorted      
    }
*/
    /// Flattened lower triangular part of a covariance matrix for a single f64 vector.
    /// Since covariance matrix is symmetric (positive semi definite), 
    /// the upper triangular part can be trivially added for all j>i by: c(j,i) = c(i,j).
    /// N.b. the indexing is always assumed to be in this order: row,column.
    /// The items of the resulting lower triangular array c[i][j] are pushed flat
    /// into a single vector in this double loop order: left to right, top to bottom 
    fn covone(self, m:&[f64]) -> Vec<f64> {
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

    /// Reconstructs the full symmetric square matrix from its lower diagonal compact form,
    /// as produced by covar, covone, wcovar
    fn symmatrix(self) -> Vec<Vec<f64>> {
        // solve quadratic equation to find the dimension of the square matrix
        let n = (((8*self.len()+1) as f64).sqrt() as usize - 1)/2;
        let mut mat = vec![vec![0_f64;n];n]; // create the square matrix 
        let mut selfindex = 0;
        for row in 0..n {     
            for column in 0..row { // excludes the diagonal  
                mat[row][column] = self[selfindex]; // just copy the value into the lower triangle
                mat[column][row] = self[selfindex]; // and into transposed upper position 
                selfindex += 1 // move to the next input value
            } // this row of lower triangle finished
            mat[row][row] = self[selfindex];  // now the diagonal element, no transpose
            selfindex += 1 // move to the next input value
        }
        mat
    }
    
}
