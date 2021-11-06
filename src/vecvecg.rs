use crate::{Stats,Vecg,MutVecf64,VecVecg,VecVec};
pub use indxvec::{merge::*,Indices};

impl<T,U> VecVecg<T,U> for &[Vec<T>] where T: Copy+PartialOrd, 
    f64: From<T>, U: Copy+PartialOrd,f64: From<U>  {

    /// Weighted Centre
    fn wacentroid(self,ws: &[U]) -> Vec<f64> where {
        let mut centre = vec![0_f64; self[0].len()];
        let mut wsum = 0_f64;
        self.iter().zip(ws).for_each(|(s,w)|
        {   let wf = f64::from(*w);
            wsum += wf;
            centre.mutvadd(&s.smult(wf))
        });
        centre.mutsmult(1.0 / wsum);
        centre
    }

    /// Trend computes the vector connecting the geometric medians of two sets of multidimensional points.
    /// This is a robust relationship between two unordered multidimensional sets.
    /// The two sets have to be in the same space but can have different numbers of points.
    fn trend(self, eps: f64, v: Vec<Vec<U>>) -> Vec<f64> {
        let mut m1 = self.gmedian(eps);       
        m1.mutvsub(&v.gmedian(eps));
        m1
    }

    /// Translates the whole set by subtracting vector m. Returns Vec of Vecs.
    /// When m is set to the geometric median, this produces the zero median form.
    /// The geometric median is invariant with respect to rotation,
    /// unlike the often used mean (`acentroid` here), or the quasi median,
    /// both of which depend on the choice of axis.
     fn translate(self, m: &[U]) -> Vec<Vec<f64>> {  
        self.iter().map(|point| point.vsub(m)).collect()   
    }

    /// Individual distances from any point v, typically not a member, to all the members of self.    
    fn dists(self, v: &[U]) -> Vec<f64> {
        self.iter().map(|p| p.vdist(v)).collect()
    }

    /// The sum of distances from any single point v, typically not a member, to all the members of self.    
    /// Geometric Median is defined as the point which minimises this function.
    fn distsum(self, v: &[U]) -> f64 {
        self.iter().map(|p| p.vdist(v)).sum::<f64>()
    }

    /// Sorted eccentricities magnitudes, w.r.t. weighted geometric median.
    /// associated cummulative probability density function in [0,1] of the weights.
    fn wsortedeccs(self, ws: &[U], eps:f64) -> ( Vec<f64>,Vec<f64>,Vec<f64> ) { 
        let mut eccs = Vec::with_capacity(self.len()); 
        let gm = self.wgmedian(ws,eps);
        for v in self { // collect true ecentricities magnitudes
            eccs.push(gm.vdist(&v)) 
        }
        // Apply linear transform
        // eccs = eccs.lintrans();
        // create sort index of the eccs
        let index = sortidx(&eccs);
        // pick the associated points weights in the order of the sorted eccs
        let mut weights = index.unindexf64(&ws,true);
        let mut sumw = 0_f64;
        // accummulate the weights
        weights.iter_mut().for_each(|w|{ sumw += *w; *w = sumw });     
        // divide by the final sum to get cummulative probabilities in [0,1]
        weights.iter_mut().for_each(|w| *w /= sumw );     
        ( gm, index.unindex(&eccs, true), weights )
    }

    /// Sorted cosines magnitudes,
    /// associated cummulative probability density function in [0,1] of the weights.
    /// Needs central median
    fn wsortedcos(self, medmed: &[U], zeromed: &[U], ws: &[U]) -> ( Vec<f64>,Vec<f64> ) { 
        let mut coses = Vec::with_capacity(self.len());  
        for p in self { // collect coses      
            coses.push(p.vsub(&medmed).vsim(&zeromed)); 
        }
        // Apply linear transform
        // coses = coses.lintrans();
        // create sort index of the coses
        let index = sortidx(&coses);
        // pick the associated points weights in the same order as the sorted coses
        let mut weights = index.unindexf64(&ws,true);
        let mut sumw = 0_f64;
        // accummulate the weights to from cpdf
        weights.iter_mut().for_each(|w|{ sumw += *w; *w = sumw });
        // divide by the sum to get cum. probabilities in [0,1]
        weights.iter_mut().for_each(|w| *w /= sumw );     
        ( index.unindex(&coses,true), weights )
    }

    /// Error vector for (usually non-member) point p, 
    /// i.e. unscaled eccentricity vector.
    /// The true geometric median would return zero vector.
    fn errorv(self, p:&[U]) -> Vec<f64> {
        let mut vsum = vec![0_f64; p.len()];
        for x in self {  vsum.mutvadd(&x.vsubunit(&p)) };
        vsum.mutsmult(1_f64 / (self.len() as f64));
        vsum
    }

    /// Next approximate weighted median, from a non member point. 
    fn wnxnonmember(self, ws:&[U], p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for i in 0..self.len() { 
            let dvmag = self[i].vdist(p);
            if !dvmag.is_normal() { continue } // zero distance, safe to ignore
            let rec = f64::from(ws[i])/dvmag; // ws[i] is weigth for this point self[i]
            vsum.mutvadd(&self[i].smult(rec)); // add weighted vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmult(1.0/recip);
        vsum
    } 

    /// Secant method with recovery from divergence
    /// for finding the weighted geometric median
    fn wgmedian(self, ws: &[U], eps: f64) -> Vec<f64> {  
        let mut p = self.wacentroid(ws);
        let mut mag1 = p.vmag();
        let mut pdif = mag1;  
        loop {  
            let mut np = self.wnxnonmember(ws,&p); 
            let e = np.vsub(&p); // new vector error, or eccentricity  
            // let e = self.errorv(&p);
            let mag2 = e.vmag(); 
            // if mag2 < eps  { return np }; 
            // overwrite np with a better secant estimate  
            np = if mag1 > mag2 {  // eccentricity magnitude decreased, good, employ secant
                p.vadd(&e.smult(pdif/(mag1-mag2)))                   
            }
            else { // recovery: probably overshot the minimum, shorten the jump 
                   // e will already be pointing moreless back
                p.vadd(&e.smult(pdif/(mag1+mag2)))                    
            };
            pdif = np.vdist(&p);
            if pdif < eps { return np };              
            mag1 = mag2; 
            p = np            
        }       
    }
    
    /// Flattened lower triangular part of a covariance matrix for f64 vectors in self.
    /// Since covariance matrix is symmetric (positive semi definite), 
    /// the upper triangular part can be trivially generated for all j>i by: c(j,i) = c(i,j).
    /// N.b. the indexing is always assumed to be in this order: row,column.
    /// The items of the resulting lower triangular array c[i][j] are here flattened
    /// into a single vector in this double loop order: left to right, top to bottom 
    fn covar(self, m:&[U]) -> Vec<f64> {
        let n = self[0].len(); // dimension of the vector(s)
        let mut cov:Vec<f64> = vec![0_f64; (n+1)*n/2]; // flat lower triangular results array  
        for thisp in self { // adding up covars for all the points
            let mut covsub = 0_usize; // subscript into the flattened array cov
            let vm = thisp.vsub(&m);  // zero mean vector
            for i in 0..n {
                let thisc = vm[i]; // ith component
                // its products up to and including the diagonal (itself)
                for j in 0..i+1 { 
                    cov[covsub] += thisc*vm[j];
                    covsub += 1
                }
            }
        }
        // now compute the means and return
        cov.mutsmult(1.0_f64/self.len()as f64);
        cov
    }

    /// Flattened lower triangular part of a comediance matrix for f64 vectors in self.
    /// Since comediance matrix is symmetric (positive semi definite), 
    /// the upper triangular part can be trivially generated for all j>i by: c(j,i) = c(i,j).
    /// N.b. the indexing is always assumed to be in this order: row,column.
    /// The items of the resulting lower triangular array c[i][j] are here flattened
    /// into a single vector in this double loop order: left to right, top to bottom. 
    /// Instead of averaging these vectors over n points, their median is returned.
    /// Warning: may run out of memory for large number of points and high dimensionality.
    fn comed(self, m:&[U], eps:f64) -> Vec<f64> { // m should be the median here
        let d = self[0].len(); // dimension of the vector(s)
        let mut coms:Vec<Vec<f64>> = Vec::with_capacity(self.len()); // vec of flat lower triangular results arrays  
        for thisp in self { // saving comeds for all the points 
            let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2);
            let vm = thisp.vsub(&m);  // zero mediance vector
            for i in 0..d {
                let thisc = vm[i]; // ith component
                // its products up to and including the diagonal (itself)
                for j in 0..i+1 { com.push(thisc*vm[j]) }
            }
            coms.push(com)
        } 
        coms.gmedian(eps) // return the median of the comeds
    }

    /// Flattened lower triangular part of a covariance matrix for weighted f64 vectors in self.
    /// Since covariance matrix is symmetric (positive semi definite), 
    /// the upper triangular part can be trivially generated for all j>i by: c(j,i) = c(i,j).
    /// N.b. the indexing is always assumed to be in this order: row,column.
    /// The items of the resulting lower triangular array c[i][j] are here flattened
    /// into a single vector in this double loop order: left to right, top to bottom 
    fn wcovar(self, ws:&[U], m:&[f64]) -> Vec<f64> {
        let n = self[0].len(); // dimension of the vector(s)
        // let mut covs:Vec<Vec<f64>> = Vec::new();
        let mut cov:Vec<f64> = vec![0_f64; (n+1)*n/2]; // flat lower triangular results array
        let mut wsum = 0_f64;
        for h in 0..self.len() { // adding up covars for all the points
            let w = f64::from(ws[h]);
            wsum += w;
            let mut covsub = 0_usize; // subscript into the flattened array cov 
            let vm = self[h].vsub(&m);  // subtract zero mean/median vector   
            for i in 0..n { // cross multiply the components of one point
                let thisc = vm[i]; // ith component of the zero mean vector               
                for j in 0..i+1 {  // its weighted products up to and including the diagonal 
                    cov[covsub] += w*thisc*vm[j];
                    covsub += 1
                }
            } 
        }
        // now compute the means and return
        cov.mutsmult(1_f64/wsum); 
        cov
    }   
    /// Flattened lower triangular part of a comediance matrix for weighted f64 vectors in self.
    /// Since comediance matrix is symmetric (positive semi definite), 
    /// the upper triangular part can be trivially generated for all j>i by: c(j,i) = c(i,j).
    /// N.b. the indexing is always assumed to be in this order: row,column.
    /// The items of the resulting lower triangular array c[i][j] are here flattened
    /// into a single vector in this double loop order: left to right, top to bottom.
    /// Instead of averaging these vectors over n points, their median is returned.
    /// Warning: may run out of memory for large number of points and high dimensionality. 
    fn wcomed(self, ws:&[U], m:&[f64], eps:f64) -> Vec<f64> {
        let d = self[0].len(); // dimension of the vector(s)
        let mut coms:Vec<Vec<f64>> = Vec::with_capacity(self.len());  
        for h in 0..self.len() { // saving comeds for all the points
            let w = f64::from(ws[h]);
            let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2);  
            let vm = self[h].vsub(&m);  // subtract zero median vector   
            for i in 0..d { // cross multiply the components of one point              
                for j in 0..i+1 {  // its weighted products up to and including the diagonal 
                    com.push(w*vm[i]*vm[j])  }
            }
            coms.push(com) 
        }
        coms.wgmedian(ws, eps) // return the median of the comeds  
    }   
}
