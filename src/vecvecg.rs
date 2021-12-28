use crate::{Vecg,Vecf64,MutVecf64,VecVecg,VecVec};
pub use indxvec::{merge::*,Indices};

impl<T,U> VecVecg<T,U> for &[Vec<T>] where T: Copy+PartialOrd+std::fmt::Display, 
    f64: From<T>, U: Copy+PartialOrd+std::fmt::Display,f64: From<U>  {

    /// Weighted Centre
    fn wacentroid(self,ws: &[U]) -> Vec<f64> where {
        let mut centre = vec![0_f64; self[0].len()];
        let mut wsum = 0_f64;
        self.iter().zip(ws).for_each(|(s,&w)|
        { 
            wsum += f64::from(w);
            centre.mutvaddf64(&s.smult(w))
        });
        centre.mutsmultf64(1.0 / wsum);
        centre
    }

    /// Trend computes the vector connecting the geometric medians of two sets of multidimensional points.
    /// This is a robust relationship between two unordered multidimensional sets.
    /// The two sets have to be in the same space but can have different numbers of points.
    fn trend(self, eps: f64, v: Vec<Vec<U>>) -> Vec<f64> {
        let mut m1 = self.gmedian(eps);       
        m1.mutvsubf64(&v.gmedian(eps));
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
    fn wsortedeccs(self, ws: &[U], gm: &[f64]) -> ( Vec<f64>,Vec<f64> ) { 
        let mut eccs = Vec::with_capacity(self.len()); 
        // collect true eccentricities magnitudes
        for v in self { eccs.push(v.vdistf64(&gm)) }
        // create sort index of the eccs
        let index = sortidx(&eccs);
        // pick the associated points weights in the order of the sorted eccs
        let mut weights = index.unindexf64(&ws,true);
        let mut sumw = 0_f64;
        // accummulate the weights
        weights.iter_mut().for_each(|w|{ sumw += *w; *w = sumw });     
        // divide by the final sum to get cummulative probabilities in [0,1]
        weights.iter_mut().for_each(|w| *w /= sumw );     
        ( index.unindex(&eccs, true), weights )
    }

    /// Sorted cosines magnitudes,
    /// associated cummulative probability density function in [0,1] of the weights.
    /// Needs central median
    fn wsortedcos(self, medmed: &[U], unitzmed: &[U], ws: &[U]) -> ( Vec<f64>,Vec<f64> ) { 
        let mut coses = Vec::with_capacity(self.len());  
        for p in self { // collect coses      
            coses.push(p.vsubunit(&medmed).dotp(&unitzmed)); 
        } 
        // create sort index of the coses
        let index = sortidx(&coses);
        // pick the associated points weights in the same order as the sorted coses
        let mut weights = index.unindexf64(&ws,true);
        let mut sumw = 0_f64;
        // accummulate the weights to form cpdf
        weights.iter_mut().for_each(|w|{ sumw += *w; *w = sumw });
        // divide by the sum to get cum. probabilities in [0,1]
        weights.iter_mut().for_each(|w| *w /= sumw );     
        ( index.unindex(&coses,true), weights )
    }

    /// Next approximate weighted median, from a non member point. 
    fn wnxnonmember(self, ws:&[U], p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for i in 0..self.len() { 
            let magsq = self[i].vdistsqf64(p);
            if !magsq.is_normal() { continue } // zero distance, safe to ignore
            let rec = f64::from(ws[i])/(magsq).sqrt(); // ws[i] is weigth for this point self[i]
            vsum.mutvaddf64(&self[i].smult(rec)); // add weighted vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmultf64(1_f64/recip);
        vsum
    }
    
    /// Estimated (computed) eccentricity vector for a non member point. 
    /// The true geometric median is as yet unknown.
    /// Returns the weighted eccentricity vector.
    /// The true geometric median would return zero vector.
    /// This function is suitable for a single non-member point. 
    fn weccnonmember(self, ws:&[U], p:&[f64]) -> Vec<f64> {
       self.wnxnonmember(ws,p).vsubf64(p)
    }

    /// Secant method with recovery from divergence
    /// for finding the weighted geometric median
    fn wgmedian(self, ws: &[U], eps: f64) -> Vec<f64> {  
        let eps2 = eps.powi(2);
        let mut point = self.wacentroid(ws); // start iterating from the Centre
        loop { // vector iteration till accuracy eps is reached
            let nextp = self.wnxnonmember(ws,&point);          
            if nextp.vdistsqf64(&point) < eps2 { return nextp }; // termination
            point = nextp
        } 
    }

    /// Secant recovery from divergence
    /// for finding the weighted geometric median.
    fn wsmedian(self, ws: &[U], eps: f64) -> Vec<f64> {  
        let eps2 = eps.powi(2);
        let mut p1 = self.wacentroid(ws);
        let mut mag1:f64 = p1.iter().map(|c| c.powi(2)).sum();
        loop {  
            let p2 = self.wnxnonmember(ws,&p1);                                       
            let e = p2.vsubf64(&p1); // new vector error, or eccentricity   
            let mag2:f64 = e.iter().map(|c| c.powi(2)).sum();
            if mag2 < eps2 { return p2 }; // termination
            p1.mutvaddf64(&e.smultf64(mag1/(mag1-mag2)));
            mag1 = mag2;
        };     
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
        cov.mutsmultf64(1.0_f64/self.len()as f64);
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
        cov.mutsmultf64(1_f64/wsum); 
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
