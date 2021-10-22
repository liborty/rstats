use std::iter::FromIterator;

use crate::{Med, MStats, MinMax, MutVecg, MutVecf64, Stats, VecVec, Vecg};
pub use indxvec::{merge::*,Indices};

impl<T> VecVec<T> for &[Vec<T>] where T: Copy+PartialOrd, 
    f64: From<T>, T: From<f64> {
    /// acentroid = simple multidimensional arithmetic mean
    /// # Example
    /// ```
    /// use rstats::{Vecg,VecVec,functions::genvec};
    /// let pts = genvec(15,15,255,30); 
    /// let dist = pts.distsum(&pts.acentroid());
    /// assert_eq!(dist, 4.14556218326653_f64);
    /// ```
    fn acentroid(self) -> Vec<f64> {
        let mut centre = vec![0_f64; self[0].len()];
        for v in self {
            centre.mvadd(&v)
        }
        centre.mutsmult(1.0 / self.len() as f64);
        centre
    }
    /// Weighted Centre
    fn wacentroid(self,ws: &[T]) -> Vec<f64> where {
        let mut centre = vec![0_f64; self[0].len()];
        let mut wsum = 0_f64;
        for i in 0..self.len() {
            let w = ws[i];
            wsum += f64::from(w);
            centre.mutvadd(&self[i].smult(w))
        }
        centre.mutsmult(1.0 / (wsum*self.len() as f64));
        centre
    }

    /// gcentroid = multidimensional geometric mean
    /// # Example
    /// ```
    /// use rstats::{Vecg,VecVec,functions::genvec};
    /// let pts = genvec(15,15,255,30);
    /// let centre = pts.gcentroid();
    /// let dist = pts.distsum(&centre);
    /// assert_eq!(dist,4.897594485332543_f64);
    /// ```
    fn gcentroid(self) -> Vec<f64> {
        let n = self.len() as f64; // number of points
        let d = self[0].len(); // dimensions
        let mut centre = vec![0_f64; d];
        let mut lnvec = vec![0_f64; d];
        for v in self { 
            for i in 0..d { lnvec[i] = f64::from(v[i]).ln() }
            centre.mutvadd(&lnvec)
        }
        centre.iter().map(|comp| (comp/n).exp()).collect()        
    }

    /// hcentroid =  multidimensional harmonic mean
    /// # Example
    /// ```
    /// use rstats::{Vecg,VecVec,functions::genvec};
    /// let pts = genvec(15,15,255,30);
    /// let centre = pts.hcentroid();
    /// let dist = pts.distsum(&centre);
    /// assert_eq!(dist, 5.623778191797538_f64);
    /// ```
    fn hcentroid(self) -> Vec<f64> {
        let mut centre = vec![0_f64; self[0].len()];
        // let t = self.translate(&self.acentroid());
        for v in self {
            centre.mutvadd(&v.vinverse())
        }
        centre.vinverse()       
    }

    /// Trend computes the vector connecting the geometric medians of two sets of multidimensional points.
    /// This is a robust relationship between two unordered multidimensional sets.
    /// The two sets have to be in the same space but can have different numbers of points.
    fn trend(self, eps: f64, v: Vec<Vec<T>>) -> Vec<f64> {
        let mut m1 = self.gmedian(eps);       
        m1.mutvsub(&v.gmedian(eps));
        m1
    }

    /// Translates the whole set by vector -m. Returns Vec of Vecs.
    /// When m is set to the geometric median, this produces the zero median form.
    /// The geometric median is invariant with respect to rotation,
    /// unlike the often misguidedly used mean (`acentroid` here), or the quasi median,
    /// both of which depend on the choice of axis.
     fn translate(self, m: &[T]) -> Vec<Vec<f64>> {  
        self.iter().map(|point| point.vsub(m)).collect()   
    }

    /// For each member point, gives its sum of distances to all other points.
    fn distsums(self) -> Vec<f64> {
        let n = self.len();
        let mut dists = vec![0_f64; n]; // distances accumulator for all points
                                        // examine all unique pairings (lower triangular part of symmetric flat matrix)
        for i in 1..n {
            let thisp = &self[i];
            for j in 0..i {
                let thatp = &self[j];
                let d = thisp.vdist(&thatp); // calculate each distance relation just once
                dists[i] += d;
                dists[j] += d; // but add it to both points
            }
        }
        dists
    }

    /// The sum of distances from one member point, given by its `indx`, to all the other points in self.
    /// For all the points, use more efficient `distsums`.    
    fn distsuminset(self, indx: usize) -> f64 {
        let n = self.len();
        let mut sum = 0_f64;
        let thisp = &self[indx];
        for i in 0..n {
            if i == indx {
                continue;
            };
            sum += self[i].vdist(&thisp)
        }
        sum
    }

    /// Individual distances from any point v, typically not a member, to all the members of self.    
    fn dists(self, v: &[T]) -> Vec<f64> {
        self.iter().map(|p| p.vdist(v)).collect()
    }

    /// The sum of distances from any single point v, typically not a member, to all the members of self.    
    /// Geometric Median is defined as the point which minimises this function.
    fn distsum(self, v: &[T]) -> f64 {
        self.iter().map(|p| p.vdist(v)).sum::<f64>()
    }

    /// Medoid is the member point (point belonging to the set of points `self`), 
    /// which has the least sum of distances to all other points.
    /// Outlier is the point with the greatest sum of distances.
    /// In other words, they are the members nearest and furthest from the median.
    /// This function returns a four-tuple:  
    /// (medoid_distance, medoid_index, outlier_distance, outlier_index).
    /// # Example
    /// ```
    /// use rstats::{Vecg,VecVec,functions::genvec};
    /// let pts = genvec(15,15,255,30);
    /// let mm = pts.medoid();
    /// assert_eq!(mm.min,4.812334638782327_f64);
    /// ```
    fn medoid(self) -> MinMax<T> where Vec<T>:FromIterator<f64> {    
        let dists:Vec<T> = self.distsums().iter().map(|&d| d).collect();
        minmax(&dists)
    }

    /// Finds approximate vectors from each member point towards the geometric median.
    /// Twice as fast as doing them individually, using symmetry.
    fn eccentricities(self) -> Vec<Vec<f64>> {
        let n = self.len();
        // allocate vectors for the results
        let mut eccs = vec![vec![0_f64; self[0].len()]; n];
        let mut recips = vec![0_f64; n];
        // ecentricities vectors accumulator for all points
        // examine all unique pairings (lower triangular part of symmetric flat matrix)
        for i in 1..n {
            let thisp = &self[i];
            for j in 0..i { 
                // calculate each unit vector between any pair of points just once
                let dvmag = self[j].vdist(&thisp);             
                if !dvmag.is_normal() { continue }
                let rec = 1.0/dvmag;
                eccs[i].mutvadd(&self[j].smultf64(rec));
                recips[i] += rec;
                // mind the vector's opposite orientations w.r.t. to the two points!
                eccs[j].mutvsub(&self[j].smultf64(rec)); 
                recips[j] += rec;
            }
        }
        for i in 0..n { 
            eccs[i].mutsmult(1.0/recips[i]); 
            eccs[i].mvsub(&self[i]) 
        }
        eccs
    }

    /// Exact eccentricity vectors from all member points by first finding the Geometric Median.
    /// Usually faster than the approximate `eccentricities` above, especially when there are many points.
    fn exacteccs(self, eps:f64) -> Vec<Vec<f64>> {        
        let mut eccs = Vec::with_capacity(self.len()); // Vectors for the results
        let gm = self.gmedian(eps);
        for v in self {
            eccs.push(gm.vsub(&v))
        }
        eccs
    }
 
    /// GM and sorted eccentricities magnitudes.
    /// Describing a set of points `self` in n dimensions
    fn sortedeccs(self, ascending:bool, eps:f64) -> ( Vec<f64>,Vec<f64> ) { 
        let mut eccs = Vec::with_capacity(self.len());       
        let gm = self.gmedian(eps);
        for v in self { // collect raw ecentricities magnitudes
            eccs.push(gm.vdist(&v)) 
        }
        ( gm, sortm(&eccs,ascending) )
    }
    
    /// Weighted geometric median, sorted eccentricities magnitudes,
    /// associated cummulative probability density function in [0,1] of the weights.
    fn wsortedeccs(self, ws: &[T], eps:f64) -> ( Vec<f64>,Vec<f64>,Vec<f64> ) { 
        let mut eccs = Vec::with_capacity(self.len()); 
        let gm = self.wgmedian(ws,eps);
        for v in self { // collect true ecentricities magnitudes
            eccs.push(gm.vdist(&v)) 
        }
        // Apply linear transform
        // eccs = eccs.lintrans();
        // create sort index of the eccs
        let index = sortidx(&eccs);
        // pick the associated points weights in the reverse order of the sorted eccs
        let mut weights = index.unindexf64(&ws,true);
        let mut sumw = 0_f64;
        // accummulate the weights 
        for i in 0..weights.len() {
            sumw += weights[i]; 
            weights[i] = sumw
        }
        // divide by the sum to get cum. probabilities in [0,1]
        weights.iter_mut().for_each(|w| *w /= sumw );     
        ( gm, index.unindex(&eccs, true), weights )
    }

    /// Sorted cosines magnitudes,
    /// associated cummulative probability density function in [0,1] of the weights.
    /// Needs central median
    fn wsortedcos(self, medmed: &[T], zeromed: &[T], ws: &[T]) -> ( Vec<f64>,Vec<f64> ) { 
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
        for i in 0..weights.len() {
            sumw += weights[i]; 
            weights[i] = sumw
        }
        // divide by the sum to get cum. probabilities in [0,1]
        weights.iter_mut().for_each(|w| *w /= sumw );     
        ( index.unindex(&coses,true), weights )
    }

    /// Vector `eccentricity` or measure of  
    /// `not being the geometric median`, added to a given member point.
    /// The true geometric median is as yet unknown.
    /// The member point in question is specified by its index `indx`.
    /// This function is suitable for a single member point. 
    fn nxmember(self, indx: usize) -> Vec<f64> {
        let n = self.len();
        let mut vsum = vec![0_f64; self[0].len()];
        let thisp = &self[indx];
        let mut recip = 0_f64;
        for i in 0..n {
            if i == indx { continue  }; // exclude this point
            let dvmag = self[i].vdist(&thisp);
            if !dvmag.is_normal() { continue } // too close to this one
            let rec = 1.0/dvmag;
            vsum.mutvadd(&self[i].smultf64(rec)); // add vector
            recip += rec // add separately the reciprocals
        }
        vsum.mutsmult(1.0/recip);
        vsum 
    }    

    /// Vector `eccentricity` or measure of  
    /// `not being the geometric median` for a member point.
    /// The true geometric median is as yet unknown.
    /// The true geometric median would return zero.
    /// The member point in question is specified by its index `indx`.
    /// This function is suitable for a single member point. 
    /// When eccentricities of all the points are needed, use `exacteccs` above.
    fn eccmember(self, indx: usize) -> Vec<f64> {
        self.nxmember(indx).vsub(&self[indx])
    }

    /// Eccentricity vector added to a non member point,
    /// while the true geometric median is as yet unknown. 
    /// This function is suitable for a single non-member point. 
    fn nxnonmember(self, p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for x in self { 
            let dvmag = x.vdistf64(p);
            if !dvmag.is_normal() { continue } // zero distance, safe to ignore
            let rec = 1.0/dvmag;
            vsum.mutvadd(&x.smultf64(rec)); // add vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmult(1.0/recip);
        vsum
    }

    /// Error vector for (usually non-member) point p, 
    /// i.e. unscaled eccentricity vector.
    /// The true geometric median would return zero vector.
    fn errorv(self, p:&[T]) -> Vec<f64> {
        let mut vsum = vec![0_f64; p.len()];
        for x in self {  vsum.mutvadd(&x.vsubunit(&p)) };
        vsum.mutsmult(1_f64 / (self.len() as f64));
        vsum
    }

    /// Eccentricity vector for a non member point,
    /// while the true geometric median is as yet unknown.
    /// Returns the eccentricity vector.
    /// The true geometric median would return zero vector.
    /// This function is suitable for a single non-member point. 
    fn eccnonmember(self, p:&[f64]) -> Vec<f64> {
        self.nxnonmember(p).vsubf64(p)
    }

    /// Next approximate weighted median, from a non member point. 
    fn wnxnonmember(self, ws:&[T], p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for i in 0..self.len() { 
            let dvmag = self[i].vdistf64(p);
            if !dvmag.is_normal() { continue } // zero distance, safe to ignore
            let rec = f64::from(ws[i])/dvmag; // ws[i] is weigth for this point self[i]
            vsum.mutvadd(&self[i].smultf64(rec)); // add weighted vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmult(1.0/recip);
        vsum
    } 

    /*  Magnitudes of all the vectors in self
    fn mags(self) -> Vec<f64> {
        let mut magsv = Vec::new();
        for v in self {
            magsv.push(v.vmag())
        }
        magsv
    } 
    */

    /// Mean and Std (in MStats struct), Median and quartiles (in Med struct), Median and Outlier (in MinMax struct) 
    /// of scalar eccentricities of points in self.
    /// These are new robust measures of a cloud of multidimensional points (or multivariate sample).  
    fn eccinfo(self, eps: f64) -> (MStats, Med, MinMax<T>) where Vec<T>:FromIterator<f64> {
        let gm = self.gmedian(eps);
        let eccs:Vec<T> = self.iter().map(|v| gm.vdist(&v)).collect();
        (eccs.ameanstd().unwrap(),eccs.median().unwrap(),minmax(&eccs))
    }

    /// Eccentricity defined Medoid and Outlier.
    /// This can give different results to `medoid` above, defined by sums of distances,
    /// especially for the outliers. See tests.rs.  
    /// Consider some point c and some other points, bunched up at a distance r from c.
    /// The sum of their distances will be nr. Now, spread those points around a circle of radius r from c.
    /// The sum of their distances from c will remain the same but the eccentricity of c will be much reduced.
    /// # Example
    /// ```
    /// use rstats::{Vecg,VecVec,functions::genvec};
    /// pub const EPS:f64 = 1e-7;
    /// let d = 6_usize;
    /// let pt = genvec(d,24,7,13); // random test data 5x20
    /// let mm = pt.emedoid(EPS);
    /// assert_eq!(mm.minindex,10); // index of e-medoid
    /// assert_eq!(mm.maxindex,20);  // index of e-outlier
    /// ```
    fn emedoid(self, eps: f64) -> MinMax<T> where Vec<T>:FromIterator<f64> {
    let gm = self.gmedian(eps);
    let eccs:Vec<T> = self.iter().map(|v| gm.vdist(&v)).collect();
    minmax(&eccs)
}

    /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    /// It has (provably) only vector iterative solutions.
    /// Search methods are slow and difficult in highly dimensional space.
    /// Weiszfeld's fixed point iteration formula had known problems with sometimes failing to converge.
    /// Especially, when the points are dense in the close proximity of the gm, or gm coincides with one of them.  
    /// However, these problems are fixed in my new algorithm here.      
    /// There will eventually be a multithreaded version.
    fn nmedian(self, eps: f64) -> Vec<f64> {
        let mut point = self.acentroid(); // start iterating from the Centre
        loop {
            let nextp = self.nxnonmember(&point);         
            if nextp.vdistf64(&point) < eps { return nextp }; // termination
            point = nextp
        } 
    }
 
    /// Possible first iteration point for geometric medians.
    /// Same as eccnonmember(origin), just saving the zero subtractions. 
    /// # Example
    /// ```
    /// use rstats::{Vecg,VecVec,functions::genvec};
    /// let pts = genvec(15,15,255,30); 
    /// let dist = pts.distsum(&pts.firstpoint());
    /// assert_eq!(dist,4.132376831171272_f64);
    /// ```
    fn firstpoint(self) -> Vec<f64> {
        let mut rsum = 0_f64;
        let mut vsum = vec![0_f64; self[0].len()];
        for thisp in self {
            let mag = thisp.vmag();
            if mag.is_normal() {  
                let invmag = 1.0_f64/mag;
                rsum += invmag;
                vsum.mutvadd(&thisp.smultf64(invmag)) // accumulate unit vectors
            }
        }
        vsum.smultf64(1.0/rsum) // scale by the sum of reciprocals
    }

    /// Secant method with recovery from divergence
    /// for finding the geometric median
    fn gmedian(self, eps: f64) -> Vec<f64> {  
        let mut p = self.acentroid();
        let mut mag1 = p.vmag();
        let mut pdif = mag1;  
        loop {  
            let mut np = self.nxnonmember(&p); 
            let e = np.vsubf64(&p); // new vector error, or eccentricity  
            // let e = self.errorv(&p);
            let mag2 = e.vmag(); 
            // if mag2 < eps  { return np }; 
            // overwrite np with a better secant estimate  
            np = if mag1 > mag2 {  // eccentricity magnitude decreased, good, employ secant
                e.smultf64(pdif/(mag1-mag2)).vaddf64(&p)
            }
            else { // recovery: probably overshot the minimum, shorten the jump 
                   // e will already be pointing moreless back
                e.smultf64(pdif/(mag1+mag2)).vaddf64(&p)                    
            };
            pdif = np.vdistf64(&p);
            if pdif < eps { return np };              
            mag1 = mag2; 
            p = np            
        }       
    }

    /// Secant method with recovery from divergence
    /// for finding the weighted geometric median
    fn wgmedian(self, ws: &[T], eps: f64) -> Vec<f64> {  
        let mut p = self.wacentroid(ws);
        let mut mag1 = p.vmag();
        let mut pdif = mag1;  
        loop {  
            let mut np = self.wnxnonmember(ws,&p); 
            let e = np.vsubf64(&p); // new vector error, or eccentricity  
            // let e = self.errorv(&p);
            let mag2 = e.vmag(); 
            // if mag2 < eps  { return np }; 
            // overwrite np with a better secant estimate  
            np = if mag1 > mag2 {  // eccentricity magnitude decreased, good, employ secant
                p.vaddf64(&e.smultf64(pdif/(mag1-mag2)))                   
            }
            else { // recovery: probably overshot the minimum, shorten the jump 
                   // e will already be pointing moreless back
                p.vaddf64(&e.smultf64(pdif/(mag1+mag2)))                    
            };
            pdif = np.vdistf64(&p);
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
    fn covar(self, m:&[f64]) -> Vec<f64> {
        let n = self[0].len(); // dimension of the vector(s)
        let mut cov:Vec<f64> = vec![0_f64; (n+1)*n/2]; // flat lower triangular results array  
        for thisp in self { // adding up covars for all the points
            let mut covsub = 0_usize; // subscript into the flattened array cov
            let vm = thisp.vsubf64(&m);  // zero mean vector
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
    fn comed(self, m:&[f64], eps:f64) -> Vec<f64> { // m should be the median here
        let d = self[0].len(); // dimension of the vector(s)
        let mut coms:Vec<Vec<f64>> = Vec::with_capacity(self.len()); // vec of flat lower triangular results arrays  
        for thisp in self { // saving comeds for all the points 
            let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2);
            let vm = thisp.vsubf64(&m);  // zero mediance vector
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
    fn wcovar(self, ws:&[f64], m:&[f64]) -> Vec<f64> {
        let n = self[0].len(); // dimension of the vector(s)
        // let mut covs:Vec<Vec<f64>> = Vec::new();
        let mut cov:Vec<f64> = vec![0_f64; (n+1)*n/2]; // flat lower triangular results array
        let mut wsum = 0_f64;
        for h in 0..self.len() { // adding up covars for all the points
            let w = ws[h];
            wsum += w;
            let mut covsub = 0_usize; // subscript into the flattened array cov 
            let vm = self[h].vsubf64(&m);  // subtract zero mean/median vector   
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
    fn wcomed(self, ws:&[f64], m:&[f64], eps:f64) -> Vec<f64> {
        let d = self[0].len(); // dimension of the vector(s)
        let mut coms:Vec<Vec<f64>> = Vec::with_capacity(self.len());  
        for h in 0..self.len() { // saving comeds for all the points
            let w = ws[h];
            let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2);  
            let vm = self[h].vsubf64(&m);  // subtract zero median vector   
            for i in 0..d { // cross multiply the components of one point              
                for j in 0..i+1 {  // its weighted products up to and including the diagonal 
                    com.push(w*vm[i]*vm[j])  }
            }
            coms.push(com) 
        }
        coms.wgmedian(ws, eps) // return the median of the comeds  
    }   
}
