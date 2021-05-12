use crate::{Med, MStats, MutVectors, Stats, VecVec, Vecf64};

impl VecVec for &[Vec<f64>] {
    /// acentroid = simple multidimensional arithmetic mean
    /// # Example
    /// ```
    /// use rstats::{Vecf64,VecVec,functions::genvec};
    /// let pts = genvec(15,15,255,30);
    /// let centre = pts.acentroid();
    /// let dist = pts.distsum(&centre);
    /// assert_eq!(dist, 4.14556218326653_f64);
    /// ```
    fn acentroid(self) -> Vec<f64> {
        let mut centre = vec![0_f64; self[0].len()];
        for v in self {
            centre.mutvadd(&v)
        }
        centre.mutsmult(1.0 / self.len() as f64);
        centre
    }

    /// hcentroid =  multidimensional harmonic mean
    /// # Example
    /// ```
    /// use rstats::{Vecf64,VecVec,functions::genvec};
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
    fn trend(self, eps: f64, v: Vec<Vec<f64>>) -> Vec<f64> {
        let m1 = self.gmedian(eps);
        let m2 = v.gmedian(eps);
        m2.vsub(&m1)
    }

    /// Translates the whole set by vector -m. Returns Vec of Vecs.
    /// When m is set to the geometric median, this produces the zero median form.
    /// The geometric median is invariant with respect to rotation,
    /// unlike the often misguidedly used mean (`acentroid` here), or the quasi median,
    /// both of which depend on the choice of axis.
     fn translate(self, m: &[f64]) -> Vec<Vec<f64>> {
        let mut result = Vec::new();
        for point in self {
            result.push(point.vsub(m))
        }
        result
    }

    /// For each member point, gives its sum of distances to all other points.
    /// This is the efficient workhorse of distance based analysis.
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
    fn dists(self, v: &[f64]) -> Vec<f64> {
        self.iter().map(|p| p.vdist(v)).collect()
    }

    /// The sum of distances from any single point v, typically not a member, to all the members of self.    
    /// Geometric Median is defined as the point which minimises this function.
    fn distsum(self, v: &[f64]) -> f64 {
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
    /// use rstats::{Vecf64,VecVec,functions::genvec};
    /// let pts = genvec(15,15,255,30);
    /// let (dm,_,_,_) = pts.medoid();
    /// assert_eq!(dm,4.812334638782327_f64);
    /// ```
    fn medoid(self) -> (f64, usize, f64, usize) {
        self.distsums().minmax()
    }

    /// Finds approximate vectors from each member point towards the geometric median.
    /// Twice as fast as doing them individually, using symmetry.
    fn eccentricities(self) -> Vec<Vec<f64>> {
        let n = self.len();
        // allocate vectors for the results
        let mut eccs = vec![vec![0_f64; self[0].len()]; n]; // can be reduced to half
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
                eccs[i].mutvadd(&self[j].smult(rec));
                recips[i] += rec;
                // mind the vector's opposite orientations w.r.t. to the two points!
                eccs[j].mutvsub(&self[j].smult(rec)); 
                recips[j] += rec;
            }
        }
        for i in 0..n { 
            eccs[i].mutsmult(1.0/recips[i]); 
            eccs[i].mutvsub(&self[i]) 
        }
        eccs
    }

    /// Vector `eccentricity` or measure of  
    /// `not being the geometric median` for a member point.
    /// The true geometric median is as yet unknown.
    /// The true geometric median would return zero.
    /// The member point in question is specified by its index `indx`.
    /// This function is suitable for a single member point. 
    /// When eccentricities of all the points are needed, use `eccentricities` above.
    fn eccmember(self, indx: usize) -> Vec<f64> {
        let n = self.len();
        let mut vsum = vec![0_f64; self[0].len()];
        let thisp = &self[indx];
        let mut recip = 0_f64;
        for i in 0..n {
            if i == indx { continue  }; // exclude this point
            let dvmag = self[i].vdist(&thisp);
            if !dvmag.is_normal() { continue } // too close to this one
            let rec = 1.0/dvmag;
            vsum.mutvadd(&self[i].smult(rec)); // add unit vector
            recip += rec // add separately the reciprocals
        }
        vsum.mutsmult(1.0/recip);
        vsum.vsub(&thisp)
    }

    /// Eccentricity vector for a non member point,
    /// while the true geometric median is as yet unknown.
    /// Returns the eccentricity vector.
    /// The true geometric median would return zero vector.
    /// This function is suitable for a single non-member point. 
    fn eccnonmember(self, p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for x in self { 
            let dvmag = x.vdist(&p);
            if !dvmag.is_normal() { continue } // zero distance, safe to ignore
            let rec = 1.0/dvmag;
            vsum.mutvadd(&x.smult(rec)); // add unit vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmult(1.0/recip);
        vsum.vsub(&p)
    }  

    /// Magnitudes of all the vectors in self
    fn mags(self) -> Vec<f64> {
        let mut magsv = Vec::new();
        for v in self {
            magsv.push(v.vmag())
        }
        magsv
    } 

    /// Mean and Std (in MStats struct), Median and quartiles (in Med struct) 
    /// of scalar eccentricities of points in self.
    /// These are new robust measures of a cloud of multidimensional points (or multivariate sample).  
    fn moe(self) -> (MStats, Med) {
        let eccs = self.eccentricities().mags();
        (eccs.ameanstd().unwrap(),eccs.median().unwrap())
    }

    /// Eccentricity defined Medoid and Outlier.
    /// This can give different results to `medoid` above, defined by sums of distances,
    /// especially for the outliers. See tests.rs.  
    /// Consider some point c and some other points, bunched up at a distance r from c.
    /// The sum of their distances will be n*r. Now, spread those points around a circle of radius r from c.
    /// The sum of their distances from c will remain the same but the eccentricity of c will be much reduced.
    /// # Example
    /// ```
    /// use rstats::{Vecf64,VecVec,functions::genvec};
    /// let d = 6_usize;
    /// let pt = genvec(d,24,7,13); // random test data 5x20
    /// let (_medoideccentricity,medei,_outlierecccentricity,outei) = pt.emedoid();
    /// assert_eq!(medei,19); // index of e-medoid
    /// assert_eq!(outei,4);  // index of e-outlier
    /// ```
    fn emedoid(self) -> (f64, usize, f64, usize) {
        self.eccentricities().mags().minmax()
    }

    /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    /// It has (provably) only vector iterative solutions.
    /// Search methods are slow and difficult in highly dimensional space.
    /// Weiszfeld's fixed point iteration formula had known problems with sometimes failing to converge.
    /// Especially, when the points are dense in the close proximity of the gm, or gm coincides with one of them.  
    /// However, these problems are fixed in my new algorithm here.      
    /// There will eventually be a multithreaded version.
    fn nmedian(self, eps: f64) -> Vec<f64> {
        let mut point = self.acentroid(); // start iterating 
        loop {
            let e = self.eccnonmember(&point);
            point.mutvadd(&e);
            if e.vmag() < eps { return point } // make the last small step anyway
        } 
    }
 
    /// Possible first iteration point for geometric medians.
    /// Same as eccnonmember(origin), just saving the zero subtractions
    /// when moving from the origin
    fn firstpoint(self) -> Vec<f64> {
        let mut rsum = 0_f64;
        let mut vsum = vec![0_f64; self[0].len()];
        for thisp in self {
            let mag = thisp.vmag();
            if mag.is_normal() {  
                let invmag = 1.0_f64/mag;
                rsum += invmag;
                vsum.mutvadd(&thisp.smult(invmag)) // accumulate unit vectors
            }
        }
        vsum.smult(1.0/rsum) // scale by the sum of reciprocals
    }

    /// Secant method with recovery from divergence
    /// for finding the geometric median
    fn gmedian(self, eps: f64) -> Vec<f64> {
        let mut p1 = self.acentroid();     
        let e1 = self.eccnonmember(&p1); // eccentricity vector1 
        let mut e1mag = e1.vmag(); 
        let mut p2 = p1.vadd(&e1);   
        loop {
            let e2 = self.eccnonmember(&p2); // eccentricity vector2
            let e2mag = e2.vmag(); 
            if e2mag < eps  { return p2 }; 
            let newp;         
            if e1mag > e2mag {  // eccentricity magnitude decreased, good, employ secant
                newp = p2.vadd(&e2.smult(p1.vsub(&p2).vmag()/(e1mag-e2mag)))                   
            }
            else { // probably overshot the minimum, go nearer, back 
               newp = p2.vadd(&e2.smult(p1.vsub(&p2).vmag()/(e1mag+e2mag)))                    
            } 
            p1 = p2;        
            p2 = newp;  
            e1mag = e2mag             
        }       
    }    

}
