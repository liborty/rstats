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
    /// Centroids of halves of points.
    /// Used to initialise the two point iterative method for finding the geometric median
    fn halfcentroid(self) -> (Vec<f64>,Vec<f64>) {
        let n = self.len();
        let nh = n/2;
        let j = self[0].len();
        let mut centre1 = vec![0_f64; j];
        let mut centre2 = vec![0_f64; j];
        for i in 0..nh {
            centre1.mutvadd(&self[i])
        }
        centre1.mutsmult(1.0 / nh as f64);
        for i in nh..n {
            centre2.mutvadd(&self[i]) 
        }
        centre2.mutsmult(1.0 / (n-nh) as f64);

        (centre1,centre2)
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

    /// `eccentricities` finds vectors from each member point towards the geometric median.
    fn eccentricities(self) -> Vec<Vec<f64>> {
        let n = self.len();
        // allocate vectors for the results
        let mut eccs = vec![vec![0_f64; self[0].len()]; n];
        // ecentricities vectors accumulator for all points
        // examine all unique pairings (lower triangular part of symmetric flat matrix)
        for i in 1..n {
            let thisp = &self[i];
            for j in 0..i { 
                // calculate each unit vector between any pair of points just once
                let e = self[j].vsub(&thisp).vunit(); 
                eccs[i].mutvadd(&e);
                // mind the vector's opposite orientations w.r.t. to the two points!
                eccs[j].mutvsub(&e); 
            }
        }
        eccs
    }

    /// Scalar positive measure of `not being the geometric median` for a member point,
    /// while the true geometric median is as yet unknown.
    /// Returns the magnitude of the eccentricity vector.
    /// The true geometric median would return zero.
    /// The member point is specified by its index `indx`.
    /// This function is suitable for a single member point. 
    /// When eccentricities of all the points are needed, use `eccentricities`.
    fn eccentrinset(self, indx: usize) -> f64 {
        let n = self.len();
        let mut vsum = vec![0_f64; self[0].len()];
        let thisp = &self[indx];
        for i in 0..n {
            if i == indx {
                continue;
            }; // exclude this point
            vsum.mutvadd(&self[i].vsub(&thisp).vunit());
        }
        vsum.vmag() / (n - 1) as f64
    }

    /// Eccentricity vector for a non member point,
    /// while the true geometric median is as yet unknown.
    /// Returns the eccentricity vector.
    /// The true geometric median would return zero vector.
    /// This function is suitable for a single non-member point. 
    fn eccnonmember(self, p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        for x in self { 
            vsum.mutvadd(&x.vsub(&p).vunit());
        }
        vsum
    }

    /// Returns (Measure, Eccentricity-Vector) of any point (typically one not belonging to the set).
    /// The first (scalar) part of the result is a positive measure of `not being a median`.
    /// The second part is the eccentricity vector, which always points towards the median.
    /// The vector is of particular value and interest.
    /// This function has no prior knowledge of the actual median.  
    /// This is suitable for a single point. When eccentricities of all the points
    /// are needed, use more efficient `eccentricities`.
    fn veccentr(self, thisp: &[f64]) -> (f64, Vec<f64>) {
        let d = thisp.len();
        let mut vsum = vec![0_f64; d];
        for thatp in self {
            let mut vdif = thatp.vsub(thisp);
            let mag = vdif.vmag();
            if mag.is_normal() {
                vdif.mutsmult(1./mag); // using already computed magnitude to find unit vdif
                vsum.mutvadd(&vdif); // add it to the sum of vector eccentricities
            } // else mag = 0, so just skip thatp, as it is the same as thisp
        }
        (vsum.vmag()/d as f64, vsum)
    }

    /// This convenience wrapper calls `veccentr` and extracts just the eccentricity (residual median error).
    /// Thus this method is the non-member equivalent of `eccentrinset`.
    /// When the eccentricity vector is needed, use `veccentr`
    fn ecc(self, v: &[f64]) -> f64 {
        let (ecc, _) = self.veccentr(v);
        ecc
    }

    /// Magnitudes of all the vectors in self
    fn mags(self) -> Vec<f64> {
        let mut magsv = Vec::new();
        for v in self {
            magsv.push(v.vmag())
        }
        magsv
    }

    /// Scalar measures of eccentricity for the whole set.
    /// The output can be typically passed to `median`
    /// or `minmax` to find the Outlier and the Medoid.
    fn scalarecc(self) -> Vec<f64> {
        let mut scecc = self.mags();
        scecc.mutsmult(1_f64 / self.len() as f64);
        scecc
    }

    /// Mean and Std (in MStats struct), Median and quartiles (in Med struct) 
    /// of scalar eccentricities of points in self.
    /// These are new robust measures of a cloud of multidimensional points (or multivariate sample).  
    fn moe(self) -> (MStats, Med) {
        let eccs = self.eccentricities().scalarecc();
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
    /// assert_eq!(medei,10); // index of e-medoid
    /// assert_eq!(outei,9);  // index of e-outlier
    /// ```
    fn emedoid(self) -> (f64, usize, f64, usize) {
        self.eccentricities().scalarecc().minmax()
    }

    /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    /// It has (provably) only vector iterative solutions.
    /// Search methods are slow and difficult in highly dimensional space.
    /// Weiszfeld's fixed point iteration formula had known problems with sometimes failing to converge.
    /// Especially, when the points are dense in the close proximity of the gm, or it coincides with one of them.  
    /// However, these problems are fixed in my new algorithm here.      
    /// There will eventually be a multithreaded version of `nmedian`.
    /// # Example
    /// ```
    /// use rstats::{VecVec,functions::genvec};
    /// let pt = genvec(15,15,255,30);
    /// let gm = pt.nmedian(1e-5);
    /// let error = pt.ecc(&gm);
    /// assert_eq!(error,0.000004826966175302838_f64);
    /// ```
    fn nmedian(self, eps: f64) -> Vec<f64> {
        let mut oldpoint = self.acentroid(); // start iterating from the centroid
        loop {
            let (rsum, mut newv) = self.betterpoint(&oldpoint);
            newv.mutsmult(rsum); // scaling the returned sum of unit vectors
            // test the magnitude of the move for termination
            if newv.vdist(&oldpoint) < eps {
                oldpoint = newv; // make the last small step anyway
                break; // from the loop
            };
            oldpoint = newv // move to the new point
        }
        oldpoint
    }

    /// Called by nmedian.
    /// Scaling by rsum is left as the final step at calling level,
    /// in order to facilitate data parallelism.
    fn betterpoint(self, v: &[f64]) -> (f64, Vec<f64>) {
        let mut rsum = 0_f64;
        let mut vsum = vec![0_f64; v.len()];
        for thatp in self {
            let dist = v.vdist(&thatp);
            if dist.is_normal() { // exclude points that are too close
                let recip = 1.0 / dist;
                rsum += recip; // accumulate reciprocal scaling weights
                vsum.mutvadd(&thatp.smult(recip)) // accumulate vectors
            }
        }
        (1.0/rsum, vsum)
    }

    /// Trend computes the vector connecting the geometric medians of two sets of multidimensional points.
    /// This is a robust relationship between two unordered multidimensional sets.
    /// The two sets have to be in the same space but can have different numbers of points.
    fn trend(self, eps: f64, v: Vec<Vec<f64>>) -> Vec<f64> {
        let m1 = self.nmedian(eps);
        let m2 = v.nmedian(eps);
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

    /// Iterative two point method for finding the geometric median
    /// without reciprocal scaling
    fn gmedian(self, eps: f64) -> Vec<f64> {
    //    let (mut op1,mut op2) = self.halfcentroid(); // start iterating from the centroids
        let vsize = self[0].len();
        let mut op1 = vec![0_f64; vsize];
        let mut op2 = self.acentroid();
        loop {
            let u = self.eccnonmember(&op1).vunit(); // eccentricity unit vectors
            let v = self.eccnonmember(&op2).vunit(); // for both points
            let uv = u.dotp(&v);
            let pd = op2.vsub(&op1);
            let udotpd = u.dotp(&pd);
            let b = (uv*udotpd-v.dotp(&pd))/(1.0-uv.powi(2));
            let a = udotpd+b*uv;

            let f1 = op1.vadd(&u.smult(a)); // parmetric vector equations 
            let f2 = op2.vadd(&v.smult(b)); // for the new points
            if (f1.vdist(&op1) < eps) || (f2.vdist(&op2) < eps) {    // termination condition, points are close 
                return f2.vadd(&f1).smult(0.5) // return their midpoint
            }
            op1 = f1;
            op2 = f2;
        }
    }

    /// Iterative annealing method for finding the geometric median
    /// Be careful, might fail to converge
    fn smedian(self, eps: f64) -> Vec<f64> {     
        let (mut anneal, mut u) = self.betterpoint(&self.acentroid());
        let mut point = u.smult(anneal);     
        loop {       
            u = self.eccnonmember(&point); // new eccentricity vector
            let mag = u.vmag();         
            let newpoint = point.vadd(&u.smult(anneal)); // move to a new point           
            if mag < eps { return newpoint } // termination 
            anneal = newpoint.vdist(&point)/mag;
            point = newpoint
        }         
    }    
}
