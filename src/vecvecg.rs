use crate::{noop,fromop,error::RError,RE,Stats,TriangMat,Vecg,MutVecg,VecVecg,VecVec};
use indxvec::Vecops;
use medians::{Median,Medianf64};
use rayon::prelude::*;

impl<T,U> VecVecg<T,U> for &[Vec<T>] 
    where T: Sync+Copy+PartialOrd+std::fmt::Display,f64:From<T>, 
    Vec<Vec<T>>: IntoParallelIterator,
    Vec<T>: IntoParallelIterator,
    U: Sync+Copy+PartialOrd+std::fmt::Display,f64:From<U>,
    Vec<Vec<U>>: IntoParallelIterator,
    Vec<U>: IntoParallelIterator {

    /// Leftmultiply (column) vector v by (rows of) matrix self
    fn leftmultv(self,v: &[U]) -> Result<Vec<f64>,RE> {
        if self[0].len() != v.len() { return Err(RError::DataError(
            "leftmultv dimensions mismatch".to_owned())); };
        Ok(self.par_iter().map(|s| s.dotp(v)).collect())
    }

    /// Rightmultiply (row) vector v by columns of matrix self
    fn rightmultv(self,v: &[U]) -> Result<Vec<f64>,RE> {
        if v.len() != self.len() { return Err(RError::DataError("rightmultv dimensions mismatch".to_owned())); }; 
        Ok((0..self[0].len()).into_par_iter().map(|colnum| v.columnp(colnum,self)).collect())
    }

    /// Rectangular Matrices multiplication: self * m.
    /// Returns DataError if lengths of rows of self: `self[0].len()` 
    /// and columns of m: `m.len()` do not match.
    /// Result dimensions are self.len() x m[0].len() 
    fn matmult(self,m: &[Vec<U>]) -> Result<Vec<Vec<f64>>,RE> {
        if self[0].len() != m.len() { return Err(RError::DataError("matmult dimensions mismatch".to_owned())); }; 
        Ok(self.par_iter().map(|srow| 
            (0..m[0].len()).map(|colnum| srow.columnp(colnum,m)).collect()
            ).collect::<Vec<Vec<f64>>>()) 
    }

    /// Weighted sum of nd points (or vectors).  
    /// Weights are associated with points, not coordinates
    fn wsumv(self,ws: &[U]) -> Vec<f64> {
        let mut resvec = vec![0_f64;self[0].len()]; 
        for (v,&w) in self.iter().zip(ws) { 
            let weight = f64::from(w);
            for (res,component) in resvec.iter_mut().zip(v) {
                *res += weight*f64::from(*component) }
        };
        resvec
    }

    /// Weighted Centre.  
    /// Weights are associated with points
    fn wacentroid(self,ws: &[U]) -> Vec<f64> {  
        let (sumvec,weightsum) = self
        .par_iter().zip(ws)
        .fold(
            || (vec![0_f64;self[0].len()], 0_f64),
            | mut pair: (Vec<f64>, f64), (p,&w) | { 
                let weight = f64::from(w); // saves converting twice 
                pair.0.mutvadd(&p.smult(weight));
                pair.1 += weight;
                pair
            }
        )
        .reduce(
            || (vec![0_f64; self[0].len()], 0_f64),
            | mut pairsum: (Vec<f64>, f64), pairin: (Vec<f64>, f64)| {
                pairsum.0.mutvadd::<f64>(&pairin.0);
                pairsum.1 += pairin.1;
                pairsum
            }
        );
        sumvec.smult(1./weightsum)
    }

    /// Trend computes the vector connecting the geometric medians of two sets of multidimensional points.
    /// This is a robust relationship between two unordered multidimensional sets.
    /// The two sets have to be in the same (dimensional) space but can have different numbers of points.
    fn trend(self, eps:f64, v:Vec<Vec<U>>) -> Result<Vec<f64>,RE> {
        if self[0].len() != v[0].len() { return Err(RError::DataError("trend dimensions mismatch".to_owned())); };
        let pair = rayon::join(||v.gmedian(eps),||self.gmedian(eps));
        Ok(pair.0.vsub::<f64>(&pair.1))
    }

    /// Translates the whole set by subtracting vector m.
    /// When m is set to the geometric median, this produces the zero median form.
    /// The geometric median is invariant with respect to rotation,
    /// unlike the often used mean (`acentroid` here), or the quasi median,
    /// both of which depend on the choice of axis.
     fn translate(self, m:&[U]) -> Result<Vec<Vec<f64>>,RE> { 
        if self[0].len() != m.len() { return Err(RError::DataError("translate dimensions mismatch".to_owned())); }; 
        Ok(self.par_iter().map(|s| s.vsub(m)).collect())   
    }

    /// Proportions of points along each +/-axis (hemisphere)
    /// Excludes points that are perpendicular to it
    /// Uses only the points specified in idx (e.g. the convex hull).
    /// Self should normally be zero mean/median vectors, 
    /// e.g. `self.translate(&median)`
    fn wtukeyvec(self, idx: &[usize], ws:&[U]) -> Result<Vec<f64>,RE> { 
        let mut wsum = 0_f64; 
        let dims = self[0].len();
        if self.len() != ws.len() { return Err(RError::DataError("wtukeyvec weights number mismatch".to_owned())); }; 
        let mut hemis = vec![0_f64; 2*dims]; 
        for &i in idx { 
            let wf = f64::from(ws[i]);
            wsum += wf;
            // let zerogm = self[i].vsub::<f64>(gm);
            for (j,&component) in self[i].iter().enumerate() {
                let cf = f64::from(component);
                if cf > 0. { hemis[j] += wf }
                else if cf < 0. { hemis[dims+j] += wf };  
            }
        }
        hemis.iter_mut().for_each(|hem| *hem /= wsum );
        Ok(hemis)
    } 

    /// Dependencies of m on each vector in self
    /// m is typically a vector of outcomes.
    /// Factors out the entropy of m to save repetition of work
    fn dependencies(self, m:&[U]) -> Result<Vec<f64>,RE> {  
        if self[0].len() != m.len() { return Err(RError::DataError("dependencies dimensions mismatch".to_owned())); }
        let entropym = m.entropy();
        return self.par_iter().map(|s| -> Result<f64,RE> {  
            Ok((entropym + s.entropy())/
            s.jointentropy(m)?-1.0)}).collect() 
    }

    /// (Median) correlations of m with each vector in self
    /// Factors out the unit vector of m to save repetition of work
    fn correlations(self, m: &[U]) -> Result<Vec<f64>,RE> {
        if self[0].len() != m.len() { return Err(RError::DataError("correlations dimensions mismatch".to_owned())); } 
        let unitzerom =  m.zeromedian( &mut fromop)?.vunit();
        self.par_iter().map(|s| -> Result<f64,RE> {
            Ok(unitzerom.dotp(&s.zeromedian(&mut fromop)?.vunit())) } )
            .collect::<Result<Vec<f64>,RE>>() 
    } 

    /// Individual distances from any point v, typically not a member, to all the members of self.    
    fn dists(self, v:&[U]) -> Result<Vec<f64>,RE> {
        if self[0].len() != v.len() { return Err(RError::DataError("dists dimensions mismatch".to_owned())); }
        Ok(self.par_iter().map(|p| p.vdist(v)).collect())
    }

    /// Sum of distances from any single point v, typically not a member, 
    /// to all members of self.    
    /// Geometric Median (gm) is defined as the point which minimises this function.
    /// This is relatively expensive measure to compute.
    /// The radius (distance) from gm is far more efficient, once gm has been found.
    fn distsum(self, v: &[U]) -> Result<f64,RE> {
        if self[0].len() != v.len() { return Err(RError::DataError("distsum dimensions mismatch".to_owned())); }
        Ok(self.par_iter().map(|p| p.vdist(v)).sum::<f64>())
    }

    /// Sorted weighted radii to all member points from the Geometric Median.
    fn wsortedrads(self, ws: &[U], gm:&[f64]) -> Result<Vec<f64>,RE> {
        if self.len() != ws.len() { 
            return Err(RError::DataError("wsortedrads self and ws lengths mismatch".to_owned())); };
        if self[0].len() != gm.len() { 
            return Err(RError::DataError("wsortedrads self and gm dimensions mismatch".to_owned())); };
        let wnorm = ws.len() as f64 / ws.iter().map(|&w|f64::from(w)).sum::<f64>(); 
        Ok (self.par_iter().zip(ws).map(|(s,&w)| wnorm*f64::from(w)*s.vdist::<f64>(gm))
            .collect::<Vec<f64>>()
            .sorth(&mut noop,true)
        )
    } 

    /// Like wgmparts, except only does one iteration from any non-member point g
    fn wnxnonmember(self, ws:&[U], g:&[f64]) -> Result<(Vec<f64>,Vec<f64>,f64),RE> {
        if self.len() != ws.len() { 
            return Err(RError::DataError("wnxnonmember and ws lengths mismatch".to_owned())); };
        if self[0].len() != g.len() { 
            return Err(RError::DataError("wnxnonmember self and gm dimensions mismatch".to_owned())); };
        // vsum is the sum vector of unit vectors towards the points
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for (x,&w) in self.iter().zip(ws) { 
            // |x-p| done in-place for speed. Could have simply called x.vdist(p)
            let mag:f64 = x.iter().zip(g).map(|(&xi,&gi)|(f64::from(xi)-gi).powi(2)).sum::<f64>(); 
            if mag.is_normal() { // ignore this point should distance be zero
                let rec = f64::from(w)/(mag.sqrt()); // reciprocal of distance (scalar)
                // vsum increments by components
                vsum.iter_mut().zip(x).for_each(|(vi,xi)| *vi += f64::from(*xi)*rec); 
                 recip += rec // add separately the reciprocals for final scaling   
            }
        }
        Ok(( vsum.iter().map(|vi| vi / recip).collect::<Vec<f64>>(),        
            vsum,
            recip ))
    }    

    /// Weighted Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    fn wgmedian(self, ws:&[U], eps: f64) -> Result<Vec<f64>,RE> { 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wgmedian and ws lengths mismatch".to_owned())); };
        let mut g = self.acentroid(); // start iterating from the Centre 
        let mut recsum = 0f64;
        loop { // vector iteration till accuracy eps is exceeded  
            let mut nextg = vec![0_f64; self[0].len()];   
            let mut nextrecsum = 0_f64;
            for (x,&w) in self.iter().zip(ws) {   
                // |x-g| done in-place for speed. Could have simply called x.vdist(g)
                //let mag:f64 = g.vdist::<f64>(&x); 
                let mag = g.iter().zip(x).map(|(&gi,&xi)|(f64::from(xi)-gi).powi(2)).sum::<f64>(); 
                if mag.is_normal() { 
                    let rec = f64::from(w)/(mag.sqrt()); // reciprocal of distance (scalar)
                    // vsum increments by components
                    nextg.iter_mut().zip(x).for_each(|(vi,&xi)| *vi += f64::from(xi)*rec); 
                    nextrecsum += rec // add separately the reciprocals for final scaling   
                } // else simply ignore this point should its distance from g be zero
            }
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);       
            // eprintln!("recsum {}, nextrecsum {} diff {}",recsum,nextrecsum,nextrecsum-recsum);
            if nextrecsum-recsum < eps { return Ok(nextg); };  // termination test
            g = nextg;
            recsum = nextrecsum;            
        }
    }

    /// Parallel (multithreaded) implementation of the weighted Geometric Median.
    /// Possibly the fastest you will find.  
    /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.  
    /// It has (provably) only vector iterative solutions.    
    /// Search methods are slow and difficult in hyper space.    
    /// Weiszfeld's fixed point iteration formula has known problems and sometimes fails to converge.  
    /// Specifically, when the points are dense in the close proximity of the gm, or gm coincides with one of them.    
    /// However, these problems are solved in my new algorithm here.     
    /// The sum of reciprocals is strictly increasing and so is used to easily evaluate the termination condition.  
    fn par_wgmedian(self, ws: &[U], eps: f64) -> Result<Vec<f64>,RE> { 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wgmedian and ws lengths mismatch".to_owned())); };
        let mut g = self.acentroid(); // start iterating from the mean  or vec![0_f64; self[0].len()];
        let mut recsum = 0_f64;
        loop {
            // vector iteration till accuracy eps is exceeded
            let (mut nextg, nextrecsum) = self
                .par_iter().zip(ws)
                .fold(
                    || (vec![0_f64; self[0].len()], 0_f64),
                    |mut pair: (Vec<f64>, f64), (p, &w)| {
                        // |p-g| done in-place for speed. Could have simply called p.vdist(g)
                        let mag: f64 = p
                            .iter()
                            .zip(&g)
                            .map(|(&vi, gi)| (f64::from(vi) - gi).powi(2))
                            .sum();
                        // let (mut vecsum, mut recsum) = pair;
                        if mag > eps {
                            let rec = f64::from(w) / (mag.sqrt()); // reciprocal of distance (scalar)
                            for (vi, gi) in p.iter().zip(&mut pair.0) {
                                *gi += f64::from(*vi) * rec
                            }
                            pair.1 += rec; // add separately the reciprocals for the final scaling
                        } // else simply ignore this point should its distance from g be zero
                        pair
                    }
                )
                // must run reduce on the partial sums produced by fold
                .reduce(
                    || (vec![0_f64; self[0].len()], 0_f64),
                    |mut pairsum: (Vec<f64>, f64), pairin: (Vec<f64>, f64)| {
                        pairsum.0.mutvadd::<f64>(&pairin.0);
                        pairsum.1 += pairin.1;
                        pairsum
                    }
                );
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);
            if nextrecsum - recsum < eps {
                return Ok(nextg);
            }; // termination test
            g = nextg;
            recsum = nextrecsum;
        }
    }
    
    /// Like `gmedian` but returns also the sum of unit vecs and the sum of reciprocals. 
    fn wgmparts(self, ws:&[U], eps: f64) -> Result<(Vec<f64>,Vec<f64>,f64),RE> { 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wgmparts and ws lengths mismatch".to_owned())); };
        let mut g = self.acentroid(); // start iterating from the Centre
        let mut recsum = 0f64; 
        loop { // vector iteration till accuracy eps is exceeded  
            let mut nextg = vec![0_f64; self[0].len()];   
            let mut nextrecsum = 0f64;
            for (x,&w) in self.iter().zip(ws) { // for all points
                // |x-g| done in-place for speed. Could have simply called x.vdist(g)
                //let mag:f64 = g.vdist::<f64>(&x); 
                let mag = g.iter().zip(x).map(|(&gi,&xi)|(f64::from(xi)-gi).powi(2)).sum::<f64>(); 
                if mag.is_normal() { 
                    let rec = f64::from(w)/(mag.sqrt()); // reciprocal of distance (scalar)
                    // vsum increments by components
                    nextg.iter_mut().zip(x).for_each(|(vi,&xi)| *vi += f64::from(xi)*rec); 
                    nextrecsum += rec // add separately the reciprocals for final scaling   
                } // else simply ignore this point should its distance from g be zero
            }
            if nextrecsum-recsum < eps { 
                return Ok((
                    nextg.iter().map(|&gi| gi/nextrecsum).collect::<Vec<f64>>(),
                    nextg,
                    nextrecsum
                )); }; // termination        
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);
            g = nextg;
            recsum = nextrecsum;            
        }
    }

    /// wmadgm median of weighted absolute deviations from weighted gm: stable nd data spread estimator.
    /// Here the weights are associated with the dimensions, not the points!
    fn wmadgm(self, ws: &[U], wgm: &[f64]) -> Result<f64,RE> { 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wgmadgm and ws lengths mismatch".to_owned())); }; 
        Ok( self
            .par_iter()
            .map(|v| v.wvdist(ws,wgm))
            .collect::<Vec<f64>>()
            .medianf64()?
        ) 
    }

    /// Symmetric covariance matrix. Becomes comediance when argument `mid`  
    /// is the geometric median instead of the centroid.
    /// The indexing is always in this order: (row,column) (left to right, top to bottom).
    /// The items are flattened into a single vector in this order.
    fn covar(self, mid:&[U]) -> Result<TriangMat,RE> {
        let d = self[0].len(); // dimension of the vector(s)
        if d != mid.len() { 
            return Err(RError::DataError("covar self and mid dimensions mismatch".to_owned())); }; 
        let mut covsum = self
            .par_iter()
            .fold(
                || vec![0_f64; (d+1)*d/2],
                | mut cov: Vec<f64>, p | {
                let mut covsub = 0_usize; // subscript into the flattened array cov
                let vm = p.vsub(mid);  // zero mean vector
                vm.iter().enumerate().for_each(|(i,thisc)| 
                    // its products up to and including the diagonal (itself)
                    vm.iter().take(i+1).for_each(|vmi| { 
                        cov[covsub] += thisc*vmi;
                        covsub += 1;
                        })); 
                cov 
                }
            )
            .reduce(
                || vec![0_f64; (d+1)*d/2],
                | mut covout: Vec<f64>, covin: Vec<f64> | {
                covout.mutvadd(&covin);
                covout
                }
            ); 
        // now compute the means and return
        let lf = self.len() as f64;
        covsum.iter_mut().for_each(|c| *c /= lf); 
        Ok(TriangMat{ kind:2,data:covsum }) // symmetric, non transposed
    }
 
    /// Weighted covariance matrix for f64 vectors in self. Becomes comediance when 
    /// argument m is the geometric median instead of the centroid.
    /// Since the matrix is symmetric, the missing upper triangular part can be trivially
    /// regenerated for all j>i by: c(j,i) = c(i,j).
    /// The indexing is always in this order: (row,column) (left to right, top to bottom).
    /// The items are flattened into a single vector in this order.
    /// The full 2D matrix can be reconstructed by `symmatrix` in the trait `Stats`.
    fn wcovar(self, ws:&[U], m:&[f64]) -> Result<TriangMat,RE> {
        let n = self[0].len(); // dimension of the vector(s)
        if n != m.len() { 
            return Err(RError::DataError("wcovar self and m dimensions mismatch".to_owned())); }; 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wcovar self and ws lengths mismatch".to_owned())); }; 
        let (mut covsum,wsum) = self
            .par_iter().zip(ws)
            .fold(
                || (vec![0_f64; (n+1)*n/2], 0_f64),
                | mut pair: (Vec<f64>, f64), (p,&w) | {
                let mut covsub = 0_usize; // subscript into the flattened array cov
                let vm = p.vsub(m);  // zero mean vector
                let wf = f64::from(w); // f64 weight for this point
                vm.iter().enumerate().for_each(|(i,thisc)| 
                    // its products up to and including the diagonal (itself)
                    vm.iter().take(i+1).for_each(|vmi| { 
                        pair.0[covsub] += wf*thisc*vmi;
                        covsub += 1;
                        }));
                pair.1 += wf; 
                pair 
                }
            )
            .reduce(
                || (vec![0_f64; (n+1)*n/2], 0_f64),
                | mut pairout: (Vec<f64>,f64), pairin: (Vec<f64>,f64) | {
                pairout.0.mutvadd(&pairin.0);
                pairout.1 += pairin.1;
                pairout 
                }
            ); 
        // now compute the means and return  
        covsum.iter_mut().for_each(|c| *c /= wsum); 
        Ok(TriangMat{ kind:2,data:covsum }) // symmetric, non transposed
        }
}
