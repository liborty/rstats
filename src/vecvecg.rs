use crate::{re_error,RError,RE,Stats,TriangMat,Vecg,MutVecg,VecVecg,VecVec};
use indxvec::Mutops;
use medians::Medianf64;
use rayon::prelude::*;

impl<T,U> VecVecg<T,U> for &[Vec<T>] 
    where T: Sync+Clone+PartialOrd+Into<f64>, 
    Vec<Vec<T>>: IntoParallelIterator,
    Vec<T>: IntoParallelIterator,
    U: Sync+Clone+PartialOrd+Into<f64>,
    Vec<Vec<U>>: IntoParallelIterator,
    Vec<U>: IntoParallelIterator {

    /// Applies scalar valued closure to all vectors in self and multiplies by their weights.
    /// Returns also the sum of weights.
    fn scalar_wfn(self, ws: &[U], f: impl Fn(&[T]) -> Result<f64,RE>)
        -> Result<(Vec<f64>,f64),RE> {
        let mut wsum = 0_f64;
        let resvec = self.iter().zip(ws).map(|(s,w)|-> Result<f64,RE> {
            let wf = w.clone().into();
            wsum += wf;
            Ok(wf*f(s)?) }).collect::<Result<Vec<f64>,RE>>()?;
        Ok((resvec,wsum))
    }

    /// Applies vector valued closure to all vectors in self and multiplies by their weights.
    /// Returns also the sum of weights
    fn vector_wfn(self, ws: &[U], f: impl Fn(&[T]) -> Result<Vec<f64>,RE>)
        -> Result<(Vec<Vec<f64>>,f64),RE> {
        let mut wsum = 0_f64;
        let resvecvec = self.iter().zip(ws).map(|(s,w)|-> Result<Vec<f64>,RE> {
            let wf = w.clone().into();
            wsum += wf;            
            Ok(f(s)?.smult(wf))}).collect::<Result<Vec<Vec<f64>>,RE>>()?;
        Ok((resvecvec,wsum))
    } 

    /// Individually time weighted time series derivative of vectors.
    /// Weighted arithmetic mean, minus the centre (geometric median). 
    fn wdvdt(self, ws: &[U], centre: &[f64]) -> Result<Vec<f64>, RE> {
        let len = self.len();
        if len < 2 {
            return re_error("empty","time series too short: {len}");
        };
        let mut weightsum:f64 = ws[0].clone().into();
        let mut sumv:Vec<f64> = self[0].smult(weightsum);
        for i in 1..len { 
            let fws = ws[i].clone().into();
            weightsum += fws;
            sumv.mutvadd(&self[i].smult(fws));
        };
        Ok(sumv.smult(1.0/weightsum).vsub(centre))
    }

    /// 1.0-dotproduct with **v**, in range [0,2] 
    fn divs(self, v: &[U]) -> Result<Vec<f64>,RE> { 
        if self.is_empty() { 
            return re_error("empty","divs given no points")?; }; 
        if self[0].len() != v.len() { 
            return re_error("size","divs dimensions mismatch")?; }; 
        let uv = v.vunit()?;
        self.scalar_fn(|p| Ok(1.0-p.vunit()?.dotp(&uv)))
    }

    /// median of weighted 1.0-dotproducts of **v**, with all in self
    fn wdivs(self, ws:&[U], v: &[f64]) -> Result<(Vec<f64>,f64),RE> { 
        if self.is_empty() { 
            return re_error("empty","wdivs given no points")?; }; 
        if self[0].len() != v.len() { 
            return re_error("size","wdivs dimensions mismatch")?; }; 
        let uv = v.vunit()?;
        self.scalar_wfn(ws,|p| Ok(1.0-p.vunit()?.dotp(&uv)))
    }

    /// median of weighted cos deviations from **v**
    fn wdivsmed(self, ws: &[U], v: &[f64]) -> Result<f64,RE> { 
        if self.is_empty() { 
            return re_error("empty","wmeddivs given no points")?; }; 
        if self[0].len() != v.len() { 
            return re_error("size","wmeddivs dimensions mismatch")?; }; 
        let (values,wsum) = self.wdivs(ws,v)?;
        Ok((self.len() as f64) * values.medf_unchecked()/wsum)
    }
 
    /// weighted radii to all points in self
    fn wradii(self, ws:&[U], gm: &[f64]) -> Result<(Vec<f64>,f64),RE> {
        if self.is_empty() { 
            return re_error("empty","wradii given no points")?; }; 
        if self[0].len() != gm.len() { 
            return re_error("size","wradii dimensions mismatch")?; }; 
        self.scalar_wfn(ws, |p| Ok(p.vdist(gm)))
    }

    /// wmadgm median of weighted deviations from (weighted) gm: stable nd data spread estimator.
    fn wmadgm(self, ws: &[U], gm: &[f64]) -> Result<f64,RE> { 
        if self.is_empty() { 
            return re_error("empty","wmadgm given no points")?; }; 
        if self[0].len() != gm.len() { 
            return re_error("size","wmadgm dimensions mismatch")?; }; 
        let (values,wsum) = self.scalar_wfn(ws,|p| Ok(p.vdist(gm)))?;
        Ok((self.len() as f64) * values.medf_unchecked()/wsum)
    }

    /// Rows of matrix self multiplying (column) vector v
    fn leftmultv(self,v: &[U]) -> Result<Vec<f64>,RE> {
        if self[0].len() != v.len() { 
            return re_error("size","leftmultv dimensions mismatch")?; };
        Ok(self.iter().map(|s| s.dotp(v)).collect())
    }

    /// Row vector v multipying columns of matrix self
    fn rightmultv(self,v: &[U]) -> Result<Vec<f64>,RE> {
        if v.len() != self.len() { 
            return re_error("size","rightmultv dimensions mismatch")?; }; 
        Ok((0..self[0].len()).map(|colnum| v.columnp(colnum,self)).collect())
    }

    /// Rectangular Matrices multiplication: self * m.
    /// Returns DataError if lengths of rows of self: `self[0].len()` 
    /// and columns of m: `m.len()` do not match.
    /// Result dimensions are self.len() x m[0].len() 
    fn matmult(self,m: &[Vec<U>]) -> Result<Vec<Vec<f64>>,RE> {
        if self[0].len() != m.len() { 
            return re_error("size","matmult dimensions mismatch")?; }; 
        Ok(self.par_iter().map(|srow| 
            (0..m[0].len()).map(|colnum| srow.columnp(colnum,m)).collect()
            ).collect::<Vec<Vec<f64>>>()) 
    }

    /// Weighted sum.  
    /// Weights are associated with vectors of self, not with coordinates
    fn wsumv(self,ws: &[U]) -> Vec<f64> {
        let mut resvec = vec![0_f64;self[0].len()]; 
        for (v,w) in self.iter().zip(ws) { 
            let weight:f64 = w.clone().into();
            for (res,component) in resvec.iter_mut().zip(v) {
                *res += weight*(component.clone().into()) }
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
            | mut pair: (Vec<f64>, f64), (p,w) | { 
                let weight:f64 = w.clone().into(); // saves converting twice 
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
        if self[0].len() != m.len() { 
            return re_error("size","translate dimensions mismatch")?; }; 
        self.vector_fn(|s| Ok(s.vsub(m)))   
    }

    /// Weighted sums of points in each hemisphere.
    /// Uses only the points specified in idx (e.g. the convex hull).
    /// Self should normally be zero median vectors, i.e. `self.translate(&median)`.  
    /// The result is normalized to unit vector.
    fn wsigvec(self, idx: &[usize], ws:&[U]) -> Result<Vec<f64>,RE> {   
        let dims = self[0].len();
        if self.len() != ws.len() { return re_error("DataError","wsigvec: weights number mismatch")?; }; 
        let mut hemis = vec![0_f64; 2*dims]; 
        for &i in idx { 
            let wf:f64 = ws[i].clone().into();           
            for (j,component) in self[i].iter().enumerate() {
                let cf:f64 = component.clone().into();
                if cf < 0. { hemis[dims+j] -= wf*cf; }
                else { hemis[j] += wf*cf; };  
            };
        };
        hemis.vunit()
    }

    /// Likelihood of zero median point **p** belonging to zero median data cloud `self`,
    /// based on the cloud's shape outside of normal plane through **p**. 
    /// Returns the weighted sum of unit vectors of its outside points, projected onto unit **p**. 
    /// Index should be in the descending order of magnitudes of self points (for efficiency).
    /// Weights ws are associated 1-1 with the points (vectors) of self.
    fn wdepth(self, descending_index: &[usize], ws:&[U], p: &[f64]) -> Result<f64,RE> {
        let p2 = p.vmagsq();
        let mut sumvec = vec![0_f64;p.len()]; 
        for &i in descending_index {
            let s = &self[i];
            let ssq = s.vmagsq();
            if ssq <= p2 { break; }; // no more outside points
            if s.dotp(p) > p2 { sumvec.mutvadd(&s.smult(ws[i].clone().into()/ssq.sqrt())) };
        };
        Ok(sumvec.dotp(&p.vunit()?))
    }

    /// Dependencies of m on each vector in self
    /// m is typically a vector of outcomes.
    /// Factors out the entropy of m to save repetition of work
    fn dependencies(self, m:&[U]) -> Result<Vec<f64>,RE> {  
        if self[0].len() != m.len() { 
            return re_error("size","dependencies: dimensions mismatch")?; };
        let entropym = m.entropy();
        return self.par_iter().map(|s| -> Result<f64,RE> {  
            Ok((entropym + s.entropy())/
            s.jointentropy(m)?-1.0)}).collect() 
    }

    /// Individual distances from any point v, typically not a member, to all the members of self.    
    fn dists(self, v:&[U]) -> Result<Vec<f64>,RE> {
        if self[0].len() != v.len() { 
            return re_error("size","dists dimensions mismatch")?; }
        self.scalar_fn(|p| Ok(p.vdist(v)))
    }

    /// Sum of distances from any single point v, typically not a member, 
    /// to all members of self.    
    /// Geometric Median (gm) is defined as the point which minimises this function.
    /// This is relatively expensive measure to compute.
    /// The radius (distance) from gm is far more efficient, once gm has been found.
    fn distsum(self, v: &[U]) -> Result<f64,RE> {
        if self[0].len() != v.len() { 
            return re_error("size","distsum dimensions mismatch")?; }
        Ok(self.iter().map(|p| p.vdist(v)).sum::<f64>())
    }

    /// Sorted weighted radii to all member points from the Geometric Median.
    fn wsortedrads(self, ws: &[U], gm:&[f64]) -> Result<Vec<f64>,RE> {
        if self.len() != ws.len() { 
            return re_error("size","wsortedrads self and ws lengths mismatch")?; };
        if self[0].len() != gm.len() { 
            return re_error("size","wsortedrads self and gm dimensions mismatch")?; };
        let wf = ws.iter().map(|x| x.clone().into()).collect::<Vec<f64>>();
        let wnorm = 1.0 / wf.iter().sum::<f64>(); 
        let mut res = self.iter().map(|s| wnorm*s.vdist::<f64>(gm))
            .collect::<Vec<f64>>();
        res.muthashsort(|f| *f);
        Ok(res)
    }

    /// Weighted Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    fn wgmedian(self, ws:&[U], eps: f64) -> Result<Vec<f64>,RE> { 
        if self.len() != ws.len() { 
            return re_error("size","wgmedian and ws lengths mismatch")?; };
        let mut g = self.wacentroid(ws); // start iterating from the weighted centre 
        let mut recsum = 0f64;
        loop { // vector iteration till accuracy eps is exceeded  
            let mut nextg = vec![0_f64; self[0].len()];   
            let mut nextrecsum = 0_f64;
            for (x,w) in self.iter().zip(ws) {   
                let mag = g.iter().zip(x).map(|(&gi,xi)|(xi.clone().into()-gi).powi(2)).sum::<f64>(); 
                if mag.is_normal() { 
                    let rec = w.clone().into()/(mag.sqrt()); // weight/distance (scalar)
                    // vsum increments by components
                    nextg.iter_mut().zip(x).for_each(|(vi,xi)| *vi += xi.clone().into()*rec); 
                    nextrecsum += rec // add separately the reciprocals for final scaling   
                } // else ignore this point should its distance from g be nearly zero
            }
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);       
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
        let mut g = self.wacentroid(ws); // start iterating from the weighted centre or from vec![0_f64; self[0].len()]
        let mut recsum = 0_f64;
        loop {
            // vector iteration till accuracy eps is exceeded
            let (mut nextg, nextrecsum) = self
                .par_iter().zip(ws)
                .fold(
                    || (vec![0_f64; self[0].len()], 0_f64),
                    |mut pair: (Vec<f64>, f64), (p, w)| {
                        // |p-g| done in-place for speed. Could have simply called p.vdist(g)
                        let mag: f64 = p
                            .iter()
                            .zip(&g)
                            .map(|(vi, gi)| (vi.clone().into() - gi).powi(2))
                            .sum();
                        // let (mut vecsum, mut recsum) = pair;
                        if mag > eps {
                            let rec = w.clone().into() / (mag.sqrt()); // reciprocal of distance (scalar)
                            for (vi, gi) in p.iter().zip(&mut pair.0) {
                                *gi += vi.clone().into() * rec
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
    fn wgmparts(self, ws:&[U], eps: f64) -> Result<(Vec<f64>,f64),RE> { 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wgmparts and ws lengths mismatch".to_owned())); };
        let mut g = self.wacentroid(ws); // start iterating from the weighted centre
        let mut recsum = 0f64; 
        loop { // vector iteration till accuracy eps is exceeded  
            let mut nextg = vec![0_f64; self[0].len()];   
            let mut nextrecsum = 0f64;
            for (x,w) in self.iter().zip(ws) { // for all points
                // |x-g| done in-place for speed. Could have simply called x.vdist(g)
                //let mag:f64 = g.vdist::<f64>(&x); 
                let mag = g.iter().zip(x).map(|(&gi,xi)|(xi.clone().into()-gi).powi(2)).sum::<f64>(); 
                if mag.is_normal() { 
                    let rec = w.clone().into()/(mag.sqrt()); // reciprocal of distance (scalar)
                    // vsum increments by components
                    nextg.iter_mut().zip(x).for_each(|(vi,xi)| *vi += xi.clone().into()*rec); 
                    nextrecsum += rec // add separately the reciprocals for final scaling   
                } // else simply ignore this point should its distance from g be zero
            }
            if nextrecsum-recsum < eps { 
                return Ok((
                    nextg.iter().map(|&gi| gi/nextrecsum).collect::<Vec<f64>>(),    
                    nextrecsum
                )); }; // termination        
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);
            g = nextg;
            recsum = nextrecsum;            
        }
    }

    /// Weighted covariance matrix for f64 vectors in self. Becomes comediance when 
    /// argument m is the geometric median instead of the centroid.
    /// Since the matrix is symmetric, the missing upper triangular part can be trivially
    /// regenerated for all j>i by: c(j,i) = c(i,j).
    /// The indexing is always in this order: (row,column) (left to right, top to bottom).
    /// The items are flattened into a single vector in this order.
    /// The full 2D matrix can be reconstructed by `symmatrix` in the trait `Stats`.
    fn wcovar(self, ws:&[U], mid:&[f64]) -> Result<TriangMat,RE> {
        let n = self[0].len(); // dimension of the vector(s)
        if n != mid.len() { 
            return re_error("data","wcovar self and m dimensions mismatch")? }; 
        if self.len() != ws.len() { 
            return re_error("data","wcovar self and ws lengths mismatch")? }; 
        let (mut covsum,wsum) = self
            .par_iter().zip(ws)
            .fold(
                || (vec![0_f64; (n+1)*n/2], 0_f64),
                | mut pair: (Vec<f64>, f64), (p,w) | {
                let mut covsub = 0_usize; // subscript into the flattened array cov
                let vm = p.vsub(mid);  // zero mean vector
                let wf = w.clone().into(); // f64 weight for this point
                vm.iter().enumerate().for_each(|(i,component)| 
                    // its products up to and including the diagonal
                    vm.iter().take(i+1).for_each(|vmi| { 
                        pair.0[covsub] += wf*component*vmi;
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
    
    /// Symmetric covariance matrix for weighted vectors.
	/// Becomes comediance when supplied argument `mid`  
    /// is the geometric median instead of the centroid.
    /// Indexing is always in this order: (row,column) (left to right, top to bottom).
    fn serial_wcovar(self, ws:&[U], mid:&[f64]) -> Result<TriangMat,RE> {
        let d = self[0].len(); // dimension of the vector(s)
        if d != mid.len() { 
            return re_error("data","serial_wcovar self and mid dimensions mismatch")? };
        if self.len() != ws.len() { 
            return re_error("data","serial_wcovar self and ws lengths mismatch")? };  
		let mut covsums = vec![0_f64; (d+1)*d/2];
        let mut wsum = 0_f64; 
 		  for (p,w) in self.iter().zip(ws) { 
            let mut covsub = 0_usize; // subscript into the flattened array cov
            let zp = p.vsub(mid);     // zero mean vector
            let wf:f64 = w.clone().into();
            wsum += wf;
            zp.iter().enumerate().for_each(|(i,thisc)| 
                  // its products up to and including the diagonal 
                  zp.iter().take(i+1).for_each(|otherc| { 
                      covsums[covsub] += wf*thisc*otherc;
                      covsub += 1;
                  }) ) 
        };
        // now compute the means and return 
        for c in covsums.iter_mut() { *c /= wsum }; 
        Ok(TriangMat{ kind:2,data:covsums }) // kind 2 = symmetric, non transposed
    }

}
