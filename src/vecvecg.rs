use crate::{error::RError,RE,Stats,TriangMat,Vecg,MutVecg,VecVecg,VecVec};
use indxvec::{Vecops};
use medians::Median;

impl<T,U> VecVecg<T,U> for &[Vec<T>] 
    where T: Copy+PartialOrd+std::fmt::Display,f64:From<T>, 
    U: Copy+PartialOrd+std::fmt::Display,f64:From<U> {

    /// Leftmultiply (column) vector v by (rows of) matrix self
    fn leftmultv(self,v: &[U]) -> Result<Vec<f64>,RE> {
        if self[0].len() != v.len() { return Err(RError::DataError("leftmultv dimensions mismatch")); };
        Ok(self.iter().map(|s| s.dotp(v)).collect())
    }

    /// Dot product of v with column c of matrix self 
    fn columnp(self,c: usize,v: &[U]) -> f64 {
        v.iter().enumerate().map(|(i,&num)| f64::from(num)*f64::from(self[i][c])).sum()        
    }

    /// Rightmultiply (row) vector v by columns of matrix self
    fn rightmultv(self,v: &[U]) -> Result<Vec<f64>,RE> {
        if v.len() != self.len() { return Err(RError::DataError("rightmultv dimensions mismatch")); }; 
        Ok(self[0].iter().enumerate().map(|(c,_)| self.columnp(c,v) ).collect())  
    }

    /// Rectangular Matrices multiplication: self * m.
    /// Returns DataError if lengths of rows of self: `self[0].len()` 
    /// and columns of m: `m.len()` do not match.
    /// Result dimensions are self.len() x m[0].len() 
    fn matmult(self,m: &[Vec<U>]) -> Result<Vec<Vec<f64>>,RE> {
        if self[0].len() != m.len() { return Err(RError::DataError("matmult dimensions mismatch")); }; 
        Ok(self.iter().map(|srow| 
            m[0].iter().enumerate().map(|(colnum,_)| m.columnp(colnum,srow) ).collect()
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

    /// Weighted Centre
    fn wacentroid(self,ws: &[U]) -> Vec<f64> { 
        let mut wsum = 0_f64;
        let mut result = vec![0_f64;self[0].len()]; 
        for (v,&w) in self.iter().zip(ws) { 
            let weight = f64::from(w); // saves converting twice
            wsum += weight;
            result.mutvadd(&v.smult(weight)); };
        result.mutsmult::<f64>(1.0/wsum); // divide by weighted sum to get the mean
        result
    }

    /// Trend computes the vector connecting the geometric medians of two sets of multidimensional points.
    /// This is a robust relationship between two unordered multidimensional sets.
    /// The two sets have to be in the same (dimensional) space but can have different numbers of points.
    fn trend(self, eps:f64, v:Vec<Vec<U>>) -> Result<Vec<f64>,RE> {
        if self[0].len() != v[0].len() { return Err(RError::DataError("trend dimensions mismatch")); };
        Ok(v.gmedian(eps).vsub::<f64>(&self.gmedian(eps)))
    }

    /// Translates the whole set by subtracting vector m.
    /// When m is set to the geometric median, this produces the zero median form.
    /// The geometric median is invariant with respect to rotation,
    /// unlike the often used mean (`acentroid` here), or the quasi median,
    /// both of which depend on the choice of axis.
     fn translate(self, m:&[U]) -> Result<Vec<Vec<f64>>,RE> { 
        if self[0].len() != m.len() { return Err(RError::DataError("translate dimensions mismatch")); }; 
        Ok(self.iter().map(|s| s.vsub(m)).collect())   
    }

    /// Proportions of points along each +/-axis (hemisphere)
    /// Excludes points that are perpendicular to it
    /// Uses only the points specified in idx (e.g. the convex hull).
    /// Self should normally be zero mean/median vectors, 
    /// e.g. `self.translate(&median)`
    fn wtukeyvec(self, idx: &[usize], ws:&[U]) -> Result<Vec<f64>,RE> { 
        let mut wsum = 0_f64; 
        let dims = self[0].len();
        if self.len() != ws.len() { return Err(RError::DataError("wtukeyvec weights number mismatch")); }; 
        let mut hemis = vec![0_f64; 2*dims]; 
        for &i in idx.iter() { 
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
        if self[0].len() != m.len() { return Err(RError::DataError("dependencies dimensions mismatch")); }
        let entropym = m.entropy();
        return self.iter().map(|s| -> Result<f64,RE> {  
            Ok((entropym + s.entropy())/
            s.jointentropy(m)?-1.0)}).collect() 
    }

    /// (Median) correlations of m with each vector in self
    /// Factors out the unit vector of m to save repetition of work
    fn correlations(self, m: &[U]) -> Result<Vec<f64>,RE> {
        if self[0].len() != m.len() { return Err(RError::DataError("correlations dimensions mismatch")); }
        let mm = m.median(); // ignore quartile fields 
        let unitzerom =  m.sadd(-mm).vunit();
        Ok (self.iter().map(|s| { 
            let ms = s.median();   
            s.sadd(-ms).vunit().dotp(&unitzerom)
        }).collect::<Vec<f64>>())
    }

    /// Individual distances from any point v, typically not a member, to all the members of self.    
    fn dists(self, v:&[U]) -> Result<Vec<f64>,RE> {
        if self[0].len() != v.len() { return Err(RError::DataError("dists dimensions mismatch")); }
        Ok(self.iter().map(|p| p.vdist(v)).collect())
    }

    /// Sum of distances from any single point v, typically not a member, 
    /// to all members of self.    
    /// Geometric Median (gm) is defined as the point which minimises this function.
    /// This is relatively expensive measure to compute.
    /// The radius (distance) from gm is far more efficient, once gm has been found.
    fn distsum(self, v: &[U]) -> Result<f64,RE> {
        if self[0].len() != v.len() { return Err(RError::DataError("distsum dimensions mismatch")); }
        Ok(self.iter().map(|p| p.vdist(v)).sum::<f64>())
    }

    /// Sorted weighted radii (eccentricity) magnitudes to all member points from the Geometric Median.
    fn wsortedrads(self, ws: &[U], gm:&[f64]) -> Result<Vec<f64>,RE> {
        if self.len() != ws.len() { 
            return Err(RError::DataError("wsortedrads self and ws lengths mismatch")); };
        if self[0].len() != gm.len() { 
            return Err(RError::DataError("wsortedrads self and gm dimensions mismatch")); };
        let wnorm = ws.len() as f64 / ws.iter().map(|&w|f64::from(w)).sum::<f64>(); 
        Ok (self.iter().zip(ws).map(|(s,&w)| wnorm*f64::from(w)*s.vdist::<f64>(gm))
            .collect::<Vec<f64>>()
            .sorth(true)
        )
    } 

    /// Like wgmparts, except only does one iteration from any non-member point g
    fn wnxnonmember(self, ws:&[U], g:&[f64]) -> Result<(Vec<f64>,Vec<f64>,f64),RE> {
        if self.len() != ws.len() { 
            return Err(RError::DataError("wnxnonmember and ws lengths mismatch")); };
        if self[0].len() != g.len() { 
            return Err(RError::DataError("wnxnonmember self and gm dimensions mismatch")); };
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
            return Err(RError::DataError("wgmedian and ws lengths mismatch")); };
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
            if nextrecsum-recsum < eps { return Ok(nextg) };  // termination test
            g = nextg;
            recsum = nextrecsum;            
        }
    }
    
    /// Like `gmedian` but returns also the sum of unit vecs and the sum of reciprocals. 
    fn wgmparts(self, ws:&[U], eps: f64) -> Result<(Vec<f64>,Vec<f64>,f64),RE> { 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wgmparts and ws lengths mismatch")); };
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
            return Err(RError::DataError("wgmadgm and ws lengths mismatch")); }; 
        Ok(self.iter().map(|v| v.wvdist(ws,wgm)).collect::<Vec<f64>>().median()) 
    }

    /// Covariance matrix for f64 vectors in self. Becomes comediance when 
    /// argument m is the geometric median instead of the centroid.
    /// Since the matrix is symmetric, the missing upper triangular part can be trivially
    /// regenerated for all j>i by: c(j,i) = c(i,j).
    /// The indexing is always in this order: (row,column) (left to right, top to bottom).
    /// The items are flattened into a single vector in this order.
    /// The full 2D matrix can be reconstructed by `symmatrix` in the trait `Stats`.
    fn covar(self, mid:&[U]) -> Result<TriangMat,RE> {
        let d = self[0].len(); // dimension of the vector(s)
        if d != mid.len() { 
            return Err(RError::DataError("covar self and mid dimensions mismatch")); }; 
        let mut cov:Vec<f64> = vec![0_f64; (d+1)*d/2]; // flat lower triangular results array  
        for thisp in self { // adding up covars for all the points
            let mut covsub = 0_usize; // subscript into the flattened array cov
            let vm = thisp.vsub(mid);  // zero mean vector
            vm.iter().enumerate().for_each(|(i,thisc)| 
                // its products up to and including the diagonal (itself)
                vm.iter().take(i+1).for_each(|vmi| { 
                    cov[covsub] += thisc*vmi;
                    covsub += 1;
                }));
            } 
        // now compute the means and return
        let lf = self.len() as f64;
        cov.iter_mut().for_each(|c| *c /= lf); 
        Ok(TriangMat{trans:false,symmetric:true,data:cov})
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
            return Err(RError::DataError("wcovar self and m dimensions mismatch")); }; 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wcovar self and ws lengths mismatch")); }; 
        let mut cov:Vec<f64> = vec![0_f64; (n+1)*n/2]; // flat lower triangular results array
        let mut wsum = 0_f64;
        self.iter().zip(ws).for_each(|(selfh,&wsh)| { // adding up covars for all the points
            let w = f64::from(wsh); wsum += w;
            let mut covsub = 0_usize; // subscript into the flattened array cov 
            let vm = selfh.vsub(m);   // subtract zero mean/median vector  
            vm.iter().enumerate().for_each(|(i,&thisc)|            
                vm.iter().take(i+1).for_each(|&vmi| {  // its weighted products up to and including the diagonal 
                    cov[covsub] += w*thisc*vmi;
                    covsub += 1;
                }));
        });
        // now compute the means and return
        cov.mutsmult::<f64>(1_f64/wsum); 
        Ok(TriangMat{trans:false,symmetric:true,data:cov})
    }

    /// Covariance matrix for f64 vectors in self. Becomes comediance when 
    /// argument m is the geometric median instead of the centroid.
    /// Since the matrix is symmetric, the missing upper triangular part can be trivially
    /// regenerated for all j>i by: c(j,i) = c(i,j).
    /// The indexing is always in this order: (row,column) (left to right, top to bottom).
    /// The items are flattened into a single vector in this order.
    /// The full 2D matrix can be reconstructed by `symmatrix` in the trait `Stats`.
    /// Similar to `covar` above but instead of averaging the covariances over n points, 
    /// their medians are returned.
    fn comed(self, m:&[U]) -> Result<TriangMat,RE> { // m should be the median here 
        let d = self[0].len(); // dimension of the vector(s)
        if d != m.len() { 
            return Err(RError::DataError("comed self and m dimensions mismatch")); }; 
        let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2); // result vec flat lower triangular array 
        let zs:Vec<Vec<f64>> = self.iter().map(|s| s.vsub(m)).collect(); // zero median vectors
        for i in 0..d { // cross multiplaying the components
            for j in 0..i+1 { // in this order so as to save memory
                let thisproduct:Vec<f64> = zs.iter().map(|v| v[i]*v[j]).collect(); 
                com.push(thisproduct.median());
            }
        }
        Ok(TriangMat{trans:false,symmetric:true,data:com})
    }

    /// Covariance matrix for weighted f64 vectors in self. Becomes comediance when 
    /// argument m is the geometric median instead of the centroid.
    /// Since the matrix is symmetric, the missing upper triangular part can be trivially
    /// regenerated for all j>i by: c(j,i) = c(i,j).
    /// The indexing is always in this order: (row,column) (left to right, top to bottom).
    /// The items are flattened into a single vector in this order.
    /// The full 2D matrix can be reconstructed by `symmatrix` in the trait `Stats`.
    /// Similar to `wcovar` above but instead of averaging the covariances over n points, 
    /// their median is returned.
    fn wcomed(self, ws:&[U], m:&[f64]) -> Result<TriangMat,RE> { // m should be the median here 
        let d = self[0].len(); // dimension of the vector(s)
        if d != m.len() { 
            return Err(RError::DataError("wcomed self and m dimensions mismatch")); }; 
        if self.len() != ws.len() { 
            return Err(RError::DataError("wcomed self and ws lengths mismatch")); }; 
        let zs:Vec<Vec<f64>> = self.iter().map(|s| s.vsub(m)).collect(); // zero median vectors
        let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2); // result vec flat lower triangular array 
        let wmean = ws.iter().map(|&w| f64::from(w)).sum::<f64>()/(self.len() as f64); 
        for i in 0..d { // cross multiplaying the components
            for j in 0..i+1 { // in this order so as to save memory
                let thisproduct:Vec<f64> = zs.iter().zip(ws).map(|(v,&w)| f64::from(w)*v[i]*v[j]).collect();  
                com.push(thisproduct.median()/wmean);
            }
        };
        Ok(TriangMat{trans:false,symmetric:true,data:com})
    }
 
}
