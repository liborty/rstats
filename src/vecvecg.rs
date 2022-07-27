use crate::{error::RError,Stats,Vecg,MutVecg,VecVecg,VecVec};
use indxvec::{F64,Vecops,Indices};
use medians::{Median};

impl<T,U> VecVecg<T,U> for &[Vec<T>] 
    where T: Copy+PartialOrd+std::fmt::Display,f64:From<T>, 
    U: Copy+PartialOrd+std::fmt::Display,f64:From<U> {

    /// Leftmultiply (column) vector v by matrix self
    fn leftmultv(self,v: &[U]) -> Result<Vec<f64>,RError> {
        if self[0].len() != v.len() { return Err(RError::DataError); };
        Ok(self.iter().map(|s| s.dotp(v)).collect())
    }

    /// Rightmultiply (row) vector v by matrix self
    fn rightmultv(self,v: &[U]) -> Result<Vec<f64>,RError> {
        if v.len() != self.len() { return Err(RError::DataError); };
        Ok(self.transpose().iter().map(|s| s.dotp(v)).collect())
    }

    /*
    /// Leftmultiply matrix m by matrix self
    fn leftmult(self,m: &[Vec<U>]) -> Result<Vec<Vec<f64>>,RError> {
        if self[0].len() != m.len() { return Err(RError::DataError); };
        self.iter().map(|s| s.dotp(
            m.iter()v).collect()
        Ok(self.iter().map(|s| s.dotp(v)).collect())
    }
    */
    
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
    /// The two sets have to be in the same space but can have different numbers of points.
    fn trend(self, eps: f64, v: Vec<Vec<U>>) -> Vec<f64> {
        let mut m1 = self.gmedian(eps);       
        m1.mutvsub::<f64>(&v.gmedian(eps));
        m1
    }

    /// Translates the whole set by subtracting vector m. Returns Vec of Vecs.
    /// When m is set to the geometric median, this produces the zero median form.
    /// The geometric median is invariant with respect to rotation,
    /// unlike the often used mean (`acentroid` here), or the quasi median,
    /// both of which depend on the choice of axis.
     fn translate(self, m: &[U]) -> Vec<Vec<f64>> {  
        self.iter().map(|s| s.vsub(m)).collect()   
    }

    /// Transforms nd data to zeromedian form
    /// essentially the same as translate but specialised to f64 gms
    fn zerogm(self, gm: &[f64]) -> Vec<Vec<f64>> {  
        self.iter().map(|s| s.vsub::<f64>(gm)).collect()
    }

    /// Dependencies of m on each vector in self
    /// m is typically a vector of outcomes.
    /// Factors out the entropy of m to save repetition of work
    fn dependencies(self, m: &[U]) -> Vec<f64> {  
        let entropym = m.entropy();
        self.iter().map(|s| (entropym + s.entropy())/s.jointentropy(m)-1.0).collect::<Vec<f64>>()
    }

    /// (Median) correlations of m with each vector in self
    /// Factors out the unit vector of m to save repetition of work
    fn correlations(self, m: &[U]) -> Vec<f64> {
        let mm = m.median(); // ignore quartile fields 
        let unitzerom =  m.sadd(-mm).vunit();
        self.iter().map(|s|{ 
            let ms = s.as_slice().median();   
            s.sadd(-ms).vunit().dotp(&unitzerom)
        }).collect::<Vec<f64>>()
    }

    /// Individual distances from any point v, typically not a member, to all the members of self.    
    fn dists(self, v: &[U]) -> Vec<f64> {
        self.iter().map(|p| p.vdist(v)).collect()
    }

    /// Sum of distances from any single point v, typically not a member, 
    /// to all members of self.    
    /// Geometric Median (gm) is defined as the point which minimises this function.
    /// This is relatively expensive measure to compute.
    /// The radius (distance) from gm is far more efficient, once gm has been found.
    fn distsum(self, v: &[U]) -> f64 {
        self.iter().map(|p| p.vdist(v)).sum::<f64>()
    }

    /// Sorted eccentricities magnitudes (radii), w.r.t. weighted geometric median.
    /// associated cummulative probability density function in [0,1] of the weights.
    fn wsortedeccs(self, ws: &[U], gm: &[f64]) -> ( Vec<f64>,Vec<f64> ) where F64:From<T> { 
        let mut eccs = Vec::with_capacity(self.len()); 
        // collect true eccentricities magnitudes
        for v in self { eccs.push(v.vdist::<f64>(gm)) }
        // create sort index of the eccs
        let index = eccs.hashsort_indexed();
        // pick the associated points weights in the order of the sorted eccs
        let weights = index.unindex(ws,true);
        let sumw = weights.iter().map(|&w|f64::from(w)).sum::<f64>();
        let mut accum = 0_f64;
        // accummulate the probabilities in [0,1]
        let probs = weights.iter().map(|&w| {
            accum += f64::from(w);
            accum / sumw }).collect::<Vec<f64>>();     
        ( index.unindex(&eccs, true), probs )
    }

    /// Like wgmparts, except only does one iteration from any non-member point g
    fn wnxnonmember(self, ws:&[U], g:&[f64]) -> (Vec<f64>,Vec<f64>,f64) {
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
        ( vsum.iter().map(|vi| vi / recip).collect::<Vec<f64>>(),        
            vsum,
            recip )
    }    

    /// Weighted Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    fn wgmedian(self, ws:&[U], eps: f64) -> Vec<f64> { 
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
            if nextrecsum-recsum < eps { return nextg };  // termination test
            g = nextg;
            recsum = nextrecsum;            
        }
    }
    
    /// Like `gmedian` but returns also the sum of unit vecs and the sum of reciprocals. 
    fn wgmparts(self, ws:&[U], eps: f64) -> (Vec<f64>,Vec<f64>,f64) { 
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
                return (
                    nextg.iter().map(|&gi| gi/nextrecsum).collect::<Vec<f64>>(),
                    nextg,
                    nextrecsum
                ); }; // termination        
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);
            g = nextg;
            recsum = nextrecsum;            
        }
    }

    /// wmadgm median of weighted absolute deviations from weighted gm: stable nd data spread estimator
    fn wmadgm(self, ws: &[U], wgm: &[f64]) -> f64 {     
        let devs:Vec<f64> = self.iter().map(|v| v.wvdist(ws,wgm)).collect();
        devs.as_slice().median()
    }

    /// Covariance matrix for f64 vectors in self. Becomes comediance when 
    /// argument m is the geometric median instead of the centroid.
    /// Since the matrix is symmetric, the missing upper triangular part can be trivially
    /// regenerated for all j>i by: c(j,i) = c(i,j).
    /// The indexing is always in this order: (row,column) (left to right, top to bottom).
    /// The items are flattened into a single vector in this order.
    /// The full 2D matrix can be reconstructed by `symmatrix` in the trait `Stats`.
    fn covar(self, m:&[U]) -> Vec<f64> {
        let d = self[0].len(); // dimension of the vector(s)
        let mut cov:Vec<f64> = vec![0_f64; (d+1)*d/2]; // flat lower triangular results array  
        for thisp in self { // adding up covars for all the points
            let mut covsub = 0_usize; // subscript into the flattened array cov
            let vm = thisp.vsub(m);  // zero mean vector
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
        cov
    }
 
    /// Weighted covariance matrix for f64 vectors in self. Becomes comediance when 
    /// argument m is the geometric median instead of the centroid.
    /// Since the matrix is symmetric, the missing upper triangular part can be trivially
    /// regenerated for all j>i by: c(j,i) = c(i,j).
    /// The indexing is always in this order: (row,column) (left to right, top to bottom).
    /// The items are flattened into a single vector in this order.
    /// The full 2D matrix can be reconstructed by `symmatrix` in the trait `Stats`.
    fn wcovar(self, ws:&[U], m:&[f64]) -> Vec<f64> {
        let n = self[0].len(); // dimension of the vector(s)
        // let mut covs:Vec<Vec<f64>> = Vec::new();
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
        cov
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
    fn comed(self, m:&[U]) -> Vec<f64> { // m should be the median here 
        let d = self[0].len(); // dimension of the vector(s)
        let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2); // result vec flat lower triangular array 
        let zs:Vec<Vec<f64>> = self.iter().map(|s| s.vsub(m)).collect(); // zero median vectors
        for i in 0..d { // cross multiplaying the components
            for j in 0..i+1 { // in this order so as to save memory
                let thisproduct:Vec<f64> = zs.iter().map(|v| v[i]*v[j]).collect(); 
                com.push(thisproduct.as_slice().median());
            }
        }
        com
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
    fn wcomed(self, ws:&[U], m:&[f64]) -> Vec<f64> { // m should be the median here 
        let d = self[0].len(); // dimension of the vector(s)
        let zs:Vec<Vec<f64>> = self.iter().map(|s| s.vsub(m)).collect(); // zero median vectors
        let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2); // result vec flat lower triangular array 
        let wmean = ws.iter().map(|&w| f64::from(w)).sum::<f64>()/(self.len() as f64); 
        for i in 0..d { // cross multiplaying the components
            for j in 0..i+1 { // in this order so as to save memory
                let thisproduct:Vec<f64> = zs.iter().zip(ws).map(|(v,&w)| f64::from(w)*v[i]*v[j]).collect();  
                com.push(thisproduct.as_slice().median()/wmean);
            }
        };
        com
    }
}
