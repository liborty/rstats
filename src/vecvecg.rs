use crate::{here,Med,Stats,Vecg,Vecf64,MutVecf64,VecVecg,VecVec};
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
        let Med{median:mm,..} = m.median() // ignore quartile fields
            .unwrap_or_else(|_|panic!("{} failed to obtain median",here!()));
        let unitzerom =  m.sadd(-mm).vunit();
        self.iter().map(|s|{ 
            let Med{median:ms,..} = s.median()            
                .unwrap_or_else(|_|panic!("{} failed to obtain median",here!()));
            s.sadd(-ms).vunit().dotp(&unitzerom)
        }).collect::<Vec<f64>>()
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
        for v in self { eccs.push(v.vdistf64(gm)) }
        // create sort index of the eccs
        let index = sortidx(&eccs);
        // pick the associated points weights in the order of the sorted eccs
        let mut weights = index.unindexf64(ws,true);
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
            coses.push(p.vsubunit(medmed).dotp(unitzmed)); 
        } 
        // create sort index of the coses
        let index = sortidx(&coses);
        // pick the associated points weights in the same order as the sorted coses
        let mut weights = index.unindexf64(ws,true);
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
        for (x,&w) in self.iter().zip(ws) {  
            let mag:f64 = x.iter().zip(p).map(|(&xi,&pi)|(f64::from(xi)-pi).powi(2)).sum::<f64>().sqrt(); 
            if mag.is_normal() {  // skip point x when its distance is zero
                let rec = f64::from(w)/mag;
                vsum.iter_mut().zip(x).for_each(|(vi,xi)| *vi += rec*f64::from(*xi)); 
                recip += rec // add separately the reciprocals    
            }
        }
        vsum.iter_mut().for_each(|vi| *vi /= recip);
        vsum
    }
    
    /// Estimated (computed) eccentricity vector for a non member point
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
        let mut p = self.wacentroid(ws); // start iterating from the Centre
        loop { // vector iteration till accuracy eps is reached
            let nextp = self.wnxnonmember(ws,&p);          
            if nextp.iter().zip(p).map(|(&xi,pi)|(xi-pi).powi(2)).sum::<f64>() < eps2 { return nextp }; // termination 
            p = nextp
        }
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
                vm.iter().take(i+1).for_each(|vmj| { 
                    cov[covsub] += thisc*vmj;
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
            let vm = selfh.vsub(m);  // subtract zero mean/median vector  
            vm.iter().enumerate().for_each(|(i,&thisc)|            
                vm.iter().take(i+1).for_each(|&vmj| {  // its weighted products up to and including the diagonal 
                    cov[covsub] += w*thisc*vmj;
                    covsub += 1;
                }));
        });
        // now compute the means and return
        cov.mutsmultf64(1_f64/wsum); 
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
    /// their median is returned.
    fn comed(self, m:&[U]) -> Vec<f64> { // m should be the median here 
        let d = self[0].len(); // dimension of the vector(s)
        let mut com:Vec<f64> = Vec::with_capacity((d+1)*d/2); // result vec flat lower triangular array 
        let zs:Vec<Vec<f64>> = self.iter().map(|s| s.vsub(m)).collect(); // zero median vectors
        for i in 0..d { // cross multiplaying the components
            for j in 0..i+1 { // in this order so as to save memory
                let thisproduct:Vec<f64> = zs.iter().map(|v| v[i]*v[j]).collect();
                let Med{median,..} = thisproduct.median().unwrap();
                com.push(median);
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
                let Med{median,..} = thisproduct.median().unwrap();
                com.push(median/wmean);
            }
        };
        com
    }
}
