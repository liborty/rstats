use std::iter::FromIterator;

use crate::{ Med, MStats, MinMax, MutVecg, MutVecf64, Stats, Vecg, Vecf64, VecVec, VecVecg};
pub use indxvec::{here,tof64,Printing,merge::*,Indices};

impl<T> VecVec<T> for &[Vec<T>] where T: Copy+PartialOrd+std::fmt::Display,
    f64: From<T> {

    /// Transpose vec of vecs like a matrix
    fn transpose(self) -> Vec<Vec<T>> {
        let n = self.len();
        let d = self[0].len();
        let mut transp:Vec<Vec<T>> = Vec::with_capacity(d);
        for i in 0..d { 
            let mut column = Vec::with_capacity(n);
            for v in self {
                column.push(v[i]);
            }
            transp.push(column); // column becomes row
        }
        transp   
    }

    /// Joint probability density function of n matched slices of the same length
    fn jointpdfn(self) -> Vec<f64> {  
        let d = self[0].len(); // their common dimensionality (length)
        for v in self.iter().skip(1) {
            if v.len() != d { panic!("{} all vectors must be of equal length!",here!()) }; 
        }
        let mut res:Vec<f64> = Vec::with_capacity(d);
        let mut tuples = self.transpose();
        let df = tuples.len() as f64; // for turning counts to probabilities
        println!("{}",df);
        // lexical sort to group together occurrences of identical tuples
        tuples.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap()); 
        let mut count = 1_usize; // running count
        let mut lastindex = 0; // initial index of the last unique tuple
        tuples.iter().enumerate().skip(1).for_each(|(i,ti)|
            if ti > &tuples[lastindex] { // new tuple ti (Vec<T>) encountered
                res.push((count as f64)/df); // save frequency count as probability
                lastindex = i; // current index becomes the new one
                count = 1_usize; // reset counter
            } 
            else { count += 1; } );        
        res.push((count as f64)/df);  // flush the rest!
        res
    } 

    /// Joint entropy of vectors of the same length
    fn jointentropyn(self) -> f64 {
        let jpdf = self.jointpdfn(); 
        jpdf.iter().map(|&x| -x*(x.ln()) ).sum() 
    }

    /// Dependence (component wise) of a set of vectors.
    /// i.e. `dependencen` returns 0 iff they are statistically independent
    /// bigger values when they are dependent
    fn dependencen(self) -> f64 { 
        self.iter().map(|v| v.entropy()).sum::<f64>()/self.jointentropyn() - 1.0
    } 

    /// Flattened lower triangular part of a symmetric matrix for column vectors in self.
    /// The upper triangular part can be trivially generated for all j>i by: c(j,i) = c(i,j).
    /// Applies closure f which computes a scalar relationship between two vectors, 
    /// that is different features stored in columns of self.
    /// The closure typically invokes one of the methods from Vecg trait (in vecg.rs),
    /// such as dependencies or correlations.
    /// Example call: `pts.transpose().crossfeatures(|v1,v2| v1.mediancorr(v2))` 
    /// computes median correlations between all column vectors (features) in pts.
 
    fn crossfeatures<F>(self,f:F) -> Vec<f64> where F: Fn(&[T],&[T]) -> f64 {
        let n = self.len(); // number of the vector(s)
        let mut codp:Vec<f64> = Vec::with_capacity((n+1)*n/2); // results 
        self.iter().enumerate().for_each(|(i,v)| 
            // its dependencies up to and including the diagonal
            self.iter().take(i+1).for_each(|vj| { 
                codp.push(f(v,vj)); 
                }));  
        codp
    }

    /// acentroid = simple multidimensional arithmetic mean
    fn acentroid(self) -> Vec<f64> {
        let mut centre = vec![0_f64; self[0].len()];
        for v in self { centre.mutvadd(v) }
        centre.mutsmultf64(1.0 / (self.len() as f64));
        centre
    }

    /// gcentroid = multidimensional geometric mean
    fn gcentroid(self) -> Vec<f64> {
        let nf = self.len() as f64; // number of points
        let d = self[0].len(); // dimensions
        let mut centre = vec![0_f64; d];
        let mut lnvec = vec![0_f64; d];
        for v in self { 
            v.iter().zip(&mut lnvec).for_each(|(&vi,lni)| *lni = f64::from(vi).ln() );
            centre.mutvaddf64(&lnvec)
        }
        centre.iter().map(|comp| (comp/nf).exp()).collect()        
    }

    /// hcentroid =  multidimensional harmonic mean
    fn hcentroid(self) -> Vec<f64> {
        let mut centre = vec![0_f64; self[0].len()]; 
        for v in self { centre.mutvaddf64(&v.vinverse().unwrap()) }
        centre.smultf64(1.0/(self.len() as f64)).vinverse().unwrap()       
    }

    /// For each member point, gives its sum of distances to all other points and their MinMax
    fn distsums(self) -> Vec<f64> {
        let n = self.len();
        let mut dists = vec![0_f64; n]; // distances accumulator for all points
        // examine all unique pairings (lower triangular part of symmetric flat matrix)
        self.iter().enumerate().for_each(|(i,thisp)| 
            self.iter().take(i).enumerate().for_each(|(j,thatp)| {
                let d = thisp.vdist(thatp); // calculate each distance relation just once
                dists[i] += d;
                dists[j] += d; // but add it to both points' sums
            }));
        dists
    } 

    /// The sum of distances from one member point, given by its `indx`,
    /// to all the other points in self.
    /// For all the points, use more efficient `distsums`.
    /// For measure of 'outlyingness', use radius from gm, in preference (is more efficient).    
    fn distsuminset(self, indx: usize) -> f64 { 
        let thisp = &self[indx];
        self.iter().enumerate().map(|(i,thatp)|
            if i == indx { 0.0 } else {thisp.vdist(thatp)}).sum()
    } 

    /// Medoid and Outlier (Medout)
    /// Medoid is the member point (point belonging to the set of points `self`), 
    /// which has the least sum of distances to all other points.
    /// Outlier is the point with the greatest sum of distances.
    /// In other words, they are the members nearest and furthest from the geometric median. 
    /// Returns struct MinMax{min,minindex,max,maxindex}
    fn medout(self) -> MinMax<f64> {  
        minmax(&self.distsums())
    }

    /// Finds approximate vectors from each member point towards the geometric median.
    /// Twice as fast using symmetry, as doing them individually.
    /// For measure of 'outlyingness' use `exacteccs` below
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
                let dvmag = self[j].vdist(thisp);             
                if !dvmag.is_normal() { continue }
                let rec = 1.0_f64/dvmag;
                eccs[i].mutvaddf64(&self[j].smultf64(rec));
                recips[i] += rec;
                // mind the vector's opposite orientations w.r.t. to the two points!
                eccs[j].mutvsubf64(&self[j].smultf64(rec)); 
                recips[j] += rec; // but scalar distances are the same
            }
        }
        for i in 0..n { 
            eccs[i].mutsmultf64(1.0/recips[i]); 
            eccs[i].mutvsub(&self[i]) 
        }
        eccs
    }

    /// Exact radii (eccentricity) vectors from all member points by first finding the Geometric Median.
    /// More accurate and usually faster as well than the approximate `eccentricities` above,
    /// especially when there are many points.
    fn exacteccs(self, eps:f64) -> Vec<Vec<f64>> {        
        let mut eccs = Vec::with_capacity(self.len()); // Vectors for the results
        let gm:Vec<f64> = self.gmedian(eps);
        for v in self {
            eccs.push(gm.vsub(v))
        }
        eccs
    }

    /// Estimated (computed) eccentricity vector for a member point.
    /// It points towards the geometric median. 
    /// The true geometric median is as yet unknown.
    /// The true geometric median would return zero vector.
    /// The member point in question is specified by its index `indx`.
    /// This function is suitable for a single member point. 
    /// When eccentricities of all the points are wanted, use `exacteccs` above.
    fn eccmember(self, indx: usize) -> Vec<f64> {
        self.nxmember(indx).vsub(&self[indx])
    }
    
    /// Estimated (computed) eccentricity vector for a non member point. 
    /// The true geometric median is as yet unknown.
    /// Returns the eccentricity vector.
    /// The true geometric median would return zero vector.
    /// This function is suitable for a single non-member point. 
    fn eccnonmember(self, p:&[f64]) -> Vec<f64> {
        self.nxnonmember(p).vsubf64(p)
    }

    /// Mean and Std (in MStats struct), Median and quartiles (in Med struct), Median and Outlier (in MinMax struct) 
    /// of scalar eccentricities of points in self.
    /// These are new robust measures of a cloud of multidimensional points (or multivariate sample).  
    fn eccinfo(self, eps: f64) -> (MStats, Med, MinMax<f64>) where Vec<f64>:FromIterator<f64> {
        let gm = self.gmedian(eps);
        let eccs:Vec<f64> = self.iter().map(|v| gm.vdist(v)).collect();
        (eccs.ameanstd().unwrap(),eccs.median().unwrap(),minmax(&eccs))
    }

    /// MADn multidimensional median of radii: stable nd data spread estimator
    fn madn(self, eps: f64) -> f64 {
        let gm = self.gmedian(eps); 
        let eccs:Vec<f64> = self.iter().map(|v| gm.vdist(v)).collect();
        let Med{median,..} = eccs.median().unwrap_or_else(|_| panic!("{},median failed\n",here!()));
        median
    }

    /// Proportions of points along each axis
    fn tukeyvec(self, gm: &[f64]) -> Vec<f64> {
        let nf = self.len() as f64; 
        let dims = self[0].len();
        let mut hemis = vec![0_f64; dims];       
        let zerogm = self.zerogm(gm);
        for v in zerogm {   
            for (i,&component) in v.iter().enumerate() {
                if component > 0. { hemis[i] += 1. }
                else if component < 0. { hemis[i] -= 1. };  
            }
        }
        hemis.iter_mut().for_each(|hem| *hem /= nf );
    hemis 
    }

    /// Mean projections of radii on each axis
    fn radvec(self, gm: &[f64]) -> Vec<f64> { 
        let nf = self.len() as f64; 
        let dims = self[0].len();
        let unitc = 1_f64.powi(-(dims as i32)); // components of unit axis vectors
        let mut projections = vec![0_f64; dims];
        let zerogm = self.zerogm(gm);
        for v in zerogm { 
            // compute projection (dot product) on all the unit axes
            // without explicitly constructing  unit axes vector
            // as all its components are the same unitc
            for (i,&component) in v.iter().enumerate() {
                projections[i] += component*unitc;    
            }
        }
        projections.iter_mut().for_each(|p| *p /= nf);
        projections      
    }
     
    /// GM and sorted eccentricities magnitudes.
    /// Describing a set of points `self` in n dimensions
    fn sortedeccs(self, ascending:bool, gm:&[f64]) -> Vec<f64> { 
        let mut eccs = Vec::with_capacity(self.len()); 
        // collect raw ecentricities magnitudes
        for v in self { eccs.push(v.vdistf64(gm)) }
        sortm(&eccs,ascending)
    }    

    /// Eccentricities of Medoid and Outlier.  
    /// Same as just the third element of a tuple returned by eccinfo
    fn emedoid(self, eps: f64) -> MinMax<f64> where Vec<f64>:FromIterator<f64> {
    let gm:Vec<f64> = self.gmedian(eps);
    let eccs:Vec<f64> = self.iter().map(|v| gm.vdist(v)).collect();
    minmax(&eccs)
} 

    /// Initial (first) point for geometric medians.
    /// Same as eccnonmember('origin') but saving the subtractions of zeroes. 
    fn firstpoint(self) -> Vec<f64> {
        let mut rsum = 0_f64;
        let mut vsum = vec![0_f64; self[0].len()];
        for p in self {
            let mag = p.vmag();
            if mag.is_normal() {  // skip if p is at the origin
                let rec = 1.0_f64/mag;
                // the sum of reciprocals of magnitudes for the final scaling  
                rsum += rec;
                // so not using simply .unitv 
                vsum.mutvaddf64(&p.smultf64(rec)) // add all unit vectors
            }
        }
        vsum.mutsmultf64(1.0/rsum); // scale by the sum of reciprocals
        vsum
    }    
    
    /// Next approximate gm computed from a member point  
    /// specified by its index `indx` to self. 
    fn nxmember(self, indx: usize) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let p = &tof64(&self[indx]);
        let mut recip = 0_f64;
        for (i,x) in self.iter().enumerate() {
            if i != indx {  // not point p
                let mag:f64 = x.iter().zip(p).map(|(&xi,&pi)|(f64::from(xi)-pi).powi(2)).sum::<f64>().sqrt(); 
                if mag.is_normal() { // ignore this point should distance be zero
                    let rec = 1.0_f64/mag;
                    vsum.iter_mut().zip(x).for_each(|(vi,xi)| *vi += rec*f64::from(*xi)); 
                    recip += rec // add separately the reciprocals    
                }
            }
        };
        vsum.iter_mut().for_each(|vi| *vi /= recip);
        vsum
    }
 
    /// Next approximate gm computed from a non-member point p
    fn nxnonmember(self, p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for x in self { 
            let mag:f64 = x.iter().zip(p).map(|(&xi,&pi)|(f64::from(xi)-pi).powi(2)).sum::<f64>().sqrt(); 
            if mag.is_normal() { // ignore this point should distance be zero
                let rec = 1.0_f64/mag;
                vsum.iter_mut().zip(x).for_each(|(vi,xi)| *vi += rec*f64::from(*xi)); 
                recip += rec // add separately the reciprocals    
            }
        }
        vsum.iter_mut().for_each(|vi| *vi /= recip);
        vsum
    }

    /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    /// It has (provably) only vector iterative solutions.
    /// Search methods are slow and difficult in highly dimensional space.
    /// Weiszfeld's fixed point iteration formula had known problems with sometimes failing to converge.
    /// Especially, when the points are dense in the close proximity of the gm, or gm coincides with one of them.  
    /// However, these problems are fixed in my new algorithm here.      
    /// There will eventually be a multithreaded version.
    fn gmedian(self, eps: f64) -> Vec<f64> {
        let eps2 = eps.powi(2);
        let mut p = self.acentroid(); // start iterating from the Centre
        loop { // vector iteration till accuracy eps2 is reached
            let nextp = self.nxnonmember(&p);          
            if nextp.iter().zip(p).map(|(&xi,pi)|(xi-pi).powi(2)).sum::<f64>() < eps2 { return nextp }; // termination 
            p = nextp;
        }; 
    }

    /// Same a gmedian but returns also the number of iterations
    fn igmedian(self, eps: f64) -> ( Vec<f64>, usize ) {  
        let eps2 = eps.powi(2);
        let mut p = self.acentroid(); 
        let mut iterations = 1_usize;
        loop {  
            let nextp = self.nxnonmember(&p);
            if nextp.vdistsqf64(&p) < eps2 { return (nextp,iterations) }; // termination
            iterations +=1;
            p = nextp;
        };     
    }
}
