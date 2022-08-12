use std::iter::FromIterator;

use crate::{ RE, RError, MStats, MinMax, MutVecg, Stats, Vecg, VecVec};
use indxvec::{Vecops};
use medians::{Med,Median};

impl<T> VecVec<T> for &[Vec<T>] 
    where T: Copy+PartialOrd+std::fmt::Display,f64:From<T> {

    /// Transpose vec of vecs matrix
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
    fn jointpdfn(self) -> Result<Vec<f64>,RE> {  
        let d = self[0].len(); // their common dimensionality (length)
        for v in self.iter().skip(1) {
            if v.len() != d { 
                return Err(RError::DataError("jointpdfn: all vectors must be of equal length!")); }; 
        }
        let mut res:Vec<f64> = Vec::with_capacity(d);
        let mut tuples = self.transpose();
        let df = tuples.len() as f64; // for turning counts to probabilities
        // println!("{}",df);
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
        Ok(res)
    } 

    /// Joint entropy of vectors of the same length
    fn jointentropyn(self) -> Result<f64,RE> {
        let jpdf = self.jointpdfn()?; 
        Ok(jpdf.iter().map(|&x| -x*(x.ln())).sum()) 
    }

    /// Dependence (component wise) of a set of vectors.
    /// i.e. `dependencen` returns 0 iff they are statistically independent
    /// bigger values when they are dependent
    fn dependencen(self) -> Result<f64,RE> { 
        Ok(self.iter().map(|v| v.entropy()).sum::<f64>()/self.jointentropyn()? - 1.0)
    } 

    /// Flattened lower triangular part of a symmetric matrix for column vectors in self.
    /// The upper triangular part can be trivially generated for all j>i by: c(j,i) = c(i,j).
    /// Applies closure f which computes a scalar relationship between two vectors, 
    /// that is different features stored in columns of self.
    /// The closure typically invokes one of the methods from Vecg trait (in vecg.rs),
    /// such as dependencies or correlations.
    /// Example call: `pts.transpose().crossfeatures(|v1,v2| v1.mediancorr(v2))` 
    /// computes median correlations between all column vectors (features) in pts. 
    fn crossfeatures(self,f:fn(&[T],&[T])->Result<f64,RE>) -> Result<Vec<f64>,RE> {
        let n = self.len(); // number of the vector(s)
        let mut codp:Vec<f64> = Vec::with_capacity((n+1)*n/2); // results 
        for (i,v) in self.iter().enumerate() {
            // its dependencies up to and including the diagonal
            for vj in self.iter().take(i+1) { 
                codp.push(f(v,vj)?); 
            };
        };
        Ok(codp)
    }

    /// Sum of nd points (or vectors)
    fn sumv(self) -> Vec<f64> {
        let mut resvec = vec![0_f64;self[0].len()]; 
        for v in self { resvec.mutvadd(v) };
        resvec
    }

    /// acentroid = multidimensional arithmetic mean
    fn acentroid(self) -> Vec<f64> {
        let mut centre = vec![0_f64; self[0].len()];
        for v in self { centre.mutvadd(v) };
        centre.mutsmult::<f64>(1.0 / (self.len() as f64));
        centre
    }

    /// gcentroid = multidimensional geometric mean
    fn gcentroid(self) -> Vec<f64> {
        let nf = self.len() as f64; // number of points
        let dim = self[0].len(); // dimensions
        let mut result = vec![0_f64; dim];
        for d in 0..dim {
            for v in self { 
                result[d] += f64::from(v[d]).ln(); 
            }
            result[d] /= nf;
            result[d] = result[d].exp() 
        }
        result
    }

    /// hcentroid =  multidimensional harmonic mean
    fn hcentroid(self) -> Vec<f64> {
        let mut centre = vec![0_f64; self[0].len()]; 
        for v in self { centre.mutvadd::<f64>(&v.vinverse().unwrap()) }
        centre.smult::<f64>(1.0/(self.len() as f64)).vinverse().unwrap()       
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
    /// For measure of 'outlyingness', use nore efficient radius from gm.    
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
    fn medout(self,gm:&[f64]) -> MinMax<f64> {  
        self.iter().map(|s| s.vdist::<f64>(gm)).collect::<Vec<f64>>().minmax() }

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
                eccs[i].mutvadd::<f64>(&self[j].smult::<f64>(rec));
                recips[i] += rec;
                // mind the vector's opposite orientations w.r.t. to the two points!
                eccs[j].mutvsub::<f64>(&self[j].smult::<f64>(rec)); 
                recips[j] += rec; // but scalar distances are the same
            }
        }
        for i in 0..n { 
            eccs[i].mutsmult::<f64>(1.0/recips[i]); 
            eccs[i].mutvsub(&self[i]) 
        }
        eccs
    }

    /// Exact radii (eccentricity) magnitudes for all member points from the Geometric Median.
    /// More accurate and usually faster as well than the approximate `eccentricities` above,
    /// especially when there are many points.
    fn radii(self, gm:&[f64]) -> Vec<f64> { 
        self.iter().map(|s| s.vdist::<f64>(gm)).collect::<Vec<f64>>()
    } 

    /// Arith mean and std (in MStats struct), Median info (in Med struct), Medoid and Outlier (in MinMax struct) 
    /// of scalar radii (eccentricities) of points in self.
    /// These are new robust measures of a cloud of multidimensional points (or multivariate sample).  
    fn eccinfo(self, gm:&[f64]) -> (MStats, Med, MinMax<f64>) where Vec<f64>:FromIterator<f64> {
        let rads:Vec<f64> = self.radii(gm);
        (rads.ameanstd().unwrap(),rads.medinfo(),rads.minmax())
    }

    /// Quasi median, recommended only for comparison purposes
    fn quasimedian(self) -> Vec<f64> {
        self.transpose()
            .iter()
            .map(|p| p.median())
            .collect()
    }

    /// Geometric median's estimated error
    fn gmerror(self,g:&[f64]) -> f64 {
        let (gm,_,_) = self.nxnonmember(g); 
        gm.vdist::<f64>(g)
    }

    /// MADGM median of absolute deviations from gm: stable nd data spread estimator
    fn madgm(self, gm: &[f64]) -> f64 {     
        let devs:Vec<f64> = self.iter().map(|v| v.vdist::<f64>(gm)).collect();
        devs.median()    
    }

    /// Initial (first) point for geometric medians.
    fn firstpoint(self) -> Vec<f64> {
        let mut rsum = 0_f64;
        let mut vsum = vec![0_f64; self[0].len()];
        for p in self {
            let mag = p.iter().map(|&pi|f64::from(pi).powi(2)).sum::<f64>(); // vmag();
            if mag.is_normal() {  // skip if p is at the origin
                let rec = 1.0_f64/(mag.sqrt());
                // the sum of reciprocals of magnitudes for the final scaling  
                rsum += rec;
                // so not using simply .unitv 
                vsum.mutvadd::<f64>(&p.smult::<f64>(rec)) // add all unit vectors
            }
        }
        vsum.mutsmult::<f64>(1.0/rsum); // scale by the sum of reciprocals
        vsum
    }    
    
    /// Next approximate gm computed from a member point  
    /// specified by its index `indx` to self. 
    fn nxmember(self, indx: usize) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let p = &self[indx].tof64();
        let mut recip = 0_f64;
        for (i,x) in self.iter().enumerate() {
            if i != indx {  // not point p
                let mag:f64 = x.iter().zip(p).map(|(&xi,&pi)|(f64::from(xi)-pi).powi(2)).sum::<f64>(); 
                if mag.is_normal() { // ignore this point should distance be zero
                    let rec = 1.0_f64/(mag.sqrt());
                    vsum.iter_mut().zip(x).for_each(|(vi,xi)| *vi += rec*f64::from(*xi)); 
                    recip += rec // add separately the reciprocals    
                }
            }
        };
        vsum.iter_mut().for_each(|vi| *vi /= recip);
        vsum
    }
 
    /// Like gmparts, except only does one iteration from any non-member point g
    fn nxnonmember(self, g:&[f64]) -> (Vec<f64>,Vec<f64>,f64) {
        // vsum is the sum vector of unit vectors towards the points
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for x in self { 
            // |x-p| done in-place for speed. Could have simply called x.vdist(p)
            let mag:f64 = x.iter().zip(g).map(|(&xi,&gi)|(f64::from(xi)-gi).powi(2)).sum::<f64>(); 
            if mag.is_normal() { // ignore this point should distance be zero
                let rec = 1.0_f64/(mag.sqrt()); // reciprocal of distance (scalar)
                // vsum increments by components
                vsum.iter_mut().zip(x).for_each(|(vi,xi)| *vi += f64::from(*xi)*rec); 
                recip += rec // add separately the reciprocals for final scaling   
            }
        }
        ( vsum.iter().map(|vi| vi / recip).collect::<Vec<f64>>(),        
          vsum,
          recip )
    }

    /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    /// It has (provably) only vector iterative solutions.
    /// Search methods are slow and difficult in highly dimensional space.
    /// Weiszfeld's fixed point iteration formula has known problems with sometimes failing to converge.
    /// Especially, when the points are dense in the close proximity of the gm, or gm coincides with one of them.  
    /// However, these problems are fixed in my new algorithm here.      
    /// There will eventually be a multithreaded version.
    /// The sum of reciprocals is strictly increasing and so is used here as
    /// easy to evaluate termination condition.
    fn gmedian(self, eps: f64) -> Vec<f64> {  
        let mut g = self.acentroid(); // start iterating from the mean  or vec![0_f64; self[0].len()];
        let mut recsum = 0f64;
        loop { // vector iteration till accuracy eps is exceeded  
            let mut nextg = vec![0_f64; self[0].len()];   
            let mut nextrecsum = 0_f64;
            for v in self {   
                // |v-g| done in-place for speed. Could have simply called x.vdist(g)
                let mag:f64 = v.iter().zip(&g).map(|(&vi,gi)|(f64::from(vi)-gi).powi(2)).sum(); 
                if mag > eps { 
                    let rec = 1.0_f64/(mag.sqrt()); // reciprocal of distance (scalar)
                    // vsum increment by components
                    for (vi,gi) in v.iter().zip(&mut nextg) { *gi += f64::from(*vi)*rec }; 
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

    /// Point by point Geometric Median (gm).
    /// Gm is the point that minimises the sum of distances to a given set of points.
    /// It has (provably) only vector iterative solutions.
    /// Search methods are slow and difficult in highly dimensional space.
    /// Weiszfeld's fixed point iteration formula has known problems with sometimes failing to converge.
    /// Especially, when the points are dense in the close proximity of the gm, or gm coincides with one of them.  
    /// However, these problems are fixed in my new algorithm here.      
    /// There will eventually be a multithreaded version.
    /// The sum of reciprocals is strictly increasing and so is used here as
    /// easy to evaluate termination condition.
    fn pmedian(self, eps: f64) -> Vec<f64> {  
        // start iterating from the centroid, alternatively from the origin: vec![0_f64; self[0].len()] 
        let mut g = self.acentroid();
        // running global sum of reciprocals
        let mut recsum = 0f64;     
        // running global sum of unit vectors
        let mut vsum = vec![0_f64; self[0].len()];
        // previous reciprocals for each point 
        let mut precs:Vec<f64> = Vec::with_capacity(self.len());
        // termination flag triggered by any one point
        let mut terminate = true; 
        
        // initial vsum,recsum and precs
        for p in self { 
            let magsq:f64 = p.iter().zip(&g).map(|(&pi,gi)|(f64::from(pi)-gi).powi(2)).sum(); 
            if magsq < eps { precs.push(0.); continue; }; // skip this point, it is too close 
            let rec = 1.0/(magsq.sqrt()); 
            // vsum incremented by components of unit vector
            for (vscomp,&pcomp) in vsum.iter_mut().zip(p) { *vscomp += rec*f64::from(pcomp) }; 
            // vsum.mutvadd::<f64>(&p.smult::<f64>(rec)); // the above, shorter but slower
            precs.push(rec); // store rec for this p
            recsum += rec; 
        }
        // first iteration done, update g
        for (gcomp,vscomp) in g.iter_mut().zip(&vsum) { *gcomp = vscomp/recsum }; 
        g = vsum.smult::<f64>(1.0/recsum); 
        loop { // vector iteration till accuracy eps is exceeded 
            for (p,rec) in self.iter().zip(&mut precs) { 
                let magsq:f64 = p.iter().zip(&g).map(|(&pi,gi)|(f64::from(pi)-gi).powi(2)).sum(); 
                if magsq < eps { *rec = 0.0; continue; }; // skip this point, it is too close
                let recip = 1.0/(magsq.sqrt()); 
                let recdelta = recip - *rec; // change in reciprocal for p
                *rec = recip; // update rec for this p for next time
                // vsum updated by components
                for (vscomp,pcomp) in vsum.iter_mut().zip(p) { *vscomp += recdelta*f64::from(*pcomp) };
                // update recsum
                recsum += recdelta; 
                // update g immediately for each point p 
                for (gcomp,vscomp) in g.iter_mut().zip(&vsum) { *gcomp = vscomp/recsum };
                // termination condition detected but do the rest of the points anyway
                if terminate && recdelta.abs() > eps { terminate = false }; 
            }
            if terminate { return g };  // termination reached 
            terminate = true
        }
    }

    /// Like `gmedian` but returns also the sum of unit vecs and the sum of reciprocals. 
    fn gmparts(self, eps: f64) -> (Vec<f64>,Vec<f64>,f64) { 
        let mut g = self.acentroid(); // start iterating from the Centre
        let mut recsum = 0f64; 
        loop { // vector iteration till accuracy eps is exceeded  
            let mut nextg = vec![0_f64; self[0].len()];   
            let mut nextrecsum = 0f64;
            for x in self { // for all points
                // |x-g| done in-place for speed. Could have simply called x.vdist(g)
                //let mag:f64 = g.vdist::<f64>(&x); 
                let mag = g.iter().zip(x).map(|(&gi,&xi)|(f64::from(xi)-gi).powi(2)).sum::<f64>(); 
                if mag.is_normal() { 
                    let rec = 1.0_f64/(mag.sqrt()); // reciprocal of distance (scalar)
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

}
