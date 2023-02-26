use crate::{sumn, MStats, MinMax, MutVecg, RError, Stats, TriangMat, VecVec, Vecg, RE};
use indxvec::{Mutops, Vecops};
use medians::{MedError, Medianf64};
use rayon::prelude::*;

impl<T> VecVec<T> for &[Vec<T>]
where
    T: Copy + PartialOrd + Sync,
    Vec<Vec<T>>: IntoParallelIterator,
    Vec<T>: IntoParallelIterator,
    f64: From<T>,
{
    /// Selects a column by number
    fn column(self, cnum: usize) -> Vec<f64> {
        self.iter().map(|row| f64::from(row[cnum])).collect()
    }

    /// Multithreaded transpose of vec of vecs matrix
    fn transpose(self) -> Vec<Vec<f64>> {
        (0..self[0].len())
            .into_par_iter()
            .map(|cnum| self.column(cnum))
            .collect()
    }

    /// Normalize columns, so that they become unit row vectors
    fn normalize(self) -> Result<Vec<Vec<f64>>, RE> {
        (0..self[0].len())
            .into_par_iter()
            .map(|cnum| -> Result<Vec<f64>, RE> { self.column(cnum).vunit() })
            .collect()
    }

    /// Householder's method returning triangular matrices (U',R), where
    /// U are the reflector generators for use by house_uapply(m).
    /// R is the upper triangular decomposition factor.
    /// Here both U and R are returned for convenience in their transposed lower triangular forms.
    /// Transposed input self for convenience, so that original columns get accessed easily as rows.
    fn house_ur(self) -> Result<(TriangMat, TriangMat),RE> {
        let n = self.len();
        let d = self[0].len();
        let min = if d <= n { d } else { n }; // minimal dimension
        let mut r = self.transpose(); // self.iter().map(|s| s.tof64()).collect::<Vec<Vec<f64>>>(); //  //
        let mut ures = vec![0.; sumn(min)];
        let mut rres = Vec::with_capacity(sumn(min));
        for j in 0..min {
            let Some(slc) = r[j].get(j..d) 
            else { return Err(RError::DataError("house_ur: failed to extract uvec slice".to_owned()));}; 
            let uvec = slc.house_reflector();
            for rlast in r.iter_mut().take(d).skip(j) {
                let rvec = uvec.house_reflect::<f64>(&rlast.drain(j..d).collect::<Vec<f64>>());
                rlast.extend(rvec);
                // drained, reflected with this uvec, and rebuilt, all remaining rows of r
            }
            // these uvecs are columns, so they must saved column-wise
            for (row, &usave) in uvec.iter().enumerate() {
                ures[sumn(row + j) + j] = usave; // using triangular index
            }
            // save completed `r[j]` components only up to and including the diagonal
            // we are not even storing the rest, so no need to set those to zero
            for &rsave in r[j].iter().take(j + 1) {
                rres.push(rsave)
            }
        }
        Ok((
            TriangMat {
                kind: 3,
                data: ures,
            }, // transposed, non symmetric kind
            TriangMat {
                kind: 3,
                data: rres,
            }, // transposed, non symmetric kind
        ))
    }

    /// Joint probability density function of n matched slices of the same length
    fn jointpdfn(self) -> Result<Vec<f64>, RE> {
        let d = self[0].len(); // their common dimensionality (length)
        for v in self.iter().skip(1) {
            if v.len() != d {
                return Err(RError::DataError(
                    "jointpdfn: all vectors must be of equal length!".to_owned(),
                ));
            };
        }
        let mut res: Vec<f64> = Vec::with_capacity(d);
        let mut tuples = self.transpose();
        let df = tuples.len() as f64; // for turning counts to probabilities
                                      // lexical sort to group together occurrences of identical tuples
        tuples.sort_unstable_by(|a, b| {
            let Some(x) = a.partial_cmp(b) 
            else { panic!("jointpdfn: comparison fail in f64 sort!"); }; x});
        let mut count = 1_usize; // running count
        let mut lastindex = 0; // initial index of the last unique tuple
        tuples.iter().enumerate().skip(1).for_each(|(i, ti)| {
            if ti > &tuples[lastindex] {
                // new tuple ti (Vec<T>) encountered
                res.push((count as f64) / df); // save frequency count as probability
                lastindex = i; // current index becomes the new one
                count = 1_usize; // reset counter
            } else {
                count += 1;
            }
        });
        res.push((count as f64) / df); // flush the rest!
        Ok(res)
    }

    /// Joint entropy of vectors of the same length
    fn jointentropyn(self) -> Result<f64, RE> {
        let jpdf = self.jointpdfn()?;
        Ok(jpdf.iter().map(|&x| -x * (x.ln())).sum())
    }

    /// Dependence (component wise) of a set of vectors.
    /// i.e. `dependencen` returns 0 iff they are statistically independent
    /// bigger values when they are dependentent
    fn dependencen(self) -> Result<f64, RE> {
        Ok((0..self.len())
            .into_par_iter()
            .map(|i| self[i].entropy())
            .sum::<f64>()
            / self.jointentropyn()?
            - 1.0)
    }

    /// Flattened lower triangular part of a symmetric matrix for vectors in self.
    /// The upper triangular part can be trivially generated for all j>i by: c(j,i) = c(i,j).
    /// Applies closure f to compute a scalar binary relation between all pairs of vector
    /// components of self.   
    /// The closure typically invokes one of the methods from Vecg trait (in vecg.rs),
    /// such as dependencies or correlations.  
    /// Example call: `pts.transpose().crossfeatures(|v1,v2| v1.mediancorrf64(v2)?)?`
    /// computes median correlations between all column vectors (features) in pts.
    fn crossfeatures(self, f: fn(&[T], &[T]) -> f64) -> Result<TriangMat, RE> {
        Ok(TriangMat {
            kind: 2, // symmetric, non transposed
            data: (0..self.len())
                .into_par_iter()
                .flat_map(|i| {
                    (0..i + 1usize)
                        .into_iter()
                        .map(|j| f(&self[i], &self[j]))
                        .collect::<Vec<f64>>()
                })
                .collect::<Vec<f64>>(),
        })
    }

    /// Sum of nd points (or vectors)
    fn sumv(self) -> Vec<f64> {
        let mut resvec = vec![0_f64; self[0].len()];
        for v in self {
            resvec.mutvadd(v)
        }
        resvec
    }

    /// acentroid = multidimensional arithmetic mean
    fn acentroid(self) -> Vec<f64> {
        self.sumv().smult::<f64>(1. / (self.len() as f64))
    }

    /// multithreaded acentroid = multidimensional arithmetic mean
    fn par_acentroid(self) -> Vec<f64> {
        let sumvec = self
            .par_iter()
            .fold(
                || vec![0_f64; self[0].len()],
                |mut vecsum: Vec<f64>, p| {
                    vecsum.mutvadd(p);
                    vecsum
                },
            )
            .reduce(
                || vec![0_f64; self[0].len()],
                |mut finalsum: Vec<f64>, partsum: Vec<f64>| {
                    finalsum.mutvadd::<f64>(&partsum);
                    finalsum
                },
            );
        sumvec.smult::<f64>(1. / (self.len() as f64))
    }

    /// gcentroid = multidimensional geometric mean
    fn gcentroid(self) -> Result<Vec<f64>,RE> { 
        let logvs = self.iter().map(|v|-> Result<Vec<f64>,RE> {
            Ok(v.vunit()?.smult::<f64>(v.vmag().ln())) })
            .collect::<Result<Vec<Vec<f64>>,RE>>()?;
        let logcentroid = logvs.acentroid();
        Ok(logcentroid.vunit()?.smult::<f64>(logcentroid.vmag().exp()))
    }

    /// hcentroid =  multidimensional harmonic mean
    fn hcentroid(self) -> Result<Vec<f64>,RE> {
        let mut centre = vec![0_f64; self[0].len()];
        for v in self {
            centre.mutvadd::<f64>(&v.vinverse()?)
        }
        Ok(centre  
            .vinverse()?.smult::<f64>(self.len() as f64)) 
    }

    /// For each member point, gives its sum of distances to all other points and their MinMax
    fn distsums(self) -> Vec<f64> {
        let n = self.len();
        let mut dists = vec![0_f64; n]; // distances accumulator for all points
                                        // examine all unique pairings (lower triangular part of symmetric flat matrix)
        self.iter().enumerate().for_each(|(i, thisp)| {
            self.iter().take(i).enumerate().for_each(|(j, thatp)| {
                let d = thisp.vdist(thatp); // calculate each distance relation just once
                dists[i] += d;
                dists[j] += d; // but add it to both points' sums
            })
        });
        dists
    }

    /// Medoid and Outlier (Medout)
    /// Medoid is the member point (point belonging to the set of points `self`),
    /// which has the least sum of distances to all other points.
    /// Outlier is the point with the greatest sum of distances.
    /// In other words, they are the members nearest and furthest from the geometric median.
    /// Returns struct MinMax{min,minindex,max,maxindex}
    fn medout(self, gm: &[f64]) -> MinMax<f64> {
        self.par_iter()
            .map(|s| s.vdist::<f64>(gm))
            .collect::<Vec<f64>>()
            .minmax()
    }

    /// Radius of a point specified by its subscript.    
    fn radius(self, i: usize, gm: &[f64]) -> Result<f64, RE> {
        if i > self.len() {
            return Err(RError::DataError("radius: invalid subscript".to_owned()));
        }
        Ok(self[i].vdist::<f64>(gm))
    }

    /// Exact radii (eccentricity) magnitudes for all member points from the Geometric Median.
    /// More accurate and usually faster as well than the approximate `eccentricities` above,
    /// especially when there are many points.
    fn radii(self, gm: &[f64]) -> Vec<f64> {
        self.iter()
            .map(|s: &Vec<T>| gm.vdist(s))
            .collect::<Vec<f64>>()
    }

    /// Arith mean and std (in MStats struct), Median and MAD (in another MStats struct), Medoid and Outlier (in MinMax struct)
    /// of scalar radii of points in self.
    /// These are new robust measures of a cloud of multidimensional points (or multivariate sample).  
    fn eccinfo(self, gm: &[f64]) -> Result<(MStats, MStats, MinMax<f64>), RE>
    where
        Vec<f64>: FromIterator<f64>,
    {
        let rads: Vec<f64> = self.radii(gm);
        Ok((rads.ameanstd()?, rads.medstatsf64()?, rads.minmax()))
    }

    /// Quasi median, recommended only for comparison purposes
    fn quasimedian(self) -> Result<Vec<f64>, RE> {
        Ok((0..self[0].len())
            .into_iter()
            .map(|colnum| self.column(colnum).medianf64())
            .collect::<Result<Vec<f64>, MedError<String>>>()?)
    }

    /// Geometric median's estimated error
    fn gmerror(self, g: &[f64]) -> f64 {
        let (gm, _, _) = self.nxnonmember(g);
        gm.vdist::<f64>(g)
    }

    /// Proportions of points along each +/-axis (hemisphere)
    /// Excludes points that are perpendicular to axis
    /// Uses only the selected points specified in idx (e.g. the hull).
    /// Self should normally be zero mean/median vectors,
    /// e.g. `self.translate(&median)`
    fn tukeyvec(self, idx: &[usize]) -> Result<Vec<f64>, RE> {
        let dims = self[0].len();
        if self.is_empty() {
            return Err(RError::NoDataError("tukeyvec given no data".to_owned()));
        };
        let mut hemis = vec![0_f64; 2 * dims];
        for &i in idx {
            for (j, &component) in self[i].iter().enumerate() {
                let cf = f64::from(component);
                if cf > 0. {
                    hemis[j] += 1.
                } else if cf < 0. {
                    hemis[dims + j] += 1.
                };
            }
        }
        hemis
            .iter_mut()
            .for_each(|count| *count /= idx.len() as f64);
        Ok(hemis)
    }

    /// madgm median of distances from gm: stable nd data spread measure
    fn madgm(self, gm: &[f64]) -> Result<f64, RE> {
        if self.is_empty() { 
            return Err(RError::NoDataError("stdgm given zero length vec!".to_owned())); };     
        Ok( self
            .iter()
            .map(|p| p.vdist::<f64>(gm))
            .collect::<Vec<f64>>()
            .medianf64()?
        )
    }

    /// stdgm mean of distances from gm: nd data spread measure, aka nd standard deviation
    fn stdgm(self, gm: &[f64]) -> Result<f64,RE> { 
        if self.is_empty() { 
            return Err(RError::NoDataError("stdgm given zero length vec!".to_owned())); };     
        Ok( self
            .iter()
            .map(|p| p.vdist::<f64>(gm))
            .sum::<f64>()/self.len() as f64
        ) 
}

    /// Collects indices of inner (or core) hull and outer hull, from zero median points in self.
    /// Defining plane of a point A goes through A and is normal to the zero median vector **a**.      
    /// B is an inner hull point, when it lies inside all other points' defining planes.  
    /// B is an outer hull point, when there is no other point beyond its own defining plane.
    /// B can belong to both hulls, as when all the points lie on a hyper-sphere around gm.   
    /// The testing is done in increasing (decreasing) radius order.  
    /// B lies outside the defining plane of **a**, when its projection onto unit **a** exceeds `|a|`:    
    /// `|b|cos(θ) > |a| => a*b > |a|^2`,  
    /// such B immediately fails as a candidate for the inner hull.
    /// Working with square magnitudes, `|a|^2` saves taking square roots and dividing the dot product by |a|.  
    /// Similarly for the outer hull, where A and B simply swap roles.
    fn hulls(self) -> (Vec<usize>, Vec<usize>) {
        let sqradii = self.par_iter().map(|s| s.vmagsq()).collect::<Vec<f64>>();
        let mut radindex = sqradii.hashsort_indexed(&mut |x| *x); // ascending square radii
        let innerindex = radindex
            .par_iter()
            .filter_map(|&b| {
                // test all points in ascending order
                for &a in &radindex {
                    // check against all points 'a' up to 'b'
                    if a == b {
                        return Some(b);
                    }; // b passed
                       // b lies outside of a => immediately reject b
                    if self[a].dotp(&self[b]) > sqradii[a] {
                        break;
                    };
                }
                None
            })
            .collect::<Vec<usize>>();
        radindex.mutrevs(); // make the order of points descending
        let mut outerindex = radindex
            .par_iter()
            .filter_map(|&b| {
                // test all points, in descending order
                for &a in &radindex {
                    // a can only be outside of b for a's of greater magnitude
                    if a == b {
                        return Some(b);
                    }; // b passed
                    let dotp = self[a].dotp(&self[b]);
                    // a lies outside of b => immediately reject b
                    if dotp > sqradii[b] {
                        break;
                    };
                }
                None
            })
            .collect::<Vec<usize>>();
        outerindex.reverse();
        (innerindex, outerindex)
    }

    /// Initial (first) point for geometric medians.
    fn firstpoint(self) -> Vec<f64> {
        let mut rsum = 0_f64;
        let mut vsum = vec![0_f64; self[0].len()];
        for p in self {
            let mag = p.iter().map(|&pi| f64::from(pi).powi(2)).sum::<f64>(); // vmag();
            if mag.is_normal() {
                // skip if p is at the origin
                let rec = 1.0_f64 / (mag.sqrt());
                // the sum of reciprocals of magnitudes for the final scaling
                rsum += rec;
                // add this unit vector to their sum
                vsum.mutvadd::<f64>(&p.smult::<f64>(rec))
            }
        }
        vsum.mutsmult::<f64>(1.0 / rsum); // scale by the sum of reciprocals
        vsum // good initial gm
    }

    /// Like gmparts, except only does one iteration from any non-member point g
    fn nxnonmember(self, g: &[f64]) -> (Vec<f64>, Vec<f64>, f64) {
        // vsum is the sum vector of unit vectors towards the points
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for x in self {
            // |x-p| done in-place for speed. Could have simply called x.vdist(p)
            let mag: f64 = x
                .iter()
                .zip(g)
                .map(|(&xi, &gi)| (f64::from(xi) - gi).powi(2))
                .sum::<f64>();
            if mag.is_normal() {
                // ignore this point should distance be zero
                let rec = 1.0_f64 / (mag.sqrt()); // reciprocal of distance (scalar)
                                                  // vsum increments by components
                vsum.iter_mut()
                    .zip(x)
                    .for_each(|(vi, xi)| *vi += f64::from(*xi) * rec);
                recip += rec // add separately the reciprocals for final scaling
            }
        }
        (
            vsum.iter().map(|vi| vi / recip).collect::<Vec<f64>>(),
            vsum,
            recip,
        )
    }

    /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
    /// It has (provably) only vector iterative solutions.
    /// Search methods are slow and difficult in highly dimensional space.
    /// Weiszfeld's fixed point iteration formula has known problems with sometimes failing to converge.
    /// Especially, when the points are dense in the close proximity of the gm, or gm coincides with one of them.  
    /// However, these problems are fixed in my new algorithm here.
    /// The sum of reciprocals is strictly increasing and so is used here as
    /// easy to evaluate termination condition.
    fn gmedian(self, eps: f64) -> Vec<f64> {
        let mut g = self.acentroid(); // start iterating from the mean  or vec![0_f64; self[0].len()];
        let mut recsum = 0f64;
        loop {
            // vector iteration till accuracy eps is exceeded
            let mut nextg = vec![0_f64; self[0].len()];
            let mut nextrecsum = 0_f64;
            for p in self {
                // |p-g| done in-place for speed. Could have simply called p.vdist(g)
                let mag: f64 = p
                    .iter()
                    .zip(&g)
                    .map(|(&vi, gi)| (f64::from(vi) - gi).powi(2))
                    .sum();
                if mag > eps {
                    let rec = 1.0_f64 / (mag.sqrt()); // reciprocal of distance (scalar)
                                                      // vsum increment by components
                    for (vi, gi) in p.iter().zip(&mut nextg) {
                        *gi += f64::from(*vi) * rec
                    }
                    nextrecsum += rec // add separately the reciprocals for final scaling
                } // else simply ignore this point v, should its distance from g be <= eps
            }
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);
            // eprintln!("recsum {}, nextrecsum {} diff {}",recsum,nextrecsum,nextrecsum-recsum);
            if nextrecsum - recsum < eps {
                return nextg;
            }; // termination test
            g = nextg;
            recsum = nextrecsum;
        }
    }

    /// Parallel (multithreaded) implementation of Geometric Median. Possibly the fastest you will find.  
    /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.  
    /// It has (provably) only vector iterative solutions.    
    /// Search methods are slow and difficult in hyper space.    
    /// Weiszfeld's fixed point iteration formula has known problems and sometimes fails to converge.  
    /// Specifically, when the points are dense in the close proximity of the gm, or gm coincides with one of them.    
    /// However, these problems are solved in my new algorithm here.     
    /// The sum of reciprocals is strictly increasing and so is used to easily evaluate the termination condition.  
    fn par_gmedian(self, eps: f64) -> Vec<f64> {
        let mut g = self.par_acentroid(); // start iterating from the mean  or vec![0_f64; self[0].len()];
        let mut recsum = 0_f64;
        loop {
            // vector iteration till accuracy eps is exceeded
            let (mut nextg, nextrecsum) = self
                .par_iter()
                .fold(
                    || (vec![0_f64; self[0].len()], 0_f64),
                    |mut pair: (Vec<f64>, f64), p: &Vec<T>| {
                        // |p-g| done in-place for speed. Could have simply called p.vdist(g)
                        let mag: f64 = p
                            .iter()
                            .zip(&g)
                            .map(|(&vi, gi)| (f64::from(vi) - gi).powi(2))
                            .sum();
                        // let (mut vecsum, mut recsum) = pair;
                        if mag > eps {
                            let rec = 1.0_f64 / (mag.sqrt()); // reciprocal of distance (scalar)
                            for (vi, gi) in p.iter().zip(&mut pair.0) {
                                *gi += f64::from(*vi) * rec
                            }
                            pair.1 += rec; // add separately the reciprocals for the final scaling
                        } // else simply ignore this point should its distance from g be zero
                        pair
                    },
                )
                // must run reduce on the partial sums produced by fold
                .reduce(
                    || (vec![0_f64; self[0].len()], 0_f64),
                    |mut pairsum: (Vec<f64>, f64), pairin: (Vec<f64>, f64)| {
                        pairsum.0.mutvadd::<f64>(&pairin.0);
                        pairsum.1 += pairin.1;
                        pairsum
                    },
                );
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);
            if nextrecsum - recsum < eps {
                return nextg;
            }; // termination test
            g = nextg;
            recsum = nextrecsum;
        }
    }

    /// Like `gmedian` but returns also the sum of unit vecs and the sum of reciprocals.
    fn gmparts(self, eps: f64) -> (Vec<f64>, Vec<f64>, f64) {
        let mut g = self.acentroid(); // start iterating from the Centre
        let mut recsum = 0f64;
        loop {
            // vector iteration till accuracy eps is exceeded
            let mut nextg = vec![0_f64; self[0].len()];
            let mut nextrecsum = 0f64;
            for x in self {
                // for all points
                // |x-g| done in-place for speed. Could have simply called x.vdist(g)
                //let mag:f64 = g.vdist::<f64>(&x);
                let mag = g
                    .iter()
                    .zip(x)
                    .map(|(&gi, &xi)| (f64::from(xi) - gi).powi(2))
                    .sum::<f64>();
                if mag.is_normal() {
                    let rec = 1.0_f64 / (mag.sqrt()); // reciprocal of distance (scalar)
                                                      // vsum increments by components
                    nextg
                        .iter_mut()
                        .zip(x)
                        .for_each(|(vi, &xi)| *vi += f64::from(xi) * rec);
                    nextrecsum += rec // add separately the reciprocals for final scaling
                } // else simply ignore this point should its distance from g be zero
            }
            if nextrecsum - recsum < eps {
                return (
                    nextg
                        .iter()
                        .map(|&gi| gi / nextrecsum)
                        .collect::<Vec<f64>>(),
                    nextg,
                    nextrecsum,
                );
            }; // termination
            nextg.iter_mut().for_each(|gi| *gi /= nextrecsum);
            g = nextg;
            recsum = nextrecsum;
        }
    }
}
