// use core::slice::SlicePattern;

use crate::{
    error::{re_error, RError, RE},
    fromop, Stats, TriangMat, Vecg
};
use indxvec::{Indices, Vecops};
use medians::Median;

impl<T> Vecg for &[T]
where
    T: Clone + PartialOrd + Into<f64>
{
    /// nd t_statistic of self against geometric median and madgm spread.     
    /// Unlike in 1d, is always positive.
    fn t_statistic(self, gm: &[f64], madgm: f64) -> Result<f64, RE> {
        Ok(self.vdist::<f64>(gm) / madgm)
    }

    /// Dot product of vector self with column c of matrix v
    fn columnp<U: Clone + Into<f64>>(self, c: usize, v: &[Vec<U>]) -> f64 {
        self.iter()
            .enumerate()
            .map(|(i, x)| x.clone().into() * v[i][c].clone().into())
            .sum::<f64>()
    }

    /// Scalar addition to a vector, creates new vec
    fn sadd<U: Into<f64>>(self, s: U) -> Vec<f64> {
        let sf: f64 = s.into();
        self.iter().map(|x| sf + x.clone().into()).collect()
    }

    /// Scalar multiplication with a vector, creates new vec
    fn smult<U: Into<f64>>(self, s: U) -> Vec<f64> {
        let sf: f64 = s.into();
        self.iter().map(|x| sf * x.clone().into()).collect()
    }

    /// Scalar product.   
    /// Must be of the same length - no error checking (for speed)
    fn dotp<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(xi, vi)| -> f64 { xi.clone().into() * vi.clone().into() })
            .sum::<f64>()
    }

    /// Product with signature vec of hemispheric frequencies.  
    /// Similar result can be obtained
    /// by projecting onto self all points but that is usually too slow.
    fn dotsig(self, sig: &[f64]) -> Result<f64, RE> {
        let dims = self.len();
        if 2 * dims != sig.len() {
            return Err(re_error(
                "size",
                "dotsig: sig vec must have double the dimensions",
            ));
        }
        let sigunit = sig.vunit()?;
        let sunit = self.vunit()?;
        let mut ressum = 0_f64;
        for (i, &scomp) in sunit.iter().enumerate() {
            if scomp > 0_f64 {
                ressum += scomp * sigunit[i];
                continue;
            };
            if scomp < 0_f64 {
                ressum -= scomp * sig[dims + i];
            };
        }
        Ok(ressum)
    }

    /// Cosine of angle between the two slices.
    /// Done in one iteration for efficiency.
    fn cosine<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        let (mut sxy, mut sy2) = (0_f64, 0_f64);
        let sx2: f64 = self
            .iter()
            .zip(v)
            .map(|(tx, uy)| {
                let x = tx.clone().into();
                let y = uy.clone().into();
                sxy += x * y;
                sy2 += y * y;
                x * x
            })
            .sum();
        sxy / (sx2 * sy2).sqrt()
    }

    /// Generalized cross product:  
    /// Sine of an angle between **self** and **v** with correct sign in any number of dimensions
    fn sine<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        let blade = self.wedge(v);
        blade.sum().signum()
            * (blade.magsq()
                / (self.vmagsq() * v.iter().map(|x| x.clone().into().powi(2)).sum::<f64>()))
            .sqrt()
    }

    /// Vector subtraction
    fn vsub<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<f64> {
        self.iter()
            .zip(v)
            .map(|(xi, vi)| xi.clone().into() - vi.clone().into())
            .collect()
    }

    /// Vectors difference unitised (done together for efficiency)
    fn vsubunit<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<f64> {
        let mut sumsq = 0_f64;
        let dif = self
            .iter()
            .zip(v)
            .map(|(xi, vi)| {
                let d = xi.clone().into() - vi.clone().into();
                sumsq += d * d;
                d
            })
            .collect::<Vec<f64>>();
        dif.smult::<f64>(1_f64 / sumsq.sqrt())
    }

    /// Vector addition
    fn vadd<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<f64> {
        self.iter()
            .zip(v)
            .map(|(xi, vi)| xi.clone().into() + vi.clone().into())
            .collect()
    }

    /// Euclidian distance   
    fn vdist<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(xi, vi)| (xi.clone().into() - vi.clone().into()).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Weighted arithmetic mean of `self:&[T]`, scaled by `ws:&[U]`
    fn wvmean<U: Clone + Into<f64>>(self, ws: &[U]) -> f64 {
        let mut wsum: f64 = 0.;
        let mut sum: f64 = 0.;
        for (s, w) in self.iter().zip(ws) {
            let fw = w.clone().into();
            sum += fw * s.clone().into();
            wsum += fw;
        }
        sum / wsum
    }

    /// Weighted distance of `self:&[T]` to `v:&[V]`, scaled by `ws:&[U]`
    /// allows all three to be of different types
    fn wvdist<U, V>(self, ws: &[U], v: &[V]) -> f64
    where
        U: Clone + Into<f64>,
        V: Clone + Into<f64>,
    {
        self.iter()
            .enumerate()
            .map(|(i, xi)| (ws[i].clone().into() * xi.clone().into() - v[i].clone().into()).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Euclidian distance squared  
    fn vdistsq<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(xi, vi)| (xi.clone().into() - vi.clone().into()).powi(2))
            .sum::<f64>()
    }

    /// cityblock distance
    fn cityblockd<U: Clone + Into<f64>>(self, v: &[U]) -> f64
    where
        U: Into<f64>,
    {
        self.iter()
            .zip(v)
            .map(|(xi, vi)| (xi.clone().into() - vi.clone().into()).abs())
            .sum::<f64>()
    }

    /// Area spanned by two vectors over their concave angle (always >= 0)
    /// |a||b||sin(theta)| == (1-cos2(theta)).sqrt()
    /// Attains maximum `|a|.|b|` when the vectors are orthogonal.
    fn varea<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> f64 {
        (self.vmagsq() * v.vmagsq() - self.dotp(v).powi(2)).sqrt()
    }

    /// Area of swept arc
    /// = |a||b|(1-cos(theta)) = 2|a||b|D
    fn varc<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> f64 {
        (v.vmagsq() * self.vmagsq()).sqrt() - self.dotp(v)
    }

    /// Positive dotp in the interval: `[0,2|a||b|]`
    /// = |a||b|(1+cos(theta)) = 2|a||b|S
    fn pdotp<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> f64 {
        (v.vmagsq() * self.vmagsq()).sqrt() + self.dotp(v)
    }

    /// We define vector similarity S in the interval [0,1] as
    /// S = (1+cos(theta))/2
    fn vsim<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        (1.0 + self.cosine(v)) / 2.0
    }

    /// We define vector dissimilarity D in the interval [0,1] as
    /// D = 1-S = (1-cos(theta))/2
    fn vdisim<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        (1.0 - self.cosine(v)) / 2.0
    }

    /// We define vector median correlation similarity in the interval [0,1] as
    fn vcorrsim(self, v:Self) -> f64 {
        (1.0 + self.mediancorr(v, fromop).expect("vcorrsim: mediancorr failed")) / 2.0
    }

    /// Lower triangular covariance matrix for a single vector.
    /// Where m is either mean or median vector (to be subtracted).
    /// Covariance matrix is symmetric (kind:2) (positive semi definite).
    fn covone<U: Clone + Into<f64>>(self, m: &[U]) -> TriangMat {
        let mut cov: Vec<f64> = Vec::new(); // flat lower triangular result array
        let vm = self.vsub(m); // zero mean/median vector
        vm.iter().enumerate().for_each(|(i,&thisc)|
            // generate its products up to and including the diagonal (itself)
            vm.iter().take(i+1).for_each(|&component| cov.push(thisc*component)));
        TriangMat { kind: 2, data: cov }
    }

    /// Kronecker product of two vectors.   
    /// The indexing is always assumed to be in this order: row,column.
    /// Flat version of outer(wedge) product
    fn kron<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<f64> {
        let vf = v.iter().map(|vi| vi.clone().into()).collect::<Vec<f64>>();
        self.iter()
            .flat_map(|s| {
                let sf: f64 = s.clone().into();
                vf.iter().map(move |&vfi| sf * vfi)
            })
            .collect::<Vec<f64>>()
    }

    /// Outer product: matrix multiplication of column self with row v.
    fn outer<U: Clone + Into<f64>>(self, v: &[U]) -> Vec<Vec<f64>> {
        let vf = v.iter().map(|vi| vi.clone().into()).collect::<Vec<f64>>();
        self.iter()
            .map(|s| {
                let sf: f64 = s.clone().into();
                vf.iter().map(|&vfi| sf * vfi).collect::<Vec<f64>>()
            })
            .collect::<Vec<Vec<f64>>>()
    }

    /// Exterior (Grassman) algebra product: produces components of 2-blade **a∧b**
    fn wedge<U: Clone + Into<f64>>(self, b: &[U]) -> TriangMat {
        let n = self.len();
        assert_eq!(n, b.len());
        let mut result: Vec<f64> = Vec::new();
        for i in 0..n {
            let ai: f64 = self[i].clone().into();
            let bi: f64 = b[i].clone().into();
            for j in 0..i + 1 {
                result.push(ai * b[j].clone().into() - bi * self[j].clone().into());
            }
        }
        TriangMat {
            kind: 1,
            data: result,
        }
    }

    /// Geometric (Clifford) algebra product: **a*b** + **a∧b** in matrix form
    /// here the elements of the dot product a*b are placed in their
    /// natural positions on the diagonal (can be easily added)
    fn geometric<U: Clone + Into<f64>>(self, b: &[U]) -> TriangMat {
        let n = self.len();
        assert_eq!(n, b.len());
        let mut result: Vec<f64> = Vec::new();
        for i in 0..n {
            let ai: f64 = self[i].clone().into();
            let bi: f64 = b[i].clone().into();
            for j in 0..i {
                result.push(ai * b[j].clone().into() - bi * self[j].clone().into());
            }
            result.push(ai * bi); // the diagonal dot product element
        }
        TriangMat {
            kind: 1,
            data: result,
        }
    }

    /// Joint probability density function of two pairwise matched slices
    fn jointpdf<U: Clone + Into<f64>>(self, v: &[U]) -> Result<Vec<f64>, RE> {
        let n = self.len();
        if v.len() != n {
            return Err(RError::DataError(
                "{jointpdf argument vectors must be of equal length!".to_owned(),
            ));
        };
        let nf = n as f64;
        let mut res: Vec<f64> = Vec::new();
        // collect successive pairs, upgrading all end types to common f64
        let mut spairs: Vec<(f64, f64)> = self
            .iter()
            .zip(v)
            .map(|(si, vi)| (si.clone().into(), vi.clone().into()))
            .collect();
        // sort them to group all same pairs together for counting
        spairs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let mut count = 1_usize; // running count
        let mut lastindex = 0;
        spairs.iter().enumerate().skip(1).for_each(|(i, si)| {
            if si > &spairs[lastindex] {
                // new pair encountered
                res.push((count as f64) / nf); // save previous probability
                lastindex = i; // current index becomes the new one
                count = 1_usize; // reset counter
            } else {
                count += 1;
            }
        });
        res.push((count as f64) / nf); // flush the rest!
        Ok(res)
    }

    /// Joint entropy of two sets of the same length
    fn jointentropy<U: Clone + Into<f64>>(self, v: &[U]) -> Result<f64, RE> {
        let jpdf = self.jointpdf(v)?;
        Ok(jpdf.iter().map(|&x| -x * (x.ln())).sum())
    }

    /// Dependence of &[T] &[U] variables in the range [0,1]
    /// returns 0 iff they are statistically component wise independent
    /// returns 1 when they are identical or all their values are unique
    fn dependence<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> Result<f64, RE> {
        Ok((self.entropy() + v.entropy()) / self.jointentropy(v)? - 1.0)
    }

    /// Independence of &[T] &[U] variables in the range [0,1]
    /// returns 1 iff they are statistically component wise independent
    fn independence<U: Clone + PartialOrd + Into<f64>>(self, v: &[U]) -> Result<f64, RE> {
        Ok(2.0 * self.jointentropy(v)? / (self.entropy() + v.entropy()) - 1.0)
    }

    /// Pearson's (most common) correlation.
    /// # Example
    /// ```
    /// use rstats::Vecg;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.correlation(&v2),-0.1076923076923077);
    /// ```
    fn correlation<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let sx: f64 = self
            .iter()
            .zip(v)
            .map(|(xt, yu)| {
                let x = xt.clone().into();
                let y = yu.clone().into();
                sy += y;
                sxy += x * y;
                sx2 += x * x;
                sy2 += y * y;
                x
            })
            .sum();
        let nf = self.len() as f64;
        (sxy - sy * sx / nf) / ((sx2 - sx * sx / nf) * (sy2 - sy * sy / nf)).sqrt()
    }
    /// Kendall Tau-B correlation.
    /// Defined by: tau = (conc - disc) / sqrt((conc + disc + tiesx) * (conc + disc + tiesy))
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vecg;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.kendalcorr(&v2),-0.07692307692307693);
    /// ```
    fn kendalcorr<U: Clone + Into<f64>>(self, v: &[U]) -> f64 {
        let (mut conc, mut disc, mut tiesx, mut tiesy) = (0_i64, 0_i64, 0_i64, 0_i64);
        for i in 1..self.len() {
            let x = self[i].clone().into();
            let y = v[i].clone().into();
            for j in 0..i {
                let xd = x - self[j].clone().into();
                let yd = y - v[j].clone().into();
                if !xd.is_normal() {
                    if !yd.is_normal() {
                        continue;
                    } else {
                        tiesx += 1;
                        continue;
                    }
                };
                if !yd.is_normal() {
                    tiesy += 1;
                    continue;
                };
                if (xd * yd).signum() > 0_f64 {
                    conc += 1
                } else {
                    disc += 1
                }
            }
        }
        (conc - disc) as f64 / (((conc + disc + tiesx) * (conc + disc + tiesy)) as f64).sqrt()
    }
    /// Spearman rho correlation.
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vecg;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.spearmancorr(&v2),-0.1076923076923077);
    /// ```
    fn spearmancorr<U: PartialOrd + Clone + Into<f64>>(self, v: &[U]) -> f64 {
        let xvec = self.rank(true);
        let yvec = v.rank(true); // rank from crate idxvec::merge
                                 // It is just Pearson's correlation of usize ranks
        xvec.ucorrelation(&yvec) // using Indices trait from idxvec
    }

    /// Delta gm that adding point self will cause
    fn contribvec_newpt(self, gm: &[f64], recips: f64) -> Result<Vec<f64>,RE> {
        let dv = self.vsub::<f64>(gm);
        let mag = dv.vmag();
        if !mag.is_normal() {
            return Err(re_error("arith","point being added is coincident with gm"));
        };
        // adding new unit vector (to approximate zero vector) and rescaling
        let recip = 1f64 / mag; 
        Ok(dv.vunit()?.smult::<f64>(recip / (recips + recip))) 
    }

    /// Normalized magnitude of change to gm that adding point self will cause
    fn contrib_newpt(self, gm: &[f64], recips: f64, nf: f64) -> Result<f64,RE> {
        let mag = self.vdist::<f64>(gm);
        if !mag.is_normal() {
            return Err(re_error("arith","point being added is coincident with gm"));
        };
        let recip = 1f64 / mag; // first had to test for division by zero
        Ok((nf + 1.0) / (recips + recip))
    }

    /// Delta gm caused by removing an existing set point self
    fn contribvec_oldpt(self, gm: &[f64], recips: f64) -> Result<Vec<f64>,RE> {
        let dv = self.vsub::<f64>(gm);
        let mag = dv.vmag();
        if !mag.is_normal() {
            return Err(re_error("arith","point being removed is coincident with gm")); 
        };
        let recip = 1f64 / mag; // first had to test for division by zero
        Ok(dv.vunit()?.smult::<f64>(recip / (recip - recips))) // scaling
    }

    /// Normalized Contribution that removing an existing set point p will make
    /// Is a negative number
    fn contrib_oldpt(self, gm: &[f64], recips: f64, nf: f64) -> Result<f64,RE> {
        let mag = self.vdist::<f64>(gm);
        if !mag.is_normal() {
            return Err(re_error("arith","point being removed is coincident with gm")); 
        };
        let recip = 1f64 / mag; // first had to test for division by zero
        Ok((nf - 1.0) / (recip - recips))
        // self.contribvec_oldpt(gm,recips,p).vmag()
    }

    /// Householder reflect
    fn house_reflect<U: Clone + PartialOrd + Into<f64>>(self, x: &[U]) -> Vec<f64> {
        x.vsub(&self.smult(x.dotp(self)))
    }
}
