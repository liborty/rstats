use crate::{Vecu8,Vecf64,MutVectors,VecVecu8};

impl Vecu8 for &[u8] {
     
    /// Scalar multiplication of a vector, creates new vec
    fn smult(self, s:f64) -> Vec<f64> {
       self.iter().map(|&x| s*x as f64).collect()
    }
    /// Scalar addition to a vector, creates new vec
    fn sadd(self, s:f64) -> Vec<f64> {
       self.iter().map(|&x| s+x as f64).collect()
    }
    /// Scalar product.   
    /// Must be of the same length - no error checking (for speed)
    fn dotp(self, v: &[f64]) -> f64 {
        self.iter().zip(v).map(|(&xi, &vi)| xi as f64 * vi).sum::<f64>()
    }
    /// Scalar product of two (positive) u8 slices.   
    /// Must be of the same length - no error checking (for speed)
    fn dotpu8(self, v: &[u8]) -> u64 {
        self.iter().zip(v).map(|(&xi, &vi)| (xi * vi)as u64).sum::<u64>()
    }
    /// Cosine between two (positive) u8 slices.
    fn cosineu8(self, v: &[u8]) -> f64 {
        let (mut sxy, mut sy2) = (0_f64, 0_f64);
        let sx2: f64 = self
            .iter()
            .zip(v)
            .map(|(&ux, &uy)| {
                let x  = ux as f64;
                let y = uy as f64;
                sxy += x * y;
                sy2 += y * y;
                x*x as f64
            })
            .sum();
        sxy / (sx2*sy2).sqrt()
    }

    /// Vector magnitude squared
    fn vmagsq(self) -> f64 {
        self.iter().map(|&x| (x as f64).powi(2)).sum::<f64>()
    } 
    /// Area of swept arc
    fn varc(self, v:&[f64]) -> f64 { 
        (self.vmagsq()*v.vmagsq()).sqrt() - self.dotp(v)
    }
    /// Euclidian distance between self u8 point and v: &[f64].  
    /// Slightly faster than vsub followed by vmag, as both are done in one loop
    fn vdist(self, v: &[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (xi as f64 - vi).powi(2))
            .sum::<f64>()
            .sqrt()
    }

}

impl VecVecu8 for  &[Vec<u8>] {

    fn acentroid(self) -> Vec<f64> {
    let mut centre = vec![0_f64; self[0].len()];
    for v in self {
        centre.mutvaddu8(&v)
    }
    centre.mutsmult(1.0 / self.len() as f64);
    centre
    }

    fn nmedian(self, eps: f64) -> Vec<f64> {
        let mut oldpoint = self.acentroid(); // start iterating from the centroid
        loop {
            let (rsum, mut newv) = self.betterpoint(&oldpoint);
            newv.mutsmult(1.0 / rsum); // scaling the returned sum of unit vectors
            // test the magnitude of the move for termination
            if newv.vdist(&oldpoint) < eps {
                oldpoint = newv; // use the last small iteration anyway
                break; // from the loop
            };
            oldpoint = newv // set up the next iteration
        }
        oldpoint
    }

    /// betterpoint is called by nmedian.
    /// Scaling by rsum is left as the final step at calling level,
    /// in order to facilitate data parallelism.
    fn betterpoint(self, v: &[f64]) -> (f64, Vec<f64>) {
        let mut rsum = 0_f64;
        let mut vsum = vec![0_f64; v.len()];
        for thatp in self {
            let dist = thatp.vdist(&v);
            if dist.is_normal() {
                let recip = 1.0 / dist;
                rsum += recip;
                vsum.mutvadd(&thatp.smult(recip))
            }
        }
        (rsum, vsum)
    }
}