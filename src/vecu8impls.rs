use crate::{Vecu8,VecVecu8,MutVectors,Vecf64,functions};
use functions::emsg;

impl Vecu8 for &[u8] {

    /// Vector magnitude squared
    fn vmagsq(self) -> f64 {
        self.iter().map(|&x| (x as f64).powi(2)).sum::<f64>()
    } 
    /// Probability density function of bytes data
    fn pdf(self) -> Vec<f64> {  
        let nf = self.len() as f64;
        let mut pdfv = vec![0_f64;256];        
        for &x in self { pdfv[x as usize] += 1_f64 };
        for i in 0..255 {  pdfv[i] /= nf };
        pdfv
    }
    /// Information (entropy) of &[u8] (in nats)
    fn entropy(self) -> f64 {
        let pdfv = self.pdf();
        let mut entr = 0_f64;
        for x in pdfv {
            if x.is_normal() // ignore zero probabilities
                { entr -= x*(x.ln()) }
        };            
        entr           
    }
    /// Joint probability density function (actually just co-occurence counts) 
    //  of two vectors of bytes of the same length.
    /// Needs n^2 x 32bits of memory. Do not use for very long vectors, 
    /// those will need hashing implementation.
    fn jointpdf(self, v:&[u8]) -> Vec<Vec<u32>> {  
        let n = self.len();
        if v.len() != n {
            panic!("{}",emsg(
            file!(),line!(),"jointpdf argument vectors must be of equal length!"))
        }
        let mut jocc = vec![vec![0_u32; 256]; 256];
        for i in 0..n { jocc[self[i] as usize][v[i] as usize] += 1 };
        jocc
    }
    /// Joint entropy of &[u8],&[u8] (in nats)
    fn jointentropy(self, v:&[u8]) -> f64 {
        let nf = self.len() as f64;
        let pdfv = self.jointpdf(v);
        let mut entr = 0_f64;
        for i in 0..255 {
            for j in 0..255 {
                let cx = pdfv[i][j];
                if cx > 0 { // ignore zero counts
                    let x = cx as f64 / nf;  // turn count into probability
                    entr -= x*(x.ln()) 
                }
            }
        };            
        entr           
    }
    /// Dependence of two &[u8] variables, the range is [0,1],
    /// i.e. it returns 0 iff they are statistically independent
    /// and 1 when they are identical
    fn dependence(self, v:&[u8]) -> f64 {
        2.0 * (1.0 - self.jointentropy(v) / (self.entropy() + v.entropy()))
    }

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
    /// Area of swept arc between self &[u8] and v:&[f64]
    fn varc(self, v:&[f64]) -> f64 { 
        (self.vmagsq()*v.vmagsq()).sqrt() - self.dotp(v)
    }
    /// Euclidian distance between self &[u8] and v:&[f64].  
    /// Faster than vsub followed by vmag, as both are done in one loop
    fn vdist(self, v: &[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (xi as f64 - vi).powi(2))
            .sum::<f64>()
            .sqrt()
    }

}

impl VecVecu8 for &[Vec<u8>] {

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