use crate::{Vecu8,VecVecu8,MutVectors,Vecf64,functions};
use functions::emsg;

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
    fn cosine(self, v: &[u8]) -> f64 {
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
    /// Vector subtraction (converts results to f64 as they can be negative)
    fn vsubu8(self, v: &[u8]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| xi as f64 - vi as f64).collect()
    }
    /// Vector addition ( converts results to f64, as they can exceed 255 )
    fn vadd(self, v: &[u8]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| (xi + vi) as f64).collect()
    }
    /// Vector magnitude
    fn vmag(self) -> f64{
        self.iter().map(|&x| (x as f64).powi(2)).sum::<f64>().sqrt()
    }
    /// Vector magnitude squared
    fn vmagsq(self) -> f64 {
        self.iter().map(|&x| (x as f64).powi(2)).sum::<f64>()
    } 
    /// Euclidian distance between self &[u8] and v:&[u8].  
    fn vdist(self, v:&[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (((xi as f64) - vi).powi(2)))
            .sum::<f64>()
            .sqrt()
    }    
    /// Euclidian distance between self &[u8] and v:&[u8].  
    /// Faster than vsub followed by vmag, as both are done in one loop
    fn vdistu8(self, v: &[u8]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (xi as f64 - vi as f64).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    ///Euclidian distance squared, the arguments are both of &[u8] type  
    fn vdistsq(self, v: &[u8]) -> u64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| {
               let x = xi as i32;
               let y = vi as i32;
                (x - y).pow(2) as u64})           
            .sum::<u64>()
    }
    /// Area proportional to the swept arc/// Area of swept arc between self &[u8] and v:&[f64]
    fn varc(self, v:&[f64]) -> f64 { 
        (self.vmagsq()*v.vmagsq()).sqrt() - self.dotp(v)
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
    /// Joint probability density function (here just co-occurence counts) 
    /// of paired values in two vectors of bytes of the same length.
    /// Needs n^2 x 32bits of memory. Do not use for very long vectors, 
    /// those need hashing implementation.
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
    /// Potentially useful clone-recast of &[u8] to &[f64] 
    fn vecu8asvecf64 (self) -> Vec<f64> {
        self.iter().map(| &x | x as f64).collect()
    }
}

impl VecVecu8 for &[Vec<u8>] {

    fn acentroid(self) -> Vec<f64> {
    let mut centre = vec![0_f64; self[0].len()];
    for v in self {
        centre.mutvaddu8(&v)
    }
    centre.mutsmult(1.0 / (self.len() as f64));
    centre
    }

    /// Eccentricity vector for a non member point.
    /// The true geometric median is as yet unknown.
    /// Returns the eccentricity vector.
    /// The true geometric median would return zero vector.
   fn eccnonmember(self, p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for x in self { 
            let dvmag = x.vdist(&p);
            if !dvmag.is_normal() { continue } // zero distance, safe to ignore
            let rec = 1.0/dvmag;
            vsum.mutvadd(&x.smult(rec)); // add unit vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmult(1.0/recip);
        vsum.vsub(&p)
    }
    /// Secant method with recovery from divergence
    /// for finding the geometric median of vectors of vectors of bytes
    fn gmedian(self, eps: f64) -> Vec<f64> {
        let mut p1 = self.acentroid();     
        let e1 = self.eccnonmember(&p1); // eccentricity vector1 
        let mut e1mag = e1.vmag(); 
        let mut p2 = p1.vadd(&e1);   
        loop {
            let e2 = self.eccnonmember(&p2); // eccentricity vector2
            let e2mag = e2.vmag(); 
            if e2mag < eps  { return p2 }; 
            let newp;         
            if e1mag > e2mag {  // eccentricity magnitude decreased, good, employ secant
                newp = p2.vadd(&e2.smult(p1.vsub(&p2).vmag()/(e1mag-e2mag)))                   
            }
            else { // probably overshot the minimum, go nearer, back 
               newp = p2.vadd(&e2.smult(p1.vsub(&p2).vmag()/(e1mag+e2mag)))                    
            } 
            p1 = p2;        
            p2 = newp;  
            e1mag = e2mag             
        }       
    }    
}