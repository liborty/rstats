use crate::{Vecu8,Vecf64,functions::emsg};

impl Vecu8 for &[u8] {

    /// Scalar multiplication of a vector, creates new vec
    fn smult(self, s:f64) -> Vec<f64> {
        self.iter().map(|&x| s*(x as f64)).collect()
     }

     /// Scalar addition to a vector, creates new vec
     fn sadd(self, s:f64) -> Vec<f64> {
        self.iter().map(|&x| s+(x as f64)).collect()
     }

    /// Scalar product.   
    /// Must be of the same length - no error checking (for speed)
    fn dotp(self, v: &[f64]) -> f64 {
        self.iter().zip(v).map(|(&xi, &vi)| (xi as f64)*vi).sum::<f64>()
    }
    /// Scalar product of two (positive) u8 slices.   
    /// Must be of the same length - no error checking (for speed)
    fn dotpu8(self, v: &[u8]) -> u64 {
        self.iter().zip(v).map(|(&xi, &vi)| (xi as u64)*(vi as u64)).sum::<u64>()
    }

    /// Cosine between &[u8] and &[f64].
    fn cosine(self, v: &[f64]) -> f64 {
        let (mut sxy, mut sy2) = (0_f64, 0_f64);
        let sx2: f64 = self
            .iter()
            .zip(v)
            .map(|(&ux, &y)| {
                let x  = ux as f64;             
                sxy += x * y;
                sy2 += y * y;
                x*x
            })
            .sum();
        sxy / (sx2*sy2).sqrt()
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
             x*x
         })
         .sum();
     sxy / (sx2*sy2).sqrt()
    }

    /// Vector subtraction 
    fn vsub(self, v: &[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| (xi as f64) - vi).collect()
    }
    /// Vector subtraction (converts results to f64 as they can be negative)
    fn vsubu8(self, v: &[u8]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| (xi as f64) - (vi as f64)).collect()
    }

    /// Vector addition
    fn vadd(self, v: &[f64]) -> Vec<f64> {
    self.iter().zip(v).map(|(&xi, &vi)| (xi as f64) + vi).collect()
    }
    /// Vector addition ( converts results to f64, as they can exceed 255 )
    fn vaddu8(self, v: &[u8]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| (xi as f64)+(vi as f64)).collect()
    }

    /// Vector magnitude
    fn vmag(self) -> f64{
        self.iter().map(|&x| (x as f64).powi(2)).sum::<f64>().sqrt()
    }

    /// Vector magnitude squared
    fn vmagsq(self) -> f64 {
        self.iter().map(|&x| (x as f64).powi(2)).sum::<f64>()
    }

    /// Euclidian distance between self &[u8] and v:&[f64].  
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
            .map(|(&xi, &vi)| ((xi as f64)-(vi as f64)).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    /// cityblock distance
    fn cityblockd(self, v:&[f64]) -> f64 {
        self.iter()
        .zip(v)
        .map(|(&xi, &vi)| (xi as f64 -vi).abs()) 
        .sum::<f64>()      
    }
    /// cityblock distance
    fn cityblockdu8(self, v:&[u8]) -> f64 {
        self.iter()
        .zip(v)
        .map(|(&xi, &vi)| { let d = xi as f64 -vi as f64; if d<0_f64 {-d} else {d} } ) 
        .sum::<f64>()      
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
    /// Area of swept arc between self &[u8] and v:&[f64]
    /// = |a||b|(1-cos(theta)) = 2|a||b|D
    fn varc(self, v:&[f64]) -> f64 { 
        (self.vmagsq()*v.vmagsq()).sqrt() - self.dotp(v)
    } 
    /// We define vector similarity S in the interval [0,1] as
    /// S = (1+cos(theta))/2
    fn vsim(self, v:&[f64]) -> f64 { (1.0+self.cosine(v))/2.0 }

    /// We define vector dissimilarity D in the interval [0,1] as
    /// D = 1-S = (1-cos(theta))/2
    fn vdisim(self, v:&[f64]) -> f64 { (1.0-self.cosine(v))/2.0 }

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
