use crate::{Vecu8,here};

impl Vecu8 for &[u8] {

    /// Scalar product of two (positive) u8 slices.   
    /// Must be of the same length - no error checking (for speed)
    fn dotpu8(self, v: &[u8]) -> u64 {
        self.iter().zip(v).map(|(&xi, &vi)| (xi as u64)*(vi as u64)).sum::<u64>()
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
 
    /// Vector subtraction (converts results to f64 as they can be negative)
    fn vsubu8(self, v: &[u8]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| (xi as f64) - (vi as f64)).collect()
    }

    /// Vector addition ( converts results to f64, as they can exceed 255 )
    fn vaddu8(self, v: &[u8]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| (xi as f64)+(vi as f64)).collect()
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
    fn cityblockdu8(self, v:&[u8]) -> f64 {
        self.iter()
        .zip(v)
        .map(|(&xi, &vi)| { let d = xi as f64 -vi as f64; if d<0_f64 {-d} else {d} } ) 
        .sum::<f64>()      
    }
    ///Euclidian distance squared, the arguments are both of &[u8] type  
    fn vdistsqu8(self, v: &[u8]) -> u64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| {
               let x = xi as i32;
               let y = vi as i32;
                (x - y).pow(2) as u64})           
            .sum::<u64>()
    }
 
    /// Probability density function of bytes data
    fn pdfu8(self) -> Vec<f64> {  
        let nf = self.len() as f64;
        let mut pdfv = vec![0_f64;256];        
        for &x in self { pdfv[x as usize] += 1_f64/nf };       
        pdfv
    }

    /// Information (entropy) of &[u8] (in nats)
    fn entropyu8(self) -> f64 {
        let pdfv = self.pdfu8();
        let mut entr = 0_f64;
        for &x in self { 
            let p = pdfv[x as usize];
            if p.is_normal() // ignore zero probabilities
                { entr -= p*(p.ln()) }
        };            
        entr           
    }
    /// Joint probability density function (here just co-occurence counts) 
    /// of paired values in two vectors of bytes of the same length.
    /// Needs n^2 x 32bits of memory. Do not use for very long vectors, 
    /// those need hashing implementation.
    fn jointpdfu8(self, v:&[u8]) -> Vec<Vec<u32>> {  
        let n = self.len();
        if v.len() != n { panic!("{} argument vectors must be of equal length!",here!()) }
        let mut jocc = vec![vec![0_u32; 256]; 256];
        for i in 0..n {
            for j in 0..v.len() { jocc[self[i] as usize][v[j] as usize] += 1 }
        }
        jocc
    }
    /// Joint entropy of &[u8],&[u8] (in nats)
    fn jointentropyu8(self, v:&[u8]) -> f64 {
        let n = self.len();
        let m = v.len();        
        let jpdf = self.jointpdfu8(v);
        let mut entr = 0_f64;
        for i in 0..n {
            for j in 0..m {
                let cx = jpdf[self[i] as usize][v[j] as usize];
                if cx > 0 { // ignore zero counts
                    // turn counts into probabilities
                    let x = (cx as f64) / ((n*m) as f64); 
                    entr -= x*(x.ln()) 
                }
            }
        };            
        entr           
    }
    /// Dependence of two &[u8] variables, the range is [0,1],
    /// i.e. it returns 0 iff they are statistically independent
    /// and 1 when they are identical
    fn dependenceu8(self, v:&[u8]) -> f64 {     
        2.0 - (self.entropyu8() + v.entropyu8())/self.jointentropyu8(v)
    }
}
