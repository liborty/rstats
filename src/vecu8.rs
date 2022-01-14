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
        for &x in self { pdfv[x as usize] += 1_f64 };
        pdfv.iter_mut().for_each(|p| if *p > 0.0 { *p /= nf });
        pdfv
    }

    /// Information (entropy) of &[u8] (in nats)
    fn entropyu8(self) -> f64 {
        let pdfv = self.pdfu8();
        pdfv.iter().map(|&p|if p > 0.0 {-p*(p.ln())} else {0.0}).sum::<f64>()    
    }

    /// Joint probability density function (here just co-occurence counts) 
    /// of successive pairs of values from two vectors of bytes 
    /// of the same lenghts n. Needs 4*256^2=262144 bytes of heap memory, 
    /// which will be sparse except for long input vectors. 
    fn jointpdfu8(self, v:&[u8]) -> Vec<Vec<u32>> {  
        let n = self.len();
        if v.len() != n { panic!("{} argument vectors must be of equal length!",here!()) }    
        let mut res:Vec<Vec<u32>> = vec![vec![0_u32;256];256]; 
        self.iter().zip(v).for_each(|(&si,&vi)| res[si as usize][vi as usize] += 1 ); 
        res
    }

    /// Joint entropy of &[u8],&[u8] (in nats)
    fn jointentropyu8(self, v:&[u8]) -> f64 {
        let n = self.len(); 
        let nf = n as f64;
        let mut entropy = 0_f64;
        // for short vecs, it is quicker to iterate through args
        if n < 65000 { 
            let mut jpdf = self.jointpdfu8(v); 
            self.iter().zip(v).for_each(|(&si,&vi)| {
              let c = jpdf[si as usize][vi as usize];
              if c > 0 {
                let p = (c as f64)/nf;                   
                entropy -= p*(p.ln()); 
                // prevent this pair's count being counted again
                jpdf[si as usize][vi as usize] = 0;  
              }
            });
            return entropy; // return value 
        } 
        // for long vecs, iterate through the counts array
        let jpdf = self.jointpdfu8(v); 
        for v in jpdf {
            for c in v {  
                if c > 0_u32 {
                    let p = (c as f64)/nf; 
                    entropy -= p*(p.ln()); 
                }
            }
        } 
        entropy              
    }

    /// Dependence in the range [0,1] of two &[u8] variables
    /// e.g. 0 is returned iff they are statistically pairwise independent
    fn dependenceu8(self, v:&[u8]) -> f64 {     
        1.0 - self.jointentropyu8(v)/(self.entropyu8() + v.entropyu8())
    }
}
