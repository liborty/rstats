use crate::{Vecu8,RE,RError};

impl Vecu8 for &[u8] {

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
    fn jointpdfu8(self, v:&[u8]) -> Result<Vec<Vec<u32>>,RE> {  
        let n = self.len();
        if v.len() != n { 
            return Err(RError::DataError("jointpdfu8 argument vectors must be of equal length!")); }    
        let mut res:Vec<Vec<u32>> = vec![vec![0_u32;256];256]; 
        self.iter().zip(v).for_each(|(&si,&vi)| res[si as usize][vi as usize] += 1 ); 
        Ok(res)
    }

    /// Joint entropy of &[u8],&[u8] (in nats)
    fn jointentropyu8(self, v:&[u8]) -> Result<f64,RE> {
        let n = self.len(); 
        let nf = n as f64;
        let mut entropy = 0_f64;
        // for short vecs, it is quicker to iterate through args
        if n < 65000 { 
            let mut jpdf = self.jointpdfu8(v)?; 
            for (&si,&vi) in self.iter().zip(v) {
              let c = jpdf[si as usize][vi as usize];
              if c > 0 {
                let p = (c as f64)/nf;                   
                entropy -= p*(p.ln()); 
                // prevent this pair's count being counted again
                jpdf[si as usize][vi as usize] = 0;  
              } 
            };
        return Ok(entropy); // return value 
        } 
        // for long vecs, iterate through the counts array
        let jpdf = self.jointpdfu8(v)?; 
        for v in jpdf {
            for c in v {  
                if c > 0_u32 {
                    let p = (c as f64)/nf; 
                    entropy -= p*(p.ln()); 
                }
            }
        } 
        Ok(entropy)              
    }

    /// Statistical pairwise dependence in range [0,1] of two &[u8] variables
    /// returns 0 iff they are statistically pairwise independent
    /// returns 1 if they are identical or all values are unique
    fn dependenceu8(self, v:&[u8]) -> Result<f64,RE> {     
        Ok((self.entropyu8() + v.entropyu8())/self.jointentropyu8(v)?-1.0)
    }

    /// Independence in the range [1,2] of two &[u8] variables
    /// e.g. 2 is returned iff they are statistically pairwise independent
    /// returns 1 if they are identical or all values are unique
    fn independenceu8(self, v:&[u8]) -> Result<f64,RE> {     
        Ok(2.0*self.jointentropyu8(v)?/(self.entropyu8() + v.entropyu8()))
    }
}
