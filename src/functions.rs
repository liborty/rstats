use std::cmp::Ordering::Equal;
use std::fmt;
use anyhow::{Result,Context};
use crate::{Vectors,MutVectors};

/// Sum of linear weights 
pub fn wsum(n: usize) -> f64 { (n*(n+1)) as f64/2. }

/// helper function for formatting error messages
pub fn emsg(file:&'static str, line:u32, msg:&'static str)-> String {
   format!("{}:{} rstats {}",file,line,msg)
}

/// betterpoint is called by nmedian. 
/// Scaling by rsum is left as the final step at calling level, 
/// in order to facilitate data parallelism. 
pub fn betterpoint(set:&[f64], d:usize, v:&[f64]) -> Result<(f64,Vec<f64>)> {
   let n = set.len()/d;
   let mut rsum = 0_f64;
   let mut vsum = vec![0_f64;d];
   for i in 0..n {
      let thatp = set.get(i*d .. (i+1)*d)
         .with_context(||emsg(file!(),line!(),"betterpoint failed to extract that point"))?; 
      let dist = v.vdist(&thatp);
      if !dist.is_normal() { continue };  
      let recip = 1.0/dist;
      rsum += recip;
      vsum.mutvadd(&thatp.smult(recip));
   }
   Ok((rsum,vsum))    
}


/// Converts the set of vectors produced by `eccentricities`
/// to their magnitudes normalised by n.
/// the output can be typically passed to `median` 
/// or `minmax` to find the Outlier and the Medoid
pub fn scalarecc(ev:Vec<Vec<f64>>) -> Vec<f64> {
   let mut scalars = Vec::new();
   let n = ev.len();
   let nf = n as f64;
   for i in 0..n { scalars.push(ev[i].vmag()/nf) };
   scalars
}

/// Sorts a mutable `Vec<f64>` in place.  
/// It is the responsibility of the user to ensure that there are no NaNs etc.
pub fn sortf(v: &mut [f64]) { 
   v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
}

/// Generates a random f64 vector of size d x n suitable for testing. It needs two seeds.  
/// Uses local closure `rand` to generate random numbers (avoids dependencies).  
/// Random numbers are in the open interval 0..1 with uniform distribution.  
pub fn genvec(d:usize, n:usize, s1:u32, s2:u32 ) -> Vec<f64> {
   let size = d*n;
   // change the seeds as desired
   let mut m_z = s1 as u32;
   let mut m_w = s2 as u32;
   let mut rand = || {
      m_z = 36969 * (m_z & 65535) + (m_z >> 16);
      m_w = 18000 * (m_w & 65535) + (m_w >> 16);
      (((m_z << 16) & m_w) as f64 + 1.0)*2.328306435454494e-10
   };
   let mut v = Vec::with_capacity(size); 
   for _i in 0..size { v.push(rand()) }; // fills the lot with random numbers
   return v
}

/// GreenIt struct facilitates printing (in green) any type
/// that has Display implemented.
pub struct GreenIt<T: fmt::Display>(pub T);
impl<T: fmt::Display> fmt::Display for GreenIt<T> {
   fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
      write!(f,"\x1B[01;92m{}\x1B[0m",self.0.to_string())  
   }
}

/// GreenVec struct facilitates printing (in green) of vectors of any type
/// that has Display implemented.
pub struct GreenVec<T: fmt::Display>(pub Vec<T>);
impl<T: fmt::Display> fmt::Display for GreenVec<T> {
   fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
      let mut s = String::from("\x1B[01;92m[");
      let n = self.0.len();
      if n > 0 {
         s.push_str(&self.0[0].to_string()); // first item
         for i in 1..n {
            s.push_str(", ");
            s.push_str(&self.0[i].to_string());
         }
      }   
      write!(f,"{}]\x1B[0m", s)
   }
}

