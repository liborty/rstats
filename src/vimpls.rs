use anyhow::{Result,Context};
use crate::{Vectors,GMedian,NDPoints,emsg};

impl Vectors for Vec<f64> { 
   
   /// Scalar multiplication of a vector, creates new vec
   fn smult(&self, s:f64) -> Vec<f64> {
      self.iter().map(|&x|s*x).collect()
   }

   /// Scalar product of two f64 slices.   
   /// Must be of the same length - no error checking for speed
   fn dotp(&self, v: &[f64]) -> f64 {
      self.iter().enumerate().map(|(i,&x)| x*v[i]).sum::<f64>()    
   }

   /// Vector subtraction, creates a new Vec result
   fn vsub(&self, v: &[f64]) -> Vec<f64> {
      self.iter().enumerate().map(|(i,&x)|x-v[i]).collect()
   }

   /// Vector addition, creates a new Vec result
   fn vadd(&self, v: &[f64]) -> Vec<f64> {
      self.iter().enumerate().map(|(i,&x)|x+v[i]).collect()
   }

   /// Euclidian distance between two n dimensional points (vectors)
   fn vdist(&self, v: &[f64]) -> f64 {
      self.iter().enumerate().map(|(i,&x)|(x-v[i]).powi(2)).sum::<f64>().sqrt()
   }

   /// Vector magnitude 
   fn vmag(&self) -> f64 { self.iter().map(|&x|x.powi(2)).sum::<f64>().sqrt() }

   /// Unit vector
   fn vunit(&self) -> Vec<f64> { self.smult(1_f64/self.vmag()) }
      
}

impl Vectors for &[f64] { 
   
   /// Scalar multiplication of a vector, creates new vec
   fn smult(&self, s:f64) -> Vec<f64> {
      self.iter().map(|&x|s*x).collect()
   }

   /// Scalar product of two f64 slices.   
   /// Must be of the same length - no error checking for speed
   fn dotp(&self, v: &[f64]) -> f64 {
      self.iter().enumerate().map(|(i,&x)| x*v[i]).sum::<f64>()    
   }

   /// Vector subtraction, creates a new Vec result
   fn vsub(&self, v: &[f64]) -> Vec<f64> {
      self.iter().enumerate().map(|(i,&x)|x-v[i]).collect()
   }

   /// Vector addition, creates a new Vec result
   fn vadd(&self, v: &[f64]) -> Vec<f64> {
      self.iter().enumerate().map(|(i,&x)|x+v[i]).collect()
   }

   /// Euclidian distance between two n dimensional points (vectors).  
   /// Slightly faster than vsub followed by vmag, as both are done in one loop
   fn vdist(&self, v: &[f64]) -> f64 {
      self.iter().enumerate().map(|(i,&x)|(x-v[i]).powi(2)).sum::<f64>().sqrt()
   }

   /// Vector magnitude 
   fn vmag(&self) -> f64 { self.iter().map(|&x|x.powi(2)).sum::<f64>().sqrt() }

   /// Unit vector
   fn vunit(&self) -> Vec<f64> { self.smult(1_f64/self.vmag()) }
   
}

impl GMedian for NDPoints<'_> {
   /// Medoid is a point in n-dimensional set of points with the least sum of distances to others.
   /// This method returns an index to the start of medoid within an NDPoints set of n-dimensional points.
   /// This computes each distance twice but it is probably faster than memoizing and looking them up,
   /// unless the dimensionality is somewhat large  
   fn medoid(&self) -> Result<(usize,f64)> {
      let n = self.buff.len()/self.dims;
      let mut minindx = 0;
      let mut mindist = f64::MAX;
      for i in 0..n {
         let thisp = self.buff.get(i*self.dims .. (i+1)*self.dims)
            .with_context(||emsg(file!(),line!(),"medoid failed to get this slice"))?;
         let mut dsum = 0_f64;
         for j in 0..n {
            if i==j { continue };
            let thatp = self.buff.get(j*self.dims .. (j+1)*self.dims)
               .with_context(||emsg(file!(),line!(),"medoid failed to get that slice"))?;
            dsum += thisp.vdist(&thatp);
            if dsum >= mindist { break } // quit adding points if minimum distance is already exceeded
            }
         // println!("Distance: {}\n",dsum);
         if dsum < mindist { mindist = dsum; minindx = i };       
      }
   Ok((minindx,mindist))
   }

   /// Distances is the sum of distances of an arbitrary point to all points in NDPoints
   fn distances(&self,v: &[f64]) -> f64 {
      let n = self.buff.len()/self.dims;
      let mut sumdist = f64::MAX;
      for i in 0..n {
         let thisp = self.buff.get(i*self.dims .. (i+1)*self.dims).unwrap();
         sumdist += v.vdist(thisp)              
      }
      sumdist
   }

}