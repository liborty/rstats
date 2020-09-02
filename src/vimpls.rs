use crate::{Vectors,GMedian,NDPoints};

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

   /// Euclidian distance between two n dimensional points (vectors)
   fn vdist(&self, v: &[f64]) -> f64 {
      self.iter().enumerate().map(|(i,&x)|(x-v[i]).powi(2)).sum::<f64>().sqrt()
   }

   /// Vector magnitude 
   fn vmag(&self) -> f64 { self.iter().map(|&x|x.powi(2)).sum::<f64>().sqrt() }
   
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

   /// Euclidian distance between two n dimensional points (vectors)
   fn vdist(&self, v: &[f64]) -> f64 {
      self.iter().enumerate().map(|(i,&x)|(x-v[i]).powi(2)).sum::<f64>().sqrt()
   }

   /// Vector magnitude 
   fn vmag(&self) -> f64 { self.iter().map(|&x|x.powi(2)).sum::<f64>().sqrt() }
   
}

impl GMedian for NDPoints<'_> {
   /// Medoid is a point in n-dimensional set of points with the least sum of distances to others.
   /// This method returns an index to the start of medoid within an NDPoints set of n-dimensional points.
   fn medoid(&self) -> usize {
      let n = self.buff.len()/self.dims;
      let mut minindx = 0;
      let mut mind = f64::MAX;
      for i in 0..n {
         let thisp = self.buff.get(i*self.dims .. (i+1)*self.dims).unwrap();
         let mut dsum = 0_f64;
         // This computes each distance twice but it is probably faster than looking them up somewhere,
         // unless the dimensionality is very large  
         for j in 0..n {
            if i==j { continue };
            let thatp = self.buff.get(j*self.dims .. (j+1)*self.dims).unwrap();
            dsum += thisp.vdist(&thatp);
            if dsum > mind { break } // quit testing points if min distance already exceeded
            }
         if dsum < mind { mind = dsum; minindx = i }
      }
      minindx
   }
}