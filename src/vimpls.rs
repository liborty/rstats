use std::cmp::Ordering::Equal;
use crate::{Vectors,GMedian,NDPoints};

impl Vectors for Vec<f64> { 
   
   /// Sorts a mutable Vec<f64> in place.  
   /// It is the responsibility of the user to ensure that there are no NaNs etc.
   fn sortf(&mut self) { 
      self.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Equal))
   }

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

   fn medoid(&self) -> &[f64] {
      self.buff
   }
}