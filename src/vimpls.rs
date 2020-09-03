use anyhow::{Result,Context,ensure};
use crate::{Vectors,emsg};

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

   /// Medoid is a point in n-dimensional set of points with the least sum of distances to all others.
   /// This method returns an index to the start of medoid within a flat vector of d-dimensional points.  
   /// `d` is the number of dimensions = length of the slices. 
   /// Set of points (slices) is held as one flat `buff:&[f64]`.  
   /// This is faster than vec of vecs but users have to handle the indices.  
   /// Note: `medoid` computes each distance twice but it is probably faster than memoizing and looking them up,  
   /// unless the dimensionality is somewhat large 
   fn medoid(&self, d:usize) -> Result<(usize,f64)> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"medoid - d must divide vector length"));
      let mut minindx = 0;
      let mut mindist = f64::MAX;
      for i in 0..n {
         let thisp = self.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"medoid failed to get this slice"))?;
         let mut dsum = 0_f64;
         for j in 0..n {
            if i==j { continue };
            let thatp = self.get(j*d .. (j+1)*d)
               .with_context(||emsg(file!(),line!(),"medoid failed to get that slice"))?;
            dsum += thisp.vdist(&thatp);
            if dsum >= mindist { break } // quit adding points if minimum distance is already exceeded
            }
         // println!("Distance: {}\n",dsum);
         if dsum < mindist { mindist = dsum; minindx = i };       
      }
   Ok((minindx,mindist))
   }

   /// The sum of distances of all points contained in &self to given point v.    
   /// v which minimises this objective function is the Geometric Median. 
   fn distsum(&self, d:usize, v:&[f64]) -> f64 {
      let n = self.len()/v.len();
      let mut sum = 0_f64;
      for i in 0..n {
         let thisp = self.get(i*d .. (i+1)*d).unwrap();
         sum += v.vdist(thisp)              
      }
      sum
   }

   /// Weiszfeld's formula for one iteration step in finding gmedian.  
   /// Has problems with choosing the starting point - may fail to converge
   fn betterpoint(&self, d:usize, v:&[f64]) -> Vec<f64> {
      let n = self.len()/d;
      let mut sum = 0_f64;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = self.get(i*d .. (i+1)*d).unwrap();
         let recip = 1.0/v.vdist(&thatp);
         sum += recip;
         vsum = vsum.as_slice().vadd(&thatp.smult(recip));
      }
      vsum.as_slice().smult(1.0/sum)
   }

   /// My innovative first steps that guarantee good convergence
   fn firstpoint(&self, d:usize, indx:usize) -> Vec<f64> {
      let n = self.len()/d;
      let mut sum = 0_f64;
      let mut vsum = vec![0_f64;d];
      let v = self.get(indx*d .. (indx+1)*d).unwrap();
      for i in 0..n {
         if i == indx { continue }; // exclude the starting point in the set 
         let thatp = self.get(i*d .. (i+1)*d).unwrap();
         let recip = 1.0/v.vdist(&thatp);
         sum += recip;
         vsum = vsum.as_slice().vadd(&thatp.smult(recip));
      }
      vsum.as_slice().smult(1.0/sum)
   }
 
   /// Geometric Median is the point that minimises the sum of distances to a given set of points.
   /// This improved algorithm has guaranteed good convergence, as it will not approach any points
   /// in the set (which caused problems to Weiszfeld).
   /// eps controls the desired accuracy (low threshold on the improvements in the sum of distances)
   fn gmedian(&self, d:usize, eps:f64) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"gmedian d must divide vector length"));
      // start from medoid
      let (indx, mut dist) = self.medoid(d)
         .with_context(||emsg(file!(),line!(),"gmedian medoid call failed"))?;
      // println!("gmedian initial medoid distance: {}",dist);
      let mut newdist: f64;
      // innovative first iteration step from the medoid, excluding the medoid
      let mut point = self.firstpoint(d,indx);
      newdist = self.distsum(d,&point);
      let mut distdif = dist-newdist;
      ensure!(distdif>0.,emsg(file!(),line!(),"gmedian failed to improve"));
      // println!("gmedian first improvement: {}",distdif);  
      dist = newdist;
      while distdif > eps { // improve the termination condition
         point = self.betterpoint(d,&point); 
         newdist = self.distsum(d,&point);
         distdif = dist - newdist;
         ensure!(distdif>0.,emsg(file!(),line!(),"gmedian failed to improve"));
         // println!("gmedian improvement: {}",distdif); 
         dist = newdist;          
      }
      // println!("gmedian final distance: {}",dist); 
      Ok((dist,point))
   }


}