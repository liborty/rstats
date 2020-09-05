use anyhow::{Result,Context,ensure};
use crate::{MutVectors,Vectors,emsg};

impl MutVectors for &mut[f64] {

   /// Scalar multiplication of a vector, mutates self
   fn mutsmult(&mut self, s:f64) {
     self.iter_mut().for_each(|x|{ *x*=s });
   }
   /// Vector subtraction, mutates self
   fn mutvsub(&mut self, v: &[f64]) {
     self.iter_mut().enumerate().for_each(|(i,x)|*x-=v[i])
   }
  /// Vector addition, mutates self
   fn mutvadd(&mut self, v: &[f64]) {
     self.iter_mut().enumerate().for_each(|(i,x)|*x+=v[i])
   }
  /// Mutate to unit vector
   fn mutvunit(&mut self) { 
      self.mutsmult(1_f64/self.iter().map(|x|x.powi(2)).sum::<f64>().sqrt())
   }
   /// Vector magnitude duplicated for mutable type 
   fn mutvmag(&mut self) -> f64 { self.iter().map(|x|x.powi(2)).sum::<f64>().sqrt() }

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

   /// Unit vector - creates a new one
   fn vunit(&self) -> Vec<f64> { self.smult(1_f64/self.vmag()) }

   /// Centroid = multidimensional arithmetic mean
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let mut pts = genvec(15,15,255,30);
   /// let centre = pts.as_slice().arcentroid(15);
   /// let dist = pts.as_slice().distsum(15,&centre);
   /// assert_eq!(dist, 4.14556218326653_f64);
   /// ```
   fn arcentroid(&self, d:usize) -> Vec<f64> {
      let n = self.len()/d;
      let mut centre = vec![0_f64;d];
      for i in 0..n {
         centre.as_mut_slice().mutvadd(self.get(i*d .. (i+1)*d).unwrap())
      }
      centre.as_mut_slice().mutsmult(1.0/n as f64);
      centre
   }

    /// Medoid is a point in n-dimensional set of points with the least sum of distances to all others.  
   /// This method returns an index to the start of medoid within a flat vector of d-dimensional points.  
   /// `d` is the number of dimensions = length of the slices. 
   /// Set of points (slices) is held as one flat `buff:&[f64]`.  
   /// This is faster than vec of vecs but users have to handle the indices.  
   /// Note: `medoid` computes each distance twice but it is probably faster than memoizing and looking them up,  
   /// unless the dimensionality is somewhat large. 
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let mut pts = genvec(15,15,255,30);
   /// let (dist,indx) = pts.as_slice().medoid(15).unwrap();
   /// assert_eq!(dist,4.812334638782327_f64);
   /// ```
   fn medoid(&self, d:usize) -> Result<(f64,usize)> {
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
   Ok((mindist,minindx))
   }

   /// The sum of distances of all points contained in &self to given point v.    
   /// v which minimises this objective function is the Geometric Median. 
   fn distsum(&self, d:usize, v:&[f64]) -> f64 {
      let n = self.len()/v.len();
      let mut sum = 0_f64;
      for i in 0..n {
         let thisp = self.get(i*d .. (i+1)*d).unwrap();
         sum += v.vdist(&thisp)              
      }
      sum
   }

   /// Weiszfeld's formula for one iteration step in finding gmedian.  
   /// It has known problems with choosing the starting point and may fail to converge.  
   /// However, nmedian below solves those problems.
   fn betterpoint(&self, d:usize, eps:f64, v:&[f64]) -> Result<(bool,Vec<f64>)> {
      let n = self.len()/d;
      let mut rsum = 0_f64;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = self.get(i*d .. (i+1)*d).unwrap();
         let dist = v.vdist(&thatp);
         ensure!(dist.is_normal(),emsg(file!(),line!(),"betterpoint encountered zero distance")); 
         let recip = 1.0/dist;
         rsum += recip;
         vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
      }
      vsum.as_mut_slice().mutsmult(1.0/rsum);
      if vsum.as_slice().vsub(&v).as_slice().vmag() < eps  
         { Ok((true,vsum)) } else { Ok((false,vsum)) }    
   }

   fn nextpoint(&self, d:usize, eps:f64, v:&[f64]) -> Result<(bool,Vec<f64>)> {
      let n = self.len()/d;
      let mut rsum = 0_f64;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = self.get(i*d .. (i+1)*d).unwrap();       
         let dist = v.vdist(&thatp);
         if dist < eps { // jump to nearby existing point
            vsum = self.firstpoint(d,i,&thatp)  // and search from there
               .with_context(||emsg(file!(),line!(),"nextpoint firstpoint call failed"))?; 
            if vsum.as_slice().vsub(&thatp).as_slice().vmag() < eps { 
               return Ok((true,vsum)) // moved less then eps from it, termination reached
               } else { return Ok((false,vsum))} // no termination, continue through the point  
         }
         let recip = 1.0/dist;
         rsum += recip;
         vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
      }
      vsum.as_mut_slice().mutsmult(1.0/rsum);
      if vsum.as_slice().vsub(&v).as_slice().vmag() < eps  
         { Ok((true,vsum)) } else { Ok((false,vsum)) }
   }

   /// My innovative step that guarantees convergence.
   fn firstpoint(&self, d:usize, indx:usize, v:&[f64]) -> Result<Vec<f64>> {
      // println!("Firstpoint used");
      let n = self.len()/d;
      let mut rsum = 0_f64;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         if i == indx { continue }; // exclude the starting medoid point  
         let thatp = self.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"firstpoint failed to extract point"))?;
         let recip = 1.0/v.vdist(&thatp);
         rsum += recip;
         vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
      }
      vsum.as_mut_slice().mutsmult(1.0/rsum);
      Ok(vsum)
   }

   /// Geometric Median is the point that minimises the sum of distances to a given set of points.
   /// This is the original Weiszfeld's algorithm for comparison.  
   /// It has problems with convergence and/or division by zero when the estimate
   /// runs too close to one of the existing points in the set.
   /// See test/tests.rs where test `gmedian` panics, whereas `nmedian` finds the correct result
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let mut pts = genvec(15,15,255,30);
   /// let (ds,gm) = pts.as_slice().gmedian(15, 1e-5).unwrap();
   /// assert_eq!(ds,4.126465898732421_f64);
   /// ```
   fn gmedian(&self, d:usize, eps:f64) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"gmedian d must divide vector length"));
      let mut oldpoint = self.arcentroid(d); // start with the centroid
      // let mut iterations = 0_usize;
      loop {
      //   iterations += 1;
         let (terminate,newpoint) = self.betterpoint(d,eps,&oldpoint)
            .with_context(||emsg(file!(),line!(),"gmedian nextpoint call failed"))?; // find new point 
         oldpoint = newpoint;
         if terminate { break }                
      }
      // println!("iterations: {}",iterations);
      Ok((self.distsum(d,&oldpoint),oldpoint))
   }
 
   /// Geometric Median is the point that minimises the sum of distances to a given set of points.  
   /// This improved  algorithm has guaranteed convergence. It will dodge any points in the set 
   /// which cause problems to Weiszfeld. It has similar running time on easy datasets 
   /// but guaranteed convergence also for the difficult cases. Eps controls the desired accuracy.
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let mut pts = genvec(15,15,255,30);
   /// let (ds,gm) = pts.as_slice().nmedian(15, 1e-5).unwrap();
   /// assert_eq!(ds,4.126465898732421_f64);
   /// ```
   fn nmedian(&self, d:usize, eps:f64) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"gmedian d must divide vector length"));
      let mut oldpoint = self.arcentroid(d); // start with the centroid
      // let mut iterations = 0_usize;
      loop {
      //   iterations += 1;
         let (terminate,newpoint) = self.nextpoint(d,eps,&oldpoint)
            .with_context(||emsg(file!(),line!(),"gmedian nextpoint call failed"))?; // find new point 
         oldpoint = newpoint;
         if terminate { break }                
      }
      // println!("iterations: {}",iterations);
      Ok((self.distsum(d,&oldpoint),oldpoint))
   }   
}