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

   
  /// Geometric Median is the point that minimises the sum of distances to a given set of points.
   /// This is the original Weiszfeld's algorithm (for comparison). It is not recommended for production use.
   /// It has problems with convergence and/or division by zero when the iterative estimate
   /// runs too close to one of the existing points in the set.
   /// See test/tests.rs, where test `gmedian` panics, whereas `nmedian` finds the correct result.
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
         let (terminate,newpoint) = betterpoint(self,d,eps,&oldpoint)
            .with_context(||emsg(file!(),line!(),"gmedian nextpoint call failed"))?; // find new point 
         oldpoint = newpoint;
         if terminate { break }                
      }
      // println!("iterations: {}",iterations);
      Ok((self.distsum(d,&oldpoint),oldpoint))
   }
 
   /// Geometric Median is the point that minimises the sum of distances to a given set of points.  
   /// This is an improved  algorithm. It dodges any points in the set 
   /// which typically cause problems to Weiszfeld. It runs slightly faster on easy datasets and maintains its
   /// convergence even on the difficult data (dense points near the geometric median).  
   /// Eps controls the desired accuracy.
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let mut pts = genvec(15,15,255,30);
   /// let (ds,gm) = pts.as_slice().nmedian(15, 1e-5).unwrap();
   /// assert_eq!(ds,4.126465898732421_f64);
   /// ```
   fn nmedian(&self, d:usize, eps:f64) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"nmedian d must divide vector length"));
      let mut oldpoint = self.arcentroid(d); // start with the centroid
      // let slc = self.get(0 .. d)
      //    .with_context(||emsg(file!(),line!(),"nmedian failed to get starting point"))?;
      //  let mut oldpoint = slc.to_vec();
      // let mut iterations = 0_usize;
      loop {
      //   iterations += 1;
         let (terminate,newpoint) = nextpoint(self,d,eps,&oldpoint)
            .with_context(||emsg(file!(),line!(),"nmedian nextpoint call failed"))?; // find new point 
         oldpoint = newpoint;
         if terminate { break }                
      }
      // println!("iterations: {}",iterations);
      Ok((self.distsum(d,&oldpoint),oldpoint))
   }   
/*
   fn pmedian(&self, d:usize, eps:f64) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d/2;
      ensure!(2*n*d == self.len(),emsg(file!(),line!(),"pmedian d must divide vector length"));
      // let mut oldpoint = self.arcentroid(d); // start with the centroid
      let slc1 = self.get(0 .. n*d)
         .with_context(||emsg(file!(),line!(),"nmedian failed to get starting point"))?;
      let (_d1,p1) = slc1.nmedian(d,eps).unwrap();
      println!("Distance of point 1: {}", self.distsum(d,&p1));
      let slc2 = self.get(n*d .. self.len())
         .with_context(||emsg(file!(),line!(),"nmedian failed to get starting point"))?;
      let (_d2,p2) = slc2.nmedian(d,eps).unwrap();
      println!("Distance of point 2: {}", self.distsum(d,&p2));
      let oldpoint = p1.as_slice().vadd(&p2).as_slice().smult(0.5_f64);
      Ok((self.distsum(d,&oldpoint),oldpoint.to_vec()))
   }   
*/
}

   /// nextpoint is called by nmedian; it checks for proximity of set points and
   /// deals with the situation by calling firstpoint.
   fn nextpoint(set:&[f64], d:usize, eps:f64, v:&[f64]) -> Result<(bool,Vec<f64>)> {
      let n = set.len()/d;
      let reps = eps / 10.0;
      let mut rsum = 0_f64;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = set.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"nextpoint failed to extract other point"))?;       
         let dist = v.vdist(&thatp);
         if dist < reps { // jump to nearby existing point
            vsum = firstpoint(set,d,i,&thatp)  // and search from there
               .with_context(||emsg(file!(),line!(),"nextpoint firstpoint call failed"))?; 
            if vsum.as_slice().vsub(&thatp).as_slice().vmag() < reps { 
               return Ok((true,vsum)) // moved less then eps/10 from it, termination reached
            } else { return Ok((false,vsum))} // no termination, continue iterating through the point  
         }
         let recip = 1.0/dist;
         rsum += recip;
         vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
      }
      vsum.as_mut_slice().mutsmult(1.0/rsum);
      if vsum.as_slice().vsub(&v).as_slice().vmag() < eps  
         { Ok((true,vsum)) } else { Ok((false,vsum)) }
   }

   /// innovative estimate from the set point that guarantees convergence.
   fn firstpoint(set:&[f64], d:usize, indx:usize, v:&[f64]) -> Result<Vec<f64>> {
      // println!("Firstpoint");
      let n = set.len()/d;
      let nf = (n as f64 - 1.)/n as f64;
      let mut rsum = 0_f64;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         if i == indx { continue }; // exclude the starting medoid point  
         let thatp = set.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"firstpoint failed to extract point"))?;
         let recip = 1.0/v.vdist(&thatp);
         rsum += recip;
         vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
      }
      vsum.as_mut_slice().mutsmult(nf/rsum);
      Ok(vsum)
   }

 
   /// Weiszfeld's formula for one iteration step in finding the geometric median (gm).
   /// It has known problems with choosing the starting point and may fail to converge.
   /// Especially in situations where the points are dense in the vicinity of the gm.
   /// However, nmedian below solves those problems.  
   /// betterpoint is called by gmedian.  
   fn betterpoint(set:&[f64], d:usize, eps:f64, v:&[f64]) -> Result<(bool,Vec<f64>)> {
      let n = set.len()/d;
      let mut rsum = 0_f64;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = set.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"betterpoint failed to extract other point"))?; 
         let dist = v.vdist(&thatp);
         ensure!(dist.is_normal(),emsg(file!(),line!(),"betterpoint collided with an existing point")); 
         let recip = 1.0/dist;
         rsum += recip;
         vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
      }
      vsum.as_mut_slice().mutsmult(1.0/rsum);
      if vsum.as_slice().vsub(&v).as_slice().vmag() < eps  
         { Ok((true,vsum)) } else { Ok((false,vsum)) }    
   }