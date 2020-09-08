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
   /// let pts = genvec(15,15,255,30);
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
   /// let pts = genvec(15,15,255,30);
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

   /// multidimensional `eccentricity` of a point within the set
   /// based on geometric median without actually having to find the median
   /// it is a measure of `not being a median`. The actual median will have eccentricity zero
   /// and the medoid will have the least ecentricity
   fn eccentr(&self, d:usize, indx:usize) -> f64 {
      let n = self.len()/d;
      let mut vsum = vec![0_f64;d];
      let thisp = self.get(indx*d .. (indx+1)*d).unwrap();
      //  .with_context(||emsg(file!(),line!(),"eccentricity failed to extract this point"))?;   
      for i in 0..n {
         if i == indx { continue }; // exclude this point  
         let thatp = self.get(i*d .. (i+1)*d).unwrap();
         //    .with_context(||emsg(file!(),line!(),"eccentricity failed to extract that point"))?;
         let unitdv = thatp.vsub(thisp).as_slice().vunit();
         vsum.as_mut_slice().mutvadd(&unitdv);   // add it to their sum
      }
      vsum.as_slice().vmag()/(n-1) as f64
   }
   /// ecentricity of a point not in the set
   fn exteccentr(&self, d:usize, thisp:&[f64]) -> f64 {
      let n = self.len()/d;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = self.get(i*d .. (i+1)*d).unwrap();
         //    .with_context(||emsg(file!(),line!(),"eccentricity failed to extract that point"))?;
         let unitdv = thatp.vsub(thisp).as_slice().vunit();
         vsum.as_mut_slice().mutvadd(&unitdv);   // add it to their sum
      }
      vsum.as_slice().vmag()/n as f64
   }

   
  /// Geometric Median is the point that minimises the sum of distances to a given set of points.
   /// This is the original Weiszfeld's algorithm (for comparison). It is not recommended for production use.
   /// It has problems with convergence and/or division by zero when the iterative estimate
   /// runs too close to one of the existing points in the set.
   /// See test/tests.rs, where test `gmedian` panics, whereas `nmedian` finds the correct result.
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let pts = genvec(15,15,255,30);
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
   /// which typically cause problems to Weiszfeld. It does so in a more sophisticated way than gmedian. 
   /// It maintains convergence even on the difficult data (dense points near the geometric median).  
   /// Use nmedian in preference to gmedian on difficult data.  
   /// Eps controls the desired accuracy.
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let pts = genvec(15,15,255,30);
   /// let (ds,gm) = pts.as_slice().nmedian(15, 1e-5).unwrap();
   /// assert_eq!(ds,4.126465898732421_f64);
   /// ```
   fn nmedian(&self, d:usize, eps:f64) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"nmedian d must divide vector length"));
      let mut oldpoint = self.arcentroid(d); // start with the centroid
      //let slc = self.get(0 .. d)
      //   .with_context(||emsg(file!(),line!(),"nmedian failed to get starting point"))?;
      //let mut oldpoint = firstpoint(self,d,0,slc)
      //   .with_context(||emsg(file!(),line!(),"nmedian firstpoint call failed"))?;
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
}

/// nextpoint is called by nmedian; it checks for proximity of set points and
/// deals with the near coincidence by calling firstpoint.
fn nextpoint(set:&[f64], d:usize, eps:f64, v:&[f64]) -> Result<(bool,Vec<f64>)> {
   let n = set.len()/d;
   let mut rsum = 0_f64;
   let mut vsum = vec![0_f64;d];
   for i in 0..n {
      let thatp = set.get(i*d .. (i+1)*d)
         .with_context(||emsg(file!(),line!(),"nextpoint failed to extract other point"))?;       
      let dist = v.vdist(&thatp);
      if dist < eps/2.0 { // jump to nearby existing set point
         vsum = firstpoint(set,d,i,thatp)  // and search from there
            .with_context(||emsg(file!(),line!(),"nextpoint firstpoint call failed"))?; 
         if vsum.as_slice().vmag() < eps { 
               return Ok((true,thatp.vadd(&vsum))) } // moved less then eps from v, termination reached
         else { return Ok((false,thatp.vadd(&vsum)))} // continue iterating through the point
      }
      let recip = 1.0/dist;
      rsum += recip;
      vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
   }
   vsum.as_mut_slice().mutsmult(1.0/rsum);
   if vsum.as_slice().vsub(&v).as_slice().vmag() < eps  
      { Ok((true,vsum)) } else { Ok((false,vsum)) }
}

/// Novel estimate from the coincident set point v. Guarantees convergence.  
/// Returns an improvement vector from point v.
fn firstpoint(set:&[f64], d:usize, indx:usize, v:&[f64]) -> Result<Vec<f64>> {
   // println!("Firstpoint");
   let n = set.len()/d;
   let mut rsum = 0_f64;
   let mut vsum = vec![0_f64;d];
   for i in 0..n {
      if i == indx { continue }; // exclude the starting medoid point  
      let thatp = set.get(i*d .. (i+1)*d)
          .with_context(||emsg(file!(),line!(),"firstpoint failed to extract point"))?;
      let mut difv = thatp.vsub(v);
      let dist = difv.as_slice().vmag();
      if !dist.is_normal() { continue }; // another coinciding point - ignore it
      let invmag = 1.0/dist;
      rsum += invmag;
      difv.as_mut_slice().mutsmult(invmag); // make difv a unit vector
      vsum.as_mut_slice().mutvadd(&difv);   // add it to their sum
   }
   vsum.as_mut_slice().mutsmult(1.0/rsum); // and scale
   Ok(vsum)
}

 
/// Weiszfeld's formula for one iteration step in finding the geometric median (gm).
/// It has known problems with choosing the starting point and may fail to converge.
/// Especially in situations, where the points are dense in the close proximity to the gm.
/// It is fixed here in the simplest way.    
/// betterpoint is called by gmedian.  
fn betterpoint(set:&[f64], d:usize, eps:f64, v:&[f64]) -> Result<(bool,Vec<f64>)> {
   let n = set.len()/d;
   let mut rsum = 0_f64;
   let mut vsum = vec![0_f64;d];
   for i in 0..n {
      let thatp = set.get(i*d .. (i+1)*d)
         .with_context(||emsg(file!(),line!(),"betterpoint failed to extract other point"))?; 
      let dist = v.vdist(&thatp);
      if !dist.is_normal() { continue };  
      let recip = 1.0/dist;
      rsum += recip;
      vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
   }
   vsum.as_mut_slice().mutsmult(1.0/rsum);
   if vsum.as_slice().vsub(&v).as_slice().vmag() < eps  
      { Ok((true,vsum)) } else { Ok((false,vsum)) }    
}