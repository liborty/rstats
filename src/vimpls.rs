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
   fn vunit(&self) ->Vec<f64> { 
      self.smult(1./self.iter().map(|x|x.powi(2)).sum::<f64>().sqrt())
   }

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

   /// `Eccentricity` of an existing d-dimensional point within the set, specified by its indx.
   /// It is a measure  between 0.0 and 1.0 of `not being a median`. It does not need the median. 
   /// The perfect median would have eccentricity zero.
   /// The medoid has the least ecentricity of the existing set points.
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

   /// Ecentricity measure and the eccentricity vector of any point.
   /// It is a measure  between 0.0 and 1.0 of `not being a median` but does not need the median.  
   /// The vector points towards the median. 
   fn veccentr(&self, d:usize, thisp:&[f64]) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d;
      let mut nf = n as f64;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = self.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"veccentr failed to extract that point"))?;
         let mut vdif = thatp.vsub(thisp);
         let mag = vdif.as_slice().vmag();
         if !mag.is_normal() { nf -= 1.0; continue }; // thisp belongs to the set
         // make vdif into a unit vector with its already known magnitude
         vdif.as_mut_slice().mutsmult(1./mag); 
         vsum.as_mut_slice().mutvadd(&vdif);   // add it to their sum
      }
      Ok((vsum.as_slice().vmag()/nf, vsum))
   }

   /// This convenience wrapper calls `veccentr` and extracts just the eccentricity (residual error for median).
   fn ecc(&self, d:usize, v:&[f64]) -> f64 {
      let (eccentricity,_) = self.veccentr(d,&v).unwrap();
      eccentricity
   }

   /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
   /// It provably only has iterative solutions over vectors. 
   /// Weiszfeld's formula had known problems with choosing the starting point and sometimes failing to converge.
   /// Especially in situations, where the points are dense in the close proximity to the gm.
   /// However, these problems are fixed here in my improved algorithm.      
   /// This is eventually going to be a multithreaded version.
   /// Results of `betterpoint` over arbitrary (data parallel) subsets will be simply added up to generate a new
   /// vector approximation.
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let pt = genvec(15,15,255,30);
   /// let pts = pt.as_slice();
   /// let gm = pts.nmedian(15, 1e-5).unwrap();
   /// let error = pts.ecc(15,&gm);
   /// assert_eq!(error,0.000004826966175302838_f64);
   /// ```
   fn nmedian(&self, d:usize, eps:f64) -> Result<Vec<f64>> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"gmedian d must divide vector length"));
      let mut oldpoint = self.arcentroid(d); // start with the centroid
      loop {
        let (rsum,mut newv) = betterpoint(self,d,&oldpoint)
            .with_context(||emsg(file!(),line!(),"nmedian betterpoint call failed"))?; // find new point 
         newv.as_mut_slice().mutsmult(1.0/rsum);
         if newv.as_slice().vdist(&oldpoint) < eps { oldpoint = newv; break };
         oldpoint = newv                       
      }
      Ok(oldpoint)
   }
}

/// betterpoint is called by nmedian.  
fn betterpoint(set:&[f64], d:usize, v:&[f64]) -> Result<(f64,Vec<f64>)> {
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
      vsum.as_mut_slice().mutvadd(&thatp.smult(recip));
   }
   Ok((rsum,vsum))    
}