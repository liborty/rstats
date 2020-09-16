use anyhow::{Result,Context,ensure};
use crate::{RStats,MutVectors,Vectors,Med};
use crate::functions::{scalarecc,betterpoint,emsg};

impl MutVectors for &mut[f64] {

   /// Scalar multiplication of a vector, mutates self
   fn mutsmult(self, s:f64) {
     self.iter_mut().for_each(|x|{ *x*=s });
   }
   /// Vector subtraction, mutates self
   fn mutvsub(self, v: &[f64]) {
     self.iter_mut().enumerate().for_each(|(i,x)|*x-=v[i])
   }
   /// Vector addition, mutates self
   fn mutvadd(self, v: &[f64]) {
     self.iter_mut().enumerate().for_each(|(i,x)|*x+=v[i])
   }
   /// Mutate to unit vector
   fn mutvunit(self) { 
      self.mutsmult(1_f64/self.iter().map(|x|x.powi(2)).sum::<f64>().sqrt())
   }
   /// Mutate a set of vetors of dimensions d to zero geometric median form.
   /// In more than one dimensions, this result is invariant with respect to rotation,
   /// unlike the often misguidedly used mean (`acentroid` here), which depends
   /// on the choice of axis. 
   /// To use separate 1-d medians for each axis is not right either.  
   /// For safety, such quasi-median is not even implemented by rstats.
   fn mutzeromd(self, d:usize, eps:f64) {
      let n = self.len()/d;
      let median = self.nmedian(d,eps,).unwrap();
      for i in 0..n {
         let point = self.get_mut(i*d .. (i+1)*d).unwrap();
         point.mutvsub(&median);
      }
   }
}

impl Vectors for &[f64] { 

   /// Retrieves from flat slice 'self' the (address of) sub-slice at index i,
   /// where the length of each sub-slice is d.
   /// Utility method for multidimensional flat vectors
   fn point(&self,d:usize,i:usize) -> &[f64] { self.get(i*d .. (i+1)*d).unwrap() }
   
   /// Scalar multiplication of a vector, creates new vec
   fn smult(self, s:f64) -> Vec<f64> {
      self.iter().map(|&x|s * x).collect()
   }
  
   /// Scalar product of two f64 slices.   
   /// Must be of the same length - no error checking for speed
   fn dotp(self, v: &[f64]) -> f64 {
      self.iter().enumerate().map(|(i,&x)| x*v[i]).sum::<f64>()    
   }

   /// Vector subtraction, creates a new Vec result
   fn vsub(self, v: &[f64]) -> Vec<f64> {
      self.iter().enumerate().map(|(i,&x)|x-v[i]).collect()
   }
   
   /// Vector addition, creates a new Vec result
   fn vadd(self, v: &[f64]) -> Vec<f64> {
      self.iter().enumerate().map(|(i,&x)|x+v[i]).collect()
   }
 
   /// Euclidian distance between two n dimensional points (vectors).  
   /// Slightly faster than vsub followed by vmag, as both are done in one loop
   fn vdist(self, v: &[f64]) -> f64 {
      self.iter().enumerate().map(|(i,&x)|(x-v[i]).powi(2)).sum::<f64>().sqrt()
   }

   /// Vector magnitude 
   fn vmag(self) -> f64 { self.iter().map(|&x|x.powi(2)).sum::<f64>().sqrt() }

   /// Unit vector - creates a new one
   fn vunit(self) ->Vec<f64> { 
      self.smult(1./self.iter().map(|x|x.powi(2)).sum::<f64>().sqrt())
   }
    
   /// Correlation coefficient of a sample of two f64 variables.
   /// # Example
   /// ```
   /// use rstats::Vectors;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// let v2 = vec![14_f64,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.];
   /// assert_eq!(v1.correlation(&v2).unwrap(),-1_f64);
   /// ```
   fn correlation(self,v:&[f64]) -> Result<f64> {
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"correlation - first sample is empty"));
      ensure!(n==v.len(),emsg(file!(),line!(),"correlation - samples are not of the same size"));
      let (mut sy,mut sxy,mut sx2,mut sy2) = (0_f64,0_f64,0_f64,0_f64);
      let sx:f64 = self.iter().enumerate().map(|(i,&x)| {
         let y = v[i]; 
         sy += y; sxy += x*y; sx2 += x*x; sy2 += y*y; x    
      }).sum();
   let nf = n as f64;
   Ok( (sxy-sx/nf*sy)/(((sx2-sx/nf*sx)*(sy2-sy/nf*sy)).sqrt()) )
   }

   /// Kendall Tau-B correlation coefficient of a sample of two f64 variables.
   /// Defined by: tau = (conc - disc) / sqrt((conc + disc + tiesx) * (conc + disc + tiesy))
   /// This is the simplest implementation with no sorting.
   /// # Example
   /// ```
   /// use rstats::Vectors;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// let v2 = vec![14_f64,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.];
   /// assert_eq!(v1.kendalcorr(&v2).unwrap(),-1_f64);
   /// ```
   fn kendalcorr(self,v:&[f64]) -> Result<f64> {
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"kendalcorr - first sample is empty"));
      ensure!(n==v.len(),emsg(file!(),line!(),"kendalcorr - samples are not of the same size"));
      let (mut conc, mut disc, mut tiesx, mut tiesy) = (0_i64,0_i64,0_i64,0_i64);
      for i in 1..n {
         let x = self[i];
         let y = v[i];
         for j in 0..i {
            let xd = x - self[j];
            let yd = y - v[j];
            if !xd.is_normal() {
               if !yd.is_normal() { continue } else { tiesx += 1; continue }
            };
            if !yd.is_normal() { tiesy += 1; continue };
            if (xd*yd).signum() > 0_f64 { conc += 1 } else { disc += 1 }               
            }
         }
      Ok((conc-disc) as f64/(((conc+disc+tiesx)*(conc+disc+tiesy)) as f64).sqrt())
   }
   /// Spearman rho correlation coefficient of a sample of two f64 variables.
   /// This is the simplest implementation with no sorting.
   /// # Example
   /// ```
   /// use rstats::Vectors;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// let v2 = vec![14_f64,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.];
   /// assert_eq!(v1.spearmancorr(&v2).unwrap(),-1_f64);
   /// ```
   fn spearmancorr(self,v:&[f64]) -> Result<f64> {
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"spearmancorr - first sample is empty"));
      ensure!(n==v.len(),emsg(file!(),line!(),"spearmancorr - samples are not of the same size"));
      let xvec = self.ranks().unwrap();
      let yvec = v.ranks().unwrap(); 
      let mx = xvec.ameanstd().unwrap();
      let my = yvec.ameanstd().unwrap();
      let mut covar = 0_f64;
      for i in 0..n {
         covar += (xvec[i]-mx.mean)*(yvec[i]-my.mean);
      }
      covar /= mx.std*my.std*(n as f64);
      // remove small truncation errors
      if covar > 1.0 { covar=1_f64 } else if covar < -1_f64 { covar=-1.0 }; 
      Ok(covar)
   }

   /// (Auto)correlation coefficient of pairs of successive values of (time series) f64 variable.
   /// # Example
   /// ```
   /// use rstats::Vectors;
   /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   /// assert_eq!(v1.autocorr().unwrap(),0.9984603532054123_f64);
   /// ```
   fn autocorr(self) -> Result<f64> {
      let n = self.len();
       ensure!(n>=2,emsg(file!(),line!(),"autocorr - sample is too small"));
       let (mut sx,mut sy,mut sxy,mut sx2,mut sy2) = (0_f64,0_f64,0_f64,0_f64,0_f64);
       for i in 0..n-1 {
          let x = self[i]; let y = self[i+1]; 
          sx += x; sy += y; sxy += x*y; sx2 += x*x; sy2 += y*y 
       }    
       let nf = n as f64;
       Ok( (sxy-sx/nf*sy)/(((sx2-sx/nf*sx)*(sy2-sy/nf*sy)).sqrt()) )
   }

   /// Centroid = simple multidimensional arithmetic mean
   /// # Example
   /// ```
   /// use rstats::{Vectors,functions::genvec};
   /// let pts = genvec(15,15,255,30);
   /// let centre = pts.acentroid(15);
   /// let dist = pts.distsum(15,&centre);
   /// assert_eq!(dist, 4.14556218326653_f64);
   /// ```
   fn acentroid(self, d:usize) -> Vec<f64> {
      let n = self.len()/d;
      let mut centre = vec![0_f64;d];
      for i in 0..n {
         centre.mutvadd(self.get(i*d .. (i+1)*d).unwrap())
      }
      centre.mutsmult(1.0/n as f64);
      centre
   }

    /// Finds minimum, minimum's index, maximum, maximum's index of &[f64]
    /// Here self is usually some data, rather than a vector
   fn minmax(self) -> (f64,usize,f64,usize) {
      let mut min = self[0]; // initialise to the first value
      let mut mini = 0;
      let mut max = self[0]; // initialised as min, allowing 'else' below
      let mut maxi = 0;
      for i in 1..self.len() {
         let x = self[i];
         if x < min { min = x; mini = i }
         else if x > max { max = x; maxi = i } 
      }
   (min,mini,max,maxi)
   }

   /// For each point, gives its sum of distances to all other points.
   /// This is the efficient workhorse of distances based analysis.
   fn distances(self, d:usize) -> Result<Vec <f64>> {  
     let n = self.len()/d;
     ensure!(n*d == self.len(),emsg(file!(),line!(),"distances - d must divide vector length"));
     let mut dists = vec![0_f64;n]; // distances accumulator for all points
      // examine all unique pairings (lower triangular part of symmetric flat matrix)
      for i in 1..n {
         let thisp = self.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"distances failed to get this slice"))?;
         for j in 0..i {
            let thatp = self.get(j*d .. (j+1)*d)
               .with_context(||emsg(file!(),line!(),"distances failed to get that slice"))?;
            let d = thisp.vdist(&thatp); // calculate each distance relation just once
            dists[i] += d; dists[j] += d;   // but add it to both points   
         }
      }
      Ok(dists)           
   }  
   
   /// The sum of distances from a set point given by its `indx` to all the other points in self.
   /// This method is suitable for a single point. For all the points, use more
   /// efficient `distances`.    
   fn distsuminset(self, d:usize, indx:usize) -> f64 {
      let n = self.len()/d;
      let mut sum = 0_f64;
      let thisp = self.get(  indx*d .. (indx+1)*d).unwrap();
      for i in 0..n {
         if i == indx { continue };
         let thatp = self.get(i*d .. (i+1)*d).unwrap();
         sum += thatp.vdist(&thisp)              
      }
      sum
   }

   /// The sum of distances from any point v (typically not in self) to all the points in self.    
   /// Geometric Median is defined as the point v which minimises this function. 
   fn distsum(self, d:usize, v:&[f64]) -> f64 {
      let n = self.len()/v.len();
      let mut sum = 0_f64;
      for i in 0..n {
         let thisp = self.get(i*d .. (i+1)*d).unwrap();
         sum += v.vdist(&thisp)              
      }
      sum
   }

   /// Medoid is the point belonging to set of points `self`,
   /// which has the least sum of distances to all other points. 
   /// Outlier is the point with the greatest sum of distances. 
   /// This function returns a four-tuple:  
   /// (medoid_distance, medoid_index, outlier_distance, outlier_index).
   /// `d` is the number of dimensions = length of the point sub-slices. 
   /// The entire set of points is held in one flat `&[f64]`.  
   /// This is faster than vec of vecs but we have to handle the indices.  
   /// # Example
   /// ```
   /// use rstats::{Vectors,functions::genvec};
   /// let pts = genvec(15,15,255,30);
   /// let (dm,_,_,_) = pts.medoid(15);
   /// assert_eq!(dm,4.812334638782327_f64);
   /// ```
   fn medoid(self, d:usize) -> (f64,usize,f64,usize) {
      self.distances(d).unwrap().minmax()
   }

   /// Eccentricity vector for each point.
   /// This is the efficient workhorse of eccentrities analysis. 
   fn eccentricities(self, d:usize) -> Result<Vec<Vec<f64>>> {  
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"distances - d must divide vector length"));
      // allocate vectors for the results
      let mut eccs = vec![vec![0_f64;d];n];
      // ecentricities vectors accumulator for all points
      // examine all unique pairings (lower triangular part of symmetric flat matrix)
      for i in 1..n {
         let thisp = self.get(i*d .. (i+1)*d).unwrap();
         for j in 0..i {
            let thatp = self.get(j*d .. (j+1)*d).unwrap();
            let e = thatp.vsub(&thisp).vunit(); // calculate each vector just once
            eccs[i].mutvadd(&e); 
            eccs[j].mutvsub(&e);  // mind the vector's orientation!   
         }
      }
   Ok(eccs)           
   }   

   /// Scalar positive measure of `not being a median` for a point belonging to the set.
   /// The point is specified by its index `indx`.
   /// The median does not have to be known. The perfect median would return zero.
   /// This is suitable for a single point. When eccentricities of all the points 
   /// are needed, use more efficient `eccentricities`. 
   fn eccentrinset(self, d:usize, indx:usize) -> f64 {
      let n = self.len()/d;
      let mut vsum = vec![0_f64;d];
      let thisp = self.get(indx*d .. (indx+1)*d).unwrap();
      for i in 0..n {
         if i == indx { continue }; // exclude this point  
         let thatp = self.get(i*d .. (i+1)*d).unwrap();       
         let unitdv = thatp.vsub(thisp).vunit();
         vsum.mutvadd(&unitdv);   // add it to their sum
      }
      vsum.vmag()/n as f64
   }

   /// Returns (Measure, Eccentricity-Vector) of any point (typically one not belonging to the set).
   /// The first (scalar) part of the result is a positive measure of `not being a median`.
   /// The second part is the eccentricity vector, which always points towards the median.
   /// The vector is of particular value and interest.
   /// This function has no prior knowledge of the actual median.  
   /// This is suitable for a single point. When eccentricities of all the points 
   /// are needed, use more efficient `eccentricities`. 
   fn veccentr(self, d:usize, thisp:&[f64]) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = self.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"veccentr failed to extract that point"))?;
         let mut vdif = thatp.vsub(thisp);
         let mag = vdif.vmag();
         if !mag.is_normal() { continue }; // thisp belongs to the set
         // use already known magnitude to compute vmag
         vdif.mutsmult(1./mag); 
         vsum.mutvadd(&vdif);   // add it to their sum
      }
      Ok((vsum.vmag()/n as f64, vsum))
   }

   /// This convenience wrapper calls `veccentr` and extracts just the eccentricity (residual error for median).
   /// Thus this method is the equivalent of `eccentr` 
   /// but suited for any explicitly given point, typically not belonging to the set.  
   /// When the eccentricity vector is needed, use `veccentr`
  fn ecc(self, d:usize, v:&[f64]) -> f64 {
      let (eccentricity,_) = self.veccentr(d,v).unwrap();
      eccentricity
   }

   /// Median of eccentricities measures (MOE).
   /// This is a new robust measure of spread of multidimensional points 
   /// (or multivariate sample).  
   fn moe(self, d:usize) -> Med {
      scalarecc(self.eccentricities(d).unwrap()) .median().unwrap()
   }

   /// Eccentricity defined Medoid and Outlier.
   /// This can give different results to `medoid` above, defined by sums of distances,
   /// especially for the outliers. See tests.rs.  
   /// Consider some point c and some other points, bunched up at a distance r from c.
   /// The sum of their distances will be n*r. Now, spread those points around a circle of radius r from c.
   /// The sum of their distances from c will remain the same but the eccentricity of c will be much reduced. 
   /// # Example
   /// ```
   /// use rstats::{Vectors,functions::genvec};
   /// let d = 6_usize;
   /// let pt = genvec(d,24,7,13); // random test data 5x20
   /// let (_medoideccentricity,medei,_outlierecccentricity,outei) = pt.emedoid(d);
   /// assert_eq!(medei,10); // index of e-medoid
   /// assert_eq!(outei,9);  // index of e-outlier
   /// ```
   fn emedoid(self, d:usize) -> (f64,usize,f64,usize) {
      scalarecc(self.eccentricities(d).unwrap()) .minmax()
   }

   /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
   /// It has (provably) only vector iterative solutions. 
   /// Search methods are slow and difficult in highly dimensional space. 
   /// Weiszfeld's fixed point iteration formula had known problems with sometimes failing to converge.
   /// Especially, when the points are dense in the close proximity of the gm, 
   /// or it coincides with one of them.  
   /// However, these problems are fixed in my new algorithm here.      
   /// There will eventually be a multithreaded version of `nmedian`.
   /// # Example
   /// ```
   /// use rstats::{Vectors,functions::genvec};
   /// let pt = genvec(15,15,255,30);
   /// let gm = pt.nmedian(15, 1e-5).unwrap();
   /// let error = pt.ecc(15,&gm);
   /// assert_eq!(error,0.000004826966175302838_f64);
   /// ```
   fn nmedian(self, d:usize, eps:f64) -> Result<Vec<f64>> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"gmedian d must divide vector length"));
      let mut oldpoint = self.acentroid(d); // start iterating from the centroid
      loop {
        let (rsum,mut newv) = betterpoint(self,d,&oldpoint)
            .with_context(||emsg(file!(),line!(),"nmedian betterpoint call failed"))?; // find new point 
         newv.mutsmult(1.0/rsum); // scaling the returned sum of unit vectors 
         if newv.vdist(&oldpoint) < eps { // test the magnitude of this move for termination
            oldpoint = newv; break // use the last small iteration anyway, as it is already computed
         };
         oldpoint = newv  // set up next iteration                     
      }
      Ok(oldpoint)
   }  

   /// Trend computes the vector connecting the geometric medians of two sets of multidimensional points.
   /// This is a robust relationship between two unordered multidimensional sets.
   /// The two sets have to be in the same space but can have different numbers of points.
   fn trend(self, d:usize, eps:f64, v:&[f64]) -> Vec<f64> {
      let m1 = self.nmedian(d,eps).unwrap();
      let m2 = v.nmedian(d,eps).unwrap();
      m2.vsub(&m1)
   }

   /// Generate a new set of vectors of dimensions d in zero (geometric) median form.
   /// Or subtract from them all any other vector `m`
   /// The geometric median is invariant with respect to rotation,
   /// unlike the often misguidedly used mean (`acentroid` here), or the quasi median,
   /// both of which depend on the choice of axis. 
   /// For safety, the quasi-median is not even implemented by rstats.
   /// Returns the zero medianized vectors as one flat vector.
   fn setsub(self, d:usize, m:&[f64]) -> Vec<f64> {
      let n = self.len()/d;
      let mut result = Vec::new();
      for i in 0..n {
         let point = self.get(i*d .. (i+1)*d).unwrap();
         let mut zp = point.vsub(m);
         result.append(&mut zp);
      }
      result
   }
}