use anyhow::{Result,Context,ensure};
use crate::{RStats,MutVectors,Vectors,Med};
use crate::f64impls::{emsg};
use std::cmp::Ordering::Equal;
use std::fmt;

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
   /// use rstats::{Vectors,vimpls::genvec};
   /// let pts = genvec(15,15,255,30);
   /// let centre = pts.acentroid(15);
   /// let dist = pts.distsum(15,&centre);
   /// assert_eq!(dist, 4.14556218326653_f64);
   /// ```
   fn acentroid(self, d:usize) -> Vec<f64> {
      let n = self.len()/d;
      let mut centre = vec![0_f64;d];
      for i in 0..n {
         centre.as_mut_slice().mutvadd(self.get(i*d .. (i+1)*d).unwrap())
      }
      centre.as_mut_slice().mutsmult(1.0/n as f64);
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
   
   /// The sum of distances from within-set point given by indx to all points in self.    
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

   /// The sum of distances from any point v to all points in self.    
   /// Geometric Median is defined as v which minimises this function. 
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
   /// use rstats::{Vectors,vimpls::genvec};
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
            eccs[i].as_mut_slice().mutvadd(&e); 
            eccs[j].as_mut_slice().mutvsub(&e);  // mind the vector's orientation!   
         }
      }
   Ok(eccs)           
   }   

   /// Scalar positive measure of `not being a median` for a point belonging to the set.
   /// The point is specified by its index indx.
   /// The median is not needed for this, however the perfect median would return zero.
   fn eccentr(self, d:usize, indx:usize) -> f64 {
      let n = self.len()/d;
      let mut vsum = vec![0_f64;d];
      let thisp = self.get(indx*d .. (indx+1)*d).unwrap();
      for i in 0..n {
         if i == indx { continue }; // exclude this point  
         let thatp = self.get(i*d .. (i+1)*d).unwrap();       
         let unitdv = thatp.vsub(thisp).vunit();
         vsum.as_mut_slice().mutvadd(&unitdv);   // add it to their sum
      }
      vsum.vmag()/n as f64
   }

   /// Returns (Measure, Eccentricity-Vector) of any point (typically one not belonging to the set).
   /// The first (scalar) part of the result is a positive measure of `not being a median`.
   /// The second part is the eccentricity vector, which points towards the median.
   /// The vector is of particular value and interest.
   /// This function has no prior knowledge of the actual median.
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
         vdif.as_mut_slice().mutsmult(1./mag); 
         vsum.as_mut_slice().mutvadd(&vdif);   // add it to their sum
      }
      Ok((vsum.vmag()/n as f64, vsum))
   }

   /// This convenience wrapper calls `veccentr` and extracts just the eccentricity (residual error for median).
   fn ecc(self, d:usize, v:&[f64]) -> f64 {
      let (eccentricity,_) = self.veccentr(d,v).unwrap();
      eccentricity
   }

   /// Median of eccentricities measures.
   /// This is a new robust measure of spread of multidimensional points 
   /// (or multivariate sample).  
   fn moe(self, d:usize) -> Med {
      scalarecc(self.eccentricities(d).unwrap()) .median().unwrap()
   }

   /// Medoid and Outlier as defined by the eccentricities
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
   /// use rstats::{Vectors,vimpls::genvec};
   /// let pt = genvec(15,15,255,30);
   /// let gm = pt.nmedian(15, 1e-5).unwrap();
   /// let error = pt.ecc(15,&gm);
   /// assert_eq!(error,0.000004826966175302838_f64);
   /// ```
   fn nmedian(self, d:usize, eps:f64) -> Result<Vec<f64>> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"gmedian d must divide vector length"));
      let mut oldpoint = self.acentroid(d); // start with the centroid
      loop {
        let (rsum,mut newv) = betterpoint(self,d,&oldpoint)
            .with_context(||emsg(file!(),line!(),"nmedian betterpoint call failed"))?; // find new point 
         newv.as_mut_slice().mutsmult(1.0/rsum); // adding vectors
         if newv.vdist(&oldpoint) < eps { 
            oldpoint = newv; break // use the last iteration anyway
         };
         oldpoint = newv                       
      }
      Ok(oldpoint)
   }  
}

/// betterpoint is called by nmedian. 
/// Scaling by rsum is left as the final step at calling level, 
/// in order to facilitate data parallelism. 
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

