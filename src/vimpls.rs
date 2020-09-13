use anyhow::{Result,Context,ensure};
use crate::{RStats,MutVectors,Vectors,Med,minmax,emsg};

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
   /// Vector magnitude duplicated for mutable type 
   fn mutvmag(self) -> f64 { self.iter().map(|x|x.powi(2)).sum::<f64>().sqrt() }
}

impl Vectors for &[f64] { 
   
   /// Scalar multiplication of a vector, creates new vec
   fn smult(self, s:f64) -> Vec<f64> {
      self.iter().map(|&x|s*x).collect()
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
   /// assert_eq!(v1.as_slice().correlation(&v2).unwrap(),-1_f64);
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
   /// assert_eq!(v1.as_slice().kendalcorr(&v2).unwrap(),-1_f64);
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
   /// assert_eq!(v1.as_slice().spearmancorr(&v2).unwrap(),-1_f64);
   /// ```
   fn spearmancorr(self,v:&[f64]) -> Result<f64> {
      let n = self.len();
      ensure!(n>0,emsg(file!(),line!(),"spearmancorr - first sample is empty"));
      ensure!(n==v.len(),emsg(file!(),line!(),"spearmancorr - samples are not of the same size"));
      let xvec = self.ranks().unwrap();
      let yvec = v.ranks().unwrap(); 
      let mx = xvec.as_slice().ameanstd().unwrap();
      let my = yvec.as_slice().ameanstd().unwrap();
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
   /// assert_eq!(v1.as_slice().autocorr().unwrap(),0.9984603532054123_f64);
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

   /// Centroid = multidimensional arithmetic mean
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let pts = genvec(15,15,255,30);
   /// let centre = pts.as_slice().acentroid(15);
   /// let dist = pts.as_slice().distsum(15,&centre);
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

   /// For each point, gives its sum of distances to all other points
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

   /// Medoid is the point belonging to set of points `self`,
   /// which has the least sum of distances to all other points. 
   /// Outlier is the point with the greatest sum of distances. 
   /// This function returns a four-tuple:  
   /// (medoid_distance, medoid_index, outlier_distance, outlier_index).
   /// `d` is the number of dimensions = length of the point sub-slices. 
   /// The entire set of points is held in one flat `&[f64]`.  
   /// This is faster than vec of vecs but we have to handle the indices.  
   /// Note: `medoid` computes each distance twice but it is probably faster than memoizing and looking them up,  
   /// unless the dimensionality is somewhat large; and it saves memory.
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let pts = genvec(15,15,255,30);
   /// let (dm,_,_,_) = pts.as_slice().medoid(15).unwrap();
   /// assert_eq!(dm,4.812334638782327_f64);
   /// ```
   fn medoid(self, d:usize) -> (f64,usize,f64,usize) {
      minmax(&self.distances(d).unwrap())
   }

   /// The sum of distances of all points in &self to given point v.    
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

   /// `Eccentricity` of a d-dimensional point belonging to the set self, specified by its indx.  
   /// Eccentricity is a measure between 0.0 and 1.0 of  a point `not being a median` of the given set. It does not need the median. 
   /// The perfect median has eccentricity zero.
   /// Of all the set points, Medoid has the lowest ecentricity and Outlier the highest.
   fn eccentr(self, d:usize, indx:usize) -> f64 {
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
      vsum.as_slice().vmag()/n as f64
   }

   /// Ecentricity measure and the eccentricity vector of any point (typically not one of the set).
   /// It is a measure  between 0.0 and 1.0 of `not being a median` but does not need the median.
   /// The eccentricity vector points towards the median and has maximum possible magnitude of n.
   fn veccentr(self, d:usize, thisp:&[f64]) -> Result<(f64,Vec<f64>)> {
      let n = self.len()/d;
      let mut vsum = vec![0_f64;d];
      for i in 0..n {
         let thatp = self.get(i*d .. (i+1)*d)
            .with_context(||emsg(file!(),line!(),"veccentr failed to extract that point"))?;
         let mut vdif = thatp.vsub(thisp);
         let mag = vdif.as_slice().vmag();
         if !mag.is_normal() { continue }; // thisp belongs to the set
         // make vdif into a unit vector with its already known magnitude
         vdif.as_mut_slice().mutsmult(1./mag); 
         vsum.as_mut_slice().mutvadd(&vdif);   // add it to their sum
      }
      Ok((vsum.as_slice().vmag()/n as f64, vsum))
   }

   /// This convenience wrapper calls `veccentr` and extracts just the eccentricity (residual error for median).
   fn ecc(self, d:usize, v:&[f64]) -> f64 {
      let (eccentricity,_) = self.veccentr(d,&v).unwrap();
      eccentricity
   }

   /// We now define MOE (median of ecentricities), a new measure of spread of multidimensional points 
   /// (or multivariate sample)  
   fn moe(self, d:usize) -> Med {
      let n = self.len()/d;
      let mut eccs = vec![0_f64;n];
      for i in 0..n { eccs[i] = self.eccentr(d, i) }
      eccs.as_slice().median().unwrap()
   }

   /// Geometric Median (gm) is the point that minimises the sum of distances to a given set of points.
   /// It has (provably) only vector iterative solutions. 
   /// Searching methods are slow and difficult in highly dimensional space. 
   /// Weiszfeld's fixed point iteration formula had known problems with sometimes failing to converge.
   /// Especially, when the points are dense in the close proximity of the gm, 
   /// or it coincides with one of them.  
   /// However, these problems are fixed in my improved algorithm here.      
   /// There is eventually going to be a multithreaded version of `nmedian`.
   /// # Example
   /// ```
   /// use rstats::{Vectors,genvec};
   /// let pt = genvec(15,15,255,30);
   /// let pts = pt.as_slice();
   /// let gm = pts.nmedian(15, 1e-5).unwrap();
   /// let error = pts.ecc(15,&gm);
   /// assert_eq!(error,0.000004826966175302838_f64);
   /// ```
   fn nmedian(self, d:usize, eps:f64) -> Result<Vec<f64>> {
      let n = self.len()/d;
      ensure!(n*d == self.len(),emsg(file!(),line!(),"gmedian d must divide vector length"));
      let mut oldpoint = self.acentroid(d); // start with the centroid
      loop {
        let (rsum,mut newv) = betterpoint(self,d,&oldpoint)
            .with_context(||emsg(file!(),line!(),"nmedian betterpoint call failed"))?; // find new point 
         newv.as_mut_slice().mutsmult(1.0/rsum); // adding unit vectors
         if newv.as_slice().vdist(&oldpoint) < eps { 
            oldpoint = newv; break // use the last iteration anyway
         };
         oldpoint = newv                       
      }
      Ok(oldpoint)
   }
}

/// betterpoint is called by nmedian. 
/// Scaling by rsum is left as the final step at higher level, 
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