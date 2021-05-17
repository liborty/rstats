use crate::{MutVectors, Vecf64, Indices};

impl Vecf64 for &[f64] {

    /// Scalar multiplication of a vector, creates new vec
    fn smult(self, s:f64) -> Vec<f64> {
        self.iter().map(|&x| s*x).collect()
    }
    
    /// Scalar addition to a vector, creates new vec
    fn sadd(self, s:f64) -> Vec<f64> {
        self.iter().map(|&x| s+x).collect()
    }    

    /// Scalar product of two f64 slices.   
    /// Must be of the same length - no error checking (for speed)
    fn dotp(self, v: &[f64]) -> f64 {
        self.iter().zip(v).map(|(&xi, &vi)| xi * vi).sum::<f64>()
    }

    fn vinverse(self) -> Vec<f64> {
        self.smult(1.0/self.vmagsq())
    }

    /// Cosine of an angle between two vectors.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.cosine(&v2),0.7517241379310344);
    /// ```
    fn cosine(self, v: &[f64]) -> f64 {
        let (mut sxy, mut sy2) = (0_f64, 0_f64);
        let sx2: f64 = self
            .iter()
            .zip(v)
            .map(|(&x, &y)| {
                sxy += x * y;
                sy2 += y * y;
                x*x
            })
            .sum();
        sxy / (sx2*sy2).sqrt()
    }


    /// Vector subtraction, creates a new Vec result
    fn vsub(self, v: &[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| xi - vi).collect()
    }

    /// Vector addition, creates a new Vec result
    fn vadd(self, v: &[f64]) -> Vec<f64> {
        self.iter().zip(v).map(|(&xi, &vi)| xi + vi).collect()
    }

    /// Euclidian distance between two n dimensional points (vectors).  
    /// Slightly faster than vsub followed by vmag, as both are done in one loop
    fn vdist(self, v: &[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (xi - vi).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    /// Euclidian distance squared between two n dimensional points (vectors).  
    /// Slightly faster than vsub followed by vmasq, as both are done in one loop
    /// Same as vdist without taking the square root
    fn vdistsq(self, v: &[f64]) -> f64 {
        self.iter()
            .zip(v)
            .map(|(&xi, &vi)| (xi - vi).powi(2))
            .sum::<f64>() 
    } 
    /// Vector magnitude
    fn vmag(self) -> f64 {
        self.iter().map(|&x| x.powi(2)).sum::<f64>().sqrt()
    }

    /// Vector magnitude squared
    fn vmagsq(self) -> f64 {
        self.iter().map(|&x| x.powi(2)).sum::<f64>()
    }

    /// Unit vector - creates a new one
    fn vunit(self) -> Vec<f64> {
        self.smult(1. / self.iter().map(|x| x.powi(2)).sum::<f64>().sqrt())
    }

    /// Area of a parallelogram between two vectors.
    /// Same as the magnitude of their cross product.
    /// Attains maximum `|a|.|b|` when the vectors are othogonal.
    fn varea(self, v:&[f64]) -> f64 {
        (self.vmagsq()*v.vmagsq() - self.dotp(v).powi(2)).sqrt()
    }

    /// Area proportional to the swept arc up to angle theta. 
    /// Attains maximum of `2|a||b|` when the vectors have opposite orientations.
    /// This is really |a||b|(1-cos(theta))
    fn varc(self, v:&[f64]) -> f64 { 
        (self.vmagsq()*v.vmagsq()).sqrt() - self.dotp(v)
    }
    
    /// Pearson's correlation coefficient of a sample of two f64 variables.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.correlation(&v2),-0.1076923076923077);
    /// ```
    fn correlation(self, v: &[f64]) -> f64 {
        let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let sx: f64 = self
            .iter()
            .zip(v)
            .map(|(&x, &y)| {
                sy += y;
                sxy += x * y;
                sx2 += x * x;
                sy2 += y * y;
                x
            })
            .sum();
        let nf = self.len() as f64;
        (sxy - sx / nf * sy) / ((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()
    }

    /// Kendall Tau-B correlation coefficient of a sample of two f64 variables.
    /// Defined by: tau = (conc - disc) / sqrt((conc + disc + tiesx) * (conc + disc + tiesy))
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.kendalcorr(&v2),-0.07692307692307693);
    /// ```
    fn kendalcorr(self, v: &[f64]) -> f64 {
        let (mut conc, mut disc, mut tiesx, mut tiesy) = (0_i64, 0_i64, 0_i64, 0_i64);
        for i in 1..self.len() {
            let x = self[i];
            let y = v[i];
            for j in 0..i {
                let xd = x - self[j];
                let yd = y - v[j];
                if !xd.is_normal() {
                    if !yd.is_normal() {
                        continue;
                    } else {
                        tiesx += 1;
                        continue;
                    }
                };
                if !yd.is_normal() {
                    tiesy += 1;
                    continue;
                };
                if (xd * yd).signum() > 0_f64 {
                    conc += 1
                } else {
                    disc += 1
                }
            }
        }
        (conc - disc) as f64 / (((conc + disc + tiesx) * (conc + disc + tiesy)) as f64).sqrt()
    }
    /// Spearman rho correlation coefficient of two f64 variables.
    /// This is the simplest implementation with no sorting.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// let v2 = vec![14_f64,1.,13.,2.,12.,3.,11.,4.,10.,5.,9.,6.,8.,7.];
    /// assert_eq!(v1.spearmancorr(&v2),-0.1076923076923077);
    /// ```
    fn spearmancorr(self, v: &[f64]) -> f64 {
        let xvec = self.mergerank();
        let yvec = v.mergerank();
        // It is just Pearson's correlation of ranks
        xvec.ucorrelation(&yvec)
    }

    /// Spearman correlation of five distances
    /// against Kazutsugi discrete outcomes [0.00,0.25,0.50,0.75,1.00], ranked as [4,3,2,1,0] 
    /// (the order is swapped to penalise distances). 
    /// The result is in the range [-1,1].
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1:Vec<f64> = vec![4.,1.,2.,0.,3.];
    /// assert_eq!(v1.kazutsugi(),0.3);
    /// ```
    fn kazutsugi(self) -> f64 {
        let xvec = self.mergerank();
        let yvec:Vec<f64> = vec![4.,3.,2.,1.,0.];
        let (mut sxy, mut sx2) = (0_f64,0_f64);
        const MY:f64 = 2.;  // y mean 
        const SY:f64 = 10.; // sum of yvec
        const SY2:f64 = 30.; // sum of y^2
        let sx:f64 = xvec.iter().zip(yvec).map(|(&ux,y)| {
            let x = ux as f64;
            sxy += x*y;
            sx2 += x*x;
            x }).sum();
        let covar = (sxy - sx*MY) / ((SY2 - SY*MY)*(sx2-sx/5.*sx)).sqrt();
        covar
    }

    /// (Auto)correlation coefficient of pairs of successive values of (time series) f64 variable.
    /// # Example
    /// ```
    /// use rstats::Vecf64;
    /// let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
    /// assert_eq!(v1.autocorr(),0.9984603532054123_f64);
    /// ```
    fn autocorr(self) -> f64 {
        let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let sx:f64 = self.windows(2).map(|w| {
            let x = w[0];
            let y = w[1];
            sy += y;
            sxy += x * y;
            sx2 += x * x;
            sy2 += y * y;
            x }).sum();
        let nf = self.len() as f64;
        (sxy - sx / nf * sy) / ((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()
    }

    /// Finds minimum, minimum's index, maximum, maximum's index of &[f64]
    /// Here self is usually some data, rather than a vector
    fn minmax(self) -> (f64, usize, f64, usize) {
        let mut min = self[0]; // initialise to the first value
        let mut mini = 0;
        let mut max = self[0]; // initialised as min, allowing 'else' below
        let mut maxi = 0;
        for i in 1..self.len() {
            let x = self[i];
            if x < min {
                min = x;
                mini = i
            } else if x > max {
                max = x;
                maxi = i
            }
        }
        (min, mini, max, maxi)
    }
    /// Linear transform to interval [0,1]
    fn lintrans(self) -> Vec<f64> {
        let (min,_,max,_) = self.minmax();
        let range = max-min;
        self.iter().map(|&x|(x-min)/range).collect()        
    }

    /// Counts how many items in sorted self are less than or equal to 
    /// the value v, using binary search. 
    fn binsearch(self, v: f64) -> usize {
        let mut lo = 0_usize; // index of the first item
        let mut hi = self.len()-1; // index of the last item
        if v < self[0] { return 0_usize }; // v is less than the first
        if v >= self[hi] { return hi+1 }; // v is at the top or over
        loop {
            let gap = hi - lo;
            if gap == 1 { return hi }
            let tryi = lo+gap/2;          
            // if value is above or equal, raise the low limit.
            // counts also repeating equal values. 
            if v >= self[tryi] { lo = tryi; continue };                 
            // else value is strictly below, reduce the high limit
            hi = tryi    
        }  
    }

    /// New sorted vector
    /// Copies self and then sorts it in place, leaving self unchanged (immutable).
    /// Calls mutsortf and that calls the standard self.sort_unstable_by
    fn sortf(self) -> Vec<f64> {
        let mut sorted:Vec<f64> = self.to_vec();
        sorted.mutsortf();
        sorted      
    }
    /// Returns new sorted vector, just as 'sortf' above
    /// but using our indexing 'mergesort' below
    fn sortm(self) -> Vec<f64> {
        self.mergesort(0,self.len()).unindex(&self)
    }

    /// Ranking of self by inverting the (merge) sort index.  
    /// Sort index is in the order of sorted items, giving their indices to the original data.
    /// Ranking is in the order of original data, giving their positions in the sort index.
    /// Very fast ranking of many f64 items, ranking `self` with only n*(log(n)+1) complexity.
    fn mergerank(self) -> Vec<usize> {
        let indx = self.mergesort(0,self.len());
        indx.revindex()    
    }    

    /// Recursive non-destructive merge sort. The data is read-only, it is not moved or mutated. 
    /// Returns vector of indices to self from i to i+n, such that the indexed values are in sort order.  
    /// Thus we are moving the index values instead of the actual values. 
    fn mergesort(self, i:usize, n:usize) -> Vec<usize> {

        if n == 1 { let res = vec![i]; return res };  // recursion termination
        if n == 2 {  // also terminate with two sorted items (for efficiency)          
            if self[i+1] < self[i] { return vec![i+1,i] } else { return vec![i,i+1] }
        }       
        let n1 = n / 2;  // the first half
        let n2 = n - n1; // the remaining second half
        let sv1 = self.mergesort(i, n1); // recursively sort the first half
        let sv2 = self.mergesort(i+n1, n2); // recursively sort the second half 

        // Now we merge the two sorted indexes into one
        let mut merged:Vec<usize> = Vec::with_capacity(n); 
        let mut k = 0_usize;   // subscript to the first sorted half
        let mut l = 0_usize;   // subscript to the second sorted half
        let mut firsthead;
        let mut secondhead;
        loop {
            // accessing the data in self only indirectly trough the sort indexes
            firsthead = self[sv1[k]];
            secondhead = self[sv2[l]];
             
            if firsthead < secondhead { // compare heads of the two sorted lists, pop the first
                merged.push(sv1[k]); k += 1 }
            else if firsthead > secondhead { // pop the second
                merged.push(sv2[l]); l += 1 } 
            else { // they are equal, so pop both, keeping their order
                merged.push(sv1[k]); k += 1;
                merged.push(sv2[l]); l += 1 }
            if k == sv1.len() {   // first one is empty, just copy the rest of the second 
                while l < sv2.len() {  merged.push(sv2[l]); l += 1  };
                break
            };
            if l == sv2.len() {   // second one is empty, just copy the rest of the first
                while k < sv1.len() {  merged.push(sv1[k]); k += 1 };
                break
            };            
        }
        return merged
    }

}