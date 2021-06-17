use crate::{MutVectors, Vecf64, Indices, here};

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

    // negated vector (components with opposite sign)
    fn negv(self) -> Vec<f64> { self.smult(-1.0) }

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

    /// cityblock distance
    fn cityblockd(self, v:&[f64]) -> f64 {
        self.iter()
        .zip(v)
        .map(|(&xi, &vi)| (xi-vi).abs()) 
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
    /// Same as the magnitude of their cross product |a ^ b| = |a||b|sin(theta).
    /// Attains maximum `|a|.|b|` when the vectors are othogonal.
    fn varea(self, v:&[f64]) -> f64 {
        (self.vmagsq()*v.vmagsq() - self.dotp(v).powi(2)).sqrt()
    }

    /// Area proportional to the swept arc up to angle theta. 
    /// Attains maximum of `2|a||b|` when the vectors have opposite orientations.
    /// This is really |a||b|(1-cos(theta)) = 2|a||b|D
    fn varc(self, v:&[f64]) -> f64 { 
         (self.vmagsq()*v.vmagsq()).sqrt() - self.dotp(v)
    }

    /// We define vector similarity S in the interval [0,1] as
    /// S = (1+cos(theta))/2
    fn vsim(self, v:&[f64]) -> f64 { (1.0+self.cosine(&v))/2.0 }

    /// We define vector dissimilarity D in the interval [0,1] as
    /// D = 1-S = (1-cos(theta))/2
    fn vdisim(self, v:&[f64]) -> f64 { (1.0-self.cosine(&v))/2.0 }

    
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
        let xvec = self.rank(true);
        let yvec = v.rank(true);
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
        let xvec = self.rank(true);
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

    /// Returns index to the first item that is strictly greater than v, 
    /// using binary search of an ascending sorted list.
    /// When none are greater, returns self.len(). 
    /// User must check for this index overflow: if the returned index == 0, then v was below the list,
    /// else use index-1 as a valid index to the last item that is less than or equal to v.
    /// This then is the right index to use for looking up cummulative probability density functions. 
    fn binsearch(self, v: f64) -> usize {
        let n = self.len();
        if n < 2 { panic!("{} list is too short!",here!()) }
        if v < self[0] { return 0_usize }; // v is smaller than the first item
        let mut hi = n-1; // valid index of the last item
        if v > self[hi] { return n }; // indicates that v is greater than the last item
        let mut lo = 0_usize; // initial index of the low limit 
    
        loop {
            let gap = hi - lo;
            if gap <= 1 { return hi }
            let tryi = lo+gap/2; 
            // if tryi index's value is above v, reduce the high index
            if self[tryi] > v { hi = tryi; continue }            
            // else indexed value is not greater than v, raise the low index;
            // jumps also repeating equal values. 
            lo = tryi
        }  
    }

    /// Merges two ascending sorted vectors &[f64]
    fn merge(self, v: &[f64]) -> Vec<f64> {
        let mut resvec:Vec<f64> = Vec::new();
        let l1 = self.len();
        let l2 = v.len();
        let mut i1 = 0;
        let mut i2 = 0;
        loop {
            if i1 == l1 { // self is now processed
                for i in i2..l2 { resvec.push(v[i]) } // copy out the rest of v
                break // and terminate
            }
            if i2 == l2 { // v is now processed
                for i in i1..l1 { resvec.push(self[i])} // copy out the rest of self
                break // and terminate
            }
            if self[i1] < v[i2] { resvec.push(self[i1]); i1 += 1; continue };
            if self[i1] > v[i2] { resvec.push(v[i2]); i2 += 1; continue }; 
            // here they are equal, so consume both
            resvec.push(self[i1]); i1 += 1;
            resvec.push(v[i2]); i2 += 1
        }
        resvec
    }
 
    /// Merges two ascending sorted vectors' indices, returns concatenated Vec<f64> and new index into it.
    /// Mostly just a wrapper for merge_indices()
    fn merge_immutable(self, idx1: &[usize], v2: &[f64], idx2: &[usize]) -> ( Vec<f64>,Vec<usize> ) {
        let resvec = [self,v2].concat(); // no sorting, just concatenation 
        let l = idx1.len();
        let idx2shifted:Vec<usize> = idx2.iter().map(|x| l+x ).collect(); // shift up the second index
        let residx = resvec.merge_indices(idx1,&idx2shifted);   
        ( resvec, residx )
    }

    /// Merges indices of two already concatenated sorted vectors: 
    /// self is untouched, only sort indices are merged.
    /// Used by `mergesort` and `merge_immutable`. 
    fn merge_indices(self, idx1:&[usize], idx2:&[usize]) -> Vec<usize> {
        let l1 = idx1.len();
        let l2 = idx2.len();
        let mut residx:Vec<usize> = Vec::new(); 
        let mut i1 = 0;  let mut i2 = 0;
        let mut head1 = self[idx1[i1]]; let mut head2 = self[idx2[i2]];
        loop {
            if head1 < head2 { 
                residx.push(idx1[i1]);
                i1 += 1;  
                if i1 == l1 { // idx1 is now fully processed
                    for i in i2..l2 { residx.push(idx2[i]) } // copy out the rest of idx2
                    break // and terminate
                }
                head1 = self[idx1[i1]]; // else move to the next value
                continue
            }
            if head1 > head2 { 
                residx.push(idx2[i2]); 
                i2 += 1; 
                if i2 == l2 { // idx2 is now processed
                    for i in i1..l1 { residx.push(idx1[i]) } // copy out the rest of idx1
                    break // and terminate
                }                    
                head2 = self[idx2[i2]]; 
                continue
            } 
            // here the heads are equal, so consume both
            residx.push(idx1[i1]); 
            i1 += 1; 
            if i1 == l1 { // idx1 is now fully processed
                for i in i2..l2 { residx.push(idx2[i]) } // copy out the rest of idx2
                break // and terminate
            }
            head1 = self[idx1[i1]];
            residx.push(idx2[i2]); 
            i2 += 1; 
            if i2 == l2 { // idx2 is now processed
                for i in i1..l1 { residx.push(idx1[i]) } // copy out the rest of idx1
                break // and terminate
            }                    
            head2 = self[idx2[i2]];            
       }
        residx
    }
    
    /// Immutable sort. Returns new sorted vector, just like 'sortf' above
    /// but using our indexing 'mergesort' below.
    /// Simply passes the boolean flag 'ascending' onto 'unindex'.
    fn sortm(self, ascending:bool) -> Vec<f64> {
        self.mergesort(0,self.len()).unindex(ascending,&self)
    }

    /// Ranking of self by inverting the (merge) sort index.  
    /// Sort index is in sorted order, giving indices to the original data positions.
    /// Ranking is in  original data order, giving their positions in the sorted order (sort index).
    /// Thus they are in an inverse relationship, easily converted by `.invindex()`
    /// Fast ranking of many f64 items, ranking `self` with only n*(log(n)+1) complexity.
    fn rank(self, ascending:bool) -> Vec<usize> {
        let n = self.len();
        let sortindex = self.mergesort(0,n);
        let mut rankvec:Vec<usize> = vec![0;n];
        if ascending { 
            for (i,&sortpos) in sortindex.iter().enumerate() {
                rankvec[sortpos] = i
            } 
        } else { // rank in the order of descending values
            for (i,&sortpos) in sortindex.iter().enumerate() {
                    rankvec[sortpos] = n-i-1 
            }
        }
        rankvec 
    }    

    /// Doubly recursive non-destructive merge sort. The data is read-only, it is not moved or mutated. 
    /// Returns vector of indices to self from i to i+n, such that the indexed values are in sort order.  
    /// Thus we are moving only the index (key) values instead of the actual values. 
    fn mergesort(self, i:usize, n:usize) -> Vec<usize> {

        if n == 1 { let res = vec![i]; return res };  // recursion termination
        if n == 2 {  // also terminate with two sorted items (for efficiency)          
            if self[i+1] < self[i] { return vec![i+1,i] } else { return vec![i,i+1] }
        }       
        let n1 = n / 2;  // the first half
        let n2 = n - n1; // the remaining second half
        let sv1 = self.mergesort(i, n1); // recursively sort the first half
        let sv2 = self.mergesort(i+n1, n2); // recursively sort the second half 
        // Now we will merge the two sorted indices into one      
        self.merge_indices(&sv1,&sv2)
    }

    /// New sorted vector. Immutable sort.
    /// Copies self and then sorts it in place, leaving self unchanged.
    /// Calls mutsortf and that calls the standard self.sort_unstable_by.
    /// Consider using our `sortm` instead.
    fn sortf(self) -> Vec<f64> {
        let mut sorted:Vec<f64> = self.to_vec();
        sorted.mutsortf();
        sorted      
    }

    /// Flattened lower triangular part of a covariance matrix for a single f64 vector.
    /// Since covariance matrix is symmetric (positive semi definite), 
    /// the upper triangular part can be trivially added for all j>i by: c(j,i) = c(i,j).
    /// N.b. the indexing is always assumed to be in this order: row,column.
    /// The items of the resulting lower triangular array c[i][j] are pushed flat
    /// into a single vector in this double loop order: left to right, top to bottom 
    fn covone(self, m:&[f64]) -> Vec<f64> {
        let n = self.len(); // dimension of the vector
        let mut cov:Vec<f64> = Vec::new(); // flat lower triangular result array
        let vm = self.vsub(&m); // zero mean vector
        for i in 0..n {
            let thisc = vm[i]; // take this component
            // generate its products up to and including the diagonal (itself)
            for j in 0..i+1 { cov.push(thisc*vm[j]) }
        }
        cov
    }

    /// Reconstructs the full symmetric square matrix from its lower diagonal compact form,
    /// as produced by covar, covone, wcovar
    fn symmatrix(self) -> Vec<Vec<f64>> {
        // solve quadratic equation to find the dimension of the square matrix
        let n = (((8*self.len()+1) as f64).sqrt() as usize - 1)/2;
        let mut mat = vec![vec![0_f64;n];n]; // create the square matrix 
        let mut selfindex = 0;
        for row in 0..n {     
            for column in 0..row { // excludes the diagonal  
                mat[row][column] = self[selfindex]; // just copy the value into the lower triangle
                mat[column][row] = self[selfindex]; // and into transposed upper position 
                selfindex += 1 // move to the next input value
            } // this row of lower triangle finished
            mat[row][row] = self[selfindex];  // now the diagonal element, no transpose
            selfindex += 1 // move to the next input value
        }
        mat
    }
    
}
