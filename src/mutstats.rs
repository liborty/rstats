use crate::{MutStats};

pub use indxvec::merge::{sortm,minmax};          

impl MutStats for &mut [f64] {  

    /// Vector with inverse magnitude
    fn minvert(self) {
        let recmag2 = 1.0/self.iter().map(|&x| f64::from(x).powi(2)).sum::<f64>();
        for c in self.iter_mut() { *c *= recmag2 }    
    }

    // negated vector (all components swap sign)
    fn mneg(self) { 
        for c in self.iter_mut() { *c *= -1_f64 }
    }

    /// Unit vector
    fn munit(self) {
        let m = (1.0 / self.iter().map(|&x| x.powi(2)).sum::<f64>()).sqrt();
        for c in self.iter_mut() { *c *= m } 
    }

    /// Linear transform to interval [0,1]
    fn mlintrans(self) {
        let mm = minmax(self);
        let range = mm.max-mm.min;
        for c in self.iter_mut() { *c = (*c-mm.min)/range }        
    }

    /// Sorts a mutable slice in place.  
    /// It is the responsibility of the user to ensure that there are no NaNs etc.
    fn msortf(self) {
        self.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap())
    }  
}
