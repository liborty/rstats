use crate::{MutVecg};
use indxvec::Vecops; 

/// Mutable vector operations on `&mut [f64]`, where the operand endtype is generic
impl MutVecg for &mut [f64] {
    /// Scalar multiplication of a vector, mutates self
    fn mutsmult<U>(self, s:U) where U: Copy+PartialOrd, f64: From<U> {
        let sf = f64::from(s);
        self.iter_mut().for_each(|x| *x *= sf);
    }

    /// Vector subtraction, mutates self
    fn mutvsub<U>(self, v: &[U]) where U: Copy+PartialOrd, f64: From<U> {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x -= f64::from(vi))
    } 

    /// Vector addition, mutates self
    fn mutvadd<U>(self, v: &[U]) where U: Copy+PartialOrd, f64: From<U> {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x += f64::from(vi))
    } 

    /// Vector with inverse magnitude
    fn minvert(self) {
        let recmag2 = 1.0/self.iter().map(|&x| x.powi(2)).sum::<f64>();
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
        let mm = self.minmax();
        let range = mm.max-mm.min;
        for c in self.iter_mut() { *c = (*c-mm.min)/range }        
    } 
}
