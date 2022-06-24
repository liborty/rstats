use crate::{MutVecg,MutVecf64};
use indxvec::Vecops; 

/// Mutable vector operations on `&mut [f64]`, where the operand endtype is generic
impl<U> MutVecg<U> for &mut [f64] where U: Copy+PartialOrd, f64: From<U> {
    /// Scalar multiplication of a vector, mutates self
    fn mutsmult(self, s:U) {
        let sf = f64::from(s);
        self.iter_mut().for_each(|x| *x *= sf);
    }

    /// Vector subtraction, mutates self
    fn mutvsub(self, v: &[U]) {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x -= f64::from(vi))
    } 

    /// Vector addition, mutates self
    fn mutvadd(self, v: &[U]) {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x += f64::from(vi))
    } 
}

/// Mutable operations on `&mut [f64]`, where the operand endtype, if any, is also f64 
impl MutVecf64 for &mut [f64] {

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
    
    /// Sorts a mutable slice in place.  
    /// It is the responsibility of the user to ensure that there are no NaNs etc.
    fn msortf(self) {
        self.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap())
    }

    /// Scalar multiplication of a vector, mutates self
    fn mutsmultf64(self, s:f64) {    
        self.iter_mut().for_each(|x| *x *= s);
    }
    /// Vector subtraction, mutates self
    fn mutvsubf64(self, v: &[f64]) {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x -= vi)
    } 
    /// Vector addition, mutates self
    fn mutvaddf64(self, v: &[f64]) {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x += vi)
    } 
}
