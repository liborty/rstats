use crate::{MutVecg,MutVecf64};
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

/// Mutable operations on `&mut [f64]`, where the operand endtype is also f64 
impl MutVecf64 for &mut [f64] {
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
