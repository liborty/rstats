use crate::{MutVecg,MutVecf64};

impl<U> MutVecg<U> for &mut [f64] where U: Copy+PartialOrd, f64: From<U> {
    /// Scalar multiplication of a vector, mutates self
    fn msmult(self, s:U) {
        let sf = f64::from(s);
        self.iter_mut().for_each(|x| *x *= sf);
    }

    /// Vector subtraction, mutates self
    fn mvsub(self, v: &[U]) {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x -= f64::from(vi))
    } 

    /// Vector addition, mutates self
    fn mvadd(self, v: &[U]) {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x += f64::from(vi))
    } 
}

impl MutVecf64 for &mut [f64] {
    /// Scalar multiplication of a vector, mutates self
    fn mutsmult(self, s:f64) {    
        self.iter_mut().for_each(|x| *x *= s);
    }
    /// Vector subtraction, mutates self
    fn mutvsub(self, v: &[f64]) {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x -= vi)
    } 
    /// Vector addition, mutates self
    fn mutvadd(self, v: &[f64]) {
        self.iter_mut().zip(v).for_each(|(x,&vi)| *x += vi)
    } 
}