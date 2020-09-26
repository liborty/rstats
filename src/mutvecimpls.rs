use crate::{MutVectors};

impl MutVectors for &mut [f64] {
    /// Scalar multiplication of a vector, mutates self
    fn mutsmult(self, s: f64) {
        self.iter_mut().for_each(|x| *x *= s);
    }
    /// Vector subtraction, mutates self
    fn mutvsub(self, v: &[f64]) {
        self.iter_mut().enumerate().for_each(|(i, x)| *x -= v[i])
    }
    /// Vector addition, mutates self
    fn mutvadd(self, v: &[f64]) {
        self.iter_mut().enumerate().for_each(|(i, x)| *x += v[i])
    }
    /// Mutate to unit vector
    fn mutvunit(self) {
        self.mutsmult(1_f64 / self.iter().map(|x| x.powi(2)).sum::<f64>().sqrt())
    }
    /// Sorts a mutable `Vec<f64>` in place.  
    /// It is the responsibility of the user to ensure that there are no NaNs etc.
    fn mutsortf(self) {
        self.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap())
    } 
    /// Vector subtraction, mutates self
    fn mutvsubu8(self, v: &[u8]) {
        self.iter_mut().enumerate().for_each(|(i, x)| *x -= v[i] as f64)
    }
    /// Vector addition, mutates self
    fn mutvaddu8(self, v: &[u8]) {
        self.iter_mut().enumerate().for_each(|(i, x)| *x += v[i] as f64)
    }  
}
