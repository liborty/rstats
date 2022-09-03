use crate::{ TriangMat }; // RE, RError, MStats, MinMax, MutVecg, Stats, Vecg, VecVec };
pub use indxvec::{Printing,printing::*};

/// Display implementation for TriangleMat<T>
impl<T> std::fmt::Display for TriangMat<T>
where
    T: Clone + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f, 
            "{YL}{} triangular matrix:{UN}\n{}",
            if self.lower {"Lower"} else {"Upper"}, self.data.gr()
        )
    }
}

/// Implementation of associated functions for struct TriangleMat<T>    
impl<T> TriangMat<T> {
/// Generates nxn identity lower TriangMat matrix
pub fn unit_lower(n:usize) -> TriangMat<f64> {
    let mut data = Vec::new();
    for i in 0..n { 
        for _ in 0..i { data.push(0_f64) };
        data.push(1_f64);
    }
    TriangMat { lower: true, rows: n, data } 
}
/// Generates nxn identity upper TriangMat matrix
pub fn unit_upper(n:usize) -> TriangMat<f64> {
    let mut data = Vec::new();
    for i in 0..n { 
        data.push(1_f64);
        for _ in i+1..n { data.push(0_f64) };
    }
    TriangMat { lower: false, rows: n, data } 
}
}