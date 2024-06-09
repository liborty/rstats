use crate::*; // MStats, MinMax, MutVecg, Stats, VecVec };
pub use indxvec::{Indices, Printing, Vecops};

/// Meanings of 'kind' field. Note that 'Upper Symmetric' would represent the same full matrix as
/// 'Lower Symmetric', so it is not used (lower symmetric matrix is never transposed)
const KINDS: [&str; 5] = [
    "Lower",
    "Lower antisymmetric",
    "Lower symmetric",
    "Upper",
    "Upper antisymmetric",
];

/// Translates single subscript to .data to a pair of  
/// (row,column) coordinates within a lower/upper triangular matrix.
/// Enables memory efficient representation of triangular matrices as one flat vector.
fn rowcol(s: usize) -> (usize, usize) {
        let row = ((((8 * s + 1) as f64).sqrt() - 1.) / 2.) as usize; // cast truncates like .floor()
        let column = s - row * (row + 1) / 2; // subtracting the last triangular number (of whole rows)
        (row, column)
}

/// Display implementation for TriangMat
impl std::fmt::Display for TriangMat {
    fn fmt<'a>(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let n = Self::dim(self);  
        write!(
            f,
            "{} ({n}x{n}) triangular matrix\n{}",
            KINDS[self.kind],
            (0..n).map(|r| self.row(r)).collect::<Vec<Vec<f64>>>().to_str()
        )
    }
}

/// Implementation of associated functions for struct TriangleMat.
/// End type is f64, as triangular matrices will be mostly computed
impl TriangMat {
    /// Length of the data vec
    pub fn len(&self) -> usize {
        self.data.len()
    }
    /// Dimension of the implied full (square) matrix
    /// from the quadratic equation: `n^2 + n - 2l = 0`
    pub fn dim(&self) -> usize {
        ((((8 * self.data.len() + 1) as f64).sqrt() - 1.) / 2.) as usize
    }
    /// Empty TriangMat test
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
    /// Squared euclidian vector magnitude (norm) of the data vector
    pub fn magsq(&self) -> f64 {
        self.data.vmagsq()
    }
    /// Sum of the elements:  
    /// when applied to the wedge product **a∧b**, returns det(**a,b**)
    pub fn sum(&self) -> f64 {
        self.data.iter().sum()
    }
    /// Diagonal elements
    pub fn diagonal(&self) -> Vec<f64> {
        let mut next = 0_usize;
        let mut skip = 1;
        let dat = &self.data;
        let mut diagonal = Vec::with_capacity(self.dim());
        while next < dat.len() {
            diagonal.push(dat[next]);
            skip += 1;
            next += skip;
        }
        diagonal
    }
    /// Determinant of C = LL' is the square of the product of the diagonal elements of L
    pub fn determinant(&self) -> f64 {
        let product = self.diagonal().iter().product::<f64>();
        product * product
    }
    /// New unit (symmetric) TriangMat matrix (data size `n*(n+1)/2`)
    pub fn unit(n: usize) -> Self {
        let mut data = Vec::new();
        for i in 0..n {
            // fill with zeros before the diagonal
            for _ in 0..i {
                data.push(0_f64)
            }
            data.push(1_f64);
        }
        TriangMat { kind: 2, data }
    }

    /// Project symmetric/antisymmetric triangmat to a smaller one of the same kind,
    /// into a subspace specified by an ascending index of dimensions.
    /// Deletes all rows and columns of the missing dimensions.
    pub fn project(&self, index: &[usize]) -> Self {
        let mut res = Vec::with_capacity(sumn(index.len()));
        for &row_idx in index {
            let row = self.row(row_idx);
            for &column_idx in index {
                if column_idx >= row.len() {
                    break;
                };
                res.push(row[column_idx]);
            }
        }
        TriangMat {
            kind: self.kind,
            data: res,
        }
    }
    /// Copy one raw data row from TriangMat
    /// To interpret the kind (plain, symmetric, assymetric, transposed),
    /// use `realrow,realcolumn,to_full`
    pub fn row(&self, r: usize) -> Vec<f64> {
        let idx = sumn(r);
        let Some(slice) = self.data.get(idx..idx + r + 1) else {
            eprintln!("row called with invalid {r}, returned empty Vec");
            return Vec::new();
        };
        slice.to_vec()
    }    
    /// Trivial implicit transposition of a mutable TriangMat.
    /// The untransposed matrix is gone.
    /// To keep the original, use `clone_transpose` below
    pub fn transpose(&mut self) {
        if self.kind != 2 {
            self.kind += 3;
            self.kind %= 6;
        }
    }
    /// Implicit transposition of a cloned TriangMat. 
    pub fn clone_transpose(&self) -> TriangMat { 
        TriangMat { 
            kind: if self.kind == 2 { self.kind } else { (self.kind + 3) % 6 }, 
            data:self.data.clone()
        } 
    } 
    /// One (short) row of a triangular matrix,  
    /// assumed to be zero filled at the end.
    /// When the matrix is transposed (kind>2),  
    /// this will be a (short) column,
    /// assumed to be zero filled upfront.
    pub fn realrow(&self, r: usize) -> Vec<f64> {
        let idx = sumn(r);
        let Some(todiag) = self.data.get(idx..idx + r + 1) else {
            eprintln!("fullrow called with invalid {r}, returned empty Vec");
            return Vec::new();
        };
        let mut rowvec = todiag.to_vec();
        // continue down from the diagonal along its column
        match self.kind % 3 {
            // symmetric
            2 => {
                for row in r + 1..self.dim() {
                    rowvec.push(self.data[sumn(row) + r]);
                }
            }
            // antisymmetric
            1 => {
                for row in r + 1..self.dim() {
                    rowvec.push(-self.data[sumn(row) + r]);
                }
            }
            // neither = plain
            _ => (), // rowvec.resize(self.dim(), 0_f64),
        };
        rowvec
    }
    /// One (short) column of a triangular matrix,
    /// assumed to be zero filled upfront.
    /// When the matrix is transposed (kind>2),  
    /// this will be a (short) row,
    /// assumed to be zero filled at the end.
    pub fn realcolumn(&self, r: usize) -> Vec<f64> {
        let idx = sumn(r);
        // reflect the corresponding row up to diagonal
        let mut columnvec = match self.kind % 3 {
            // symmetric
            2 => self
                .data
                .iter()
                .skip(idx)
                .take(r)
                .copied()
                .collect::<Vec<f64>>(),
            // antisymmetric
            1 => self
                .data
                .iter()
                .skip(idx)
                .take(r)
                .map(|&dataitem| -dataitem)
                .collect::<Vec<f64>>(),
            // neither = plain, fill with zeroes
            _ => Vec::with_capacity(self.dim()), // vec![0_f64; r]
        };
        // now add the column starting below the diagonal
        for row in r..self.dim() {
            columnvec.push(self.data[sumn(row) + r]);
        }
        columnvec
    }
    /// Unpacks all kinds of TriangMat to equivalent full matrix form
    /// For multiplications, use `rmultv,lmultv,mult` instead, to save this unpacking.
    pub fn to_full(&self) -> Vec<Vec<f64>> {
        let n = self.dim();
        if self.kind > 2 {
            // transpose
            (0..self.dim())
                .map(|rownum| {
                    let mut column = vec![0_f64; rownum]; // fill zeroes
                    column.append(&mut self.realcolumn(rownum));
                    column
                })
                .collect::<Vec<Vec<f64>>>()
        } else {
            (0..self.dim())
                .map(|rownum| {
                    let mut shortrow = self.realrow(rownum);
                    shortrow.resize(n, 0_f64); // fill zeroes
                    shortrow
                })
                .collect::<Vec<Vec<f64>>>()
        }
    }
    /// Postmultiply row vector v by triangular matrix `self`.
    /// When a column of self is shorter, it is as if padded with zeroes upfront.
    /// When v is shorter, it is as if padded with zeroes at the end.
    pub fn rmultv<U>(&self, v: &[U]) -> Vec<f64>
    where
        U: Copy + PartialOrd + std::fmt::Display,
        f64: From<U>,
    {
        if self.kind > 2 {
            // transpose
            (0..self.dim())
                .map(|rownum| self.realrow(rownum).dotp(v))
                .collect::<Vec<f64>>()
        } else {
            (0..self.dim())
                .map(|rownum| v.dotp(&self.realcolumn(rownum)))
                .collect::<Vec<f64>>()
        }
    }
    /// Premultiply column vector v by triangular matrix `self`.
    /// When a row of self is shorter, it is as if padded with zeroes at the end.
    /// When v is shorter, it is as if padded with zeroes upfront.
    /// The output is (assumed to be) a column.
    pub fn lmultv<U>(&self, v: &[U]) -> Vec<f64>
    where
        U: Copy + PartialOrd + std::fmt::Display,
        f64: From<U>,
    {
        if self.kind > 2 {
            // transpose
            (0..self.dim())
                .map(|rownum| v.dotp(&self.realcolumn(rownum)))
                .collect::<Vec<f64>>()
        } else {
            (0..self.dim())
                .map(|rownum| self.realrow(rownum).dotp(v))
                .collect::<Vec<f64>>()
        }
    }
    /// One element of a product matrix, used by `mult`
    /// given its precomputed (short) row/column vectors
    /// self is used here only to test its `kind`
    fn dotmult(&self, selfvec: &[f64], otvec: &[f64], otherkind: usize) -> f64 {
        if self.kind > 2 {
            if otherkind > 2 {
                otvec.dotp(selfvec)
            } else if selfvec.len() > otvec.len() {
                selfvec.dotp(otvec)
            } else {
                otvec.dotp(selfvec)
            }
        } else if otherkind > 2 {
            if selfvec.len() > otvec.len() {
                otvec.dotp(selfvec)
            } else {
                selfvec.dotp(otvec)
            }
        } else {
            selfvec.dotp(otvec)
        }
    }
    /// General multiplication of two triangular matrices (of any kind).
    /// The triangular matrices are not expanded and
    /// incomplete rows/columns are not even padded (very effient).
    pub fn mult(&self, other: &Self) -> Vec<Vec<f64>> {
        (0..self.dim())
            .map(|rownum| {
                let selfvec = if self.kind > 2 {
                    self.realcolumn(rownum)
                } else {
                    self.realrow(rownum)
                };
                (0..other.dim())
                    .map(|colnum| {
                        let otvec = if other.kind > 2 {
                            other.realrow(colnum)
                        } else {
                            other.realcolumn(colnum)
                        };
                        self.dotmult(&selfvec, &otvec, other.kind)
                    })
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<Vec<f64>>>()
    }
    /// Efficient Cholesky-Banachiewicz matrix decomposition into `LL'`,
    /// where L is the returned lower triangular matrix and L' its upper triangular transpose.
    /// Expects as input a positive definite matrix
    /// in TriangMat compact form, such as a covariance matrix produced by `covar`.
    /// The computations are all done on the compact form,
    /// making this implementation memory efficient for large (symmetric) matrices.
    /// Reports errors if the above conditions are not satisfied.
    pub fn cholesky(&self) -> Result<Self, RE> {
        let sl = self.data.len();
        // input not long enough to compute anything
        if sl < 3 {
            return nodata_error("cholesky needs at least 3x3 TriangMat");
        };
        // n is the dimension of the implied square matrix.
        // Not needed as an extra argument. We compute it
        // by solving a quadratic equation in seqtosubs()
        let (n, c) = rowcol(sl);
        // input is not a triangular number, is of wrong size
        if c != 0 {
            return data_error("cholesky needs a triangular matrix");
        };
        let mut res = vec![0.0; sl]; // result L is of the same size as the input
        for i in 0..n {
            let isub = i * (i + 1) / 2; // matrix row index to the compact vector index
            for j in 0..(i + 1) {
                // i+1 to include the diagonal
                let jsub = j * (j + 1) / 2; // matrix column index to the compact vector index
                let mut sum = 0.0;
                for k in 0..j {
                    sum += res[isub + k] * res[jsub + k];
                }
                let dif = self.data[isub + j] - sum;
                res[isub + j] = if i == j {
                    // diagonal elements
                    // dif <= 0 means that the input matrix is not positive definite,
                    // or is ill-conditioned, so we return ArithError
                    if dif <= 0_f64 {
                        return arith_error("cholesky matrix is not positive definite");
                    };
                    dif.sqrt()
                }
                // passed, so enter real non-zero square root
                else {
                    dif / res[jsub + j]
                };
            }
        }
        Ok(TriangMat { kind: 0, data: res })
    }

    /// Mahalanobis scaled magnitude m(d) of a (column) vector d.
    /// Self is the decomposed lower triangular matrix L, as returned by `cholesky`
    /// decomposition of covariance/comediance positive definite matrix: C = LL',
    /// where ' denotes transpose. Mahalanobis distance is defined as:    
    /// `m(d) = sqrt(d'inv(C)d) = sqrt(d'inv(LL')d) = sqrt(d'inv(L')inv(L)d)`,
    /// where `inv()` denotes matrix inverse, which is never explicitly computed.  
    /// Let  `x = inv(L)d` ( and therefore also  `x' = d'inv(L')` ).
    /// Substituting x into the above definition: `m(d) = sqrt(x'x) = |x|.  
    /// We obtain x by setting Lx = d and solving by forward substitution.
    /// All the calculations are done in the compact triangular form.
    pub fn mahalanobis<U>(&self, d: &[U]) -> Result<f64, RE>
    where
        U: Copy + PartialOrd + std::fmt::Display,
        f64: From<U>,
    {
        Ok(self.forward_substitute(d)?.vmag())
    }

    /// Solves for x the system of linear equations Lx = b,
    /// where L (self) is a lower triangular matrix.   
    pub fn forward_substitute<U>(&self, b: &[U]) -> Result<Vec<f64>, RE>
    where
        U: Copy + PartialOrd + std::fmt::Display,
        f64: From<U>
    {
        if self.kind != 0 { return data_error("forward-substitute expects plain lower kind"); };
        let data = &self.data;
        if data.len() < 3 {
            return Err(RError::NoDataError(
                "forward-substitute needs at least three items".to_owned(),
            ));
        };
        // 2d matrix dimension
        let n = self.dim();
        // dimensions/lengths mismatch
        if n != b.len() {
            return Err(RError::DataError(
                "forward_substitute mismatch of self and b dimension".to_owned(),
            ));
        };
        let mut res: Vec<f64> = Vec::with_capacity(n); // result of the same size and shape as b
        res.push(f64::from(b[0]) / self.data[0]);
        for (row, &b_component) in b.iter().enumerate().take(n).skip(1) {
            let rowoffset = sumn(row);
            let mut sumtodiag = 0_f64;
            for (column, res_component) in res.iter().enumerate() {
                sumtodiag += self.data[rowoffset + column] * res_component;
            }
            res.push((f64::from(b_component) - sumtodiag) / self.data[rowoffset + row]);
        }
        // println!("Forward substitution: {}",res.gr());
        Ok(res)
    }

    /// Householder's Q*M matrix product without explicitly computing Q
    pub fn house_uapply<T>(&self, m: &[Vec<T>]) -> Vec<Vec<f64>>
    where
        T: Copy + PartialOrd + std::fmt::Display,
        f64: From<T>,
    {
        let u = self.to_full();
        let mut qm = m.iter().map(|mvec| mvec.tof64()).collect::<Vec<Vec<f64>>>();
        for uvec in u.iter().take(self.dim()) {
            qm.iter_mut()
                .for_each(|qvec| *qvec = uvec.house_reflect::<f64>(qvec))
        }
        qm
    }
}
