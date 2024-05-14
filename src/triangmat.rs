use crate::*; // MStats, MinMax, MutVecg, Stats, VecVec };
pub use indxvec::{printing::*, Indices, Printing, Vecops};

/// Meanings of 'kind' field. Note that 'Upper Symmetric' would represent the same full matrix as
/// 'Lower Symmetric', so it is not used (lower symmetric matrix is never transposed)
const KINDS: [&str; 5] = [
    "Lower",
    "Lower antisymmetric",
    "Lower symmetric",
    "Upper",
    "Upper antisymmetric",
];

/// Display implementation for TriangMat
impl std::fmt::Display for TriangMat {
    fn fmt<'a>(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let dim = Self::dim(self);
        write!(
            f,
            "{YL}{} ({dim}x{dim}) triangular matrix:\n{}",
            KINDS[self.kind],
            self.to_triangle().gr()
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
        };
        diagonal
    }
    /// Determinant of A = LL' is the square of the product of the diagonal elements.
    /// However, the diagonal elements squared are not the eigenvalues of A!
    pub fn determinant(&self) -> f64 {
        let mut next = 0_usize;
        let mut skip = 1;
        let mut product = 1_f64;
        let dat = &self.data; 
        while next < dat.len() {
            product *= dat[next];
            skip += 1;
            next += skip;  
        };
        product*product
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
    /// Translates subscripts to a 1d vector, i.e. natural numbers, to a pair of
    /// (row,column) coordinates within a lower/upper triangular matrix.
    /// Enables memory efficient representation of triangular matrices as one flat vector.
    pub fn rowcol(s: usize) -> (usize, usize) {
        let row = ((((8 * s + 1) as f64).sqrt() - 1.) / 2.) as usize; // cast truncates like .floor()
        let column = s - row * (row + 1) / 2; // subtracting the last triangular number (of whole rows)
        (row, column)
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
    /// Extract one row from TriangMat
    pub fn row(&self, r: usize) -> Vec<f64> {
        let idx = sumn(r);
        self.data.get(idx..idx + r + 1).unwrap().to_vec()
    }
    /// Unpacks flat TriangMat Vec to triangular Vec<Vec> form
    pub fn to_triangle(&self) -> Vec<Vec<f64>> {
        let n = self.dim();
        let mut res = Vec::with_capacity(n);
        for r in 0..n {
            res.push(self.row(r));
        }
        res
    }
    /// TriangMat trivial implicit transposition
    pub fn transpose(&mut self) {
        if self.kind != 2 {
            self.kind += 3;
            self.kind %= 6;
        }
    }
    /// Unpacks all kinds of TriangMat to equivalent full matrix form
    pub fn to_full(&self) -> Vec<Vec<f64>> {
        // full matrix dimension(s)
        let n = self.dim();
        let mut res = vec![vec!(0_f64; n); n];
        // function pointer for primitive filling actions, depending on the matrix kind
        let fill: fn(usize, usize, &mut Vec<Vec<f64>>, f64) = match self.kind % 3 {
            // symmetric
            2 => |row: usize, col: usize, res: &mut Vec<Vec<f64>>, item: f64| {
                res[row][col] = item;
                if row != col {
                    res[col][row] = item;
                };
            },
            // antisymmetric
            1 => |row: usize, col: usize, res: &mut Vec<Vec<f64>>, item: f64| {
                res[row][col] = item;
                if row != col {
                    res[col][row] = -item;
                };
            },
            // plain
            _ => |row: usize, col: usize, res: &mut Vec<Vec<f64>>, item: f64| {
                res[row][col] = item;
                if row != col {
                    res[col][row] = 0_f64;
                };
            },
        };
        if self.kind > 2 {
            // transpose
            for (i, &item) in self.data.iter().enumerate() {
                let (row, col) = Self::rowcol(i);
                fill(col, row, &mut res, item);
            }
        } else {
            // do not transpose
            for (i, &item) in self.data.iter().enumerate() {
                let (row, col) = Self::rowcol(i);
                fill(row, col, &mut res, item);
            }
        };
        res
    }

    /// Efficient Cholesky-Banachiewicz matrix decomposition into `LL'`,
    /// where L is the returned lower triangular matrix and L' its upper triangular transpose.
    /// Expects as input a symmetric positive definite matrix
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
        let (n, c) = TriangMat::rowcol(sl);
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
    fn forward_substitute<U>(&self, b: &[U]) -> Result<Vec<f64>, RE>
    where
        U: Copy + PartialOrd + std::fmt::Display,
        f64: From<U>,
    {
        let sl = self.data.len();
        if sl < 3 {
            return Err(RError::NoDataError(
                "forward-substitute needs at least three items".to_owned(),
            ));
        };
        // 2d matrix dimensions
        let (n, c) = TriangMat::rowcol(sl);
        if c != 0 {
            return Err(RError::DataError(
                "forward_substitute needs a triangular matrix".to_owned(),
            ));
        };
        // dimensions/lengths mismatch
        if n != b.len() {
            return Err(RError::DataError(
                "forward_substitute mismatch of self and b dimension".to_owned(),
            ));
        };
        let mut res: Vec<f64> = Vec::with_capacity(n); // result of the same size and shape as b
        res.push(f64::from(b[0]) / self.data[0]);
        for (row, &bitem) in b.iter().enumerate().take(n).skip(1) {
            let mut sumtodiag = 0_f64;
            let rowoffset = sumn(row);
            for (column, resc) in res.iter().enumerate().take(row) {
                sumtodiag += self.data[rowoffset + column] * resc;
            }
            res.push((f64::from(bitem) - sumtodiag) / self.data[rowoffset + row]);
        }
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

    /* Leftmultiply (column) vector v by upper triangular matrix self
    fn utriangmultv<U>(self,v: &[U]) -> Result<Vec<f64>,RE>
        where U: Copy+PartialOrd+std::fmt::Display, f64:From<U> {
        let sl = self.data.len();
        if sl < 1 { return Err(RError::NoDataError("utriangmultv needs at least one item"));};
        // 2d matrix dimensions
        let (n,c) = TriangMat::rowcol(sl);
        if c != sl { return Err(RError::DataError("utriangmultv expects a triangular matrix"));};
        if n != v.len() { return Err(RError::DataError("utriangmultv dimensions mismatch")); };
        let mut res:Vec<f64> = vec![0_f64;n];
        for row in 0..n {
            for j in row..n {
                res[row] += self.data[sumn(row)+j]*f64::from(v[j])
            };
        };
        Ok(res)
    }
    */
}
