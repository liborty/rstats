use crate::{sumn, RError, Stats, TriangMat, RE}; // RE, RError, MStats, MinMax, MutVecg, Stats, Vecg, VecVec };
pub use indxvec::{printing::*, Printing};

/// Display implementation for TriangMat
/// Converts to, and displays, the full matrix form
impl std::fmt::Display for TriangMat {
    fn fmt<'a>(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let dim = Self::rowcol(self.data.len()).0;
        write!(
            f,
            "{YL}{} triangular {} matrix {dim}x{dim}:\n{}",
            if self.lower { "Lower" } else { "Upper" },
            if self.symmetric {
                "symmetric"
            } else {
                "non-symmetric"
            },
            self.to_full().gr()
        )
    }
}

/// Implementation of associated functions for struct TriangleMat.
/// End type is f64 here, as triangular matrices will be mostly computed,
/// so we may as well keep it simple.    
impl TriangMat {
    /// Generates new unit TriangMat matrix of size (n+1)*n/2
    pub fn unit(n: usize, lower: bool) -> TriangMat {
        let mut data = Vec::new();
        for i in 0..n {
            // fill with zeros before the diagonal
            for _ in 0..i {
                data.push(0_f64)
            }
            data.push(1_f64);
        }
        TriangMat {
            lower,
            symmetric: true,
            data,
        }
    }
    /// Translates subscripts to a 1d vector, i.e. natural numbers, to a pair of
    /// full coordinates within a lower/upper triangular matrix.
    /// Enables memory efficient representation of triangular matrices as a flat vector.
    pub fn rowcol(s: usize) -> (usize, usize) {
        let row = ((((8 * s + 1) as f64).sqrt() - 1.) / 2.) as usize; // cast truncates, like .floor()
        let column = s - row * (row + 1) / 2; // subtracting the previous triangular number
        (row, column)
    }

    fn to_full(&self) -> Vec<Vec<f64>> {
        let (n, _) = TriangMat::rowcol(self.data.len());
        let mut res = vec![vec!(0_f64; n); n];
        for (i, &item) in self.data.iter().enumerate() {
            let (row, col) = TriangMat::rowcol(i);
            if row == col {
                res[row][col] = item;
                continue;
            }; //diagonal
            if self.lower || self.symmetric {
                res[row][col] = item;
            };
            if !self.lower || self.symmetric {
                res[col][row] = item;
            };
        }
        res
    }

    /// Efficient Cholesky-Banachiewicz matrix decomposition into `LL'`,
    /// where L is the returned lower triangular matrix and L' its upper triangular transpose.
    /// Expects as input a symmetric positive definite matrix
    /// in TriangMat compact form, such as a covariance matrix produced by `covar`.
    /// The computations are all done on the compact form,
    /// making this implementation memory efficient for large (symmetric) matrices.
    /// Reports errors if the above conditions are not satisfied.
    pub fn cholesky(self) -> Result<TriangMat, RE> {
        let sl = self.data.len();
        // input not long enough to compute anything
        if sl < 3 {
            return Err(RError::NoDataError(
                "cholesky needs at least three Vec items",
            ));
        };
        // n is the dimension of the implied square matrix.
        // Not needed as an extra argument. We compute it
        // by solving a quadratic equation in seqtosubs()
        let (n, c) = TriangMat::rowcol(sl);
        // input is not a triangular number, is of wrong size
        if c != 0 {
            return Err(RError::DataError("cholesky needs a triangular matrix"));
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
                        return Err(RError::ArithError(
                            "cholesky needs a positive definite matrix",
                        ));
                    };
                    dif.sqrt()
                }
                // passed, so enter real non-zero square root
                else {
                    dif / res[jsub + j]
                };
            }
        }
        Ok(TriangMat {
            lower: true,
            symmetric: false,
            data: res,
        })
    }

    /// Mahalanobis scaled magnitude M(d) of vector d:
    /// self is a precomputed lower triagonal matrix L, as returned by `cholesky`
    /// from covariance/comediance positive definite matrix C = LL'.
    /// M(d) = sqrt(d'inv(C)d) = sqrt(d'inv(LL')d) = sqrt(d'inv(L')inv(L)d),
    /// where ' denotes transpose and inv denotes inverse.
    /// Putting Lx = d and solving for x by forward substitution, we obtain x = inv(L)d
    ///  => M(d) = sqrt(x'x) = |x|.
    /// We stay in the compact triangular, from C to M(d).
    pub fn mahalanobis<U>(self, d: &[U]) -> Result<f64, RE>
    where
        U: Copy + PartialOrd + std::fmt::Display,
        f64: From<U>,
    {
        Ok(self.forward_substitute(d)?.vmag())
    }

    /// Solves the system of linear equations Lx = b,
    /// where L (self) is a lower triangular matrix in left to right 1d scan form   
    fn forward_substitute<U>(self, b: &[U]) -> Result<Vec<f64>, RError<&'static str>>
    where
        U: Copy + PartialOrd + std::fmt::Display,
        f64: From<U>,
    {
        let sl = self.data.len();
        if sl < 3 {
            return Err(RError::NoDataError(
                "forward-substitute needs at least three items",
            ));
        };
        // 2d matrix dimensions
        let (n, c) = TriangMat::rowcol(sl);
        if c != 0 {
            return Err(RError::DataError(
                "forward_substitute needs a triangular matrix",
            ));
        };
        // dimensions/lengths mismatch
        if n != b.len() {
            return Err(RError::DataError(
                "forward_substitute mismatch of self and b dimension",
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
