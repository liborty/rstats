# Rstats [<img alt="crates.io" src="https://img.shields.io/crates/v/Rstats?logo=rust">](https://crates.io/crates/rstats) [<img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/liborty/Rstats/HEAD?logo=github">](https://github.com/liborty/Rstats) [![Actions Status](https://github.com/liborty/rstats/actions/workflows/tests.yml/badge.svg)](https://github.com/liborty/rstats/actions)

Statistics, Information Measures, Vector Algebra, Linear Algebra, Cholesky Matrix Decomposition, Mahalanobis Distance, Householder QR Decomposition, Multidimensional Data Analysis, Geometric Median, Hulls, Machine Learning ...

## Usage

Insert `Rstats = "^1"` in the `Cargo.toml` file, under `[dependencies]`.

Use in your source files any of the following structs, as and when needed:

```rust  
use Rstats::{RE,TriangMat,Mstats,MinMax,Med};
```

and any of the following traits:

```rust 
use Rstats::{Stats,Vecg,Vecu8,MutVecg,VecVec,VecVecg};
```
and any of the following auxiliary functions:

```rust
use Rstats::{noop,fromop,sumn,st_error,unit_matrix};
```

The latest (nightly) version is always available in the github repository [Rstats](https://github.com/liborty/Rstats). Sometimes it may be only in some details a little ahead of the crates.io release versions.

It is highly recommended to read and run [tests.rs](https://github.com/liborty/Rstats/blob/master/tests/tests.rs) for examples of usage. To run all the tests, use a single thread in order to print the results in the right order:

```bash  
cargo test --release -- --test-threads=1 --nocapture --color always
```

Alternatively, just to get a quick idea of the methods provided and their usage, you can now read the output produced by an [automated test run](https://github.com/liborty/rstats/actions). There are test logs generated for each new push to the github repository. Click the latest (top) one, then `Rstats` and then `Run cargo test` ... The badge at the top of this document lights up green when all the tests have passed and clicking it gets you to these logs as well.

Any compilation errors arising out of `rstats` crate indicate most likely that some of the dependencies have become out of date. Issuing `cargo update` command will usually fix this.

## Introduction

`Rstats` has a small footprint. Only the best methods are implemented, primarily with Data Analysis and Machine Learning in mind. They include multidimensional `nd` analysis, i.e. characterising sets of n points in space of d dimensions.

Several branches of mathematics: statistics, information theory, set theory and linear algebra are combined in this one consistent crate, based on the abstraction that they all operate on the same data objects (here Rust Vecs). The only difference being that an ordering of their components is sometimes assumed (in linear algebra, set theory) and sometimes it is not (in statistics, information theory, set theory).

`Rstats` begins with basic statistical measures, information measures, vector algebra and linear algebra. These provide self-contained tools for the multidimensional algorithms but are also useful in their own right.

`Non analytical (non parametric) statistics` is preferred, whereby the 'random variables' are replaced by vectors of real data. Probabilities densities and other parameters are always obtained from the real data, not from some assumed distributions.

`Linear algebra` uses by default `Vec<Vec<T>>`, generic data structure capable of representing irregular matrices. Also, `struct TriangMat` is defined and used for symmetric and triangular matrices (for memory efficiency reasons).

Our treatment of multidimensional sets of points is constructed from the first principles. Some original concepts, not found elsewhere, are introduced and implemented here:

* `median correlation`- in one dimension, our `mediancorr` method is to replace `Pearson's correlation`. We define `median correlation` as the cosine of an angle between two zero median samples (instead of Pearson's zero mean samples). This conceptual clarity is one of the benefits of interpreting a data sample of length d as a single point in d dimensional space (or vector).

* `gmedian and pmedian` - fast multidimensional `geometric median (gm)` algorithms.

* `madgm` - our generalization to n dimensions of a robust data spread estimator known as `MAD` (median of absolute deviations from median).

* `standard error` - also generalized to `nd`. Here the role of the central tendency is taken by the `geometric median` and the spread by `madgm`. Thus a single scalar standard error is obtained in any number of dimensions.

* `contribution` - of a point to an `nd` set. Defined as a magnitude of `gm` adjustment caused by adding the point to the set. It is related to the point's radius (distance from the `gm`) but is not the same, as it depends on the radii of all the other points as well.

* `comediance` - instead of covariance (triangular matrix). It is obtained by supplying `covar` with the geometric median instead of the usual centroid. Thus `zero median vectors` are replacing `zero mean vectors` in covariance calculations. The results are similar but more stable with respect to the outliers.

*Zero median vectors are generally preferable to the commonly used zero mean vectors.*

In n dimensions, many authors 'cheat' by using `quasi medians` (1-d medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. In a highly dimensional space, they are also slower to compute than is our `gm`.

*Specifically, all such 1d measures are sensitive to the choice of axis and thus are affected by their rotation.*

In contrast, analyses based on the true geometric median (`gm`) are axis (rotation) independent. Also, they are more stable, as medians have a 50% breakdown point (the maximum possible). They are computed here by methods `gmedian` and its weighted version `wgmedian`, in traits `vecvec` and `vecvecg` respectively.

## Additional Documentation

For more detailed comments, plus some examples, see [docs.rs](https://docs.rs/rstats/latest/rstats). You may have to go directly to the modules source. These traits are implemented for  'out of this crate' rust `Vec` type and rust docs do not display 'implementations on foreign types' very well.

## Terminology

#### Including some new definitions for sets of nd points, i.e. n points in d dimensional space

* `Median correlation` between
 two samples. We define it analogously to Pearson, as cosine of an angle between two 'normalised' vectors. Pearson 'normalises' by subtracting the mean from all components, we subtract the median.

* `Centroid/Centre/Mean` of an 'nd' set. It is the (generally non member) point that minimises its sum of *squares* of distances to all member points. Thus it is susceptible to outliers. Specifically, it is the d-dimensional arithmetic mean. It is sometimes called 'the centre of mass'. Centroid can also sometimes mean the member of the set which is the nearest to the Centre. Here we follow the common usage: Centroid = Centre = Arithmetic Mean.

* `Quasi/Marginal Median` is the point minimising sums of distances separately in each dimension (its coordinates are  medians along each axis). It is a mistaken concept which we do not recommend using.

* `Tukey Median` is the point maximising `Tukey's Depth`, which is the minimum number of (outlying) points found in a hemisphere in any direction. Potentially useful concept but only partially implemented here by `tukeyvec`, as its advantages over the geometric median are not clear.

* `Median or the true geometric median (gm)`, is the point (generally non member), which minimises the sum of distances to all member points. This is the one we want. It is much less susceptible to outliers than the centroid. In addition, unlike quasi median, `gm` is rotation independent.

* `Medoid` is the member of the set with the least sum of distances to all other members. Equivalently, the member which is the nearest to the `gm`.

* `Outlier` is the member of the set with the greatest sum of distances to all other members. Equivalently, it is the point furthest from the `gm`.

* `Outer Hull` is a subset containing zero median member points p, such that no other points lie outside the normal plane through p. The points that do not satisfy this condition are the `internal` points.

* `Inner Hull or Core` is a subset containing zero median member points p, such that none of the points lie outside (the normal planes through) any other.
Note that in a highly dimensional space up to all points may belong to both hulls.

* `Zero median vectors` are obtained by subtracting the `gm` (placing the origin of the coordinate system at the `gm`). This is a proposed  alternative to the commonly used `zero mean vectors`, obtained by subtracting the centroid.

* `MADGM` (median of distances from `gm`). This is our generalisation of `MAD` (median of absolute differences) measure from one dimension to any number of dimensions (d>1). It is a robust measure of `nd` data spread.

* `Comediance` is similar to `covariance`, except that zero median vectors are used to compute it (instead of zero mean vectors).

* `Mahalanobis Distance` is a weighted distace, where the weights are derived from the axis of variation of the `nd` data points cloud. Thus distances in the directions in which there are few points are penalised (increased) and vice versa. Efficient Cholesky singular (eigen) value decomposition is used. Cholesky method decomposes the covariance or comediance positive definite matrix S into a product of two triangular matrices: S = LL'. For more details see the comments in the source code.

* `Contribution` One of the key questions of Machine Learning (ML) is how to quantify the contribution that each example point (typically a member of some large `nd` set) makes to the recognition concept, or class, represented by that set. In answer to this, we define the `contribution` of a point as the magnitude of adjustment to `gm` caused by adding that point. Generally, outlying points make greater contributions to the `gm` but not as much as they would to the centroid. The contribution depends not only on the radius of the example point in question but also on the radii of all other existing example points.

* `Tukey Vector` Proportions of points in each hemisphere from gm. This is a useful 'signature' of a data cloud. For a new point (that typically needs to be classified) we can then quickly determine whether it lies in a well populated direction. This could also be done by projecting all the existing points on its unit radius vector but that would be much slower, as there are many points. Also, in keeping with the stability properties of medians, we are only using counts of points in the hemispheres, not their distances.


## Implementation

The main constituent parts of Rstats are its traits. The selection of traits (to `use`) is primarily determined by the types of objects handled. These are mostly vectors of arbitrary length/dimensionality (`d`). The main traits are implementing methods applicable to:

* `Stats`: a single vector (of numbers),
* `Vecg`: methods (of vector algebra) operating on two vectors, e.g. scalar product
* `Vecu8`: some specialized methods for end-type `u8`
* `MutVecg`: some of the above methods, mutating self
* `VecVec`: methods operating on n vectors, 
* `VecVecg`: methods for n vectors, plus another generic argument, e.g. vector of weights.

In other words, the traits and their methods operate on arguments of their required categories. In classical statistical terminology, the main categories correspond to the number of 'random variables'.

`Vec<Vec<T>>` type is used for full rectangular matrices (could also be irregular), whereas `TriangMat` struct is used specifically for symmetric and triangular matrices (to save memory).

The vectors' end types (of the actual data) are mostly generic: usually some numeric type. End type `f64` is mostly used for the computed results.

## Errors

`Rstats` crate produces custom errors `RError`:

```rust
pub enum RError<T> where T:Sized+Debug {
    /// Insufficient data
    NoDataError(T),
    /// Wrong kind/size of data
    DataError(T),
    /// Invalid result, such as prevented division by zero
    ArithError(T),
    /// Other error converted to RError
    OtherError(T)
}
```
Each of its enum variants also carries a generic payload `T`. Most commonly this will be a `String` message giving more helpful explanation, e.g.:

```rust 
if dif <= 0_f64 {
    return Err(RError::ArithError(format!(
        "cholesky needs a positive definite matrix {}", dif )));
};
```

`format!(...)` is used to insert values of variables to the payload String, as shown. These potential errors are returned and can then be automatically converted (with `?`) to users' own errors. Some such conversions are implemented at the bottom of `errors.rs` file and used in `tests.rs`.

 There is a type alias shortening return declarations to, e.g.: `Result<Vec<f64>,RE>`, where

 ```rust
pub type RE = RError<String>;
```

## Structs

### `struct MStats` 
holds the central tendency of `1d` data, e.g. some kind of mean or median, and its dispersion measure, e.g. standard deviation or MAD.

### `struct TriangMat` 
holds lower/upper triangular symmetric/non-symmetric matrix in compact form that avoids zeros and duplications. Beyond the usual conversion to full matrix form, a number of (the best) Linear Algebra methods are implemented directly on `TriangMat`, in module `triangmat.rs`, such as:

* **Cholesky-Banachiewicz** matrix decomposition: M = LL' (where ' denotes a transpose). This decomposition is used by `mahalanobis`.
* **Mahalanobis Distance**
* **Householder UR** (M = QR) matrix decomposition

Also, some methods implemented for `VecVecg` produce `TriangMat` matrices, specifically the covariance/comedience calculations: `covar`,`wcovar`,`comed` and `wcomed`. Their results are positive definite. Whenever this condition is satisfied, then the most efficient Cholesky-Banachiewics decomposition is applicable.

##  Quantify Functions

Most methods in `medians::Median` trait and hashort methods in `indxvec` crate require explicit closure to tell them how to quantify into f64 user data of any end type T. Variety of different quantifying methods can then be dynamically employed. For example, in analysis of words (&str type), it can be the word length, or the numerical value of its first few bytes, etc. Then we can sort them or find their means/medians/dispersions under these different measures. We do not necessarily want to explicitly store all such quantifications, as data can be voluminous. Rather, we want to be able to compute them on demand.

### `noop`

is a shorthand dummy function to supply to these methods, when the data is already of f64 end type. The second line is the full equivalent version that can be used instead:

```rust
&mut noop
&mut |f:&f64| *f
```

When T is a primitive type, such as i64, u64, usize, that can only be converted to f64 by explicit truncation, use:

```rust
&mut |f:&T| *f as f64
```

### `fromop`

When T is a type convertible by an existing `From` implementation and `f64:From<T>` has been duly added everywhere as a trait bound, then you can pass in one of these: 

```rust
&mut fromop
&mut |f:&T| f64::from(*f)
```

This also works for 'smaller' primitive types.

All other cases were previously only possible with manual implementation written for the (global) From trait for each type T and each different conversion method, whereby the different conversions would conflict. Now the user can simply pass in a custom 'quantify' closure. This generality is obtained at the price of one small inconvenience: using the above closures for the simple cases.

## Auxiliary Functions

* `sumn`: sum of the sequence `1..n = n*(n+1)/2`. It is also the size of a lower/upper triangular matrix.
* `st_error`: standard error of a value, on the basis of any central tendency and dispersion measures.
* `unit_matrix`: - generates full unit matrix.

## Trait Stats

One dimensional statistical measures implemented for all numeric end types.

Its methods operate on one slice of generic data and take no arguments.
For example, `s.amean()` returns the arithmetic mean of the data in slice `s`.
These methods are checked and will report RError(s), such as an empty input. This means you have to apply `?` to their results to pass the errors up, or explicitly match them to take recovery actions, depending on the error variant.

Included in this trait are:

* means (arithmetic, geometric and harmonic),
* standard deviations,
* linearly weighted means (useful for time analysis),
* probability density function (pdf)
* autocorrelation, entropy
* linear transformation to [0,1],
* other measures and vector algebra operators

Note that fast implementation of 1d medians is, as of version 1.1.0, performed by crate `medians`.  


## Trait Vecg

Generic vector algebra operations between two slices `&[T]`, `&[U]` of any (common) length  (dimensions). Note that it may be necessary to invoke some using the 'turbofish' `::<type>` syntax to indicate the type U of the supplied argument, e.g.: 
`datavec.methodname::<f64>(arg)`. 
This is because Rust is currently incapable of inferring its type ('the inference bug').

Methods implemented by this trait:

* Vector additions, subtractions and products (scalar, kronecker, outer),
* Other relationships and measures of difference,
* Pearson's, Spearman's and Kendall's correlations,
* `Median correlation`, which we define analogously to Pearson's, as cosine of an angle between two zero median vectors (instead of his zero mean vectors).
* Joint pdf, joint entropy, statistical independence (based on mutual information).
* `Contribution` measure of a point to geometric median

The simpler methods of this trait are sometimes unchecked (for speed), so some caution with data is advisable.

## Trait MutVecg

A select few of the `Stats` and `Vecg` methods (e.g. mutable vector addition, subtraction and multiplication) are reimplemented under this trait, so that they can mutate `self` in-place. This is more efficient and convenient in some circumstances, such as in vector iterative methods.

## Trait Vecu8

Some vector algebra as above that can be more efficient when the end type happens to be u8 (bytes). These methods have u8 appended to their names to avoid confusion with Vecg methods. These specific algorithms are different to their generic equivalents in Vecg.

* Frequency count of bytes by their values (histogram, pdf, jointpdf)
* Entropy, jointentropy, independence.

## Trait VecVec

Relationships between n vectors (in d dimensions).
This general data domain is denoted here as (nd). It is in nd where the main original contribution of this library lies. True geometric median (gm) is found by fast and stable iteration, using improved Weiszfeld's algorithm `gmedian`. This algorithm solves Weiszfeld's convergence and stability problems in the neighbourhoods of existing set points. Its variant `pmedian` iterates point-by-point, which gives even better convergence.

* centroid, medoid, outliers, gm
* sums of distances, radius of a point (as its distance from gm)
* characterisation of a set of multidimensional points by the mean, standard deviation, median and MAD of its points' radii. These are useful recognition measures for the set.
* transformation to zero geometric median data,
* multivariate trend (regression) between two sets of nd points,
* covariance and comediance matrices.
* inner and outer hulls

## Trait VecVecg

Methods which take an additional generic vector argument, such as a vector of weights for computing weighted geometric medians (where each point has its own weight). Matrices multiplications.

## Appendix: Recent Releases

* **Version 1.2.24** - added `st_error` method to trait Vecg. It is a generalization of standard error to 'nd'. The central tendency is (usually) the geometric median and the spread is (usually) MADGM. Also tidied up `hulls`.

* **Version 1.2.23** - `convex_hull => hulls`. Now computes both inner and outer hulls. See above for definitions. Also, added `st_error` to auxiliary functions.

* **Version 1.2.22** - Improved Display of TriangMat - it now prints just the actual triangular form. Other minor cosmetic improvements.

* **Version 1.2.21** - Updated dependency `medians` to v 2.0.2 and made the necessary compatibility changes (see Quantify Functions above). Moved all remaining methods to do with 1d medians from here to crate `medians`. Removed auxiliary function i64tof64, as it was a trivial mapping of `as f64`. Made `dfdt` smoothed and median based.

* **Version 1.2.20** - Added `dfdt` to `Stats` trait (approximate weighted time series derivative at the last point). Added automatic conversions (with `?`) of any potential errors returned from crates `ran`, `medians` and `times`. Now demonstrated in `tests.rs`.

* **Version 1.2.19** - Presentation only: github actions now run automatically the full battery of `cargo test`. Detailed and informative tests output can be seen in the github actions log and overall success is indicated by the green badge at the head of this readme file.

* **Version 1.2.18** - Updated dependency `ran v1.0.4`. Added github action `cargo check`.

* **Version 1.2.17** - Rectangular matrices (as `Vec<Vec<T>>`): multiplications made more efficient. Added `mat` test to tests.rs.

* **Version 1.2.16** - TriangMat developed. Methods working with triangular matrices are now implemented for this struct.

* **Version 1.2.15** - Introducing `struct TriangMat`: better representation of triangular matrices.

* **Version 1.2.14** - Householder's UR matrix decomposition and orthogonalization.

* **Version 1.2.13** - Updated dependency to `indxvec v1.4.2`. Added `normalize` and `wvmean`. Added `radius` to `VecVec`. Added `wtukeyvec` a `dottukey`. Removed bulky test/tests.rs from the crate, get them from the github repository. Householder decomposition/orthogonalization is 'to do'.

* **Version 1.2.12** - Updated dependency `indxvec v1.3.4'.

* **Version 1.2.11** - Added `convex_hull` to trait VecVec. Added more error checking: VecVecg trait is now fully checked, be prepared to append `?` after most method calls.

* **Version 1.2.10** - Minor: corrected some examples, removed all unnecessary `.as_slice()` conversions.

* **Version 1.2.9** - More RError forwarding. Removed all deliberate panics.

* **Version 1.2.8** - Fixed a silly bug in `symmatrix` and made it return Result.

* **Version 1.2.7** - Added efficient `mahalanobis` distance and its test.

* **Version 1.2.6** - Added test `matrices` specifically for matrix operations. Added type alias `Rstats::RE` to shorten method headings returning `RErrors` that carry `&str` payloads (see subsection Errors above). 

* **Version 1.2.5** - Added some more matrix algebra. Added generic payload `T` to RError: `RError<T>` to allow it to carry more information. 

* **Version 1.2.4** - Added Cholesky–Banachiewicz algorithm `cholesky` to trait `Statsg` for efficient matrix decomposition.

* **Version 1.2.3** - Fixed `hwmeanstd`. Some more tidying up using RError. `Autocorr` and `lintrans` now also check their data and return `Result`. 

* **Version 1.2.2** - Introduced custom error RError, potentially returned by some methods of trait `Statsg`. Removed the dependency on crate `anyhow`.

* **Version 1.2.1** - Code pruning - removed `wsortedcos` of questionable utility from trait `VecVecg`.
