# Rstats

[<img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github">](https://github.com/liborty/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/v/rstats?logo=rust">](https://crates.io/crates/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/d/rstats?logo=rust">](https://crates.io/crates/rstats)
[<img alt="docs.rs" src="https://img.shields.io/docsrs/rstats?logo=rust">](https://docs.rs/rstats)
[![Actions Status](https://github.com/liborty/rstats/workflows/compilation/badge.svg)](https://github.com/liborty/rstats/actions)

Statistics, Linear Algebra, Cholesky Matrix Decomposition, Mahalanobis Distance, Householder Matrix Decomposition, Information Measures, Machine Learning, Data Analysis, Geometric Median, Convex Hull and more ...

## Usage

Insert `rstats = "^1"` in the `Cargo.toml` file, under `[dependencies]`.

Use in your source files any of the following structs, when needed:

```rust  
use rstats::{RE,TriangMat,Mstats,MinMax,Med};
```

and any of the following rstats defined traits:

```rust 
use rstats::{Stats,Vecg,Vecu8,MutVecg,VecVec,VecVecg};
```

The latest (nightly) version is always available in the github repository [rstats](https://github.com/liborty/rstats). Sometimes it may be a little ahead of the crates.io release versions.

It is highly recommended to read and run [tests.rs](https://github.com/liborty/rstats/blob/master/tests/tests.rs) from the github repository as examples of usage.
To run all the tests, use single thread in order to produce the results in the right order:

```bash  
cargo test --release -- --test-threads=1 --nocapture --color always
```

## Introduction

`Rstats` is primarily about characterising sets of n points, each represented as a `Vec`, in space of `d` dimensions (where d is `vec.len()`). It has applications mostly in Machine Learning and Data Analysis. 

Several branches of mathematics: statistics, linear algebra, information theory and set theory are combined in this one consistent crate, based on the fact that they all operate on these same objects. The only difference being that an ordering of their components is sometimes assumed (linear algebra, set theory) and sometimes it is not (statistics, information theory, set theory).

`RStats` begins with basic statistical measures and vector algebra, which provide self-contained tools for the machine learning (ML) multidimensional algorithms but can also be used in their own right. 

`Non analytical statistics` is preferred, whereby the `random variables` are replaced by vectors of real data. Probabilities densities and other parameters are always obtained from the data, not from some assumed distributions.

Our treatment of multidimensional sets of points is constructed from the first principles. Some original concepts, not found elsewhere, are introduced and implemented here:

* `median correlation`- in one dimension, our `mediancorr` method is to replace `Pearson's correlation`. We define `median correlation` as cosine of an angle between two zero median samples (instead of Pearson's zero mean samples). This conceptual clarity is one of the benefits of our thinking of a sample as a vector in d dimensional space. 

* `gmedian and pmedian` - fast multidimensional `geometric median (gm)` algorithms.

* `madgm` - generalisation to `nd` of a robust data spread estimator known as `MAD` (median of absolute deviations from median). 

* `contribution` - of a point to an nd set. Defined as a magnitude of gm adjustment, when the point is added to the set. It is related to the point's radius (distance from the gm) but not the same, as it depends on the radii of all the other points as well.

* `comediance` - instead of covariance (matrix). It is obtained by supplying `covar` with the geometric median instead of the usual centroid. Thus `zero median vectors` are replacing `zero mean vectors` in covariance calculations.

*Zero median vectors are generally preferable to the commonly used zero mean vectors.*

In n dimensions, many authors 'cheat' by using `quasi medians` (1-d medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. In a highly dimensional space, they are also much slower to compute than our gm.

*Specifically, all such 1d measures are sensitive to the choice of axis and thus are affected by their rotation.*

In contrast, analyses based on the true geometric median (gm) are axis (rotation) independent. Also, they are more stable, as medians have a 50% breakdown point (the maximum possible). They are computed here by methods `gmedian` and its weighted version `wgmedian`, in traits `vecvec` and `vecvecg` respectively.

## Terminology

#### Including some new definitions for sets of nd points, i.e. n points in d dimensional space

* `Median correlation` between
 two samples. We define it analogously to Pearson, as cosine of an angle between two 'normalised' vectors. Pearson 'normalises' by subtracting the mean from all components, we subtract the median.

* `Centroid/Centre/Mean` of an nd set. It is the (generally non member) point that minimises its sum of *squares* of distances to all member points. Thus it is susceptible to outliers. Specifically, it is the d-dimensional arithmetic mean. It is sometimes called 'the centre of mass'. Centroid can also sometimes mean the member of the set which is the nearest to the Centre. Here we follow the common usage: Centroid = Centre = Arithmetic Mean.

* `Quasi/Marginal Median` is the point minimising sums of distances separately in each dimension (its coordinates are  medians along each axis). It is a mistaken concept which we do not recommend using.

* `Tukey Median` is the point maximising `Tukey's Depth`, which is the minimum number of (outlying) points found in a hemisphere in any direction. Potentially useful concept but only partially implemented here by `tukeyvec`, as its advantages over the geometric median are not clear.

* `Median or the true geometric median (gm)`, is the point (generally non member), which minimises the sum of distances to all member points. This is the one we want. It is much less susceptible to outliers than the centroid. In addition, unlike quasi median, `gm` is rotation independent.

* `Medoid` is the member of the set with the least sum of distances to all other members. Equivalently, the member which is the nearest to the `gm`.

* `Outlier` is the member of the set with the greatest sum of distances to all other members. Equivalently, it is the point furthest from the `gm`.

* `Convex Hull` is the subset consisting of selected points p, such that no other member points lie outside the plane through p and normal to its radius vector. The remaining points are the `internal` points.

* `Zero median vectors` are obtained by subtracting the `gm` (placing the origin of the coordinate system at the `gm`). This is a proposed  alternative to the commonly used `zero mean vectors`, obtained by subtracting the centroid.

* `MADGM` (median of distances from gm). This is our generalisation of `MAD` (median of absolute differences) measure from one dimension to any number of dimensions (d>1). It is a robust measure of nd data spread.

* `Comediance` is similar to `covariance`, except that zero median vectors are used to compute it,  instead of zero mean vectors.

* `Mahalanobis Distance` is a weighted distace, where the weights are derived from the axis of variation of the nd data points cloud. Thus distances in the directions in which there are few points are penalised (increased) and vice versa. Efficient Cholesky singular (eigen) value decomposition is used. Cholesky method decomposes the covariance/comediance positive definite matrix S into a product of two triangular matrices: S = LL'. See more details in the code comments.

* `Contribution` One of the key questions of Machine Learning (ML) is how to quantify the contribution that each example point (typically a member of some large nd set) makes to the recognition concept, or class, represented by that set. In answer to this, we define the `contribution` of a point as the magnitude of adjustment to gm caused by adding that point. Generally, outlying points make greater contributions to the gm but not as much as they would to the centroid. The contribution depends not only on the radius of the example point in question but also on the radii of all other (existing) examples.

* `Tukey Vector` Proportions of points in each hemisphere from gm. This is a useful 'signature' of a data cloud. For a new point (that typically needs to be classified) we can then quickly determine whether it lies in a well populated direction. This could be done by projecting all the existing points on it but that is much slower, as there are many. Also, in keeping with the stability properties of medians, we are only using counts of points, not their distances.


## Implementation

The main constituent parts of Rstats are its traits. The selection of traits (to import) is primarily determined by the types of objects to be handled. These are mostly vectors of arbitrary length (dimensionality). The main traits are implementing methods applicable to:

* `Stats`: a single vector (of numbers),
* `Vecg`: methods (of vector algebra) operating on two vectors, e.g. scalar product
* `Vecu8`: some special methods for end-type u8
* `MutVecg`: some of the above methods, mutating self
* `VecVec`: methods operating on n vectors, 
* `VecVecg`: methods for n vectors, plus another generic argument, e.g. vector of weights.

In other words, the traits and their methods operate on arguments of their required categories. In classical statistical terminology, the main categories correspond to the number of 'random variables'.

`Vec<Vec<T>>` type is used for full rectangular matrices, whereas `TriangMat` struct is used for symmetric and triangular matrices (to save memory).

The vectors' end types (for the actual data) are mostly generic: usually some numeric type. There are also some traits specialised for input end type `u8` and some that take mutable self. End type `f64` is most commonly used for the results.

## Errors

RStats crate produces custom errors `RError`:

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
Each of its enum variants also carries a generic payload `T`. Most commonly this will be simply a `&'static str` message giving more helpful explanation, e.g.:

```rust 
return Err(RError::ArithError("cholesky needs a positive definite matrix"));
```

 There is a type alias shortening return declarations to, e.g.: `Result<Vec<f64>,RE>`, where

 ```rust
pub type RE = RError<&'static str>;
```

More error checking will be added in later versions, where it makes sense. 

## Documentation

For more detailed comments, plus some examples, see the source. You may have to unclick the 'implementations on foreign types' somewhere near the bottom of the page in the rust docs to get to it.  (Since these traits are implemented over the pre-existing Rust Vec type).

## Structs

### `struct MStats` 
holds the central tendency, e.g. some kind of mean or median, and dispersion, e.g. standard deviation or MAD.

### `struct TriangMat` 
holds lower/upper triangular symmetric/non-symmetric matrix in compact form that avoids zeros and duplications. Beyond the usual conversion to full matrix form, a number of (the best) Linear Algebra methods are implemented directly on TriangMag, in module `triangmag.rs`, such as:

* **Cholesky-Banachiewicz** matrix decomposition: M = LL' (where ' denotes a transpose), used by:
* **Mahalanobis Distance**
* **Householder UR** (M = QR) matrix decomposition

Also, there are some methods implemented for `VecVecg` that produce `TriangMat`, specifically the covariance/comedience calculations: `covar`,`wcovar`,`comed` and `wcomed`. Their results will be typically used by `mahalanobis`.


##  Auxiliary Functions

* `i64tof64`: converts an i64 vector to f64, 
* `sumn`: sum of a sequence 1..n, also the size of a lower/upper triangular matrix below/above the diagonal (n*(n+1)/2.),
* `unit_matrix` full unit matrix

## Trait Stats

One dimensional statistical measures implemented for all numeric end types.

Its methods operate on one slice of generic data and take no arguments.
For example, `s.amean()` returns the arithmetic mean of the data in slice `s`.
These methods are checked and will report RError(s), such as an empty input. This means you have to apply `?` to their results to pass the errors up, or explicitly match them to take recovery actions, depending on the variant.

Included in this trait are:

* means (arithmetic, geometric and harmonic),
* standard deviations,
* linearly weighted means (useful for time dependent data analysis),
* probability density function (pdf)
* autocorrelation, entropy
* linear transformation to [0,1],
* other measures and vector algebra operators

Note that fast implementation of 1d medians is as of version 1.1.0 in crate `medians`.  


## Trait Vecg

Generic vector algebra operations between two slices `&[T]`, `&[U]` of any length (dimensionality). It may be necessary to invoke some using the 'turbofish' `::<type>` syntax to indicate the type U of the supplied argument, e.g.:  
`datavec.methodname::<f64>(arg)`  
This is because Rust is currently incapable of inferring the type ('the inference bug').

* Vector additions, subtractions and products (scalar, kronecker, outer),
* Other relationships and measures of difference,
* Pearson's, Spearman's and Kendall's correlations,
* `Median correlation`, which we define analogously to Pearson's, as cosine of an angle between two zero median vectors (instead of his zero mean vectors).
* Joint pdf, joint entropy, statistical independence (based on mutual information).
* `Contribution` measure of a point w.r.t gm

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
* convex hull points

## Trait VecVecg

Methods which take an additional generic vector argument, such as a vector of weights for computing weighted geometric medians (where each point has its own weight). Matrices multiplications.

## Appendix: Recent Releases

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

* **Version 1.2.6** - Added test `matrices` specifically for matrix operations. Added type alias `RStats::RE` to shorten method headings returning `RErrors` that carry `&str` payloads (see subsection Errors above). 

* **Version 1.2.5** - Added some more matrix algebra. Added generic payload `T` to RError: `RError<T>` to allow it to carry more information. 

* **Version 1.2.4** - Added Cholesky–Banachiewicz algorithm `cholesky` to trait `Statsg` for efficient matrix decomposition.

* **Version 1.2.3** - Fixed `hwmeanstd`. Some more tidying up using RError. `Autocorr` and `lintrans` now also check their data and return `Result`. 

* **Version 1.2.2** - Introduced custom error RError, potentially returned by some methods of trait `Statsg`. Removed the dependency on crate `anyhow`.

* **Version 1.2.1** - Code pruning - removed `wsortedcos` of questionable utility from trait `VecVecg`.
