# Rstats [![crates.io](https://img.shields.io/crates/v/Rstats?logo=rust)](https://crates.io/crates/rstats) [![crates.io](https://img.shields.io/crates/d/rstats?logo=rust)](https://crates.io/crates/rstats) [![GitHub last commit](https://img.shields.io/github/last-commit/liborty/Rstats/HEAD?logo=github)](https://github.com/liborty/Rstats) [![Actions Status](https://github.com/liborty/rstats/actions/workflows/tests.yml/badge.svg)](https://github.com/liborty/rstats/actions)

## Author: Libor Spacek

## Usage

This crate is written in 100% safe Rust.

Use in your source files any of the following structs, as and when needed:

```rust  
use Rstats::{RE,RError,Params,TriangMat,MinMax};
```

and any of the following traits:

```rust
use Rstats::{Stats,Vecg,Vecu8,MutVecg,VecVec,VecVecg};
```

and any of the following auxiliary functions:

```rust
use Rstats::{
    fromop,sumn,tm_stat,unit_matrix,nodata_error,data_error,
    arith_error,other_error };
```

or just simply use everything:

```rust
use Rstats::*;
```

The latest (nightly) version is always available in the github repository [Rstats](https://github.com/liborty/Rstats). Sometimes it may be (only in some details) a little ahead of the `crates.io` release versions.

It is highly recommended to read and run [tests.rs](https://github.com/liborty/Rstats/blob/master/tests/tests.rs) for examples of usage. To run all the tests, use a single thread in order not to print the results in confusing mixed-up order:

```bash  
cargo test --release -- --test-threads=1 --nocapture
```

However, `geometric_medians`, which compares multithreading performance, should be run separately in multiple threads, as follows:

```bash
cargo test -r geometric_medians -- --nocapture
```

Alternatively, just to get a quick idea of the methods provided and their usage, read the output produced by an [automated test run](https://github.com/liborty/rstats/actions). There are test logs generated for each new push to the github repository. Click the latest (top) one, then `Rstats` and then `Run cargo test` ... The badge at the top of this document lights up green when all the tests have passed and clicking it gets you to these logs as well.

Any compilation errors arising out of `rstats` crate indicate most likely that some of the dependencies are out of date. Issuing `cargo update` command will usually fix this.

## Introduction

`Rstats` has a small footprint. Only the best methods are implemented, primarily with Data Analysis and Machine Learning in mind. They include multidimensional (`nd` or 'hyperspace') analysis, i.e. characterising clouds of n points in space of d dimensions.

Several branches of mathematics: statistics, information theory, set theory and linear algebra are combined in this one consistent crate, based on the abstraction that they all operate on the same data objects (here Rust Vecs). The only difference being that an ordering of their components is sometimes assumed (in linear algebra, set theory) and sometimes it is not (in statistics, information theory, set theory).

`Rstats` begins with basic statistical measures, information measures, vector algebra and linear algebra. These provide self-contained tools for the multidimensional algorithms but they are also useful in their own right.

`Non analytical (non parametric) statistics` is preferred, whereby the 'random variables' are replaced by vectors of real data. Probabilities densities and other parameters are in preference obtained from the real data (pivotal quantity), not from some assumed distributions.

`Linear algebra` uses generic data structure `Vec<Vec<T>>` capable of representing irregular matrices.

`Struct TriangMat` is defined and used for symmetric, anti-symmetric, and triangular matrices, and their transposed versions, saving memory.

Our treatment of multidimensional sets of points is constructed from the first principles. Some original concepts, not found elsewhere, are defined and implemented here (see the next section).

*Zero median vectors are generally preferred to commonly used zero mean vectors.*

In n dimensions, many authors 'cheat' by using `quasi medians` (one dimensional (`1d`) medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. Also, they are actually slower to compute than our **gm** ( `geometric median`), as soon as the number of dimensions exceeds trivial numbers.

*Specifically, all such 1d measures are sensitive to the choice of axis and thus are affected by their rotation.*

In contrast, our methods based on **gm** are axis (rotation) independent. Also, they are more stable, as medians have a maximum possible breakdown point.

We compute geometric medians by our method `gmedian` and its parallel version `par_gmedian` in trait `VecVec` and their weighted versions `wgmedian` and `par_wgmedian` in trait `VecVecg`. It is mostly these efficient algorithms that make our new concepts described below practical.

### Additional Documentation

For more detailed comments, plus some examples, see [rstats in docs.rs](https://docs.rs/rstats/latest/rstats). You may have to go directly to the modules source. These traits are implemented for existing 'out of this crate' rust `Vec` type and unfortunately rust docs do not display 'implementations on foreign types' very well.

## New Concepts and their Definitions

* `zero median points` (or vectors) are obtained by moving the origin of the coordinate system to the median (in `1d`), or to the **gm** (in `nd`). They are our alternative to the commonly used `zero mean points`, obtained by moving the origin to the arithmetic mean (in 1d) or to the arithmetic centroid (in `nd`).

* `median correlation` between two 1d sets of the same length.  
We define this correlation similarly to Pearson, as cosine of an angle between two normalised samples, interpreted as coordinate vectors. Pearson first normalises each set by subtracting its mean from all components. Whereas we subtract the median (cf. zero median points above). This conceptual clarity is one of the benefits of interpreting a data sample of length `d` as a single vector in `d` dimensional space.

* `gmedian, par_gmedian, wgmedian and par_wgmedian`  
our fast multidimensional `geometric median` (**gm**) algorithms.

* `madgm` (median of distances from `gm`)  
is our generalisation of `mad` (**m**edian of **a**bsolute **d**ifferences from median), to n dimensions. `1d` median is replaced in `nd` by **gm**. Where `mad` was a robust measure of 1d data spread, `madgm` becomes a robust measure of `nd` data spread. We define it as: `median(`**|pi-gm|**,`i=1..n)`, where **p1`..`pn** are a sample of n data points, each of which is now a vector.

* `tm_stat`  
 We define our generalized `tm_stat` of a single scalar observation x as: `(x-centre)/spread`, with the recommendation to replace mean by median and `std` by `mad`, whenever possible. Compare with common `t-stat`, defined as `(x-mean)/std`, where `std` is the standard deviation.  
 These are similar to the well known `standard z-score`, except that the central tendency and spread are obtained from the sample (pivotal quantity), rather than from any old assumed population distribution.

* `tm_statistic`  
we now generalize `tm_stat` from scalar domain to vector domain of any number of dimensions, defining `tm_statistic` as |**p-gm**|`/madgm`, where **p** is a single observation point in `nd` space. For sample central tendency now serves the `geometric median` **gm** vector and the spread is the `madgm` scalar (see above). The error distance of observation **p** from the median: **|p-gm|**, is also a scalar. Thus the co-domain of `tm_statistic` is a simple positive scalar, regardless of the dimensionality of the vector space in question.

* `contribution`  
one of the key questions of Machine Learning is how to quantify the contribution that each example (typically represented as a member of some large `nd` set) makes to the recognition concept, or outcome class, represented by that set. In answer to this, we define the `contribution` of a point **p** as the magnitude of displacement of `gm`, caused by adding **p** to the set. Generally, outlying points make greater contributions to the `gm` but not as much as to the `centroid`. The contribution depends not only on the radius of **p** but also on the radii of all other existing set points and on their number.

* `comediance`  
is similar to `covariance`. It is a triangular symmetric matrix, obtained by supplying method `covar` with the geometric median instead of the usual centroid. Thus `zero mean vectors` are replaced by `zero median vectors` in the covariance calculations. The results are similar but more stable with respect to the outliers.

* `outer_hull` is a subset of all zero median points **p**, such that no other points lie outside the normal plane through **p**. The points that do not satisfy this condition are called the `internal` points.

* `inner_hull` is a subset of all zero median points **p**, that do not lie outside the normal plane of any other point. Note that in a highly dimensional space up to all points may belong to both the inner and the outer hulls, as, for example, when they all lie on the same hypersphere.

* `depth` is a measure of likelihood of a zero median point **p** belonging to a data cloud. More specifically, it is the projection onto unit **p** of a sum of unit vectors that lie outside the normal through **p**. For example, all outer hull points have by their definition `depth = 0`, whereas the inner hull points have high values of depth. This is intended as an improvement on Mahalanobis distance which has a similar goal but says nothing about how well enclosed **p** is. Whereas `tm_statistic` only informs about the probability pertaining to the whole cloud, not to its local shape near **p**.

* `sigvec (signature vector)`  
Proportional projections of a cloud of zero median vectors on all hemispheric axis. When a new zero median point **p** needs to be classified, we can quickly estimate how well populated is its direction from **gm**. Similar could be done by projecting all the points directly onto **p** but this is usually impractically slow, as there are typically very many such points. However, `signature_vector` only needs to be precomputed once and is then the only vector to be projected onto **p**.

## Previously Known Concepts and Terminology

* `centroid/centre/mean` of an `nd` set.  
Is the point, generally non member, that minimises its sum of *squares* of distances to all member points. The squaring makes it susceptible to outliers. Specifically, it is the d-dimensional arithmetic mean. It is sometimes called 'the centre of mass'. Centroid can also sometimes mean the member of the set which is the nearest to the Centre. Here we follow the common usage: Centroid = Centre = Arithmetic Mean.

* `quasi/marginal median`  
is the point minimising sums of distances separately in each dimension (its coordinates are medians along each axis). It is a mistaken concept which we do not recommend using.

* `Tukey median`  
is the point maximising `Tukey's Depth`, which is the minimum number of (outlying) points found in a hemisphere in any direction. Potentially useful concept but its advantages over the geometric median are not clear.

* `true geometric median` (**gm**)  
is the point (generally non member), which minimises the sum of distances to all member points. This is the one we want. It is less susceptible to outliers than the centroid. In addition, unlike quasi median, **gm** is rotation independent.

* `medoid`  
is the member of the set with the least sum of distances to all other members. Equivalently, the member which is the nearest to the **gm** (has the minimum radius).

* `outlier`  
is the member of the set with the greatest sum of distances to all other members. Equivalently, it is the point furthest from the **gm** (has the maximum radius).

* `Mahalanobis distance`  
is a scaled distance, whereby the scaling is derived from the axes of covariances / `comediances` of the data points cloud. Distances in the directions in which there are few points are increased and distances in the directions of significant covariances / `comediances` are decreased. Requires matrix decomposition.

* `Cholesky-Banachiewicz matrix decomposition`  
decomposes any positive definite matrix S (often covariance or comediance matrix) into a product of two triangular matrices: `S = LL'`. The eigenvalues and the determinant are easily obtained from the diagonal of L. We implemented it on `TriangMat` for maximum efficiency. It is used mainly by `eigenvalues`, `eigenvectors`, `mahalanobis` and `pca_reduction`.

* `Householder's decomposition`  
in cases where the precondition (positive definite matrix S) for the Cholesky-Banachiewicz decomposition is not satisfied, Householder's (UR) decomposition is often used as the next best method. It is implemented here on our efficient `struct TriangMat`.

* `wedge product, geometric product`  
products of the Grassman and Clifford algebras, respectively. Wedge product is used here to generalize the cross product of two vectors into any number of dimensions, determining the correct sign (sidedness of their common plane).

## Implementation Notes

The main constituent parts of Rstats are its traits. The different traits are determined by the types of objects to be handled. The objects are mostly vectors of arbitrary length/dimensionality (`d`). The main traits are implementing methods applicable to:

* `Stats`: a single vector (of numbers),
* `Vecg`: methods operating on two vectors, e.g. scalar product,
* `Vecu8`: some methods specialized for end-type `u8`,
* `MutVecg`: some of the above methods, mutating self,
* `VecVec`: methods operating on n vectors (rows of numbers),
* `VecVecg`: methods for n vectors, plus another generic argument, typically a vector of n weights, expressing the relative significance of the vectors.

The traits and their methods operate on arguments of their required categories. In classical statistical parlance, the main categories correspond to the number of 'random variables'.

**`Vec<Vec<T>>`** type is used for rectangular matrices (could also have irregular rows).
  
**`struct TriangMat`** is used for symmetric / antisymmetric / transposed / triangular matrices and wedge and geometric products. All instances of `TriangMat` store only `n*(n+1)/2` items in a single flat vector, instead of `n*n`, thus almost halving the memory requirements. Their transposed versions only set up a flag `kind >=3` that is interpreted by software, instead of unnecessarily rewriting the whole matrix. Thus saving processing of all transposes (a common operation). All this is put to a good use in our implementation of the matrix decomposition methods.

The vectors' end types (of the actual data) are mostly generic: usually some numeric type. `Copy` trait bounds on these generic input types have been relaxed to `Clone`, to allow cloning user's own end data types in any way desired. There is no difference for primitive  types.

The computed results end types are usually `f64`.

## Errors

`Rstats` crate produces custom error `RError`:

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

Each of its enum variants also carries a generic payload `T`. Most commonly this will be a `String` message, giving more helpful explanation, e.g.:

```rust
if dif <= 0_f64 {
    return Err(RError::ArithError("cholesky needs a positive definite matrix".to_owned())));
};
```

`format!(...)` can be used to insert (debugging) run-time values to the payload String. These errors are returned and can then be automatically converted (with `?`) to users' own errors. Some such error conversions are implemented at the bottom of `errors.rs` file and used in `tests.rs`.

 There is a type alias shortening return declarations to, e.g.: `Result<Vec<f64>,RE>`, where

 ```rust
pub type RE = RError<String>;
```

Convenience functions `nodata_error, data_error, arith_error, other_error` are used to construct and return these errors. Their message argument can be either literal `&str`, or `String` (e.g. constructed by `format!`). They return `ReError<String>` already wrapped up as an `Err` variant of `Result`. cf.:

```rust
if dif <= 0_f64 {
    return arith_error("cholesky needs a positive definite matrix");
};
```

## Structs

### `struct Params`

holds the central tendency of `1d` data, e.g. some kind of mean or median, and its spread measure, e.g. standard deviation or 'mad'.

### `struct TriangMat`

holds triangular matrices of all kinds, as described in Implementation section above. Beyond the usual conversion to full matrix form, a number of (the best) Linear Algebra methods are implemented directly on `TriangMat`, in module `triangmat.rs`, such as:

* **Cholesky-Banachiewicz** matrix decomposition: `S = LL'` (where ' denotes the transpose). This decomposition is used by `mahalanobis`, `eigenvalues`, etc.
* **Eigenvectors, eigenvalues, determinant**
* **Mahalanobis Distance** for ML recognition tasks
* **Householder UR** (M = QR) more general matrix decomposition

Some methods, specifically the covariance/comedience calculations in `VecVec` and `VecVecg` return `TriangMat` matrices. These are positive definite, which makes the most efficient Cholesky-Banachiewicz decomposition applicable to them.

## Quantify Functions (Dependency Injection)

Most methods in `medians` and some in `indxvec` crates, e.g. `find_any` and `find_all`, require explicit closure passed to them, usually to tell them how to quantify input data of any type T into f64. Variety of different quantifying methods can then be dynamically employed.

For example, in text analysis (`&str` end type), it can be the word length, or the numerical value of its first few letters, or the numerical value of its consonants, etc. Then we can sort them or find their means / medians / spreads under all these different measures. We do not necessarily want to explicitly store all such different quantifications values, as input data can be voluminous. It is preferable to be able to compute any of them on demand, using these closure arguments.

When data is already of the required end-type, use the 'dummy' closure:

```rust
|&f| f
```

When T is a primitive type, such as i64, u64, usize, that can be converted to f64, possibly with some loss of accuracy, use:

```rust
|&f| f as f64
```

### `fromop`

When T is convertible by an existing custom `From` implementation (and `f64:From<T>, T:Clone` have been duly added everywhere as trait bounds), then simply pass in `fromop`, defined as:

```rust
/// Convenience From quantification invocation
pub fn fromop<T: Clone + Into<f64>>(f: &T) -> f64 {
    f.clone().into()
}|
```

The remaining general cases previously required new manual implementations to be written for the (global) `From` trait for each new type and for each different quantification method, and adding the trait bounds everywhere. Even then, the different implementations of `From` would conflict with each other. Now we can simply implement all the custom quantifications within the closures. This generality is obtained at the price of a small inconvenience: having to supply one of the above closures argument for the primitive types as well.

## Auxiliary Functions

* `fromop`: see above.

* `sumn`: the sum of the sequence `1..n = n*(n+1)/2`. It is also the size of a lower/upper triangular matrix.

* `tm_stat`: (x-centre)/dispersion. Generalised t-statistic in one dimension.

* `unit_matrix`: - generates full square unit matrix.

* `nodata_error, data_error, arith_error, other_error` - construct custom RE errors (see section Errors above).

## Trait Stats

One dimensional statistical measures implemented for all numeric end types.

Its methods operate on one slice of generic data and take no arguments.
For example, `s.amean()?` returns the arithmetic mean of the data in slice `s`.
These methods are checked and will report RError(s), such as an empty input. This means you have to apply `?` to their results to pass the errors up, or explicitly match them to take recovery actions, depending on the error variant.

Included in this trait are:

* 1d medians (classic, geometric and harmonic) and their spreads
* 1d means (arithmetic, geometric and harmonic) and their spreads
* linearly weighted means (useful for time analysis),
* probability density function (pdf)
* autocorrelation, entropy
* linear transformation to [0,1],
* other measures and basic vector algebra operators

Note that fast implementations of 1d 'classic' medians are, as of version 1.1.0, provided in a separate crate `medians`.

## Trait Vecg

Generic vector algebra operations between two slices `&[T]`, `&[U]` of any (common) length  (dimensions). Note that it may be necessary to invoke some using the 'turbofish' `::<type>` syntax to indicate the type U of the supplied argument, e.g.:

```rust
datavec.somemethod::<f64>(arg)
```

Methods implemented by this trait:

* Vector additions, subtractions and products (scalar, kronecker, outer),
* Other relationships and measures of difference,
* Pearson's, Spearman's and Kendall's correlations,
* Joint pdf, joint entropy, statistical independence (based on mutual information).
* `Contribution` measure of a point's impact on the geometric median

Note that our `median correlation` is implemented in a separate crate `medians`.

Some simpler methods of this trait may be unchecked (for speed), so some caution with data is advisable.

## Trait MutVecg

A select few of the `Stats` and `Vecg` methods (e.g. mutable vector addition, subtraction and multiplication) are reimplemented under this trait, so that they can mutate `self` in-place. This is more efficient and convenient in some circumstances, such as in vector iterative methods.

However, these methods do not fit in with the functional programming style, as they do not explicitly return anything (their calls are statements with side effects, rather than expressions).

## Trait Vecu8

Some vector algebra as above that can be more efficient when the end type happens to be u8 (bytes). These methods have u8 appended to their names to avoid confusion with Vecg methods. These specific algorithms are different to their generic equivalents in Vecg.

* Frequency count of bytes by their values (histogram, pdf, jointpdf)
* Entropy, jointentropy, independence.

## Trait VecVec

Relationships between n vectors in d dimensions.
This (hyper-dimensional) data domain is denoted here as (`nd`). It is in `nd` where the main original contribution of this library lies. True geometric median (gm) is found by fast and stable iteration, using improved Weiszfeld's algorithm `gmedian`. This algorithm solves Weiszfeld's convergence and stability problems in the neighbourhoods of existing set points. Its variant, `par_gmedian`, employs multithreading for faster execution and gives otherwise  the same result.

* centroid, medoid, outliers, gm
* sums of distances, radius of a point (as its distance from gm)
* characterisation of a set of multidimensional points by the mean, standard deviation, median and MAD of its points' radii. These are useful recognition measures for the set.
* transformation to zero geometric median data,
* multivariate trend (regression) between two sets of `nd` points,
* covariance and comediance matrices.
* inner and outer hulls

## Trait VecVecg

Methods which take an additional generic vector argument, such as a vector of weights for computing weighted geometric medians (where each point has its own significance weight). Matrices multiplications.

## Appendix: Recent Releases

* **Version 2.1.8** - Improved `TriangMat::diagonal()`, restored `TriangMat::determinant()`, tidied up `triangmat` test.

* **Version 2.1.7** - Removed suspect eigen values/vectors computations. Improved 'householder' test.

* **Version 2.1.5** - Added `projection` to trait `VecVecg` to project all self vectors to a new basis. This can be used e.g. for Principal Components Analysis data reduction, using some of the eigenvectors as the new basis.

* **Version 2.1.4** - Tidied up some error processing.

* **Version 2.1.3** - Added `normalize` (normalize columns of a matrix and transpose them to rows).

* **Version 2.1.2** - Added function `project` to project a `TriangMat` to a lower dimensional space of selected dimensions. Removed `rows` which was a duplicate of `dim`.

* **Version 2.1.0** - Changed the type of `mid` argument to covariance methods from U -> f64, making the normal expectation for the type of precise geometric medians explicit. Accordingly, moved `covar` and `serial_covar` from trait `VecVecg` to `VecVec`. This might potentially require changing some `use` declarations in your code.

* **Version 2.0.12** - added `depth_ratio`

* **Version 2.0.11** - removed not so useful `variances`. Tidied up error processing in `vecvecg.rs`. Added to it `serial_covar` and `serial_wcovar` for when heavy loading of all the cores may not be wanted.

* **Version 2.0.9** - Pruned some rarely used methods, simplified `gmparts` and `gmerror`, updated dependencies.

* **Version 2.0.8**' - Changed initial guess in iterative weighted gm methods to weighted mean. This, being more accurate than plain mean, leads to fewer iterations. Updated some dependencies.

* **Version 2.0.7** - Updated to `ran 2.0`.

* **Version 2.0.6** - Added convenience method `medmad` to Stats trait. It packs median and mad into `struct Params`, similarly to `ameanstd` and others. Consequently simplified the printouts in some tests.

* **Version 2.0.5** - Corrected `wsigvec` to also return normalized result. Updated dependency `Medians` to faster version 3.0.1.

* **Version 2.0.4** - Made a corresponding change: `winsideness` -> `wdepth`.

* **Version 2.0.3** - Improved `insideness` to be projection of a sum of unit vectors instead of just a simple count. Renamed it to `depth` to avoid confusion. Also some fixes to `hulls`.

* **Version 2.0.2** - Significantly speeded up `insideness` and added weighted version `winsideness` to `VecVecg` trait.

* **Version 2.0.1** - Added `TriangMat::dim()` and tidied up some comments.

* **Version 2.0.0** - Renamed `MStats` -> `Params` and its variant `dispersion` -> `spread`. This may cause some backwards incompatibilities, hence the new major version. Added 'centre' as an argument to `dfdt`,`dvdt`,`wdvdt`, so that it does not have to be recomputed.
