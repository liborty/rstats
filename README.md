# Rstats

[<img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github">](https://github.com/liborty/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/v/rstats?logo=rust">](https://crates.io/crates/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/d/rstats?logo=rust">](https://crates.io/crates/rstats)
[<img alt="docs.rs" src="https://img.shields.io/docsrs/rstats?logo=rust">](https://docs.rs/rstats)

## Usage

Insert `rstats = "^1"` in the `Cargo.toml` file, under `[dependencies]`.

Use in your source files any of the following structs, as needed:  
`use rstats::{MinMax,Med,Mstats};`  
and any of the following traits:  
`use rstats::{Stats,MutStats,Vecu8,Vecf64,Vecg,MutVecg,VecVec,VecVecg};`  
and any of the following helper functions:  
`use rstats::{i64tof64,wsum};`

The latest (nightly) version of this readme file and everything, is always available in the github repository [rstats](https://github.com/liborty/rstats). Sometimes it may be a little ahead of the crates.io release versions.

It is highly recommended to read and run `tests/tests.rs`, which shows examples of usage.

To run all the tests, use single thread in order to produce the results in the right order:  
`cargo test --release -- --test-threads=1 --nocapture --color always`

## Introduction

`Rstats` is primarily about characterising multidimensional sets of points, with applications to Machine Learning and Big Data Analysis. It uses `non analytical statistics`, where the 'random variables' are replaced by vectors of real data. Probabilities densities and other parameters are always obtained from the data, not from some assumed distributions.

This crate begins with basic statistical measures and vector algebra, which provide self-contained tools for the multidimensional algorithms but can also be used in their own right.

Our treatment of multidimensional sets of points (vectors) is constructed from the first principles. Some original concepts, not found elsewhere, are introduced and implemented here:

* `gmedian` - fast multidimensional (geometric) median algorithm.

* `madgm` - generalisation of robust data spread estimator known as 'MAD' in 1d (median of absolute deviations from median), to multiple dimensions (nd).    

* `comediance` - instead of covariance (matrix). It is obtained by supplying `covar` with the geometric median instead of the usual centroid. Thus *zero median vectors* are replacing *zero mean vectors* in covariance calculations.

* `median correlation`- in one dimension, our `mediancorr` method is to replace *Pearson's correlation*. We define *median correlation*  as cosine of an angle between two zero median vectors (instead of Pearson's zero mean vectors).

*Zero median vectors are generally preferable to the commonly used zero mean vectors.*

In n dimensions, many authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. In a highly dimensional space, they are not even any faster to compute.

*Specifically, all such 1-d measures are sensitive to the choice of axis and thus are affected by rotation.*

In contrast, analyses based on the true geometric median (gm) are axis (rotation) independent. Also, they are more stable, as medians have a 50% breakdown point (the maximum possible). They are computed here by methods `gmedian` and its weighted version `wgmedian`, in traits `vecvec` and `vecvecg` respectively.

## Implementation

The main constituent parts of Rstats are its traits. The selection of traits (to import) is primarily determined by the types of objects to be handled. These are mostly vectors of arbitrary length (dimensionality). The main traits are implementing methods applicable to:

* `Stats`: a single vector (of numbers),
* `Vecg`: methods (of vector algebra) operating on two vectors, e.g. scalar product
* `VecVec`: methods operating on n vectors, 
* `VecVecg`: methods for n vectors, plus another generic argument, e.g. vector of weights.

In other words, the traits and their methods operate on arguments of their required categories. In classical statistical parlance, the main categories correspond to the number of 'random variables'. However, the vectors' end types (for the actual data) are mostly generic: usually some numeric type. There are also some traits specialised for input end types `f64` and `u8` and some that take mutable self. End type `f64` is most commonly used for the results.

### Documentation

For more detailed comments, plus some examples, see the source. You may have to unclick the 'implementations on foreign types' somewhere near the bottom of the page in the rust docs to get to it.  (Since these traits are implemented over the pre-existing Rust Vec type).

## Structs

* `struct Med` holds the median and quartiles

* `struct MStats` holds the mean and standard deviation

* `struct MinMax`, re-exported from crate `indxvec`. It holds min and max values of a vector and their indices. It is returned by function `indxvec::merge::minmax`.

##  Auxiliary Functions

* `i64tof64`: converts an i64 vector to f64, 
* `wsum`: sum of a sequence 1..n, also the size of a lower/upper triangular matrix below/above the diagonal (n*(n+1)/2.).

## Trait Stats

One dimensional statistical measures implemented for all numeric end types.

Its methods operate on one slice of generic data and take no arguments.
For example, `s.amean()` returns the arithmetic mean of the data in slice `s`.
Some of these methods are checked and will report all kinds of errors, such as an empty input. This means you have to apply to their results `?`, `.unwrap()` or something better.

Included in this trait are:

* means (arithmetic, geometric and harmonic),
* standard deviations,
* linearly weighted means (useful for time dependent data analysis),
* median and quartiles,
* probability density function (pdf)
* autocorrelation, entropy
* linear transformation to [0,1],
* other measures and vector algebra operators

## Trait MutStats

A few of the `Stats` methods are reimplemented under this trait
(only for f64), so that they mutate `self` in-place.
This is more efficient and convenient in some circumstances, such as in
vector iterative methods.

## Trait Vecg

Vector algebra operations between two slices `&[T]`, `&[U]` of any length (dimensionality):

* Vector additions, subtractions and products (scalar, kronecker, outer),
* Other relationships and measures of difference,
* Pearson's, Spearman's and Kendall's correlations,
* `Median correlation`, which we define analogously to Pearson's, as cosine of an angle between two zero median vectors (instead of his zero mean vectors).
* Joint pdf, joint entropy, statistical independence (based on mutual information).

This trait is unchecked (for speed), so some caution with data is advisable.

## Trait Vecf64

A handful of methods from `Vecg`, specialised to an argument of known end type `f64` or `&[f64]`.

## Traits MutVecg & MutVecf64

Mutable vector addition, subtraction and multiplication.  
Mutate `self` in-place.
This is for efficiency and convenience. Specifically, in
vector iterative methods.

`MutVecf64` is to be used in preference, when the end type of `self` is known to be `f64`. Beware that these methods work by side-effect and do not return anything, so they can not be functionally chained.

## Trait Vecu8

Some vector algebra as above that can be more efficient when the end type happens to be u8 (bytes). They have u8 appended to their names to avoid confusion with Vecg methods.

* Relationships between two vectors (of bytes)
* Frequency count of bytes by their values (histogram, pdf, jointpdf)
* Entropy, jointentropy, independence (different algorithms to those in Vecg)

## Trait VecVec

Relationships between n vectors (in d dimensions).
This general data domain is denoted here as (nd). It is in nd where the main original contribution of this library lies. True geometric median (gm) is found by fast and stable iteration, using improved Weiszfeld's algorithm `gmedian`. This algorithm solves Weiszfeld's convergence and stability problems in the neighbourhoods of existing set points.

* centroid, medoid, outliers, gm
* sums of distances, radius of a point (as its distance from gm)
* characterisation of a set of multidimensional points by the mean, standard deviation, median of its points' radii. These are useful recognition measures for the set.
* transformation to zero geometric median data,
* multivariate trend (regression) between two sets of nd points,
* covariance and comediance matrices (weighted and unweighted).

Warning: trait VecVec is entirely unchecked, so check your data upfront.

## Trait VecVecg

Methods which take an additional generic vector argument, such as a vector of weights for computing weighted geometric medians.

## Appendix I: Terminology

#### Including some new definitions for sets of nd points, i.e. n points in d dimensional space

* `Centroid/Centre/Mean` is the (generally non member) point that minimises the sum of *squares* of distances to all member points. Thus it is susceptible to outliers. Specifically, it is the n-dimensional arithmetic mean. By drawing physical analogy with gravity, it is sometimes called 'the centre of mass'. Centroid can also sometimes mean the member of the set which is the nearest to the Centre. Here we follow the common (if somewhat confusing) usage: Centroid = Centre = Arithmetic Mean.

* `Quasi/Marginal Median` is the point minimising sums of distances separately in each dimension (its coordinates are 1-d medians along each axis). It is a mistaken concept which we do not use here.

* `Tukey Median` is the point maximising `Tukey's Depth`, which is the minimum number of (outlying) points found in a hemisphere in any direction. Potentially useful concept but only partially implemented here by `tukeyvec`, as its advantages over the geometric median are not clear.

* `Median or the true geometric median (gm)`, is the point (generally non member), which minimises the sum of distances to all members. This is the one we want. It is much less susceptible to outliers than centroid. In addition, unlike quasi median, `gm` is rotation independent.

* `Medoid` is the member of the set with the least sum of distances to all other members. Equivalently, the member which is the nearest to the `gm`.

* `Outlier` is the member of the set with the greatest sum of distances to all other members. Equivalently, it is the point furthest from the `gm`.

* `Zero median vectors` are obtained by subtracting the `gm` (placing the origin of the coordinate system at the `gm`). This is a proposed  alternative to the commonly used `zero mean vectors`, obtained by subtracting the centroid.

* `MADGM` (median of distances from gm). This is a generalisation of `MAD` (median of absolute differences) measure from 1d to nd. It is a robust measure of data spread.

* `Comediance` is similar to `covariance`, except that zero median vectors are used to compute it,  instead of zero mean vectors.

* `Median correlation` between
 two vectors. We define it analogously to Pearson, as cosine of an angle between two 'normalised' vectors. Pearson 'normalises' by subtracting the mean from all components, we subtract the median.

## Appendix II: Recent Releases

* **Version 1.0.18** - Renamed `madn` to `madgm` (median of absolute deviations, i.e. radii, from gm). Added its weighted version `wmadgm`. They now take `gm` or `wgm` respectively as an argument, to avoid recomputation. Removed `radvec`, as it was a simple difference of `gm` and `centroid`.

* **Version 1.0.16** - Added `tukeyvec` and test of tukeyvec. Also changed usage of `ran` crate to its generic methods within `vecvec` test.

* **Version 1.0.14** - Some improvements of README.md.

* **Version 1.0.13** - Updated `ran` dev-dependency to "^0.3".

* **Version 1.0.12** - New random number generators are now in their own crate `ran`. It has been added here to development dependencies, where it properly belongs. `tests.rs` have been changed accordingly. No other changes.

* **Version 1.0.11** - The random number generators in indxvec have been been moved to their own module `random`. To keep `tests.rs` compatible, the import is now changed accordingly, to: `use indxvec::random::*;`

* **Version 1.0.10** Now using new random numbers generators from `indxvec` for testing.

* **Version 1.0.9** Removed `genvec` and `genvecu8` as non-essential. Removed some examples in the code that were using them. Changed the printing of vecs to utilise the new trait `Printing` from `indxvec`. See `testing.rs` for usage.

* **Version 1.0.8** Pruned some non-essential code, such as `smedian`. Gmedian now performs consistently a bit better.

* **Version 1.0.7** Achieved further 20% speedup of `gmedian` by optimising some inner loops.

* **Version 1.0.6** Added `crossfeatures` - computes relationships between all pairs of vectors in self. Returns flattened lower triangular (symmetric) matrix.

    Dependence of two vectors is now normalised to the range [0,1], e.g. the dependence of two identical vectors without repetitions is 1. Same for vectors of any real values that are all unique. In these cases it is better to fall back to correlations. N-dependence of a whole set of vectors will often be more than one.    

* **Version 1.0.5** Added 1D median correlation `medaincorr`. This is a more robust measure. Added `dependencies` and `correlations` which efficiently map these relationships of a single given vector (typically of outcomes), to a set of vectors (typically feature vectors).

* **Version 1.0.4** Added joint pdf, joint entropy and dependence for a set of n vectors.

* **Version 1.0.3** Better implementations of joint probability and joint entropy. Code style and testing improvements.

* **Version 1.0.2** Updated the dependency `indxvec` to version 1. A few minor changes to this document.

* **Version 1.0.1** Minor change: `sortedeccs` and `wsortedeccs` now take gm as an argument for more efficient repeated use. Vecvec test improved.

* **Version 1.0.0** Rstats reaches stability (of sorts)!
