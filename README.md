# Rstats

[<img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github">](https://github.com/liborty/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/v/rstats?logo=rust">](https://crates.io/crates/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/d/rstats?logo=rust">](https://crates.io/crates/rstats)
[<img alt="docs.rs" src="https://img.shields.io/docsrs/rstats?logo=rust">](https://docs.rs/rstats)

## Usage

Insert `rstats = "^1"` in the `Cargo.toml` file, under `[dependencies]`.

Use in your source files any of the following structs, if needed:

```rust  
use rstats::{Mstats};  
use indxvec::{MinMax};  
use medians::Med;
```

and any of the following rstats traits:
```rust 
use rstats::{ Stats, Vecg, Vecu8, MutVecg, VecVec, VecVecg };
```

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

* `median correlation`- in one dimension (1d), our `mediancorr` method is to replace *Pearson's correlation*. We define *median correlation*  as cosine of an angle between two zero median vectors (instead of Pearson's zero mean vectors).

* `gmedian` - fast multidimensional geometric median (gm) algorithm.

* `madgm` - generalisation of robust data spread estimator known as 'MAD' (median of absolute deviations from median),  from 1d to nd. 

* `contribution` - of a point to an nd set. Defined as gm displacement if the point was added/removed. Related to the point's radius but not the same, as it depends on all the points.

* `comediance` - instead of covariance (matrix). It is obtained by supplying `covar` with the geometric median instead of the usual centroid. Thus *zero median vectors* are replacing *zero mean vectors* in covariance calculations.

*Zero median vectors are generally preferable to the commonly used zero mean vectors.*

In n dimensions (nd), many authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. In a highly dimensional space, they are not even any faster to compute.

*Specifically, all such 1d measures are sensitive to the choice of axis and thus are affected by rotation.*

In contrast, analyses based on the true geometric median (gm) are axis (rotation) independent. Also, they are more stable, as medians have a 50% breakdown point (the maximum possible). They are computed here by methods `gmedian` and its weighted version `wgmedian`, in traits `vecvec` and `vecvecg` respectively.

## Implementation

The main constituent parts of Rstats are its traits. The selection of traits (to import) is primarily determined by the types of objects to be handled. These are mostly vectors of arbitrary length (dimensionality). The main traits are implementing methods applicable to:

* `Stats`: a single vector (of numbers),
* `Vecg`: methods (of vector algebra) operating on two vectors, e.g. scalar product
* `MutVecg`: some of the above methods, mutating self
* `VecVec`: methods operating on n vectors, 
* `VecVecg`: methods for n vectors, plus another generic argument, e.g. vector of weights.

In other words, the traits and their methods operate on arguments of their required categories. In classical statistical parlance, the main categories correspond to the number of 'random variables'. However, the vectors' end types (for the actual data) are mostly generic: usually some numeric type. There are also some traits specialised for input end type `u8` and some that take mutable self. End type `f64` is most commonly used for the results.

### Documentation

For more detailed comments, plus some examples, see the source. You may have to unclick the 'implementations on foreign types' somewhere near the bottom of the page in the rust docs to get to it.  (Since these traits are implemented over the pre-existing Rust Vec type).

## Struct

* `struct MStats` holds the central tendency, e.g. mean, and spread, e.g. standard deviation.

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
* probability density function (pdf)
* autocorrelation, entropy
* linear transformation to [0,1],
* other measures and vector algebra operators

Note that fast implementation of 1d medians is as of version 1.1.0 in crate `medians`:  
`use medians::{Med,Median};`

## Trait Vecg

Generic vector algebra operations between two slices `&[T]`, `&[U]` of any length (dimensionality). It may be necessary to invoke some using the 'turbofish' `::<type>` syntax to indicate the type U of the supplied argument, e.g.:  
`datavec.as_slice().methodname::<f64>(arg)`  
This is because Rust is currently incapable of inferring the type ('the inference bug').

* Vector additions, subtractions and products (scalar, kronecker, outer),
* Other relationships and measures of difference,
* Pearson's, Spearman's and Kendall's correlations,
* `Median correlation`, which we define analogously to Pearson's, as cosine of an angle between two zero median vectors (instead of his zero mean vectors).
* Joint pdf, joint entropy, statistical independence (based on mutual information).

This trait is unchecked (for speed), so some caution with data is advisable.

## Trait MutVecg

A select few of the `Stats` and `Vecg` methods (e.g. mutable vector addition, subtraction and multiplication) are reimplemented under this trait, so that they can mutate `self` in-place. This is more efficient and convenient in some circumstances, such as in vector iterative methods.

## Trait Vecu8

Some vector algebra as above that can be more efficient when the end type happens to be u8 (bytes). These methods have u8 appended to their names to avoid confusion with Vecg methods. These specific algorithms are different to their generic equivalents in Vecg.

* Frequency count of bytes by their values (histogram, pdf, jointpdf)
* Entropy, jointentropy, independence.

## Trait VecVec

Relationships between n vectors (in d dimensions).
This general data domain is denoted here as (nd). It is in nd where the main original contribution of this library lies. True geometric median (gm) is found by fast and stable iteration, using improved Weiszfeld's algorithm `gmedian`. This algorithm solves Weiszfeld's convergence and stability problems in the neighbourhoods of existing set points.

* centroid, medoid, outliers, gm
* sums of distances, radius of a point (as its distance from gm)
* characterisation of a set of multidimensional points by the mean, standard deviation, median and MAD of its points' radii. These are useful recognition measures for the set.
* transformation to zero geometric median data,
* multivariate trend (regression) between two sets of nd points,
* covariance and comediance matrices.

Warning: trait VecVec is entirely unchecked, so check your data upfront.

## Trait VecVecg

Methods which take an additional generic vector argument, such as a vector of weights for computing weighted geometric medians (where each point has its own weight).

## Appendix I: Terminology

#### Including some new definitions for sets of nd points, i.e. n points in d dimensional space

* `Median correlation` between
 two 1d vectors. We define it analogously to Pearson, as cosine of an angle between two 'normalised' vectors. Pearson 'normalises' by subtracting the mean from all components, we subtract the median.

* `Centroid/Centre/Mean` is the (generally non member) nd point that minimises the sum of *squares* of distances to all other member points. Thus it is susceptible to outliers. Specifically, it is the n-dimensional arithmetic mean. It is sometimes called 'the centre of mass'. Centroid can also sometimes mean the member of the set which is the nearest to the Centre. Here we follow the common usage: Centroid = Centre = Arithmetic Mean.

* `Quasi/Marginal Median` is the point minimising sums of distances separately in each dimension (its coordinates are 1-d medians along each axis). It is a mistaken concept which we do not recommend using.

* `Tukey Median` is the point maximising `Tukey's Depth`, which is the minimum number of (outlying) points found in a hemisphere in any direction. Potentially useful concept but only partially implemented here by `tukeyvec`, as its advantages over the geometric median are not clear.

* `Median or the true geometric median (gm)`, is the point (generally non member), which minimises the sum of distances to all members. This is the one we want. It is much less susceptible to outliers than centroid. In addition, unlike quasi median, `gm` is rotation independent.

* `Medoid` is the member of the set with the least sum of distances to all other members. Equivalently, the member which is the nearest to the `gm`.

* `Outlier` is the member of the set with the greatest sum of distances to all other members. Equivalently, it is the point furthest from the `gm`.

* `Zero median vectors` are obtained by subtracting the `gm` (placing the origin of the coordinate system at the `gm`). This is a proposed  alternative to the commonly used `zero mean vectors`, obtained by subtracting the centroid.

* `MADGM` (median of distances from gm). This is a generalisation of `MAD` (median of absolute differences) measure from 1d to nd. It is a robust measure of nd data spread.

* `Comediance` is similar to `covariance`, except that zero median vectors are used to compute it,  instead of zero mean vectors. By Cholesky singular value decomposition of this positive definite matrix, it is possible to calculate *Mahalanobis distance* (weighted distace, where the weights are derived from the shape of the data points cloud). 

* `Contribution`: one of the questions of interest to Machine Learning (ML) is how to quantify the significance of the contribution that each example point (typically a member of some large nd set) makes to the recognition concept, or class, represented by that set. In answer to this, we define `the contribution` of a point as the change to gm caused by adding/deleting that point. Generally more outlying points make greater contributions but not as much as is the case with means. The contribution depends on the arrangement of other set points as well.

## Appendix II: Recent Releases

* **Version 1.1.6** - Moved `contrib_vec_newp`, `contrib_newpt`, `contribvec.oldpt`, `contrib.oldpt` to more appropriate trait Vecg.

* **Version 1.1.5** - Added `contribution` methods to `VecVec`: `contribvec_newp`, `contrib_newpt`, `contribvec.oldpt`, `contrib.oldpt`. Named according to whether they return a vector change to gm, or just a magnitude of change measure and whether it is a new point, or an existing point.

* **Version 1.1.4** - Speedup: invented faster termination condition for gm algorithm! Generality: invented `contribution` measure for nd set points. Tidied up functions associated with gm. Deleted and added some.

* **Version 1.1.3** - Added `quasimedian`, `gmerror`. Some more tidying up.

* **Version 1.1.2** - Tidying up of methods in trait VecVec. Removed some, added some. Some preparations for parallelism.

* **Version 1.1.1** - Crate size reduction. Some more code pruning, this time of methods in trait `Vecu8`, which are performed perfectly adequately by their generic counterparts in `Vecg`. Leaving only those with different u8 specific algorithms. 

* **Version 1.1.0** - Big release. Added dependency on crate `medians` for fast 1D medians. Simplifications: subsumed module `mutstats.rs` into `mutvec.rs`. Removed traits `Mutstats` and `MutVecf64` and added their few methods to trait `MutVecg`. Added some more doc comments here. Generalisations: methods in `Vecg` and `MutVecg` now work on any type T of self and a potentially different type U for their argument. They should be called with the 'turbofish' syntax.
