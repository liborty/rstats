# Rstats

[<img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github">](https://github.com/liborty/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/v/rstats?logo=rust">](https://crates.io/crates/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/d/rstats?logo=rust">](https://crates.io/crates/rstats)
[<img alt="docs.rs" src="https://img.shields.io/docsrs/rstats?logo=rust">](https://docs.rs/rstats)

## Usage

Insert `rstats = "^1"` in the `Cargo.toml` file, under `[dependencies]`.

Use in source files any of the following structs, as needed:  
`use rstats::{MinMax,Med,Mstats};`  
and any of the following helper functions:  
`use rstats::{i64tof64,tof64,here,wi,wv,wsum,printvv,genvec,genvecu8};`  
and any of the following traits:  
`use rstats::{Stats,MutStats,Vecu8,Vecf64,Vecg,MutVecg,VecVec,VecVecg};`

It is highly recommended to read and run `tests/tests.rs`, which shows examples of usage.

To run all the tests, use single thread in order to produce the results in the right order:  
`cargo test --release -- --test-threads=1 --nocapture --color always`

## Introduction

`Rstats` is primarily about characterising multidimensional sets of points, with applications to Machine Learning and Big Data Analysis. It uses `non analytical statistics`, where the 'random variables' are replaced by vectors of real data. Probabilities densities and other parameters are always obtained from the data, not from some assumed distributions.

This crate begins with basic statistical measures and vector algebra, which provide self-contained tools for the multidimensional algorithms but can also be used in their own right.

Our treatment of multidimensional sets of points (vectors) is constructed from the first principles. Some original concepts, not found elsewhere, are introduced and implemented here:

* `gmedian` - fast multidimensional (geometric) median algorithm. 

* `comediance`, is a suggested replacement for covariance (matrix). It is obtained simply by supplying `covar` with the geometric median instead of the usual centroid. Thus zero median vectors are replacing zero mean vectors.

* similarly, in just one dimension, our `mediancorr` is to replace Pearson's correlation. We define median correlation  as cosine of an angle between two zero median vectors (instead of zero mean vectors as per Pearson). 

*Zero median vectors are generally preferable to the commonly used zero mean vectors.*

In n dimensions, many authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. In a highly dimensional space, they are not even any faster to compute.

*Specifically, all 1-d measures are sensitive to the choice of axis and thus are affected by rotation.*

In contrast, analyses based on the true geometric median (gm) are axis (rotation) independent. Also, they are more stable, as medians have a 50% breakdown point (the maximum possible). They are computed here by  `gmedian` and its weighted version `wgmedian`.

## Implementation

The main constituent parts of Rstats are its traits. The selection of traits (to import) is primarily determined by the types of objects to be handled. These are mostly vectors of arbitrary length (dimensionality). The main traits are implementing methods applicable to a single vector (of numbers) - `Stats`, methods (of vector algebra) for two vectors - `Vecg`, methods for n vectors - `VecVec`, and methods for n vectors with another generic argument - `VecVecg`.

In other words, the traits and their methods operate on arguments of their required categories. In classical statistical parlance, the main categories correspond to the number of 'random variables'. However, the vectors' end types (for the actual data) are mostly generic: usually some numeric type. There are also some traits specialised for input end types `f64` and `u8` and some that take mutable self. End type `f64` is most commonly used for the results.

### Documentation

For more detailed comments, plus some examples, see the source. You may have to unclick the 'implementations on foreign types' somewhere near the bottom of the page in the rust docs to get to it.

## Structs and auxiliary functions

* `struct Med` to hold median and quartiles

* `struct MStats` to hold mean and standard deviation

* `struct MinMax` re exported from crate `indxvec` to hold min and max values of a vector and their indices. It is returned by function `indxvec::merge::minmax`.

* auxiliary functions: `i64tof64, tof64, here, wsum, wi, wv, printvv, genvec, genvecu8`.

## Trait Stats

One dimensional statistical measures implemented for all numeric end types.

Its methods operate on one slice of generic data and take no arguments.
For example, `s.amean()` returns the arithmetic mean of the data in slice `s`.
Some of these methods are checked and will report all kinds of errors, such as an empty input. This means you have to call `.unwrap()` or something better on their results.

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
* `Median correlation`, which we define analogously to Pearson's, as cosine of an angle between two zero median vectors (instead of zero mean vectors).
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

Relationships between n vectors (in d dimensions). This is the main original contribution of this library. True geometric median is found by fast and stable iteration, using improved Weiszfeld's algorithm `gmedian`. This algorithm solves Weiszfeld's convergence and stability problems in the neighbourhood of existing set points.

* sums of distances, eccentricity (radius) measure,
* centroid, medoid, outliers, true geometric median,
* characterisations of sets of multidimensional points (of d random variables): means, stds, medians 
* transformation to zero geometric median data,
* multivariate trend (regression) between two sets of nd points,
* covariance and comediance matrices (weighted and unweighted).

Trait VecVec is entirely unchecked, so check your data upfront.

## Trait VecVecg

Methods which take an additional generic vector argument, such as a vector of weights for computing the weighted geometric medians.

## Appendix I: Terminology (and some new definitions) for sets of nD points

* `Centroid\Centre\Mean` is the (generally non member) point that minimises the sum of *squares* of distances to all member points. Thus it is susceptible to outliers. Specifically, it is the n-dimensional arithmetic mean. By drawing physical analogy with gravity, it is sometimes called 'the centre of mass'. Centroid can also sometimes mean the member of the set which is the nearest to the Centre. Here we follow the common (if somewhat confusing) usage: Centroid = Centre = Arithmetic Mean.

* `Quasi\Marginal Median` is the point minimising sums of distances separately in each dimension (its coordinates are 1-d medians along each axis). It is a mistaken concept which we do not use here.

* `Tukey Median` is the point maximising `Tukey's Depth`, which is the minimum number of (outlying) points found in a hemisphere in any direction. Potentially useful concept but not yet implemented here, as its advantages over `gm` are not clear.

* `Medoid` is the member of the set with the least sum of distances to all other members.

* `Outlier` is the member of the set with the greatest sum of distances to all other members.

* `Median or the true geometric median (gm)`, is the point (generally non member), which minimises the sum of distances to all members. This is the one we want. It is much less susceptible to outliers than centroid. In addition, unlike quasi median, `gm` is rotation independent.

* `Zero median vector` is obtained by subtracting the geometric median. This is a proposed  alternative to the commonly used `zero mean vector`, obtained by subtracting the centroid.

* `Comediance` is similar to `covariance`, except zero median vectors are used to compute it,  instead of zero mean vectors.

* `Median correlation` we define analogously to Pearson, as cosine of an angle between two vectors 'normalised' by subtracting their 1d medians from all components, instead of subtracting their means. 

## Appendix II: Recent Releases

* **Version 1.0.8** Pruned some non-essential code.

* **Version 1.0.7** Achieved further 20% speedup of `gmedian` by optimising some inner loops. 

* **Version 1.0.6** Added `crossfeatures` - computes relationships between all pairs of vectors in self. Returns flattened lower triangular (symmetric) matrix.

    Dependence of two vectors is now normalised to the range [0,1], e.g. the dependence of two identical vectors without repetitions is 1. Same for vectors of any real values that are all unique. In these cases it is better to fall back to correlations. N-dependence of a whole set of vectors will often be more than one.    

* **Version 1.0.5** Added 1D median correlation `medaincorr`. This is a more robust measure. Added `dependencies` and `correlations` which efficiently map these relationships of a single given vector (typically of outcomes), to a set of vectors (typically feature vectors).

* **Version 1.0.4** Added joint pdf, joint entropy and dependence for a set of n vectors.

* **Version 1.0.3** Better implementations of joint probability and joint entropy. Code style and testing improvements.

* **Version 1.0.2** Updated the dependency `indxvec` to version 1. A few minor changes to this document.

* **Version 1.0.1** Minor change: `sortedeccs` and `wsortedeccs` now take gm as an argument for more efficient repeated use. Vecvec test improved.

* **Version 1.0.0** Rstats reaches stability (of sorts)!

* **Version 0.9.4** Organisation improvements. Added trait `Vecf64` and moved into it relevant methods from `Vecg`. Added a few functions to MutVecf64 trait. Simplified `gmedian`.

* **Version 0.9.3** Added `hwmeanstd` - harmonic weighted mean and standard deviation. Tidied up readme badges and some tests. Simplified random number generation. Weights for the weighted means are now ascending (more intuitive).

* **Version 0.9.2** Fixed some tests broken by moved functions. Added harmonic standard deviation, harmonic centroid  and more tests.

* **Version 0.9.1** Made the auxiliary functions more visible by moving them to `lib.rs` (the top level of the crate).

* **Version 0.9.0** Added `kron` and `outer` products to `Vecg` trait.
