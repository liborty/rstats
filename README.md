# Rstats

[<img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github" height="20">](https://github.com/liborty/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/v/rstats.svg?style=for-the-badge&color=green&logo=rust" height="20">](https://crates.io/crates/rstats)
[<img alt="crates.io" src="https://img.shields.io/crates/d/rstats?logo=rust" height="20">](https://crates.io/crates/rstats)
[<img alt="docs.rs" src="https://img.shields.io/badge/docs.rs-rstats-green?style=for-the-badge&labelColor=555555&logoColor=white&logo=data:image/svg+xml;base64,PHN2ZyByb2xlPSJpbWciIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgdmlld0JveD0iMCAwIDUxMiA1MTIiPjxwYXRoIGZpbGw9IiNmNWY1ZjUiIGQ9Ik00ODguNiAyNTAuMkwzOTIgMjE0VjEwNS41YzAtMTUtOS4zLTI4LjQtMjMuNC0zMy43bC0xMDAtMzcuNWMtOC4xLTMuMS0xNy4xLTMuMS0yNS4zIDBsLTEwMCAzNy41Yy0xNC4xIDUuMy0yMy40IDE4LjctMjMuNCAzMy43VjIxNGwtOTYuNiAzNi4yQzkuMyAyNTUuNSAwIDI2OC45IDAgMjgzLjlWMzk0YzAgMTMuNiA3LjcgMjYuMSAxOS45IDMyLjJsMTAwIDUwYzEwLjEgNS4xIDIyLjEgNS4xIDMyLjIgMGwxMDMuOS01MiAxMDMuOSA1MmMxMC4xIDUuMSAyMi4xIDUuMSAzMi4yIDBsMTAwLTUwYzEyLjItNi4xIDE5LjktMTguNiAxOS45LTMyLjJWMjgzLjljMC0xNS05LjMtMjguNC0yMy40LTMzLjd6TTM1OCAyMTQuOGwtODUgMzEuOXYtNjguMmw4NS0zN3Y3My4zek0xNTQgMTA0LjFsMTAyLTM4LjIgMTAyIDM4LjJ2LjZsLTEwMiA0MS40LTEwMi00MS40di0uNnptODQgMjkxLjFsLTg1IDQyLjV2LTc5LjFsODUtMzguOHY3NS40em0wLTExMmwtMTAyIDQxLjQtMTAyLTQxLjR2LS42bDEwMi0zOC4yIDEwMiAzOC4ydi42em0yNDAgMTEybC04NSA0Mi41di03OS4xbDg1LTM4Ljh2NzUuNHptMC0xMTJsLTEwMiA0MS40LTEwMi00MS40di0uNmwxMDItMzguMiAxMDIgMzguMnYuNnoiPjwvcGF0aD48L3N2Zz4K" height="20">](https://docs.rs/rstats/)

## Usage

Insert `rstats = "^0.9"` in the `Cargo.toml` file under `[dependencies]`.

Use any of the following structs that you need in your source files:  
`use rstats::{MinMax,Med,Mstats};`  
Use any of the following functions that you need in your source files:  
`use rstats::{i64tof64,tof64,here,wi,wv,wsum,printvv,genvec,genvecu8};`  
Use any of the following traits that you need in your source files:  
`use rstats::{Stats,MutStats,Vecu8,Vecg,MutVecg,VecVec,VecVecg};`

## Introduction

`Rstats` is primarily about characterising multidimensional sets of points, with applications to Machine Learning and Big Data Analysis. It begins with basic statistical measures and vector algebra, which provide self-contained tools for the multidimensional algorithms but can also be used in their own right.

Our treatment of multidimensional sets of points is constructed from the first principles. Some original concepts, not found elsewhere, are introduced and implemented here. Specifically, the new multidimensional (geometric) median algorithm. Also, the `comediance matrix`  as a replacement for the covariance matrix. It is obtained simply by supplying `covar` with the geometric median instead of the usual centroid (mean vector).

*Zero median vectors are generally preferable to the commonly used zero mean vectors.*

Most authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. In a highly dimensional space, they are not even any easier to compute.

*Specifically, all 1-d measures are sensitive to the choice of axis and thus are affected by rotation.*

In contrast, analyses based on the true geometric median (gm), computed here by the novel methods `gmedian` and `wgmedian`, are axis (rotation) independent.

### Implementation

The main constituent parts of Rstats are Rust traits grouping together methods applicable to a single vector (of numbers) - `Stats`, two vectors - `Vecg`, or n vectors - `VecVec` and `VecVecg`. End type `f64` is most commonly used for the results, whereas the inputs to the generic methods can be vectors (or slices) of any numeric end types.

### Documentation

To see more detailed comments, plus some examples, see the source.  
It is highly recommended to read and run `tests/tests.rs`, which shows examples of usage.

To run all the tests, use single thread in order to produce the results in the right order:  
`cargo test --release -- --test-threads=1 --nocapture --color always`

## Structs and auxiliary functions

* `struct Med` to hold median and quartiles

* `struct MStats` to hold mean and standard deviation

* `struct MinMax` reexported from crate `indxvec` to hold min and max values of a vector and their indices. It is returned by function `indxvec::merge::minmax`.

* functions: `i64tof64,tof64,here,wsum,wi,wv,printvv,genvec,genvecu8`

## Traits

### Stats

One dimensional statistical measures implemented for all numeric end types.

Its methods operate on one slice of generic data and take no arguments.
For example, `s.amean()` returns the arithmetic mean of the data in slice `s`.
Some of these methods are checked and will report all kinds of errors, such as an empty input. This means you have to call `.unwrap()` or something similar on their  results.

Included in this trait are:

* means (arithmetic, geometric and harmonic),
* standard deviations,
* linearly weighted means (useful for time dependent data analysis),
* median and quartiles,
* autocorrelation, entropy
* linear transformation to [0,1],
* other measures and vector algebra operators

### MutStats

A few of the `Stats` methods are reimplemented under this trait
(only for f64), so that they mutate `self` in-place.
This is more efficient and convenient in some circumstances, such as in
vector iterative methods.

### Vecg

Vector algebra operations between two slices `&[T]`, `&[U]` of any length (dimensionality):

* Vector additions, subtractions and products (scalar, kronecker, outer),
* Other relationships and measures,
* Pearson's, Spearman's and Kendall's correlations.

This trait is unchecked (for speed), so some caution with data is advisable.

### MutVecg & MutVecf64

Mutable vector addition, subtraction and multiplication.  
Mutate `self` in-place.
This is for efficiency and convenience. Specifically, in
vector iterative methods. `MutVecf64` is to be used in preference, when the end type of `self` is known to be `f64`. Beware that these methods work by side-effect and do not return anything, so they can not be functionally chained.

### Vecu8

* Some vector algebra as above that can be more efficient when the end type happens to be u8 (bytes).
* Frequency count of bytes by their values (Histogram or Probability Density Function).
* Entropy measures in units of e (using natural logarithms).

### VecVec

Relationships between n vectors (nD). This is the original contribution of this library. True geometric median is found by fast and stable iteration, using improved Weiszfeld's algorithm boosted by multidimensional secant method.

* sums of distances, eccentricity measure for nD points,
* centroid, medoid, outliers, true geometric median,
* transformation to zero (geometric) median data,
* relationship between two sets of multidimensional vectors: trend,
* covariance and comediance matrices (weighted and unweighted).

Trait VecVec is entirely unchecked, so check your data upfront.

### VecVecg

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

## Appendix II: Recent Releases

* **Version 0.9.1** Made the auxiliary functions more visible by moving them to `lib.rs` (the top level of the crate).

* **Version 0.9.0** Added `kron` and `outer` products to `Vecg` trait. Added `printvv` utility for pretty printing generic vectors of vectors.

* **Version 0.8.9** Minor improvements to `readme` and `vecvecg`.

* **Version 0.8.8** More generics: added `VecVecg` trait and removed `VecVecu8` as its functionality is now subsumed by this addition. Removed `benches/benchmark.rs` as it was not really needed. There are timings of the geometric medians computation in `tests/tests.rs`.

* **Version 0.8.7** Some simplification of reporting, using struct MinMax from crate `indxvec`.

* **Version 0.8.6** Added `comed` and `wcomed` methods to `VecVec` trait.
