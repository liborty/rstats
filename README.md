# Rstats - Rust Stats

![Crates.io](https://img.shields.io/crates/v/rstats?logo=rust) ![GitHub last commit (branch)](https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github)  

## Usage

Insert in your `Cargo.toml` file under `[dependencies]` `rstats = "^0.8"`
and in your source file(s) `use rstats::` followed by any of these functions and/or traits that you need: `{functions, Stats, MutStats, Vecg, MutVecg, VecVec, Vecu8, VecVecu8};`

## Introduction

Rstats is primarily about characterising multidimensional sets of points, with applications to Machine Learning and Data Analysis. It begins with statistical measures and vector algebra, which provide some basic self-contained tools for the more interesting algorithms but can also be used in their own right.

Our treatment of multidimensional sets of points is constructed from the first principles. Some original concepts, not found elsewhere, are introduced and implemented here. Specifically, the new multidimensional (geometric) median algorithm. Also, the `comediance matrix`; a replacement for the covariance matrix. It is obtained simply by supplying `covar` with the geometric median instead of the centroid.

*Zero median vectors are generally preferable to the commonly used zero mean vectors.*

Most authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. In a highly dimensional space, they are not even any easier to compute.

*Specifically, all such 1-d measures are sensitive to the choice of axis.*

Our methods based on the True Geometric Median, computed here by `gmedian`, are axis (rotation) independent from the first step.

### Implementation

The constituent parts of Rstats are Rust traits grouping together methods applicable to a single vector (`Stats`), two vectors (`Vecg`), or n vectors of data. End type f64 is most commonly used for the results.

### Documentation

Follow the documentation link. Then select a trait of interest to see the skeletal comments on the prototype function declarations in lib.rs. To see more detailed comments, plus some examples from the implementation files, scroll to the bottom of the trait and unclick [+] to the left of the `implementations` of the trait. To see the tests, consult `tests/tests.rs`.

To run the tests, use single thread. It will be slower but will produce the results in the right order:  
`cargo test --release -- --test-threads=1 --nocapture --color always`

## Structs and functions

* `pub struct Med` to hold median and quartiles

* `pub struct MStats` to hold a mean and standard deviation

* functions: `tof64, i64tof64, wsum, genvec, genvecu8` (see documentation for the module `functions.rs`).

## Traits

### Stats

One dimensional statistical measures implemented for all 'numeric' types.

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

* Vector additions, subtractions, products and other relationships and measures.
* Pearson's, Spearman's and Kendall's correlations.

This trait is unchecked (for speed), so some caution with data is advisable.

### MutVectors

Mutable vector operations that take one generic argument.  
A few of the essential `Vecg` methods are reimplemented here 
to mutate `self` (only `&[f64]`) in-place. 
This is for efficiency and convenience. For example, in
vector iterative methods.
Clearly, they can only be applied to a mutable variable. Beware that they work by side-effect and do not return anything, so they can not be chained.

### Vecu8

* Some vector algebra as above that can be made more efficient when the end type happens to be u8 (bytes).
* Frequency count of bytes by their values (Histogram or Probability Density Function).
* Entropy measures in units of e (using natural logarithms).

### VecVec

Relationships of one vector to a set of vectors:

* sums of distances, eccentricity,
* centroid, medoid, true geometric median,
* transformation to zero (geometric) median data,
* relationship between sets of multidimensional vectors: trend,
* covariance and comediance matrices (weighted and unweighted).

Trait VecVec is entirely unchecked, so check your data upfront. This is the more sophisticated part of the library. The true geometric median is found iteratively.

### VecVecu8

Some of the above for vectors of vectors of bytes.

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

* **Version 0.8.6** Added `comed` and `wcomed` methods to `VecVec` trait.

* **Version 0.8.5** Split `MutVectors` trait into `MutStats` (with no arguments) and `MutVecg` (with one generic argument). They are both still implemented only for f64 and will remain so. However, it is now possible, for example, to mutably subtract a slice of any end type. This allowed the deletion of `mutvaddu8` and `mutvsubu8` as special cases. Fixed some in-code tests that were not yet using the new `Vecg` trait. Trait `VecVecf64` generalised and renamed to `VecVec`.

* **Version 0.8.4** Significant reorganisation. `Vecf64` trait and its source module `vecf64.rs` have been replaced by `Vecg` generic trait and `vecg.rs` module respectively. Numerous methods have been sorted more carefully into `Vecg` trait or `Stats` trait, according to whether or not they take an argument. Some methods have also been moved out of `Vecu8` trait and generalised in the process. Methods remaining in `Vecu8` now all have names ending in `u8` for clarity and to avoid confusion with their generic versions. Some bugs in entropy methods have been fixed.

* **Version 0.8.3** Simplification of generic `Stats`. `GSlice` is no longer needed. The only restriction remaining is the necessity to explicitly convert `&[i64] -> &[f64]`, using function `statsg::i64tof64(s: &[i64])`. All other end types are fine. This made possible the removal of two modules, `statsf64.rs` and `stasi64.rs`. They are now superceded by a single generic `statsg.rs`. This rationalisation work will continue with the remaining traits as well.

* **Version 0.8.2** Added `statsgen.rs` (generic) module to add the capability of applying the trait `Stats` to all numeric end types, as long as their slices are wrapped in `GSlice(&s)`. This is a step towards more generality, as `Stats` methods can now work on all primitive numeric types.  f64 and i64 remain as previously, so they should not be wrapped.

* **Version 0.8.0** Simplified, more stable version. Moved auxiliary macro `here` and functions `wv,wi` to crate `indxvec`. Tidied up the tests accordingly.

* **Version 0.7.17** Updated Cargo.toml dependency to `indxvec = "^0.2"`.
