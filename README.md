# Rstats - Rust Stats

![Crates.io](https://img.shields.io/crates/v/rstats?logo=rust) ![GitHub last commit (branch)](https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github)  

## Usage

Insert into your Cargo.toml file [dependencies] section:

```rust
rstats = "^0.8" 
```

and import into your source file(s) any of these functions and/or traits that you want:

```rust
use rstats::{functions,Stats,Vecf64,Vecu8,VecVecf64,VecVecu8,Mutvectors};
```

## Introduction

Rstats is primarily about characterising multidimensional sets of points, with applications to Machine Learning and Data Analysis. It begins with statistical measures and vector algebra, which provide some basic self-contained tools for the more interesting algorithms but can also be used in their own right.

Our treatment of multidimensional sets of points is constructed from the first principles. Some original concepts, not found elsewhere, are introduced and implemented here. Specifically, the new multidimensional (geometric) median algorithm. Also, the `comediance matrix`; a replacement for the covariance matrix. It is obtained simply by supplying `covar` with the geometric median instead of the centroid.

*Zero median vectors are generally preferable to the commonly used zero mean vectors.*

Most authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians are a poor start to stable characterisation of multidimensional data. In a highly dimensional space, they are not even any easier to compute.

*Specifically, all such 1-d measures are sensitive to the choice of axis.*

Our methods based on the True Geometric Median, computed here by `gmedian`, are axis (rotation) independent from the first step.

### Implementation

The constituent parts of Rstats are Rust traits grouping together functions applicable to vectors of data of relevant end types.End type f64 is most commonly used. Facilities for other end types are limited. For lots of data of other end types, it is always possible to clone to f64, see for example the included utility function `vecu8asvecf64`.

### Documentation

Follow the documentation link. Then select a trait of interest to see the skeletal comments on the prototype function declarations in lib.rs. To see more detailed comments, plus some examples from the implementation files, scroll to the bottom of the trait and unclick [+] to the left of the `implementations` of the trait. To see the tests, consult `tests.rs`.

To run the tests, use single thread. It will be slower but will produce the results in the right order:

```rust
cargo test --release -- --test-threads=1 --nocapture --color always
```

## Structs and functions

* pub struct Med to hold median and quartiles

* pub struct MStats to hold a mean and standard deviation

* functions wsum, genvec, genvecu8 (see documentation for the module `functions`).

## Traits

### Stats

One dimensional statistical measures implemented for `&[i64]` and `&[f64]`.

All these methods operate on one vector of data and take no arguments.
For example, `s.amean()` returns the arithmetic mean of slice `s` of either type.
This is the only attempt at genericity.  
This trait is carefully checked and will report all kinds of errors, such as empty input.
This means you have to call `.unwrap()` or something similar on its results.

Included in this trait are:

* means (arithmetic, geometric and harmonic),
* standard deviations,
* linearly weighted means (useful for time dependent data analysis),
* median and quartiles.

### Vecf64

Vector algebra implemented on one or two `&[f64]` slices of any length (dimensionality):

* Vector algebra
* Autocorrelation, Pearson's, Spearman's and Kendall's correlations.
* Linear transformation to [0,1], etc.

This trait is sometimes unchecked (for speed), so some caution with data is advisable.

### Vecu8

* Some vector algebra as above for vectors of u8 (bytes).
* Frequency count of bytes by their values (Histogram or Probability Density Function).
* Entropy measures in units of e (using natural logarithms).

### MutVectors

Some of the above functions are for memory efficiency reasons reimplemented in this trait so that they mutate `self` in place, instead of creating a new Vec. Clearly, they can only be applied to a mutable variable. They are useful in vector iterative methods. Beware that they work by side-effect and do not return anything, so they can not be chained.

### VecVec

Relationships of one vector to a set of vectors (of `&[f64]` end types):

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

* `Tukey Median` is the point maximising `Tukey's Depth`, which is the minimum number of (outlying) points found in a hemisphere in any direction. Potentially useful concept but not yet implemented here, as its advantages over GM are not clear.

* `Medoid` is the member of the set with the least sum of distances to all other members.

* `Outlier` is the member of the set with the greatest sum of distances to all other members.

* `Median or the true geometric median (gm)`, is the point (generally non member), which minimises the sum of distances to all members. This is the one we want. It is much less susceptible to outliers and is rotation independent.

* `Zero median vector` is obtained by subtracting the geometric median. This is a proposed  alternative to the commonly used `zero mean vector`, obtained by subtracting the centroid.

* `Comediance` is similar to covariance, except zero median vectors are used to compute it  instead of zero mean vectors.

## Appendix II: Recent Releases

* **Version 0.8.0** Simplified, more stable version. Moved auxiliary macro `here` and functions `wv,wi` to crate `indxvec`. Tidied up the tests accordingly.

* **Version 0.7.17** Updated Cargo.toml dependency to `indxvec = "^0.2"`.

* **Version 0.7.16** Simplification & generalisation: changed `sortf`, `ranks`, `iranks`, to `sortm`,`rank` from `indxvec`. Made `minmax` into generic function and moved it to `indxvec::merge`.

* **Version 0.7.15** Added `wcovar` also for `Vec<f64>`.

* **Version 0.7.14** Added weighted centroid `wacentroid`. Updated weighted geometric median `wgmedian`. Updated `gmedian` and `wgmedian` in `vecvecu8` to include the optimisations of v. 0.7.11. Ensured compatibility with the latest release of crate `indxvec`.

* **Version 0.7.12** Split off Index trait and associated functions into a new crate `indxvec`.

* **Version 0.7.11** Removed Kazutsugi (too specialised). Added `gcentroid` (geometric centroid). Further optimisations to `gmedian`.

* **Version 0.7.10** Added `symmatrix` to reconstruct full symmetric matrix from its lower triangular part (for compatibility with crates which duplicate data). Renamed mergerank to plain `rank` and added boolean argument to facilitate ranking in ascending or descending order. Expanded vecf64() tests (see it for instructive example usage).

* **Version 0.7.9** Added `wcovar` of weighted points. Improved struct GV and tests. Replaced `emsg` with macro `here!()` for easier diagnostics. Moved all structs into lib.rs.

* **Version 0.7.8** Added `covar` = covariance or comediance matrix computations. Some changes to this text (Readme.md).

* **Version 0.7.7** Fixed `merge_immutable` and added a test. Added `cityblockd` and `vaddu8`.

* **Version 0.7.6** Added `merge_immutable` and `merge_indices`. Simplified `mergesort`.

* **Version 0.7.5** Renamed VecVec trait to VecVecf64 to make the naming consistent. Added `unindexu8`. Removed `wsortedeccs` and `wsortedcos` for being too application specific.

* **Version 0.7.4** Added merge of two sorted &[f64]. Added `ascending` boolean flag to `unindex`, `sortm` and functions that call them, to facilitate easy sorting in ascending or descending order. Added `genvecu8` to `functions` to generate sets of random u8 vectors. Normalised cummulative probability density functions  to [0,1].

* **Version 0.7.3** Replaced varc with vector similarity and dissimilarity in [0,1] in terms of their cosines. Similar to unstandardised Pearson's correlation.

* **Version 0.7.2** Added weighted `wgmedian` and `wsortedeccs` to VecVecu8. Created new source file for VecVecu8 trait.

* **Version 0.7.1** Ported the improved gmedian also to VecVecu8.

* **Version 0.7.0** Made gmedian slightly more accurate. Added Weighted Geometric Median and supporing functions. Added vecu8asvecf64 utility conversion.
