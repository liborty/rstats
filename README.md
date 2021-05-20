README

# Rstats - Rust Stats

![Crates.io](https://img.shields.io/crates/v/rstats?logo=rust) ![GitHub last commit (branch)](https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github)  
Now forward compliant with Rust 2021 Edition!

## Usage

Insert into your Cargo.toml file [dependencies] section:

```rust
rstats = "^0" 
```

## Introduction

Rstats is primarily about characterising multidimensional sets of points, with applications to Machine Learning and Data Analysis. It begins with statistical measures and vector algebra, which provide some basic self-contained tools for the more interesting algorithms but can also be used in their own right. Other general tools included are efficient ranking, sorting and searching. 

Our treatment of multidimensional sets of points is constructed from the first principles. Some original concepts, not to be found elsewhere, are introduced and implemented here. Specifically, the new multidimensional (geometric) median algorithm.

Going beyond one dimension, most authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians may be easy to compute but they are a poor start to stable characterisation of multidimensional data.

*Specifically, all such 1-d measures are sensitive to the choice of axis.* 

Such dependence has to be later removed by Principle Components Analysis or similar methods. In contradistinction to this, our methods based on the True Geometric Median, computed here by `gmedian`, are axis (rotation) independent from the first step.

### Terminology for sets of points in n dimensions

* `Centroid\Centre\Mean` is the (generally non member) point that minimises the sum of *squares* of distances to all member points. Thus it is susceptible to outliers. In other words, it is the n-dimensional arithmetic mean. By drawing physical analogy with gravity, it is sometimes called 'the centre of mass'. Centroid can also sometimes mean the member of the set which is the nearest to the Centre. Here we follow the common (if confusing) usage and always mean the actual Centre.

* `Quasi Median` is the point minimising sums of distances separately in each dimension (its coordinates are 1-d medians along each axis). It is a mistaken concept which we do not use here. The only good thing about it is that it is easy to compute.

* `Medoid` is the member of the set with the least sum of distances to all other members.

* `Outlier` is the member of the set with the greatest sum of distances to all other members.

* `Median or the true geometric median (gm)`, is the point (generally non member), which minimises the sum of distances to all other members. This is the one we want. It is much less susceptible to outliers and it is rotation independent.

### Implementation

Rstats is a lean minimalistic library that only depends on *anyhow* (for its simple error handling).

The constituent parts of Rstats are Rust traits grouping together functions applicable to vectors of data of relevant end types. This division is necessary because generic vectors are problematic in Rust. 

End type f64 is most commonly used. Facilities for other end types are limited. For lots of data of other end types, it is always possible to clone to f64, see for example the included utility function `vecu8asvecf64`.

### Documentation

Follow the documentation link. Then select a trait of interest to see the skeletal comments on the prototype function declarations in lib.rs. To see more detailed comments, plus some examples from the implementation files, scroll to the bottom of the trait and unclick [+] to the left of the `implementations` of the trait. To see the tests, consult `tests.rs`.

To run the tests, use single thread. It will be slower but will produce the results in the right order:

```rust
cargo test --release -- --test-threads=1 --nocapture --color always 
```

## Trait Stats

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

## Trait Vecf64

Vector algebra implemented on one or two `&[f64]` slices of any length (dimensionality):

* Autocorrelation, Pearson's, Spearman's and Kendall's correlations.
* Finding minimum and maximum, linear transformation to [0,1].
* Indirect merge sort, binary search. 

This trait is sometimes unchecked (for speed), so some caution with data is advisable.

## Trait Vecu8

* Some vector algebra as above for vectors of u8 (bytes).
* Frequency count of bytes by their values (Histogram or Probability Density Function).
* Entropy measures in units of e (using natural logarithms).

## Trait MutVectors

Some of the above functions are for memory efficiency reasons reimplemented in this trait so that they mutate `self` in place, instead of creating a new Vec. Clearly, they can only be applied to a mutable variable. They are useful in vector iterative methods. Beware that they work by side-effect and do not return anything, so they can not be chained.

## Trait VecVec

Relationships of one vector to a set of vectors (of `&[f64]` end types):

* sums of distances, eccentricity,
* centroid, medoid, true geometric median,
* transformation to zero (geometric) median data,
* relationship between sets of multidimensional vectors: trend.

Trait VecVec is entirely unchecked, so check your data upfront. This is the more sophisticated part of the library. The true geometric median is found iteratively.

## Trait VecVecu8

Some of the above for vectors of vectors of bytes.

## Trait Index

The functions of this trait are implemented for vectors of subscripts, i.e. `&[usize]`.

* `ucorrelation`(self, v: &[usize]) -> f64; Pearson's correlation coefficient of two slices, typically containing the ranks.  
* `revindex`(self) -> Vec\<usize\>; method for reversing an index, e.g. given a sort index, returns ranks and vice versa.
* `unindex`(self, v:&[f64]) -> Vec\<f64\>; collects values from v in the order given by self index.

## Recent Releases

* **Version 0.7.2** Added weighted `wgmedian` and `wsortedeccs` to VecVecu8. Created new source file for VecVecu8 trait.

* **Version 0.7.1** Ported the improved gmedian also to VecVecu8. Added to wsortedeccs outputs associated cummulative probability density function of the weights.

* **Version 0.7.0** Made gmedian slightly more accurate. Added Weighted Geometric Median and supporing functions. Added vecu8asvecf64 utility conversion.

* **Version 0.6.10** Added to vecu8 for completeness. 

* **Version 0.6.9** Added `sortedeccs` : good descriptive measure for a set of points in nD. Added `binsearch`.

* **Version 0.6.8** Added `exacteccs` for obtaining eccentricities after the gm has been found. Added some tests.

* **Version 0.6.7** Eccentricities and geometric medians optimisations. Gmedian ported also to points defined by &[u8].

* **Version 0.6.6** Deleted some defunct functions.

* **Version 0.6.5** Simplified eccentricities functions. The new geometric median algorithm `gmedian` consistently beats on time even the improved Weizsfeld `nmedian`.

* **Version 0.6.4** Added twopoint and secant GM algorithms (experimental).

* **Version 0.6.3** Fixed dependence measure interval to [0,1].

* **Version 0.6.2** Fixed entropy bug, added jointpdf, joint entropy and dependence.

* **Version 0.6.1** Improved documentation and tests.