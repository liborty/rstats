README

# Rstats - Rust Stats

![Crates.io](https://img.shields.io/crates/v/rstats?logo=rust) ![GitHub last commit (branch)](https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github)  
Now forward compliant with Rust 2021 Edition!

## Introduction
Rstats is primarily about characterising multidimensional sets of points, with applications to Machine Learning and Data Analysis. It begins with statistical measures and vector algebra, which provide some basic self-contained tools for the more interesting algorithms but can also be used in their own right.

Our treatment of multidimensional sets of points is constructed from the first principles. Some original concepts, not to be found elsewhere, are introduced and implemented here. Specifically, new multidimensional median algorithm.

Going beyond one dimension, most authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians may be easy to compute but they are a poor start to stable characterisation of multidimensional data.  
*Specifically, all such 1-d measures are not invariant with respect to the choice of axis.*  
Such dependence has to be later removed by Principle Components Analysis or similar methods. In contradistinction to this, our methods based on the True Geometric Median, computed here by `nmedian`, are axis (rotation) independent from the first step.

### Terminology for sets of points in n dimensions
* `Centroid\Centre\Mean` is the (generally non member) point that minimises the sum of `squares` of distances to all member points. Thus it is susceptible to outliers. In other words, it is the n-dimensional arithmetic mean. By drawing physical analogy with gravity, is is sometimes called 'the centre of mass'. Centroid can also sometimes mean the member of the set which is the nearest to the Centre. Here we follow the common (if confusing) usage and always mean the actual Centre.

* `Quasi Median` is the point minimising sums of distances separately in each dimension (its coordinates are 1-d medians along each axis). It is a mistaken concept which we do not use here. The only good thing about it is that it is dead easy to compute.

* `Medoid` is the member of the set with the least sum of distances to all other members.

* `Outlier` is the member of the set with the greatest sum of distances to all other members.

* `Median or the true geometric median`, is the point (generally non member), which minimises the sum of distances to all other members. This is the one we want. It is much less susceptible to outliers and it is rotation independent.

### Implementation 
Rstats is a lean minimalistic library that only depends on *anyhow* (for its error handling).

The constituent parts of Rstats are Rust traits grouping together functions applicable to vectors of data of relevant types. This division is necessary because generic vectors are problematic in Rust.

### Documentation 
To see the documentation, click the link on the right. Then, to see just the skeletal comments, select a trait of interest. To see more deailed comments plus some examples, scroll to the bottom of the trait and unclick [+] to the left of the `implementations` of the trait. To see tests, consult `tests.rs`.

To run the tests, use single thread. It will be slower but will produce the results in the right order:   
`cargo test --release -- --nocapture --color always --test-threads=1` 

## Trait Stats

One dimensional statistical measures implemented for `&[i64]` and `&[f64]`. 
All these methods operate on one vector of data and take no arguments.
For example, `s.amean()` returns the arithmetic mean of slice `s` of either type. This is the only attempt at genericity.  
Trait Stats is carefully checked and will report all kinds of errors, such as empty input.

Included in this trait are:

* means (arithmetic, geometric and harmonic),
* standard deviations,
* linearly weighted means (useful for time dependent data analysis),
* median and quartiles.

## Trait Vecf64

Vector algebra implemented on one or two `&[f64]` slices of any length (dimensionality):
* Autocorrelation, Pearson's, Spearman's and Kendall's correlations.
* Finding minimum and maximum, linear transformation.

This trait is sometimes unchecked (for speed), so some caution with data is advisable.

## Trait Vecu8

* Some vector algebra as above for vectors of u8 (bytes).
* Frequency count of bytes by their values (Histogram or Probability Density Function).
* Entropy measure in units of e (using natural logarithms).

## Trait MutVectors

Some of the methods are for memory efficiency reasons reimplemented in this trait so that they mutate `self` in place, instead of creating a new Vec. They are useful in vector iterative methods. Beware that they work by side-effect and do not return anything, so they can not be chained.

## Trait VecVec

Relationships of one vector to a set of vectors (of `&[f64]` end types):
* Sums of distances, eccentricity,
* centroid, medoid, true geometric median, 
* transformation to zero (geometric) median data,
* relationship between sets of multidimensional vectors: trend.

Trait VecVec is entirely unchecked, so check your data upfront. This is the more sophisticated part of the library. The true geometric median is found iteratively.

## Trait VecVecu8

Some of the above for sets of vectors of bytes.

## Trait Index
 
* `ucorrelation`(self, v: &[usize]) -> f64; Pearson's correlation coefficient of two slices, typically containing the ranks.  
* `revindex`(self) -> Vec<usize>; method for reversing an index, e.g. given a sort index, returns ranks and vice versa.
* `unindex`(self, v:&[f64]) -> Vec<f64>; collects values from v in the order given by self index. 

The methods of this trait are implemented for `&[usize]`.

## Recent Releases
* **Version 0.6.2** 

* **Version 0.6.1** Improved documentation and tests.

* **Version 0.5.10** Added information theory measure.

* **Version 0.5.9** Added a few more u8 implementations, plus 'sortm' convenience wrapper.

* **Version 0.5.8** Improved some comments. Implemented `varc` for u8.

* **Version 0.5.7** Minor version - added `dists`.

* **Version 0.5.6** Added minimal support also for vectors of bytes (of u8 end type) for vector algebra over files and images.

* **Version 0.5.5** Introduced `revindex`, `mergerank` and `mergesort`. Made 1-d quartiles more accurate. Changed all correlations to required (unchecked) methods.