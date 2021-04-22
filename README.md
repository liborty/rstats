README

# Rstats - Rust Stats

![Crates.io](https://img.shields.io/crates/v/rstats?logo=rust) ![GitHub last commit (branch)](https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github)  
22/04/2021 23:46

Rstats is primarily about characterising multidimensional sets of points, with applications to Machine Learning and Data Analysis. It begins with statistical measures and vector algebra, which provide some basic self-contained tools for the more interesting algorithms but can also be used in their own right.

Our treatment of multidimensional sets of points is constructed from the first principles. Thus some original concepts, unlikely to be found elsewhere, are introduced and implemented here.

Going beyond one dimension, most authors  'cheat' by using *quasi medians* (1-d medians along each axis). Quasi medians are quick to compute but they are a poor start to characterising multidimensional clouds of points reliably. 
*Specifically, all such 1-d measures depend on the choice of axis.*
Such dependence has to be later removed by Principle Components Analysis or similar methods. In contradistinction to this, our methods based on true Geometric Median, computed here by `nmedian`, are axis (rotation) independent.

Rstats is a lean minimalistic library that only depends on `anyhow` (for its error handling).

The constituent parts of Rstats are Rust traits grouping together various functions applicable to vectors of data of relevant types: 

## Trait Stats

One dimensional statistical measures implemented for `&[i64]` and `&[f64]`.
All these methods operate on one vector of data and take no arguments.
For example, `s.amean()` returns the arithmetic mean of slice `s` of either type.
Trait Stats is carefully checked and will report all kinds of errors, such as empty input.

Included in this trait are:

* means (arithmetic, geometric and harmonic),
* standard deviations,
* linearly weighted means (useful for time dependent data analysis),
* median and quartiles.

## Trait Vecf64

Vector algebra implemented on one or two `&[f64]` slices of any length (vector dimensionality):
* Autocorrelation, Pearson's, Spearman's and Kendall's correlations.
* Finding minimum and maximum, linear transformation.

This trait is sometimes unchecked for speed, so some caution with data is advisable.

## Trait Vecu8

Some but not all vector algebra as above, for vectors of u8.

## Trait MutVectors

Some of the methods are for memory efficiency reasons reimplemented in this trait so that they mutate `self` in place, instead of creating a new Vec. They are useful in vector iterative methods. Beware that they work by side-effect and do not return anything, so they can not be chained.

## Trait VecVec

Relationships of one vector to a set of vectors (of `&[f64]` end types):
* Sums of distances, eccentricity,
* centroid, medoid, geometric median, 
* transformation to zero (geometric) median data,
* relationship between sets of multidimensional vectors: trend.

Trait VecVec is entirely unchecked, so check your data upfront.

## Trait VecVecu8

Some of the above for sets of vectors of u8 end type.

## Trait Index
 
* `ucorrelation`(self, v: &[usize]) -> f64; Pearson's correlation coefficient of two slices, typically containing the ranks.  
* `revindex`(self) -> Vec<usize>; method for reversing an index, e.g. given a sort index, returns ranks.
* `unindex`(self, v:&[f64]) -> Vec<f64>; collects values from v in the order given by self index. 

The methods of this trait are implemented for `&[usize]`.

## Recent Releases

* **Version 0.5.10** Added information theory measure.

* **Version 0.5.9** Added a few more u8 implementations, plus 'sortm' convenience wrapper.

* **Version 0.5.8** Improved some comments. Implemented `varc` for u8.

* **Version 0.5.7** Minor version - added `dists`.

* **Version 0.5.6** Added minimal support also for vectors of bytes (of u8 end type) for vector algebra over files and images.

* **Version 0.5.5** Introduced `revindex`, `mergerank` and `mergesort`. Made 1-d quartiles more accurate. Changed all correlations to required (unchecked) methods.

* **Version 0.5.4** Added `irank,varc,kazutsugi`.

* **Version 0.5.3** Added `varea` =  magnitude of the cross product. Changed status of some methods from 'required' to 'provided'.

* **Version 0.5.2** Renamed trait RStats to Stats, to avoid naming confusion. Separated MutVecs implementations to their own module `mutvecimpls.rs`. Added some more tests. Expanded `moe` to include mean and std of eccentricities.

* **Version 0.5.1** Added scalar addition `sadd` and linear transformation `lintrans` to `Vectors`.

* **Version 0.5.0** Introduces *VecVec* trait for all multi-point methods, now implemented for type `&[Vec<f64>]`. This is a breaking change but it did allow streamlining of the code and a clean separation of the traits. Main benefit to the user is in no longer having to explicitly pass around the dimensionality d.
