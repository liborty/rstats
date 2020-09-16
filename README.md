# Rstats - Rust Stats

![Crates.io](https://img.shields.io/crates/v/rstats?logo=rust) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/liborty/rstats?logo=github) ![GitHub last commit (branch)](https://img.shields.io/github/last-commit/liborty/rstats/HEAD?logo=github)

Rstats is particularly useful for analysis of multidimensional sets of points, with applications to Machine Learning and Data Analysis. It begins with basic statistical measures and vector algebra, which provide self-contained tools for the more interesting algorithms but can also be used in their own right.

Our treatment of multidimensional sets of points is constructed from the first principles. Thus some original concepts, unlikely to be found anywhere else, are introduced and implemented here. Going beyond one dimension, other people mostly cheat by using centroids or 'quasi medians' (1-d medians along each axis). They may be quicker to compute but they are a bad start to characterising multidimensional clouds of points. Specifically, all such 1-d measures depend on the choice of axis. Typically, such dependence has to be later removed by Principle Components Analysis or similarly laborious methods.

This is a lean minimalistic library that only depends on `anyhow` (for its error handling).
Trait RStats is carefully checked and will report all kinds of errors, such as empty input.
Trait Vectors is sometimes unchecked for speed, so some caution is advisable.

## Trait RStats

Statistical measures implemented for `&[i64]` and `&[f64]`.
All these methods operate on one vector of data and take no arguments.
For example, `s.amean()` computes the arithmetic mean of slice `s` of either type.

Included are:

* means (arithmetic, geometric and harmonic),
* standard deviations,
* linearly weighted means (useful for time dependent data analysis),
* median and quartiles.

## Trait Vectors

* Vector algebra implemented for one or two `&[f64]` slices of any length (or space dimensionality).
* Pearson's, Spearman's and Kendall's correlations.
* Relationships of one vector to a set of vectors (geometric median, eccentricity).
* Relationships between sets of multidimensional vectors.

## Trait MutVectors

Some of the methods are for memory efficiency reasons reimplemented in this trait so that they mutate `self` in place instead of creating a new Vec. They are useful in vector iterative methods. Beware that they do not return anything, so they can not be chained.

## Releases

* **Version 0.4.14** Added `mutzeromd` - transforms mutable self to zero-median form.

* **Version 0.4.13** Added `trend` between two sets of points. More comments, tests and examples.

* **Version 0.4.12** Some more utilities.

* **Version 0.4.10**  Moved unimportant helper functions out of the main module.

* **Version 0.4.9** Streamlining, introduced `distances` and `eccentricities`, with speedups up to 50%.

* **Version 0.4.8** Added generics to emphasise print vectors and items of various types.

* **Version 0.4.7** Added Spearman's Rho Correlation. Removed some spurious indirections.

* **Version 0.4.6** Made eccentricity measure continuous. Added Kendall's Tau (rank) correlation. Moved all correlations to Vectors trait. Improved readme and doc comments.

* **Version 0.4.4** Medoid now finds the Outlier as well. Improved tests.
Defined and added MOE = median of eccentricities (spread) of multivariate data.

* **Version 0.4.3** Introduced computation of `nmedian` residual errors. Improved tests and benchmarks.

* **Version 0.4.2** Added `eccentricity` measures.

* **Version 0.4.1** Tidier testing and benchmarks.

* **Version 0.4.0** Cleanup. Changed the implementation types from Vecs to slices everywhere for consistency. You may need more .as_slice() conversions here and there. Made some subsidiary functions private.
