# Rstats - Rust Stats

Rstats is particularly useful for analysis of multidimensional sets of points, with applications to Machine Learning.

This is a lean minimalistic library that only depends on `anyhow` (for its error handling).
Trait RStats is carefully checked and will report all kinds of errors, such as empty input.
Trait Vectors is sometimes unchecked for speed, so some caution is advisable.

Basic statistical measures and vector algebra are implemented here for self-containment of more sophisticated algorithms.

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

Vector algebra implemented for `&[f64]` slices of any length (dimensions of space). First there are functions applying to just one and two vectors. Then there are functions expressing  the relationships of one vector to a whole set of vectors. See doc examples.

## Trait MutVectors

Some of the above basic Vector methods are for memory efficiency reasons reimplemented so that they mutate `self` in place instead of creating a new Vec. They are useful in vector iterative methods. Beware that they do not return anything, so they can not be chained.

## Releases

* **Version 0.4.12** Some more utilities

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
