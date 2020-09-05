# Rstats - Rust Stats
![GitHub last commit](https://img.shields.io/github/last-commit/liborty/rstats)
[![Crates.io](https://img.shields.io/crates/v/rstats)](https://docs.rs/rstats)

Rstats is a lean minimalistic library that only depends on `anyhow` (for its error handling).
Trait RStats is carefully checked and will report all kinds of errors, such as empty input.
Trait Vectors is unchecked to achieve speed, so some caution is advisable.

## Trait RStats 

has statistical methods implemented for i64 and f64 vectors (and their slices).
For example, `v.amean()` computes the arithmetic mean of vector `v` of either `Vec<i64>` or `Vec<f64>` type (or their slices).

Included are:

* means (arithmetic, geometric and harmonic), 

* standard deviations,

* linearly weighted means useful for time dependent data analysis,

* correlation and autocorrelation,

* median and quartiles.

## Trait Vectors

has basic vector algebra implemented for `&[f64]` slices.
Should you get errors when applying them to `Vec<f64>`, just convert it using `.as_slice()`.

## Trait MutVectors

Some of the above Vector methods are for efficiency reasons reimplemented here so that they mutate `self` in place instead of creating a new Vec.

## Releases
* **Version 0.3.3** Added `nmedian` as the definitive algorithm for finding n-dimensional medians; `gmedian` is now the defunct Weiszfeld's algorithm which will panic and/or infinitely loop on some data. Also added benchmarks and tidied up the tests.

* **Version 0.3.2** Added `arcentroid` = n-dimensional arithmetic mean. Added some more doc examples.

* **Version 0.3.1** Geometric Median speeded up. Added trait MutVectors.

* **Version 0.3.0** completed the Geometric Median. Removed duplicated implementations of Vector for `Vec<f64>`.
