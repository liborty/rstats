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

* **Version 0.3.2** Added `arcentroid` = n-dimensional arithmetic mean. Added some more doc examples.

* **Version 0.3.1** Geometric Median speeded up. Added trait MutVectors.

* **Version 0.3.0** completed the Geometric Median. Removed duplicated implementations of Vector for `Vec<f64>`.

* **Version 0.2.5** breaking change: removed the struct NDPoints and  its trait GMedian for more simplicity and generality. N-dimensional vector algebra is now included in Vectors.

* **Version 0.2.4** completed basic vector algebra.

* **Version 0.2.3** added Medoid.

* **Version 0.2.2** added simple vector calculations via trait `Vectors`. Renamed `correlation` -> `icorrelation` and `fcorrelation` -> `correlation`. As normal expectation here is for the arguments to be f64 vectors rather than i64s.

* **Version 0.2.1** added (slower) median for f64 vectors. When your data is discretisable to i64 and speed is important, use the integer version, as that one avoids sorting. Nota bene: data with fixed number of decimal places is easily discretised just by deleting the decimal points and finally dividing the result by the relevant scaling factor.

* **Version 0.2.0** introduced a breaking change from function calls: `amean(&v)` to trait method invocations: `v.amean()`.
