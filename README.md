# Rstats - Rust Stats
![GitHub last commit](https://img.shields.io/github/last-commit/liborty/rstats)
[![Crates.io](https://img.shields.io/crates/v/rstats)](https://docs.rs/rstats)
## Synopsis

Rstats is a lean minimalistic library that only depends on `anyhow` (for its error handling).

The implemented trait methods compute means (arithmetic, geometric and harmonic) of i64 and f64 vectors.

For example, `v.amean()` computes the arithmetic mean of vector `v` of type either `Vec<i64>` or `Vec<f64>` (or their slices).

Also included are:

* linearly weighted means useful for time dependent data analysis,

* correlation and autocorrelation,

* median and quartiles.

## Releases

* **Version 0.2.1** added (slower) median for f64 vectors. When your data is discretisable to i64 and speed is important, use the integer version, as that one avoids sorting. Nota bene: data with fixed number of decimal places is easily discretised just by deleting the decimal points and finally dividing the result by the relevant scaling factor.

* **Version 0.2.0** introduced a breaking change from function calls: `amean(&v)` to trait method invocations: `v.amean()`.

## Future Work
More methods will be added in future versions; in particular the geometric median of n-dimensional data, using an improved Weiszfeld's algorithm.
