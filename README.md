# Rstats - Rust Stats
![GitHub last commit](https://img.shields.io/github/last-commit/liborty/rstats)
[![Crates.io](https://img.shields.io/crates/v/rstats)](https://docs.rs/rstats)
### Synopsis

Rstats is a lean minimalistic library that only depends on `anyhow` (for its error handling).

The functions supplied compute means (arithmetic, geometric and harmonic) of i64 and f64 vectors. 

They are implemented as trait methods. For example, `v.amean()` computes the arithmetic mean of vector `v` of type either `Vec<i64>` or `Vec<f64>` (or their slices).

Also included are:

* linearly weighted means useful for time dependent data analysis,

* correlation and autocorrelation,

* median and quartiles.

### Releases 
Version 0.2.0 marks breaking change from function calls: `amean(&v)` to method invocations: `v.amean()`

More methods will be added in future versions.