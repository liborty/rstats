#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use anyhow::{Result};
use crate::{RStats,correlation,autocorr,median,awmean,hmean,hwmean,ameanstd,awmeanstd,gmeanstd,gwmeanstd};
// constant testing vectors used by tests
const VEC1:[i64;14] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
const VEC2:[i64;14] = [14,13,12,11,10,9,8,7,6,5,4,3,2,1];

/*
#[test]
fn mean() -> Result<()> { 
let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   println!("\nThe arithmetic mean of\n{:?} is \x1B[01;92m{}\x1B[0m",&v1,&v1.amean().unwrap());
   Ok(())
}

#[test]
fn fmean() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
      println!("\nThe arithmetic mean of\n{:?} is \x1B[01;92m{}\x1B[0m",&v1,&v1.amean().unwrap());
      Ok(())
   }
  
#[test]
fn corr() -> Result<()> { 
   println!("\nCorrelation of\n{:?}\n{:?}\n\x1B[01;92m{}\x1B[0m",
   &VEC1,&VEC2,correlation(&VEC1,&VEC2).unwrap());
   Ok(())
}

#[test]
fn acorr() -> Result<()> { 
   println!("\nAuto correlation of\n{:?}\n\x1B[01;92m{}\x1B[0m",
   &VEC1,autocorr(&VEC1).unwrap());
   Ok(())
}

#[test]
fn med() -> Result<()> { 
   println!("\nMedian of\n{:?}\n\x1B[01;92m{}\x1B[0m",
   &VEC1,median(&VEC1).unwrap());
   Ok(())
}

#[test]
fn gmv() -> Result<()> { 
   println!("\nArithmetic mean of\n{:?}\n\x1B[01;92m{}\x1B[0m\n",&VEC1,ameanstd(&VEC1).unwrap());
   println!("\nWeighted arithmetic mean of\n{:?}\n\x1B[01;92m{}\x1B[0m\n",&VEC1,awmeanstd(&VEC1).unwrap());
   println!("\nGeometric mean of\n{:?}\n\x1B[01;92m{}\x1B[0m\n",&VEC1,gmeanstd(&VEC1).unwrap());
   println!("\nWeighted geometric mean of\n{:?}\n\x1B[01;92m{}\x1B[0m\n",&VEC1,gwmeanstd(&VEC1).unwrap());
   println!("\nHarmonic mean of\n{:?}\n\x1B[01;92m{}\x1B[0m\n",&VEC1,hmean(&VEC1).unwrap());
   println!("\nWeighted harmonic mean of\n{:?}\n\x1B[01;92m{}\x1B[0m\n",&VEC1,hwmean(&VEC1).unwrap());   
   Ok(())
}
*/