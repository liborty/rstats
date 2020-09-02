#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use anyhow::{Result};
use crate::{RStats,Vectors};

#[test]
fn vectors() -> Result<()> { 
   let v = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   let mut v2 = v.clone();
   v2.reverse();
   println!("\nThe magnitude of\n{:?} is \x1B[01;92m{}\x1B[0m",&v,&v.vmag());
   println!("\nThe scalar product of\n{:?}\n{:?} is \x1B[01;92m{}\x1B[0m",&v,&v,v.dotp(&v));
   println!("\nThe difference of\n{:?}\n{:?} is\n\x1B[01;92m{:?}\x1B[0m",&v,&v,v.vsub(&v));
   println!("\nThe distance between\n{:?}\n{:?} is\n\x1B[01;92m{:?}\x1B[0m",&v,&v2,v.vdist(&v2));
   Ok(())
}

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
*/  
