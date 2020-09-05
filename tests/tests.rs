#![allow(unused_imports)]
#![allow(dead_code)]
#[cfg(test)]

use anyhow::{Result};
use rstats::{RStats,MutVectors,Vectors,genvec};

#[test]
#[ignore]
fn fstats() -> Result<()> { 
   let v1 = vec![1_f64,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.];
   println!("\n{:?}",&v1);
   println!("The arithmetic mean is \x1B[01;92m{}\x1B[0m",v1.amean().unwrap());
   println!("Arithmetic\n\x1B[01;92m{}\x1B[0m",v1.ameanstd().unwrap());
   println!("The geometric mean is \x1B[01;92m{}\x1B[0m",v1.gmean().unwrap());
   println!("The harmonic mean is \x1B[01;92m{}\x1B[0m",v1.hmean().unwrap());
   println!("The magnitude is \x1B[01;92m{}\x1B[0m",v1.as_slice().vmag());
   let mut v2 = v1.clone();
   v2.reverse();
   println!("{:?}",v2);
   println!("The scalar product is \x1B[01;92m{}\x1B[0m",v1.as_slice().dotp(&v2));
   println!("The difference is \x1B[01;92m{:?}\x1B[0m",v1.as_slice().vsub(&v2));
   println!("The distance is \x1B[01;92m{:?}\x1B[0m",v1.as_slice().vdist(&v2));
   Ok(())
}

#[test]
#[ignore]
fn intstats() -> Result<()> { 
   let v1 = vec![1_i64,2,3,4,5,6,7,8,9,10,11,12,13,14];
   println!("\n{:?}",&v1);
   println!("The arithmetic mean is \x1B[01;92m{}\x1B[0m",&v1.amean().unwrap());
   println!("Arithmetic\n\x1B[01;92m{}\x1B[0m",&v1.ameanstd().unwrap());
   println!("The geometric mean is \x1B[01;92m{}\x1B[0m",&v1.gmean().unwrap());
   println!("The harmonic mean is \x1B[01;92m{}\x1B[0m",&v1.hmean().unwrap());
   Ok(())
}

#[test]
fn mc() -> Result<()> { 
   let pts = genvec(20,50,1,2);
   let (dist,indx) = pts.as_slice().medoid(20).unwrap();
   println!("Sum of Medoid distances:\x1B[01;92m{}\x1B[0m Index: {}",dist,indx);
   let centroid = pts.as_slice().arcentroid(20);
   println!("Sum of Centroid distances:\x1B[01;92m{}\x1B[0m",pts.as_slice().distsum(20,&centroid));
   Ok(())
}
#[test]
#[should_panic]
fn gmedian() {
   let pts:Vec<f64> = vec![0.,0.,0.,0., 1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
      -1.,0.,0.,0., 0.,-1.,0.,0., 0.,0.,-1.,0., 0.,0.,0.,-1.];
   let (ds,_gm) = pts.as_slice().gmedian(4, 1e-5).unwrap();
   println!("Sum of Median distances:  \x1B[01;92m{}\x1B[0m",ds); 
}

#[test]
fn nmedian() -> Result<()> {
   let pts:Vec<f64> = vec![0.,0.,0.,0., 1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.,
      -1.,0.,0.,0., 0.,-1.,0.,0., 0.,0.,-1.,0., 0.,0.,0.,-1.];
   let (ds,_gm) = pts.as_slice().nmedian(4, 1e-5).unwrap();
   println!("Sum of Nmedian distances:  \x1B[01;92m{}\x1B[0m",ds);
   Ok(())
}
