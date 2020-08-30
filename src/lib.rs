mod tests;
use anyhow::{Result,Context,ensure};

/// Computed quartiles and median or time weighted versions of one column
#[derive(Default)]
pub struct Med {
    pub lquartile: f64,
    pub median: f64,
    pub uquartile: f64
}
impl std::fmt::Display for Med {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "(LQ: {}, M: {}, UQ: {})", self.lquartile, self.median, self.uquartile)
    }
}
/// Geometric mean and standard deviations ratios
#[derive(Default)]
pub struct MStats {
    pub mean: f64,
    pub std: f64
}
impl std::fmt::Display for MStats {
   fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
      write!(f, "Mean:\t{}\nStd:\t{}", self.mean, self.std)
   }
}

fn cmsg(file:&'static str, line:u32, msg:&'static str)-> String {
   format!("{}:{} stats {}",file,line,msg)
}

/// The sum of linear weights is given by this simple formula
fn wsum(n: usize) -> f64 { (n*(n+1)) as f64/2. }

/// Arithmetic mean of an i64 slice
pub fn amean(dvec: &[i64]) -> Result<f64> { 
   let n = dvec.len();
   ensure!(n > 0, "{}:{} amean - supplied sample is empty!",file!(),line!() );
   Ok( dvec.iter().sum::<i64>() as f64 / (n as f64) )
}

pub fn ameanstd(dvec: &[i64]) -> Result<MStats> {
   let n = dvec.len();
   ensure!(n > 0,"{}:{} ameanstd - supplied sample is empty!",file!(),line!());
   let mut sx2:i64 = 0;
   let mean = dvec.iter().map(|&x|{ sx2+=x*x; x}).sum::<i64>() as f64 / (n as f64);
   Ok( MStats { 
      mean : mean, 
      std : (sx2 as f64/(n as f64) - mean.powi(2)).sqrt() } )
}

/// Liearly weighted arithmetic mean of an i64 slice
/// Linearly descending weights from the number of items down to one
pub fn awmean(dvec: &[i64]) -> Result<f64> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} awmean - supplied sample is empty!",file!(),line!());
	let mut iw = dvec.len() as i64 + 1; // descending linear weights
	Ok( dvec.iter().map(|&x| { iw -= 1; iw*x }).sum::<i64>() as f64 / wsum(n))
}

pub fn awmeanstd(dvec: &[i64]) -> Result<MStats> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} awmeanstd - supplied sample is empty!",file!(),line!());
   let mut sx2 = 0f64;
   let mut iw = n as f64; // descending linear weights
   let mean = dvec.iter().map( |&x| { 
      let wx = iw*x as f64;
      sx2 += wx*x as f64;
      iw -= 1.; 
      wx } ).sum::<f64>() as f64 / wsum(n);
   Ok( MStats { 
      mean : mean, 
      std : (sx2 as f64/wsum(n) - mean.powi(2)).sqrt() } )  
}

/// Harmonic mean of an i64 slice
pub fn hmean(dvec: &[i64]) -> Result<f64> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} hmean - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   for &x in dvec {
      ensure!(x != 0i64,"{}:{} hmean does not accept zero valued data!",file!(),line!());  
      sum += 1.0/(x as f64) 
   }
   Ok ( n as f64 / sum )
}

/// Linearly weighted harmonic mean of an f64 slice
/// Linearly descending weights from the number of items down to one
pub fn hwmean(dvec: &[i64]) -> Result<f64> {
   let mut n = dvec.len();
   ensure!(n>0,"{}:{} hwmean - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   for &x in dvec {
      ensure!(x!=0i64,
         "{}:{} hwmean does not accept zero valued data!",file!(),line!());  
      sum += n as f64/x as f64;
      n -= 1; 
   }
   Ok( wsum(dvec.len()) / sum )
}

/// Geometric mean of an i64 slice.
/// The geometric mean is just an exponential of an arithmetic mean
/// of log data (natural logarithms of the data items).
/// The geometric mean is less sensitive to outliers near maximal value.
/// However, we must normalise the data to avoid taking ln of zero,
/// as minus infinity does not play nice on computers.

pub fn gmean(dvec: &[i64]) -> Result<f64> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} gmean - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   for &x in dvec {   
      ensure!(x!=0i64,
         "{}:{} gmean does not accept zero valued data!",file!(),line!()); 
      sum += (x as f64).ln()
   }
   Ok( (sum/(n as f64)).exp() )
}

pub fn gmeanstd(dvec: &[i64]) -> Result<MStats> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} gmeanstd - supplied sample is empty!",file!(),line!());
   let mut sum = 0f64;
   let mut sx2 = 0f64;
   for &x in dvec { 
      ensure!(x!=0i64,
         "{}:{} gmeanstd does not accept zero valued data!",file!(),line!());   
      let lx = (x as f64).ln();
      sum += lx;
      sx2 += lx*lx    
   }
   sum /= n as f64;
   Ok( MStats { 
      mean: sum.exp(), 
      std: (sx2/(n as f64) - sum.powi(2)).sqrt().exp() }
    )
}

/// Time linearly weighted geometric mean of an i64 slice
/// Linearly descending weights from the number of items down to one
/// The geometric mean is just an exponential of an arithmetic mean
/// of log data (natural logarithms of the data items)
/// The geometric mean is less sensitive to outliers near maximal value
pub fn gwmean(dvec: &[i64]) -> Result<f64> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} gwmean - supplied sample is empty!",file!(),line!());
   let mut iw = n as i64; // descending weights
   let mut sum = 0f64;
   for &x in dvec {  
      ensure!(x!=0i64,
         "{}:{} gwmean does not accept zero valued data!",file!(),line!()); 
      sum += (iw as f64)*(x as f64).ln();
      iw -= 1;
   }
   Ok( (sum/wsum(n)).exp() )
}	

pub fn gwmeanstd(dvec: &[i64]) -> Result<MStats> {
   let n = dvec.len();
   ensure!(n>0,"{}:{} gwmeanstd - supplied sample is empty!",file!(),line!());
   let mut iw = n as i64; // descending weights
   let mut sum = 0f64;
   let mut sx2 = 0f64;
   for &x in dvec { 
      ensure!(x!=0i64,
         "{}:{} gwmeanstd does not accept zero valued data!",file!(),line!());  
      let lx = (x as f64).ln();
      sum += (iw as f64)*lx;
      sx2 += (iw as f64)*lx*lx;
      iw -= 1;
   }
   sum /= wsum(n);
   Ok( MStats { 
      mean : sum.exp(),
      std : (sx2 as f64/wsum(n) - sum.powi(2)).sqrt().exp() }
    )
}	

/// Fast median (avoids sorting)
/// Works on  data discretised into ints [u32] with minimum and maximum [low, high].
/// high-low must be a moderate range, not exceeding u16size (= 65535)
pub fn median(data: &[i64]) -> Result<Med> {
   let max = *data.iter().max().with_context(||cmsg(file!(),line!(),"median failed to find maximum"))?;
   let min = *data.iter().min().with_context(||cmsg(file!(),line!(),"median failed to find minimum"))?;
   let range =  (max-min+1) as usize;
   ensure!(range <= u16::max_value() as usize, // range too big to use as subscripts
      "{}:{} median range {} of values exceeds u16",file!(),line!(),range);
	let mut acc = vec![0_u16; range]; // min max values inclusive
   for &item in data { acc[(item-min) as usize] += 1_u16 } // frequency distribution
   let mut result: Med = Default::default();
   let rowlength = data.len();
   let mut cumm = 0_usize;
   let mut i2;

   for i in 0..range {
      // find the lower quartile
      cumm += (acc[i]) as usize; // accummulate frequencies
      if 4 * cumm >= rowlength {
         result.lquartile = (i as i64 + min) as f64; // restore min value
         break;
      }
   }
   cumm = 0usize;
   for i in (0..range).rev() {
      // find the upper quartile
      cumm += (acc[i]) as usize; // accummulate frequencies
      if 4 * cumm >= rowlength {
         result.uquartile = (i as i64 + min) as f64;
         break;
      }
   }
   cumm = 0usize;
   for i in 0..range {
   // find the midpoint of the frequency distribution
      cumm += (acc[i]) as usize; // accummulate frequencies
      if 2 * cumm == rowlength {
         // even, the other half must have the same value
         i2 = i + 1;
         while acc[i2] == 0 { i2 += 1 }
         // first next non-zero acc[i2] must represent the other half
         result.median = ((i + i2) as i64 + 2*min) as f64 / 2.;
         break;
      }
      if 2 * cumm > rowlength {
         result.median = (i as i64 + min) as f64;
         break;
      }
      // first over the half, this must be the odd midpoint
   }
   Ok(result)
}
/*
/// Fast weighted median (avoids sorting)
/// First discretise data into ints [u32] with minimum and maximum [low, high].
/// high - low must be a moderate range, not exceeding u16size (= 65535).
pub fn wmedian(max: usize, min: usize, data: &[u32]) -> Result<Med> {
    let range = max - min + 1usize; // range too big to use as subscripts
    if range > u16::max_value() as usize {
        bail!("{}:{} wmedian range {} of values exceeds u16",file!(),line!(),range);
    }
    let mut acc = vec![0u16; range]; // create subscripts vector
    let rowlength = data.len();
    for i in 0..rowlength {
        // construct weighted frequency distribution
        acc[data[i] as usize - min] += (rowlength - i) as u16;
    }
    let mut result: Med = Default::default();
    let mut cumm = 0usize;
    let mut i2;
    // median will be at half of the sum of linear weights 1 to rowlength
    let sumw = rowlength * (rowlength + 1);
    // linearly (time) weighted observations, older get lower weights
    // weights start at the number of observations for the most recent and
    // end with 1 for the oldest. Most suitable for time series data.
    for i in 0..range {
        // find the lower quartile
        cumm += (acc[i]) as usize; // accummulate frequencies
        if 8 * cumm >= sumw {
            result.lquartile = 10usize * (i + min);
            break;
        }
    }
    cumm = 0usize;
    for i in (0..range).rev() {
        // find the upper quartile
        cumm += (acc[i]) as usize; // accummulate frequencies
        if 8 * cumm >= sumw {
            result.uquartile = 10usize * (i + min);
            break;
        }
    }
    cumm = 0usize;
    for i in 0..range {
        // find the midpoint of the frequency distribution
        cumm += (acc[i]) as usize; // accummulate frequencies
        if 4 * cumm == sumw {
            // even, the other half must have the same value
            i2 = i + 1;
            while acc[i2] == 0 {
                i2 += 1;
            }
            // first next non-zero acc[i2] must represent the other half
            result.median = 5usize * (i + i2 + 2 * min);
            break;
        }
        if 4 * cumm > sumw {
            result.median = 10usize * (i + min);
            break;
        }
    }
    Ok(result) // to get here would be an error
}
*/
/// Correlation coefficient of a sample of two integer variables
pub fn correlation(v1:&[i64],v2:&[i64]) -> Result<f64> {
   let n = v1.len();
   ensure!(n>0,cmsg(file!(),line!(),"correlation - first sample is empty"));
   ensure!(n==v2.len(),cmsg(file!(),line!(),"correlation - samples are not of the same size"));
   let (mut sy,mut sxy,mut sx2,mut sy2) = (0,0,0,0);
   let sx:i64 = v1.iter().enumerate().map(|(i,&x)| {
      let y = v2[i]; sy += y; sxy += x*y; sx2 += x*x; sy2 += y*y; x    
      }).sum();
   let (sxf,syf,sxyf,sx2f,sy2f,nf) = 
       (sx as f64,sy as f64,sxy as f64,sx2 as f64,sy2 as f64,n as f64);
   Ok( (sxyf-sxf/nf*syf)/(((sx2f-sxf/nf*sxf)*(sy2f-syf/nf*syf)).sqrt()) )
}

/// Correlation coefficient of pairs of successive values of time series integer variable
pub fn autocorr(v1:&[i64]) -> Result<f64> {
   let n = v1.len();
   ensure!(n>=2,cmsg(file!(),line!(),"autocorr - sample is too small"));
   let (mut sx,mut sy,mut sxy,mut sx2,mut sy2) = (0,0,0,0,0);
   for i in 0..n-1 {
      let x = v1[i]; let y = v1[i+1]; 
      sx += x; sy += y; sxy += x*y; sx2 += x*x; sy2 += y*y }    
   let (sxf,syf,sxyf,sx2f,sy2f,nf) = 
       (sx as f64,sy as f64,sxy as f64,sx2 as f64,sy2 as f64,n as f64);
   Ok( (sxyf-sxf/nf*syf)/(((sx2f-sxf/nf*sxf)*(sy2f-syf/nf*syf)).sqrt()) )
}
