use std::fmt;
/// generate random data for testing, plus some auxilliary pretty printing functions

/// Sum of linear weights
pub fn wsum(n: usize) -> f64 {
    (n * (n + 1)) as f64 / 2.
}

/// helper function for formatting error messages
pub fn emsg(file: &'static str, line: u32, msg: &'static str) -> String {
    format!("{}:{} rstats {}", file, line, msg)
}

/// Generates a vector of n vectors, each of length d, all filled with random numbers for testing.
/// It needs two seeds s1 and s2. Same seeds will produce the same random sequence.  
/// Uses local closure `rand` to generate random numbers (avoids dependencies).  
/// Random numbers are in the open interval 0..1 with uniform distribution.  
pub fn genvec(d: usize, n: usize, s1: u32, s2: u32) -> Vec<Vec<f64>> {
    if n * d < 1 {
        panic!("{}",emsg(
            file!(),
            line!(),
            "genvec given zero or wrong dimensions"
        ))
    }
    // random numbers generating closure with captured seeds
    let mut m_z = s1 as u32;
    let mut m_w = s2 as u32;
    let mut rand = || {
        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
        (((m_z << 16) & m_w) as f64 + 1.0) * 2.328306435454494e-10
    };
    let mut v: Vec<Vec<f64>> = Vec::with_capacity(n);
    for _i in 0..n {
        let mut pt = Vec::with_capacity(d);
        for _j in 0..d {
            pt.push(rand())
        }
        v.push(pt)
    } // fills the lot with random numbers
    return v;
}

/// Generates a vector of n vectors, each of length d, all filled with random numbers for testing.
/// It needs two seeds s1 and s2. Same seeds will produce the same random sequence.  
/// Uses local closure `rand` to generate random numbers (avoids dependencies).  
/// Random numbers are in the closed interval 0..255 with uniform distribution.  
pub fn genvecu8(d: usize, n: usize, s1: u32, s2: u32) -> Vec<Vec<u8>> {
    if n * d < 1 {
        panic!("{}",emsg(
            file!(),
            line!(),
            "genvecu8 given zero or wrong dimensions"
        ))
    }
    // random numbers generating closure with captured seeds
    let mut m_z = s1 as u32;
    let mut m_w = s2 as u32;
    let mut rand = || {
        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
        (256.0*(((m_z << 16) & m_w) as f32 + 1.0)*2.328306435454494e-10).floor() as u8
    };
    let mut v: Vec<Vec<u8>> = Vec::with_capacity(n);
    for _i in 0..n {
        let mut pt = Vec::with_capacity(d);
        for _j in 0..d {
            pt.push(rand())
        }
        v.push(pt)
    } // fills the lot with random numbers
    return v;
}


/// GreenIt (GI) struct facilitates printing (in green) any type
/// that has Display implemented.
pub struct GI<T: fmt::Display>(pub T);
impl<T: fmt::Display> fmt::Display for GI<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "\x1B[01;92m{}\x1B[0m", self.0.to_string())
    }
}

/// GreenVec (GV) struct facilitates printing (in green) vectors of any type
/// that has Display implemented.
pub struct GV<T: fmt::Display>(pub Vec<T>);
impl<T: fmt::Display> fmt::Display for GV<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::from("\x1B[01;92m[");
        let n = self.0.len();
        if n > 0 {
            s.push_str(&self.0[0].to_string()); // first item
            for i in 1..n {
                s.push_str(", ");
                s.push_str(&self.0[i].to_string());
            }
        }
        write!(f, "{}]\x1B[0m", s)
    }
}
