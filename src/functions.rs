use crate::here;

/// Sum of linear weights 1..n.
/// Also the size of an upper or lower triangular array (including the diagonal)
pub fn wsum(n: usize) -> f64 {
    (n * (n + 1)) as f64 / 2.
}

/// Generates a vector of n vectors, each of length d, all filled with random numbers for testing.
/// It needs two seeds s1 and s2. Same seeds will produce the same random sequence.  
/// Uses local closure `rand` to generate random numbers (avoids dependencies).  
/// Random numbers are in the open interval 0..1 with uniform distribution.  
pub fn genvec(d: usize, n: usize, s1: u32, s2: u32) -> Vec<Vec<f64>> {
    if n * d < 1 { panic!("{}\n\tzero or wrong dimensions",here!()) }
    // random numbers generating closure with captured seeds m_z m_w
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
    if n * d < 1 { panic!("{}\n\tzero or wrong dimensions",here!()) } 
    // random numbers generating closure with captured seeds m_z m_w
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


