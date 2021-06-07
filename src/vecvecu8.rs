use crate::{Vecu8,VecVecu8,MutVectors,Vecf64};

impl VecVecu8 for &[Vec<u8>] {
    
    fn acentroid(self) -> Vec<f64> {
    let mut centre = vec![0_f64; self[0].len()];
    for v in self {
        centre.mutvaddu8(&v)
    }
    centre.mutsmult(1.0 / (self.len() as f64));
    centre
    }

    /// Eccentricity vector added to a non member point,
    /// while the true geometric median is as yet unknown. 
    /// This function is suitable for a single non-member point. 
    fn nxnonmember(self, p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for x in self { 
            let dvmag = x.vdist(p);
            if !dvmag.is_normal() { continue } // zero distance, safe to ignore
            let rec = 1.0/dvmag;
            vsum.mutvadd(&x.smult(rec)); // add vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmult(1.0/recip);
        vsum
    }

    /// Next approximate weighted median, from a non member point. 
    fn wnxnonmember(self, ws:&[u8], p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for i in 0..self.len() { 
            let dvmag = self[i].vdist(&p);
            if !dvmag.is_normal() { continue } // zero distance, safe to ignore
            let rec = ws[i] as f64/dvmag; // ws[i] is weigth for this point self[i]
            vsum.mutvadd(&self[i].smult(rec)); // add weighted vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmult(1.0/recip);
        vsum
    } 

    /// Secant method with recovery from divergence
    /// for finding the geometric median
    fn gmedian(self, eps: f64) -> Vec<f64> {
        let mut p1 = self.acentroid();     
        let mut p2 = self.nxnonmember(&p1); 
        let mut e1mag = p2.vdist(&p1);    
        loop {
            // will not use this np as does nmedian, using secant instead
            let mut np = self.nxnonmember(&p2); 
            let e2 = np.vsub(&p2); // new vetor error, or eccentricity
            let e2mag = e2.vmag(); 
            if e2mag < eps  { return np };  
            if e1mag > e2mag {  // eccentricity magnitude decreased, good, employ secant
                np = p2.vadd(&e2.smult(p1.vsub(&p2).vmag()/(e1mag-e2mag)));                   
            }
            else { // recovery: probably overshot the minimum, shorten the jump 
                   // e2 will already be pointing moreless back
                np = p2.vadd(&e2.smult(p1.vsub(&p2).vmag()/(e1mag+e2mag)));                    
            } 
            p1 = p2;        
            p2 = np;  
            e1mag = e2mag             
        }       
    }

    /// Secant method with recovery
    /// for finding the weighted geometric median
    fn wgmedian(self, ws: &[u8],  eps: f64) -> Vec<f64> {
        let mut p1 = self.acentroid();     
        let mut p2 = self.wnxnonmember(ws,&p1); 
        let mut e1mag = p2.vdist(&p1);    
        loop {
            // will not use this np directly as does nmedian, using secant instead
            let mut np = self.wnxnonmember(ws,&p2); 
            let e2 = np.vsub(&p2); // new vector error, or eccentricity
            let e2mag = e2.vmag(); 
            if e2mag < eps  { return np };  
            if e1mag > e2mag {  // eccentricity magnitude decreased, good, employ secant
                np = p2.vadd(&e2.smult(p1.vsub(&p2).vmag()/(e1mag-e2mag)));                   
            }
            else { // recovery: probably overshot the minimum, shorten the jump 
                   // e2 will already be pointing moreless back
                np = p2.vadd(&e2.smult(p1.vsub(&p2).vmag()/(e1mag+e2mag)));                    
            } 
            p1 = p2;        
            p2 = np;  
            e1mag = e2mag             
        }       
    }
    
    /// Flattened lower triangular part of a covariance matrix for u8 vectors in self.
    /// Since covariance matrix is symmetric (positive semi definite), 
    /// the upper triangular part can be trivially generated for all j>i by: c(j,i) = c(i,j).
    /// N.b. the indexing is always assumed to be in this order: row,column.
    /// The items of the resulting lower triangular array c[i][j] are here flattened
    /// into a single vector in this double loop order: left to right, top to bottom 
    fn covar(self, m:&[f64]) -> Vec<f64> {
        let n = self[0].len(); // dimension of the vector(s)
        let mut cov:Vec<f64> = vec![0_f64; (n+1)*n/2]; // flat lower triangular results array
  
        for thisp in self { // adding up covars for all the points
            let mut covsub = 0_usize; // subscript into the flattened array cov
            let vm = thisp.vsub(&m);  // zero mean vector
            for i in 0..n {
                let thisc = vm[i]; // ith component
                // its products up to and including the diagonal (itself)
                for j in 0..i+1 { 
                    cov[covsub] += thisc*vm[j];
                    covsub += 1
                }
            }
        }
        // now compute the means and return
        cov.mutsmult(1.0_f64/self.len()as f64);
        cov
    }  

}
