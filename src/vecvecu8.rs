use crate::{Vecu8,VecVecu8,MutVectors,Vecf64,Indices};

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
    fn wnxnonmember(self, ws:&[f64], p:&[f64]) -> Vec<f64> {
        let mut vsum = vec![0_f64; self[0].len()];
        let mut recip = 0_f64;
        for i in 0..self.len() { 
            let dvmag = self[i].vdist(&p);
            if !dvmag.is_normal() { continue } // zero distance, safe to ignore
            let rec = ws[i]/dvmag; // ws[i] is weigth for this point self[i]
            vsum.mutvadd(&self[i].smult(rec)); // add weighted vector
            recip += rec // add separately the reciprocals    
        }
        vsum.mutsmult(1.0/recip);
        vsum
    } 

    /// Weighted geometric median, sorted eccentricities magnitudes,
    /// associated reversed cummulative probability density function of the weights.
    /// So that small ecc. magnitude means high probability of being in set represented by self.
    fn wsortedeccs(self, ws: &[f64], eps:f64) -> ( Vec<f64>,Vec<f64>,Vec<f64> ) { 
        let mut eccs = Vec::with_capacity(self.len()); 
        let gm = self.wgmedian(ws,eps);
        for v in self { // collect ecentricities magnitudes
            eccs.push(v.vdist(&gm)) 
        }
        // create sort index of the eccs
        let index = eccs.mergesort(0,self.len());
        // pick the associated points weights in the same order as the sorted eccs
        let mut weights = index.unindex(true,&ws);
        let mut sumw = 0_f64;
        // accummulate the weights 
        for i in (0..weights.len()).rev() {
            sumw += weights[i]; 
            weights[i] = sumw
        }
        // divide by the sumw to get cummulative probabilities in [0,1]
        for i in 0..weights.len() { weights[i] /= sumw }; 
        ( gm, index.unindex(true, &eccs), weights )
    }
    /// Sorted cosines magnitudes,
    /// associated cummulative probability density function in [0,1] of the weights.
    /// Needs central median. 
    /// Small cos means low probability 
    fn wsortedcos(self, medmed: &[f64], zeromed: &[f64], ws: &[f64]) -> ( Vec<f64>,Vec<f64> ) { 
        let mut coses = Vec::with_capacity(self.len());  
        for p in self { // collect coses      
            coses.push(p.vsub(&medmed).vdisim(&zeromed)); 
        }
        // create sort index of the coses
        let index = coses.mergesort(0,self.len());
        // pick the associated points weights in the same order as the sorted coses
        let mut weights = index.unindex(true,&ws);
        let mut sumw = 0_f64;
        // accummulate the weights to from cpdf
        for i in 0..weights.len() {
            sumw += weights[i]; 
            weights[i] = sumw
        }
        // divide by the sum to get cum. probabilities in [0,1]
        for i in 0..weights.len() { weights[i] /= sumw };
        ( index.unindex(true,&coses), weights )
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
    fn wgmedian(self, ws: &[f64],  eps: f64) -> Vec<f64> {
        let mut p1 = self.acentroid();     
        let mut p2 = self.wnxnonmember(ws,&p1); 
        let mut e1mag = p2.vdist(&p1);    
        loop {
            // will not use this np as does nmedian, using secant instead
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

}
