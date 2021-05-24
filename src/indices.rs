use crate::Indices;

impl Indices for &[usize] {

    /// Constructs reversed (inverted) index, eg. from sort index to data ranks
    /// This is a symmetric operation, i.e. any even number of applications 
    /// leads back to the original form.
    fn revindex(self) -> Vec<usize> {
        let n = self.len();
        let mut index = vec![0_usize;n];
        for i in 0..self.len() { index[self[i]] = i };
        index
    }

    /// Collects values from v in the order given by self index 
    /// when ascending is true, otherwise in reverse (descending) order.
    /// Used by msort for ascending or descending sort.
    /// Good for efficient sorting of any vectors by their separate keys.    
    fn unindex(self, ascending: bool, v:&[f64]) -> Vec<f64> {
        if ascending { self.iter().map(|&i| v[i]).collect() }
        else { self.iter().rev().map(|&i| v[i]).collect()   } 
    }

    /// Pearson's correlation coefficient of a two $[usize] slices,
    /// typically the ranks. In which case this is the Spearman's correlation, where the ranks
    /// have been computed previously.
    fn ucorrelation(self, v: &[usize]) -> f64 {
        let (mut sy, mut sxy, mut sx2, mut sy2) = (0_f64, 0_f64, 0_f64, 0_f64);
        let sx: f64 = self
            .iter()
            .zip(v)
            .map(|(&ux, &uy)| {
                let x = ux as f64;
                let y = uy as f64;
                sy += y;
                sxy += x * y;
                sx2 += x * x;
                sy2 += y * y;
                x
            })
            .sum();
        let nf = self.len() as f64;
        (sxy - sx / nf * sy) / ((sx2 - sx / nf * sx) * (sy2 - sy / nf * sy)).sqrt()
    }
}
