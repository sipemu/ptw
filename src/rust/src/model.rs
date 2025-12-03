use crate::optim::optimize;

#[derive(Clone, Debug)]
pub struct PtwModel {
    pub coeffs: Vec<Vec<f64>>, // Normalized coefficients (b)
    pub crit_values: Vec<f64>, // Scores
    pub warp_type: String,
    pub optim_crit: String,
    pub trwdth: usize,
    pub _smooth_param: f64,
}

impl PtwModel {
    pub fn new(warp_type: String, optim_crit: String, trwdth: usize, smooth_param: f64) -> Self {
        Self {
            coeffs: vec![],
            crit_values: vec![],
            warp_type,
            optim_crit,
            trwdth,
            _smooth_param: smooth_param,
        }
    }

    pub fn fit(
        &mut self, 
        refs: &[Vec<f64>], 
        samps: &[Vec<f64>], 
        init_coeffs: Option<Vec<f64>>,
        try_restart: bool
    ) {
        if samps.is_empty() { return; }
        let n_samp = samps.len();
        let n_ref = refs.len();
        let n_points = samps[0].len();
        
        // Default 'a' coeffs (unscaled)
        let default_coeffs_a = init_coeffs.unwrap_or_else(|| vec![0.0, 1.0, 0.0]); 
        
        // Convert to 'b' coeffs (scaled)
        let init_b = scale_coeffs(&default_coeffs_a, n_points);
        
        if self.warp_type == "global" {
            let (optimized_b, score) = self.run_optimization(&init_b, refs, samps, try_restart);
            self.coeffs = vec![optimized_b];
            self.crit_values = vec![score];
        } else {
            self.coeffs = Vec::with_capacity(n_samp);
            self.crit_values = Vec::with_capacity(n_samp);
            for i in 0..n_samp {
                let ref_sig = if n_ref == 1 { &refs[0] } else { &refs[i] };
                let samp_sig = &samps[i];
                
                let (optimized_b, score) = self.run_optimization_single(&init_b, ref_sig, samp_sig, try_restart);
                self.coeffs.push(optimized_b);
                self.crit_values.push(score);
            }
        }
    }

    fn run_optimization(
        &self, 
        init_b: &[f64], 
        refs: &[Vec<f64>], 
        samps: &[Vec<f64>],
        try_restart: bool
    ) -> (Vec<f64>, f64) {
        let objective = |c: &[f64]| -> f64 {
            let mut total_err = 0.0;
            for i in 0..samps.len() {
                let ref_sig = if refs.len() == 1 { &refs[0] } else { &refs[i] };
                let samp_sig = &samps[i];
                total_err += self.calculate_error(ref_sig, samp_sig, c);
            }
            total_err
        };

        let (mut best_coeffs, mut best_score) = optimize(init_b, objective);
        
        if try_restart {
             // Perturbations in 'b' space. Identity is [0, 1, 0].
             let perturbations = vec![
                 vec![0.0, 1.0, 0.0], 
                 vec![0.0, 1.0, 0.0, 0.0], 
             ];
             
             for p in perturbations {
                 if p.len() == init_b.len() {
                     let (c, s) = optimize(&p, objective);
                     if s < best_score {
                         best_score = s;
                         best_coeffs = c;
                     }
                 }
             }
        }
        (best_coeffs, best_score)
    }
    
    fn run_optimization_single(
        &self, 
        init_b: &[f64], 
        ref_sig: &[f64], 
        samp_sig: &[f64],
        try_restart: bool
    ) -> (Vec<f64>, f64) {
        let objective = |c: &[f64]| -> f64 {
            self.calculate_error(ref_sig, samp_sig, c)
        };
        
        let (mut best_coeffs, mut best_score) = optimize(init_b, objective);
        
        if try_restart {
             let alt_init = vec![0.0, 1.0, 0.0]; 
             if alt_init.len() == init_b.len() {
                 let (c, s) = optimize(&alt_init, objective);
                 if s < best_score {
                     best_coeffs = c;
                     best_score = s;
                 }
             }
        }
        (best_coeffs, best_score)
    }

    fn calculate_error(&self, ref_sig: &[f64], samp_sig: &[f64], coeffs_b: &[f64]) -> f64 {
        let warped = warp_signal(samp_sig, coeffs_b);
        
        if self.optim_crit == "WCC" {
             let val = wcc(ref_sig, &warped, self.trwdth);
             1.0 - val
        } else {
            let mut sum_sq = 0.0;
            let n = std::cmp::min(ref_sig.len(), warped.len());
            for i in 0..n {
                sum_sq += (ref_sig[i] - warped[i]).powi(2);
            }
            (sum_sq / n as f64).sqrt()
        }
    }

    pub fn predict(&self, samples: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let mut results = Vec::with_capacity(samples.len());
        for (i, samp) in samples.iter().enumerate() {
            let coeffs = if self.warp_type == "global" {
                &self.coeffs[0]
            } else {
                &self.coeffs[i]
            };
            results.push(warp_signal(samp, coeffs));
        }
        results
    }
}

pub fn scale_coeffs(coeffs: &[f64], n: usize) -> Vec<f64> {
    if n <= 1 { return coeffs.to_vec(); }
    let m = (n - 1) as f64;
    coeffs.iter().enumerate().map(|(k, &a)| {
        if k == 0 {
            a / m
        } else {
            a * m.powi((k as i32) - 1)
        }
    }).collect()
}

pub fn unscale_coeffs(coeffs: &[f64], n: usize) -> Vec<f64> {
    if n <= 1 { return coeffs.to_vec(); }
    let m = (n - 1) as f64;
    coeffs.iter().enumerate().map(|(k, &b)| {
        if k == 0 {
            b * m
        } else {
            b * m.powi(1 - (k as i32))
        }
    }).collect()
}

pub fn warp_signal(signal: &[f64], coeffs: &[f64]) -> Vec<f64> {
    let n = signal.len();
    let w_indices = warp_time(n, coeffs);
    let mut warped = Vec::with_capacity(n);
    
    for &w in &w_indices {
        warped.push(interpolate(signal, w));
    }
    warped
}

// Uses normalized time t in [0, 1]
fn warp_time(n: usize, coeffs: &[f64]) -> Vec<f64> {
    let mut w = Vec::with_capacity(n);
    let m = if n > 1 { (n - 1) as f64 } else { 1.0 };
    for i in 0..n {
        let t = i as f64 / m;
        let mut val = 0.0;
        for (p, c) in coeffs.iter().enumerate() {
            val += c * t.powi(p as i32);
        }
        // Map back to index
        w.push(val * m);
    }
    w
}

fn interpolate(signal: &[f64], index: f64) -> f64 {
    if index <= 0.0 {
        return signal[0];
    }
    if index >= (signal.len() - 1) as f64 {
        return signal[signal.len() - 1];
    }
    
    let i = index.floor() as usize;
    let frac = index - i as f64;
    
    let y0 = signal[i];
    let y1 = signal[i + 1];
    
    y0 + frac * (y1 - y0)
}

pub fn wcc(ref_sig: &[f64], samp_sig: &[f64], width: usize) -> f64 {
    if width == 0 {
        return pearson_corr(ref_sig, samp_sig);
    }
    
    let ref_smooth = triangle_smooth(ref_sig, width);
    let samp_smooth = triangle_smooth(samp_sig, width);
    
    pearson_corr(&ref_smooth, &samp_smooth)
}

fn triangle_smooth(sig: &[f64], width: usize) -> Vec<f64> {
    let n = sig.len();
    let mut smoothed = vec![0.0; n];
    let half_w = width as isize;
    
    for i in 0..n {
        let mut sum = 0.0;
        let mut w_sum = 0.0;
        for k in -half_w..=half_w {
            let idx = i as isize + k;
            if idx >= 0 && idx < n as isize {
                let weight = 1.0 - (k.abs() as f64 / (half_w as f64 + 1.0));
                if weight > 0.0 {
                    sum += sig[idx as usize] * weight;
                    w_sum += weight;
                }
            }
        }
        smoothed[i] = if w_sum > 0.0 { sum / w_sum } else { sig[i] };
    }
    smoothed
}

fn pearson_corr(x: &[f64], y: &[f64]) -> f64 {
    let n = std::cmp::min(x.len(), y.len());
    if n == 0 { return 0.0; }
    
    let mut sum_xy = 0.0;
    let mut sum_sq_x = 0.0;
    let mut sum_sq_y = 0.0;
    
    for i in 0..n {
        sum_xy += x[i] * y[i];
        sum_sq_x += x[i] * x[i];
        sum_sq_y += y[i] * y[i];
    }
    
    let denom = (sum_sq_x * sum_sq_y).sqrt();
    if denom == 0.0 { 0.0 } else { sum_xy / denom }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wcc_basic() {
        let s1 = vec![1.0, 2.0, 3.0, 2.0, 1.0];
        let s2 = vec![1.0, 2.0, 3.0, 2.0, 1.0];
        assert!((wcc(&s1, &s2, 2) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_global_fit() {
        let mut ref_sig = vec![0.0; 100];
        for i in 0..100 {
            let x = i as f64;
            ref_sig[i] = (-0.5 * (x - 50.0).powi(2) / 10.0).exp();
        }
        
        let mut samp_sig_small = vec![0.0; 100];
        for i in 0..100 {
            let x = i as f64;
            samp_sig_small[i] = (-0.5 * (x - 52.0).powi(2) / 10.0).exp();
        }
        let samps_small = vec![samp_sig_small];
        let refs = vec![ref_sig];
        
        let mut model = PtwModel::new("global".to_string(), "RMS".to_string(), 0, 0.0);
        
        model.fit(&refs, &samps_small, Some(vec![0.0, 1.0, 0.0]), false);
        
        let a = unscale_coeffs(&model.coeffs[0], 100);
        println!("Coeffs (unscaled): {:?}", a);
        assert!((a[0] - 2.0).abs() < 0.5);
    }

    #[test]
    fn test_individual_fit() {
         let mut ref_sig = vec![0.0; 50];
        for i in 0..50 {
            let x = i as f64;
            ref_sig[i] = (-0.5 * (x - 25.0).powi(2) / 5.0).exp();
        }
        
        let mut s1 = vec![0.0; 50];
        for i in 0..50 {
            let x = i as f64;
            s1[i] = (-0.5 * (x - 27.0).powi(2) / 5.0).exp();
        }
        
        let samps = vec![s1];
        let refs = vec![ref_sig];
        
        let mut model = PtwModel::new("individual".to_string(), "RMS".to_string(), 0, 0.0);
        model.fit(&refs, &samps, Some(vec![0.0, 1.0, 0.0]), false);
        
        let a = unscale_coeffs(&model.coeffs[0], 50);
        assert!((a[0] - 2.0).abs() < 0.5);
    }
}
