
#[cfg(test)]
mod tests {
    use crate::model::{PtwModel, unscale_coeffs};

    #[test]
    fn test_large_shift_repro() {
        // Reproduce the user's case
        // time <- seq(0, 100, length.out = 10000)
        // ref <- exp(-0.5 * (time - 50)^2 / 10)       # Reference peak at 50
        // samp <- exp(-0.5 * (time - 60)^2 / 10)      # Sample peak at 60 (shifted)
        
        let n = 10000;
        let mut ref_sig = vec![0.0; n];
        let mut samp_sig = vec![0.0; n];
        
        for i in 0..n {
            // map index i (0..9999) to time (0..100)
            let t = i as f64 * 100.0 / (n as f64 - 1.0);
            
            ref_sig[i] = (-0.5 * (t - 50.0).powi(2) / 10.0).exp();
            samp_sig[i] = (-0.5 * (t - 60.0).powi(2) / 10.0).exp();
        }
        
        // Shift calculation:
        // Peak at t=50 (index 5000) needs to map to t=60 (index 6000).
        // w(idx) = a0 + a1*idx.
        // w(5000) = 6000.
        // If a1=1, a0 = 1000.
        
        let mut model = PtwModel::new("global".to_string(), "RMS".to_string(), 0, 0.0);
        
        // 1. Try Linear Fit [0.0, 1.0] - Should be exact
        model.fit(&vec![ref_sig.clone()], &vec![samp_sig.clone()], Some(vec![0.0, 1.0]), false);
        let a_lin = unscale_coeffs(&model.coeffs[0], n);
        println!("Linear Coeffs: {:?}", a_lin);
        assert!((a_lin[0] - 1000.0).abs() < 50.0, "Linear fit failed. Got {}", a_lin[0]);

        // 2. Try Quadratic Fit [0.0, 1.0, 0.0] - Might degenerate but should be close
        // Reset model? No, fit overwrites.
        model.fit(&vec![ref_sig.clone()], &vec![samp_sig.clone()], Some(vec![0.0, 1.0, 0.0]), false);
        let a_quad = unscale_coeffs(&model.coeffs[0], n);
        println!("Quadratic Coeffs: {:?}", a_quad);
        
        // Check if it warps correctly at peak
        // w(5000) should be ~6000
        let w_peak = a_quad[0] + a_quad[1] * 5000.0 + a_quad[2] * 5000.0 * 5000.0;
        println!("Warped peak index: {}", w_peak);
        assert!((w_peak - 6000.0).abs() < 100.0, "Quadratic warp failed at peak");
    }
}
