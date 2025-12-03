use extendr_api::prelude::*;

mod model;
mod optim;

use model::PtwModel;

#[extendr]
fn ptw_fit_r(
    ref_vec: Vec<f64>, 
    ref_n_samples: i32,
    samp_vec: Vec<f64>,
    samp_n_samples: i32,
    init_coef: Nullable<Vec<f64>>,
    warp_type: String,
    optim_crit: String,
    trwdth: i32,
    smooth_param: f64,
    try_restart: bool
) -> List {
    
    // Helper to reconstruct Vec<Vec<f64>> from flat vector
    let reconstruct = |flat: &[f64], n_sigs: usize| -> Vec<Vec<f64>> {
        if n_sigs == 0 { return vec![]; }
        let len = flat.len() / n_sigs;
        flat.chunks(len).map(|chk| chk.to_vec()).collect()
    };

    let refs = reconstruct(&ref_vec, ref_n_samples as usize);
    let samps = reconstruct(&samp_vec, samp_n_samples as usize);

    let init = match init_coef {
        Nullable::NotNull(v) => Some(v),
        Nullable::Null => None,
    };

    let mut model = PtwModel::new(
        warp_type, 
        optim_crit, 
        trwdth as usize, 
        smooth_param
    );
    
    model.fit(&refs, &samps, init, try_restart);

    // Flatten coefficients for return
    // If global, coeffs is 1 vector. R expects matrix (1 row).
    // If individual, coeffs is N vectors. R expects matrix (N rows).
    // We return flat vector + dimensions, let R shape it.
    // Actually, simpler to just return flat vector and let R wrap it.
    
    let coeffs_flat: Vec<f64> = model.coeffs.iter().flatten().cloned().collect();
    
    // Generate warped samples using predict logic (needed for return value)
    let warped_samps = model.predict(&samps);
    let warped_flat: Vec<f64> = warped_samps.iter().flatten().cloned().collect();

    list!(
        coeffs = coeffs_flat,
        warped = warped_flat,
        n_coeffs = model.coeffs[0].len() as i32
    )
}

#[extendr]
fn ptw_predict_r(
    samp_vec: Vec<f64>, 
    samp_n_samples: i32,
    coeffs_flat: Vec<f64>, 
    coeffs_n_samples: i32
) -> Vec<f64> {
    // Reconstruct
    let samps = if samp_n_samples > 0 {
        let len = samp_vec.len() / samp_n_samples as usize;
        samp_vec.chunks(len).map(|c| c.to_vec()).collect()
    } else {
        vec![]
    };
    
    let coeffs = if coeffs_n_samples > 0 {
        let len = coeffs_flat.len() / coeffs_n_samples as usize;
        coeffs_flat.chunks(len).map(|c| c.to_vec()).collect()
    } else {
        vec![]
    };

    // Temporary model for prediction
    let mut model = PtwModel::new("".to_string(), "".to_string(), 0, 0.0);
    model.coeffs = coeffs;
    
    // Decide warp type based on coeffs vs samps
    // If coeffs has 1 set, global. If N sets, individual (assuming N matches samps)
    model.warp_type = if model.coeffs.len() == 1 { 
        "global".to_string() 
    } else { 
        "individual".to_string() 
    };

    let predicted = model.predict(&samps);
    predicted.iter().flatten().cloned().collect()
}

extendr_module! {
    mod ptw;
    fn ptw_fit_r;
    fn ptw_predict_r;
}
