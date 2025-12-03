# Internal wrappers for Rust functions

ptw_fit_r <- function(ref_vec, ref_n, samp_vec, samp_n, init_coef, warp_type, optim_crit, trwdth, smooth_param, try_restart) {
  .Call("wrap__ptw_fit_r", ref_vec, ref_n, samp_vec, samp_n, init_coef, warp_type, optim_crit, trwdth, smooth_param, try_restart, PACKAGE = "ptwRust")
}

ptw_predict_r <- function(samp_vec, samp_n, coeffs_flat, coeffs_n) {
  .Call("wrap__ptw_predict_r", samp_vec, samp_n, coeffs_flat, coeffs_n, PACKAGE = "ptwRust")
}
