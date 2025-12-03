#' Parametric Time Warping
#'
#' Aligns a sample signal to a reference signal using Parametric Time Warping.
#'
#' @param ref Reference signal (numeric vector or matrix).
#' @param samp Sample signal (numeric vector or matrix).
#' @param selected.traces Optional vector of indices indicating which traces to use.
#' @param init.coef Initial coefficients (vector or matrix).
#' @param try Logical. If TRUE, tries multiple initializations.
#' @param warp.type "global" or "individual".
#' @param optim.crit "RMS" or "WCC".
#' @param trwdth Width of the triangle window for WCC (integer).
#' @param smooth.param Smoothing parameter (numeric).
#' @return A list containing the warped sample, coefficients, and parameters.
#' @export
ptw <- function(ref, samp, selected.traces = NULL, init.coef = NULL, try = FALSE,
                warp.type = c("individual", "global"),
                optim.crit = c("RMS", "WCC"),
                smooth.param = 0, trwdth = 20) {
  
  warp.type <- match.arg(warp.type)
  optim.crit <- match.arg(optim.crit)
  
  # Ensure ref and samp are matrices for uniform handling
  # If vector, convert to 1-row matrix (1 sample)
  # Standard: rows are samples, cols are time points
  if (is.vector(ref)) ref <- matrix(ref, nrow = 1)
  if (is.vector(samp)) samp <- matrix(samp, nrow = 1)
  
  n_samp <- nrow(samp)
  n_ref <- nrow(ref)
  n_time <- ncol(samp)
  
  # Transpose for Rust (so samples are contiguous in memory)
  # t() makes it n_time x n_samp. Flattening gives Sample 1, Sample 2...
  ref_vec <- as.vector(t(ref))
  samp_vec <- as.vector(t(samp))
  
  # Handle init.coef
  if (!is.null(init.coef)) {
     if (is.vector(init.coef)) {
         # If global or 1 sample, fine. If individual and vector, replicate?
         # Original ptw handles this. We pass flat vector.
     }
     # Ensure doubles
     init.coef <- as.double(as.vector(t(init.coef))) # Transpose if matrix?
  }
  
  # Call Rust
  res <- ptw_fit_r(
      as.double(ref_vec), as.integer(n_ref),
      as.double(samp_vec), as.integer(n_samp),
      if(is.null(init.coef)) NULL else init.coef,
      as.character(warp.type),
      as.character(optim.crit),
      as.integer(trwdth),
      as.double(smooth.param),
      as.logical(try)
  )
  
  # Process result
  coeffs_flat <- res$coeffs
  warped_flat <- res$warped
  n_coeffs_per_samp <- res$n_coeffs
  
  # Reshape coeffs
  # If global, 1 row. If individual, n_samp rows.
  n_coeff_rows <- if (warp.type == "global") 1 else n_samp
  coeffs <- matrix(coeffs_flat, nrow = n_coeff_rows, byrow = TRUE)
  
  # Reshape warped
  warped_sample <- matrix(warped_flat, nrow = n_samp, byrow = TRUE) # Rust returned contiguous samples
  
  structure(list(
    reference = ref,
    sample = samp,
    warped.sample = warped_sample,
    warp.coef = coeffs,
    warp.fun = coeffs, # Legacy support
    crit.value = NA, # TODO: return scores
    optim.crit = optim.crit,
    warp.type = warp.type
  ), class = "ptw")
}
