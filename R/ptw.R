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
#' @param mode Warp mode (e.g., "forward" or "backward"). Currently ignored in Rust implementation.
#' @param verbose Logical. If TRUE, prints progress.
#' @param ... Additional arguments.
#' @return A list containing the warped sample, coefficients, and parameters.
#' @examples
#' \dontrun{
#' # Load example data if available
#' data(gaschrom, package = "ptw")
#' ref <- gaschrom[1,]
#' samp <- gaschrom[16,]
#' gaschrom.ptw <- ptw(ref, samp)
#' summary(gaschrom.ptw)
#' 
#' ## Not run: 
#' ## comparison between backward and forward warping
#' gaschrom.ptw <- ptw(ref, samp, init.coef = c(0, 1, 0, 0), mode = "backward")
#' summary(gaschrom.ptw)
#' gaschrom.ptw <- ptw(ref, samp, init.coef = c(-10, 1, 0, 0), mode = "forward")
#' summary(gaschrom.ptw)
#' 
#' ## #############################
#' ## many samples warped on one reference
#' ref <- gaschrom[1,]
#' samp <- gaschrom[2:16,]
#' gaschrom.ptw <-
#'   ptw(ref, samp, warp.type = "individual", verbose = TRUE,
#'       optim.crit = "WCC",  trwdth = 100, init.coef = c(0, 1, 0))
#' summary(gaschrom.ptw)
#' 
#' ## #############################
#' ## several samples on several references individually
#' ref <- gaschrom[1:8,]
#' samp <- gaschrom[9:16,]
#' gaschrom.ptw <-
#'   ptw(ref, samp, warp.type = "individual",
#'       optim.crit = "WCC", trwdth = 100, init.coef = c(0, 1, 0))
#' summary(gaschrom.ptw)
#' 
#' ## #############################
#' ## several samples on several references: one, global warping
#' gaschrom.ptw <- ptw(ref, samp, warp.type = "global",
#'                     optim.crit = "WCC", init.coef = c(0, 1, 0))
#' summary(gaschrom.ptw)
#' }
#' @export
ptw <- function(ref, samp, selected.traces = NULL, init.coef = NULL, try = FALSE,
                warp.type = c("individual", "global"),
                optim.crit = c("RMS", "WCC"),
                smooth.param = 0, trwdth = 20,
                mode = "forward", verbose = FALSE, ...) {
  
  warp.type <- match.arg(warp.type)
  optim.crit <- match.arg(optim.crit)
  
  # Ensure ref and samp are matrices for uniform handling
  if (is.vector(ref)) ref <- matrix(ref, nrow = 1)
  if (is.vector(samp)) samp <- matrix(samp, nrow = 1)
  
  n_samp <- nrow(samp)
  n_ref <- nrow(ref)
  
  # Transpose for Rust (so samples are contiguous in memory)
  ref_vec <- as.vector(t(ref))
  samp_vec <- as.vector(t(samp))
  
  # Handle init.coef
  if (!is.null(init.coef)) {
     init.coef <- as.double(as.vector(t(init.coef)))
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
  
  coeffs_flat <- res$coeffs
  warped_flat <- res$warped
  n_coeffs_per_samp <- res$n_coeffs
  
  # Reshape coeffs
  n_coeff_rows <- if (warp.type == "global") 1 else n_samp
  coeffs <- matrix(coeffs_flat, nrow = n_coeff_rows, byrow = TRUE)
  
  # Reshape warped
  warped_sample <- matrix(warped_flat, nrow = n_samp, byrow = TRUE)
  
  structure(list(
    reference = ref,
    sample = samp,
    warped.sample = warped_sample,
    warp.coef = coeffs,
    warp.fun = coeffs, # Legacy support
    crit.value = res$crit_value,
    optim.crit = optim.crit,
    warp.type = warp.type,
    mode = mode # Store mode
  ), class = "ptwRust")
}
