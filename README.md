# ptwRust: Parametric Time Warping in Rust

> **Experimental Port**: This project is an experiment in porting the R package [`ptw`](https://github.com/rwehrens/ptw) to Rust. It leverages Rust's type safety and performance features (including parallel execution via `rayon`) while maintaining full API compatibility with the original R package.

## Overview

This package implements Parametric Time Warping (PTW), a technique used to align signal patterns (e.g., chromatograms, spectra) by non-linearly warping the time axis. The core computational logic, including the warping algorithms and optimization routines (Nelder-Mead), has been rewritten in Rust to assess the feasibility and benefits of a hybrid R/Rust architecture.

## Features

*   **Full Feature Parity**: Supports all standard `ptw` options, including:
    *   **Warping Modes**: `global` (one warp for all samples) and `individual` (per-sample warping).
    *   **Optimization Criteria**: `RMS` (Root Mean Square) and `WCC` (Weighted Cross Correlation).
    *   **Initialization**: Support for custom initial coefficients and automatic restart strategies (`try = TRUE`).
*   **High Performance**:
    *   **Parallelization**: Multithreaded processing for multiple samples using `rayon`, automatically utilizing available CPU cores.
    *   **Numerical Stability**: Implements advanced parameter scaling to ensure robust optimization convergence even for large datasets (N > 10,000) where polynomial warping typically faces conditioning issues.
*   **Rust Backend**:
    *   Custom **Nelder-Mead Simplex Optimizer** implemented in pure Rust.
    *   Seamless R integration via `extendr`.

## Installation

### Binary Installation (Recommended for non-Rust users)

If you do not have the Rust toolchain installed, you can install the pre-compiled binary packages from the [GitHub Releases](https://github.com/sipemu/ptw/releases) page.

1.  Download the appropriate file for your OS:
    *   **Windows**: `ptwRust_*.zip`
    *   **macOS**: `ptwRust_*.tgz`
2.  Install in R:

```r
# Replace path with the actual location of the downloaded file
install.packages("~/Downloads/ptwRust_0.1.0.zip", repos = NULL)
```

### From Source (Requires Rust)

Ensure you have the [Rust toolchain](https://rustup.rs/) installed (including `cargo`).

```r
# Install directly from the repository or local source
devtools::install()
```

## Usage

The API mirrors the original `ptw` package, but the package name is `ptwRust`.

```r
library(ptwRust)

# Generate synthetic data
time <- seq(0, 100, length.out = 100)
ref <- exp(-0.5 * (time - 50)^2 / 10)       # Reference peak at 50
samp <- exp(-0.5 * (time - 60)^2 / 10)      # Sample peak at 60 (shifted)

# Align sample to reference using Global Warping and RMS
model <- ptw(ref, samp, warp.type = "global", optim.crit = "RMS")

# Inspect coefficients (Expect intercept ~10.0)
print(model$warp.coef)

# Access warped signal
plot(time, ref, type = "l", col = "black", main = "PTW Alignment")
lines(time, samp, col = "red", lty = 2)
lines(time, model$warped.sample, col = "blue", lwd = 2)
legend("topright", legend = c("Reference", "Original", "Warped"), 
       col = c("black", "red", "blue"), lty = c(1, 2, 1))
```

## Architecture

*   **`R/ptw.R`**: R wrapper function that prepares data matrices and calls the Rust backend.
*   **`src/rust/`**: The Rust crate containing the logic.
    *   `src/lib.rs`: `extendr` bindings exposing `ptw_fit` and `ptw_predict`.
    *   `src/model.rs`: `PtwModel` struct, warping logic, and WCC implementation. Uses `rayon` for parallelizing individual warping tasks.
    *   `src/optim.rs`: Generic Nelder-Mead optimizer with robust parameter handling.

## Benchmarks & Verification

Unit tests are implemented in both Rust (via `cargo test`) and R (via `testthat`) to ensure numerical correctness against the behavior of the original package.

## License

MIT
