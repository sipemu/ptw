# ptw: Parametric Time Warping in Rust

> **Experimental Port**: This project is an experiment in porting the R package [`ptw`](https://github.com/rwehrens/ptw) to Rust, leveraging the `trueno` numeric crate and `extendr` for R bindings. The goal is to explore performance improvements and type safety in scientific computing workflows.

## Overview

This package implements Parametric Time Warping (PTW), a technique used to align signal patterns (e.g., chromatograms, spectra) by non-linearly warping the time axis. The core computational logic, including the warping algorithms and optimization routines (Nelder-Mead), has been rewritten in Rust to assess the feasibility and benefits of a hybrid R/Rust architecture.

## Features

*   **Full Feature Parity**: Supports all standard `ptw` options, including:
    *   **Warping Modes**: `global` (one warp for all samples) and `individual` (per-sample warping).
    *   **Optimization Criteria**: `RMS` (Root Mean Square) and `WCC` (Weighted Cross Correlation).
    *   **Initialization**: Support for custom initial coefficients and automatic restart strategies (`try = TRUE`).
*   **Rust Backend**:
    *   Custom **Nelder-Mead Simplex Optimizer** implemented in pure Rust.
    *   High-performance signal interpolation and array manipulation.
    *   Seamless R integration via `extendr`.

## Installation

Ensure you have the [Rust toolchain](https://rustup.rs/) installed (including `cargo`).

### From Source

```r
# Install directly from the repository or local source
devtools::install()
```

## Usage

The API mirrors the original `ptw` package.

```r
library(ptw)

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
    *   `src/model.rs`: `PtwModel` struct, warping logic, and WCC implementation.
    *   `src/optim.rs`: Generic Nelder-Mead optimizer.

## Benchmarks & Verification

Unit tests are implemented in both Rust (via `cargo test`) and R (via `testthat`) to ensure numerical correctness against the behavior of the original package.

## License

MIT

