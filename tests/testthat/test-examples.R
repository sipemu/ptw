
test_that("Original package examples (Mock Data)", {
  # Since we might not have 'gaschrom' data without original 'ptw',
  # we create mock data simulating 'gaschrom' structure.
  
  # Simulate chromatograms (Gaussian peaks)
  time <- seq(0, 100, length.out = 100)
  
  # Ref: 2 peaks
  ref <- exp(-0.5 * (time - 30)^2 / 2) + exp(-0.5 * (time - 70)^2 / 2)
  
  # Samp: 2 peaks shifted by +2
  samp <- exp(-0.5 * (time - 32)^2 / 2) + exp(-0.5 * (time - 72)^2 / 2)
  
  # Example 1: Simple ptw
  # gaschrom.ptw <- ptw(ref, samp)
  # summary(gaschrom.ptw)
  gaschrom.ptw <- ptw(ref, samp, init.coef = c(0, 1, 0))
  expect_s3_class(gaschrom.ptw, "ptwRust")
  expect_equal(nrow(gaschrom.ptw$warp.coef), 1)
  
  # Example 2: Backward/Forward warping
  # gaschrom.ptw <- ptw(ref, samp, init.coef = c(0, 1, 0, 0), mode = "backward")
  # Note: we added 'mode' to function signature but implementation might ignore it currently.
  # Just verifying API compatibility.
  ptw_back <- ptw(ref, samp, init.coef = c(0, 1, 0), mode = "backward")
  expect_s3_class(ptw_back, "ptwRust")
  expect_equal(ptw_back$mode, "backward")
  
  ptw_fwd <- ptw(ref, samp, init.coef = c(-2, 1, 0), mode = "forward")
  expect_s3_class(ptw_fwd, "ptwRust")
  
  # Example 3: Many samples warped on one reference
  # ref <- gaschrom[1,]
  # samp <- gaschrom[2:16,]
  # gaschrom.ptw <- ptw(ref, samp, warp.type = "individual", verbose = TRUE, optim.crit = "WCC", trwdth = 100, init.coef = c(0, 1, 0))
  
  # Mock: 3 samples
  samp_mat <- rbind(samp, samp, samp)
  ptw_ind <- ptw(ref, samp_mat, warp.type = "individual", optim.crit = "RMS", trwdth = 20, init.coef = c(0, 1, 0))
  expect_equal(nrow(ptw_ind$warp.coef), 3)
  
  # Example 4: Global warping
  # gaschrom.ptw <- ptw(ref, samp, warp.type = "global", optim.crit = "WCC", init.coef = c(0, 1, 0))
  ptw_glob <- ptw(ref, samp_mat, warp.type = "global", optim.crit = "RMS", init.coef = c(0, 1, 0))
  expect_equal(nrow(ptw_glob$warp.coef), 1)
})

test_that("Original package examples (Real Data if available)", {
  skip_if_not_installed("ptw")
  
  # If 'ptw' is installed, we can load its data
  # We need to be careful not to mask our 'ptw' function.
  # data() loads into global env usually.
  
  # Load data in a safe way
  e <- new.env()
  data("gaschrom", package = "ptw", envir = e)
  
  ref <- e$gaschrom[1,]
  samp <- e$gaschrom[16,]
  
  # Run ptw from OUR package
  res <- ptw(ref, samp, init.coef = c(0, 1, 0))
  expect_s3_class(res, "ptwRust")
  
  # Run individual on subset
  samp_set <- e$gaschrom[2:5,]
  res_ind <- ptw(ref, samp_set, warp.type = "individual", init.coef = c(0, 1, 0))
  expect_equal(nrow(res_ind$warp.coef), 4)
  
  # Run global
  res_glob <- ptw(ref, samp_set, warp.type = "global", init.coef = c(0, 1, 0))
  expect_equal(nrow(res_glob$warp.coef), 1)
})


