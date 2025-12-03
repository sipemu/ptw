
test_that("crit.value is returned", {
  time <- seq(0, 100, length.out = 100)
  ref <- exp(-0.5 * (time - 50)^2 / 10)
  samp <- exp(-0.5 * (time - 60)^2 / 10)
  
  # Global
  res <- ptw(ref, samp, warp.type = "global", optim.crit = "RMS")
  expect_false(any(is.na(res$crit.value)))
  expect_true(is.numeric(res$crit.value))
  # RMS should be small for perfect alignment
  expect_lt(res$crit.value[1], 0.1)
  
  # Individual
  res_ind <- ptw(ref, rbind(samp, samp), warp.type = "individual", optim.crit = "RMS")
  expect_false(any(is.na(res_ind$crit.value)))
  expect_length(res_ind$crit.value, 2)
})

