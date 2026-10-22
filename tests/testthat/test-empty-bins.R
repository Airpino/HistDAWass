test_that("empty bins get a weight of 1e-8 at creation", {
  d <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.5, 0.5, 1))
  w <- diff(d@p)
  expect_true(all(w > 0))
  expect_equal(w[2], 1e-8, tolerance = 1e-6)
  expect_equal(d@p[1], 0)
  expect_equal(d@p[4], 1)
  # mean and std are computed on the fixed cdf
  expect_equal(d@m, meanH(d))
  expect_equal(d@s, stdH(d))
})

test_that("checkEmptyBins is idempotent and works on MatH", {
  d <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.5, 0.5, 1))
  expect_identical(checkEmptyBins(d), d)
  # an object with an empty bin created by slot assignment (as in old versions)
  d2 <- d
  d2@p <- c(0, 0.5, 0.5, 1)
  d3 <- checkEmptyBins(d2)
  expect_true(all(diff(d3@p) > 0))
  M <- checkEmptyBins(MatH(list(d2, d2), nrows = 2, ncols = 1))
  expect_true(all(sapply(M@M, function(o) all(diff(o@p) > 0))))
})

test_that("data2hist produces no empty bins", {
  set.seed(1)
  x <- c(rnorm(200), rnorm(200, 30))
  for (a in c("base", "histogram")) {
    h <- data2hist(x, algo = a)
    expect_true(all(diff(h@p) > 0))
  }
  h <- data2hist(c(1:10, 50), algo = "ManualBreaks", breaks = c(0, 10, 20, 30, 51))
  expect_true(all(diff(h@p) > 0))
})

test_that("the shipped datasets have no empty bins", {
  for (o in list(China_Seas, China_Month, Age_Pyramids_2014)) {
    expect_true(all(vapply(o@M, function(d) all(diff(d@p) > 0), logical(1))))
  }
})

test_that("invalid distributions are rejected", {
  expect_error(distributionH(c(1, 3, 2), c(0, 0.5, 1)))
  expect_error(distributionH(c(1, 2, 3), c(0, 0.6, 0.5)))
  expect_error(distributionH(c(1, 2, 3), c(0, 1)))
  expect_equal(distributionH(5)@x, c(5, 5))
})

test_that("distances are correct also for objects with empty bins", {
  # brute force: integral of squared difference of quantile functions
  qf <- function(d, u) stats::approx(d@p, d@x, xout = u, ties = "ordered")$y
  u <- (seq_len(1e5) - 0.5) / 1e5
  a <- distributionH(c(0, 1, 2, 3, 4), c(0, 0.2, 0.2, 0.7, 1))
  b <- distributionH(c(1, 2, 4, 5), c(0, 0.3, 0.9, 1))
  # an "old style" object with a true zero-weight bin
  a_old <- a
  a_old@p <- c(0, 0.2, 0.2, 0.7, 1)
  bf <- mean((qf(a_old, u) - qf(b, u))^2)
  expect_equal(WassSqDistH(a, b), bf, tolerance = 1e-4)
  expect_equal(WassSqDistH(a_old, b), bf, tolerance = 1e-4)
  X <- MatH(list(a_old, b), nrows = 2, ncols = 1)
  expect_equal(WH_MAT_DIST(X)[1, 2]^2, bf, tolerance = 1e-4)
})
