test_that("moving averages refer to the central element of the window", {
  H <- subsetHTS(RetHTS, 1, 10)
  ma <- HTS.moving.averages(H, k = 3)$AveragedHTS
  expect_length(ma@data, 8)
  expect_equal(ma@data[[1]]@tstamp, H@data[[2]]@tstamp)
  expect_equal(ma@data[[8]]@tstamp, H@data[[9]]@tstamp)
  ma5 <- HTS.moving.averages(H, k = 5)$AveragedHTS
  expect_equal(ma5@data[[1]]@tstamp, H@data[[3]]@tstamp)
  expect_error(HTS.moving.averages(H, k = 3, weights = c(1, 1)))
})

test_that("k-NN forecast handles ties and uses the k nearest neighbours", {
  d1 <- new("TdistributionH", tstamp = 1, x = c(0, 1), p = c(0, 1))
  d2 <- new("TdistributionH", tstamp = 2, x = c(5, 6), p = c(0, 1))
  H <- new("HTS", epocs = 6, ListOfTimedElements = list(d1, d2))
  # the last element (d2) is equally distant from positions 2 and 4: a tie
  expect_no_warning(pr <- HTS.predict.knn(H, position = 6, k = 2))
  expect_s4_class(pr, "distributionH")
  expect_equal(pr@m, 0.5) # the elements following the two neighbours are d1
  expect_error(HTS.predict.knn(H, position = 3, k = 3))
})

test_that("arithmetic with numbers and registered distributions", {
  a <- distributionH(c(0, 2, 5), c(0, 0.3, 1))
  b <- distributionH(c(0, 1, 2), c(0, 0.3, 1))
  # same cdf: failed up to 1.0.9
  expect_equal((a - b)@x, c(0, 1, 3))
  # negative products and number minus distribution are valid distributions
  for (o in list(-1 * a, a * (-2), 10 - a)) {
    expect_true(all(diff(o@x) >= 0))
    expect_true(all(diff(o@p) >= 0))
    expect_equal(o@m, meanH(o))
    expect_equal(o@s, stdH(o), tolerance = 1e-12)
  }
  expect_equal((10 - a)@m, 10 - a@m)
  expect_equal((-2 * a)@s, 2 * a@s)
  expect_equal(compQ(-1 * a, 0.3), -compQ(a, 0.7))
  # adding a number shifts the distribution
  expect_equal((a + 3)@x, a@x + 3)
  expect_equal((a + 3)@m, a@m + 3)
})

test_that("fast computations give the same results as the definitions", {
  X <- BLOOD
  # sum of squares: matrix form vs definition with distances from the mean
  m <- WH.vec.mean(X[, 1])
  ssq <- sum(sapply(1:nrow(X@M), function(i) WassSqDistH(X@M[[i, 1]], m)))
  expect_equal(WH.SSQ(X[, 1])[1, 1], ssq, tolerance = 1e-8)
  # registerMH keeps the distributions
  R <- registerMH(X[1:4, 1])
  for (i in 1:4) {
    expect_equal(WassSqDistH(R@M[[i, 1]], X@M[[i, 1]]), 0, tolerance = 1e-10)
  }
  expect_true(is.registeredMH(R))
})

test_that("Kohonen maps run with proportional initial weights", {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })
  set.seed(1)
  for (sch in c(1, 2)) {
    expect_no_error(suppressWarnings(utils::capture.output(
      WH_2d_Adaptive_Kohonen_maps(BLOOD, net = list(xdim = 2, ydim = 2, topo = "rectangular"),
        repetitions = 1, schema = sch, init.weights = "PROP", simplify = TRUE, qua = 10)
    )))
  }
})
