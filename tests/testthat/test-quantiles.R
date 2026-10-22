test_that("compQ handles repeated values of x", {
  d <- distributionH(c(1, 2, 2, 3), c(0, 0.3, 0.6, 1))
  expect_equal(compQ(d, 0.7), 2.25)
  expect_equal(compQ(d, 0.5), 2)
  expect_equal(compQ(d, 0), 1)
  expect_equal(compQ(d, 1), 3)
  vp <- c(0.1, 0.35, 0.7, 0.95)
  expect_equal(sapply(vp, function(p) compQ(d, p)),
               HistDAWass:::COMP_Q_VECT(d@x, d@p, vp))
})

test_that("compP is the inverse of compQ", {
  d <- distributionH(c(1, 2, 3, 4), c(0, 0.5, 0.5, 1))
  expect_equal(compP(d, 3.5), 0.75, tolerance = 1e-7)
  d2 <- distributionH(c(0, 1, 5, 6), c(0, 0.1, 0.8, 1))
  for (p in c(0.05, 0.3, 0.9)) expect_equal(compP(d2, compQ(d2, p)), p)
})
