test_that("WH_adaptive.kmeans works when k is lower than the number of variables", {
  # see Rcpp PR #1513: up to 1.0.8 a vector of length k was assigned into a
  # row of length vars
  set.seed(1)
  expect_no_error(res <- WH_adaptive.kmeans(BLOOD, k = 2, rep = 2, simplify = TRUE,
                                            qua = 10, standardize = TRUE))
  expect_true(res$quality > 0 && res$quality < 1)
  set.seed(1)
  expect_no_error(WH_adaptive.kmeans(BLOOD, k = 4, rep = 2, simplify = TRUE,
                                     qua = 10, standardize = TRUE))
})

test_that("distance matrix matches the pairwise Wasserstein distances", {
  X <- BLOOD[1:6, 1:2]
  D <- WH_MAT_DIST(X)
  for (i in 1:5) for (j in (i + 1):6) {
    d2 <- WassSqDistH(X@M[[i, 1]], X@M[[j, 1]]) + WassSqDistH(X@M[[i, 2]], X@M[[j, 2]])
    expect_equal(D[i, j]^2, d2, tolerance = 1e-10)
  }
})

test_that("clustering methods run and return a partition", {
  set.seed(2)
  res <- WH_kmeans(BLOOD, k = 3, rep = 2)
  expect_length(res$solution$IDX, nrow(BLOOD@M))
  expect_true(all(res$solution$IDX %in% 1:3))
  h <- WH_hclust(BLOOD)
  expect_s3_class(h, "hclust")
  set.seed(3)
  # the Kohonen maps draw the final map: use a temporary device
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })
  expect_no_error(WH_2d_Kohonen_maps(BLOOD, net = list(xdim = 2, ydim = 2, topo = "rectangular"),
                                     rep = 1, simplify = TRUE, qua = 10, standardize = TRUE))
})
