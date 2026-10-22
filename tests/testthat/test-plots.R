test_that("plot methods draw and return an invisible ggplot without warnings", {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })
  d <- BLOOD[1, 1]@M[[1]]
  for (ty in c("HISTO", "CDF", "QF", "DENS", "HBOXPLOT", "VBOXPLOT")) {
    expect_no_warning(v <- withVisible(plot(d, type = ty)))
    expect_false(v$visible)
    expect_s3_class(v$value, "ggplot")
  }
  for (ty in c("HISTO", "DENS", "BOXPLOT")) {
    expect_no_warning(suppressMessages(p <- plot(BLOOD[1:3, ], type = ty)))
    expect_s3_class(p, "ggplot")
  }
  for (ty in c("BOXPLOT", "VIOLIN")) {
    expect_no_warning(p <- plot(subsetHTS(RetHTS, 1, 20), type = ty, maxno.perplot = 10))
    expect_s3_class(p, "ggplot")
  }
})

test_that("summaryHTS returns one row per timestamp", {
  s <- summaryHTS(subsetHTS(RetHTS, 1, 10))
  expect_equal(nrow(s), 10)
  expect_equal(s$MED[1], compQ(RetHTS@data[[1]], 0.5))
})

test_that("histograms on the factorial plane have the requested size", {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })
  r <- WH.MultiplePCA(BLOOD, 1:3, graph = FALSE)
  co <- r$ind$coord[, 1:2]
  for (H in c(TRUE, FALSE)) {
    pl <- WH.plot_multiple_indivs(BLOOD, r, var = 1, HISTO = H)
    d <- pl[[1]]$data
    w <- tapply(d$x_o, d$ID_o, function(v) diff(range(v)))
    h <- tapply(d$y_o, d$ID_o, max) - co[, 2]
    expect_equal(mean(w) / diff(range(co[, 1])), 0.1)
    expect_equal(max(h) / diff(range(co[, 2])), 1 / 8)
  }
  pl <- WH.plot_multiple_indivs(BLOOD, r, var = 1, strx = 0.2, stry = 0.3)
  d <- pl[[1]]$data
  expect_equal(mean(tapply(d$x_o, d$ID_o, function(v) diff(range(v)))) /
    diff(range(co[, 1])), 0.2)
})
