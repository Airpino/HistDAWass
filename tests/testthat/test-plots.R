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
