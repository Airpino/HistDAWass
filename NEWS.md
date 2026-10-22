<!-- README.md is generated from README.Rmd. Please edit that file -->

# HistDAWass 1.0.9

## Bug fixes

* `c_WH_ADPT_KMEANS_TOTALSSQ()` (used by `WH_adaptive.kmeans()` and the
  adaptive Kohonen maps): the vector of cluster weights assigned to each row
  of the weight matrix was allocated with `k` (number of clusters) elements
  instead of `vars` (number of variables). When `k < vars` the remaining
  variables read uninitialised memory. This made the examples fail with the
  forthcoming Rcpp release (Rcpp PR #1513, which checks the length of vectors
  assigned into matrix rows). Thanks to the Rcpp team for the report and the
  patch.
* Empty bins (bins with zero weight) are now handled consistently: a weight of
  `1e-8` is assigned to each empty bin and the cdf is renormalised.
  This happens when a `distributionH` (or `TdistributionH`) object is created,
  in `data2hist()` (that used a weight of 0.001 counts), in `checkEmptyBins()`
  and in the internal C++ registration of distributions (that used `1e-14`).
  Up to 1.0.8 a distribution with an empty bin produced misaligned quantiles
  when several distributions were put on a common grid of cdf values, with
  errors in the computed distances (up to 7% in the `China_Seas` data).
* The `China_Seas`, `China_Month` and `Age_Pyramids_2014` datasets contained
  683 distributions with empty bins: they have been fixed accordingly.
* `compQ()` returned wrong quantiles for distributions with repeated values
  in `x` (e.g. `compQ(distributionH(c(1,2,2,3), c(0,.3,.6,1)), 0.7)` returned
  2 instead of 2.25).
* `compP()` returned wrong values for distributions with empty bins.
* The validity check of `distributionH` objects never checked that `p` is
  not decreasing; `distributionH(x)` with a single value now creates an
  impulse distribution as documented.
* `whichLowerEQThanX()` (C++) could read one element beyond the end of the
  vector; `c_DISTA_M2()` (C++) accessed the second element of a list of length
  one when data had a single variable.

## Graphics

* All the plot methods (`distributionH`, `TdistributionH`, `MatH`, `HTS`) and
  `plotPredVsObs()` now draw the plot on the current graphic device and return
  the `ggplot` object invisibly. Up to 1.0.8 the `distributionH` and `MatH`
  methods returned the object without drawing it, so nothing was shown inside
  loops, functions or `source()`.
* `WH.1d.PCA()` no longer opens new graphic devices with `dev.new()`: the four
  plots of each factorial plane are drawn on the current device (2 x 2
  panels), restoring the graphical parameters on exit. In the plot coloured by
  the standard deviations the densities were centred on the standard
  deviation instead of the mean.
* `WH.MultiplePCA()` gains a `graph` argument (default `TRUE`) passed to
  `FactoMineR::MFA()`.
* `WH.plot_multiple_Spanish.funs()` and `WH.plot_multiple_indivs()` no longer
  print their plots twice; `WH.plot_multiple_indivs()` failed with
  `stat` different from `"mean"`, and its cut of the highest densities was
  applied to the wrong column.
* The `VIOLIN` plot of `HTS` objects showed the data of the first block of
  timestamps in all the rows of the plot: fixed. `HTS` plots are now a single
  `ggplot` with one facet for each block of `maxno.perplot` timestamps.
* `plotPredVsObs()` returns a single `ggplot` (faceted by object) instead of
  drawing several plots in a grid of viewports.
* Updated to ggplot2 >= 3.4.0 (no more deprecation warnings with ggplot2 4.x):
  `linewidth` instead of `size` for lines, `guides(fill = "none")`, `aes()`
  with the `.data` pronoun instead of `aes_all()` and `with()`.
* The data of the plots are built in a vectorised way (quantiles computed in
  C++ for all the probabilities at once, no `rbind()` in loops):
  `summaryHTS()` and the `HTS` violin plot are about twice as fast.
* `data2hist(algo = "base")` no longer draws a histogram as a side effect.
* Removed leftover `browser()` calls from the code of the adaptive Kohonen
  maps and of the adaptive weights.

## Documentation

* The documentation has been revised and regenerated with roxygen2 8
  (`Config/roxygen2/version` replaces `RoxygenNote`): descriptive titles, a
  single usage for each method (no duplicated generic/method usage, no
  `initialize` methods in the usage of the classes), accurate return values,
  documented slots of the classes, dataset formats and examples, and an
  overview of the package in `?HistDAWass`.
* Related functions are documented together: the arithmetic operators
  (`?"Arith-distributionH"`), `WH.SSQ`/`WH.var.covar`/`WH.correlation` and
  their two-matrix versions, `WH.vec.sum`/`WH.vec.mean`, `WH.bind*`,
  `registerMH`/`is.registeredMH`, `get.cell.MatH`/`set.cell.MatH` and the
  `get.MatH.*` accessors.
* `\dontrun{}` examples replaced by `\donttest{}` (they all run).

## New features

* `checkEmptyBins()` has a method for `MatH` objects, useful for objects
  created with previous versions of the package.
* Added unit tests (testthat).

## Performance

* The C++ functions computing L2 Wasserstein distances between quantile
  functions (`c_MM_L2_SQ_WASS_D()`, `c_DISTA_M()`, `c_Fast_D_Mat()`,
  `c_STEP_2_ADA_KMEANS()`, `c_ComputeFast_L2_SQ_WASS_DMAT()`) no longer
  allocate temporary vectors nor access S4 slots in the inner loops, and
  `c_WH_ADPT_KMEANS_TOTALSSQ()` computes cluster-invariant distances once.
  Results are numerically identical to 1.0.8; on a 300 x 4 matrix of
  histograms `WH_MAT_DIST()` and `WH_hclust()` are about 20 times faster,
  `WH_kmeans()` about 6 times, `WH_adaptive.kmeans()` about 4 times and the
  Kohonen maps about 3 times faster.

# HistDAWass 1.0.7 What’s new

Minor bug fixes fo compatibility isseues with ggridges.

# HistDAWass 1.0.6 What’s new

Added the WH_WASS_MAT_DIST for computing a distance matrix

# HistDAWass 1.0.5 What’s new

Fixed some minor bugs. Added the names of columns and the rows in
weights matrices of adaptive distances based methods.

# HistDAWass 1.0.4 What’s new

Fixed some minor bugs.

# HistDAWass 1.0.3 What’s new

## Histogram-valued Data analysis using Wasserstein metric

Fixed some visualizations accordingly to the new ggplot2 release.

#HistDAWass 1.0.2 What’s new ##(Histogram-valued Data analysis using
Wasserstein metric)

Soelf Organizing maps have been improved for speed and aptions.

#HistDAWass 1.0.1 What’s new ##(Histogram-valued Data analysis using
Wasserstein metric)

Some bugs have been resolved when processing data with a single
variable.

#HistDAWass 1.0.0 What’s new ##(Histogram-valued Data analysis using
Wasserstein metric)

Now, a set of methods for MatH (WH.SSQ, WH.vec.mean, and all the methods
using them) are much faster due to RCPP implementation.

#HistDAWass 0.1.7 What’s new ##(Histogram-valued Data analysis using
Wasserstein metric)

Some general improvements have been added for reducing elaboration
times.

> unsupervised learning

WH_kmeans now returns more informative outputs

> principal components analysis

Multiple factor analysis has new functions for plotting results

#HistDAWass 0.1.6 What’s new ##(Histogram-valued Data analysis using
Wasserstein metric)

> unsupervised learning

fuzzy c-means and kohonen maps are 10 times faster

#HistDAWass 0.1.5 What’s new ##(Histogram-valued Data analysis using
Wasserstein metric)

> some minor bug fixed

#HistDAWass 0.1.4 What’s new ##(Histogram-valued Data analysis using
Wasserstein metric)

> some minor bug fixed

#HistDAWass 0.1.3 What’s new Some fixes to the worksapce

#HistDAWass 0.1.2 What’s new ##(Histogram-valued Data analysis using
Wasserstein metric) \> From raw data to histograms

A new function was added for generating a histogram from a set of raw
data

> Visualization

Now for plotting, `ggplot2` package is used generally.

A new plot for visualizing Histogram Time Series was added.

A new plot for comparing obversed vs predicted histogram was added.

A new plot for showing the differences between an observed and a
predicted histogram was added.

> Regression and predictive techniques

A new function for computing Goodness of Fit of a model was added
