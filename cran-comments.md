## Submission of HistDAWass 1.0.10

This is a bug-fix release following 1.0.9, submitted shortly after it because
it fixes errors found after the 1.0.9 submission:

* `HTS.moving.averages()` assigned wrong timestamps to the averages;
* `HTS.predict.knn()` could select wrong neighbours in case of ties;
* the adaptive Kohonen maps failed with `init.weights = "PROP"`;
* the difference of two histograms sharing the same cdf failed, and products
  by negative numbers returned invalid objects.

It also makes the computations with histograms having different bins much
faster (up to 25 times), with identical results, and adds five vignettes
(see NEWS.md).

## Test environments

* Ubuntu 24.04, R 4.3.3, released Rcpp and Rcpp development branch of
  PR #1513 (1.1.2.4, with the sugar length checks enabled), ggplot2 3.4.4 and
  4.0.3
* win-builder: R-devel and R-release
* macOS (mac-builder), R-release

## R CMD check results

0 errors | 0 warnings | NOTEs related only to the local environment
(installed size of the compiled code with debug symbols, compilation flags
set by the local R installation).

## Reverse dependencies

There are two reverse dependencies: dataSDA (Imports) and ggInterval
(Suggests). The test suite and the examples of dataSDA give the same results
with HistDAWass 1.0.8 and 1.0.10; ggInterval only uses the `BLOOD` dataset
and the `get.MatH.*` accessors, which are unchanged.
