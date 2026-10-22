## Submission of HistDAWass 1.0.10

This is a bug-fix release following 1.0.9, submitted shortly after it because
it fixes errors found after the 1.0.9 submission:

* `HTS.moving.averages()` assigned wrong timestamps to the averages;
* `HTS.predict.knn()` could select wrong neighbours in case of ties;
* the adaptive Kohonen maps failed with `init.weights = "PROP"`;
* the difference of two histograms sharing the same cdf failed, and products
  by negative numbers returned invalid objects.

It also makes the computations with histograms having different bins much
faster (up to 25 times), with identical results (see NEWS.md).

## Test environments

* Ubuntu 24.04, R 4.3.3, released Rcpp and Rcpp development branch of
  PR #1513 (1.1.2.4, with the sugar length checks enabled), ggplot2 3.4.4 and
  4.0.3
* (to add: win-builder R-devel, macOS via R-hub / mac-builder)

## R CMD check results

0 errors | 0 warnings | NOTEs related only to the local environment
(installed size of the compiled code with debug symbols, compilation flags
set by the local R installation).

## Reverse dependencies

(to verify before submission with
`tools::package_dependencies("HistDAWass", reverse = TRUE, which = "all")`)
