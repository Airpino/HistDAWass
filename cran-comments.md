## Submission of HistDAWass 1.0.9

This is a maintenance release, requested by the Rcpp team ahead of the next
Rcpp release (Rcpp PR #1513 checks that a vector assigned into a matrix row
has the row's length). The `WH_adaptive.kmeans()` example failed because
`c_WH_ADPT_KMEANS_TOTALSSQ()` assigned a vector of length k into a row of
length `vars`. The bug has been fixed (it also read uninitialised memory with
the released Rcpp when k < vars).

The release also fixes the handling of zero-weight bins of histograms and a
few other bugs, speeds up the C++ distance computations and updates the
graphics for ggplot2 >= 3.4.0 (see NEWS.md).

## Test environments

* Ubuntu 24.04, R 4.3.3, released Rcpp and Rcpp development branch of
  PR #1513 (1.1.2.4, with the sugar length checks enabled)
* (to add: win-builder R-devel, macOS via R-hub / mac-builder)

## R CMD check results

0 errors | 0 warnings | 2 NOTEs (environment-related)

* installed package size (libs 5.4Mb, debug symbols of the compiled code)
* compilation flags (`-mno-omit-leaf-frame-pointer`), set by the local R
  installation, not by the package

## Reverse dependencies

(to verify before submission with
`tools::package_dependencies("HistDAWass", reverse = TRUE, which = "all")`)
