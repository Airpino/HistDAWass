# A class for histogram-valued data ---------------------------------------
#' Histogram objects: the distributionH class
#'
#' The \code{distributionH} class describes a histogram (a distribution) by means
#' of its cumulative distribution function (cdf): the bounds of the bins and the
#' values of the cdf at the bounds. All the methods of the package are based on
#' the L2 Wasserstein distance between distributions.
#'
#' @slot x A non-decreasing numeric vector: the bounds of the bins (the domain of
#' the distribution).
#' @slot p A non-decreasing numeric vector of the same length of \code{x}, from 0
#' to 1: the cdf at the values in \code{x}.
#' @slot m The mean of the distribution.
#' @slot s The standard deviation of the distribution.
#' @param x A numeric vector, the bounds of the bins, in non-decreasing order.
#' A single value defines an impulse distribution.
#' @param p A numeric vector of the same length of \code{x}, the values of the
#' cdf at \code{x}, non-decreasing from 0 to 1.
#' @section Creating objects:
#' Objects are created with \code{distributionH(x, p)} or, equivalently, with
#' \code{new("distributionH", x, p, m, s)}, where the optional \code{m} and
#' \code{s} (mean and standard deviation) are computed when not given.
#'
#' Empty bins, i.e. bins with zero weight (two consecutive equal values in
#' \code{p}), may produce numerical problems. Since version 1.0.9 each empty bin
#' gets a weight of \code{1e-8} and the cdf is renormalised when the object is
#' created (see also \code{\link{checkEmptyBins}}).
#' @return \code{distributionH()} returns a \code{distributionH} object.
#' @references Irpino, A., Verde, R. (2015) \emph{Basic statistics for
#' distributional symbolic variables: a new metric-based approach}, Advances in
#' Data Analysis and Classification, 9(2), 143-175.
#' \doi{10.1007/s11634-014-0176-4}
#' @seealso \code{\link{data2hist}} to build a histogram from raw data;
#' \code{\link{meanH}}, \code{\link{stdH}}, \code{\link{plot,distributionH-method}}.
#' @examples
#' # the histogram
#' # | Bins  | Prob | cdf |
#' # | [1,2) | 0.4  | 0.4 |
#' # | [2,3] | 0.6  | 1.0 |
#' mydist <- distributionH(x = c(1, 2, 3), p = c(0, 0.4, 1))
#' mydist
#' str(mydist)
#' @name distributionH-class
#' @rdname distributionH-class
#' @aliases distributionH-class
#' @exportClass distributionH
setClass(
  Class = "distributionH",
  representation = representation(
    x = "numeric",
    p = "numeric",
    m = "numeric",
    s = "numeric"
  ),
  validity = function(object) {
    if (length(object@x) == 0) {
      return(TRUE)
    }
    .check_xp(object@x, object@p)
  }
)
#' @rdname distributionH-class
#' @aliases initialize,distributionH-method
#' @usage NULL
# showClass("distributionH")
setMethod("initialize", "distributionH",
  definition = function(.Object, x = numeric(0), p = numeric(0), m = numeric(0), s = numeric(0)) {
    # a single value is a degenerate distribution (an impulse)
    if (length(x) == 1L) {
      x <- c(x, x)
      p <- c(0, 1)
    }
    .Object@x <- as.numeric(x)
    .Object@p <- as.numeric(p)
    if (length(x) > 0) {
      validObject(.Object)
      # zero-weight (empty) bins get a weight of 1e-8 (since v. 1.0.9)
      fx <- .fix_empty_bins_p(.Object@p)
      if (fx$fixed) {
        .Object@p <- fx$p
        m <- numeric(0) # mean and std must be recomputed
        s <- numeric(0)
      }
      if (length(m) == 0) .Object@m <- meanH(.Object) else .Object@m <- m
      if (length(s) == 0) .Object@s <- stdH(.Object) else .Object@s <- s
    }

    return(.Object)
  }
)

# A class of a matrix of histogram-valued data ----------------------------
#' Matrices of histograms: the MatH class
#'
#' The \code{MatH} class describes a matrix of histograms: each cell contains a
#' \code{distributionH} object, rows are the individuals and columns the
#' (histogram-valued) variables.
#'
#' @slot M A matrix of mode list, each cell containing a \code{distributionH}
#' object.
#' @section Creating objects:
#' Objects are created with \code{MatH()} or with
#' \code{new("MatH", nrows, ncols, ListOfDist, names.rows, names.cols, by.row)},
#' where \code{ListOfDist} is a list of \code{distributionH} objects and
#' \code{names.rows}, \code{names.cols} are the names of rows and columns.
#' @references Irpino, A., Verde, R. (2015) \emph{Basic statistics for
#' distributional symbolic variables: a new metric-based approach}, Advances in
#' Data Analysis and Classification, 9(2), 143-175.
#' \doi{10.1007/s11634-014-0176-4}
#' @seealso \code{\link{get.MatH.main.info}}, \code{\link{plot,MatH-method}},
#' \code{\link{WH.vec.mean}}, \code{\link{WH.var.covar}}.
#' @examples
#' ## a list of six distributionH objects
#' ListOfDist <- vector("list", 6)
#' ListOfDist[[1]] <- distributionH(c(1, 2, 3), c(0, 0.4, 1))
#' ListOfDist[[2]] <- distributionH(c(7, 8, 10, 15), c(0, 0.2, 0.7, 1))
#' ListOfDist[[3]] <- distributionH(c(9, 11, 20), c(0, 0.5, 1))
#' ListOfDist[[4]] <- distributionH(c(2, 5, 8), c(0, 0.3, 1))
#' ListOfDist[[5]] <- distributionH(c(8, 10, 15), c(0, 0.75, 1))
#' ListOfDist[[6]] <- distributionH(c(20, 22, 24), c(0, 0.12, 1))
#' ## a 3 x 2 MatH filled by columns
#' MyMAT <- MatH(
#'   x = ListOfDist, nrows = 3, ncols = 2,
#'   rownames = c("I1", "I2", "I3"), varnames = c("Var1", "Var2")
#' )
#' MyMAT
#' ## an empty 10 x 4 matrix of histograms
#' MAT <- MatH(nrows = 10, ncols = 4)
#' @name MatH-class
#' @rdname MatH-class
#' @aliases MatH-class
#' @exportClass MatH

setClass(
  Class = "MatH",
  representation = representation(M = "matrix"),
)
#' @rdname MatH-class
#' @aliases initialize,MatH-method
#' @usage NULL
setMethod("initialize", "MatH",
  definition = function(.Object, nrows = 1, ncols = 1, ListOfDist = NULL, names.rows = NULL, names.cols = NULL, by.row = FALSE) {
    tt <- list(new("distributionH"))
    .Object@M <- matrix(tt, nrows, ncols)

    if (length(ListOfDist) > 0) {
      nOBJ <- length(ListOfDist)
      if (by.row) {
        count <- 0
        for (i in 1:nrows) {
          for (j in 1:ncols) {
            count <- count + 1
            if (count == nOBJ) count <- 1
            .Object@M[i, j][[1]] <- ListOfDist[[count]]
          }
        }
      }
      else {
        count <- 0
        for (j in 1:ncols) {
          for (i in 1:nrows) {
            count <- count + 1
            if (count > nOBJ) count <- 1
            .Object@M[i, j][[1]] <- ListOfDist[[count]]
          }
        }
      }
    }
 
    if (length(names.rows) > 0) {
      count <- 0
      rnames <- vector("list", nrows)

      for (i in 1:nrows) {
        count <- count + 1

        if (count > length(names.rows)) {
          rnames[[count]] <- paste("I", count, sep = "")
        }
        else {
          rnames[[count]] <- names.rows[[count]]
        }
      }
      rownames(.Object@M) <- rnames
    }
    else {
      rownames(.Object@M) <- paste("I", 1:nrows, sep = "")
    }

    if (length(names.cols) > 0) {
      count <- 0
      cnames <- vector("list", ncols)
      for (i in 1:ncols) {
        count <- count + 1
        if (count > length(names.cols)) {
          cnames[[count]] <- paste("X", count, sep = "")
        }
        else {
          cnames[[count]] <- names.cols[[count]]
        }
        colnames(.Object@M) <- cnames
      }
    } else {
      colnames(.Object@M) <- paste("X", 1:ncols, sep = "")
    }

    return(.Object)
  }
)
## Classes for Histogram Time Series ------
## A single distribution with a time stamp ------
#' Histograms with a time reference: the TdistributionH class
#'
#' The \code{TdistributionH} class extends \code{distributionH} with a time
#' reference: a timestamp or a time period. It is the element of a histogram time
#' series (\code{\link{HTS-class}}).
#'
#' @slot tstamp A number, the timestamp.
#' @slot period A list with two elements, \code{start} and \code{end}, the bounds
#' of the time period (alternative to \code{tstamp}).
#' @slot x The bounds of the bins (see \code{\link{distributionH-class}}).
#' @slot p The cdf at \code{x}.
#' @slot m The mean.
#' @slot s The standard deviation.
#' @section Creating objects:
#' Objects are created with
#' \code{new("TdistributionH", tstamp, period, x, p, m, s)}, where \code{x},
#' \code{p}, \code{m} and \code{s} are as in \code{\link{distributionH}}.
#' @seealso \code{\link{HTS-class}}, \code{\link{plot,TdistributionH-method}}
#' @examples
#' td <- new("TdistributionH", tstamp = 1, x = c(1, 2, 3), p = c(0, 0.4, 1))
#' td@tstamp
#' @name TdistributionH-class
#' @rdname TdistributionH-class
#' @aliases TdistributionH-class
#' @exportClass TdistributionH
setClass(
  Class = "TdistributionH",
  representation = representation(
    tstamp = "numeric",
    period = "list"
  ),
  contains = "distributionH",
  validity = function(object) {
    if (length(object@tstamp) > 1) {
      stop("time stamp slot may contain a single value")
    }
    if (length(object@period) > 2) {
      stop("period slot may contain no more that two ordered values")
    }
    if ((length(object@period) == 2) && (object@period[[1]] > object@period[[2]])) {
      stop("period slot error: the end of period is lower than the starting one")
    }
  }
)
## Initialize ------
#' @rdname TdistributionH-class
#' @aliases initialize,TdistributionH-method
#' @usage NULL
setMethod("initialize", "TdistributionH",
  definition = function(.Object,
                        tstamp = numeric(0),
                        period = list(start = -Inf, end = -Inf),
                        x = numeric(0),
                        p = numeric(0),
                        m = numeric(0),
                        s = numeric(0)) {
    if (length(x) == 1L) {
      x <- c(x, x)
      p <- c(0, 1)
    }
    .Object@x <- as.numeric(x)
    .Object@p <- as.numeric(p)
    .Object@m <- m
    .Object@s <- s
    .Object@tstamp <- tstamp
    .Object@period <- period
    if (length(x) > 0) {
      validObject(.Object)
      # zero-weight (empty) bins get a weight of 1e-8 (since v. 1.0.9)
      fx <- .fix_empty_bins_p(.Object@p)
      if (fx$fixed) {
        .Object@p <- fx$p
        m <- numeric(0)
        s <- numeric(0)
      }
      if (length(m) == 0) .Object@m <- meanH(.Object) else .Object@m <- m
      if (length(s) == 0) .Object@s <- stdH(.Object) else .Object@s <- s
    }
    return(.Object)
  }
)
# Coerce a TdistributionH into a distributionH
setAs(
  from = "TdistributionH", to = "distributionH",
  function(from, to) {
    to <- new("distributionH",
      x = from@x,
      p = from@p,
      m = from@m,
      s = from@s
    )
    return(to)
  }
)

## A mulvariate MatH with a time stamp  ----
#' Matrices of histograms with a time reference: the TMatH class
#'
#' The \code{TMatH} class extends \code{MatH} with a time reference (timestamps or
#' time periods).
#'
#' @slot tstamp A numeric vector of timestamps.
#' @slot period A list of time periods (alternative to \code{tstamp}).
#' @slot M The matrix of histograms, as in \code{\link{MatH-class}}.
#' @section Creating objects:
#' Objects are created with \code{new("TMatH", mat, tstamp, period)}, where
#' \code{mat} is a \code{MatH} object.
#' @seealso \code{\link{MatH-class}}, \code{\link{HTS-class}}
#' @name TMatH-class
#' @rdname TMatH-class
#' @aliases TMatH-class
#' @exportClass TMatH
setClass(
  Class = "TMatH",
  representation = representation(
    tstamp = "numeric",
    period = "list"
  ),
  contains = "MatH",
  validity = function(object) {
    if (length(tstamp) > 1) {
      stop("time stamp slot may contain a single value")
    }
    if (length(period) > 2) {
      stop("period slot may contain no more that two ordered values")
    }
    if ((length(period) == 2) && (period[1] > period[2])) {
      stop("period slot error: the end of period is lower than the starting one")
    }
  }
)
## Initialize TMatH ------
#' @rdname TMatH-class
#' @aliases initialize,TMatH-method
#' @usage NULL
setMethod("initialize", "TMatH",
  definition = function(.Object,
                        tstamp = numeric(0),
                        period = list(start = -Inf, end = -Inf),
                        mat = new("MatH")) {
    .Object@M <- mat
    .Object@tstamp <- tstamp
    .Object@period <- period

    return(.Object)
  }
)

## A Histogram Time Series HTS-----
#' Histogram time series: the HTS class
#'
#' The \code{HTS} class describes a histogram time series, i.e. a sequence of
#' histograms (\code{TdistributionH} objects) observed along time.
#'
#' @slot data A list of \code{TdistributionH} objects.
#' @section Creating objects:
#' Objects are created with \code{new("HTS", epocs, ListOfTimedElements)}, where
#' \code{epocs} is the number of time points and \code{ListOfTimedElements} a
#' list of \code{TdistributionH} objects.
#' @seealso \code{\link{RetHTS}} for an example dataset;
#' \code{\link{subsetHTS}}, \code{\link{plot,HTS-method}},
#' \code{\link{HTS.moving.averages}}, \code{\link{HTS.exponential.smoothing}},
#' \code{\link{HTS.predict.knn}}.
#' @examples
#' RetHTS
#' subsetHTS(RetHTS, from = 1, to = 5)
#' @name HTS-class
#' @rdname HTS-class
#' @aliases HTS-class
#' @exportClass HTS
setClass(
  Class = "HTS",
  representation = representation(data = "vector"),
  validity = function(object) {
    if (length(object@data) > 0) {
      count <- 0
      while (count < length(object@data)) {
        count <- count + 1
        if (count == 1) {
          type <- is(object@data[[1]])[1]
          if (!(type == "TdistributionH") || (type == "TMatH")) {
            stop("data must be  TdistributionH or TMatH objects")
          }
        }
        else {
          if (is(object@data[[count]])[1] != type) {
            stop("all data must be of the same type")
          }
          if (type == "TMatH") {
            if ((nrow(object@data[[i]]@M) != nrow(object@data[[i - 1]]@M)) ||
              (ncol(object@data[[i]]@M) != ncol(object@data[[i - 1]]@M))) {
              stop("all TMatH objects must have the same dimensions")
            }
          }
        }
      }
    }
  }
)
## Initialize HTS ------
#' @rdname HTS-class
#' @aliases initialize,HTS-method
#' @usage NULL
setMethod("initialize", "HTS",
  definition = function(.Object, epocs = 1, ListOfTimedElements = c(new("TdistributionH"))) {
    counts <- 0
    for (i in 1:epocs) {
      counts <- counts + 1
      if (counts > length(ListOfTimedElements)) {
        counts <- 1
      }
      .Object@data <- c(.Object@data, ListOfTimedElements[[counts]])
    }
    validObject(.Object)
    return(.Object)
  }
)
