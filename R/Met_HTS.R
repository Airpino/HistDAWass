## ---- Show overridding for Time series of TdistributionH and TMatH ----
setMethod("show",
  signature(object = "TdistributionH"),
  definition = function(object) {
    newdist <- distributionH(object@x, object@p)
    show(newdist)
  }
)
setMethod("show",
  signature(object = "TMatH"),
  definition = function(object) {
    cat("to be implemented")
  }
)
setMethod("show",
  signature(object = "HTS"),
  definition = function(object) {
    x <- summaryHTS(object)
    show(x[, 1:7])
  }
)

## --- Plot overloading for Timed TdistributionH  TMatH and HTS----
#' Plot a histogram with a time reference
#'
#' Plots a \code{TdistributionH} object: see \code{\link{plot,distributionH-method}}
#' for the types of plot.
#'
#' @param x A \code{TdistributionH} object.
#' @inheritParams plot-distributionH
#' @return The plot is drawn on the current device and the \code{ggplot} object is
#' returned invisibly.
#' @examples
#' td <- new("TdistributionH", tstamp = 1, x = c(1, 2, 3), p = c(0, 0.4, 1))
#' plot(td)
#' @name plot-TdistributionH
#' @aliases plot,TdistributionH-method
#' @export
setMethod(
  "plot",
  signature(x = "TdistributionH"),
  function(x, type = "HISTO", col = "green", border = "black") {
    newdist <- distributionH(x@x, x@p)
    plot(newdist, type = type, col = col, border = border)
  }
)
# setMethod("plot",
#           signature(x = "TMatH"),
#           function (x, y="missing", type="HISTO",col="green",border="black")
#           {
#             stop("to be implemented")
#           }
# )
## --- Plot overloading
#' Plot a histogram time series
#'
#' Plots a \code{HTS} object as a sequence of violin plots or boxplots, one for
#' each timestamp.
#'
#' @param x A \code{HTS} object.
#' @param y Not used.
#' @param type The type of plot: \code{"VIOLIN"} (default) or \code{"BOXPLOT"}.
#' @param border The colour of the borders (default \code{"black"}).
#' @param maxno.perplot The maximum number of timestamps in each row of the plot
#' (default 30): longer series are split into several rows.
#' @return The plot is drawn on the current device and the \code{ggplot} object is
#' returned invisibly, so that it can be further modified.
#' @examples
#' plot(subsetHTS(RetHTS, from = 1, to = 10))
#' \donttest{
#' plot(RetHTS, type = "BOXPLOT", border = "blue", maxno.perplot = 20)
#' plot(RetHTS, type = "VIOLIN", maxno.perplot = 20)
#' }
#' @name plot-HTS
#' @rdname plot-HTS
#' @aliases plot,HTS-method
#' @export
setMethod(
  "plot",
  signature(x = "HTS"),
  function(x, y = "missing", type = "VIOLIN", border = "black", maxno.perplot = 30) {
    .draw(.plot_HTS(x, type = type, border = border, maxno.perplot = maxno.perplot))
  }
)


#' Subset of a histogram time series
#'
#' Extracts the elements of a histogram time series between two positions.
#'
#' @param object A \code{HTS} object.
#' @param from,to Integers, the positions of the first and of the last element.
#' @return A \code{HTS} object.
#' @examples
#' SUB_RetHTS <- subsetHTS(RetHTS, from = 1, to = 20) # the first 20 elements
#' @name subsetHTS
#' @rdname subsetHTS-methods
#' @exportMethod subsetHTS
setGeneric("subsetHTS", function(object, from, to) standardGeneric("subsetHTS"))
#' @rdname subsetHTS-methods
#' @aliases subsetHTS,HTS,numeric,numeric-method
#' @usage NULL
setMethod(
  f = "subsetHTS",
  signature = c(object = "HTS", from = "numeric", to = "numeric"),
  function(object, from, to) {
    nel <- to - from + 1
    SubHTS <- new("HTS")
    cc <- 0
    for (i in from:to) {
      cc <- cc + 1
      SubHTS@data[[cc]] <- object@data[[i]]
    }
    return(SubHTS)
  }
)
