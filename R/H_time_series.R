## functions for analysis and plotting
## smoothing techniques ----
## moving averages -----
#' Moving averages of a histogram time series
#'
#' (Beta version) Smoothing of a histogram time series by (weighted) moving
#' averages, based on the L2 Wasserstein distance.
#'
#' @param HTS A \code{HTS} object.
#' @param k An odd integer, the number of elements of each average (default 3).
#' @param weights A vector of \code{k} non-negative weights (default: equal
#' weights).
#' @return A list with \code{k}, the (normalised) \code{weights} and
#' \code{AveragedHTS}, the smoothed \code{HTS} object.
#' @seealso \code{\link{HTS.exponential.smoothing}}
#' @examples
#' mov.av.smoothed <- HTS.moving.averages(HTS = RetHTS, k = 5)
#' mov.av.smoothed$AveragedHTS
#' @export
HTS.moving.averages <- function(HTS, k = 3, weights = rep(1, k)) {
  # check and standardize weights
  if (min(weights) < 0) {
    stop("Weights must be positive")
  }
  weights <- weights / sum(weights)
  if (k %% 2 == 0) { # the number of elements is even
    stop("for k even it is not yet implemented")
  }
  else {
    epocs <- length(HTS@data) - k + 1
    AveragedHTS <- new("HTS", epocs = epocs)
    counts <- 0
    shift <- (k + 1) / 2
    for (i in 1:epocs) {
      AveragedHTS@data[[i]]@tstamp <- HTS@data[[i + shift]]@tstamp
      AveragedHTS@data[[i]]@period <- HTS@data[[i + shift]]@period
      tmp <- weights[1] * HTS@data[[i]]
      for (j in 2:k) {
        tmp <- tmp + weights[j] * HTS@data[[i + j - 1]]
      }
      # write in the new HTS
      AveragedHTS@data[[i]]@x <- tmp@x
      AveragedHTS@data[[i]]@p <- tmp@p
      AveragedHTS@data[[i]]@m <- tmp@m
      AveragedHTS@data[[i]]@s <- tmp@s
    }
  }
  return(list(k = k, weights = weights, AveragedHTS = AveragedHTS))
}
## exponential smoothing ----
#' Exponential smoothing of a histogram time series
#'
#' (Beta version) Exponential smoothing of a histogram time series, based on the
#' L2 Wasserstein distance.
#'
#' @param HTS A \code{HTS} object.
#' @param alpha A number between 0 and 1, the smoothing parameter (default 0.9).
#' @return A list with \code{smoothing.alpha}, the smoothing parameter, and
#' \code{SmoothedHTS}, the smoothed \code{HTS} object.
#' @seealso \code{\link{HTS.moving.averages}}
#' @examples
#' mov.expo.smooth <- HTS.exponential.smoothing(HTS = RetHTS, alpha = 0.8)
#' mov.expo.smooth$SmoothedHTS
#' @export
HTS.exponential.smoothing <- function(HTS, alpha = 0.9) {
  if ((alpha < 0) || (alpha > 1)) {
    stop("The smoothing parameter alpha must be between 0 and 1")
  }
  epocs <- length(HTS@data)
  SmoothHTS <- new("HTS", epocs = epocs)
  s <- HTS@data[[1]]
  SmoothHTS@data[[1]]@tstamp <- HTS@data[[1]]@tstamp
  SmoothHTS@data[[1]]@period <- HTS@data[[1]]@period
  SmoothHTS@data[[1]]@x <- s@x
  SmoothHTS@data[[1]]@p <- s@p
  SmoothHTS@data[[1]]@m <- s@m
  SmoothHTS@data[[1]]@s <- s@s
  for (i in 2:epocs) {
    SmoothHTS@data[[i]]@tstamp <- HTS@data[[i]]@tstamp
    SmoothHTS@data[[i]]@period <- HTS@data[[i]]@period
    s <- (1 - alpha) * s + alpha * HTS@data[[i]]
    SmoothHTS@data[[i]]@x <- s@x
    SmoothHTS@data[[i]]@p <- s@p
    SmoothHTS@data[[i]]@m <- s@m
    SmoothHTS@data[[i]]@s <- s@s
  }
  return(list(smoothing.alpha = alpha, SmoothedHTS = SmoothHTS))
}
## prediction techniques ----
## predicting via k-nn ----
#' K-NN forecast of a histogram time series
#'
#' (Beta version) Forecasts a histogram of a time series with the k-nearest
#' neighbours method, based on the L2 Wasserstein distance: the forecast is the
#' mean of the histograms following the \code{k} histograms most similar to the
#' one at \code{position}.
#'
#' @param HTS A \code{HTS} object.
#' @param position An integer, the position of the last observed histogram
#' (default: the last one).
#' @param k The number of neighbours (default 3).
#' @return A \code{distributionH} object, the forecast.
#' @references Arroyo, J., Mate, C. (2009). Forecasting histogram time series
#' with k-nearest neighbours methods. International Journal of Forecasting,
#' 25(1), 192-207. \doi{10.1016/j.ijforecast.2008.07.003}
#' @examples
#' prediction <- HTS.predict.knn(HTS = RetHTS, position = 108, k = 3)
#' @export
HTS.predict.knn <- function(HTS, position = length(HTS@data), k = 3) {
  if (length(HTS@data) > position) {
    HTS@data <- HTS@data[1:position]
  }
  # compute distances
  dista <- vector(mode = "numeric", position - 1)
  for (i in 1:(position - 1)) {
    dista[i] <- WassSqDistH(HTS@data[[position]], HTS@data[[i]])
  }
  firstksorted <- sort(dista)[1:k]
  positions <- vector(mode = "numeric", k)
  vectorofdistr <- new(Class = "MatH", nrows = k, ncols = 1)
  for (v in 1:k) {
    positions[v] <- which(dista == (firstksorted[v])) + 1
    vectorofdistr@M[v, 1][[1]] <- new(
      "distributionH",
      HTS@data[[positions[v]]]@x,
      HTS@data[[positions[v]]]@p,
      HTS@data[[positions[v]]]@m,
      HTS@data[[positions[v]]]@s
    )
  }
  # compute the mean of histograms at position positions
  pred <- WH.vec.mean(vectorofdistr)
  return(pred)
}
