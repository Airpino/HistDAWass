# Constructor
#' @rdname distributionH-class
#' @export
#' @import methods
distributionH <- function(x = numeric(0), p = numeric(0)) {
  object <- new("distributionH", x = x, p = p)
  return(object)
}

# Get methods -----

#' Mean of a histogram (stored value)
#'
#' Returns the mean stored in the \code{m} slot of a \code{distributionH} object.
#' To compute it from the distribution use \code{\link{meanH}}.
#'
#' @param object A \code{distributionH} object.
#' @return A number, the mean of the distribution.
#' @seealso \code{\link{meanH}}, \code{\link{get.s}}
#' @examples
#' D <- distributionH(x = c(1, 2, 3, 4), p = c(0, 0.2, 0.6, 1))
#' get.m(D)
#' @name get.m
#' @rdname get.m-methods
#' @exportMethod get.m
setGeneric("get.m", function(object) standardGeneric("get.m"))
#' @rdname get.m-methods
#' @usage NULL
setMethod(
  "get.m", "distributionH",
  function(object) {
    if (!is.null(object)) {
      return(object@m)
    } else {
      return(NA)
    }
  }
)

#' Standard deviation of a histogram (stored value)
#'
#' Returns the standard deviation stored in the \code{s} slot of a
#' \code{distributionH} object. To compute it from the distribution use
#' \code{\link{stdH}}.
#'
#' @param object A \code{distributionH} object.
#' @return A non-negative number, the standard deviation of the distribution.
#' @seealso \code{\link{stdH}}, \code{\link{get.m}}
#' @examples
#' D <- distributionH(x = c(1, 2, 3, 4), p = c(0, 0.2, 0.6, 1))
#' get.s(D)
#' @name get.s
#' @rdname get.s-methods
#' @exportMethod get.s
setGeneric("get.s", function(object) standardGeneric("get.s"))

#' @rdname get.s-methods
#' @usage NULL
setMethod(
  "get.s", "distributionH",
  function(object) {
    if (!is.null(object)) {
      return(object@s)
    } else {
      return(NA)
    }
  }
)
#' Cumulative distribution function of a histogram
#'
#' Returns the points of the cumulative distribution function (cdf) of a
#' \code{distributionH} object.
#'
#' @param object A \code{distributionH} object.
#' @return A \code{data.frame} with two columns: \code{x}, the bounds of the
#' bins, and \code{p}, the values of the cdf at \code{x}.
#' @seealso \code{\link{get.histo}}
#' @examples
#' D <- distributionH(x = c(1, 2, 3, 4), p = c(0, 0.2, 0.6, 1))
#' get.distr(D)
#' @name get.distr
#' @rdname get.distr-methods
#' @exportMethod get.distr
setGeneric("get.distr", function(object) standardGeneric("get.distr"))

#' @rdname get.distr-methods
#' @usage NULL
setMethod(
  "get.distr", "distributionH",
  function(object) {
    if (!is.null(object)) {
      MAT <- cbind(object@x, object@p)
      colnames(MAT) <- c("x", "p")
      return(MAT = as.data.frame(MAT))
    } else {
      return(NA)
    }
  }
)
#' Bins of a histogram
#'
#' Returns the bins of a \code{distributionH} object and their probabilities.
#'
#' @param object A \code{distributionH} object.
#' @return A \code{data.frame} with one row for each bin and three columns:
#' \code{min.x} and \code{max.x}, the bounds of the bin, and \code{p}, its
#' probability (relative frequency).
#' @seealso \code{\link{get.distr}}
#' @examples
#' D <- distributionH(x = c(1, 2, 3, 4), p = c(0, 0.2, 0.6, 1))
#' get.histo(D)
#' @name get.histo
#' @rdname get.histo-methods
#' @exportMethod get.histo
setGeneric("get.histo", function(object) standardGeneric("get.histo"))

# Get the distribution of a histogram
#' @rdname get.histo-methods
#' @usage NULL
setMethod(
  "get.histo", "distributionH",
  function(object) {
    if (!is.null(object)) {
      M <- get.distr(object)
      MAT <- cbind(
        M$x[1:(length(M$x) - 1)],
        M$x[2:length(M$x)],
        M$p[2:length(M$p)] - M$p[1:(length(M$p) - 1)]
      )
      colnames(MAT) <- c("min.x", "max.x", "p")
      return(MAT = as.data.frame(MAT))
    } else {
      return(NA)
    }
  }
)

# Basic statistics of distributions --------
#' Mean of a histogram
#'
#' Computes the mean (first moment) of a \code{distributionH} object.
#'
#' @param object A \code{distributionH} object.
#' @return A number, the mean of the distribution.
#' @seealso \code{\link{stdH}}, \code{\link{skewH}}, \code{\link{kurtH}}
#' @examples
#' mydist <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' meanH(mydist) # 4.4
#' @name meanH
#' @rdname meanH-methods
#' @exportMethod meanH
setGeneric("meanH", function(object) standardGeneric("meanH"))
#' Standard deviation of a histogram
#'
#' Computes the standard deviation (square root of the second centred moment)
#' of a \code{distributionH} object.
#'
#' @param object A \code{distributionH} object.
#' @return A non-negative number, the standard deviation of the distribution.
#' @seealso \code{\link{meanH}}, \code{\link{skewH}}, \code{\link{kurtH}}
#' @examples
#' mydist <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' stdH(mydist) # 2.563851
#' @name stdH
#' @rdname stdH-methods
#' @exportMethod stdH
setGeneric("stdH", function(object) standardGeneric("stdH"))
#' Skewness of a histogram
#'
#' Computes the skewness (third standardised moment) of a \code{distributionH}
#' object.
#'
#' @param object A \code{distributionH} object.
#' @return A number, the skewness index.
#' @seealso \code{\link{meanH}}, \code{\link{stdH}}, \code{\link{kurtH}}
#' @examples
#' mydist <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' skewH(mydist) # 0.7148855
#' @name skewH
#' @rdname skewH-methods
#' @exportMethod skewH
setGeneric("skewH", function(object) standardGeneric("skewH"))
#' Kurtosis of a histogram
#'
#' Computes the kurtosis (fourth standardised moment) of a \code{distributionH}
#' object.
#'
#' @param object A \code{distributionH} object.
#' @return A number, the kurtosis index (3 for a Gaussian distribution).
#' @seealso \code{\link{meanH}}, \code{\link{stdH}}, \code{\link{skewH}}
#' @examples
#' mydist <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' kurtH(mydist) # 2.133533
#' @name kurtH
#' @rdname kurtH-methods
#' @exportMethod kurtH
setGeneric("kurtH", function(object) standardGeneric("kurtH"))
#' Centres, radii and weights of the bins of a histogram
#'
#' Computes the centre (midpoint), the radius (half-width) and the weight
#' (probability) of each bin of a \code{distributionH} object. This
#' representation allows a fast computation of the statistics based on the L2
#' Wasserstein distance.
#'
#' @param object A \code{distributionH} object.
#' @return A list with three numeric vectors: \code{Centers}, the midpoints of
#' the bins; \code{Radii}, the half-widths of the bins; \code{Weights}, the
#' probabilities of the bins (summing to one).
#' @references Irpino, A., Verde, R., Lechevallier, Y. (2006) \emph{Dynamic
#' clustering of histograms using Wasserstein metric}, In: Proceedings of
#' COMPSTAT 2006, Physica-Verlag, 869-876.
#' @examples
#' mydist <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' crwtransform(mydist)
#' @name crwtransform
#' @rdname crwtransform-methods
#' @exportMethod crwtransform
setGeneric("crwtransform", function(object) standardGeneric("crwtransform"))

#' @rdname meanH-methods
#' @usage NULL
setMethod(
  "meanH", "distributionH", # Mean of a distr --------
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      m <- M_STD_H(object)
      # resu=crwtransform(object)
      # c=resu[[1]]
      # w=resu[[3]]
      # m=t(w)%*%c
      # m=as.numeric(m)
      return(m[1])
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)
#' @rdname stdH-methods
#' @usage NULL
setMethod(
  "stdH", "distributionH", # Std of a distr --------
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      m <- M_STD_H(object)
      #   #resu=crwtransform(object)
      #   c=0.5*(object@x[2:length(object@x)]+object@x[1:(length(object@x)-1)])
      #   r=0.5*(object@x[2:length(object@x)]-object@x[1:(length(object@x)-1)]) # resu[[2]]
      #   w=(object@p[2:length(object@p)]-object@p[1:(length(object@p)-1)])#resu[[3]]
      #   std=sqrt(abs(sum(w*c^2+1/3*w*r^2)-(sum(w*c))^2))
      # #  std=as.numeric(std)
      return(m[2])
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)
#' @rdname skewH-methods
#' @usage NULL
setMethod(
  "skewH", "distributionH", # skewness of a distr --------
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      resu <- crwtransform(object)
      cs <- (resu[[1]] - object@m) / (object@s)
      rs <- resu[[2]] / (object@s)
      w <- resu[[3]]
      sk <- t(w) %*% (cs * (rs^2 + cs^2))
      sk <- as.numeric(sk)
      return(sk)
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)
#' @rdname kurtH-methods
#' @usage NULL
setMethod(
  "kurtH", "distributionH", # Kurtosis of a distr --------
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      resu <- crwtransform(object)
      cs <- (resu[[1]] - object@m) / (object@s)
      rs <- resu[[2]] / (object@s)
      w <- resu[[3]]
      ku <- 0.2 * t(w) %*% (5 * cs^4 + 10 * (cs^2) * rs^2 + rs^4)
      ku <- as.numeric(ku)
      return(ku)
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)
#' @rdname crwtransform-methods
#' @usage NULL
setMethod(
  "crwtransform", "distributionH", # crw transform -----
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      nv <- length(object@x)
      c <- (object@x[2:nv] + object@x[1:(nv - 1)]) / 2
      r <- (object@x[2:nv] - object@x[1:(nv - 1)]) / 2
      w <- object@p[2:nv] - object@p[1:(nv - 1)]
      resu <- list(Centers = c, Radii = r, Weights = w)
      return(resu)
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)
# Overloading of the sum of two distribution according to the L2 w --------
#' Arithmetic operators for histograms
#'
#' Sum, difference and product of \code{distributionH} objects and numbers,
#' consistent with the L2 Wasserstein metric: the operations are applied to the
#' quantile functions of the distributions.
#' \itemize{
#' \item \code{e1 + e2}: the sum of the quantile functions. A number is treated as
#' an impulse distribution, so \code{x + a} shifts \code{x} by \code{a}.
#' \item \code{e1 - e2}: the difference of the quantile functions. The difference
#' of two distributions is a distribution only if \code{e1} is more variable
#' than \code{e2} (the result must be a non-decreasing quantile function);
#' otherwise it may not be valid. \code{x - a} shifts \code{x} by \code{-a}.
#' \item \code{a * x} or \code{x * a}: the product of a distribution and a number,
#' i.e. the distribution of \code{a} times the random variable (for \code{a < 0}
#' the distribution is reflected). Similarly, \code{a - x} is the distribution
#' of \code{a} minus the random variable.
#' }
#'
#' @param e1,e2 Two \code{distributionH} objects, or a \code{distributionH} and a
#' number.
#' @return A \code{distributionH} object.
#' @examples
#' d1 <- distributionH(x = c(1, 2, 3), p = c(0, 0.4, 1))
#' d2 <- distributionH(x = c(7, 8, 10, 15), p = c(0, 0.2, 0.7, 1))
#' d1 + d2
#' d2 + 10
#' 2 * d1
#' d2 - d1
#' @name Arith-distributionH
#' @rdname Arith-distributionH
#' @aliases +,distributionH,distributionH-method
#' @export
setMethod(
  "+",
  signature(e1 = "distributionH", e2 = "distributionH"),
  function(e1, e2) {
    if (!identical(e1@p, e2@p)) {
      tmp <- register(e1, e2)
      x <- callGeneric(tmp[[1]]@x, tmp[[2]]@x)
      e1@p <- tmp[[1]]@p
    } else {
      x <- callGeneric(e1@x, e2@x)
    }

    e1@x <- x

    e1@m <- e1@m + e2@m
    e1@s <- stdH(e1)
    return(e1)
    #            OBJ_NEW=new("distributionH",x,tmp[[1]]@p,(tmp[[1]]@m+tmp[[2]]@m))
  }
)
#' @rdname Arith-distributionH
#' @aliases +,numeric,distributionH-method
#' @usage NULL
#' @export

setMethod(
  "+",
  signature(e1 = "numeric", e2 = "distributionH"),
  function(e1, e2) {
    x <- callGeneric(rep(e1, length(e2@x)), e2@x)
    OBJ_NEW <- .distH(x, e2@p, (e1 + e2@m), e2@s)
  }
)
#' @rdname Arith-distributionH
#' @aliases +,distributionH,numeric-method
#' @usage NULL
#' @export
setMethod(
  "+",
  signature(e1 = "distributionH", e2 = "numeric"),
  function(e1, e2) {
    x <- callGeneric(rep(e2, length(e1@x)), e1@x)
    OBJ_NEW <- .distH(x, e1@p, (e2 + e1@m), e1@s)
  }
)
# Overloading of the difference of two distributions according to the L2 wasserstein --------
#' @rdname Arith-distributionH
#' @aliases -,distributionH,distributionH-method
#' @export
setMethod(
  "-",
  signature(e1 = "distributionH", e2 = "distributionH"),
  function(e1, e2) {
    if (!identical(e1@p, e2@p)) {
      tmp <- register(e1, e2)
      x <- callGeneric(tmp[[1]]@x, tmp[[2]]@x)
      e1@p <- tmp[[1]]@p
    } else {
      x <- callGeneric(e1@x, e2@x)
    }
    # (up to 1.0.9 x was recomputed here from tmp, that does not exist when the
    # two distributions share the same cdf: d1 - d2 failed for registered
    # distributions)
    OBJ_NEW <- new("distributionH", x, e1@p)
  }
)
#' @rdname Arith-distributionH
#' @aliases -,numeric,distributionH-method
#' @usage NULL
setMethod(
  "-",
  signature(e1 = "numeric", e2 = "distributionH"),
  function(e1, e2) {
    # the distribution of e1 - X: the quantile function is e1 - Q(1 - p)
    # (up to 1.0.9 it failed, x was not increasing)
    x <- rev(e1 - e2@x)
    OBJ_NEW <- .distH(x, rev(1 - e2@p), e1 - e2@m, e2@s)
  }
)
#' @rdname Arith-distributionH
#' @aliases -,distributionH,numeric-method
#' @usage NULL
setMethod(
  "-",
  signature(e1 = "distributionH", e2 = "numeric"),
  function(e1, e2) {
    x <- callGeneric(e1@x, rep(e2, length(e1@x)))
    OBJ_NEW <- .distH(x, e1@p, (e1@m - e2), e1@s)
  }
)

# Overloading of the product of a number by a distribution according to the L2 w --------
#' @rdname Arith-distributionH
#' @aliases *,distributionH,distributionH-method
#' @export
setMethod(
  "*",
  signature(e1 = "distributionH", e2 = "distributionH"),
  function(e1, e2) {
    stop("please use dotpW function product between distributions")
  }
)
#' @rdname Arith-distributionH
#' @aliases *,numeric,distributionH-method
#' @usage NULL
setMethod(
  "*",
  signature(e1 = "numeric", e2 = "distributionH"),
  function(e1, e2) {
    x <- callGeneric(rep(e1, length(e2@x)), e2@x)
    if (e1 < 0) {
      # the distribution of e1 * X, with e1 < 0 (up to 1.0.9 the result had a
      # decreasing x, i.e. it was not a valid distribution)
      x <- rev(x)
      e2@p <- rev(1 - e2@p)
    }
    e2@x <- x
    e2@m <- e1 * e2@m
    e2@s <- abs(e1) * e2@s
    return(e2)
  }
)
#' @rdname Arith-distributionH
#' @aliases *,distributionH,numeric-method
#' @usage NULL
setMethod(
  "*",
  signature(e1 = "distributionH", e2 = "numeric"),
  function(e1, e2) {
    x <- callGeneric(rep(e2, length(e1@x)), e1@x)
    if (e2 < 0) {
      x <- rev(x)
      e1@p <- rev(1 - e1@p)
    }
    e1@x <- x
    e1@m <- e2 * e1@m
    e1@s <- abs(e2) * e1@s
    return(e1)
  }
)

# Utilities for single or couples of distributionH --------------------------------
#' Fix the empty bins of histograms
#'
#' Checks for empty bins, i.e. bins with zero weight (two consecutive equal values
#' of the cdf). Each empty bin gets a weight of \code{1e-8} and the cdf is
#' renormalised. Since version 1.0.9 empty bins are fixed when a
#' \code{distributionH} object is created (see \code{\link{distributionH}}), so
#' this method is useful for objects created with previous versions of the
#' package or modified by assigning their slots directly. The \code{MatH} method
#' applies the check to all the cells of the matrix.
#'
#' @param object A \code{distributionH} or a \code{MatH} object.
#' @return An object of the same class of \code{object}, without empty bins.
#' @examples
#' ## an empty bin in an object modified by slot assignment
#' mydist <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.5, 0.6, 1))
#' mydist@p <- c(0, 0.5, 0.5, 1)
#' diff(checkEmptyBins(mydist)@p)
#' ## all the distributions of a MatH
#' newMat <- checkEmptyBins(China_Seas[1:5, 1:3])
#' @name checkEmptyBins
#' @rdname checkEmptyBins-methods
#' @exportMethod checkEmptyBins
setGeneric("checkEmptyBins", function(object) standardGeneric("checkEmptyBins"))
#' Quantile of a histogram
#'
#' Computes the quantile of a \code{distributionH} object at a given level,
#' i.e. the value \eqn{Q(p) = F^{-1}(p)} of its quantile function.
#'
#' @param object A \code{distributionH} object.
#' @param p A number between 0 and 1.
#' @return A number, the quantile of order \code{p}.
#' @seealso \code{\link{compP}}
#' @examples
#' mydist <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' compQ(mydist, 0.5) # the median
#' compQ(mydist, 0) # the minimum
#' compQ(mydist, 1) # the maximum
#' compQ(mydist, 0.25) # the first quartile
#' @name compQ
#' @rdname compQ-methods
#' @exportMethod compQ
setGeneric("compQ", function(object, p) standardGeneric("compQ"))
#' Cumulative distribution function of a histogram at a value
#'
#' Computes the value \eqn{F(q)} of the cumulative distribution function of a
#' \code{distributionH} object.
#'
#' @param object A \code{distributionH} object.
#' @param q A number.
#' @return A number between 0 and 1.
#' @seealso \code{\link{compQ}}
#' @examples
#' mydist <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' compP(mydist, 5)
#' @name compP
#' @rdname compP-methods
#' @exportMethod compP
setGeneric("compP", function(object, q) standardGeneric("compP"))
#' Register two histograms on a common cdf
#'
#' Given two \code{distributionH} objects, returns two equivalent
#' distributions sharing the same cdf values (the union of the two \code{p}
#' vectors). This is the basis for computing the statistics based on the L2
#' Wasserstein distance.
#'
#' @param object1,object2 Two \code{distributionH} objects.
#' @return A list: the two \code{distributionH} objects, with the same \code{p}
#' slot, and a matrix with their quantiles and the common cdf.
#' @references Irpino, A., Lechevallier, Y. and Verde, R. (2006): \emph{Dynamic
#' clustering of histograms using Wasserstein metric} In: Rizzi, A., Vichi, M.
#' (eds.) COMPSTAT 2006. Physica-Verlag, Berlin, 869-876.
#' @examples
#' mydist1 <- distributionH(c(1, 2, 3), c(0, 0.4, 1))
#' mydist2 <- distributionH(c(7, 8, 10, 15), c(0, 0.2, 0.7, 1))
#' regDist <- register(mydist1, mydist2)
#' regDist[[1]]@p
#' regDist[[1]]@x
#' regDist[[2]]@x
#' @name register
#' @rdname register-methods
#' @exportMethod register
setGeneric("register", function(object1, object2) standardGeneric("register"))

#' @rdname register-methods
#' @usage NULL
setMethod(
  f = "register", signature = c(object1 = "distributionH", object2 = "distributionH"),
  function(object1, object2) {
    res <- REGISTER2(object1, object2)
    return(res)
    # DIG=14
    # #if (!identical(object1@p,object2@p)){
    #   #object1@p=round(object1@p,digits = DIG)
    #   if (min(
    #       object1@p[2:length(object1@p)]-object1@p[1:(length(object1@p)-1)]
    #       )<((1e-1)^DIG)
    #       ){
    #     diffs=object1@p[2:length(object1@p)]-object1@p[1:(length(object1@p)-1)]
    #     diffs[which(diffs<((1e-1)^DIG))]=(1e-1)^DIG
    #     p=c(0,cumsum(diffs))/sum(diffs)
    #     object1@p=p
    #   }
    #   #object2@p=round(object2@p,digits = DIG)
    #   if (
    #     min(
    #       object2@p[2:length(object2@p)]-object2@p[1:(length(object2@p)-1)]
    #     )<((1e-1)^DIG)
    #
    #   ){
    #    # browser()
    #     diffs=object2@p[2:length(object2@p)]-object2@p[1:(length(object2@p)-1)]
    #     diffs[which(diffs<((1e-1)^DIG))]=(1e-1)^DIG
    #     p=c(0,cumsum(diffs))/sum(diffs)
    #     object2@p=p
    #   }
    #
    #   commoncdf=sort(unique(c(object1@p,object2@p)))
    #   nr=length(commoncdf)
    #   commoncdf[1]=0
    #   commoncdf[nr]=1
    #
    #   tmp00=setdiff(commoncdf,object1@p)
    #   if (length(tmp00)>0){
    #     tmp1=compQ_vect(object1,vp=tmp00)
    #     object1@x=sort(c(object1@x,tmp1))
    #   }
    #   tmp01=setdiff(commoncdf,object2@p)
    #   if (length(tmp01)>0){
    #     tmp2=compQ_vect(object2,vp=tmp01)
    #   object2@x=sort(c(object2@x,tmp2))
    #   }
    #  # if(length(object1@x)!=length(object2@x)){browser()}
    #   object1@p=commoncdf
    #   object2@p=commoncdf
    #
    #   return(c(object1,object2))
    # # }
    # # else{
    # #   return(c(object1,object2))
    # # }
  }
)
#' @rdname checkEmptyBins-methods
#' @usage NULL
setMethod(
  f = "checkEmptyBins", signature = "distributionH",
  function(object) {
    fx <- .fix_empty_bins_p(object@p)
    if (fx$fixed) {
      object@p <- fx$p
      object@m <- meanH(object)
      object@s <- stdH(object)
    }
    return(object)
  }
)
#' @rdname checkEmptyBins-methods
#' @usage NULL
setMethod(
  f = "checkEmptyBins", signature = "MatH",
  function(object) {
    object@M[] <- lapply(object@M, function(d) {
      if (length(d@x) > 0) checkEmptyBins(d) else d
    })
    return(object)
  }
)
#' @rdname compQ-methods
#' @usage NULL
setMethod(
  f = "compQ", signature = c(object = "distributionH", p = "numeric"),
  function(object, p) {
    # %Computes the p-th quantile p=[0,1] of the distribution o1
    # %INPUT  - p  a value in [0,1]
    # %           - o1 a distribution
    # %OUTPUT-  res the computed quantile
    # %example:q=compQ(o1,0.5) returns the median of the
    # % distribution o1

    # Check for errors
    if (p < 0 || p > 1) stop("p must be a value between 0 and 1")

    if (p <= 0) {
      return(q = object@x[1])
    }
    if (p >= 1) {
      return(q = object@x[length(object@x)])
    }

    # last position with cdf <= p (up to 1.0.8 which.max(x[p<=p]) was used, that
    # returned a wrong position when x has repeated values)
    pos1 <- max(which(object@p <= p))
    ini <- object@x[pos1]
    fin <- object@x[pos1 + 1]
    if (ini == fin) {
      q <- ini
    }
    else {
      q <- ini + (fin - ini) * (p - object@p[pos1]) / (object@p[pos1 + 1] - object@p[pos1])
    }
    return(q)
  }
)
#' @rdname compP-methods
#' @usage NULL
setMethod(
  f = "compP", signature = c(object = "distributionH", q = "numeric"),
  function(object, q) {
    domain <- object@x
    cumul <- object@p
    if (q <= domain[1]) {
      return(p = 0)
    }
    if (q >= domain[length(domain)]) {
      return(p = 1)
    }

    # last position with domain <= q (up to 1.0.8 which.max(cumul[domain<=q])
    # was used, that returned a wrong position in presence of empty bins)
    pos1 <- max(which(domain <= q))
    ini <- cumul[pos1]
    pos2 <- pos1 + 1
    fin <- cumul[pos2]

    if (ini == fin) {
      p <- fin
    }
    else {
      p <- ini + (cumul[pos2] - cumul[pos1]) * (q - domain[pos1]) / (domain[pos2] - domain[pos1])
    }
    return(as.numeric(p))
  }
)
# L2 Wasserstein distance between two distributions and related results ----
#' Squared L2 Wasserstein distance between two histograms
#'
#' Computes the squared L2 Wasserstein distance between two \code{distributionH}
#' objects, i.e. the integral of the squared difference of their quantile
#' functions. Optionally it returns its decomposition into a position, a size and
#' a shape component.
#'
#' @param object1,object2 Two \code{distributionH} objects.
#' @param details A logical value (default \code{FALSE}). If \code{TRUE}, the
#' decomposition of the distance is returned.
#' @return If \code{details = FALSE}, a number, the squared distance. If
#' \code{details = TRUE}, a named vector with the squared distance
#' (\code{SQ_W_dist}), its position (\code{POSITION}), size (\code{SIZE}) and
#' shape (\code{SHAPE}) components, and the correlation between the quantile
#' functions (\code{rQQ}).
#' @references
#' Irpino, A. and Romano, E. (2007): \emph{Optimal histogram representation of large data sets:
#' Fisher vs piecewise linear approximations}. RNTI E-9, 99-110.
#'
#' Irpino, A., Verde, R. (2015) \emph{Basic statistics for distributional
#' symbolic variables: a new metric-based approach}, Advances in Data Analysis
#' and Classification, 9(2), 143-175. \doi{10.1007/s11634-014-0176-4}
#' @seealso \code{\link{WH_MAT_DIST}} for the distance matrix of a \code{MatH}.
#' @examples
#' mydist1 <- distributionH(x = c(1, 2, 3), p = c(0, 0.4, 1))
#' mydist2 <- distributionH(x = c(7, 8, 10, 15), p = c(0, 0.2, 0.7, 1))
#' WassSqDistH(mydist1, mydist2)
#' WassSqDistH(mydist1, mydist2, details = TRUE)
#' @name WassSqDistH
#' @rdname WassSqDistH-methods
#' @usage NULL
#' @exportMethod WassSqDistH
setGeneric("WassSqDistH", function(object1, object2, ...) standardGeneric("WassSqDistH")) # Wasserstein distance between two distributions
#' Correlation between the quantile functions of two histograms
#'
#' Computes the correlation between the quantile functions of two
#' \code{distributionH} objects.
#'
#' @param e1,e2 Two \code{distributionH} objects.
#' @return A number, the correlation between the quantile functions.
#' @references Irpino, A., Verde, R. (2015) \emph{Basic statistics for
#' distributional symbolic variables: a new metric-based approach}, Advances in
#' Data Analysis and Classification, 9(2), 143-175.
#' \doi{10.1007/s11634-014-0176-4}
#' @examples
#' mydist1 <- distributionH(x = c(1, 2, 3), p = c(0, 0.4, 1))
#' mydist2 <- distributionH(x = c(7, 8, 10, 15), p = c(0, 0.2, 0.7, 1))
#' rQQ(mydist1, mydist2) # 0.916894
#' @name rQQ
#' @rdname rQQ-methods
#' @exportMethod rQQ
setGeneric("rQQ", function(e1, e2) standardGeneric("rQQ")) # Quantile-Quantile correlation between twi distributions

#' @rdname WassSqDistH-methods
setMethod(
  f = "WassSqDistH", signature = c(object1 = "distributionH", object2 = "distributionH"),
  # Computes the L2 Wasserstein squared distance between two distributions
  # INPUT: object1 and object2 - two distributionH objects
  # OUTPUT: A list containing the distance and its decomposition in three parts (position, size and shape)
  function(object1 = object1, object2 = object2, details = FALSE) {
    cp <- c_dotpW(object1, object2)
    x2 <- object1@s^2 + object1@m^2
    y2 <- object2@s^2 + object2@m^2
    D <- x2 + y2 - 2 * cp


    if (details) {
      DC <- (object1@m - object2@m)^2
      DS <- (object1@s - object2@s)^2

      DR <- ifelse(abs(D - DC - DS) < 1e-30, yes = 0, no = abs(D - DC - DS))
      rho <- (cp - object1@m * object2@m) / (object1@s * object2@s)
      if (rho < 0) rho <- 0
      resu <- c(D, DC, DS, DR, rho)
      names(resu) <- c("SQ_W_dist", "POSITION", "SIZE", "SHAPE", "rQQ")
      return(resu)
    }
    else {
      return(as.numeric(D))
    }
  }
)

#' Dot product of two histograms
#'
#' Computes the dot product, inducing the L2 Wasserstein metric, between two
#' \code{distributionH} objects, or between a \code{distributionH} and a number
#' (considered as an impulse distribution).
#'
#' @param e1,e2 Two \code{distributionH} objects, or a \code{distributionH} and a
#' number.
#' @return A number.
#' @references Irpino, A., Verde, R. (2015) \emph{Basic statistics for
#' distributional symbolic variables: a new metric-based approach}, Advances in
#' Data Analysis and Classification, 9(2), 143-175.
#' \doi{10.1007/s11634-014-0176-4}
#' @examples
#' mydist1 <- distributionH(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' mydist2 <- distributionH(x = c(5, 7, 15), p = c(0, 0.7, 1))
#' dotpW(mydist1, mydist2) # 39.51429
#' dotpW(mydist1, 3) # 13.2
#' dotpW(3, mydist1) # 13.2
#' @name dotpW
#' @rdname dotpW-methods
#' @exportMethod dotpW
setGeneric("dotpW", function(e1, e2) standardGeneric("dotpW")) # dotproduct from L2 Wasserstein
#' @rdname dotpW-methods
#' @aliases dotpW,distributionH,distributionH-method
#' @usage NULL
setMethod("dotpW",
  signature(e1 = "distributionH", e2 = "distributionH"),
  definition = function(e1, e2) {
    return(c_dotpW(e1, e2))
  }
)
#' @rdname dotpW-methods
#' @aliases dotpW,numeric,distributionH-method
#' @usage NULL
setMethod(
  "dotpW",
  signature(e1 = "numeric", e2 = "distributionH"),
  function(e1, e2) {
    dprod <- e1 * e2@m
    return(dprod)
  }
)
#' @rdname dotpW-methods
#' @aliases dotpW,distributionH,numeric-method
#' @usage NULL
setMethod(
  "dotpW",
  signature(e1 = "distributionH", e2 = "numeric"),
  function(e1, e2) {
    dprod <- e2 * e1@m
    return(dprod)
  }
)
#' @rdname rQQ-methods
#' @usage NULL
setMethod("rQQ",
  signature(e1 = "distributionH", e2 = "distributionH"),
  definition = function(e1, e2) {
    rQQ <- (dotpW(e1, e2) - e1@m * e2@m) / (e1@s * e2@s)
    return(rQQ)
  }
)
## ---- Show overridding for distributionH and MatH ----
#' Print a histogram
#'
#' Prints a \code{distributionH} object as a table of bins and probabilities.
#' If the histogram has many bins, only the first and the last five are shown.
#'
#' @param object A \code{distributionH} object.
#' @return Invisibly \code{NULL}; called for its side effect.
#' @examples
#' mydist <- distributionH(x = c(7, 8, 10, 15), p = c(0, 0.2, 0.7, 1))
#' mydist
#' @name show-distributionH
#' @rdname show-distributionH-methods
#' @aliases show,distributionH-method
setMethod("show",
  signature(object = "distributionH"),
  definition = function(object) {
    if (length(object@p) > 1) {
      mymat <- matrix(0, length(object@p) - 1, 2)
      if (length(object@p) > 11) {
        cat("Output shows the first five and the last five bins due to eccesive length \n")
        mymat <- matrix(0, 12, 2)
        mymat[1, 1] <- "X"
        mymat[1, 2] <- "p"
        count <- 0
        for (i in 2:6) {
          count <- count + 1
          mymat[count + 1, 1] <- paste("[", format(object@x[(i - 1)], digits = 5), "-", format(object@x[i], digits = 5), ")", sep = "")
          mymat[count + 1, 2] <- paste(format(object@p[i] - object@p[i - 1], digits = 4))
        }
        count <- count + 1
        mymat[count + 1, 1] <- paste("...")
        mymat[count + 1, 2] <- paste("...")
        for (i in (length(object@p) - 4):length(object@p)) {
          count <- count + 1
          mymat[count + 1, 1] <- paste("[", format(object@x[(i - 1)], digits = 5), " ; ", format(object@x[i], digits = 5), ")", sep = "")
          mymat[count + 1, 2] <- paste(format(object@p[i] - object@p[i - 1], digits = 4))
        }
        rownames(mymat) <- c(
          "", paste("Bin", 1:5, sep = "_"), "...",
          paste("Bin", (length(object@p) - 5):(length(object@p) - 1), sep = "_")
        )
        write.table(format(mymat, justify = "right", digits = 5),
          row.names = T, col.names = F, quote = F
        )
      }
      else {
        mymat <- matrix(0, length(object@p), 2)
        mymat[1, 1] <- "X"
        mymat[1, 2] <- "p"
        if (length(object@p) > 2) {
          for (i in 2:(length(object@p) - 1)) {
            mymat[i, 1] <- paste("[ ", format(object@x[(i - 1)], digits = 5), " ; ", format(object@x[i], digits = 5), " )", sep = "")
            mymat[i, 2] <- paste(format(object@p[i] - object@p[i - 1], digits = 4))
          }
        }
        mymat[length(object@p), 1] <- paste("[ ", format(object@x[(length(object@p) - 1)], digits = 5),
          " ; ", format(object@x[length(object@p)], digits = 5), " ]",
          sep = ""
        )
        mymat[length(object@p), 2] <- paste(format(object@p[length(object@p)] - object@p[length(object@p) - 1], digits = 4))

        rownames(mymat) <- c(" ", paste("Bin", 1:(length(object@p) - 1), sep = "_"))
        write.table(format(mymat, justify = "right"),
          row.names = T, col.names = F, quote = F
        )
      }
      cat(paste("\n mean = ", object@m, "  std  = ", object@s, "\n "))
    } else {
      (cat("Empty distributionH\n"))
    }
  }
)




## --- Plot  for distributionH  ----
#' Plot a histogram
#'
#' Plots a \code{distributionH} object as a histogram, a cumulative distribution
#' function, a quantile function, a density or a boxplot.
#'
#' @param x A \code{distributionH} object.
#' @param type The type of plot: \code{"HISTO"} (default, histogram), \code{"CDF"}
#' (cumulative distribution function), \code{"QF"} (quantile function),
#' \code{"DENS"} (kernel density approximation), \code{"HBOXPLOT"} (horizontal
#' boxplot) or \code{"VBOXPLOT"} (vertical boxplot).
#' @param col The fill colour (default \code{"green"}).
#' @param border The colour of the borders (default \code{"black"}).
#' @return The plot is drawn on the current device and the \code{ggplot} object is
#' returned invisibly, so that it can be further modified, e.g.
#' \code{p <- plot(mydist); p + ggplot2::theme_bw()}.
#' @examples
#' mydist <- distributionH(x = c(7, 8, 10, 15), p = c(0, 0.2, 0.7, 1))
#' plot(mydist)
#' plot(mydist, type = "HISTO", col = "red", border = "blue")
#' plot(mydist, type = "DENS", col = "red", border = "blue")
#' plot(mydist, type = "HBOXPLOT")
#' plot(mydist, type = "VBOXPLOT")
#' plot(mydist, type = "CDF")
#' plot(mydist, type = "QF")
#' @name plot-distributionH
#' @aliases plot,distributionH-method
#' @importFrom utils write.table
#' @export
setMethod(
  "plot",
  signature(x = "distributionH"),
  function(x, type = "HISTO", col = "green", border = "black") {
    .draw(.plot_distributionH(x, type = type, col = col, border = border))
  }
)
