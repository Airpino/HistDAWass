#' @rdname MatH-class
#' @param x A list of \code{distributionH} objects (default \code{NULL}, an empty
#' matrix).
#' @param nrows,ncols Integers, the number of rows (individuals) and of columns
#' (variables). Default 1.
#' @param rownames,varnames Character vectors (or lists of strings), the names
#' of the rows and of the columns (variables). Default \code{NULL}.
#' @param by.row A logical value (default \code{FALSE}): if \code{TRUE} the matrix
#' is filled by rows, otherwise by columns.
#' @return \code{MatH()} returns a \code{MatH} object.
#' @export
MatH <- function(x = NULL, nrows = 1, ncols = 1, rownames = NULL, varnames = NULL, by.row = FALSE) {
  MAT <- new("MatH",
    nrows = nrows,
    ncols = ncols,
    ListOfDist = x,
    names.rows = rownames,
    names.cols = varnames, by.row = by.row
  )
  return(MAT)
}

# overriding of "[" operator for MatH object ----
#' Extract a sub-matrix of histograms
#'
#' Extracts rows and/or columns of a \code{MatH} object, as for an ordinary
#' matrix.
#'
#' @param x A \code{MatH} object.
#' @param i,j Row and column indices (integers, names or logical vectors).
#' @param ... Not used.
#' @param drop Not used (the result is always a \code{MatH} object).
#' @return A \code{MatH} object.
#' @examples
#' SUB_D <- BLOOD[c(1, 2, 5), c(1, 2)]
#' SUB_D
#' @name extract-MatH
#' @rdname extract-methods
#' @aliases [,MatH,ANY,ANY,ANY-method [,MatH-method
#' @importFrom stats variable.names
#' @export
setMethod(
  "[",
  signature(x = "MatH"),
  function(x, i, j, ..., drop = TRUE) {
    if (missing(i) && missing(j)) {
      i <- c(1:nrow(x@M))
      j <- c(1:ncol(x@M))
    }
    else {
      if (missing(i)) i <- c(1:nrow(x@M))
      if (missing(j)) j <- c(1:ncol(x@M))
    }
    # consider negative indexes!TO BE DONE!!
    if (min(i) <= 0 | min(j) <= 0) {
      stop("negative indexes are not allowed in subsetting [,] a MatH object")
    }


    x@M <- matrix(x@M[i, j],
      nrow = length(i), ncol = length(j),
      dimnames = list(row.names(x@M)[i], colnames(x@M)[j])
    )

    return(x)
  }
)
# methods for getting information from a MatH

setGeneric("get.MatH.nrows", function(object) standardGeneric("get.MatH.nrows"))
#' Dimensions and names of a matrix of histograms
#'
#' Functions returning the number of rows and of columns, the names of the rows
#' (individuals) and of the columns (variables) of a \code{MatH} object.
#'
#' @param object A \code{MatH} object.
#' @return \code{get.MatH.nrows} and \code{get.MatH.ncols} return an integer;
#' \code{get.MatH.rownames} and \code{get.MatH.varnames} return a character
#' vector; \code{get.MatH.main.info} returns a list with all of them
#' (\code{nrows}, \code{ncols}, \code{rownames}, \code{varnames}).
#' @examples
#' get.MatH.nrows(BLOOD)
#' get.MatH.ncols(BLOOD)
#' get.MatH.rownames(BLOOD)
#' get.MatH.varnames(BLOOD)
#' get.MatH.main.info(BLOOD)
#' @name get.MatH.main.info
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.nrows get.MatH.nrows,MatH-method
#' @usage get.MatH.nrows(object)
#' @exportMethod get.MatH.nrows

setMethod(
  f = "get.MatH.nrows", signature = c(object = "MatH"),
  function(object) {
    return(nrow(object@M))
  }
)
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.ncols
#' @exportMethod get.MatH.ncols

setGeneric("get.MatH.ncols", function(object) standardGeneric("get.MatH.ncols"))
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.ncols,MatH-method
#' @usage NULL
setMethod(
  f = "get.MatH.ncols", signature = c(object = "MatH"),
  function(object) {
    return(ncol(object@M))
  }
)
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.rownames
#' @exportMethod get.MatH.rownames
setGeneric("get.MatH.rownames", function(object) standardGeneric("get.MatH.rownames"))
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.rownames,MatH-method
#' @usage NULL
setMethod(
  f = "get.MatH.rownames", signature = c(object = "MatH"),
  function(object) {
    return(rownames(object@M))
  }
)
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.varnames
#' @exportMethod get.MatH.varnames
setGeneric("get.MatH.varnames", function(object) standardGeneric("get.MatH.varnames"))
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.varnames,MatH-method
#' @usage NULL
setMethod(
  f = "get.MatH.varnames", signature = c(object = "MatH"),
  function(object) {
    return(colnames(object@M))
  }
)
#' @rdname get.MatH.main.info-methods
#' @exportMethod get.MatH.main.info
setGeneric("get.MatH.main.info", function(object) standardGeneric("get.MatH.main.info"))
#' @rdname get.MatH.main.info-methods
#' @aliases get.MatH.main.info,MatH-method
#' @usage NULL
setMethod(
  f = "get.MatH.main.info", signature = c(object = "MatH"),
  function(object) {
    return(list(
      nrows = get.MatH.nrows(object), ncols = get.MatH.ncols(object),
      rownames = get.MatH.rownames(object), varnames = get.MatH.varnames(object)
    ))
  }
)
#' Statistics of the histograms of a matrix
#'
#' Computes a statistic (mean, median, standard deviation, ...) for each
#' histogram of a \code{MatH} object.
#'
#' @param object A \code{MatH} object.
#' @param stat The statistic: \code{"mean"} (default), \code{"median"},
#' \code{"min"}, \code{"max"}, \code{"std"}, \code{"skewness"},
#' \code{"kurtosis"} or \code{"quantile"}.
#' @param prob A number between 0 and 1, the level of the quantile when
#' \code{stat = "quantile"}.
#' @return A list with the name of the statistic (\code{stat}), the level
#' \code{prob} (only when \code{stat = "quantile"}) and the matrix of values
#' (\code{mat}).
#' @examples
#' get.MatH.stats(BLOOD) # means
#' get.MatH.stats(BLOOD, stat = "median")
#' get.MatH.stats(BLOOD, stat = "quantile", prob = 0.05) # 5th percentiles
#' get.MatH.stats(BLOOD, stat = "std")
#' @name get.MatH.stats
#' @rdname get.MatH.stats-methods
#' @usage NULL
#' @exportMethod get.MatH.stats
setGeneric("get.MatH.stats", function(object, ...) standardGeneric("get.MatH.stats"))

#' @rdname get.MatH.stats-methods
#' @aliases get.MatH.stats,MatH-method
setMethod(
  f = "get.MatH.stats", signature = c(object = "MatH"),
  function(object, stat = "mean", prob = 0.5) {
    r <- get.MatH.nrows(object)
    c <- get.MatH.ncols(object)
    MAT <- matrix(NA, get.MatH.nrows(object), get.MatH.ncols(object))
    rownames(MAT) <- get.MatH.rownames(object)
    colnames(MAT) <- get.MatH.varnames(object)
    for (i in 1:r) {
      for (j in 1:c) {
        if (length(object@M[i, j][[1]]@x) > 0) {
          if (stat == "mean") {
            MAT[i, j] <- object@M[i, j][[1]]@m
          }
          if (stat == "std") {
            MAT[i, j] <- object@M[i, j][[1]]@s
          }
          if (stat == "skewness") {
            MAT[i, j] <- skewH(object@M[i, j][[1]])
          }
          if (stat == "kurtosis") {
            MAT[i, j] <- kurtH(object@M[i, j][[1]])
          }
          if (stat == "median") {
            MAT[i, j] <- compQ(object = object@M[i, j][[1]], p = 0.5)
          }
          if (stat == "quantile") {
            MAT[i, j] <- compQ(object = object@M[i, j][[1]], p = prob)
          }
          if (stat == "min") {
            MAT[i, j] <- compQ(object = object@M[i, j][[1]], p = 0)
          }
          if (stat == "max") {
            MAT[i, j] <- compQ(object = object@M[i, j][[1]], p = 1)
          }
        }
      }
    }
    if (stat == "quantile") {
      return(list(stat = stat, prob = prob, mat = MAT))
    } else {
      return(list(stat = stat, mat = MAT))
    }
  }
)




# methods for collating by row or by column two MatHs ----
#' Bind matrices of histograms
#'
#' \code{WH.bind.row} binds by rows two \code{MatH} objects with the same number
#' of columns, \code{WH.bind.col} binds by columns two \code{MatH} objects with
#' the same number of rows. \code{WH.bind} does one of the two according to
#' \code{byrow}.
#'
#' @param object1,object2 Two \code{MatH} objects.
#' @param byrow A logical value: if \code{TRUE} the objects are bound by rows,
#' otherwise by columns.
#' @return A \code{MatH} object.
#' @examples
#' WH.bind.row(BLOOD[1:3, ], BLOOD[5:8, ])
#' WH.bind.col(BLOOD[1:10, 1], BLOOD[1:10, 3])
#' WH.bind(BLOOD[1:10, 1], BLOOD[1:10, 3], byrow = FALSE)
#' @name WH.bind
#' @rdname WH.bind-methods
#' @aliases WH.bind.row
#' @exportMethod WH.bind.row
setGeneric("WH.bind.row", function(object1, object2) standardGeneric("WH.bind.row")) #
#' @rdname WH.bind-methods
#' @aliases WH.bind.col
#' @exportMethod WH.bind.col
setGeneric("WH.bind.col", function(object1, object2) standardGeneric("WH.bind.col")) #
#' @rdname WH.bind-methods
#' @exportMethod WH.bind
setGeneric("WH.bind", function(object1, object2, byrow) standardGeneric("WH.bind")) #
#' @rdname WH.bind-methods
#' @aliases WH.bind.row,MatH-method
#' @usage NULL
setMethod(
  f = "WH.bind.row", signature = c(object1 = "MatH", object2 = "MatH"),
  function(object1, object2) {
    ncol1 <- ncol(object1@M)
    ncol2 <- ncol(object2@M)
    nrow1 <- nrow(object1@M)
    nrow2 <- nrow(object2@M)
    if (ncol1 != ncol2) {
      stop("The two matrix must have the same number of columns")
    }
    object1@M <- rbind(object1@M, object2@M)
    return(object1)
  }
)
#' @rdname WH.bind-methods
#' @aliases WH.bind.col,MatH-method
#' @usage NULL
setMethod(
  f = "WH.bind.col", signature = c(object1 = "MatH", object2 = "MatH"),
  function(object1, object2) {
    ncol1 <- ncol(object1@M)
    ncol2 <- ncol(object2@M)
    nrow1 <- nrow(object1@M)
    nrow2 <- nrow(object2@M)
    if (nrow1 != nrow2) {
      stop("The two matrix must have the same number of rows")
    }
    # NewMat=new("MatH", nrows=nrow1,ncols=ncol1+ncol2)
    object1@M <- cbind(object1@M, object2@M)
    return(object1)
  }
)
#' @rdname WH.bind-methods
#' @aliases WH.bind,MatH-method
#' @usage NULL
setMethod(
  f = "WH.bind", signature = c(object1 = "MatH", object2 = "MatH"),
  function(object1, object2, byrow = TRUE) {
    ncol1 <- ncol(object1@M)
    ncol2 <- ncol(object2@M)
    nrow1 <- nrow(object1@M)
    nrow2 <- nrow(object2@M)
    if (byrow == TRUE) {
      NewMat <- WH.bind.row(object1, object2)
    }
    else {
      NewMat <- WH.bind.col(object1, object2)
    }
    return(NewMat)
  }
)

# methods for MatH based on the L2 Wasserstein distance between distributions ----
#' Sum of two matrices of histograms
#'
#' Sums two \code{MatH} objects of the same dimensions, cell by cell, by summing
#' the quantile functions of the histograms (a sum consistent with the L2
#' Wasserstein metric).
#'
#' @param object1,object2 Two \code{MatH} objects with the same dimensions.
#' @return A \code{MatH} object.
#' @examples
#' MAT <- WH.mat.sum(BLOOD[1:5, ], BLOOD[6:10, ])
#' @name WH.mat.sum
#' @rdname WH.mat.sum-methods
#' @exportMethod WH.mat.sum
setGeneric("WH.mat.sum", function(object1, object2) standardGeneric("WH.mat.sum")) # ok matrix sum
#' Product of two matrices of histograms
#'
#' Matrix product of two \code{MatH} objects, where the product of two cells is
#' the dot product of the histograms inducing the L2 Wasserstein metric (see
#' \code{\link{dotpW}}).
#'
#' @param object1,object2 Two \code{MatH} objects with conformable dimensions.
#' @param traspose1,traspose2 Logical values (default \code{FALSE}): if
#' \code{TRUE} the corresponding matrix is transposed.
#' @return A numeric matrix.
#' @examples
#' MAT <- WH.mat.prod(BLOOD[1:5, ], BLOOD[6:10, ], traspose1 = TRUE)
#' @name WH.mat.prod
#' @rdname WH.mat.prod-methods
#' @usage NULL
#' @exportMethod WH.mat.prod
setGeneric("WH.mat.prod", function(object1, object2, ...) standardGeneric("WH.mat.prod")) # ok matrix product
#' @rdname WH.mat.sum-methods
#' @aliases WH.mat.sum,MatH-method
#' @usage NULL
setMethod(
  f = "WH.mat.sum", signature = c(object1 = "MatH", object2 = "MatH"),
  # sums two MatH, i.e.  two matrices of distributionsH
  # INPUT:
  # OUTPUT:

  function(object1, object2) {
    nrows1 <- nrow(object1@M)
    ncols1 <- ncol(object1@M)
    nrows2 <- nrow(object1@M)
    ncols2 <- ncol(object1@M)
    if (!identical(dim(object1@M), dim(object2@M))) {
      stop("the two matrices must be of the same dimension")
    }
    else {
      MATS <- object1
      TMP <- new("MatH", 1, 2)
      for (r in 1:nrows1) {
        for (c in 1:ncols1) {
          TMP@M[1, 1][[1]] <- object1@M[r, c][[1]]
          TMP@M[1, 2][[1]] <- object2@M[r, c][[1]]
          TMP <- registerMH(TMP)
          MATS@M[r, c][[1]] <- new(
            "distributionH",
            (TMP@M[1, 1][[1]]@x + TMP@M[1, 2][[1]]@x),
            TMP@M[1, 1][[1]]@p,
            (TMP@M[1, 1][[1]]@m + TMP@M[1, 2][[1]]@m)
          )
        }
      }
    }
    return(MATS)
  }
)
#' @rdname WH.mat.prod-methods
#' @aliases WH.mat.prod,MatH-method
setMethod(
  f = "WH.mat.prod", signature = c(object1 = "MatH", object2 = "MatH"),
  # sums two MatH, i.e.  two matrics of distributionsH
  # INPUT:
  # OUTPUT:

  function(object1, object2, traspose1 = FALSE, traspose2 = FALSE) {
    if (traspose1 == TRUE) {
      # trasposing the first matrix
      object1@M <- t(object1@M)
    }
    if (traspose2 == TRUE) {
      # trasposing the second matrix
      object2@M <- t(object2@M)
    }
    nrows1 <- nrow(object1@M)
    ncols1 <- ncol(object1@M)
    nrows2 <- nrow(object2@M)
    ncols2 <- ncol(object2@M)
    if (ncols1 != nrows2) {
      cat(
        "Fisrt matrix dimensions ", nrow(object1@M), "x", ncol(object1@M), "\n",
        "Second matrix dimensions ", nrow(object2@M), "x", ncol(object2@M), "\n"
      )
      stop("Dimensions of matrices are not compatible")
    }



    MAT <- matrix(0, nrows1, ncols2)
    #             cat("Fisrt matrix dimensions ", nrow(object1@M), "x", ncol(object1@M), "\n",
    #                 "Second matrix dimensions ", nrow(object2@M), "x", ncol(object2@M), "\n")
    for (r in 1:nrows1) {
      for (c in 1:ncols2) {
        for (els in 1:ncols1) {
          MAT[r, c] <- MAT[r, c] + dotpW(object1@M[r, els][[1]], object2@M[els, c][[1]])
        }
      }
    }
    return(MAT)
  }
)

# L2 Wasserstein basic operations and basic statistics for matrices of distributionH ----
#' Weighted sum and mean of the histograms of a matrix
#'
#' \code{WH.vec.sum} and \code{WH.vec.mean} compute the (weighted) sum and the
#' (weighted) mean of all the histograms of a \code{MatH} object, consistently
#' with the L2 Wasserstein metric: the result is the histogram whose quantile
#' function is the (weighted) sum or mean of the quantile functions. The mean is
#' the Wasserstein barycenter of the histograms.
#'
#' @param object A \code{MatH} object.
#' @param w A vector of non-negative weights, one for each cell of the matrix.
#' Default: equal weights.
#' @return A \code{distributionH} object.
#' @examples
#' WH.vec.mean(BLOOD)
#' WH.vec.sum(BLOOD)
#' # weighted mean with random weights
#' RN <- runif(get.MatH.nrows(BLOOD) * get.MatH.ncols(BLOOD))
#' WH.vec.mean(BLOOD, w = RN)
#' @name WH.vec.mean
#' @rdname WH.vec.mean-methods
#' @aliases WH.vec.sum
#' @usage NULL
#' @exportMethod WH.vec.sum
setGeneric("WH.vec.sum", function(object, ...) standardGeneric("WH.vec.sum")) # OK weighted sum of a vector of distributionH
#' @rdname WH.vec.mean-methods
#' @usage NULL
#' @exportMethod WH.vec.mean
setGeneric("WH.vec.mean", function(object, ...) standardGeneric("WH.vec.mean")) # OK weighted mean of a vector of distributionH
#' Sum of squares, covariance and correlation of histogram variables
#'
#' Compute the matrix of the sums of squared deviations from the mean
#' (\code{WH.SSQ}), the variance-covariance matrix (\code{WH.var.covar}) and the
#' correlation matrix (\code{WH.correlation}) of the variables of a \code{MatH}
#' object, according to the L2 Wasserstein metric. The functions with suffix 2
#' compute the same (rectangular) matrices between the variables of two
#' \code{MatH} objects with the same number of rows.
#'
#' @param object A \code{MatH} object.
#' @param object1,object2 Two \code{MatH} objects with the same number of rows.
#' @param w A vector of non-negative weights, one for each row (individual).
#' Default: equal weights.
#' @return A numeric matrix.
#' @references Irpino, A., Verde, R. (2015) \emph{Basic statistics for
#' distributional symbolic variables: a new metric-based approach}, Advances in
#' Data Analysis and Classification, 9(2), 143-175.
#' \doi{10.1007/s11634-014-0176-4}
#' @examples
#' WH.SSQ(BLOOD)
#' WH.var.covar(BLOOD)
#' WH.correlation(BLOOD)
#' # weighted, with random weights
#' RN <- runif(get.MatH.nrows(BLOOD))
#' WH.var.covar(BLOOD, w = RN)
#' # between two sets of variables
#' WH.correlation2(BLOOD[, 1], BLOOD[, 2:3])
#' @name WH.var.covar
#' @rdname WH.var.covar-methods
#' @aliases WH.SSQ
#' @usage NULL
#' @exportMethod WH.SSQ
setGeneric("WH.SSQ", function(object, ...) standardGeneric("WH.SSQ")) # weighted de-codeviance matrix
#' @rdname WH.var.covar-methods
#' @usage NULL
#' @exportMethod WH.var.covar
setGeneric("WH.var.covar", function(object, ...) standardGeneric("WH.var.covar")) # weighted variance variance matrix
#' @rdname WH.var.covar-methods
#' @aliases WH.correlation
#' @usage NULL
#' @exportMethod WH.correlation
setGeneric("WH.correlation", function(object, ...) standardGeneric("WH.correlation")) # weighted corelation matrix
#' @rdname WH.var.covar-methods
#' @aliases WH.SSQ2
#' @usage NULL
#' @exportMethod WH.SSQ2
setGeneric("WH.SSQ2", function(object1, object2, ...) standardGeneric("WH.SSQ2")) # weighted de-codeviance matrix
#' @rdname WH.var.covar-methods
#' @aliases WH.var.covar2
#' @usage NULL
#' @exportMethod WH.var.covar2
setGeneric("WH.var.covar2", function(object1, object2, ...) standardGeneric("WH.var.covar2")) # weighted variance variance matrix
#' @rdname WH.var.covar-methods
#' @aliases WH.correlation2
#' @usage NULL
#' @exportMethod WH.correlation2
setGeneric("WH.correlation2", function(object1, object2, ...) standardGeneric("WH.correlation2")) # weighted corelation matrix

#' @rdname WH.vec.mean-methods
#' @aliases WH.vec.sum,MatH-method
setMethod(
  f = "WH.vec.sum", signature = c(object = "MatH"),
  function(object, w = numeric(0)) {
    nrows <- nrow(object@M)
    ncols <- ncol(object@M)
    nelem <- nrows * ncols
    if (missing(w)) {
      w <- rep(1, nelem)
    }
    else {
      if (length(object@M) != length(w)) {
        stop("Wheights must have the same dimensions of the input matrix of distributions")
      }
      if (min(w) < 0) {
        stop("Weights must be positive!!")
      }
    }
    w <- matrix(w, nrows, ncols)

    SUM <- new("distributionH", c(0, 0), c(0, 1))
    for (c in 1:ncols) {
      for (r in 1:nrows) {
        SUM <- SUM + w[r, c] * object@M[r, c][[1]]
      }
    }
    return(SUM)
  }
)
#' @rdname WH.vec.mean-methods
#' @aliases WH.vec.mean,MatH-method
setMethod(
  f = "WH.vec.mean", signature = c(object = "MatH"),
  function(object, w = numeric(0)) {
    # if (length(object@M)==1) return(object)
    # WH MEAN H qua si puo migliorare -----
    nrows <- nrow(object@M)
    ncols <- ncol(object@M)
    nelem <- nrows * ncols
    if (missing(w)) {
      w <- rep(1 / nelem, nelem)
    }
    else {
      if (length(object@M) != length(w)) {
        stop("Wheights must have the same dimensions of the input matrix of distributions")
      }
      if (min(w) < 0) {
        stop("Weights must be positive!!")
      }
    }
    w <- matrix(w, nrows, ncols)
    w <- w / sum(w)
    if (ncols == 1) {
      MEAN <- MEAN_VA(object, w)
    }
    else {
      w2 <- colSums(w)
      w2 <- w2 / sum(w2)
      MEAN <- w2[1] * MEAN_VA(object[, 1], w[, 1])
      for (c in 2:ncols) {
        MEAN <- MEAN + w2[c] * MEAN_VA(object[, c], w[, c])
      }
    }
    return(MEAN)
  }
)
#' @rdname WH.var.covar-methods
#' @aliases WH.SSQ,MatH-method
setMethod(
  f = "WH.SSQ", signature = c(object = "MatH"),
  function(object, w = numeric(0)) {
    nrows <- nrow(object@M)
    ncols <- ncol(object@M)
    nelem <- nrows * ncols
    if (missing(w)) {
      w <- rep(1, nrows)
    }
    else {
      if (nrows != length(w)) {
        stop("Wheights must have the same length of rows of the input matrix of distributions")
      }
      if (min(w) < 0) {
        stop("Weights must be positive!!")
      }
    }
    w <- matrix(w, nrows, 1)
    DEV_MAT <- SSQ_RCPP(object, w)
    # w=w/sum(w)
    # DEV_MAT=matrix(0,ncols,ncols)
    colnames(DEV_MAT) <- colnames(object@M)
    rownames(DEV_MAT) <- colnames(object@M)
    # compute the means
    # MEANS=new("MatH",1,ncols)
    # for (v1 in 1:ncols){
    #   MEANS@M[1,v1][[1]]=WH.vec.mean(object[,v1],w)
    # }
    # for (v1 in 1:ncols){
    #   for (v2 in v1:ncols){
    #     for (indiv in 1:nrows){
    #       if (v1==v2){
    #         DEV_MAT[v1,v2]=DEV_MAT[v1,v2]+
    #           w[indiv,1]*((object@M[indiv,v1][[1]]@s)^2+(object@M[indiv,v1][[1]]@m)^2)
    #       }else{
    #         DEV_MAT[v1,v2]=DEV_MAT[v1,v2]+
    #           w[indiv,1]*dotpW(object@M[indiv,v1][[1]],object@M[indiv,v2][[1]])
    #       }
    #     }
    #     if (v2>v1){
    #       DEV_MAT[v1,v2]=DEV_MAT[v1,v2]-sum(w)*dotpW(MEANS@M[1,v1][[1]],MEANS@M[1,v2][[1]])
    #       DEV_MAT[v2,v1]=DEV_MAT[v1,v2]
    #     }else{
    #       DEV_MAT[v1,v1]=DEV_MAT[v1,v1]-sum(w)*(MEANS@M[1,v1][[1]]@s^2+MEANS@M[1,v1][[1]]@m^2)
    #     }
    #   }
    # }
    # if(ncols==1){
    #   return(as.vector(DEV_MAT))
    # }
    # else
    return(DEV_MAT)
  }
)
#' @rdname WH.var.covar-methods
#' @aliases WH.var.covar,MatH-method
setMethod(
  f = "WH.var.covar", signature = c(object = "MatH"),
  function(object, w = numeric(0)) {
    nrows <- nrow(object@M)
    ncols <- ncol(object@M)
    nelem <- nrows * ncols
    if (missing(w)) {
      w <- rep(1, nrows)
    }
    else {
      if (nrows != length(w)) {
        stop("Weights must have the same length of rows of the input matrix of distributions")
      }
      if (min(w) < 0) {
        stop("Weights must be positive!!")
      }
    }
    w <- matrix(w, nrows, 1)
    w <- w / sum(w)
    COV_MAT <- WH.SSQ(object, w)
    return(COV_MAT)
  }
)
#' @rdname WH.var.covar-methods
#' @aliases WH.correlation,MatH-method
setMethod(
  f = "WH.correlation", signature = c(object = "MatH"),
  function(object, w = numeric(0)) {
    nrows <- nrow(object@M)
    ncols <- ncol(object@M)
    nelem <- nrows * ncols
    if (missing(w)) {
      w <- rep(1, nrows)
    }
    else {
      if (nrows != length(w)) {
        stop("Wheights must have the same length of rows of the input matrix of distributions")
      }
      if (min(w) < 0) {
        stop("Weights must be positive!!")
      }
    }
    w <- matrix(w, nrows, 1)
    w <- w / sum(w)
    COV_MAT <- WH.var.covar(object, w)
    CORR_MAT <- as.matrix(COV_MAT)
    # browser()
    # a=Sys.time()
    CORR_MAT <- COV_MAT / (t(t(sqrt(diag(COV_MAT)))) %*% sqrt(diag(COV_MAT)))
    # b=Sys.time()
    # print(b-a)
    #
    #  for (v1 in 1:ncols){
    #    for (v2 in v1:ncols){
    #      CORR_MAT[v1,v2]= COV_MAT[v1,v2]/sqrt((COV_MAT[v1,v1]*COV_MAT[v2,v2]))
    #      CORR_MAT[v2,v1]=CORR_MAT[v1,v2]
    #    }
    #  }
    #  c=Sys.time()
    #  print(c-b)
    # #
    # browser()
    return(CORR_MAT)
  }
)
#' @rdname WH.var.covar-methods
#' @aliases WH.SSQ2,MatH-method
setMethod(
  f = "WH.SSQ2", signature = c(object1 = "MatH", object2 = "MatH"),
  function(object1, object2, w = numeric(0)) {
    nrows1 <- nrow(object1@M)
    ncols1 <- ncol(object1@M)

    nrows2 <- nrow(object2@M)
    ncols2 <- ncol(object2@M)

    if (nrows1 != nrows2) {
      stop("The two matrices have a different number of rows")
    }

    if (missing(w)) {
      w <- rep(1, nrows1)
    }
    else {
      if (nrows1 != length(w)) {
        stop("Wheights must have the same length of rows of the input matrix of distributions")
      }
      if (min(w) < 0) {
        stop("Weights must be positive!!")
      }
    }
    w <- matrix(w, nrows1, 1)
    # w=w/sum(w)
    DEV_MAT <- matrix(0, ncols1, ncols2)
    rownames(DEV_MAT) <- colnames(object1@M)
    colnames(DEV_MAT) <- colnames(object2@M)

    # compute the means
    MEANS1 <- new("MatH", 1, ncols1)
    for (v1 in 1:ncols1) {
      MEANS1@M[1, v1][[1]] <- WH.vec.mean(object1[, v1], w)
    }
    MEANS2 <- new("MatH", 1, ncols2)
    for (v2 in 1:ncols2) {
      MEANS2@M[1, v2][[1]] <- WH.vec.mean(object2[, v2], w)
    }


    for (v1 in 1:ncols1) {
      for (v2 in 1:ncols2) {
        for (indiv in 1:nrows1) {
          DEV_MAT[v1, v2] <- DEV_MAT[v1, v2] + w[indiv, 1] * dotpW(object1@M[indiv, v1][[1]], object2@M[indiv, v2][[1]])
        }

        DEV_MAT[v1, v2] <- DEV_MAT[v1, v2] - sum(w) * dotpW(MEANS1@M[1, v1][[1]], MEANS2@M[1, v2][[1]])
      }
    }
    if (ncols1 == 1 && ncols2 == 1) {
      return(as.vector(DEV_MAT))
    }
    else {
      return(DEV_MAT)
    }
  }
)
#' @rdname WH.var.covar-methods
#' @aliases WH.var.covar2,MatH-method
setMethod(
  f = "WH.var.covar2", signature = c(object1 = "MatH", object2 = "MatH"),
  function(object1, object2, w = numeric(0)) {
    nrows1 <- nrow(object1@M)
    ncols1 <- ncol(object1@M)

    nrows2 <- nrow(object2@M)
    ncols2 <- ncol(object2@M)

    if (nrows1 != nrows2) {
      stop("The two matrices have a different number of rows")
    }

    if (missing(w)) {
      w <- rep(1, nrows1)
    }
    else {
      if (nrows1 != length(w)) {
        stop("Wheights must have the same length of rows of the input matrix of distributions")
      }
      if (min(w) < 0) {
        stop("Weights must be positive!!")
      }
    }
    w <- matrix(w, nrows1, 1)
    w <- w / sum(w)
    VAR_MAT <- matrix(0, ncols1, ncols2)
    rownames(VAR_MAT) <- colnames(object1@M)
    colnames(VAR_MAT) <- colnames(object2@M)

    # compute the means
    MEANS1 <- new("MatH", 1, ncols1)
    for (v1 in 1:ncols1) {
      MEANS1@M[1, v1][[1]] <- WH.vec.mean(object1[, v1], w)
    }
    MEANS2 <- new("MatH", 1, ncols2)
    for (v2 in 1:ncols2) {
      MEANS2@M[1, v2][[1]] <- WH.vec.mean(object2[, v2], w)
    }


    for (v1 in 1:ncols1) {
      for (v2 in 1:ncols2) {
        for (indiv in 1:nrows1) {
          VAR_MAT[v1, v2] <- VAR_MAT[v1, v2] + w[indiv, 1] * dotpW(object1@M[indiv, v1][[1]], object2@M[indiv, v2][[1]])
        }

        VAR_MAT[v1, v2] <- VAR_MAT[v1, v2] - sum(w) * dotpW(MEANS1@M[1, v1][[1]], MEANS2@M[1, v2][[1]])
      }
    }
    if (ncols1 == 1 && ncols2 == 1) {
      return(as.vector(VAR_MAT))
    }
    else {
      return(VAR_MAT)
    }
  }
)
#' @rdname WH.var.covar-methods
#' @aliases WH.correlation2,MatH-method
setMethod(
  f = "WH.correlation2", signature = c(object1 = "MatH", object2 = "MatH"),
  function(object1, object2, w = numeric(0)) {
    nrows1 <- nrow(object1@M)
    ncols1 <- ncol(object1@M)

    nrows2 <- nrow(object2@M)
    ncols2 <- ncol(object2@M)

    if (nrows1 != nrows2) {
      stop("The two matrices have a different number of rows")
    }

    if (missing(w)) {
      w <- rep(1, nrows1)
    }
    else {
      if (nrows1 != length(w)) {
        stop("Wheights must have the same length of rows of the input matrix of distributions")
      }
      if (min(w) < 0) {
        stop("Weights must be positive!!")
      }
    }
    w <- matrix(w, nrows1, 1)
    w <- w / sum(w)
    COV_MAT <- WH.var.covar2(object1, object2, w)
    CORR_MAT <- as.matrix(COV_MAT)
    # qua perde tempo
    for (v1 in 1:ncols1) {
      for (v2 in 1:ncols2) {
        CORR_MAT[v1, v2] <- COV_MAT[v1, v2] / sqrt(WH.var.covar(object1[, v1], w) * WH.var.covar(object2[, v2], w))
      }
    }
    if (length(CORR_MAT) == 1) {
      return(as.vector(CORR_MAT))
    } else {
      return(CORR_MAT)
    }
  }
)

# Utility methods for registration of distributions ----
#' Register the histograms of a matrix on a common cdf
#'
#' \code{registerMH} recomputes all the histograms of a \code{MatH} object so that
#' they share the same cdf values (the same \code{p} slot), without changing the
#' distributions. This allows a fast computation of all the methods based on the
#' L2 Wasserstein metric. \code{is.registeredMH} checks whether the histograms
#' already share the same cdf.
#'
#' @param object A \code{MatH} object.
#' @return \code{registerMH} returns a \code{MatH} object whose histograms share
#' the same cdf; \code{is.registeredMH} returns \code{TRUE} or \code{FALSE}.
#' @references Irpino, A., Lechevallier, Y. and Verde, R. (2006): \emph{Dynamic
#' clustering of histograms using Wasserstein metric} In: Rizzi, A., Vichi, M.
#' (eds.) COMPSTAT 2006. Physica-Verlag, Berlin, 869-876.
#' @examples
#' mydist1 <- distributionH(c(1, 2, 3), c(0, 0.4, 1))
#' mydist2 <- distributionH(c(7, 8, 10, 15), c(0, 0.2, 0.7, 1))
#' mydist3 <- distributionH(c(9, 11, 20), c(0, 0.8, 1))
#' MyMAT <- MatH(x = list(mydist1, mydist2, mydist3), nrows = 1, ncols = 3)
#' is.registeredMH(MyMAT) # FALSE
#' MATregistered <- registerMH(MyMAT)
#' is.registeredMH(MATregistered) # TRUE
#' MATregistered@M[[1, 1]]@p
#' @name registerMH
#' @rdname registerMH-methods
#' @aliases is.registeredMH
#' @exportMethod is.registeredMH
setGeneric("is.registeredMH", function(object) standardGeneric("is.registeredMH")) # OK
#' @rdname registerMH-methods
#' @aliases is.registeredMH,MatH-method
#' @usage NULL
setMethod(
  f = "is.registeredMH", signature = c(object = "MatH"),
  # check if all the distributions share the same cdf
  # INPUT: object11  - a vector or a matrix two distributions
  # OUTPUT: resu - a matrix of distributionH objects with
  # recomputed quantiles on a common cdf
  function(object) {
    nrows <- nrow(object@M)
    ncols <- ncol(object@M)
    ndis <- nrows * ncols
    # Check if the distribution are registered
    OK <- 1
    count <- 1
    r <- 1
    tmpcdf <- object@M[1, 1][[1]]@p
    while (OK == 1) {
      count <- count + 1
      if (count <= ndis) {
        if (!identical(tmpcdf, object@M[count][[1]]@p)) {
          OK <- 0
          return(FALSE)
        }
      }
      else {
        OK <- 0
        return(TRUE)
      }
    }
  }
)
#' @rdname registerMH-methods
#' @exportMethod registerMH
setGeneric("registerMH", function(object) standardGeneric("registerMH")) # OK
#' @rdname registerMH-methods
#' @aliases registerMH,MatH-method
#' @usage NULL
setMethod(
  f = "registerMH", signature = c(object = "MatH"),
  # register a row or a column vector of qfs of distributionH:
  # if the cdf are different a a matrix resu is returned with the quantiles of the two
  # distribution computed at the same levels of a common vector of cdfs.
  # INPUT: object11  - a vector or a matrix two distributions
  # OUTPUT: resu - a matrix of distributionH objects with
  # recomputed quantiles on a common cdf
  function(object) {
    nrows <- nrow(object@M)
    ncols <- ncol(object@M)
    ndis <- nrows * ncols
    # Check if the distributions are registered
    if (is.registeredMH(object)) {
      return(object)
    }
    # all the cdf values (up to 1.0.8 collected with rbind() in a loop)
    commoncdf <- unlist(lapply(object@M, function(o) o@p), use.names = FALSE)
    commoncdf <- sort(unique(round(commoncdf, digits = 10)))
    commoncdf[1] <- 0
    commoncdf[length(commoncdf)] <- 1
    # check for tiny bins and for very long vectors of wheights

    # end of check
    nr <- length(commoncdf)
    result <- matrix(0, nr, (ndis + 1))
    result[, (ndis + 1)] <- commoncdf
    NEWMAT <- new("MatH", nrows, ncols)
    for (r in 1:nrows) {
      for (c in 1:ncols) {
        x <- compQ_vect(object@M[r, c][[1]], vp = commoncdf)
        # x=numeric(0)
        # for (rr in 1:nr){
        #   x=c(x,compQ(object@M[r,c][[1]],commoncdf[rr]))
        # }
        NEWMAT@M[r, c][[1]] <- .distH(x, commoncdf)
      }
    }
    return(NEWMAT)
  }
)




#' Centre the histograms of a matrix
#'
#' Shifts each histogram of a \code{MatH} object so that its mean is zero.
#'
#' @param object A \code{MatH} object.
#' @return A \code{MatH} object whose histograms have zero mean.
#' @examples
#' CEN_BLOOD <- Center.cell.MatH(BLOOD)
#' get.MatH.stats(CEN_BLOOD, stat = "mean")
#' @name Center.cell.MatH
#' @rdname Center.cell.MatH-methods
#' @exportMethod Center.cell.MatH
setGeneric("Center.cell.MatH", function(object) standardGeneric("Center.cell.MatH")) # OK
#' @rdname Center.cell.MatH-methods
#' @aliases Center.cell.MatH,MatH-method
#' @usage NULL
setMethod(
  f = "Center.cell.MatH", signature = c(object = "MatH"),
  function(object) {
    nr <- get.MatH.nrows(object)
    nc <- get.MatH.ncols(object)
    NM <- object
    for (i in 1:nr) {
      for (j in 1:nc) {
        NM@M[i, j][[1]]@x <- NM@M[i, j][[1]]@x - NM@M[i, j][[1]]@m
        NM@M[i, j][[1]]@m <- 0
      }
    }
    return(NM)
  }
)
## Show overridding ----
#' Print a matrix of histograms
#'
#' Prints a \code{MatH} object showing, for each histogram, its mean and its
#' standard deviation.
#'
#' @param object A \code{MatH} object.
#' @return Invisibly \code{NULL}; called for its side effect.
#' @examples
#' BLOOD
#' @name show-MatH
#' @rdname show-MatH-methods
#' @aliases show,MatH-method
setMethod("show",
  signature(object = "MatH"),
  definition = function(object) {
    cat("a matrix of distributions \n", paste(
      ncol(object@M), " variables ",
      nrow(object@M), " rows \n"
    ), "each distibution in the cell is represented by the mean and the standard deviation \n ")
    mymat <- matrix(0, nrow(object@M) + 1, ncol(object@M))
    for (i in 1:ncol(object@M)) {
      mymat[1, i] <- colnames(object@M)[i]
    }
    for (i in 1:nrow(object@M)) {
      for (j in 1:ncol(object@M)) {
        if (length(object@M[i, j][[1]]@x) == 0) {
          mymat[i + 1, j] <- paste("Empty distribution")
        }
        else {
          if ((abs(object@M[i, j][[1]]@m) > 1e5 || abs(object@M[i, j][[1]]@m) < 1e-5) &&
            (object@M[i, j][[1]]@s > 1e5 || object@M[i, j][[1]]@s < 1e-5)) {
            mymat[i + 1, j] <- paste(
              "[m=", format(object@M[i, j][[1]]@m, digits = 5, scientific = TRUE),
              " ,s=", format(object@M[i, j][[1]]@s, digits = 5, scientific = TRUE), "]"
            )
          }
          if ((abs(object@M[i, j][[1]]@m) <= 1e5 && abs(object@M[i, j][[1]]@m) >= 1e-5) &&
            (object@M[i, j][[1]]@s <= 1e5 || object@M[i, j][[1]]@s >= 1e-5)) {
            mymat[i + 1, j] <- paste(
              "[m=", format(object@M[i, j][[1]]@m, digits = 5),
              " ,s=", format(object@M[i, j][[1]]@s, digits = 5), "]"
            )
          }
          if ((abs(object@M[i, j][[1]]@m) > 1e5 || abs(object@M[i, j][[1]]@m) < 1e-5) &&
            (object@M[i, j][[1]]@s <= 1e5 && object@M[i, j][[1]]@s >= 1e-5)) {
            mymat[i + 1, j] <- paste(
              "[m=", format(object@M[i, j][[1]]@m, digits = 5, scientific = TRUE),
              " ,s=", format(object@M[i, j][[1]]@s, digits = 5), "]"
            )
          }
          if ((abs(object@M[i, j][[1]]@m) <= 1e5 && abs(object@M[i, j][[1]]@m) >= 1e-5) &&
            (object@M[i, j][[1]]@s > 1e5 || object@M[i, j][[1]]@s < 1e-5)) {
            mymat[i + 1, j] <- paste(
              "[m=", format(object@M[i, j][[1]]@m, digits = 5),
              " ,s=", format(object@M[i, j][[1]]@s, digits = 5, scientific = TRUE), "]"
            )
          }
        }
      }
    }

    rownames(mymat) <- c(
      paste(rep(" ", nchar(rownames(object@M)[1])), collapse = ""),
      row.names(object@M)
    )
    write.table(format(mymat, justify = "centre"), row.names = T, col.names = F, quote = F)
  }
)
if (!isGeneric("plot")) {
  setGeneric(
    "plot",
    function(x, y, ...) standardGeneric("plot")
  )
}
#   Plot overloading ----
#' Plot a matrix of histograms
#'
#' Plots all the histograms of a \code{MatH} object: one panel for each variable,
#' one row for each individual.
#'
#' @param x A \code{MatH} object.
#' @param y Not used.
#' @param type The type of plot: \code{"HISTO"} (default, histograms),
#' \code{"DENS"} (kernel density approximations) or \code{"BOXPLOT"} (boxplots).
#' @param border The colour of the borders (default \code{"black"}).
#' @param angL Not used (kept for compatibility).
#' @return The plot is drawn on the current device and the \code{ggplot} object is
#' returned invisibly, so that it can be further modified.
#' @examples
#' plot(BLOOD)
#' \donttest{
#' plot(BLOOD, type = "DENS", border = "blue")
#' plot(BLOOD, type = "BOXPLOT")
#' }
#' @name plot-MatH
#' @rdname plot-MatH
#' @aliases plot,MatH-method
#' @importFrom utils write.table
#' @export

setMethod(
  "plot",
  signature(x = "MatH"),
  function(x, y = "missing", type = "HISTO", border = "black", angL = 330) {
    .draw(.plot_MatH(x, type = type, border = border, angL = angL))
  }
)
#' Get or set a cell of a matrix of histograms
#'
#' \code{get.cell.MatH} returns the histogram in a cell of a \code{MatH} object;
#' \code{set.cell.MatH} assigns a histogram to a cell.
#'
#' @param object For \code{get.cell.MatH}, a \code{MatH} object; for
#' \code{set.cell.MatH}, the \code{distributionH} object to assign.
#' @param mat A \code{MatH} object.
#' @param r,c Integers, the row and the column of the cell.
#' @return \code{get.cell.MatH} returns a \code{distributionH} object,
#' \code{set.cell.MatH} the modified \code{MatH} object.
#' @examples
#' get.cell.MatH(BLOOD, r = 1, c = 1)
#' mydist <- distributionH(x = c(0, 1, 2, 3, 4), p = c(0, 0.1, 0.6, 0.9, 1))
#' MAT <- set.cell.MatH(mydist, BLOOD, r = 1, c = 1)
#' @name get.cell.MatH
#' @rdname get.cell.MatH-methods
#' @exportMethod get.cell.MatH
setGeneric("get.cell.MatH", function(object, r, c) standardGeneric("get.cell.MatH")) # OK
#' @rdname get.cell.MatH-methods
#' @aliases get.cell.MatH,MatH,numeric,numeric-method
#' @usage NULL
setMethod(
  f = "get.cell.MatH", signature = c(object = "MatH", r = "numeric", c = "numeric"),
  function(object, r, c) {
    nr <- get.MatH.nrows(object)
    nc <- get.MatH.ncols(object)
    r <- as.integer(r)
    c <- as.integer(c)
    if (r > nr | r < 1 | c < 1 | c > nc) {
      print("Indices out of range")
      return(NULL)
    } else {
      Dist <- object@M[r, c][[1]]
    }
    return(Dist)
  }
)


#' @rdname get.cell.MatH-methods
#' @aliases set.cell.MatH
#' @exportMethod set.cell.MatH
setGeneric("set.cell.MatH", function(object, mat, r, c) standardGeneric("set.cell.MatH")) # OK
#' @rdname get.cell.MatH-methods
#' @aliases set.cell.MatH,distributionH,MatH,numeric,numeric-method
#' @usage NULL
setMethod(
  f = "set.cell.MatH", signature = c(object = "distributionH", mat = "MatH", r = "numeric", c = "numeric"),
  function(object, mat, r, c) {
    nr <- get.MatH.nrows(mat)
    nc <- get.MatH.ncols(mat)
    r <- as.integer(r)
    c <- as.integer(c)
    if (r > nr | r < 1 | c < 1 | c > nc) {
      print("Indices out of range")
      return(NULL)
    } else {
      mat@M[r, c][[1]] <- object
    }
    return(mat)
  }
)
