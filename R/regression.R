## Regression analysis ----
# Multiple regression analysis of histogram variable based on Wasserstein distance----
#' Two-component regression of histogram variables
#'
#' Linear regression of a histogram variable on a set of histogram variables,
#' based on the decomposition of the L2 Wasserstein distance into a component
#' for the means and a component for the centred distributions. The parameters
#' are estimated by least squares (non-negative least squares for the
#' centred component).
#'
#' @param data A \code{MatH} object.
#' @param Yvar An integer, the column of the dependent variable.
#' @param Xvars A vector of integers, the columns of the explanatory variables.
#' @param simplify A logical value (default \code{FALSE}). If \code{TRUE}, the
#' histograms are approximated by \code{qua} equally spaced quantiles, for faster
#' computations.
#' @param qua An integer, the number of quantiles when \code{simplify = TRUE}
#' (default 20).
#' @return A named vector with the estimated parameters: the intercept and the
#' coefficients of the means (prefix \code{AV_}) and the coefficients of the
#' centred distributions (prefix \code{CEN_}).
#' @references Irpino, A., Verde, R. (2015). Linear regression for numeric
#' symbolic variables: a least squares approach based on Wasserstein Distance.
#' Advances in Data Analysis and Classification, 9, 81-106.
#' \doi{10.1007/s11634-015-0197-7}
#' @seealso \code{\link{WH.regression.two.components.predict}},
#' \code{\link{WH.regression.GOF}}
#' @examples
#' model.parameters <- WH.regression.two.components(data = BLOOD, Yvar = 1, Xvars = c(2:3))
#' model.parameters
#' @importFrom stats as.formula lm
#' @export
WH.regression.two.components <- function(data, Yvar, Xvars, simplify = FALSE, qua = 20) {
  # check input
  if (is(data) != "MatH") {
    stop("HW.regression.two.components accepts only MatH objects as input data table")
  }
  if (length(Yvar) > 1) {
    stop("Multiple choice not allowed for Y variable")
  }
  # check for missing values and set up working data
  selected <- c(1:nrow(data@M))

  Y <- data[selected, Yvar]
  X <- data[selected, Xvars]

  n <- length(selected)
  d <- ncol(X@M)
  # simplify, i.e. use a fixed number of quantiles for a rapid computation

  if (simplify) {
    pr <- c(0:qua) / qua
    for (i in 1:n) {
      tmpx <- numeric(0)
      for (q in 0:qua) {
        tmpx <- c(tmpx, compQ(Y@M[i, 1][[1]], pr[q + 1]))
      }
      Y@M[i, 1][[1]] <- new("distributionH", x = tmpx, p = pr)
      for (j in 1:d) {
        tmpx <- numeric(0)
        for (q in 0:qua) {
          tmpx <- c(tmpx, compQ(X@M[i, j][[1]], pr[q + 1]))
        }
        X@M[i, j][[1]] <- new("distributionH", x = tmpx, p = pr)
      }
    }
  }
  #
  # extract means and do multiple regression on means
  MatAver <- matrix(0, n, (d + 1)) # the matrix of means of the distributions
  for (i in 1:n) {
    MatAver[i, 1] <- Y@M[i, 1][[1]]@m
    for (j in 1:d) {
      MatAver[i, (j + 1)] <- X@M[i, j][[1]]@m
    }
  }

  colnames(MatAver) <- paste0("AV_", c(colnames(Y@M), colnames(X@M)))
  rownames(MatAver) <- rownames(Y@M)
  MatAver <- as.data.frame(MatAver)
  xnam <- colnames(MatAver)[2:(d + 1)]
  fmla <- as.formula(paste(colnames(MatAver)[1], paste(xnam, collapse = "+"), sep = "~"))
  fit <- lm(fmla, data = MatAver)
  AveCoeff <- fit$coefficient
  names(AveCoeff)[1] <- "(AV_Intercept)"

  # center data
  CENDATA <- new("MatH", nrows = n, ncols = d + 1)
  for (i in 1:n) {
    CENDATA@M[i, 1][[1]] <- Y@M[i, 1][[1]] - Y@M[i, 1][[1]]@m
    for (j in 1:d) {
      CENDATA@M[i, (j + 1)][[1]] <- X@M[i, j][[1]] - X@M[i, j][[1]]@m
    }
  }
  if (!simplify) {
    CENDATA <- registerMH(CENDATA)

    # register all data
  }
  YC <- CENDATA[, 1]
  XC <- CENDATA[, 2:(d + 1)]
  # do Non Negative Least Squares modifying
  # Lawson, Charles L.; Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
  gammas <- WH.NNLS(XC, YC)
  gammas <- as.vector(gammas)
  names(gammas) <- paste0("CEN_", colnames(X@M))
  return(parameters = c(AveCoeff, gammas))
}
## NNLS for histogram variables -----

WH.NNLS <- function(X, Y) {
  # Non Negative Least Squares for histogram variables modifying
  # Lawson, Charles L.; Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
  tol <- 1e-14
  # Check input
  if ((is(Y)[1] != is(X)[1]) && (is(X)[1] != "MatH")) {
    stop("X and Y maust be a MatH objects")
  }
  if (nrow(Y@M) != nrow(X@M)) {
    stop("X and Y must have the same number of rows")
  }
  # Lawson and Hanson NNLS algorithm
  # step 1
  P <- integer(0)
  Z <- c(1:ncol(X@M))
  gamma <- matrix(0, nrow = ncol(X@M), ncol = 1)
  stp <- 0
  # step 2
  while (stp == 0) {
    w <- WH.mat.prod(X, Y, traspose1 = TRUE) - WH.mat.prod(X, X, traspose1 = TRUE) %*% gamma

    # step 3: if Z is empty or if all w are less or equal to 0 end
    if ((length(Z) == 0) || (sum(w <= 0) == length(w))) {
      stp <- 1
    }
    else {
      P <- c(P, Z[which.max(w[Z, 1])])
      Z <- Z[-which.max(w[Z, 1])]
      stp2 <- 0
      while (stp2 == 0) {
        P <- sort(P)
        Z <- sort(Z)
        tmpX <- X[, P]
        GG <- solve(WH.mat.prod(tmpX, tmpX, traspose1 = TRUE)) %*% 
          WH.mat.prod(tmpX, Y, traspose1 = TRUE)
        zvect <- matrix(0, ncol(X@M), 1)
        zvect[P, 1] <- GG
        # step 7
        if (length(zvect[P, 1]) == sum(zvect[P, 1] > 0)) {
          gamma <- zvect
          stp <- 0
          stp2 <- 1
          # go to step 2
        }
        else { # step 8 and 9
          tmp3 <- gamma / (gamma - zvect)
          alpha <- min(tmp3[zvect[P, 1] <= 0, 1])
          # step 10
          gamma <- gamma + alpha * (zvect - gamma)

          # step 11
          # remove from P all the indices for wich gamma_j is 0 and move them to Z
          tmp4 <- P[abs(gamma[P, 1]) < tol]
          P <- P[-tmp4]
          Z <- c(Z, tmp4)
        }
      }
    }
  }
  return(gamma)
}
## diagnostics measures to be implemented
#' Predictions of a two-component regression of histogram variables
#'
#' Predicts the histograms of the dependent variable using the parameters
#' estimated by \code{\link{WH.regression.two.components}}.
#'
#' @param data A \code{MatH} object with the explanatory variables (in the same
#' order used for the estimation).
#' @param parameters The named vector of parameters returned by
#' \code{\link{WH.regression.two.components}}.
#' @return A \code{MatH} object with one column, the predicted histograms.
#' @seealso \code{\link{WH.regression.GOF}}
#' @examples
#' model.parameters <- WH.regression.two.components(data = BLOOD, Yvar = 1, Xvars = c(2:3))
#' Predicted.BLOOD <- WH.regression.two.components.predict(
#'   data = BLOOD[, 2:3],
#'   parameters = model.parameters
#' )
#' @export
WH.regression.two.components.predict <- function(data, parameters) {
  if (is(data) != "MatH") {
    stop("HW.regression.two.components.predict accepts only MatH objects as input data table")
  }
  if (ncol(data@M) != (length(parameters) - 1) / 2) {
    stop("Predictors (colums of data) are not consistent with the number of parameters (that must be (num.of colums of data)x2+1)")
  }
  npred <- ncol(data@M)
  if (sum(abs(parameters[(1 + npred + 1):length(parameters)]) - (parameters[(1 + npred + 1):length(parameters)])) > 0) {
    stop(cat("HW.regression.two.components.predict accepts only last ", npred, " parameters positives"))
  }
  indiv <- nrow(data@M)
  predictions <- new("MatH", nrows = indiv, ncols = 1, names.cols = c("Predicted"), names.rows = rownames(data@M))
  parameters <- as.numeric(parameters)
  for (i in 1:indiv) {
    tmp.pred <- new("distributionH", x = c(parameters[1], parameters[1]), p = c(0, 1))
    for (j in 1:npred) {
      tmpdist <- data@M[i, j][[1]]
      tmp.pred <- tmp.pred + (parameters[1 + j] * tmpdist@m) + (parameters[1 + npred + j] * (tmpdist - tmpdist@m))
    }
    predictions@M[i, 1][[1]] <- tmp.pred
  }
  return(predictions)
}
#' Goodness of fit of a two-component regression of histogram variables
#'
#' Computes goodness of fit indices of the predictions of a
#' \code{\link{WH.regression.two.components}} model.
#'
#' @param observed A \code{MatH} object with one column, the observed histograms.
#' @param predicted A \code{MatH} object with one column, the predicted
#' histograms.
#' @return A list with the root mean squared error based on the L2 Wasserstein
#' distance (\code{RMSE_W}), the \code{OMEGA} index and the pseudo-R2 index
#' (\code{PSEUDOR2}, with the details of the decomposition of the sum of
#' squares).
#' @references Irpino, A., Verde, R. (2015). Linear regression for numeric
#' symbolic variables: a least squares approach based on Wasserstein Distance.
#' Advances in Data Analysis and Classification, 9, 81-106.
#' \doi{10.1007/s11634-015-0197-7}
#' @examples
#' model.parameters <- WH.regression.two.components(data = BLOOD, Yvar = 1, Xvars = c(2:3))
#' Predicted.BLOOD <- WH.regression.two.components.predict(
#'   data = BLOOD[, 2:3],
#'   parameters = model.parameters
#' )
#' WH.regression.GOF(observed = BLOOD[, 1], predicted = Predicted.BLOOD)
#' @export
WH.regression.GOF <- function(observed, predicted) {
  indiv <- nrow(observed@M)
  WassD2 <- matrix(0, nrow = indiv, ncol = 1)
  TOTSSQ <- WH.SSQ(observed)
  MO <- WH.vec.mean(observed)
  SSQR <- 0
  RMSE_W <- 0
  NUM_OMEGA <- 0
  DEN_OMEGA <- 0
  for (i in 1:indiv) {
    WassD2[i, 1] <- WassSqDistH(observed@M[i, 1][[1]], predicted@M[i, 1][[1]])
    SSQR <- SSQR + WassSqDistH(predicted@M[i, 1][[1]], MO)
    NUM_OMEGA <- NUM_OMEGA + (predicted@M[i, 1][[1]]@m - MO@m)^2 + (predicted@M[i, 1][[1]]@s)^2
    DEN_OMEGA <- DEN_OMEGA + (observed@M[i, 1][[1]]@m - MO@m)^2 + (observed@M[i, 1][[1]]@s)^2
  }

  return(indices = list(
    RMSE_W = sqrt(sum(WassD2) / indiv),
    OMEGA = NUM_OMEGA / DEN_OMEGA,
    PSEUDOR2 = list(
      index = max(0, min(SSQR / TOTSSQ, 1 - sum(WassD2) / TOTSSQ)),
      details = c(
        TotSSQ = TOTSSQ, SSQ.R = SSQR, SSQ.E = sum(WassD2),
        Bias = TOTSSQ - SSQR - sum(WassD2), SSQ.R.rel = SSQR / TOTSSQ,
        SSQ.E.rel = sum(WassD2) / TOTSSQ,
        SSQ.bias.rel = (TOTSSQ - SSQR - sum(WassD2)) / TOTSSQ
      )
    )
  ))
}

WH.regr.two.comp.Bootstrap <- function(data, Yvar, Xvars, simplify = FALSE,
                                       qua = 20, rep = 100, GOF = FALSE) {
  nobj <- get.MatH.nrows(data)
  ind <- sample(nobj, nobj, replace = TRUE)
  pars <- WH.regression.two.components(data[ind, ], Yvar, Xvars, simplify = simplify, qua = qua)

  if (GOF == TRUE) {
    OBS <- data[ind, Yvar]
    PRED <- WH.regression.two.components.predict(data[ind, Xvars], pars)
    idxs <- WH.regression.GOF(OBS, PRED)
    MP <- c(pars, RMSEW = idxs$RMSE_W, OMEGA = idxs$OMEGA, PSEUDOR2 = idxs$PSEUDOR2$index)
  } else {
    MP <- pars
  }
  for (i in 1:rep) {
    ind <- sample(nobj, nobj, replace = TRUE)
    pars <- WH.regression.two.components(data[ind, ], Yvar, Xvars, simplify = simplify, qua = qua)

    if (GOF == TRUE) {
      OBS <- data[ind, Yvar]
      PRED <- WH.regression.two.components.predict(data[ind, Xvars], pars)
      idxs <- WH.regression.GOF(OBS, PRED)
      TMP <- c(pars, RMSEW = idxs$RMSE_W, OMEGA = idxs$OMEGA, PSEUDOR2 = idxs$PSEUDOR2$index)
      MP <- rbind(MP, TMP)
    } else {
      MP <- rbind(MP, pars)
    }
  }
  print(summary(MP))
  return(MP)
}
