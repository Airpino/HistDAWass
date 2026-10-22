## 1v PCA for quantiles -------
# Principal components analysis  of histogram variable based on Wasserstein distance----
#' Principal component analysis of a histogram variable
#'
#' Principal component analysis (PCA) of one histogram variable based on the L2
#' Wasserstein distance: a centred (not standardised) PCA of the quantiles of the
#' histograms. It provides a dimensional reduction and a visualization of the
#' distributions.
#'
#' @param data A \code{MatH} object.
#' @param var An integer, the column of the variable to analyse.
#' @param quantiles An integer, the number of quantiles used in the analysis
#' (default 10).
#' @param plots A logical value (default \code{TRUE}). If \code{TRUE}, for each
#' factorial plane four plots (the variables and the individuals, plain and
#' coloured by mean and by standard deviation) are drawn on the current graphic
#' device, in a 2 x 2 layout.
#' @param listaxes A vector of integers, the axes of the factorial planes to
#' plot (default \code{1:4}, i.e. the planes 1-2 and 3-4).
#' @param axisequal A logical value (default \code{FALSE}). If \code{TRUE}, the x
#' and y axes of the plots have the same scale.
#' @param qcut A number between 0.5 and 1 (default 1): densities higher than
#' their \code{qcut} quantile are cut in the plots, to avoid very peaked
#' densities.
#' @param outl A number between 0 (default) and 0.5, the probability mass
#' removed from each tail of the distributions (e.g. 0.1 removes 10\% from each
#' tail).
#' @return A list with \code{PCAout}, the results of
#' \code{\link[FactoMineR]{PCA}}, \code{WASSVARIANCE} (not computed, always 0)
#' and \code{INERTIA}, the total inertia.
#' @references Verde, R., Irpino, A., Balzanella, A. (2016). Dimension reduction
#' techniques for distributional symbolic data. IEEE Transactions on
#' Cybernetics, 46(2), 344-355. \doi{10.1109/TCYB.2015.2389653}
#' @seealso \code{\link{WH.MultiplePCA}}, \code{\link{WH.plot_multiple_Spanish.funs}}
#' @examples
#' results <- WH.1d.PCA(data = BLOOD, var = 1, listaxes = c(1:2))
#' @importFrom FactoMineR PCA
#' @importFrom graphics plot abline axis hist plot.new plot.window polygon segments text title par
#' @importFrom grDevices rgb
#' @importFrom stats density quantile
#' @export

WH.1d.PCA <- function(data, var, quantiles = 10, plots = TRUE, listaxes = c(1:4), axisequal = FALSE, qcut = 1, outl = 0) {
  if (is(data)[1] != "MatH") {
    stop("Input data must be a MatH object (A matrix of histograms)")
  }
  # Build the matrix of Quantiles
  VARS <- ncol(data@M)
  INDIV <- nrow(data@M)
  if ((qcut < 0.5) || (qcut > 1)) {
    qcut <- 1
  }
  if ((outl < 0) || (outl > 0.5)) {
    outl <- 0
  }
  if ((outl == 0.5)) {
    outl <- 0.49
  }
  if (missing(var)) {
    var <- 1
    varname <- colnames(data@M)[1]
    cat(paste("Var is missing, We do a PCA on variable ", varname, "\n"))
  } else {
    if (length(var) > 1) {
      varname <- colnames(data@M)[var[1]]
      cat(paste("Var has several values, We do a PCA on variable -->", var[1], "\n"))
      var <- var[1]
    }
    else {
      if (var > VARS) {
        stop(paste("The variables are less than ", var))
      }
      else {
        varname <- colnames(data@M)[var]
        cat(paste("We do a PCA on variable ---> ", varname, "\n"))
      }
    }
  }
  data <- data[, var]
  # check indiv and remove empty rows
  toremove <- numeric(0)
  for (i in 1:INDIV) {
    if (length(data@M[i, 1][[1]]@x) < 2) {
      toremove <- c(toremove, i)
    }
  }
  if (length(toremove) > 0) {
    data@M <- as.matrix(data@M[-toremove, 1])
  }
  colnames(data@M) <- varname
  INDIV <- nrow(data@M)
  if (INDIV < 2) {
    stop("At least two individuals are necessary, please check data")
  }

  # Create the quantile matrix
  MATQ <- matrix(0, nrow = INDIV, ncol = (quantiles + 1))
  if (outl == 0) {
    p <- c(0, 1:quantiles) / quantiles
    namec <- list("Min")
    for (j in 2:(quantiles)) {
      namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
    }
    namec <- c(namec, "Max")
  }
  else {
    p <- c(0, 1:quantiles) / quantiles * (1 - 2 * outl) + outl
    namec <- list(paste0("Q(", format(p[1], digits = 2, nsmall = 2), ")"))
    for (j in 2:length(p)) {
      namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
    }
  }

  for (i in 1:INDIV) {
    MATQ[i, ] <- .qfun(data@M[i, 1][[1]], p)
  }
  rownames(MATQ) <- rownames(data@M)
  colnames(MATQ) <- namec

  # do a PCA
  vmeans <- numeric(0)
  vstds <- numeric(0)

  for (ind in 1:INDIV) {
    vmeans <- c(vmeans, data@M[ind, 1][[1]]@m)
    vstds <- c(vstds, data@M[ind, 1][[1]]@s)
  }
  minM <- min(vmeans)
  maxM <- max(vmeans)
  minS <- min(vstds)
  maxS <- max(vstds)

  MATF <- cbind(MATQ / sqrt((quantiles + 1)), vmeans, vstds)
  res.pca <- PCA(MATF, quanti.sup = c((ncol(MATF) - 1), ncol(MATF)), scale.unit = FALSE, graph = F)

  VARW <- 0 # WH.var.covar(data)
  TOTINE <- sum(res.pca$eig[, 1])
  ## plotting PCA results ----------------
  if (plots) {
    # since 1.0.9: the plots are drawn on the current device, four panels
    # for each plane (up to 1.0.8 a new device was opened for each plot)
    oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(oldpar), add = TRUE)
    planes <- matrix(0, 1, 2)
    for (cc in 1:floor(length(listaxes) / 2)) {
      planes <- rbind(planes, c(listaxes[(cc * 2 - 1)], listaxes[(cc * 2)]))
    }
    planes <- as.matrix(planes[-1, ])
    dim(planes) <- c(floor(length(listaxes) / 2), 2)
    if ((length(listaxes) %% 2) == 1) {
      planes <- rbind(planes, c(
        listaxes[(length(listaxes) - 1)],
        listaxes[length(listaxes)]
      ))
    }

    # plot Spanish-fan plot
    for (pl in 1:nrow(planes)) {
      axe1 <- planes[pl, 1]
      axe2 <- planes[pl, 2]

      labX <- paste("Axis ", axe1, " (", format(res.pca$eig[axe1, 2], digits = 2, nsmall = 2), "%)")
      labY <- paste("Axis ", axe2, " (", format(res.pca$eig[axe2, 2], digits = 2, nsmall = 2), "%)")
      CVAR <- res.pca$var$coord[, c(axe1, axe2)]
      plot.new()

      if (axisequal) {
        xrange <- c(min(0, min(cbind(CVAR[, 1], CVAR[, 2]))), max(0, max(cbind(CVAR[, 1], CVAR[, 2]))))
        yrange <- xrange
        plot.window(xrange, yrange)
        title(xlab = labX)
        title(ylab = labY)

        #         plot(type="n",
        #              xlab=labX,ylab=labY)
        segments(xrange[1], 0, xrange[2], 0)
        segments(yrange[1], 0, yrange[2], 0)
      }
      else {
        plot.window(
          c(min(0, CVAR[, 1]), max(0, CVAR[, 1])),
          c(min(0, CVAR[, 2]), max(0, CVAR[, 2]))
        )
        title(xlab = labX)
        title(ylab = labY)
        #       plot(c(min(0,CVAR[,1]),max(0,CVAR[,1])), c(min(0,CVAR[,2]),max(0,CVAR[,2])),type="n",
        #            xlab=labX,ylab=labY)Y)
        segments(min(0, CVAR[, 1]), 0, max(0, CVAR[, 1]), 0)
        segments(0, min(0, CVAR[, 2]), 0, max(0, CVAR[, 2]))
      }
      title("PCA Variable plot (Spanish-fan plot)")
      if ((nrow(CVAR) %% 2) == 0) {
        centr <- nrow(CVAR) / 2
      } else {
        centr <- (nrow(CVAR) - 1) / 2
      }
      centr <- nrow(CVAR) / 2
      cxl <- 0.8
      for (tr in 2:nrow(CVAR)) {
        centrality <- (abs(tr - centr - 1) / (centr - 1))

        # cat(centrality,"\n")
        x <- c(0, CVAR[(tr - 1), 1], CVAR[tr, 1], 0)
        y <- c(0, CVAR[(tr - 1), 2], CVAR[tr, 2], 0)
        red <- 1
        green <- centrality
        # cat(red,green,"\n")
        polygon(x, y, col = rgb(red, green, 0, 0.7), lty = 1, lwd = 1, border = "black")
      }
      text(x = CVAR[1, 1], y = CVAR[1, 2], labels = colnames(MATQ)[1], cex = cxl)
      if (nrow(CVAR) < 8) {
        for (tr in 2:(nrow(CVAR) - 1)) {
          text(x = CVAR[tr, 1], y = CVAR[tr, 2], labels = colnames(MATQ)[tr], cex = cxl)
        }
      }
      else {
        text(
          x = CVAR[ceiling(nrow(CVAR) / 8), 1], y = CVAR[ceiling(nrow(CVAR) / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) / 8)], cex = cxl
        )
        text(
          x = CVAR[ceiling(nrow(CVAR) * 2 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 2 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 2 / 8)], cex = cxl
        )
        text(
          x = CVAR[ceiling(nrow(CVAR) * 3 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 3 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 3 / 8)], cex = cxl
        )
        text(x = CVAR[centr + 1, 1], y = CVAR[centr + 1, 2], labels = colnames(MATQ)[centr + 1], cex = cxl)
        text(
          x = CVAR[ceiling(nrow(CVAR) * 5 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 5 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 5 / 8)], cex = cxl
        )
        text(
          x = CVAR[ceiling(nrow(CVAR) * 6 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 6 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 6 / 8)], cex = cxl
        )
        text(
          x = CVAR[ceiling(nrow(CVAR) * 7 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 7 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 7 / 8)], cex = cxl
        )
      }
      text(
        x = CVAR[nrow(CVAR), 1], y = CVAR[nrow(CVAR), 2],
        labels = colnames(MATQ)[ncol(MATQ)], cex = cxl
      )
      axis(1, pos = 0, cex.axis = 0.7) # , at=xM,labels=labX)
      axis(2, pos = 0, cex.axis = 0.7) # , at=yM,labels=labY)

      # plot individuals

      # layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

      xm <- min(res.pca$ind$coord[, axe1])
      xM <- max(res.pca$ind$coord[, axe1])
      ym <- min(res.pca$ind$coord[, axe2])
      yM <- max(res.pca$ind$coord[, axe2])
      Xlen <- xM - xm
      Ylen <- yM - ym

      # Compute scaling factor for domain
      matX <- matrix(0, 128, INDIV)
      matD <- matrix(0, 128, INDIV)
      for (ind in 1:INDIV) {
        distr <- data@M[ind, 1][[1]]
        # generate 200 random points according to the QF
        rn <- 200
        xn <- .qfun(distr, c(0:(rn - 1)) / rn)
        d <- density(xn, n = 128)
        matX[, ind] <- d$x
        matD[, ind] <- d$y
      }
      MinX <- min(matX)
      MaxX <- max(matX)
      RX <- MaxX - MinX
      q95D <- quantile(matD, probs = qcut)
      matD[matD > q95D] <- as.numeric(q95D)
      MaxD <- max(matD)

      meanXRange <- diff(range(apply(matX, 2, FUN = range)))
      xfact <- 0.3
      yfact <- 0.1
      if (axisequal) {
        xm1 <- min(xm, ym)
        xM1 <- max(xM, yM)
        ym1 <- xm1
        yM1 <- xM1
        Xlen <- xM1 - xm1
        Ylen <- yM1 - ym1
        xm <- (xM1 + xm1) / 2 - Xlen / 2 * (1 + 1.5 * xfact)
        xM <- (xM1 + xm1) / 2 + Xlen / 2 * (1 + 1.5 * xfact)
        ym <- (yM1 + ym1) / 2 - Ylen / 2 * (1 + 0.5 * yfact)
        yM <- (yM1 + ym1) / 2 + Ylen / 2 * (1 + 1.5 * yfact)
      } else {
        xm1 <- xm
        xM1 <- xM
        ym1 <- ym
        yM1 <- yM
        xm <- (xM1 + xm1) / 2 - Xlen / 2 * (1 + 1.5 * xfact)
        xM <- (xM1 + xm1) / 2 + Xlen / 2 * (1 + 1.5 * xfact)
        ym <- (yM1 + ym1) / 2 - Ylen / 2 * (1 + 0.5 * yfact)
        yM <- (yM1 + ym1) / 2 + Ylen / 2 * (1 + 1.5 * yfact)
      }
      plot.new()
      plot.window(c(xm, xM), c(ym, yM))
      title(main = "PCA plot of distributions", sub = "Smoothed distributions")
      for (ind in 1:INDIV) {
        x <- matX[, ind]
        y <- matD[, ind]
        x <- (x - data@M[ind, 1][[1]]@m) / meanXRange * xfact * Xlen + res.pca$ind$coord[ind, axe1]
        x <- c(x, x[length(x)], x[1])
        y <- c(y, 0, 0) / q95D * yfact * Ylen + res.pca$ind$coord[ind, axe2]
        polygon(x, y, col = rgb(1, 1, 0, 0.6))
      }
      for (ind in 1:INDIV) {
        text(res.pca$ind$coord[ind, axe1], res.pca$ind$coord[ind, axe2],
          label = rownames(MATQ)[ind], pos = 1, cex = 0.7
        )
      }
      axis(1, pos = 0, cex.axis = 0.7) # , at=xM,labels=labX)
      axis(2, pos = 0, cex.axis = 0.7) # , at=yM,labels=labY)
      segments(xm, 0, xM, 0)
      segments(0, ym, 0, yM)

      title(xlab = labX)
      title(ylab = labY)
      plot.new()
      plot.window(c(xm, xM), c(ym, yM))
      title(main = "Coloured using the mean values")
      trasp <- 0.6

      for (ind in 1:INDIV) {
        x <- matX[, ind]
        y <- matD[, ind]
        x <- (x - data@M[ind, 1][[1]]@m) / meanXRange * xfact * Xlen + res.pca$ind$coord[ind, axe1]
        x <- c(x, x[length(x)], x[1])
        y <- c(y, 0, 0) / q95D * yfact * Ylen + res.pca$ind$coord[ind, axe2]
        red <- (data@M[ind, 1][[1]]@m - minM) / (maxM - minM)
        if (red > 0.5) {
          green <- 1 - red
          blue <- 0
        }
        else {
          green <- 1 - red
          red <- 0
          blue <- 1 - green
        }
        polygon(x, y, col = rgb(red, green, blue, trasp))
      }
      for (ind in 1:INDIV) {
        text(res.pca$ind$coord[ind, axe1], res.pca$ind$coord[ind, axe2],
          label = rownames(MATQ)[ind], pos = 1, cex = 0.7
        )
      }
      axis(1, pos = 0, cex.axis = 0.7) # , at=xM,labels=labX)
      axis(2, pos = 0, cex.axis = 0.7) # , at=yM,labels=labY)
      segments(xm, 0, xM, 0)
      segments(0, ym, 0, yM)

      abline(0, res.pca$quanti.sup$coord[1, 2] / res.pca$quanti.sup$coord[1, 1], lty = 3, lwd = 2, col = "red")
      title(xlab = labX)
      title(ylab = labY)
      plot.new()
      plot.window(c(xm, xM), c(ym, yM))
      title(main = "Coloured using std values")
      for (ind in 1:INDIV) {
        x <- matX[, ind]
        y <- matD[, ind]
        x <- (x - data@M[ind, 1][[1]]@m) / meanXRange * xfact * Xlen + res.pca$ind$coord[ind, axe1]
        x <- c(x, x[length(x)], x[1])
        y <- c(y, 0, 0) / q95D * yfact * Ylen + res.pca$ind$coord[ind, axe2]
        red <- (data@M[ind, 1][[1]]@s - minS) / (maxS - minS)
        if (red > 0.5) {
          green <- 1 - red
          blue <- 0
        }
        else {
          green <- 1 - red
          red <- 0
          blue <- 1 - green
        }
        polygon(x, y, col = rgb(red, green, blue, trasp))
      }
      for (ind in 1:INDIV) {
        text(res.pca$ind$coord[ind, axe1], res.pca$ind$coord[ind, axe2],
          label = rownames(MATQ)[ind], pos = 1, cex = 0.7
        )
      }
      axis(1, pos = 0, cex.axis = 0.7) # , at=xM,labels=labX)
      axis(2, pos = 0, cex.axis = 0.7) # , at=yM,labels=labY)
      segments(xm, 0, xM, 0)
      segments(0, ym, 0, yM)

      abline(0, res.pca$quanti.sup$coord[2, 2] / res.pca$quanti.sup$coord[2, 1], lty = 3, lwd = 2, col = "red")

      title(xlab = labX)
      title(ylab = labY)
    }
  }
  return(list(PCAout = res.pca, WASSVARIANCE = VARW, INERTIA = TOTINE))
}
# Principal components analysis  of a set of histogram variables based on Wasserstein distance----
#' Principal component analysis of a set of histogram variables
#'
#' (Beta version) Principal component analysis of several histogram variables
#' based on the L2 Wasserstein distance: a multiple factor analysis
#' (\code{\link[FactoMineR]{MFA}}) of the quantiles of the histograms, where each
#' variable is a group of quantiles. It extends \code{\link{WH.1d.PCA}}.
#'
#' @inheritParams WH.1d.PCA
#' @param list.of.vars A vector of integers, the columns of the variables to
#' analyse.
#' @param graph A logical value (default \code{TRUE}): if \code{TRUE} the plots of
#' \code{\link[FactoMineR]{MFA}} are drawn.
#' @return The results of \code{\link[FactoMineR]{MFA}}.
#' @seealso \code{\link{WH.plot_multiple_Spanish.funs}},
#' \code{\link{WH.plot_multiple_indivs}}
#' @examples
#' res <- WH.MultiplePCA(BLOOD, list.of.vars = c(1:3), graph = FALSE)
#' @importFrom FactoMineR MFA
#' @importFrom stats rnorm
#' @export
## MFA PCA for quantiles -------
WH.MultiplePCA <- function(data, list.of.vars, quantiles = 10, outl = 0, graph = TRUE) {
  data <- data[, list.of.vars]
  VARS <- ncol(data@M)
  INDIV <- nrow(data@M)
  MATQ <- matrix(0, nrow = INDIV, ncol = VARS * (quantiles + 1))
  if ((outl < 0) || (outl > 0.5)) {
    outl <- 0
  }
  if ((outl == 0.5)) {
    outl <- 0.49
  }
  # browser()
  # create the names of the variables
  #################################################
  COLnames <- list()
  if (outl == 0) {
    p <- c(0, 1:quantiles) / quantiles
    for (v in 1:VARS) {
      namec <- list("Min")
      for (j in 2:(quantiles)) {
        namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
      }
      namec <- c(namec, "Max")
      COLnames <- c(COLnames, paste(substr(colnames(data@M)[v], 1, 4), namec, sep = "."))
      for (i in 1:INDIV) {
        MATQ[i, ((v - 1) * (quantiles + 1) + 1:(quantiles + 1))] <- .qfun(data@M[i, v][[1]], p)
      }
    }
  }
  else {
    p <- c(0, 1:quantiles) / quantiles * (1 - 2 * outl) + outl
    namec <- list(paste0("Q(", format(p[1], digits = 2, nsmall = 2), ")"))
    for (j in 2:length(p)) {
      namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
    }
    for (v in 1:VARS) {
      namec <- list(paste0("Q(", format(p[1], digits = 2, nsmall = 2), ")"))
      for (j in 2:length(p)) {
        namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
      }
      COLnames <- c(COLnames, paste(substr(colnames(data@M)[v], 1, 4), namec, sep = "."))
      for (i in 1:INDIV) {
        MATQ[i, ((v - 1) * (quantiles + 1) + 1:(quantiles + 1))] <- .qfun(data@M[i, v][[1]], p)
      }
    }
  }
  #################################################
  # COLnames=list();
  # for (v in 1:VARS){
  #   namec=list("Min")
  #   for (j in 2:(quantiles)){
  #     namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
  #   }
  #   namec=c(namec, "Max")
  #   COLnames=c(COLnames, paste(substr(colnames(data@M)[v],1,4),namec,sep="."))
  #   for (i in 1:INDIV)
  #   {
  #     for (j in 1:(quantiles+1)){
  #       MATQ[i,((v-1)*(quantiles+1)+j)]=compQ(data@M[i,v][[1]],p[j])
  #     }
  #   }
  # }

  rownames(MATQ) <- rownames(data@M)
  colnames(MATQ) <- COLnames

  # do a MFA
  MATF <- cbind(MATQ / sqrt((quantiles + 1)))
  # check all columns
  for (cc in 1:ncol(MATF)) {
    if (sd(MATF[, cc]) < 1e-100) {
      MATF[, cc] <- MATF[, cc] + rnorm(nrow(MATF), 0, 1e-20)
    }
  }
  MFA.res <- MFA(MATF,
    group = c(rep(quantiles + 1, length(list.of.vars))), type = c(rep("c", length(list.of.vars))),
    ncp = 3, name.group = colnames(data@M), graph = graph
  )
  # do Plots!!!!
  # plot Spanish-fan plot
  # plot individuals
  return(MFA.res)
}

# Plotting Spanish fun plots for Multiple factor analysis of Histogram Variables ----
#' Plot the variables of a PCA of histogram variables (Spanish-fan plot)
#'
#' Plots the quantiles of the histogram variables on a factorial plane after
#' \code{\link{WH.MultiplePCA}} (or \code{\link{WH.1d.PCA}}): the quantiles of
#' each variable form a fan.
#'
#' @param res The results of \code{\link{WH.MultiplePCA}} or of
#' \code{\link{WH.1d.PCA}}.
#' @param axes A vector of two integers, the factorial axes (default
#' \code{c(1, 2)}).
#' @param var A vector of integers, the variables to plot.
#' @param LABS A logical value (default \code{TRUE}): if \code{TRUE} the quantiles
#' are labelled.
#' @param multi A logical value: \code{TRUE} (default) if \code{res} comes from
#' \code{WH.MultiplePCA}, \code{FALSE} if it comes from \code{WH.1d.PCA}.
#' @param corplot A logical value: if \code{TRUE} (default) the correlations are
#' plotted (circle of correlations), otherwise the coordinates.
#' @return The plot is drawn on the current device and the \code{ggplot} object
#' is returned invisibly.
#' @examples
#' res <- WH.MultiplePCA(BLOOD, list.of.vars = c(1:3), graph = FALSE)
#' WH.plot_multiple_Spanish.funs(res, axes = c(1, 2), var = c(1:3))
#' @export
WH.plot_multiple_Spanish.funs <- function(res,
                                          axes = c(1, 2), var = 1, LABS = TRUE,
                                          multi = TRUE, corplot = TRUE) {
  # require(ggplot2)
  if (multi) {
    labX <- paste("Axis ", axes[1], " (", format(res$eig[axes[1], 2],
      digits = 2, nsmall = 2
    ), "%)")
    labY <- paste("Axis ", axes[2], " (", format(res$eig[axes[2], 2],
      digits = 2, nsmall = 2
    ), "%)")
  }
  else {
    labX <- paste("Axis ", axes[1], " (", format(res$PCAout$eig[axes[1], 2],
      digits = 2, nsmall = 2
    ), "%)")
    labY <- paste("Axis ", axes[2], " (", format(res$PCAout$eig[axes[2], 2],
      digits = 2, nsmall = 2
    ), "%)")
  }
  v1 <- numeric()
  v2 <- numeric()
  tt <- character()
  tt2 <- numeric()
  VA <- numeric()
  x1 <- numeric()
  x2 <- numeric()
  hj <- numeric()
  vj <- numeric()
  if (!multi) {
    var <- 1
  }
  for (j in var) {

    # plot Spanish-fan plot
    if (multi) {
      els <- which(res$summary.quanti$group == j)
    }
    else {
      els <- c(1:nrow(res$PCAout$var$coord))
    }
    tmp_v1 <- numeric()
    tmp_v2 <- numeric()
    tmp_tt <- character()
    tmp_tt2 <- numeric()
    tmp_VA <- numeric()
    if (multi) {
      if (corplot) {
        CVAR <- res$global.pca$var$cor[els, c(axes[1], axes[2])]
      }
      else {
        CVAR <- res$global.pca$var$coord[els, c(axes[1], axes[2])]
      }
    }
    else {
      if (corplot) {
        CVAR <- res$PCAout$var$cor[els, c(axes[1], axes[2])]
      }
      else {
        CVAR <- res$PCAout$var$coord[els, c(axes[1], axes[2])]
      }
    }

    for (pp in 1:(nrow(CVAR) - 1)) {
      tmp_v1 <- c(tmp_v1, CVAR[pp, 1], 0, CVAR[pp + 1, 1], CVAR[pp, 1])
      tmp_v2 <- c(tmp_v2, CVAR[pp, 2], 0, CVAR[pp + 1, 2], CVAR[pp, 2])
      tmp_tt <- c(tmp_tt, paste(as.character(j), rep(as.character(pp), 4), sep = "_"))
      tmp_tt2 <- c(tmp_tt2, rep(pp, 4))
      tmp_VA <- c(tmp_VA, rep(j, 4))
    }
    tmp_tt2 <- abs(tmp_tt2 - mean(tmp_tt2))
    tmp_tt2 <- 1 - tmp_tt2 / max(tmp_tt2)
    v1 <- c(v1, tmp_v1, 0, 0, 0, 0)
    v2 <- c(v2, tmp_v2, 0, 0, 0, 0)
    tt <- c(tt, tmp_tt, paste(as.character(j), rep(as.character(pp), 4), sep = "_"))
    tt2 <- c(tt2, tmp_tt2, 1, 1, 1, 1)
    tmp_VA <- c(tmp_VA, rep(j, 4))
    VA <- c(VA, tmp_VA)
    x1 <- c(x1, CVAR[, 1])
    x2 <- c(x2, CVAR[, 2])
  }

  d <- data.frame(v1 = v1, v2 = v2, tt = tt, tt2 = tt2, VA = VA)
  # CVAR=res$global.pca$var$cor[,c(axes[1],axes[2])]
  txt_lab <- data.frame(v1 = x1, v2 = x2, hj = (-sign(x1) + 1) / 2, vj = (-sign(x2) + 1) / 2)
  ### the circle of correlation
  angle <- seq(-pi, pi, length = 250)
  x_c <- sin(angle)
  y_c <- cos(angle)
  df <- data.frame(x_c = x_c, y_c = y_c)
  ##
  l1 <- data.frame(x = c(-1.2, 1.2), y = c(0, 0))
  l2 <- data.frame(x = c(0, 0), y = c(-1.2, 1.2))


  txt_lab$lab <- rownames(txt_lab)
  p <- ggplot(d, aes(x = .data$v1, y = .data$v2)) +
    geom_polygon(aes(group = .data$tt, fill = as.factor(.data$VA), alpha = .data$tt2), colour = "grey30")
  if (LABS == TRUE) {
    p <- p + geom_text(data = txt_lab, aes(
      x = .data$v1, y = .data$v2,
      label = .data$lab,
      hjust = .data$hj, vjust = .data$vj
    ))
  }
  if (corplot) {
    p <- p + xlim(-1.2, 1.2) + ylim(-1.2, 1.2) +
      ggtitle("Correlation plot of variables (Spanish-fan plot)")
  }
  else {
    XL <- range(d$v1)
    XL <- (XL - mean(XL)) * 1.2 + mean(XL)
    YL <- range(d$v2)
    YL <- (YL - mean(YL)) * 1.2 + mean(YL)
    p <- p + xlim(XL) + ylim(YL) +
      ggtitle("Plot of variables (Spanish-fan plot)")
  }

  p <- p + theme_bw() +
    xlab(labX) + ylab(labY) +
    theme(legend.position = "none")
  if (corplot) {
    p <- p + coord_fixed(ratio = 1)
  }
  p <- p + theme(plot.title = element_text(hjust = 0.5))
  if (corplot) {
    p <- p + geom_path(data = df, aes(x = .data$x_c, y = .data$y_c), inherit.aes = FALSE, linetype = 2)
  }
  p <- p + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)


  print(p)

  return(invisible(p))
}

# Plot histograms of individuals after a  Multiple factor analysis of Histogram Variables ----
#' Plot the histograms of the individuals on a factorial plane
#'
#' (Beta version) Plots the histograms (or densities) of one variable for each
#' individual, placed at the coordinates of the individuals on a factorial plane
#' of \code{\link{WH.MultiplePCA}}.
#'
#' @param data A \code{MatH} object, the data of the analysis.
#' @param res The results of \code{\link{WH.MultiplePCA}}.
#' @param axes A vector of two integers, the factorial axes (default
#' \code{c(1, 2)}).
#' @param indiv A vector of integers, the individuals (rows) to plot. Default 0:
#' all the individuals.
#' @param var An integer, the variable to plot.
#' @param strx A number, the average width of the drawn distributions as a
#' fraction of the range of the coordinates of the individuals on the x axis
#' (default 0.1: on average, the domain of a distribution takes one tenth of
#' the x axis). All the distributions share the same scale, so their spreads
#' remain comparable.
#' @param stry A number, the height of the highest density as a fraction of the
#' range of the coordinates of the individuals on the y axis (default 1/8).
#' All the densities share the same scale; spikes (bins of almost zero
#' width) are avoided by drawing histograms with very narrow bins on 30
#' equal-width bins.
#' @param HISTO A logical value: if \code{TRUE} (default) histograms are drawn,
#' otherwise smoothed densities.
#' @param coor A matrix of coordinates of the individuals. Default 0: the
#' coordinates in \code{res}.
#' @param stat The statistic used to label and colour the histograms in the
#' second plot: \code{"mean"} (default), \code{"std"}, \code{"skewness"} or
#' \code{"kurtosis"}.
#' @details The distributions are measured on the scale of the original
#' variable, while the factorial plane has the scale of the principal
#' components: \code{strx} and \code{stry} set the size of the drawn
#' distributions automatically, as fractions of the ranges of the coordinates
#' of the individuals. Each distribution is drawn with its mean at the
#' coordinates of the individual.
#' @return The two plots are drawn on the current device and a list with the two
#' \code{ggplot} objects is returned invisibly.
#' @examples
#' \donttest{
#' results <- WH.MultiplePCA(BLOOD, list.of.vars = c(1:3), graph = FALSE)
#' WH.plot_multiple_indivs(BLOOD, results, axes = c(1, 2), var = 1, HISTO = TRUE)
#' WH.plot_multiple_indivs(BLOOD, results, axes = c(1, 2), var = 1, HISTO = FALSE)
#' }
#' @export
WH.plot_multiple_indivs <- function(data, res, axes = c(1, 2), indiv = 0, var = 1,
                                    strx = 0.1, stry = 1 / 8, HISTO = TRUE, coor = 0,
                                    stat = "mean") {
  # require(ggplot2)

  if (indiv[1] == 0) {
    el_ind <- c(1:get.MatH.nrows(data))
  } else {
    el_ind <- indiv
  }
  if (length(coor) == 1) {
    coord <- res$ind$coord[el_ind, axes]
  } else {
    coord <- coor
  }
  # browser()

  if (HISTO) {
    x <- numeric()
    y <- numeric()
    ID_O <- numeric()
    minx <- numeric()
    mea <- numeric()
    rmax <- -Inf
    for (i in el_ind) {
      tmp_r <- max(data@M[i, var][[1]]@x) - data@M[i, var][[1]]@x[1] + 1e-100
      mm <- data@M[i, var][[1]]@m
      if (tmp_r > rmax) {
        rmax <- tmp_r
      }
      vv <- .poly_hist(data@M[i, var][[1]])
      x <- c(x, min(vv$x), vv$x, max(vv$x))
      y <- c(y, 0, vv$y, 0)
      ID_O <- c(ID_O, rep(i, length(vv$x) + 2))
      mea <- c(mea, rep(mm, length(vv$x) + 2))
      minx <- c(minx, rep(min(vv$x), length(vv$x) + 2))
    }
  } else {
    # compute densities
    pt <- 300
    x <- numeric()
    y <- numeric()
    ID_O <- numeric()
    minx <- numeric()
    mea <- numeric()
    rmax <- -Inf

    for (i in el_ind) {
      vals <- .qfun(data@M[i, var][[1]], c(0:pt) / pt)
      mm <- data@M[i, var][[1]]@m
      tmp_r <- max(vals) - min(vals)
      if (tmp_r > rmax) {
        rmax <- tmp_r
      }
      vv <- density(vals, n = 128, from = min(vals), to = max(vals))
      x <- c(x, min(vv$x), vv$x, max(vv$x))
      y <- c(y, 0, vv$y, 0)
      ID_O <- c(ID_O, rep(i, length(vv$x) + 2))
      mea <- c(mea, rep(mm, length(vv$x) + 2))
      minx <- c(minx, rep(min(vals), length(vv$x) + 2))
    }

    ### end of compute density
  }
  ID_o <- ID_O
  x_o <- x
  y_o <- y
  ALL_OBJ <- data.frame(x_o = x_o, y_o = y_o, ID_o = ID_o, minx = minx, mea = mea)
  # Automatic scaling (since 1.0.10). The distributions are measured on the
  # scale of the variable, the plane on the scale of the components:
  # - x: a common factor such that the average width of the domains is
  #   strx times the range of the coordinates on the first axis;
  # - y: a common factor such that the highest density is stry times the
  #   range of the coordinates on the second axis (histograms with very narrow
  #   bins are re-binned for drawing, see .poly_hist).
  # Each distribution is centred on its mean. (Up to 1.0.9 the original values
  # were used, and the histograms were drawn as flat lines.)
  coord <- as.matrix(coord)
  rX <- diff(range(coord[, 1]))
  rY <- diff(range(coord[, 2]))
  if (!is.finite(rX) || rX <= 0) rX <- 1
  if (!is.finite(rY) || rY <= 0) rY <- 1
  widths <- tapply(ALL_OBJ$x_o, ALL_OBJ$ID_o, function(v) diff(range(v)))
  avg_width <- mean(widths)
  if (!is.finite(avg_width) || avg_width <= 0) avg_width <- 1
  sx <- strx * rX / avg_width
  limY <- max(ALL_OBJ$y_o[is.finite(ALL_OBJ$y_o)])
  if (!is.finite(limY) || limY <= 0) limY <- 1
  sy <- stry * rY / limY
  ALL_OBJ$y_o <- pmin(ALL_OBJ$y_o, limY)
  cc <- 0
  for (i in el_ind) {
    cc <- cc + 1
    xc <- coord[cc, 1]
    yc <- coord[cc, 2]
    sel <- which(ALL_OBJ$ID_o == i)
    tmp_x <- (ALL_OBJ$x_o[sel] - ALL_OBJ$mea[sel]) * sx + xc
    tmp_y <- ALL_OBJ$y_o[sel] * sy + yc
    ALL_OBJ$x_o[which(ALL_OBJ$ID_o == i)] <- tmp_x
    ALL_OBJ$y_o[which(ALL_OBJ$ID_o == i)] <- tmp_y
  }


  title <- paste(
    "Plot of individuals: distributions for variable",
    colnames(data@M)[var]
  )
  labX <- paste("Axis ", axes[1], " (", format(res$eig[axes[1], 2],
    digits = 2, nsmall = 2
  ), "%)")
  labY <- paste("Axis ", axes[2], " (", format(res$eig[axes[2], 2],
    digits = 2, nsmall = 2
  ), "%)")

  colnames(coord) <- c("x", "y")
  coord <- as.data.frame(coord)
  coord$lab <- rownames(coord)
  p <- ggplot(data = ALL_OBJ, aes(x = .data$x_o, y = .data$y_o)) +
    geom_polygon(aes(group = .data$ID_o, fill = as.factor(.data$ID_o)), alpha = 0.4, colour = "black") +
    geom_point(data = coord, aes(x = .data$x, y = .data$y)) +
    geom_text(
      data = coord,
      aes(x = .data$x, y = .data$y, label = .data$lab),
      hjust = 0.5, vjust = 1
    ) +
    theme_bw() +
    xlab(labX) +
    ylab(labY) +
    ggtitle(title) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)
  ### add a line for mean snd sd
  Ms <- get.MatH.stats(data[el_ind, var], stat = stat)$mat
  # SD=get.MatH.stats(data[el_ind,var],stat="std")$mat
  # Skew=get.MatH.stats(data[el_ind,var],stat="std")$mat
  # kurtH()=get.MatH.stats(data[el_ind,var],stat="std")$mat
  ### print a plot with means

  # the label is the requested statistic (up to 1.0.8 it failed for stat != "mean")
  coord$statval <- format(as.vector(Ms), digits = 1, nsmall = 2)

  p2 <- ggplot(ALL_OBJ, aes(x = .data$x_o, y = .data$y_o)) +
    scale_fill_gradient(low = "yellow", high = "red") +
    geom_polygon(aes(group = .data$ID_o, fill = .data$mea), alpha = 0.4, colour = "black") +
    geom_point(data = coord, aes(x = .data$x, y = .data$y)) +
    theme_bw() +
    geom_text(
      data = coord,
      aes(x = .data$x, y = .data$y, label = .data$statval),
      hjust = 0.5, vjust = 1
    ) +
    xlab(labX) +
    ylab(labY) +
    ggtitle(paste(title, "each label is the ", stat, ".")) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)

  #+ geom_abline(intercept = 0, slope = cor(coord)[1,3]*sd(coord$x),linetype=2)
  # browser()
  print(p)
  print(p2)

  return(invisible(list(p, p2)))
}

# Polygon of a histogram for drawing. If some bins are narrower than 1/100 of
# the domain (e.g. histograms of quantiles with ties, whose density is huge or
# infinite), the histogram is re-binned on max_bins equal-width bins, keeping
# the probabilities (the cdf at the new breaks); otherwise it is drawn as is.
.poly_hist <- function(object, max_bins = 30) {
  x <- object@x
  rng <- diff(range(x))
  if (length(x) > 2 && rng > 0 && min(diff(x)) < rng / 100) {
    br <- seq(x[1], x[length(x)], length.out = max_bins + 1)
    p <- .cdffun(object, br)
    p[1] <- 0
    p[length(p)] <- 1
    object <- .distH(br, p)
  }
  COMP_POLY_H(object)
}

COMP_POLY_H <- function(object) {
  # browser()
  x <- object@x[1]
  y <- 0
  if (length(object@x) < 2) {
    x <- c(x, x, x)
    y <- c(y, 1, 0)
  }
  else {
    for (i in 2:length(object@x)) {
      tmp_dens <- (object@p[i] - object@p[i - 1]) / (object@x[i] - object@x[i - 1] + 1e-100)
      x <- c(x, object@x[i - 1], object@x[i])
      y <- c(y, tmp_dens, tmp_dens)
    }
    x <- c(x, object@x[length(object@x)])
    y <- c(y, 0)
  }

  return(res = data.frame(x = x, y = y))
}
