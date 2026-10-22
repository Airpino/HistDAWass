#' @importFrom ggridges geom_density_ridges
NULL

# Internal helpers for plotting (since v. 1.0.9) ---------------------------
# All the plotting functions draw the plot on the current graphic device and
# return the ggplot object invisibly, so that it can be further modified.
# Data for the plots are built in a vectorised way.

# quantiles of a distribution for a vector of probabilities (in any order)
.qfun <- function(d, p) {
  o <- order(p)
  q <- numeric(length(p))
  # COMP_Q_VECT sorts its third argument in place: pass a sorted copy
  q[o] <- COMP_Q_VECT(d@x, d@p, as.numeric(p[o]) + 0)
  q
}

# cdf of a distribution for a vector of values (vectorised compP)
.cdffun <- function(d, q) {
  x <- d@x
  p <- d@p
  n <- length(x)
  res <- numeric(length(q))
  res[q >= x[n]] <- 1
  inside <- (q > x[1]) & (q < x[n])
  if (any(inside)) {
    qi <- q[inside]
    pos1 <- findInterval(qi, x) # last position with x <= q
    pos2 <- pos1 + 1
    res[inside] <- p[pos1] + (p[pos2] - p[pos1]) * (qi - x[pos1]) / (x[pos2] - x[pos1])
  }
  res
}

# density of the bins, capped for very narrow bins
.bins_density <- function(d) {
  n <- length(d@x)
  lowers <- d@x[-n]
  uppers <- d@x[-1]
  maxdens <- 10 / (max(uppers) - lowers[1])
  dens <- diff(d@p) / (uppers - lowers)
  dens[!is.finite(dens) | dens > maxdens] <- maxdens
  data.frame(xm = lowers, xM = uppers, ym = 0, yM = dens)
}

# draw a ggplot and return it invisibly
.draw <- function(p) {
  print(p)
  invisible(p)
}

#' Summary statistics of a histogram time series
#'
#' Computes basic statistics for each histogram of a \code{HTS} object.
#'
#' @param x A \code{HTS} object.
#' @return A data frame with one row for each histogram: timestamp and period,
#' mean, standard deviation, skewness, kurtosis, minimum, quartiles and maximum.
#' @examples
#' summaryHTS(subsetHTS(RetHTS, from = 1, to = 10))
#' @export
summaryHTS <- function(x) {
  first_or_na <- function(v) if (length(v) > 0) as.numeric(v[[1]]) else NA_real_
  D <- x@data
  Q <- t(vapply(D, function(o) .qfun(o, c(0.25, 0.5, 0.75)), numeric(3)))
  data.frame(
    Tstamp = vapply(D, function(o) first_or_na(o@tstamp), numeric(1)),
    T_period_start = vapply(D, function(o) first_or_na(o@period$start), numeric(1)),
    T_period_end = vapply(D, function(o) first_or_na(o@period$end), numeric(1)),
    mean = vapply(D, function(o) o@m, numeric(1)),
    std = vapply(D, function(o) o@s, numeric(1)),
    skew = vapply(D, skewH, numeric(1)),
    kurt = vapply(D, kurtH, numeric(1)),
    min = vapply(D, function(o) o@x[1], numeric(1)),
    Q1 = Q[, 1],
    MED = Q[, 2],
    Q3 = Q[, 3],
    Max = vapply(D, function(o) o@x[length(o@x)], numeric(1))
  )
}

# plot of a single distribution ------------------------------------------
.plot_distributionH <- function(x, type = "HISTO", col = "green", border = "black") {
  xlabel <- paste("m= ", format(x@m, nsmall = 4), " std= ", format(x@s, nsmall = 4))
  type <- toupper(type)
  p <- switch(type,
    "HISTO" = {
      df <- .bins_density(x)
      ggplot(df, aes(xmin = .data$xm, xmax = .data$xM, ymin = .data$ym, ymax = .data$yM)) +
        geom_rect(alpha = 0.7, fill = col, colour = border) +
        geom_vline(xintercept = x@m, colour = col, linetype = "dashed", linewidth = 0.5) +
        xlab(xlabel) +
        ylab("density") +
        ggtitle("Histogram")
    },
    "CDF" = {
      df <- data.frame(x = x@x, cdf = x@p)
      ggplot(df, aes(x = .data$x, y = .data$cdf)) +
        geom_line(colour = border) +
        xlab("x") +
        ylab("p") +
        ggtitle("Cumulative Distribution Function")
    },
    "QF" = {
      df <- data.frame(x = x@x, cdf = x@p)
      ggplot(df, aes(x = .data$cdf, y = .data$x)) +
        geom_line(colour = border) +
        xlab("p") +
        ylab("x") +
        ggtitle("Quantile Function")
    },
    "DENS" = {
      # 201 points according to the quantile function
      df <- data.frame(x = .qfun(x, c(0:200) / 200))
      ggplot(df, aes(x = .data$x)) +
        geom_density(alpha = 0.7, fill = col, colour = border) +
        geom_vline(xintercept = x@m, colour = "black", linetype = "dashed", linewidth = 0.5) +
        xlab(xlabel) +
        ylab("density") +
        ggtitle("Density plot (KDE)")
    },
    "HBOXPLOT" = ,
    "VBOXPLOT" = {
      q <- .qfun(x, c(0, 0.25, 0.5, 0.75, 1))
      df <- data.frame(
        x = factor(0), ymin = q[1], lower = q[2], middle = q[3],
        upper = q[4], ymax = q[5]
      )
      pp <- ggplot(df, aes(
        x = .data$x, ymin = .data$ymin, lower = .data$lower,
        middle = .data$middle, upper = .data$upper, ymax = .data$ymax
      )) +
        geom_boxplot(stat = "identity", fill = col, linewidth = 1, colour = border)
      if (type == "HBOXPLOT") {
        pp + ylab(xlabel) + xlab("") + ggtitle("Horizontal Boxplot") + coord_flip()
      } else {
        pp + ylab("x") + xlab(xlabel) + ggtitle("Vertical Boxplot")
      }
    },
    stop("type must be one of HISTO, CDF, QF, DENS, HBOXPLOT, VBOXPLOT")
  )
  return(p)
}

# data for the ridgeline plots of a MatH ---------------------------------
# quantiles of each distribution (DENS)
.joy_dens_data <- function(DATA, vars) {
  vp <- c(0:1000) / 1000
  rn <- get.MatH.rownames(DATA)
  vn <- get.MatH.varnames(DATA)
  pieces <- lapply(vars, function(v) {
    q <- unlist(lapply(seq_along(rn), function(i) .qfun(DATA@M[i, v][[1]], vp)))
    data.frame(x = q, y = rep(rn, each = length(vp)), z = vn[v])
  })
  DF <- do.call(rbind, pieces)
  DF$y <- factor(DF$y, levels = rev(rn))
  DF
}
# profile of the histograms (HISTO)
.joy_histo_data <- function(DATA, vars, tr = 1) {
  rn <- get.MatH.rownames(DATA)
  vn <- get.MatH.varnames(DATA)
  pieces <- lapply(vars, function(v) {
    one <- lapply(seq_along(rn), function(i) {
      d <- DATA@M[i, v][[1]]
      n <- length(d@x)
      widths <- diff(d@x)
      widths[widths == 0] <- 1e-16
      densit <- diff(d@p) / widths
      data.frame(
        X = as.vector(rbind(d@x[-n], d@x[-n], d@x[-1], d@x[-1])),
        Y = i,
        H = as.vector(rbind(0, densit, densit, 0)),
        y = rn[i]
      )
    })
    DF <- do.call(rbind, one)
    if (tr == 1) {
      treshold <- 2
    } else {
      treshold <- 1.5 * diff(stats::quantile(DF$H, probs = c(0.25, 0.75))) +
        stats::quantile(DF$H, probs = 0.75)
    }
    DF$H[DF$H > treshold] <- treshold
    DF$z <- vn[v]
    DF
  })
  DF <- do.call(rbind, pieces)
  DF$z <- factor(DF$z, levels = vn[vars])
  DF$y <- factor(DF$y, levels = rev(rn))
  DF
}

# plot of a MatH ---------------------------------------------------------
.plot_MatH <- function(x, type = "HISTO", border = "black", angL = 330) {
  varsno <- ncol(x@M)
  indno <- nrow(x@M)
  type <- toupper(type)
  if (varsno == 1 && indno == 1) {
    return(.plot_distributionH(x@M[1, 1][[1]], type = type, col = "green", border = border))
  }
  ridge_theme <- theme(
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(face = "bold", vjust = 0),
    legend.position = "none",
    axis.title.x = element_blank(), axis.title.y = element_blank()
  )
  p <- switch(type,
    "HISTO" = {
      df <- .joy_histo_data(x, seq_len(varsno), tr = 1)
      ggplot(df, aes(
        x = .data$X, y = .data$y, height = .data$H,
        group = .data$y, fill = .data$z
      )) +
        geom_density_ridges(stat = "identity", scale = 1.1, alpha = 0.7, color = border) +
        facet_wrap(~z, scales = "free_x", ncol = varsno) +
        ridge_theme
    },
    "DENS" = {
      df <- .joy_dens_data(x, seq_len(varsno))
      ggplot(df, aes(x = .data$x, y = .data$y, fill = .data$z)) +
        geom_density_ridges(scale = 1.3, rel_min_height = 0.01, alpha = 0.7) +
        facet_wrap(~z, scales = "free_x", ncol = varsno) +
        ridge_theme
    },
    "BOXPLOT" = {
      cells <- expand.grid(i = seq_len(indno), j = seq_len(varsno))
      Q <- t(mapply(function(i, j) .qfun(x@M[i, j][[1]], c(0, 0.25, 0.5, 0.75, 1)),
        cells$i, cells$j
      ))
      df <- data.frame(
        ymin = Q[, 1], yQ1 = Q[, 2], yME = Q[, 3], yQ3 = Q[, 4], ymax = Q[, 5],
        indiv = factor(rownames(x@M)[cells$i], levels = rownames(x@M)),
        vars = factor(colnames(x@M)[cells$j], levels = colnames(x@M))
      )
      ggplot(df, aes(
        x = factor(0), ymin = .data$ymin, lower = .data$yQ1, middle = .data$yME,
        upper = .data$yQ3, ymax = .data$ymax, fill = .data$vars
      )) +
        geom_boxplot(stat = "identity", linewidth = 0.5, colour = border) +
        facet_grid(vars ~ indiv, scales = "free_y") +
        scale_x_discrete(breaks = NULL) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 12, face = "bold", angle = 30),
          strip.text.y = element_text(size = 12, face = "bold")
        )
    },
    stop("type must be one of HISTO, DENS, BOXPLOT")
  )
  return(p)
}

# plot of a histogram time series -----------------------------------------
# The timestamps are split in blocks of maxno.perplot elements, each block is
# a row (a facet) of the plot.
.plot_HTS <- function(x, type = "BOXPLOT", border = "black", maxno.perplot = 15) {
  type <- toupper(type)
  numofp <- length(x@data)
  df1 <- summaryHTS(x)
  df1$block <- factor((seq_len(numofp) - 1) %/% maxno.perplot + 1)
  nblocks <- nlevels(df1$block)
  hts_theme <- theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )
  p <- switch(type,
    "BOXPLOT" = {
      ggplot(df1, aes(
        x = factor(.data$Tstamp), ymin = .data$min, lower = .data$Q1, middle = .data$MED,
        upper = .data$Q3, ymax = .data$Max, fill = factor(.data$Tstamp)
      )) +
        geom_boxplot(stat = "identity", colour = border) +
        facet_wrap(~block, ncol = 1, scales = "free_x") +
        hts_theme +
        theme(axis.text.x = element_text(vjust = 0.5, size = 8))
    },
    "VIOLIN" = {
      TT <- max(df1$Tstamp) - min(df1$Tstamp)
      Tini <- min(df1$Tstamp)
      grad <- if (TT > 0) (df1$Tstamp - Tini) / TT else rep(0.5, numofp)
      # 100 points from the quantile function of each distribution
      vp <- c(0:99) / 100
      df <- data.frame(
        x = unlist(lapply(x@data, .qfun, p = vp)),
        t = rep(df1$Tstamp, each = length(vp)),
        grad = rep(grad, each = length(vp)),
        block = rep(df1$block, each = length(vp))
      )
      df2 <- data.frame(x = df1$mean, t = df1$Tstamp, grad = grad, block = df1$block)
      ggplot(df, aes(factor(.data$t), .data$x, fill = .data$grad)) +
        geom_violin(colour = border) +
        geom_point(data = df2, aes(factor(.data$t), .data$x)) +
        geom_line(
          data = df2, aes(x = factor(.data$t), y = .data$x, group = .data$block),
          linewidth = 0.5, alpha = 0.7, linetype = "dashed"
        ) +
        scale_fill_gradient2(limits = c(0, 1), low = "red", mid = "white", high = "green", midpoint = 0.5) +
        facet_wrap(~block, ncol = 1, scales = "free_x") +
        hts_theme +
        theme(
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          axis.text.x = element_text(angle = 330, vjust = 0.5, size = 10)
        )
    },
    stop("type must be one of BOXPLOT, VIOLIN")
  )
  return(p)
}

#' Plot observed versus predicted histograms
#'
#' Compares observed and predicted histograms (e.g. from
#' \code{\link{WH.regression.two.components.predict}}), one panel for each
#' individual.
#'
#' @param PRED A \code{MatH} object with one column, the predicted histograms.
#' @param OBS A \code{MatH} object with one column, the observed histograms.
#' @param type The type of plot: \code{"HISTO"} (default, histograms), \code{"CDF"}
#' (cumulative distribution functions) or \code{"DENS"} (kernel density
#' approximations).
#' @param ncolu The number of columns of panels (default 2).
#' @return The plot is drawn on the current device and the \code{ggplot} object
#' is returned invisibly.
#' @seealso \code{\link{plot_errors}}
#' @examples
#' pars <- WH.regression.two.components(BLOOD, Yvar = 1, Xvars = c(2:3))
#' PRED <- WH.regression.two.components.predict(data = BLOOD[, 2:3], parameters = pars)
#' OBS <- BLOOD[, 1]
#' plotPredVsObs(PRED, OBS, "HISTO")
#' \donttest{
#' plotPredVsObs(PRED, OBS, "CDF")
#' plotPredVsObs(PRED, OBS, "DENS")
#' }
#' @export
plotPredVsObs <- function(PRED, OBS, type = "HISTO", ncolu = 2) {
  type <- toupper(type)
  nobj <- get.MatH.nrows(PRED)
  names_o <- rownames(OBS@M)
  pieces <- lapply(seq_len(nobj), function(ind) {
    o <- OBS@M[ind, 1][[1]]
    pr <- PRED@M[ind, 1][[1]]
    two <- switch(type,
      "HISTO" = list(OBS = .bins_density(o), PRED = .bins_density(pr)),
      "CDF" = list(
        OBS = data.frame(x = o@x, cdf = o@p),
        PRED = data.frame(x = pr@x, cdf = pr@p)
      ),
      "DENS" = {
        vp <- c(0:199) / 200
        list(OBS = data.frame(x = .qfun(o, vp)), PRED = data.frame(x = .qfun(pr, vp)))
      },
      stop("type must be one of HISTO, CDF, DENS")
    )
    two$OBS$Type <- "OBS"
    two$PRED$Type <- "PRED"
    res <- rbind(two$OBS, two$PRED)
    res$NameO <- names_o[ind]
    res
  })
  M <- do.call(rbind, pieces)
  M$Type <- factor(M$Type, levels = c("OBS", "PRED"))
  M$NameO <- factor(M$NameO, levels = names_o)
  bottom <- theme(legend.position = "bottom", legend.title = element_blank())
  p <- switch(type,
    "HISTO" = ggplot(M, aes(
      xmin = .data$xm, xmax = .data$xM, ymin = .data$ym,
      ymax = .data$yM, fill = .data$Type
    )) +
      geom_rect(alpha = 0.5) +
      bottom,
    "CDF" = ggplot(M, aes(x = .data$x, y = .data$cdf)) +
      geom_line(aes(linetype = .data$Type)) +
      bottom + theme(axis.title.x = element_blank()),
    "DENS" = ggplot(M, aes(x = .data$x, fill = .data$Type)) +
      geom_density(alpha = 0.6, colour = "gray") +
      bottom + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  )
  p <- p + facet_wrap(~NameO, ncol = ncolu, dir = "v", strip.position = "right")
  .draw(p)
}

#' Plot the errors of predicted histograms
#'
#' Plots the differences between observed and predicted histograms (e.g. from
#' \code{\link{WH.regression.two.components.predict}}) and prints the root mean
#' squared error based on the L2 Wasserstein distance.
#'
#' @param PRED A \code{MatH} object with one column, the predicted histograms.
#' @param OBS A \code{MatH} object with one column, the observed histograms.
#' @param type The type of plot: \code{"HISTO_QUA"} (default, differences of the
#' quantile functions), \code{"HISTO_DEN"} (differences of the densities) or
#' \code{"DENS_KDE"} (differences of kernel density approximations).
#' @param np The number of points used to compute the differences (default 200).
#' @return The root mean squared error, invisibly. The plot is drawn on the
#' current device.
#' @seealso \code{\link{plotPredVsObs}}, \code{\link{WH.regression.GOF}}
#' @examples
#' pars <- WH.regression.two.components(BLOOD, Yvar = 1, Xvars = c(2:3))
#' PRED <- WH.regression.two.components.predict(data = BLOOD[, 2:3], parameters = pars)
#' OBS <- BLOOD[, 1]
#' plot_errors(PRED, OBS, "HISTO_QUA")
#' plot_errors(PRED, OBS, "HISTO_DEN")
#' plot_errors(PRED, OBS, "DENS_KDE")
#' @export
plot_errors <- function(PRED, OBS, type = "HISTO_QUA", np = 200) {
  type <- toupper(type)
  nobj <- get.MatH.nrows(PRED)
  minimum <- min(min((get.MatH.stats(OBS, stat = "min"))$mat), min((get.MatH.stats(PRED, stat = "min"))$mat))
  maximum <- max(max((get.MatH.stats(OBS, stat = "max"))$mat), max((get.MatH.stats(PRED, stat = "max"))$mat))
  vals <- minimum + (maximum - minimum) * c(0:np) / np
  NM <- WH.bind.col(OBS, PRED)
  colnames(NM@M) <- c("OBS", "PRED")
  NM <- registerMH(NM)
  names_o <- rownames(OBS@M)
  dists <- vapply(seq_len(nobj), function(i) WassSqDistH(NM@M[i, 1][[1]], NM@M[i, 2][[1]]), numeric(1))
  RMS_W <- sqrt(mean(dists))
  pieces <- lapply(seq_len(nobj), function(i) {
    xo <- NM@M[i, 1][[1]]
    xp <- NM@M[i, 2][[1]]
    df <- switch(type,
      "HISTO_QUA" = {
        p <- c(0:np) / np
        data.frame(quantile = .qfun(xo, p) - .qfun(xp, p), p = p)
      },
      "HISTO_DEN" = {
        v <- vals[1:np]
        mp <- 0.5 * (v[-1] + v[-np])
        ran <- diff(v)
        data.frame(x = mp, density = diff(.cdffun(xo, v)) / ran - diff(.cdffun(xp, v)) / ran)
      },
      "DENS_KDE" = {
        p <- c(1:np) / np
        do <- stats::density(.qfun(xo, p), from = minimum, to = maximum)
        dp <- stats::density(.qfun(xp, p), from = minimum, to = maximum)
        data.frame(x = do$x, density = do$y - dp$y)
      },
      stop("type must be one of HISTO_QUA, HISTO_DEN, DENS_KDE")
    )
    df$NameO <- names_o[i]
    df
  })
  M <- do.call(rbind, pieces)
  M$NameO <- factor(M$NameO, levels = names_o)
  p <- switch(type,
    "HISTO_QUA" = ggplot(M, aes(x = .data$p, y = .data$quantile, group = .data$NameO, colour = .data$NameO)) +
      geom_line() +
      geom_point() +
      ggtitle("Plot of quantile differences"),
    "HISTO_DEN" = ggplot(M, aes(x = .data$x, y = .data$density, group = .data$NameO, colour = .data$NameO)) +
      geom_line() +
      ggtitle("Plot of density differences"),
    "DENS_KDE" = ggplot(M, aes(x = .data$x, y = .data$density, group = .data$NameO, colour = .data$NameO)) +
      geom_line() +
      ggtitle("Plot of smoothed (KDE) \n density differences")
  )
  print(p)
  print(RMS_W)
  invisible(RMS_W)
}
