# Batch SOM Kohonen Maps of HD using adaptive distances -----
#' Batch Kohonen self-organising maps of histogram data with adaptive distances
#'
#' Batch Kohonen self-organising 2d map of the individuals of a \code{MatH}
#' object, based on adaptive L2 Wasserstein distances.
#'
#' @param x A \code{MatH} object.
#' @param net A list describing the topology of the map:
#' \code{list(xdim = <rows>, ydim = <columns>, topo = "rectangular" or "hexagonal")},
#' see \code{\link[class]{somgrid}}. Default: a 4 x 3 rectangular map.
#' @param kern.param The parameter of the RBF kernel (default 2).
#' @param TMAX,Tmin The initial and final values of the neighbourhood radius.
#' When negative (default), they are computed from the size of the map.
#' @param niter The maximum number of iterations (default 30).
#' @param repetitions The number of repetitions of the algorithm (default 5);
#' the best solution is returned.
#' @inheritParams WH_kmeans
#' @param schema An integer, the system of weights: 1 = a weight for each
#' variable; 2 = a weight for the position and one for the variability component
#' of each variable; 3 and 4 = as 1 and 2, with a different set of weights for
#' each neuron (product of the weights of each variable equal to one); 5 and 6 =
#' as 1 and 2, with a different set of weights for each neuron (product of all
#' the weights equal to one, default 6).
#' @param init.weights The initial weights: \code{"EQUAL"} (default),
#' \code{"PROP"} or \code{"RANDOM"}.
#' @param weight.sys The constraint on the weights: \code{"PROD"} (default) or
#' \code{"SUM"}.
#' @param theta A number, the exponent of the weights when
#' \code{weight.sys = "SUM"} (default 2).
#' @param Wfix A logical value (default \code{FALSE}). If \code{TRUE} the weights
#' are not updated (non adaptive distances).
#' @param atleast An integer (default 2): solutions where fewer than
#' \code{atleast} neurons attract individuals are discarded.
#'
#' @return A list with:
#' \item{solution}{The best solution among the repetitions: a list with
#' \code{MAP} (the grid of the map), \code{IDX} (the neuron of each individual),
#' \code{cardinality} (the number of individuals for each neuron), \code{proto}
#' (a \code{MatH} object with the prototypes of the neurons), \code{Crit} (the
#' criterion), \code{weights.comp} (the weights of the components of the
#' variables) and \code{Weight.sys} (the system of weights).}
#' \item{quality}{The proportion of the total sum of squares explained by the
#' map (the higher the better).}
#' and further details on the parameters of the algorithm.
#' The final map is plotted on the current graphic device.
#' @references Irpino A., Verde R., De Carvalho F.A.T. (2012). Batch self
#' organizing maps for interval and histogram data. In: Proceedings of COMPSTAT
#' 2012, ISI/IASC, 143-154.
#' @details Batch Self-Organising Maps (BSOM) for histogram data, with a cost
#' function based on the L2 Wasserstein distance, i.e. the Euclidean distance
#' between quantile functions. The adaptive version introduces an automatic
#' system of weights in the cost function, to take into account the different
#' effect of the variables on the map.
#'
#' @seealso \code{\link{WH_2d_Kohonen_maps}}
#' @importFrom class somgrid
#' @importFrom stats lm dist runif
#' @examples
#' \donttest{
#' results <- WH_2d_Adaptive_Kohonen_maps(
#'   x = BLOOD,
#'   net = list(xdim = 2, ydim = 3, topo = c("rectangular")),
#'   repetitions = 2, simplify = TRUE,
#'   qua = 10, standardize = TRUE
#' )
#' }
#' @export
WH_2d_Adaptive_Kohonen_maps <- function(x, net = list(xdim = 4, ydim = 3, topo = c("rectangular")),
                                        kern.param = 2, TMAX = -9999, Tmin = -9999,
                                        niter = 30, repetitions,
                                        simplify = FALSE,
                                        qua = 10,
                                        standardize = FALSE, schema = 6,
                                        init.weights = "EQUAL", weight.sys = "PROD",
                                        theta = 2, Wfix = FALSE, verbose = FALSE, atleast = 2) {
  tol <- 1e-10
  if (weight.sys == "PROD") {
    weight_sys <- 1
    theta <- 1
  }
  if (weight.sys == "SUM") {
    weight_sys <- 2
  }
  # remember to fix TMAX e Tmin
  # require(class)#for somgrid function
  # Check initial conditions and passed parameters

  ind <- nrow(x@M)
  vars <- ncol(x@M)
  ## we homogeneize data for speeding up the code if required
  tmp <- Prepare(x, simplify, qua, standardize)

  MM <- tmp$MM
  x <- tmp$x
  remove(tmp)
  # end of the preprocessing step

  TOTSSQ <- 0
  SSQ_comp <- rep(0, vars * 2)
  tmpc <- as.matrix(get.MatH.stats(x)$mat)

  colm <- colMeans(tmpc)
  for (v in 1:vars) {
    tmpSSQ <- c_ComputeFASTSSQ(MM[[v]])$SSQ
    TOTSSQ <- TOTSSQ + tmpSSQ
    C_comp <- sum((tmpc[, v] - colm[v])^2)
    V_comp <- tmpSSQ - C_comp
    SSQ_comp[v * 2 - 1] <- C_comp
    SSQ_comp[v * 2] <- V_comp
  }
  ## ------------------------------------------------batchKOHONENMAP-----------------------------------------------------
  solutions <- list()
  criteria <- numeric(0)
  MAP <- somgrid(net$xdim, net$ydim, net$topo)
  k <- nrow(MAP$pts)
  if (missing(repetitions)) {
    repetitions <- 5
  }

  restart <- 1
  repet <- 0


  while (repet < repetitions) {
    repet <- repet + 1
    # random selection of prototypes of neurons
    if (k <= ind) {
      proto <- x[sample(1L:ind, k, replace = FALSE), ]
    } else {
      proto <- x[sample(1L:ind, k, replace = TRUE), ]
    }

    rownames(proto@M) <- paste0("Prot. ", 1:k)
    nhbrdist <- as.matrix(dist(MAP$pts))
    dmax <- (max(MAP$pts[, 1]) - min(MAP$pts[, 1]))^2 + (max(MAP$pts[, 2]) - min(MAP$pts[, 2]))^2
    dmax <- (max(as.matrix(dist(MAP$pts))))^2

    if (TMAX < 0) {
      TMAX <- sqrt(-(dmax / 4) / (2 * log(0.99)))
    }
    if (Tmin < 0) {
      Tmin <- sqrt(-(1 / 4) / (2 * log(0.01)))
    }

    TT <- TMAX
    KT <- exp(-(nhbrdist^2) / (2 * TT^2))
    KT[which(KT < 1e-10)] <- 0
    # initialize weights
    if (Wfix) {
      lambdas <- matrix(1, 2 * vars, k)
    } else {
      if (init.weights == "EQUAL") {
        cat("Weights initialization === EQUAL  \n")
        if (weight.sys == "PROD") {
          lambdas <- matrix(1, 2 * vars, k)
        }
        else {
          if (schema == 1 || schema == 2 || schema == 3 || schema == 4) {
            lambdas <- matrix(1 / vars, 2 * vars, k)
          }
          if (schema == 5 || schema == 6) {
            lambdas <- matrix(1 / (2 * vars), 2 * vars, k)
          }
        }
      }
      if (init.weights == "PROP") {
        if (weight.sys == "PROD") {
          # 1=A weight for each variable (default)
          # 2=A weight for the average and the dispersion component of each variable
          # 3=Same as 1 but a different set of weights for each cluster (indipprod=1)
          # 4=Same as 2 but a different set of weights for each cluster (indipprod=1)
          # 5=Same as 1 but a different set of weights for each cluster (allprod=1)
          # 6=Same as 2 but a different set of weights for each cluster (allprod=1)
          if (schema == 1 || schema == 3 || schema == 5) {
            tmpvals <- SSQ_comp[c(1:v) * 2 - 1] + SSQ_comp[c(1:v) * 2]
            mg <- (prod(tmpvals))^(1 / v)
            tmpw <- mg / tmpvals
            # one weight for each variable, repeated for its two components
            # (up to 1.0.9 rep(tmpw, 2, each = 2): twice the needed length,
            # unnoticed because this branch was then overwritten)
            elew <- rep(tmpw, each = 2)
            lambdas <- matrix(rep(elew, k), nrow = 2 * vars, ncol = k, byrow = FALSE)
          }
          if (schema == 2 || schema == 4) { # (up to 1.0.9: schema == 2 || 4, always TRUE)
            tmpvalsc <- SSQ_comp[c(1:v) * 2 - 1]
            tmpvalsv <- SSQ_comp[c(1:v) * 2]
            mgc <- (prod(tmpvalsc))^(1 / v)
            mgv <- (prod(tmpvalsv))^(1 / v)

            tmpwc <- mgc / tmpvalsc
            tmpwv <- mgv / tmpvalsv
            elew <- rep(0, 2 * v)
            elew[c(1:v) * 2 - 1] <- tmpwc
            elew[c(1:v) * 2] <- tmpwv

            lambdas <- matrix(rep(elew, k), nrow = 2 * vars, ncol = k, byrow = FALSE)
          }
          if (schema == 6) {
            tmpvals <- SSQ_comp
            mg <- (prod(tmpvals))^(1 / (2 * v))
            tmpw <- mg / tmpvals
            lambdas <- matrix(rep(tmpw, k), nrow = 2 * vars, ncol = k, byrow = FALSE)
          }
        } else {
          if ((schema == 1 || schema == 2 || schema == 3 || schema == 4)) {
            lambdas <- matrix(1 / v, nrow = 2 * vars, ncol = k, byrow = FALSE)
          }
          if (schema == 5 || schema == 6) {
            lambdas <- matrix(2 / (2 * v), nrow = 2 * vars, ncol = k, byrow = FALSE)
          }
        }
      }
      if ((init.weights != "EQUAL") & (init.weights != "PROP")) { # Random initialization
        cat("Weights initialization === RANDOM  \n")
        m1 <- matrix(runif((vars * k), 0.01, 0.99), vars, k)
        m2 <- matrix(runif((vars * k), 0.01, 0.99), vars, k)
        m1 <- m1 / matrix(rep(apply(m1, 2, sum)), vars, k, byrow = TRUE)
        m2 <- m2 / matrix(rep(apply(m2, 2, sum)), vars, k, byrow = TRUE)
        if (weight.sys == "PROD") {
          m1 <- exp(m1 * vars - 1)
          m2 <- exp(m2 * vars - 1)
        }

        if (schema == 1) {
          m1 <- matrix(m1[, 1], nrow = vars, ncol = k)
          m2 <- m1
        }
        if (schema == 2) {
          m1 <- matrix(rep(m1[, 1], k), vars, k)
          m2 <- matrix(rep(m2[, 1], k), vars, k)
        }
        if (schema == 3) {
          m2 <- m1
        }

        lambdas <- matrix(0, 2 * vars, k)
        colnames(lambdas) <- paste("Clust", c(1:k), sep = "_")

        n1 <- paste("M_Var.", c(1:vars), sep = "_")
        n2 <- paste("C_Var.", c(1:vars), sep = "_")
        nr <- list()
        for (nn in 1:vars) {
          nr[[nn * 2 - 1]] <- n1[[nn]]
          nr[[nn * 2]] <- n2[[nn]]
        }
        rownames(lambdas) <- nr
        lambdas[(c(1:vars) * 2 - 1), ] <- m1
        lambdas[(c(1:vars) * 2), ] <- m2
      }
    }



    # assign objects to closest neuron
    # compute adaptive distances to prototypes
    distances <- array(0, dim = c(vars, k, 2))
    diINDtoPROT <- array(0, dim = c(ind, vars, k, 2))
    # compute distances

    MD <- (c_ComputeFast_L2_SQ_WASS_DMAT(MM, proto))$Dist
    KT[KT < 1e-50] <- 0
    Dind2Neu <- MD %*% KT

    # initialize matrix of memberships
    IDX <- max.col(-Dind2Neu, ties.method = "first")
    memb <- matrix(0, ind, k)
    for (i in (1:ind)) {
      memb[i, IDX[i]] <- 1
    }
    cat(paste("---------> repetitions  ", repet, "\n"))
    ## initialize clusters and prototypes

    SSQ1 <- matrix(0, k, vars)
    GenCrit <- Inf
    ##  compute initial criterion
    tmpIDX <- Table2(IDX, k, FALSE)
    SSQ <- sum(memb * Dind2Neu)
    GenCrit <- SSQ
    OK <- 1
    maxit <- 100

    t <- 0
    ft <- 10
    fts <- 0
    while ((TT > Tmin) || (fts <= ft)) {
      if (t < (niter - 1)) {
        t <- t + 1
        TT <- TMAX * (Tmin / TMAX)^(t / (niter - 1))
        KT <- exp(-(nhbrdist^2) / (2 * TT^2))
        KT[which(KT < 1e-20)] <- 0
      }
      if ((t == (niter - 1)) & (fts <= ft)) {
        fts <- fts + 1
      }
      # computing prototypes
      CardinalityOfNeuron <- colSums(memb)
      proto <- c_PROTO_KOHONEN(proto, k, ind, MM, vars, KT, IDX)
      #########################################################################################
      # first compute distances using kernel
      #  computing distances

      ###########    NEW part
      tmpD <- c_ComputeFast_L2_SQ_WASS_DMAT(MM, proto)

      GRKT <- KT[IDX, 1:k]
      distances <- array(0, dim = c(vars, k, 2))
      for (j in 1:vars) {
        distances[j, , 1] <- distances[j, , 1] + colSums(tmpD$DET[[j]]$DM * GRKT)
        distances[j, , 2] <- distances[j, , 2] + colSums(tmpD$DET[[j]]$DV * GRKT)
      }

      tmp2 <- matrix(0, ind, k)
      if (Wfix == FALSE) { # Weights computation
        # S2.2) Weights computation
        # NEW
        if (t > 0) {
          oldlambda <- lambdas
        }

        lambdas <- c_ADA_F_WHEIGHT(
          list(as.matrix(distances[, , 1]), as.matrix(distances[, , 2])),
          k, vars, ind, schema, weight_sys, theta, 1
        )
        if (t > 0) {
          for (cluster in 1:k) {
            if (weight.sys == "PROD") {
              if ((schema >= 1 & schema <= 4) &
                ((abs(prod(lambdas[, cluster]) - 1)) > 0.0001)) {
                print(cat("ini ", prod(lambdas[, cluster])))
                lambdas[, cluster] <- oldlambda[, cluster]
                print(cat("t", t, "CLU", cluster, "new", prod(lambdas[, cluster])))
              }
              if ((schema == 5 || schema == 6) &
                ((abs(prod(lambdas[, cluster]) - 1)) > 0.0001)) {
                print(cat("ini ", prod(lambdas[, cluster])))
                lambdas[, cluster] <- oldlambda[, cluster]
                print(cat("t", t, "CLU", cluster, "new", prod(lambdas[, cluster])))
              }
            }
            if (weight.sys == "SUM") {
              if ((schema >= 1 & schema <= 4) &
                ((abs(sum(lambdas[, cluster]) - 2)) > 0.0001)) {
                print(cat("ini ", sum(lambdas[, cluster])))
                lambdas[, cluster] <- oldlambda[, cluster]
                print(cat("t", t, "CLU", cluster, "new", sum(lambdas[, cluster])))

                # browser()
              }
              if ((schema == 5 || schema == 6) &
                ((abs(sum(lambdas[, cluster]) - 1)) > 0.0001)) {
                print(cat("ini ", sum(lambdas[, cluster])))
                lambdas[, cluster] <- oldlambda[, cluster]
                print(cat("t", t, "CLU", cluster, "new", sum(lambdas[, cluster])))
              }
            }
          }
        }
        wM <- (lambdas[((1:vars) * 2 - 1), ])^theta
        wV <- (lambdas[(1:vars) * 2, ])^theta


        # assign data to neurons
        ####################
        if (schema == 1 | schema == 3) {
          ####################
          for (variables in (1:vars)) {
            tmpM <- matrix(wM[variables, ], ind, k, byrow = TRUE)
            tmp2 <- tmp2 + tmpD$DET[[variables]]$D * tmpM
          }
          ###############
        } else {
          for (variables in (1:vars)) {
            tmpM <- matrix(wM[variables, ], ind, k, byrow = TRUE)
            tmp2 <- tmp2 + tmpD$DET[[variables]]$DM * tmpM
            tmpM <- matrix(wV[variables, ], ind, k, byrow = TRUE)
            tmp2 <- tmp2 + tmpD$DET[[variables]]$DV * tmpM
          }
        }
      }
      if (Wfix == TRUE) {

        ####################
        if (schema == 1 | schema == 3) { # one weight for one variable
          for (variables in (1:vars)) {
            tmp2 <- tmp2 + tmpD$DET[[variables]]$D
          }
        } else {
          for (variables in (1:vars)) {
            tmp2 <- tmp2 + tmpD$DET[[variables]]$DM
            tmp2 <- tmp2 + tmpD$DET[[variables]]$DV
          }
        }
      }
      remove(tmpD)
      Dind2Neu <- tmp2 %*% KT # matrix(Inf,ind,k)

      # recompute criterion
      IDX <- max.col(-Dind2Neu, ties.method = "first")
      old_NUM <- tmpIDX
      tmpIDX <- Table2(IDX, k, FALSE)
      OldCrit <- GenCrit
      GenCrit <- 0
      for (individual in 1:ind) {
        GenCrit <- GenCrit + sum(Dind2Neu[individual, IDX[individual]])
      }
      if (verbose) {
        cat(paste(t, GenCrit, "\n", sep = "---->"))
      }
      if (is.na(GenCrit)) {
        print("A GREAT PROBLEM HERE!")
      }
      if (verbose) {
        print(print_map(MAP, IDX))
      }
    }
    # crisp assignment

    cardinality <- Table2(IDX, k, FALSE)
    names(IDX) <- get.MatH.rownames(x)
    colnames(lambdas) <- get.MatH.rownames(proto)
    rownames(lambdas) <- paste(rep(get.MatH.varnames(x), each = 2),c(" (Cen.)", " (Var.)"), sep = ".")

    if (length(cardinality) >= atleast) {
      solutions <- c(solutions, list(solution = list(
        MAP = MAP, IDX = IDX, cardinality = cardinality, proto = proto,
        Crit = GenCrit, weights.comp = lambdas, Weight.sys = weight.sys
      )))
      criteria <- c(criteria, GenCrit)
    } else {
      repet <- repet - 1
    }
  }

  FM <- solutions[[which.min(criteria)]]$MAP
  FIDX <- solutions[[which.min(criteria)]]$IDX
  print(print_map(FM, FIDX, gr = TRUE))

  best.solution <- list(
    solution = solutions[[which.min(criteria)]],
    repetitions = which.min(criteria),
    quality = 1 - min(criteria) / TOTSSQ, dmax = dmax, TMAX = TMAX, Tmin = Tmin, DATA = x
  )
  memb <- matrix(0, ind, k)
  memb[(best.solution$solution$IDX - 1) * ind + c(1:ind)] <- 1
  resTOTSQ <- c_WH_ADPT_KMEANS_TOTALSSQ(x, memb,
    m = 1,
    lambdas = best.solution$solution$weights.comp,
    proto = best.solution$solution$proto
  )
  TTSQ2 <- sum(resTOTSQ$TSQ2_m + resTOTSQ$TSQ2_v)
  best.solution$quality2 <- 1 - min(criteria) / TTSQ2
  return(best.solution)
}

# Batch SOM of HD -----
#' Batch Kohonen self-organising maps of histogram data
#'
#' Batch Kohonen self-organising 2d map of the individuals of a \code{MatH}
#' object, based on the L2 Wasserstein distance.
#'
#' @inheritParams WH_2d_Adaptive_Kohonen_maps
#' @param TMAX,Tmin The initial and final values of the neighbourhood radius
#' (defaults 2 and 0.2).
#'
#' @return A list with:
#' \item{solution}{The best solution among the repetitions: a list with
#' \code{MAP} (the grid of the map), \code{IDX} (the neuron of each individual),
#' \code{cardinality} (the number of individuals for each neuron), \code{proto}
#' (a \code{MatH} object with the prototypes of the neurons), \code{Crit} (the
#' criterion), \code{weights.comp} (the weights of the components of the
#' variables) and \code{Weight.sys} (the system of weights).}
#' \item{quality}{The proportion of the total sum of squares explained by the
#' map (the higher the better).}
#' and further details on the parameters of the algorithm.
#' The final map is plotted on the current graphic device.
#' @references Irpino A., Verde R., De Carvalho F.A.T. (2012). Batch self
#' organizing maps for interval and histogram data. In: Proceedings of COMPSTAT
#' 2012, ISI/IASC, 143-154.
#' @details Batch Self-Organising Maps (BSOM) for histogram data, with a cost
#' function based on the L2 Wasserstein distance, i.e. the Euclidean distance
#' between quantile functions. The adaptive version introduces an automatic
#' system of weights in the cost function, to take into account the different
#' effect of the variables on the map.
#'
#' @seealso \code{\link{WH_2d_Adaptive_Kohonen_maps}}
#' @importFrom class somgrid
#' @examples
#' \donttest{
#' results <- WH_2d_Kohonen_maps(
#'   x = BLOOD,
#'   net = list(xdim = 2, ydim = 3, topo = c("rectangular")),
#'   repetitions = 2, simplify = TRUE,
#'   qua = 10, standardize = TRUE
#' )
#' }
#' @export
WH_2d_Kohonen_maps <- function(x, net = list(xdim = 4, ydim = 3, topo = c("rectangular")),
                               kern.param = 2, TMAX = 2, Tmin = 0.2,
                               niter = 30, repetitions = 5,
                               simplify = FALSE,
                               qua = 10,
                               standardize = FALSE, verbose = FALSE) {
  SOL <- WH_2d_Adaptive_Kohonen_maps(x, net, kern.param, TMAX, Tmin,
    niter, repetitions,
    simplify,
    qua,
    standardize,
    schema = 1,
    init.weights = "EQUAL",
    weight.sys = "PROD", theta = 2,
    Wfix = TRUE, verbose = verbose
  )

  return(SOL)
}

print_map <- function(MAP, IDX, gr = FALSE) {
  pp <- table(IDX)
  print(pp)
  indici <- as.numeric(names(pp))
  mappa <- MAP$pts
  countv <- rep(0, nrow(mappa))
  for (h in 1:length(indici)) {
    countv[indici[h]] <- pp[h]
  }
  mappa <- cbind(mappa, count = countv)
  mappa[, 1] <- as.factor(mappa[, 1])
  mappa[, 2] <- as.factor(mappa[, 2])
  MAPPA <- matrix(0, nrow = max(mappa[, 1]), ncol = max(mappa[, 2]))
  for (i in 1:nrow(mappa)) {
    MAPPA[mappa[i, 1], mappa[i, 2]] <- mappa[i, 3]
  }
  if (gr) {
    plot(MAP)
    text(MAP$pts[, 1], MAP$pts[, 2], as.character(c(1:nrow(MAP$pts))), pos = 2, cex = 0.7)
    text(MAP$pts[which(mappa[, 3] > 0), 1], MAP$pts[which(mappa[, 3] > 0), 2], as.character(mappa[which(mappa[, 3] > 0), 3]), pos = 1, cex = 0.9, col = "BLUE")
  }
  return(MAPPA)
}
