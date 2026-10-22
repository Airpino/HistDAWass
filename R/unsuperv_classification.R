# K-means of histogram data based on Wasserstein ----
#' K-means clustering of histogram data
#'
#' K-means (dynamic clustering) of the individuals of a \code{MatH} object, based
#' on the L2 Wasserstein distance between histograms.
#'
#' @param x A \code{MatH} object.
#' @param k An integer, the number of clusters.
#' @param rep An integer, the number of repetitions of the algorithm (with
#' random starting partitions); the best solution is returned. Default 5.
#' @param simplify A logical value (default \code{FALSE}). If \code{TRUE}, each
#' histogram is approximated by \code{qua} equally spaced quantiles: the
#' computations are much faster for large datasets, with a small approximation
#' error.
#' @param qua An integer, the number of quantiles used when \code{simplify = TRUE}
#' (default 10).
#' @param standardize A logical value (default \code{FALSE}). If \code{TRUE}, the
#' variables are standardised by their Wasserstein-based standard deviation.
#' @param verbose A logical value (default \code{FALSE}). If \code{TRUE}, details
#' of the computations are printed.
#' @return A list with:
#' \item{solution}{The best solution among the repetitions (minimum criterion): a
#' list with \code{IDX} (the cluster of each individual), \code{cardinality}
#' (the size of each cluster), \code{centers} (a \code{MatH} object with the
#' prototypes) and \code{Crit} (the within-cluster sum of squares).}
#' \item{quality}{The proportion of the total sum of squares explained by the
#' partition (the higher the better).}
#' @references Irpino A., Verde R., Lechevallier Y. (2006). Dynamic clustering of
#' histograms using Wasserstein metric. In: Rizzi A., Vichi M. (eds.) COMPSTAT
#' 2006 - Advances in computational statistics, Physica-Verlag, 869-876.
#' @seealso \code{\link{WH_adaptive.kmeans}}, \code{\link{WH_fcmeans}},
#' \code{\link{WH_hclust}}
#' @examples
#' results <- WH_kmeans(
#'   x = BLOOD, k = 2, rep = 10, simplify = TRUE,
#'   qua = 10, standardize = TRUE
#' )
#' results$solution$IDX
#' results$quality
#' @export
WH_kmeans <- function(x, k, rep = 5,
                      simplify = FALSE,
                      qua = 10,
                      standardize = FALSE, verbose = FALSE) {
  if ((k < 2) || (k >= nrow(x@M))) {
    str <- paste(
      "The number of clusters is not appropriate:\n 2<= k<",
      nrow(x@M), "(no. of individuals in x)\n"
    )
    stop(str)
  }

  init <- "RPART"
  ind <- nrow(x@M)
  vars <- ncol(x@M)
  NAMER <- get.MatH.rownames(x)
  NAMEV <- get.MatH.varnames(x)
  ## we homogeneize data for speeding up the code if required
  tmp <- Prepare(x, simplify, qua, standardize)
  MM <- tmp$MM
  x <- tmp$x

  ## compute total sum of
  TOTSSQ <- 0

  tmpd <- distributionH()
  for (v in 1:vars) {
    #    tmp=ComputeFastSSQ(MM[[v]])
    tmp <- c_ComputeFASTSSQ(MM[[v]])
    TOTSSQ <- TOTSSQ + tmp$SSQ
  }
  solutions <- list()
  criteria <- numeric(0)
  for (repet in 1:rep) {
    if (verbose) {
      cat(paste("---------> rep  ", repet, "\n"))
    }
    ## initialize clusters and prototypes
    proto <- new("MatH", nrows = k, ncols = vars)
    cards <- rep(0, k)
    SSQ <- matrix(0, k, vars)
    GenCrit <- Inf
    if (init == "RPART") {
      rperm <- sample(1:ind)
      tmp <- rep(c(1:k), ceiling(ind / k))[1:ind]
      IDX <- tmp[rperm]
      ##  compute prototypes
      for (clu in 1:k) {
        sel <- (1:ind)[IDX == clu]
        cards <- length(sel)
        for (j in 1:vars) {
          # tmp=ComputeFastSSQ(MM[[j]][,c(sel,(ind+1))])
          tmp <- c_ComputeFASTSSQ(MM[[j]][, c(sel, (ind + 1))])
          proto@M[clu, j][[1]]@x <- tmp$mx
          proto@M[clu, j][[1]]@p <- tmp$mp
          tt <- M_STD_H(proto@M[clu, j][[1]])
          proto@M[clu, j][[1]]@m <- tt[1]
          proto@M[clu, j][[1]]@s <- tt[2]
          #
          # tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          # proto@M[clu,j][[1]]=tmpH
          SSQ[clu, j] <- tmp$SSQ
        }
      }
    }
    GenCrit <- sum(SSQ)
    OK <- 1
    maxit <- 100
    treshold <- GenCrit * 1e-10
    itcount <- 0

    while (OK == 1) {
      itcount <- itcount + 1
      # assign
      tmp <- c_DISTA_M(MM, proto)

      MD <- tmp$dist
      IDX <- tmp$wm
      ##########################
      # for (i in 1:ind){
      #     IDX[i]=which.min(MD[i,])
      # }
      # browser()
      ################################
      # recompute
      OldCrit <- GenCrit

      # check for empty clusters
      moved <- numeric(0)
      for (i in 1:k) {
        sel <- (1:ind)[IDX == i]
        # show(sel)
        if (length(sel) == 0) {
          cat("empty cluster\n")
          which.max(apply(MD, 1, max))
          tomove <- which.max(apply(MD, 1, max))
          IDX[tomove] <- i
          moved <- c(moved, tomove)
          MD[tomove, ] <- rep(0, k)
          # show(MD)
        }
      }
      for (clu in 1:k) {
        sel <- (1:ind)[IDX == clu]

        cards <- length(sel)
        for (j in 1:vars) {
          # tmp=ComputeFastSSQ(MM[[j]][,c(sel,(ind+1))])
          tmp <- c_ComputeFASTSSQ(MM[[j]][, c(sel, (ind + 1))])
          proto@M[clu, j][[1]]@x <- tmp$mx
          proto@M[clu, j][[1]]@p <- tmp$mp
          tt <- M_STD_H(proto@M[clu, j][[1]])
          proto@M[clu, j][[1]]@m <- tt[1]
          proto@M[clu, j][[1]]@s <- tt[2]
          # tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          # proto@M[i,j][[1]]=tmpH
          SSQ[clu, j] <- tmp$SSQ
        }
      }
      GenCrit <- sum(SSQ)
      if (verbose) {
        cat(paste(itcount, GenCrit, "\n", sep = "---->"))
      }
      # check criterion
      if (abs(GenCrit - OldCrit) < treshold) {
        OK <- 0
      }
    }
    cardinality <- table(IDX)
    names(IDX) <- NAMER
    dimnames(cardinality)$IDX <- paste("Cl", dimnames(cardinality)$IDX, sep = ".")
    row.names(proto@M) <- paste("Cl", c(1:k), sep = ".")
    colnames(proto@M) <- NAMEV
    solutions <- c(solutions, list(solution = list(
      IDX = IDX, cardinality = cardinality,
      centers = proto, Crit = GenCrit
    )))
    criteria <- c(criteria, GenCrit)
  }
  # plot(solutions[[which.min(criteria)]]$centers,type="DENS")
  return(best.solution = list(
    solution = solutions[[which.min(criteria)]],
    quality = 1 - min(criteria) / TOTSSQ
  ))
}


# Hierarchical clustering of histogram data based on Wasserstein ----
#' Hierarchical clustering of histogram data
#'
#' Hierarchical clustering of the individuals of a \code{MatH} object, based on
#' the L2 Wasserstein distance. It is a wrapper of \code{\link[stats]{hclust}}.
#'
#' @inheritParams WH_kmeans
#' @param distance The distance: only \code{"WDIST"}, the L2 Wasserstein
#' distance, is available.
#' @param method The agglomeration method, as in \code{\link[stats]{hclust}}:
#' one of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}
#' (default), \code{"average"}, \code{"mcquitty"}, \code{"median"} or
#' \code{"centroid"}.
#' @return An object of class \code{hclust}.
#' @references Irpino A., Verde R. (2006). A new Wasserstein based distance for
#' the hierarchical clustering of histogram symbolic data. In: Batanjeli et al.
#' (eds.) Data Science and Classification, IFCS 2006, Springer, 185-192.
#' @seealso \code{\link{WH_MAT_DIST}}
#' @examples
#' results <- WH_hclust(x = BLOOD, simplify = TRUE, method = "complete")
#' plot(results)
#' cutree(results, k = 5)
#' @importFrom stats hclust quantile as.dist
#' @export
WH_hclust <- function(x,
                      simplify = FALSE,
                      qua = 10,
                      standardize = FALSE,
                      distance = "WDIST",
                      method = "complete") {
  ind <- nrow(x@M)
  vars <- ncol(x@M)

  ## we homogeneize data for speeding up the code if required
  tmp <- Prepare(x, simplify, qua, standardize)
  MM <- tmp$MM
  x <- tmp$x

  ## compute distance matrix
  if (distance == "WDIST") {
    d <- sqrt(c_Fast_D_Mat(MM))
  }

  rownames(d) <- rownames(x@M)
  colnames(d) <- rownames(x@M)
  hc <- hclust(as.dist(d), method = method)
  return(hc)
}


# Adaptive-distance based k-means ----
#' K-means clustering of histogram data with adaptive distances
#'
#' K-means of the individuals of a \code{MatH} object based on adaptive L2
#' Wasserstein distances: the algorithm also estimates a weight for each variable
#' (or for each component of the variables, and/or for each cluster).
#'
#' @inheritParams WH_kmeans
#' @param schema An integer, the system of weights: 1 = a weight for each
#' variable (default); 2 = a weight for the position and one for the variability
#' component of each variable; 3 = as 1, with a different set of weights for
#' each cluster; 4 = as 2, with a different set of weights for each cluster.
#' @param init Not used.
#' @param weight.sys The constraint on the weights: \code{"PROD"} (default), the
#' product of the weights is one; \code{"SUM"}, the weights sum to one.
#' @param theta A number, the exponent of the weights when
#' \code{weight.sys = "SUM"} (default 2).
#' @param init.weights The initial weights: \code{"EQUAL"} (default) or
#' \code{"RANDOM"}.
#' @return A list with \code{IDX} (the cluster of each individual),
#' \code{cardinality} (the size of each cluster), \code{proto} (a \code{MatH}
#' object with the prototypes), \code{weights} (the weights of the variables or
#' components), \code{Crit} (the criterion), \code{TOTSSQ}, \code{BSQ},
#' \code{WSQ} (total, between and within sums of squares) and \code{quality} (the
#' proportion of the total sum of squares explained by the partition).
#' @references Irpino A., Verde R., De Carvalho F.A.T. (2014). Dynamic clustering
#' of histogram data based on adaptive squared Wasserstein distances. Expert
#' Systems with Applications, 41, 3351-3366. \doi{10.1016/j.eswa.2013.12.001}
#' @seealso \code{\link{WH_kmeans}}, \code{\link{WH_adaptive_fcmeans}}
#' @examples
#' results <- WH_adaptive.kmeans(x = BLOOD, k = 2, rep = 10,
#'             simplify = TRUE, qua = 10, standardize = TRUE)
#' results$weights
#' @importFrom stats runif
#' @export

WH_adaptive.kmeans <- function(x, k,
                               schema = 1, # 1=VariableGLOBAL 2=componentGLOBAL 3=Variable x Cluster 4=Components x cluster
                               init, rep,
                               simplify = FALSE,
                               qua = 10,
                               # empty.action="singleton",
                               standardize = FALSE,
                               weight.sys = "PROD",
                               theta = 2,
                               init.weights = "EQUAL",
                               verbose = FALSE) {
  if ((k < 2) || (k >= nrow(x@M))) {
    str <- paste(
      "The number of clusters is not appropriate:\n 2<= k<",
      nrow(x@M), "(no. of individuals in x)\n"
    )
    stop(str)
  }
  if (missing(init)) {
    init <- "RPART"
  }
  if (missing(rep)) {
    rep <- 5
  }
  ind <- nrow(x@M)
  vars <- ncol(x@M)
  if (weight.sys == "PROD") {
    theta <- 1
  }
  ## we homogeneize data for speeding up the code if required
  tmp <- Prepare(x, simplify, qua, standardize)
  MM <- tmp$MM
  x <- tmp$x

  ## compute total sum of
  TOTSSQ <- 0
  for (v in 1:vars) {
    # tmp=ComputeFastSSQ(MM[[v]])
    tmp <- c_ComputeFASTSSQ(MM[[v]])
    TOTSSQ <- TOTSSQ + tmp$SSQ
  }

  solutions <- list()
  criteria <- numeric(0)
  repet <- 0
  while (repet < rep) {
    repet <- repet + 1
    if (verbose) cat(paste("---------> rep  ", repet, "\n"))
    # initialize matrix of variables' weights
    if (init.weights == "EQUAL") {
      if (verbose) cat("Weights initialization === EQUAL  \n")
      if (weight.sys == "PROD") {
        lambdas <- matrix(1, 2 * vars, k)
      }
      else {
        lambdas <- matrix(1 / vars, 2 * vars, k)
      }
    }
    else { # Random initialization
      if (verbose) cat("Weights initialization === RANDOM  \n")
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
    ## initialize clusters and prototypes
    proto <- new("MatH", nrows = k, ncols = vars)
    cards <- rep(0, k)
    SSQ <- matrix(0, k, vars)
    GenCrit <- Inf
    if (init == "RPART") {
      rperm <- sample(1:ind)
      tmp <- rep(c(1:k), ceiling(ind / k))[1:ind]
      IDX <- tmp[rperm]
    }
    memb <- matrix(0, ind, k)
    for (indiv in 1:ind) {
      memb[indiv, IDX[indiv]] <- 1
    }
    GenCrit <- Inf
    OK <- 1
    maxit <- 100
    treshold <- 1e-5
    itcount <- 0

    while (OK == 1) {
      itcount <- itcount + 1
      # STEP 1 ##  compute prototypes
      for (cluster in 1:k) {
        sel <- (1:ind)[IDX == cluster]

        cards <- length(sel)
        for (variables in 1:vars) {
          tmp <- c_ComputeFASTSSQ(MM[[variables]][, c(sel, (ind + 1))])
          proto@M[cluster, variables][[1]]@x <- tmp$mx
          proto@M[cluster, variables][[1]]@p <- tmp$mp
          tt <- M_STD_H(proto@M[cluster, variables][[1]])
          proto@M[cluster, variables][[1]]@m <- tt[1]
          proto@M[cluster, variables][[1]]@s <- tt[2]
          SSQ[cluster, variables] <- tmp$SSQ
        }
      }
      #       TMP_SSQ=WH.ADPT.FCMEANS.SSQ(x,memb,1,lambdas,proto)
      #       cat('\n after prototypes SSQ:', TMP_SSQ,'\n')
      # STEP 2 ##  compute weights (Lambda) Fixed Prototypes and partition
      #######################################################################################

      TMPRESU <- c_STEP_2_ADA_KMEANS(MM, proto, schema, memb)
      distances <- array(0, dim = c(vars, k, 2))
      distances[, , 1] <- TMPRESU$distances1
      distances[, , 2] <- TMPRESU$distances2
      diINDtoPROT <- array(0, dim = c(ind, vars, k, 2))
      for (cluster in 1:k) {
        diINDtoPROT[, , cluster, 1] <- TMPRESU$dIpro_m[[cluster]]
        diINDtoPROT[, , cluster, 2] <- TMPRESU$dIpro_v[[cluster]]
      }



      ##  compute weights
      # S2.2) Weights computation
      if (weight.sys == "PROD") {
        PROSUM <- 1
      } else {
        PROSUM <- 2
      }
      lambdas <- c_STEP_2_2_WEIGHTS_ADA_KMEANS(
        TMPRESU$distances1,
        TMPRESU$distances2,
        PROSUM, schema, theta
      )

      ############################################################################
      # recompute
      OldCrit <- GenCrit
      # S3) affectation (prototype, weights are fixed)
      resu <- c_STEP_3_AFFECT_ADA_KMEANS(
        lambdas,
        TMPRESU$dIpro_m, TMPRESU$dIpro_v, ind, k, vars
      )
      IDX <- resu$IDX

      # check for empty clusters and restart
      empty <- FALSE
      for (i in 1:k) {
        if (sum(IDX == i) == 0) {
          cat("empty cluster\n")
          empty <- TRUE
          OK <- 0
        }
      }
      if (empty) {
        repet <- repet - 1
      }
      memb <- matrix(0, ind, k)
      for (indiv in 1:ind) {
        memb[indiv, IDX[indiv]] <- 1
      }
      TMP_SSQ <- resu$SSQ
      GenCrit <- TMP_SSQ
      if (is.na(GenCrit)) {
        cat("isNAN")
      }

      if (verbose) cat(paste(itcount, GenCrit, "\n", sep = "---->"))
      # check criterion
      if (abs(GenCrit - OldCrit) < treshold) {
        OK <- 0
      }
    }
    if (!empty) {
      cardinality <- table(IDX)
      dimnames(cardinality)$IDX <- paste("Cl", dimnames(cardinality)$IDX, sep = ".")
      solutions <- c(solutions, list(solution = list(IDX = IDX, 
                                                     cardinality = cardinality, 
                                                     proto = proto, weights = lambdas, 
                                                     Crit = GenCrit)))
      criteria <- c(criteria, GenCrit)
    }
  }
  best.solution <- c(solutions[[which.min(criteria)]])
  # plot(best.solution$proto,type="DENS")
  memb <- matrix(0, ind, k)
  for (i in 1:ind) {
    memb[i, best.solution$IDX[i]] <- 1
  }
  colnames(best.solution$weights) <- paste0("Cl.", c(1:k))
  row.names(best.solution$weights) <- paste0(rep(get.MatH.varnames(x), each = 2), c("(Cen.)", "(Var.)"))
  resTOTSQ <- c_WH_ADPT_KMEANS_TOTALSSQ(x, memb,
    m = 1,
    lambdas = best.solution$weights, proto = best.solution$proto
  )
  GEN_proto <- resTOTSQ$gen
  DET_TSQ <- list(TSQ_M = resTOTSQ$TSQ_m, TSQ_V = resTOTSQ$TSQ_v)
  TOTFSSQ <- sum(resTOTSQ$TSQ_m + resTOTSQ$TSQ_v)
  DET_TSQ2 <- list(TSQ2_M = resTOTSQ$TSQ2_m, TSQ2_V = resTOTSQ$TSQ2_v)
  TOTFSSQ2 <- sum(resTOTSQ$TSQ2_m + resTOTSQ$TSQ2_v)
  BSQ <- sum(resTOTSQ$BSQ_m + resTOTSQ$BSQ_v)

  quality1 <- 1 - min(criteria) / TOTFSSQ
  quality2 <- 1 - min(criteria) / TOTFSSQ2
  best.solution <- c(best.solution, TOTSSQ = TOTFSSQ2, BSQ = BSQ, WSQ = best.solution$Crit, quality = BSQ / TOTFSSQ2)
  return(best.solution)
}

# L2  Wasserstein distance matrix----
#' L2 Wasserstein distance matrix
#'
#' Computes the matrix of the L2 Wasserstein distances between the individuals
#' (rows) of a \code{MatH} object, summing the squared distances over the
#' variables.
#'
#' @inheritParams WH_kmeans
#' @return A symmetric matrix of L2 Wasserstein distances.
#' @references Irpino A., Verde R. (2006). A new Wasserstein based distance for
#' the hierarchical clustering of histogram symbolic data. In: Batanjeli et al.
#' (eds.) Data Science and Classification, IFCS 2006, Springer, 185-192.
#' @seealso \code{\link{WassSqDistH}}, \code{\link{WH_hclust}}
#' @examples
#' DMAT <- WH_MAT_DIST(x = BLOOD, simplify = TRUE)
#' @importFrom stats quantile
#' @export
WH_MAT_DIST<- function(x,
                      simplify = FALSE,
                      qua = 10,
                      standardize = FALSE) {
  ind <- nrow(x@M)
  vars <- ncol(x@M)
  
  ## we homogeneize data for speeding up the code if required
  tmp <- Prepare(x, simplify, qua, standardize)
  MM <- tmp$MM
  x <- tmp$x
  
  ## compute distance matrix
  d <- sqrt(c_Fast_D_Mat(MM))
  
  rownames(d) <- rownames(x@M)
  colnames(d) <- rownames(x@M)
  return(d)
}
