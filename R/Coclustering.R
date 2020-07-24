# Coclustering of histogram data based on Wasserstein  Distance ----
#' Coclustering of a dataset of histogram-valued data
#' @description The function implements coclustering for a set of histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param kR An integer, the number of row clusters.
#' @param kC An integer, the number of column clusters.
#' @param rep An integer, maximum number of repetitions of the algorithm (default \code{rep}=5).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @param verbose A logic value (default is FALSE). If TRUE, details on computations are shown.   
#' @return a list with the results of the k-means of the set of Histogram-valued data \code{x} into  \code{k} cluster.
#' @slot solution A list.Returns the best solution among the \code{rep}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$IDXR A vector. The row-clusters at which the objects are assigned.
#' @slot solution$IDXC A vector. The column-clusters at which the columns are assigned.
#' @slot solution$cardinalityR A vector. The cardinality of each final row-cluster.
#' @slot solution$cardinalityC A vector. The cardinality of each final column-cluster.
#' @slot solution$proto A \code{MatH} object with the description of column-clusters.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references De Carvalho, F.A.T., Verde, R., Balzanella A., Irpino A., (2020 to appear). Coclustering of histogram data.
#' submitted to Information Science
#' @examples
#' results=WH_coclustering(x = China_Seas,kR = 3, kC=4, rep = 5,simplify = TRUE,
#' qua = 10,standardize = TRUE,verbose=TRUE)
#' @export

WH_coclustering =function (x,kR,kC, rep=5, 
                     simplify=FALSE,
                     qua=10,
                     standardize=FALSE,verbose=FALSE){
  if ((kR<2) || (kR>=nrow(x@M))){
    str=paste("The number of ROW clusters is not appropriate:\n 2<= kR<",
              nrow(x@M), "(no. of rows in x)\n")
    stop(str)
  }
  if ((kC<2) || (kC>=nrow(x@M))){
    str=paste("The number of clusters is not appropriate:\n 2<= kC<",
              ncol(x@M), "(no. of columns in x)\n")
    stop(str)
  }
  
  init="RPART"
  ind=nrow(x@M)
  vars=ncol(x@M)
  NAMER=get.MatH.rownames(x)
  NAMEV=get.MatH.varnames(x)
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  
  ## compute total sum of 
  TOTSSQ=0
  
  tmpd=distributionH()
  for (v in 1:vars){
    #    tmp=ComputeFastSSQ(MM[[v]])
    tmp=c_ComputeFASTSSQ(MM[[v]])
    TOTSSQ=TOTSSQ+tmp$SSQ
    
  }
  solutions=list()
  criteria=numeric(0)
  for (repet in 1:rep){
    if(verbose){
      cat(paste("---------> rep  ",repet,"\n"))}
    ## initialize clusters and prototypes
    proto=new("MatH",nrows=k,ncols=vars)
    cards=rep(0,k)
    SSQ=matrix(0,k,vars)
    GenCrit=Inf;
    if (init=="RPART"){
      rperm=sample(1:ind)
      tmp=rep(c(1:k),ceiling(ind/k))[1:ind]
      IDX=tmp[rperm]
      ##  compute prototypes
      for (clu in 1:k){
        sel=(1:ind)[IDX==clu]
        cards=length(sel)
        for (j in 1:vars){
          #tmp=ComputeFastSSQ(MM[[j]][,c(sel,(ind+1))])
          tmp=c_ComputeFASTSSQ(MM[[j]][,c(sel,(ind+1))])
          proto@M[clu,j][[1]]@x=tmp$mx
          proto@M[clu,j][[1]]@p=tmp$mp
          tt=M_STD_H(proto@M[clu,j][[1]])
          proto@M[clu,j][[1]]@m=tt[1]
          proto@M[clu,j][[1]]@s=tt[2]
          # 
          # tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          # proto@M[clu,j][[1]]=tmpH
          SSQ[clu,j]=tmp$SSQ
        }
      }
    }
    GenCrit=sum(SSQ)
    OK=1 
    maxit=100
    treshold=GenCrit*1e-10
    itcount=0
    
    while (OK==1){
      itcount=itcount+1
      #assign
      tmp=c_DISTA_M(MM,proto)
      
      MD=tmp$dist
      IDX=tmp$wm
      ##########################
      # for (i in 1:ind){
      #     IDX[i]=which.min(MD[i,])
      # }
      # browser()
      ################################
      #recompute
      OldCrit=GenCrit
      
      #check for empty clusters
      moved=numeric(0)
      for (i in 1:k){
        sel=(1:ind)[IDX==i]
        #show(sel)
        if (length(sel)==0){
          cat("empty cluster\n")
          which.max(apply(MD,1,max))
          tomove=which.max(apply(MD,1,max))
          IDX[tomove]=i
          moved=c(moved,tomove)
          MD[tomove,]=rep(0,k)
          #show(MD)
        }
      }
      for (clu in 1:k){
        sel=(1:ind)[IDX==clu]
        
        cards=length(sel)
        for (j in 1:vars){
          #tmp=ComputeFastSSQ(MM[[j]][,c(sel,(ind+1))])
          tmp=c_ComputeFASTSSQ(MM[[j]][,c(sel,(ind+1))])
          proto@M[clu,j][[1]]@x=tmp$mx
          proto@M[clu,j][[1]]@p=tmp$mp
          tt=M_STD_H(proto@M[clu,j][[1]])
          proto@M[clu,j][[1]]@m=tt[1]
          proto@M[clu,j][[1]]@s=tt[2]
          # tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          # proto@M[i,j][[1]]=tmpH
          SSQ[clu,j]=tmp$SSQ
        }
      }
      GenCrit=sum(SSQ)
      if (verbose){
        cat(paste(itcount,GenCrit, "\n", sep="---->"))}
      #check criterion
      if (abs(GenCrit-OldCrit)<treshold){
        OK=0
      }
    }
    cardinality=table(IDX)
    names(IDX)=NAMER
    dimnames(cardinality)$IDX=paste("Cl",dimnames(cardinality)$IDX,sep=".")
    row.names(proto@M)=paste("Cl",c(1:k),sep=".")
    colnames(proto@M)=NAMEV
    solutions=c(solutions,list(solution=list(IDX=IDX,cardinality=cardinality,
                                             centers=proto,Crit=GenCrit)))
    criteria=c(criteria,GenCrit)
  }
  #plot(solutions[[which.min(criteria)]]$centers,type="DENS")
  return(best.solution=list(solution=solutions[[which.min(criteria)]],
                            quality=1-min(criteria)/TOTSSQ))
}