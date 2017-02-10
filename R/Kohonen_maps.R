# Batch SOM Kohonen Maps of HD using adaptive distances -----
#' Batch Kohonen self-organizing 2d maps using adaptive distances for  histogram-valued data
#' @description The function implements a Batch Kohonen self-organizing 2d maps algorithm for  histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param net a list describing the topology of the net \code{list(xdim=number of rows,
#' ydim=numbers of columns,topo=c('rectangular' or 'hexagonal'))}, see \code{somgrid} sintax in package\pkg{class}
#' default \code{net=list(xdim=4,ydim=3,topo=c('rectangular'))}
#' @param kern.param (default =2) the kernel parameter for the RBF kernel used in the algorithm
#' @param TMAX a parameter useful for the iterations (default=2)
#' @param Tmin a parameter useful for the iterations (default=0.2)
#' @param niter maximum number of iterations (default=30)
#' @param repetitions number of repetion of the algorithm (default=5), beacuase each launch may generate a local optimum
#' @param simplify a logical parameter for speeding up computations (default=FALSE). If true data are recoded in order to have fast computations
#' @param qua if \code{simplify=TRUE} number of equally spaced quantiles for recodify the histograms (default=10)
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one. 
#' @param schema a number from 1 to 4 \cr
#' 1=A weight for each variable (default) \cr
#' 2=A weight for the average and the dispersion component of each variable\cr
#' 3=Same as 1 but a different set of weights for each cluster\cr
#' 4=Same as 2 but a different set of weights for each cluster 
#' @param init.weights a string how to initialize weights: 'EQUAL' (default), all weights are the same, 
#' @param weight.sys a string. Weights may add to one ('SUM') or their product is equal to 1 ('PROD', default).
#' @param theta a number. A parameter if \code{weight.sys='SUM'}, default is 2.  
#' @param Wfix a logical parameter (default=FALSE). If TRUE the algorithm does not use adaptive distances
#'         
#' @return a list with the results of the Batch Kohonen map
#' @slot solution A list.Returns the best solution among the \code{repetitions}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$MAP The map topology.
#' @slot solution$IDX A vector. The clusters at which the objects are assigned.
#' @slot solution$cardinality A vector. The cardinality of each final cluster.
#' @slot solution$proto A \code{MatH} object with the description of centers.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot solution$Weights.comp the final weights assigned to each component of the histogram variables
#' @slot solution$Weight.sys a string the type of weighting system ('SUM' or 'PRODUCT')
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references 
#' Irpino A, Verde R, De Carvalho FAT (2012). Batch self organizing maps for interval and histogram data.
#' In: Proceedings of COMPSTAT 2012. p. 143-154, ISI/IASC, ISBN: 978-90-73592-32-2
#' @details An extension of Batch Self Organised Map (BSOM) is here proposed for  histogram data.
#'  These kind of data have been defined in the context of symbolic data analysis.
#'   The BSOM cost function is then based on a distance 
#'   function: the L2 Wasserstein distance. This distance has been widely proposed in several
#'    techniques of analysis (clustering, regression) when input data are expressed by distributions 
#'    (empirical by histograms or theoretical by probability distributions).
#'    The peculiarity of such distance is to be an Euclidean distance between quantile functions so 
#'    that all the properties proved for L2 distances are verified again. An adaptative versions of 
#'    BSOM is also introduced considering an automatic system of weights in the cost function in 
#'    order to take into account the different effect of the several variables in the Self-Organised Map 
#'    grid. 
#' 
#' @importFrom class somgrid
#' @importFrom stats lm dist runif
#' @examples
#' \dontrun{
#' results=WH_2d_Adaptive_Kohonen_maps(x = BLOOD,k = 2,
#'                                    net=list(xdim=2,ydim=3,topo=c('rectangular')), 
#'                                    repetitions = 2,simplify = TRUE,
#'                                    qua = 10,standardize = TRUE)
#'                                    }
#' @export
WH_2d_Adaptive_Kohonen_maps =function (x,net=list(xdim=4,ydim=3,topo=c('rectangular')), 
                                       kern.param=2, TMAX=-9999, Tmin=-9999, 
                                       niter=30,repetitions ,
                                       simplify=FALSE,
                                       qua=10,
                                       standardize=FALSE, schema=4,
                                       init.weights='EQUAL',weight.sys='PROD',
                                       theta=2,Wfix=FALSE){
  tol=1e-10
  # remember to fix TMAX e Tmin
  #require(class)#for somgrid function
  # Check initial conditions and passed parameters
  
  ind=nrow(x@M)
  vars=ncol(x@M)
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  # end of the preprocessing step
  
  TOTSSQ=0
  for (v in 1:vars){
    TOTSSQ=TOTSSQ+(WH.SSQ(x[,v]))
  }
  ##------------------------------------------------batchKOHONENMAP-----------------------------------------------------
  solutions=list()
  criteria=numeric(0)
  MAP=somgrid(net$xdim,net$ydim,net$topo)
  k=nrow(MAP$pts)
  if (missing(repetitions)){
    repetitions=5
  }
  print(repetitions)
  
  for (repet in (1:repetitions)){
    
    
    data=x@M
    nd=nrow(data)
    #random selection of prototypes of neurons
    init=data[sample(1L:nd,k, replace=FALSE), ,drop=FALSE]
    proto=new("MatH")
    proto@M=init
    nhbrdist=as.matrix(dist(MAP$pts))
    dmax=(max(MAP$pts[,1])-min(MAP$pts[,1]))^2+(max(MAP$pts[,2])-min(MAP$pts[,2]))^2
    dmax=(max(as.matrix(dist(MAP$pts))))^2

    if (TMAX<0){
      TMAX=sqrt(-(dmax/4)/(2*log(0.99)))
    }
    if (Tmin<0){
      Tmin=sqrt(-(1/4)/(2*log(0.01)))
    }
    TT=TMAX
    KT=exp(-(nhbrdist^2)/(2*TT^2))
    #initialize weights
    if (Wfix){lambdas=matrix(1,2*vars,k)}else{
      if (init.weights=='EQUAL'){
        cat("Weights initialization === EQUAL  \n")
        if (weight.sys=='PROD'){
          lambdas=matrix(1,2*vars,k)}
        else{lambdas=matrix(1/vars,2*vars,k)}
      }
      else{#Random initialization
        cat("Weights initialization === RANDOM  \n")
        m1=matrix(runif((vars*k),0.01,0.99),vars,k)
        m2=matrix(runif((vars*k),0.01,0.99),vars,k)
        m1=m1/matrix(rep(apply(m1,2,sum)),vars,k,byrow = TRUE)
        m2=m2/matrix(rep(apply(m2,2,sum)),vars,k,byrow = TRUE)
        if (weight.sys=='PROD'){
          m1=exp(m1*vars-1)
          m2=exp(m2*vars-1)
          
        }
        
        if (schema==1){
          m1=matrix(m1[,1],nrow = vars,ncol = k)
          m2=m1
          
        }
        if (schema==2){
          m1=matrix(rep(m1[,1],k),vars,k)
          m2=matrix(rep(m2[,1],k),vars,k)
        }
        if (schema==3){m2=m1}
        
        lambdas=matrix(0,2*vars,k)
        colnames(lambdas)=paste('Clust',c(1:k),sep="_")
        
        n1=paste('M_Var.',c(1:vars),sep="_")
        n2=paste('C_Var.',c(1:vars),sep="_")
        nr=list()
        for (nn in 1:vars){
          nr[[nn*2-1]]=n1[[nn]]
          nr[[nn*2]]=n2[[nn]]
        }
        rownames(lambdas)=nr
        lambdas[(c(1:vars)*2-1),]=m1
        lambdas[(c(1:vars)*2),]=m2
      }}
    #assign objects to closest neuron
    #compute adaptive distances to prototypes
    distances=array(0,dim=c(vars,k,2))
    diINDtoPROT=array(0,dim=c(ind,vars,k,2))
    # compute distances
    ###############################
    #computing distances
    Dmat=array(0,c(ind,vars,k))
    MD=matrix(0,ind,k)
    for (cl in 1:k){
      tmp=ComputeFast_L2_SQ_WASS_DMAT(MM,proto[cl,])
      Dmat[,,cl]=tmp
      MD[,cl]=MD[,cl]+apply(tmp,1,sum)
    }
   
    Dind2Neu=MD%*%KT
    # for (indiv in 1:ind){
    #   for (cluster in 1:k){#s
    #     if (Dind2Neu[indiv,cluster]==Inf) Dind2Neu[indiv,cluster]=0
    #     for (otherclust in 1:k){#l
    #       Dind2Neu[indiv,cluster]=Dind2Neu[indiv,cluster]+KT[otherclust,cluster]*MD[indiv,otherclust]
    #     }
    #   }
    # }
    #initialize matrix of memberships
    IDX=apply(Dind2Neu,1,which.min)
    memb=matrix(0,ind,k)
    for (i in (1:ind)){
      memb[i,IDX[i]]=1
    }
    cat(paste("---------> repetitions  ",repet,"\n"))
    ## initialize clusters and prototypes
    proto=new("MatH",list(new("distributionH")),nrows=k,ncols=vars)
    SSQ1=matrix(0,k,vars)
    GenCrit=Inf
    ##  compute initial criterion
    
    SSQ=0
    for (cluster in 1:k){
      for (variables in (1:vars)){
        for (indiv in 1:ind){
          if (memb[indiv,cluster]>0)
            SSQ=SSQ+Dind2Neu[indiv,IDX[indiv]]
        }
      }
    }
    GenCrit=SSQ
    OK=1 
    maxit=100
    
    t=0
    while (TT>=Tmin){
      t=t+1
      TT=TMAX*(Tmin/TMAX)^(t/(niter-1))
      KT=exp(-(nhbrdist^2)/(2*TT^2))
      KT[which(KT<1e-100)]=0
      #computing prototypes
      CardinalityOfNeuron=apply(memb,2,sum)
      
      for (cluster in 1:k){
        kern=rep(0,ind)
        for (individuals in 1:ind){
          kern[individuals]=KT[IDX[individuals],cluster]
        }
        for (variables in 1:vars){
          if(max(kern)<1e-200){
            kern[which(kern<1e-200)]=1e-200
           }
          M_tmp=matrix(kern/sum(kern),nrow=nrow(MM[[variables]]),ncol=ind,byrow=TRUE)
          M_tmp=M_tmp*MM[[variables]][,1:ind]
          if(length(which(is.nan(M_tmp)))>0){browser()}
          x_tmp=apply(M_tmp,1,sum)
          p_tmp=MM[[variables]][,(ind+1)]
          proto@M[cluster,variables][[1]]@x=x_tmp
          proto@M[cluster,variables][[1]]@p=p_tmp
          proto@M[cluster,variables][[1]]@m=meanH(proto@M[cluster,variables][[1]])
          proto@M[cluster,variables][[1]]@s=stdH(proto@M[cluster,variables][[1]])
          #### RICORDATI DI CALCOLARE M E S
          #proto@M[cluster,variables][[1]]=distributionH(x_tmp,p_tmp)
        }
      }
      
      #compute weights
      if (Wfix==FALSE){
        #first compute distances using kernel
        #  computing distances
        Dmat=array(0,c(ind,vars,k))
        Dmat_M=array(0,c(ind,vars,k))
        Dmat_V=array(0,c(ind,vars,k))
        MD_TOT=matrix(0,vars,k)
        MD_M=matrix(0,vars,k)
        MD_V=matrix(0,vars,k)
        for (cl in 1:k){
          tmp =ComputeFast_L2_SQ_WASS_DMAT(MM,proto[cl,],DETA=TRUE)
          if(length(which(is.nan(tmp$Dist)))>0){browser()}
          Dmat[,,cl]=tmp$Dist
          Dmat_M[,,cl]=tmp$DM
          Dmat_V[,,cl]=tmp$DV
        }
        GRKT=0*Dmat
        for (i1 in 1:ind){
          for (cl in 1:k){
            GRKT[i1,,cl]=KT[IDX[i1],cl]
          }
        }
        # 
        
        distances=array(0,dim=c(vars,k,2))
        distances[,,1]=apply(Dmat_M*GRKT,c(2,3),sum)
        distances[,,2]=apply(Dmat_V*GRKT,c(2,3),sum)
        diINDtoPROT=array(0,dim=c(ind,vars,k,2))
        for (cluster in 1:k){
          diINDtoPROT[,,cluster,1]=Dmat_M[,,cluster]
          diINDtoPROT[,,cluster,2]=Dmat_V[,,cluster]
        }
        #################################################################### DONE 
        #Weights computation
        
        lambdas=ADA_F_DIST(distances,k,vars,ind,schema,weight.sys,theta,m=1)
      }
      #assign data to neurons
      #assign objects to closest neuron
      #compute adaptive distances to prototypes
      distances=array(0,dim=c(vars,k,2))
      diINDtoPROT=array(0,dim=c(ind,vars,k,2))
      ####################
      for (cluster in 1:k){
        for (variables in (1:vars)){
          for (indiv in 1:ind){
            if (schema==1){#one weight for one variable
              tmpD=Dmat_M[indiv,variables,cluster]
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*lambdas[(variables*2-1),cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD
            }
            if (schema==2|schema==5){#two weights for the two components for each variable
              tmpD=Dmat_M[indiv,variables,cluster]
              tmpD_mean=Dmat_M[indiv,variables,cluster]
              tmpD_centered=Dmat_V[indiv,variables,cluster]
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*lambdas[(variables*2-1),cluster]
              distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*lambdas[(variables*2),cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
              diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
            }
            if (schema==3){#a weight for one variable and for each cluster
              tmpD=Dmat_M[indiv,variables,cluster]
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*lambdas[(variables*2-1),cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD
            }
            if (schema==4|schema==6){#two weights for the two components for each variable and each cluster
              tmpD=Dmat_M[indiv,variables,cluster]
              tmpD_mean=Dmat_M[indiv,variables,cluster]
              tmpD_centered=Dmat_V[indiv,variables,cluster]
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*lambdas[(variables*2-1),cluster]
              distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*lambdas[(variables*2),cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
              diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
            }
          }
        }
      }
      tmp=apply(diINDtoPROT,c(1,3),sum)
     # tmp[which(is.na(tmp))]=1e100
      Dind2Neu=tmp%*%KT#matrix(Inf,ind,k)
      # for (indiv in 1:ind){
      #   for (cluster in 1:k){#s
      #     if (Dind2Neu[indiv,cluster]==Inf) Dind2Neu[indiv,cluster]=0
      #     for (otherclust in 1:k){#l
      #       Dind2Neu[indiv,cluster]=Dind2Neu[indiv,cluster]+KT[otherclust,cluster]*sum(diINDtoPROT[indiv,,otherclust,])
      #     }
      #   }
      # }
      #recompute criterion
      IDX=apply(Dind2Neu,1,which.min)
      OldCrit=GenCrit
      GenCrit=0
      for (individual in 1:ind){
        GenCrit=GenCrit+sum(Dind2Neu[individual, IDX[individual]])
      }
      
      cat(paste(t,GenCrit, "\n", sep="---->"))
      if (is.na(GenCrit)){
        browser()
        print('A GREAT PROBLEM HERE!')
      }
      print(print_map(MAP,IDX))
      #print(TT)
      #print(max(KT))
      #print(min(KT[which(KT>1e-16)]))
    }
    #crisp assignment
    
    cardinality=table(IDX)
    #dimnames(cardinality)$IDX=paste("Cl",1:k,sep=".")
    solutions=c(solutions,list(solution=list(MAP=MAP, IDX=IDX,cardinality=cardinality,proto=proto,
                                             Crit=GenCrit,weights.comp=lambdas,Weight.sys=weight.sys)))
    criteria=c(criteria,GenCrit)
    #plot(proto,type="DENS")
  }
  FM=solutions[[which.min(criteria)]]$MAP
    FIDX=solutions[[which.min(criteria)]]$IDX
  print_map(FM,FIDX,gr=TRUE)
  #print(dmax)
  #print(TMAX)
  #print(Tmin)
  return(best.solution=list(solution=solutions[[which.min(criteria)]],
                            repetitions=which.min(criteria),
                            quality=1-min(criteria)/TOTSSQ,dmax=dmax,TMAX=TMAX,Tmin=Tmin))
}

# Batch SOM of HD -----
#' Batch Kohonen self-organizing 2d maps for  histogram-valued data
#' @description The function implements a Batch Kohonen self-organizing 2d maps algorithm for  histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param net a list describing the topology of the net \code{list(xdim=number of rows,
#' ydim=numbers of columns,topo=c('rectangular' or 'hexagonal'))}, see \code{somgrid} sintax in package\pkg{class} 
#' default \code{net=list(xdim=4,ydim=3,topo=c('rectangular'))}
#' @param kern.param (default =2) the kernel parameter for the RBF kernel used in the algorithm
#' @param TMAX a parameter useful for the iterations (default=2)
#' @param Tmin a parameter useful for the iterations (default=0.2)
#' @param niter maximum number of iterations (default=30)
#' @param repetitions number of repetion of the algorithm (default=5), beacuase each launch may generate a local optimum
#' @param simplify a logical parameter for speeding up computations (default=FALSE). If true data are recoded in order to have fast computations
#' @param qua if \code{simplify=TRUE} number of equally spaced quantiles for recodify the histograms (default=10)
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @return a list with the results of the Batch Kohonen map
#'@slot solution A list.Returns the best solution among the \code{repetitions}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$MAP The map topology.
#' @slot solution$IDX A vector. The clusters at which the objects are assigned.
#' @slot solution$cardinality A vector. The cardinality of each final cluster.
#' @slot solution$proto A \code{MatH} object with the description of centers.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references 
#' Irpino A, Verde R, De Carvalho FAT (2012). Batch self organizing maps for interval and histogram data.
#' In: Proceedings of COMPSTAT 2012. p. 143-154, ISI/IASC, ISBN: 978-90-73592-32-2
#' @details An extension of Batch Self Organised Map (BSOM) is here proposed for  histogram data.
#'  These kind of data have been defined in the context of symbolic data analysis.
#'   The BSOM cost function is then based on a distance 
#'   function: the L2 Wasserstein distance. This distance has been widely proposed in several
#'    techniques of analysis (clustering, regression) when input data are expressed by distributions 
#'    (empirical by histograms or theoretical by probability distributions).
#'    The peculiarity of such distance is to be an Euclidean distance between quantile functions so 
#'    that all the properties proved for L2 distances are verified again. An adaptative versions of 
#'    BSOM is also introduced considering an automatic system of weights in the cost function in 
#'    order to take into account the different effect of the several variables in the Self-Organised Map 
#'    grid. 
#' 
#' @importFrom class somgrid
#' @examples
#' \dontrun{
#' results=WH_2d_Kohonen_maps(x = BLOOD,k = 2,
#'                                    net=list(xdim=2,ydim=3,topo=c('rectangular')), 
#'                                    repetitions = 2,simplify = TRUE,
#'                                    qua = 10,standardize = TRUE)
#'                                    }
#' @export
WH_2d_Kohonen_maps =function (x,net=list(xdim=4,ydim=3,topo=c('rectangular')),
                              kern.param=2, TMAX=2, Tmin=0.2, 
                              niter=30,repetitions=5 ,      
                              simplify=FALSE,
                              qua=10,
                              standardize=FALSE){
  SOL=WH_2d_Adaptive_Kohonen_maps(x,net, kern.param, TMAX, Tmin, 
                                  niter,repetitions ,
                                  simplify,
                                  qua,
                                  standardize, schema=1,
                                  init.weights='EQUAL',weight.sys='PROD',theta=2,Wfix=TRUE)
  
  return(SOL)
}

WH.ADPT.KOHONEN.SSQ=function(x,memb,m,lambdas,proto){
  vars=ncol(x@M)
  ind=nrow(x@M)
  k=ncol(memb)
  lambdas[is.na(lambdas)]=0
  SSQ=0
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
        tmpD_mean=as.numeric(tmpD[2])
        tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
        SSQ=SSQ+((memb[indiv,cluster])^m)*(lambdas[(variables*2-1),cluster]*tmpD_mean+lambdas[(variables*2),cluster]*tmpD_centered)
        
      }
    }
  }
  return(SSQ)
}

print_map=function(MAP,IDX,gr=FALSE){
  pp=table(IDX)
  print(pp)
  indici=as.numeric(names(pp))
  mappa=MAP$pts
  countv=rep(0,nrow(mappa))
  for (h in 1:length(indici)){
    countv[indici[h]]=pp[h]
  }
  mappa=cbind(mappa,count=countv)
  mappa[,1]=as.factor(mappa[,1])
  mappa[,2]=as.factor(mappa[,2])
  MAPPA=matrix(0,nrow=max(mappa[,1]),ncol=max(mappa[,2]))
  for (i in 1:nrow(mappa)){
    MAPPA[mappa[i,1],mappa[i,2]]=mappa[i,3]
  }
  if (gr){
  plot(MAP)
  text(MAP$pts[,1],MAP$pts[,2],as.character(c(1:nrow(MAP$pts))),pos=2,cex=0.7)
  text(MAP$pts[which(mappa[,3]>0),1],MAP$pts[which(mappa[,3]>0),2],as.character(mappa[which(mappa[,3]>0),3]),pos=1,cex=0.9,col='BLUE')
  }
  return(MAPPA)
}