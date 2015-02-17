
# res=Prepare(x,simplify,qua,standardize)
Prepare =function (x, simplify, qua, standardize){
  ind=nrow(x@M)
  vars=ncol(x@M)
  MM=list()
  if (simplify){
    p=(0:qua)/qua
    for (j in 1:vars){
      MM[[j]]=matrix(0,(qua+1),(ind+1))
      MM[[j]][,(ind+1)]=p
      for (i in 1:ind){
        dom=numeric(0)
        for (q in 0:qua){
          dom=c(dom, compQ(x@M[i,j][[1]],q/qua))
        }
        
        tmp=new("distributionH",dom,p) 
        tmp=(tmp-tmp@m)*(x@M[i,j][[1]]@s/tmp@s)+x@M[i,j][[1]]@m #transformation with invariance with respect mean and std
        x@M[i,j][[1]]=tmp
        MM[[j]][,i]=tmp@x
      }
    } 
  }
  else{
    
    for (j in 1:vars){
      tmp=registerMH(x[,j])
      MM[[j]]=matrix(0,length(tmp@M[1,1][[1]]@x),ind+1)
      MM[[j]][,ind+1]=tmp@M[1,1][[1]]@p
      for (i in 1:ind){
        x@M[i,j][[1]]=tmp@M[i,1][[1]]
        MM[[j]][,i]=tmp@M[i,1][[1]]@x
      }
    }
  }
  ## standardize data if required
  if (standardize){
    cat("Standardizing data...\n")
    STAND=rep(0,vars)
    Mc=rep(0,vars)
    # compute varianaces
    for (v in 1:vars){
      STAND[v]=sqrt(WH.var.covar(x[,v]))
      Mc[v]=(WH.vec.mean(x[,v]))@m
      for (i in 1:ind){
        if (STAND[v]>0){
          x@M[i,v][[1]]=new("distributionH",x=(x@M[i,v][[1]]@x-Mc[v])/STAND[v],p=x@M[i,v][[1]]@p)
          MM[[v]][,i]=(MM[[v]][,i]-Mc[v])/STAND[v]
          #x@M[i,v][[1]]@x=(x@M[i,v][[1]]@x-Mc[v])/STAND[v]
        }
      }
    }
    
  }
  return(list(MM=MM,x=x))
}

# compute fast SSQ
ComputeFastSSQ=function(subMM){
  ind=ncol(subMM)-1
  rr=nrow(subMM)
  SSQ=0
  m1=apply(as.matrix(subMM[,1:ind]),1,mean)
  p1=subMM[2:rr,ind+1]-subMM[1:(rr-1),ind+1]
  cm=(m1[2:rr]+m1[1:(rr-1)])/2
  rm=(m1[2:rr]-m1[1:(rr-1)])/2
  for (indiv in 1:ind){
    ci=(subMM[2:rr,indiv]+subMM[1:(rr-1),indiv])/2
    ri=(subMM[2:rr,indiv]-subMM[1:(rr-1),indiv])/2
    SSQ=SSQ+sum(p1*((ci-cm)^2)+ p1/3*((ri-rm)^2))
  }
  return(list(SSQ=SSQ,mx=m1, mp=subMM[,ind+1]))
}

# compute fast DIST
ComputeFast_L2_SQ_WASS_D=function(subMM){
  ind=ncol(subMM)-1
  rr=nrow(subMM)
  p1=subMM[2:rr,ind+1]-subMM[1:(rr-1),ind+1]
  
  c1=(subMM[2:rr,1]+subMM[1:(rr-1),1])/2
  r1=(subMM[2:rr,1]-subMM[1:(rr-1),1])/2
  c2=(subMM[2:rr,2]+subMM[1:(rr-1),2])/2
  r2=(subMM[2:rr,2]-subMM[1:(rr-1),2])/2
  
  Dist=sum(p1*((c1-c2)^2)+ p1/3*((r1-r2)^2))
  
  return(Dist)
}

WH.ADPT.KMEANS.TOTALSSQ=function(x,memb,m,lambdas,proto){
  vars=ncol(x@M)
  ind=nrow(x@M)
  k=ncol(memb)
  lambdas[is.na(lambdas)]=0
  #Compute general weights for the global PROTOTYPE
  mu=apply(memb^m,2,sum)
  protoGEN=x[1,]
  for (variables in (1:vars)){
    protoGEN@M[1,variables][[1]]=WH.vec.mean(proto[,variables],mu)
  }
  #Compute general weights for the global PROTOTYPE2
  Mat_of_means=matrix(0,ind,vars)
  CentredD=x
  for (i in 1:ind){
    for (j in 1:vars){
      Mat_of_means[i,j]=x@M[i,j][[1]]@m
      CentredD@M[i,j][[1]]=x@M[i,j][[1]]-x@M[i,j][[1]]@m
    }
  }
  W_centers=matrix(0,ind,vars)
  W_dist_cent=matrix(0,ind,vars)
  for (i in 1:ind){
    for (j in 1:vars){
      W_centers[i,j]=lambdas[(j*2-1),which.max(memb[i,])]
      W_dist_cent[i,j]=lambdas[(j*2),which.max(memb[i,])]
    }
  }
  protoGEN2=x[1,]
  mprot=matrix(0,1,vars)
  for (variables in (1:vars)){
    mprot[1,variables]=sum(Mat_of_means[,variables]*W_centers[,variables])/sum(W_centers[,variables])
    protoGEN2@M[1,variables][[1]]=WH.vec.mean(CentredD[,variables],W_dist_cent[,variables])+mprot[1,variables]
  }
  #Compute BETWEEN
  BSQ_clu=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (variables in (1:vars)){
    GENPROTm=protoGEN2@M[1,variables][[1]]@m
    GENPROTcent=protoGEN2@M[1,variables][[1]]-GENPROTm
    for (clu in 1:k){
      LOCPROTm=proto@M[clu,variables][[1]]@m
      LOCPROTcent=proto@M[clu,variables][[1]]-LOCPROTm
      BSQ_clu[1,variables,clu]=sum(memb[,clu])*lambdas[(variables*2-1),clu]*(GENPROTm-LOCPROTm)^2
      BSQ_clu[2,variables,clu]=sum(memb[,clu])*lambdas[(variables*2),clu]*WassSqDistH(GENPROTcent,LOCPROTcent)
    }
  }
  #Compute the total fuzzy sum of SQUARES
  TSQ_clu=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        tmpD=WassSqDistH(x@M[indiv,variables][[1]],protoGEN@M[1,variables][[1]],details=T)
        tmpD_mean=as.numeric(tmpD[2])
        tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
        TSQ_clu[1,variables,cluster]=TSQ_clu[1,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2-1),cluster]*tmpD_mean
        TSQ_clu[2,variables,cluster]=TSQ_clu[2,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2),cluster]*tmpD_centered
        #         TSQ[1,variables]=TSQ[1,variables]+((memb[indiv,cluster])^m)*lambdas[(variables*2-1),cluster]*tmpD_mean
        #         TSQ[2,variables]=TSQ[2,variables]+((memb[indiv,cluster])^m)*lambdas[(variables*2),cluster]*tmpD_centered
        #         TSQ_clu[1,variables,cluster]=TSQ_clu[1,variables,cluster]+((memb[indiv,cluster])^m)*tmpD_mean
        #         TSQ_clu[2,variables,cluster]=TSQ_clu[2,variables,cluster]+((memb[indiv,cluster])^m)*tmpD_centered
        
      }
    }
  }
  TSQ_clu2=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        tmpD=WassSqDistH(x@M[indiv,variables][[1]],protoGEN2@M[1,variables][[1]],details=T)
        tmpD_mean=as.numeric(tmpD[2])
        tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
        TSQ_clu2[1,variables,cluster]=TSQ_clu2[1,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2-1),cluster]*tmpD_mean
        TSQ_clu2[2,variables,cluster]=TSQ_clu2[2,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2),cluster]*tmpD_centered
        #         TSQ[1,variables]=TSQ[1,variables]+((memb[indiv,cluster])^m)*lambdas[(variables*2-1),cluster]*tmpD_mean
        #         TSQ[2,variables]=TSQ[2,variables]+((memb[indiv,cluster])^m)*lambdas[(variables*2),cluster]*tmpD_centered
        
      }
    }
  }
  return(RES=list(protoGEN=proto,TSQdetailed=TSQ_clu, TSQ=sum(TSQ_clu),TSQdetailed2=TSQ_clu2, TSQ2=sum(TSQ_clu2),
                  BSQdetailed=BSQ_clu, BSQ=sum(BSQ_clu)))
}

ComputeFastSSQ_Fuzzy=function(subMM,memb,m){
  ind=ncol(subMM)-1
  rr=nrow(subMM)
  SSQ=0
  W=memb^m
  W2=W/sum(W)
  m1=subMM[,1:ind]%*%as.matrix(W2)
  
  m1=as.vector(m1)
  p1=subMM[2:rr,ind+1]-subMM[1:(rr-1),ind+1]
  cm=(m1[2:rr]+m1[1:(rr-1)])/2
  rm=(m1[2:rr]-m1[1:(rr-1)])/2
  for (indiv in 1:ind){
    ci=(subMM[2:rr,indiv]+subMM[1:(rr-1),indiv])/2
    ri=(subMM[2:rr,indiv]-subMM[1:(rr-1),indiv])/2
    SSQ=SSQ+(sum(p1*((ci-cm)^2)+ p1/3*((ri-rm)^2)))*W[indiv]
  }
  return(list(SSQ=SSQ,mx=m1, mp=subMM[,ind+1]))
}

ComputeFast_Fuzzy_TOT_SSQ=function(MM,proto,memb,m){
  ind=ncol(MM[[1]])-1
  var=length(MM)
  k=ncol(memb)
#computeGenProt
mus=apply(memb^m,2,sum)
mus=mus/sum(mus)
GP=new('MatH',1,var)
for (variable in 1:var){
  dG=MM[[variable]][,1]*0
  for (clu in 1:k){
    dG=dG+proto@M[clu,variable][[1]]@x*mus[clu]
  }
  tmp=new('distributionH',x=dG,p=MM[[variable]][,(ind+1)])
  GP@M[1,variable][[1]]=tmp
}
#compute TOTALSSQ
SSQ_det=array(0,dim = c(2,var,clu), dimnames=list(c('Mean','Variability'),colnames(proto@M),rownames(proto@M)))
for (variable in 1:var){
  rr=nrow(MM[[variable]])
  p1=MM[[variable]][2:rr,ind+1]-MM[[variable]][1:(rr-1),ind+1]
  cm=(GP@M[1,variable][[1]]@x[2:rr]+GP@M[1,variable][[1]]@x[1:(rr-1)])/2
  rm=(GP@M[1,variable][[1]]@p[2:rr]-GP@M[1,variable][[1]]@p[1:(rr-1)])/2
  for (indiv in 1:ind){
    ci=(MM[[variable]][2:rr,indiv]+MM[[variable]][1:(rr-1),indiv])/2
    ri=(MM[[variable]][2:rr,indiv]-MM[[variable]][1:(rr-1),indiv])/2
    dist=sum(p1*((ci-cm)^2+1/3*(ri-rm)^2))
    dc=(sum(p1*ci)-GP@M[1,variable][[1]]@m)^2
      dv=dist-dc
    for (clu in 1:k){
      SSQ_det[1,variable,clu]=SSQ_det[1,variable,clu]+memb[indiv,clu]^m*dc
      SSQ_det[2,variable,clu]=SSQ_det[2,variable,clu]+memb[indiv,clu]^m*dv
    }
  }
}
return(list(SSQ=sum(SSQ_det),SSQ_det=SSQ_det,ProtoGEN=GP))
}

ComputeFast_Fuzzy_Adaptive_TOT_SSQ=function(MM,proto,memb,m,lambdas){
  ind=ncol(MM[[1]])-1
  var=length(MM)
  k=ncol(memb)
  #computeGenProt
  
  GP=new('MatH',1,var)
  for (variable in 1:var){
    rr=nrow(MM[[variable]])
    ci=(MM[[variable]][2:rr,1:ind]+MM[[variable]][1:(rr-1),1:ind])/2
    ri=(MM[[variable]][2:rr,1:ind]-MM[[variable]][1:(rr-1),1:ind])/2
    wi=(MM[[variable]][2:rr,(ind+1)]-MM[[variable]][1:(rr-1),(ind+1)])
    LMc=memb^m*matrix(lambdas[(variable*2-1),],ind,k, byrow=TRUE)
    LMc=LMc/sum(LMc)
    LMv=memb^m*matrix(lambdas[(variable*2),],ind,k, byrow=TRUE)
    LMv=LMv/sum(LMv)
    Cen=MM[[variable]][,1]*0
    MG=0
    McG=ci[,1]*0
    RcG=ri[,1]*0
    for (clu in 1:k){
      for (i in 1:ind){
      mui=sum(ci[,i]*wi)
      MG=MG+LMc[i,clu]*mui
      Cen=Cen+(MM[[variable]][,i]-mui)*LMv[i,clu]
      }
    }
    Cen=Cen+MG
    tmp=new('distributionH',x=Cen,p=MM[[variable]][,(ind+1)])
    GP@M[1,variable][[1]]=tmp
  }
  #compute G by protos
#   GP2=new('MatH',1,var)
#   for (v in 1:var){
#     LMc=memb^m*matrix(lambdas[(v*2-1),],ind,k, byrow=TRUE)
#     LMc=LMc/sum(LMc)
#     LWc=apply(LMc,2,sum)
#     LMv=memb^m*matrix(lambdas[(v*2),],ind,k, byrow=TRUE)
#     LMv=LMv/sum(LMv)
#     LWv=apply(LMv,2,sum)
#     CG=proto@M[1,v][[1]]@m*LWc[1]
#     CEN=(proto@M[1,v][[1]]-proto@M[1,v][[1]]@m)*LWv[1]
#     for (clu in 2:k){
#       CG=CG+proto@M[clu,v][[1]]@m*LWc[clu]
#       CEN=CEN+(proto@M[clu,v][[1]]-proto@M[clu,v][[1]]@m)*LWv[clu]
#     }
#     GP2@M[1,v][[1]]=CEN+CG
# 
#   
#   }
  
  
  #compute TOTALSSQ
  SSQ_det=array(0,dim = c(2,var,clu), dimnames=list(c('Mean','Variability'),colnames(proto@M),rownames(proto@M)))
  for (variable in 1:var){
    rr=nrow(MM[[variable]])
    p1=MM[[variable]][2:rr,ind+1]-MM[[variable]][1:(rr-1),ind+1]
    cm=(GP@M[1,variable][[1]]@x[2:rr]+GP@M[1,variable][[1]]@x[1:(rr-1)])/2
    rm=(GP@M[1,variable][[1]]@p[2:rr]-GP@M[1,variable][[1]]@p[1:(rr-1)])/2
    for (indiv in 1:ind){
      ci=(MM[[variable]][2:rr,indiv]+MM[[variable]][1:(rr-1),indiv])/2
      ri=(MM[[variable]][2:rr,indiv]-MM[[variable]][1:(rr-1),indiv])/2
      dist=sum(p1*((ci-cm)^2+1/3*(ri-rm)^2))
      dc=(sum(p1*ci)-GP@M[1,variable][[1]]@m)^2
      dv=dist-dc
      for (clu in 1:k){
        SSQ_det[1,variable,clu]=SSQ_det[1,variable,clu]+lambdas[(variable*2-1),clu]*(memb[indiv,clu]^m)*dc
        SSQ_det[2,variable,clu]=SSQ_det[2,variable,clu]+lambdas[(variable*2),clu]*(memb[indiv,clu]^m)*dv
      }
    }
  }
  return(list(SSQ=sum(SSQ_det),SSQ_det=SSQ_det,ProtoGEN=GP))
}