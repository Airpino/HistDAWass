//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int whichLowerEQThanX(NumericVector x, double y) {
  // Rcpp supports STL-style iterators
  int elem=0;
  int i=0;
  int stop=0;
  while ((i<=x.size())&&(stop==0)){
    if (x[i]>y){
      stop=1; 
      elem=i;
    }
    i++;
  }
  
  return elem;
}

// [[Rcpp::export]]
NumericVector COMP_Q_VECT(NumericVector x,NumericVector p, NumericVector vp) {
  
  std::sort(vp.begin(), vp.end());
  int vals;
  
  double p_val;
  NumericVector qua(vp.size());
  
  for (vals=0;vals<vp.size();vals++){
    p_val=vp[vals];
    double value=0;
    if (p_val<=0) value=x[0];
    if (p_val>=1) value=x[x.size()-1];  
    if ((p_val>0) && (p_val<1)){
      int pos1;
      double ini,fin;
      pos1=whichLowerEQThanX(p,p_val);
      ini=x[pos1-1];
      fin=x[pos1];
      value=ini+(fin-ini)*(p_val-p[pos1-1])/(p[pos1]-p[pos1-1]);
      
    }
    qua[vals]=value;
  }
  return qua;
}


// [[Rcpp::export]]
NumericVector concatenate_and_sort(NumericVector a,NumericVector b){
  // int DIG=14;
  std::vector<double> ttmp1=as<std::vector<double> >(a);
  std::vector<double> ttmp2=as<std::vector<double> >(b);
  
  ttmp1.insert(ttmp1.end(),ttmp2.begin(),ttmp2.end());
  NumericVector ttmp3;
  ttmp3=wrap(ttmp1);
  // NumericVector ttmp4;
  ttmp3=sort_unique(ttmp3);
  return ttmp3;
}
NumericVector concatenate_and_sort_not_unique(NumericVector a,NumericVector b){
  // int DIG=14;
  std::vector<double> ttmp1=as<std::vector<double> >(a);
  std::vector<double> ttmp2=as<std::vector<double> >(b);
  
  ttmp1.insert(ttmp1.end(),ttmp2.begin(),ttmp2.end());
  NumericVector ttmp3;
  ttmp3=wrap(ttmp1);
  // NumericVector ttmp4;
  std::sort(ttmp3.begin(),ttmp3.end());
  return ttmp3;
}

// [[Rcpp::export]]
Rcpp::List REGISTER2(S4 a, S4 b) {
  S4 x;
  S4 y;
  x=clone(a);
  y=clone(b);
  int DIG=14;
  int i;
  List My(2);
  NumericVector tmp, tmp_p;
  ///////////////////////////////////////////////////////////
  // check if there are empty bins in x
  tmp_p=x.slot("p");
  tmp_p=diff(tmp_p);
  
  if (min(tmp_p)<powf(0.1,DIG)){
    Rcout<<"\n  zero weighted bin!! \n";
    tmp_p[tmp_p<powf(0.1,DIG)]=powf(0.1,DIG);
    tmp_p=tmp_p/(sum(as<NumericVector>(tmp_p)));
    NumericVector tt=cumsum(tmp_p);
    NumericVector tmp_p2(tt.size()+1);
    tmp_p2[0]=0;
    for (i=1;i<(tmp_p2.size()-1);i++)
      tmp_p2[i]=tt[i-1];
    tmp_p2[tmp_p2.size()-1]=1;
    x.slot("p")=tmp_p2;
  }
  //////////////////////////////////////////////////////////
  // check if there are empty bins in y
  tmp_p=y.slot("p");
  tmp_p=diff(tmp_p);
  
  if (min(tmp_p)<powf(0.1,DIG)){
    Rcout<<"\n  zero weighted bin!! \n";
    tmp_p[tmp_p<powf(0.1,DIG)]=powf(0.1,DIG);
    tmp_p=tmp_p/(sum(as<NumericVector>(tmp_p)));
    NumericVector tt=cumsum(tmp_p);
    NumericVector tmp_p2(tt.size()+1);
    tmp_p2[0]=0;
    for (i=1;i<(tmp_p2.size()-1);i++)
      tmp_p2[i]=tt[i-1];
    tmp_p2[tmp_p2.size()-1]=1;
    y.slot("p")=tmp_p2;
  }
  
  ///////////////////////////////////////////////////////
  //create the commoncdf
  NumericVector ttmp3;
  
  ttmp3=concatenate_and_sort(x.slot("p"),y.slot("p"));
  tmp_p=y.slot("p");
  NumericVector ciccio;
  NumericVector t1,t2;
  t1=x.slot("p");
  t2=y.slot("p");
  
  ciccio=setdiff(ttmp3,t1);
  
  x.slot("p")=ttmp3;
  
  if (ciccio.size()>0){
    NumericVector tres,fintrex;
    tres=COMP_Q_VECT(x.slot("x"),t1,ciccio);
    fintrex=concatenate_and_sort_not_unique(x.slot("x"),tres);
    x.slot("x")=fintrex;
  }
  
  y.slot("p")=ttmp3;
  ciccio=setdiff(ttmp3,t2);
  
  if (ciccio.size()>0){
    NumericVector tres,fintrex;
    tres=COMP_Q_VECT(y.slot("x"),t2,ciccio);
    fintrex=concatenate_and_sort_not_unique(y.slot("x"),tres);
    y.slot("x")=fintrex;
  }
  My[0]=x;
  My[1]=y;
  return My;
}


// [[Rcpp::export]]
List PREPARE_A_VEC_MAT(S4 MAT){
  //int DIG=14;
  List My(2);
  NumericVector a;
  List b=MAT.slot("M");
  int i,j;
  std::vector<double> tmp1;
  
  
  for (i=0;i<b.size();i++){
    S4 o=b[i];
    NumericVector t=o.slot("p");
    std::vector<double> tmp2=as<std::vector<double> >(t);
    tmp1.insert(tmp1.end(),tmp2.begin(),tmp2.end());
  }
  NumericVector tmp3;
  tmp3=wrap(tmp1);
  tmp3=sort_unique(tmp3);
  
  NumericMatrix M(b.size(),tmp3.size());
  //rows are individuals, columns are quantiles
  for (i=0;i<b.size();i++){
    S4 o=b[i];
    NumericVector x=o.slot("x");
    NumericVector t2=o.slot("p");
    NumericVector xf(clone(x));
    NumericVector t2f(clone(t2));  
    //std::vector<double> tmp2=as<std::vector<double> >(t);
    //tmp1.insert(tmp1.end(),tmp2.begin(),tmp2.end());
    NumericVector ciccio=setdiff(tmp3,t2f);
    
    if (ciccio.size()>0){
      NumericVector tres,fintrex;
      tres=COMP_Q_VECT(xf,t2f,ciccio);
      fintrex=concatenate_and_sort_not_unique(xf,tres);
      xf=fintrex;
    }
    //Rcout<<"\n STOASSEGNANDO \n";
    for (j=0;j<tmp3.size();j++){
      M(i,j)=xf(j);
    }
    
    
  }
  
  //Rcout<<"\n"<<tmp3<<"\n";
  My(0)=M;
  My(1)=tmp3;
  return My;
}
// [[Rcpp::export]]
S4 MEDIA_V(S4 MAT, NumericVector wei){ //compute the average distribution
  List TMP;
  TMP=PREPARE_A_VEC_MAT(MAT);
  NumericVector p;
  //assuming wei>0 and of size n
  //MAT has n rows and a column  
  p=TMP(1);
  //rows are individuals, columns are quantiles
  NumericMatrix A=TMP(0);
  int ind=A.nrow();
  int qua=A.ncol();
  NumericVector x(qua);
  double SumWei=0;
  int i,j;
  for (i=0;i<ind;i++){
    SumWei=SumWei+wei(i);
    for (j=0;j<qua;j++){
      //   Rcout<<"\n A(i,j) "<<A(i,j)<< " wei i " <<wei(i)<<" x(j) "<<x(j)<<"\n";
      x(j)=x(j)+A(i,j)*wei(i);
      
    }
  }
  
  NumericMatrix resu(qua,2);
  for (j=0;j<qua;j++){
    resu(j,0)=x(j)/SumWei;
    resu(j,1)=p(j);
  }
  
  S4 o("distributionH");
  NumericVector xx=resu(_,0);
  NumericVector pp=resu(_,1);
  
  o.slot("x")=xx;
  o.slot("p")=pp;
  NumericVector rr=diff(xx);
  NumericVector cc(rr.size());
  NumericVector pp2=diff(pp);
  double m=0;
  double s=0;
  
  for (i=0;i<(xx.size()-1);i++){
    cc(i)=(xx[i]+xx[i+1])/2;
    m=m+cc[i]*pp2[i];
    s=s+pp2[i]*(cc[i]*cc[i]+rr[i]*rr[i]/3);
  }
  s=sqrt(s-m*m);
  
  o.slot("m")=m;
  o.slot("s")=s;
  return o;
}

// [[Rcpp::export]]
NumericMatrix SSQ_RCPP(S4 MAT,NumericVector wei){
  List resu;
  ListMatrix MM=MAT.slot("M");
  float Sumwei=sum(wei);
  arma::mat MIXP;
  MIXP.zeros(MM.ncol(),MM.ncol());
  int i;
  for (i=0;i<MM.nrow();i++){
    
    S4 TMPL("MatH");
    TMPL.slot("M")=MM(i,_);
    resu=PREPARE_A_VEC_MAT(TMPL);
    arma::mat M=as<arma::mat>(resu[0]);
    arma::vec p=resu[1];
    arma::mat Mc;
    arma::mat Mr;
    
    
    Mc=(M(span(0,(M.n_rows-1)),span(1,(M.n_cols-1)))+M(span(0,(M.n_rows-1)),span(0,(M.n_cols-2))))/2;
    Mr=(M(span(0,(M.n_rows-1)),span(1,(M.n_cols-1)))-M(span(0,(M.n_rows-1)),span(0,(M.n_cols-2))))/2;
    arma::mat W=diagmat(diff(p));
    
    
    MIXP=MIXP+Mc*W*Mc.t()*wei(i)+Mr*W*Mr.t()*wei(i)/3;
  }
  //now compute the second part
  
  S4 MEANS("MatH");
  ListMatrix TML(1,MM.ncol());
  int j;
  for (j=0;j<MM.ncol();j++){
    S4 TMPL("MatH");
    TMPL.slot("M")=MM(_,j);
    S4 tmp("distributionH");
    tmp=MEDIA_V(TMPL, wei);
    TML(0,j)=tmp;
  }
  MEANS.slot("M")=TML;
  List resu2;
  resu2=PREPARE_A_VEC_MAT(MEANS);
  
  ////////////////////////////////////////////////////////////////////////
  arma::mat M=as<arma::mat>(resu2[0]);
  arma::vec p=resu2[1];
  arma::mat Mc;
  arma::mat Mr;
  arma::mat W=diagmat(diff(p));
  Mc=(M(span(0,(M.n_rows-1)),span(1,(M.n_cols-1)))+M(span(0,(M.n_rows-1)),span(0,(M.n_cols-2))))/2;
  Mr=(M(span(0,(M.n_rows-1)),span(1,(M.n_cols-1)))-M(span(0,(M.n_rows-1)),span(0,(M.n_cols-2))))/2;
  arma::mat MIXP2;
  MIXP=MIXP-Mc*W*Mc.t()*Sumwei-Mr*W*Mr.t()*Sumwei/3;
  
  ////////////////////////////////////////////////////////////////////////
  NumericMatrix fin_resu;
  return fin_resu=wrap(MIXP);
}
// [[Rcpp::export]]
NumericMatrix COV_RCPP(S4 MAT,NumericVector wei){
  NumericMatrix cov;
  cov=SSQ_RCPP(MAT,wei);
  cov=cov/sum(wei);
  return cov;
}
// [[Rcpp::export]]
NumericMatrix CORR_RCPP(S4 MAT,NumericVector wei){
  arma::mat corr=as<arma::mat>(COV_RCPP(MAT,wei));
  arma::vec dia=sqrt(corr.diag());
  corr=corr/(dia*dia.t());
  NumericMatrix fin_resu;
  return fin_resu=wrap(corr);
}


// [[Rcpp::export]]
NumericVector M_STD_H(S4 o){ //compute mean and std of a distrbution 
  S4 o2;
  o2=clone(o);
  NumericVector xx=o2.slot("x");
  NumericVector pp=o2.slot("p");
  
  NumericVector rr=diff(xx);
  rr=rr/2;
  NumericVector cc(rr.size());
  NumericVector pp2=diff(pp);
  double m=0;
  int i;
  double s=0;
  
  for (i=0;i<(xx.size()-1);i++){
    cc(i)=(xx[i]+xx[i+1])/2;
    m=m+cc[i]*pp2[i];
    s=s+pp2[i]*(cc[i]*cc[i]+rr[i]*rr[i]/3);
  }
  
  s=sqrt(s-m*m);
  NumericVector resu(2);
  resu[0]=m;
  resu[1]=s;
  return resu;    
}