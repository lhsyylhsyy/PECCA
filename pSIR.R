#---------------------------------------------------------------------------------------------# 
# R functions for sliced inverse regression with L1 / L2 regularizations                      #
# Author: Lexin Li                                                                            # 
# Reference: Li, L., and Yin, X. (2008). Sliced Inverse Regression with Regularizations.      #
#            Biometrics. 64, 124-131.                                                         #
#---------------------------------------------------------------------------------------------#

library(lasso2)
library(dr)

# estimate the structural dimension using the BIC type criterion
ssir.d<-function(X, y, nslices=5, pen.term=1) 
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  
  # center X and slice y
  Xc<-scale(X)
  if(is.null(dim(y))) {
    sy<-dr.slices(y, nslices)
  } else {
    sy<-dr.slices.double(y, nslices)
  }
  h<-sy$nslices
  
  # kernel matrix in x-scale
  M.sir<-matrix(0, nrow=p, ncol=p)
  for(s in 1:h) {
    Xc.s<-Xc[sy$slice.indicator == s, ]
    if(sy$slice.sizes[s] == 1) Xc.s<-matrix(Xc.s, nrow=1)
    X.hs<-as.vector(apply(Xc.s, 2, mean))
    M.sir<-M.sir + (sy$slice.sizes[s]/n) * X.hs %*% t(X.hs)
  }
  
  # eigen decomposition
  dv<-eigen(M.sir)$values
  dv<-(dv+abs(dv))/2
  
  # penalty term
  if(pen.term == 1) C<-log(n) * h/n 
  if(pen.term == 2) C<-(0.5*log(n) + 0.5*n^(1/3)) * h/n
  
  # BIC-type criterion
  dv2<-eigen(M.sir + diag(p))$values
  dv2<-(dv2+abs(dv2))/2
  crit<-NULL
  for(m in 0:(p-1)) {
    crit<-c(crit, 0.5*n*sum((log(dv2)+1-dv2)[(1+min(sum(dv2>1), m)):p]) - 0.5*C*m*(2*p-m+1))
  }
  d<-order(crit, decreasing=T)[1] - 1
  
  # return
  ans<-list(d=d, crit=crit, dv=dv)
  return(ans)  
}


# L2 regularized SIR for a given lambda
ssir.L2<-function(X, y, d, lambda, nslices=5, max.iter=100, eps.conv=1e-6, B.init=NULL)
{
  # parameters 
  n<-nrow(X)
  p<-ncol(X)
  Sigma.x<-cov(X)
  
  # center X and slice y
  Xc<-scale(X)
  if(is.null(dim(y))) {
    sy<-dr.slices(y, nslices)
  } else {
    sy<-dr.slices.double(y, nslices)
  }
  h<-sy$nslices
  
  # compute components for each slice 
  X.h<-f.h<-NULL
  for(s in 1:h) {
    Xc.s<-Xc[sy$slice.indicator == s, ]
    if(sy$slice.sizes[s] == 1) Xc.s<-matrix(Xc.s, nrow=1)
    X.hs<-as.vector(apply(Xc.s, 2, mean))
    X.h<-cbind(X.h, X.hs)
    f.h<-c(f.h, sy$slice.sizes[s]/n)
  }
  
  # initial estimate of B
  if(is.null(B.init)) {
    B.init<-diag(1, p)[,1:d]; if(d == 1) B.init<-matrix(B.init, ncol=1)
  } 
  
  # estimate B and C iteratively
  iter<-0
  B<-B.init
  G<-99
  G.diff<-1
  while((iter < max.iter) & (G.diff > eps.conv)){
    comp1<-solve(t(B) %*% Sigma.x %*% Sigma.x %*% B) %*% t(B) %*% Sigma.x
    C<-NULL; for(s in 1:h) { C<-cbind(C, comp1 %*% X.h[,s]) }
    comp2<-solve((C %*% diag(f.h) %*% t(C)) %x% (Sigma.x %*% Sigma.x) + lambda * diag(1, nrow=p*d, ncol=p*d))
    comp3<-comp2 %*% ((C %*% diag(f.h)) %x% Sigma.x) 
    vecB<-comp3 %*% as.vector(X.h)
    B.new<-matrix(vecB, nrow=p, byrow=FALSE) 
    #B.new<-apply(B.new, 2, norm)
    
    S.lambda<-((diag(f.h^0.5) %*% t(C)) %x% Sigma.x) %*% comp2 %*% ((C %*% diag(f.h^0.5)) %x% Sigma.x) 
    
    comp4<-(diag(1, nrow=p*h, ncol=p*h) - S.lambda) %*% (diag(f.h^0.5) %x% diag(1, p, p)) %*% as.vector(X.h)
    G.new<-as.vector(t(comp4) %*% comp4)
    G.diff<-G[iter+1] - G.new
    
    B<-B.new
    G<-c(G, G.new)
    
    iter<-iter + 1
  }
  G<-G[-1]
  
  # return 
  rss<-G[length(G)]
  n.e<-p * h
  p.e<-sum(diag(S.lambda))
  aic<-n.e * log(rss / n.e) + 2 * p.e
  bic<-n.e * log(rss / n.e) + log(n.e) * p.e
  ric<-(n.e - p.e) * log(rss / (n.e - p.e)) + p.e * (log(n.e) - 1) + 4 / (n.e - p.e - 2)
  gcv<-rss / (n.e * (1 - p.e/n.e)^2)
  
  ans<-list(B=B, G=G, G.min=rss, p.e=p.e, iter=iter, aic=aic, bic=bic, ric=ric, gcv=gcv)
  return(ans)
}


# L1 regularized SIR for a given bound
ssir.L1<-function(X, y, B, lambda, A_constrain, nslices=5, lasso_weight_A,K)
{
  #K is the number of constrain
  # parameters 
  n<-nrow(X)
  p<-ncol(X)
  Sigma.x<-cov(X)
  
  # center X and slice y
  Xc<-scale(X)
  if(is.null(dim(y))) {
    sy<-dr.slices(y, nslices)
  } else {
    sy<-dr.slices.double(y, nslices)
  }
  h<-sy$nslices
  
  # form data for lasso
  y.star<-x.star<-NULL
  for(s in 1:h) {
    Xc.s<-Xc[sy$slice.indicator == s, ]; if(sy$slice.sizes[s] == 1) Xc.s<-matrix(Xc.s, nrow=1)
    X.hs<-as.vector(apply(Xc.s, 2, mean))
    y.star<-c(y.star, X.hs * sqrt(sy$slice.sizes[s]/n))
    
    comp4<-B %*% solve(t(B) %*% Sigma.x %*% Sigma.x %*% B) %*% t(B) %*% Sigma.x * sqrt(sy$slice.sizes[s]/n) 
    x.star<-rbind(x.star, Sigma.x %*% diag(as.vector(comp4 %*% X.hs)))
  }
  
  # estimate alpha by ADMM
  result=ADMM(x.star%*%solve(diag(c(B))),y.star,A_constrain,10^6,lambda,K,lasso_weight_A);
  rss<-result[[2]];
  
  # shrinkage estimate
  B.s<-result[[1]];
  
  # information criteria
  n.e<-nrow(x.star)
  p.e<-sum(abs(result[[1]]) >10^(-6))
  aic<-n.e * log(rss / n.e) + 2 * p.e
  bic<-n.e * log(rss / n.e) + log(n.e) * p.e
  ric<-(n.e - p.e) * log(rss / (n.e - p.e)) + p.e * (log(n.e) - 1) + 4 / (n.e - p.e - 2)
  gcv<-rss / (n.e * (1 - p.e/n.e)^2)
  
  # return 
  ans<-list(B.s=B.s, B=B, rss=rss, p.e=p.e, aic=aic, bic=bic, ric=ric, gcv=gcv, result=result)
  return(ans)
}


# L1+L2 regularized SIR for a series values of lambdas (lds) and a serious values of bounds (bds)
# the best L2 solution is chosen based on GCV; the best L1 solutions for AIC/BIC/RIC/GCV are all given
est.rssir<-function(X, y, d, lds, tuning, nslices=5, max.iter=100, eps.conv=1e-3, A_constrain,lasso_weight_A,K)
{
  # parameters
  n<-nrow(X)
  p<-ncol(X)
  
  # ridge estimation
  rss.L2<-p.e.L2<-aic.L2<-bic.L2<-ric.L2<-gcv.L2<-NULL
  out.L2s<-list()
  for(j in 1:length(lds)) {
    lambda<-lds[j]
    
    out<-ssir.L2(X, y, d=d, lambda=lambda, nslices=nslices, max.iter=max.iter, eps.conv=eps.conv)
    
    out.L2s[[j]]<-out
    rss.L2<-c(rss.L2, out$G.min)
    p.e.L2<-c(p.e.L2, out$p.e)
    aic.L2<-c(aic.L2, out$aic)
    bic.L2<-c(bic.L2, out$bic)
    ric.L2<-c(ric.L2, out$ric)
    gcv.L2<-c(gcv.L2, out$gcv)
  }
  out.L2<-out.L2s[[order(gcv.L2)[1]]]
  crit.L2<-rbind(rss.L2, p.e.L2, aic.L2, bic.L2, ric.L2, gcv.L2)
  rownames(crit.L2)<-c("rss", "p.e", "aic", "bic", "ric", "gcv")
  colnames(crit.L2)<-as.character(lds)
  main=diag(p)
  diag(main)=out.L2$B
  
  
  # lasso estimation
  rss.L1<-p.e.L1<-aic.L1<-bic.L1<-ric.L1<-gcv.L1<-NULL
  for(j in 1:length(tuning)) {
    lambda<-tuning[j]
    
    out<-ssir.L1(X, y, out.L2$B, lambda, t(A_constrain), nslices=5, lasso_weight_A,K)
    if(class(out) == "try-error") out<-list(rss=NA, p.e=NA, aic=NA, bic=NA, ric=NA, gcv=NA)
    
    rss.L1<-c(rss.L1, out$rss)
    p.e.L1<-c(p.e.L1, out$p.e)
    aic.L1<-c(aic.L1, out$aic)
    bic.L1<-c(bic.L1, out$bic)
    ric.L1<-c(ric.L1, out$ric)
    gcv.L1<-c(gcv.L1, out$gcv)
  }
  out.L1.aic<-ssir.L1(X, y, out.L2$B, lambda=tuning[order(aic.L1)[1]], t(A_constrain), nslices=nslices, lasso_weight_A,K)
  out.L1.bic<-ssir.L1(X, y, out.L2$B, lambda=tuning[order(bic.L1)[1]], t(A_constrain), nslices=nslices, lasso_weight_A,K)
  out.L1.ric<-ssir.L1(X, y, out.L2$B, lambda=tuning[order(ric.L1)[1]], t(A_constrain), nslices=nslices, lasso_weight_A,K)
  out.L1.gcv<-ssir.L1(X, y, out.L2$B, lambda=tuning[order(gcv.L1)[1]], t(A_constrain), nslices=nslices, lasso_weight_A,K)
  crit.L1<-rbind(rss.L1, p.e.L1, aic.L1, bic.L1, ric.L1, gcv.L1)
  rownames(crit.L1)<-c("rss", "p.e", "aic", "bic", "ric", "gcv")
  colnames(crit.L1)<-as.character(tuning)
  
  # return
  ans<-list(out.L1.aic=out.L1.aic, out.L1.bic=out.L1.bic, out.L1.ric=out.L1.ric, out.L1.gcv=out.L1.gcv, 
            out.L2s=out.L2s, crit.L2=crit.L2, crit.L1=crit.L1)
  return(ans)
}