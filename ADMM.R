############################################################################
##Use ADMM to optimize parameter with linear constrain on stiefel manifold##
############################################################################
ADMM<-function(X,Y,A,u,lambda,K,ada_w){
  ##X and Y are input data
  ##A is the constrian A%*%beta=0
  ##K is the number of constrain
  ##ada_w is the adaptive lasso weight
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(A,nrow=K,ncol=p)
  ZERO=matrix(0,nrow=K,ncol=p)
  C1=rbind(A,diag(p))
  C2=rbind(ZERO,-diag(p))
  C3=rbind(ZERO,diag(p))
  beta=rep(0,p)
  delta=rep(0,p)
  gamma=rep(0,p)
  a=rep(0,p+K)
  eps=10
  iter=1
  while(iter<20 & eps>10^(-9)){
    ##Optimize beta and project it on manifold
    beta=solve(2*t(X)%*%X/n+t(C1)%*%C1)%*%(2*t(X)%*%Y/n+t(C1)%*%a-t(C1)%*%C2%*%gamma-t(C1)%*%C3%*%delta)
    beta=beta/sqrt(sum(beta*beta))
    ##Solve gamma with regular lasso process
    gamma=Opt_gamma(C1,C2,C3,beta,delta,a,lambda,p,ada_w)
    delta=solve(u*diag(p)+t(C3)%*%C3)%*%(t(C3)%*%a-t(C3)%*%(C1%*%beta+C2%*%gamma))
    a=a-(C1%*%beta+C2%*%gamma+C3%*%delta)
    iter=iter+1
    eps=sum(abs(A%*%beta))
  }
  Loss=sum((Y-X%*%gamma)*(Y-X%*%gamma))
  return(list(gamma,Loss,beta,eps,iter))
}
##############################################
##Optimize gamma given regular lasso process##
##############################################
Opt_gamma<-function(C1,C2,C3,beta,delta,a,lambda,p,ada_w){
  ##ada_w is the adaptive lasso weight
  library(glmnet)
  library(expm)
  ##The design matrix D and first order term
  D=t(C2)%*%C2/2
  Res=t(a)%*%C2-t(C1%*%beta+C3%*%delta)%*%C2
  X=sqrtm(D)
  Y=solve(X)%*%c(Res/2)
  result=glmnet(X,Y,family="gaussian",lambda=lambda,intercept=FALSE,penalty.factor=ada_w)
  return(result$beta[1:p])
}