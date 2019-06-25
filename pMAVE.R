############################################################################
##The function for finding bandwidth h and kernel weight matrix given MAVE##
############################################################################
Kernel_H<-function(x,y,B){
  ##The dimension for input matrix x##
  p=ncol(x); 
  n=nrow(x);
  ##Find the optimal bandwidth h and kernel weight given x, y and projector B##
  m=p;
  V=x%*%B;
  h=1.2*mean(apply(V,2,"sd"))/n^(1/(m+4));
  h2=2*h*h;
  v1n=matrix(1,n,1);
  v1p=matrix(1,p,1);
  x=x-v1n%*%apply(x,2,mean);
  invX=MASS::ginv(t(x)%*%x/n + diag(p)/n^3);
  eig=eigen(invX); eigL = eig$va; eigV = eig$ve;
  ss=eigV%*%diag(sqrt(eigL),length(eigL),length(eigL))%*%solve(eigV);
  x=x%*%ss;
  C=V%*%t(V)
  dC=diag(C)
  dxij=dC%*%t(v1n)+v1n%*%dC-2*C;
  ds=apply(dxij,2,sort);
  h2j=ifelse(h2>ds[2*m,],h2,ds[2*m,]);
  ker=exp(-dxij/(v1n%*%h2j));
  ##ker is the optimal kernel weight w_{ij}##
  return(ker);
}


###############################################################################################################
##The function for update B given a, b, ker, tuning parameter lambda and linear constrain B_constrain in MAVE##
###############################################################################################################
Update_B<-function(x,y,a,b,ker,lambda,B_constrain,lasso_weight_B,K){
  library(glmnet)
  library(ncvreg)
  ##K is the number of constrain
  ##The dimension for input matrix x##
  p=ncol(x); 
  n=nrow(x);
  ##Define Lasso input X and Y##
  X=matrix(0,n*n,p);
  Y=rep(0,n);
  for(i in 1:n){
    for(j in 1:n){
      X[((i-1)*n+j),]=ker[i,j]^(1/2)*b[j]*(x[i,]-x[j,]);
      Y[((i-1)*n+j)]=ker[i,j]^(1/2)*(y[i]-a[j]);
    }
  }
  ##Solve the lasso estimator B given X,Y and tuning parameter lambda##
  result=ADMM(X,Y,t(B_constrain),10^6,lambda,K,lasso_weight_B)
  B=result[[1]];
  ##B is the estimation for projector B##
  ##Find the corresponding BIC value##
  N=n^2;
  l=log(sum((Y-X%*%B)*(Y-X%*%B))/N/2)
  BIC=log(sum((Y-X%*%B)*(Y-X%*%B))/N/2)+log(n)*length(which(abs(B)>0))/n;
  return(list(B,BIC,l));
}


#####################################################
##The function for update a and b given projector B##
#####################################################
Update_ab<-function(x,y,B,ker){
  ##The dimension for input matrix x##
  p=ncol(x); 
  n=nrow(x);
  ##Estimate a and b##
  a=rep(0,n);
  b=rep(0,n);
  ##Estimate a and b element by element##
  for(j in 1:n){
    Attr=matrix(0,2,2);
    Res=rep(0,2);
    for(i in 1:n){
      S_ij=rbind(1,sum((x[i,]-x[j,])*B));
      Attr=Attr+S_ij%*%t(S_ij)*ker[i,j];
      Res=Res+S_ij*y[i]*ker[i,j];
    }
    ##Solve the j-th element in a and b##
    c=solve(Attr)%*%Res;
    a[j]=c[1];
    b[j]=c[2];
  }
  return(list(a,b));
}


###############################################################################################################################################
##The main function for estimate B given x,y,initilizer B0, linear constrain B_constrain and A_constrain and tuning parameter lambda for MAVE##
###############################################################################################################################################
Update_MAVE<-function(x,y,lambda,B0,B_constrain,lasso_weight_B,K){
  ##K is the number of constrain
  ##The dimension for input matrix x##
  p=ncol(x); 
  n=nrow(x);
  ##Set the initial value for projector matrix B##
  ##Here B has dimension p*1##
  BIC=100;
  l=100;
  B=B0;
  eps=10;
  iter=1;
  iter_max=5;
  while(eps>0.001 & iter<=iter_max){
    B0=B;
    ##Find the kernel weight
    ker=Kernel_H(x,y,B);
    ##Update a and b, the intercept and gradient for g(Bx)##
    ab=Update_ab(x,y,B,ker);
    a=ab[[1]];
    b=ab[[2]];
    ##Update projector matrix B##
    result=Update_B(x,y,a,b,ker,lambda,B_constrain,lasso_weight_B,K);
    B=result[[1]];
    ##If the all the element in projector B is too small, we will stop and treat it as 0
    if(sum(abs(B))>10^(-3)){
      BIC=result[[2]];
      l=result[[3]];
      B=c(B);
      iter=iter+1;
      eps=sum((B-B0)*(B-B0));
    }
    if(sum(abs(B))<10^(-3)){
      break;
    }
  }
  return(list(B,BIC,a,l));
}


#######################################################################################################
##The joint SIR and MAVE method to optimize A and B given x,y,lambda and linear constrain B_constrain##
#######################################################################################################
pECCA_sub<-function(x,y,lambda,B_constrain,A_constrain,lasso_weight_B,K){
  library(CCA)
  library(dr)
  q=dim(y)[2]
  ##Use CCA to find the initial value##
  B=CCA_initial(x,y,B_constrain,A_constrain)[[1]];
  A=CCA_initial(x,y,B_constrain,A_constrain)[[2]];
  eps=10;
  iter=1;
  iter_max=5;
  while(eps>0.001 & iter<=iter_max){
    A0=A;
    B0=B;
    ##Given A, optimize B in g(Bx)##
    result_B=Update_MAVE(x,y%*%A,lambda,B0,B_constrain,lasso_weight_B,K);
    B=result_B[[1]];
    ##If all the elements in B is too small, then stop the computation##
    if(sum(abs(B))<10^(-3)){
      break;
    }
    ##Given B, optimize A in Ay##
    lds=seq(0.05,1,0.05);
    tuning=c(0);
    result_A=est.rssir(y,c(x%*%B),d=1,lds,tuning,nslices=5,max.iter=100,eps.conv=1e-3,A_constrain,rep(1,q),K);
    A=result_A$out.L1.bic[[1]];
    if(sum(abs(A))>0) A=A/sqrt(sum(A*A));
    eps=sum((B-B0)*(B-B0))+sum((A-A0)*(A-A0));
    iter=iter+1;
  }
  ##If B is exactly 0, BIC is defined as a large number##
  if(sum(abs(B))>10^(-3)){
    B=B/sqrt(sum(B*B));
    BIC=result_B[[2]];
    l=result_B[[4]]
  }
  if(sum(abs(B))<10^(-3)){
    BIC=100;
    l=100;
  }
  return(list(A,B,BIC,l));
}


#####################################################################################################
##Estimate A and B given tuning parameter sequence and linear constrain B_constrain and A_constrain##
#####################################################################################################
pECCA<-function(x,y,B_constrain,A_constrain,K,lam_seq){
  ##Denote a sequence to record A and B##
  total=length(lam_seq)
  p=ncol(x); 
  n=nrow(x);
  q=ncol(y);
  AA=matrix(0,total,q);
  BB=matrix(0,total,p);
  BIC_sequence=rep(0,total);
  ##Define the adaptive lasso weight
  lasso_weight_B=rep(1,p);
  result=pECCA_sub(x,y,0,B_constrain,A_constrain,lasso_weight_B,K);
  lasso_weight_B=1/abs(result[[2]]+10^(-3));
  ##Try grid search tuning parameter##
  for(t in 1:total){
    lambda=lam_seq[t]
    result=pECCA_sub(x,y,lambda,B_constrain,A_constrain,lasso_weight_B,K);
    AA[t,]=result[[1]];
    BB[t,]=result[[2]];
    BIC_sequence[t]=result[[3]];
  }
  ##Find the tuning parameter that minimize the BIC value##
  index=min(which(BIC_sequence==min(BIC_sequence)));
  ##Return optimal A and B##
  A=AA[index,];
  B=BB[index,];
  return(list(A,B,BIC_sequence,index));
}