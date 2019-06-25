#############################
##Sparse ECCA main function##
#############################
####################################
##Estimate paris of sparse A and B##
####################################
Sparse_ECCA<-function(x,y,num,lam_seq_1,lam_seq_2,tun_seq_1,tun_seq_2){
  ##num is the total pairs
  ##lam_seq is the tuning parameters for B and tun_seq is the tuning parameters for A
  p=ncol(x);
  q=ncol(y);
  B_constrain=rep(0,p);
  A_constrain=rep(0,q);
  ##Estimate first pair of A and B
  ##Optimization for first pair of B
  result_B=pECCA(x,y,B_constrain,A_constrain,0,lam_seq_1);
  B=result_B[[2]];
  lasso_weight_A=1/abs(result_B[[1]]+10^(-3));
  ##Optimization for first pair of A
  lds=10^seq(-6,0,0.1);
  ##tuning=tun_seq;
  result_A=est.rssir(y,c(x%*%B),d=1,lds,tun_seq_1,nslices=5,max.iter=100,eps.conv=1e-3,A_constrain,lasso_weight_A,0);
  A=result_A$out.L1.bic[[1]];
  if(sum(abs(A))>0) A=A/sqrt(sum(A*A));
  A_constrain=A;
  B_constrain=B;
  ##Estimate next num-1 pairs
  if(num>1){
    for(time in 1:(num-1)){
      ##Optimization for time-th pair of B
      result_B=pECCA(x,y,B_constrain,A_constrain,time,lam_seq_2);
      B=result_B[[2]];
      lasso_weight_A=1/abs(result_B[[1]]+10^(-3));
      ##Optimization for time-th pair of A
      lds=10^seq(-6,0,0.1);
      ##tuning=tun_seq
      result_A=est.rssir(y,c(x%*%B),d=1,lds,tun_seq_2,nslices=5,max.iter=100,eps.conv=1e-3,A_constrain,lasso_weight_A,time);
      A=result_A$out.L1.bic[[1]];
      if(sum(abs(A))>0) A=A/sqrt(sum(A*A));
      A_constrain=cbind(A_constrain,A);
      B_constrain=cbind(B_constrain,B);
    }
  }
  ##Output A_constrain and B_constrain are the coefficient estimation
  return(list(A_constrain,B_constrain));
}