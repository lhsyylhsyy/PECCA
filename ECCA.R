
#############################
##ECCA without penalty term##
#############################
ECCA<-function(x,y,num){
  library(MAVE)
  library(dr)
  library(CCA)
  ##num is the total pairs
  A_seq=NULL
  B_seq=NULL
  my_cor=rep(0,num)
  result=AB(x,y)
  # switch x and y
  result.switch = AB(y, x)
  
  # keep whichever result gives the larger correlation
  if(result[[3]] >= result.switch[3]){
    A_seq=cbind(A_seq,result[[1]])
  B_seq=cbind(B_seq,result[[2]])
  my_cor[1]=result[[3]]
  } else{
    A_seq=cbind(A_seq,result.switch[[2]])
    B_seq=cbind(B_seq,result.switch[[1]])
    my_cor[1]=result.switch[[3]]
  }
  
  # subsequent directions
  if(num>=2){
    for(i in 2:num){
      result=AB_next(x,y,A_seq,B_seq)
      # switch x and y
      result.switch = AB_next(y, x, B_seq, A_seq)
      # keep the result with larger correlation
      if(result[[3]] >= result.switch[3]){
        A_seq=cbind(A_seq,result[[1]])
        B_seq=cbind(B_seq,result[[2]])
        my_cor[i]=result[[3]]
      } else{
        A_seq=cbind(A_seq,result.switch[[2]])
        B_seq=cbind(B_seq,result.switch[[1]])
        my_cor[i]=result.switch[[3]]
      }
    }
  }
  #08/19/2017. I added this line of warning message - Ning
  ##if(iter>50) warning('algorithm failed to converge')
  return(list(A_seq,B_seq,my_cor))
}
##########################################################
##Given input x and y, find optimal A and B (first pair)##
##########################################################
AB<-function(x,y){
  eps=10
  iter=1
  ##CCA provides initializer##
  result=cc(x,y);
  A=result$ycoef[,1];
  B=result$xcoef[,1];
  ##Update A and B by SIR and MAVE alternately##
  while(eps>10^(-6) & iter<=50){
    A0=A
    B0=B
    result_B=mave(y%*%A~x,method='meanopg')
    B=result_B$dir[[1]]
    result_A=dr(x%*%B~y,numdir=1,method="sir",nslices=5)
    A=result_A$evectors[,1]
    eps=sum(abs(A-A0))+sum(abs(B-B0))
    iter=iter+1
  }
  ##Provide correlation estimates##
  my_cor=result_A$evalues[1]
  ##Make the solution unique by adjusting the first informative dimension##
  i=which(abs(A)/sum(abs(A))>0.1)[1]
  j=which(abs(B)/sum(abs(B))>0.1)[1]
  if(A[1]<0) A=-A
  if(B[1]<0) B=-B
  return(list(A,B,my_cor))
}
#########################################################
##Given input x and y, find optimal A and B (i-th pair)##
#########################################################
AB_next<-function(x,y,A_seq,B_seq){
  eps=10
  iter=1
  p=dim(x)[2]
  q=dim(y)[2]
  ##CCA provides initialier##
  A_seq=as.matrix(A_seq)
  B_seq=as.matrix(B_seq)
  time=ncol(B_seq);
  P=diag(p)-B_seq%*%solve(t(B_seq)%*%B_seq)%*%t(B_seq);
  Q=diag(q)-A_seq%*%solve(t(A_seq)%*%A_seq)%*%t(A_seq);
  U=svd(P)$u;
  V=svd(Q)$u;
  tilde_U=U[,1:(p-time)];
  tilde_V=V[,1:(q-time)];
  x_sub=x%*%tilde_U;
  y_sub=y%*%tilde_V;
  result=cc(x_sub,y_sub);
  A=result$ycoef[,1];
  B=result$xcoef[,1];
  ##Update A and B by SIR and MAVE alternately##
  while(eps>10^(-6) & iter<=50){
    A0=A
    B0=B
    result_B=mave(y_sub%*%A~x_sub,method='meanopg')
    B=result_B$dir[[1]]
    result_A=dr(x_sub%*%B~y_sub,numdir=1,method="sir",nslices=5)
    A=result_A$evectors[,1]
    eps=sum(abs(A-A0))+sum(abs(B-B0))
    iter=iter+1
  }
  B=tilde_U%*%B;
  A=tilde_V%*%A;
  ##Provide correlation estimates##
  my_cor=result_A$evalues[1]
  ##Make the solution unique by adjusting the first informative dimension##
  i=which(abs(A)/sum(abs(A))>0.1)[1]
  j=which(abs(B)/sum(abs(B))>0.1)[1]
  if(A[1]<0) A=-A
  if(B[1]<0) B=-B
  return(list(A,B,my_cor))
}