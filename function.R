############ Codes for our algorithm 
#dim:dimensions of Y, W and X.
dim=c(p,q,r)
Hier=function(data,ob,lambda1,lambda2,lambda3,dim,maxiter=50){
  #######Input:
  #data: a list of data whose elements are matrix
  #ob: number of observation
  #lambda1,lambda2,lambda3: tuning parameters
  #dim: a vector whose elements are variable dimensions of 
  #each data matrix in the list 'data'
  #maxiter: maximum number of cycles
  #######Output:
  # A list including: 
  #three precision matrices
  #a vector called 'iter' that records the number of iterations in each regression
  #bic_value:recording BIC value
  res_1=res_2=res_3=matrix(0,p,p)
  Xp=matrix(0,ob,p)
  diag(res_1)=diag(res_2)=diag(res_3)=1
  iter=rep(0,p)
  Bic=vector()
  for (i in 1:p) {
    Xp[,p]=data[[1]][,i]
    Xp[,-p]=data[[1]][,-i]
    Z1=data[[2]]
    Z2=data[[3]]
    Z1Z2=cbind(dataZ1,dataZ2)
    
    #Setting initial value
    B1=rep(0,(dim[1]-1))
    R1=rep(0,dim[2])
    B2=rep(0,(dim[1]-1))
    N1=rep(0,dim[3])
    Delta=rep(0,(dim[1]-1))
    R2N2=rep(0,(dim[2]+dim[3]))
    
    P1=mean((Xp[,p]-Xp[,-p]%*%B1-Z1%*%R1)^2)+mean((Xp[,p]-Xp[,-p]%*%B2-Z2%*%N1)^2)+
      mean((Xp[,p]-Xp[,-p]%*%(B1*B2*Delta)-Z1Z2%*%R2N2)^2)
    P2=lambda1*(sum(abs(B1))+sum(abs(B2)))+lambda2*sum(abs(Delta))+
      lambda3*(sum(abs(R1))+sum(abs(N1))+sum(abs(R2N2)))
    P=P1+P2
    for (j in 1:maxiter) {# j:number of current iterations
      Delta_old=Delta
      Beta1_old=B1
      Beta2_old=B2
      R1_old=R1
      N1_old=N1
      R2N2_old=R2N2
      #update Delta
      for(k in 1:(p-1)){
        if(B1[k]!=0&B2[k]!=0){
          Coef=Xp[,k]*B1[k]*B2[k]
          mean_Coef=mean(Coef^2)
          judge=mean((Xp[,p]-Xp[,-c(k,p)]%*%(B1[-k]*B2[-k]*Delta[-k])-Z1Z2%*%R2N2)*Coef)
          if((abs(judge)-lambda2)>0){
            Delta[k]=sign(judge)*(abs(judge)-lambda2)/mean_Coef
          }
          else{
            Delta[k]=0
          }
        }
        else{
          Delta[k]=0
        }
      }
      #update R2 N2
      for (k in 1:(2*p)) {
        Coef=Z1Z2[,k]
        mean_Coef=mean(Coef^2)
        judge=mean((Xp[,p]-Xp[,-p]%*%(B1*B2*Delta)-Z1Z2[,-k]%*%R2N2[-k])*Coef)
        if((abs(judge)-lambda3)>0){
          R2N2[k]=sign(judge)*(abs(judge)-lambda3)/mean_Coef
        }
        else{
          R2N2[k]=0
        }
      }
      #update Beta1
      for (k in 1:(p-1)) {
        judge=mean((Xp[,p]-Xp[,-c(k,p)]%*%B1[-k]-Z1%*%R1)*Xp[,k]+
                     (Xp[,p]-Xp[,-c(k,p)]%*%(B1[-k]*B2[-k]*Delta[-k])-Z1Z2%*%R2N2)*
                     (Xp[,k]*B2[k]*Delta[k]))
        if((abs(judge)-lambda1)>0){
          B1[k]=sign(judge)*(abs(judge)-lambda1)/mean(Xp[,k]^2+(Xp[,k]*B2[k]*Delta[k])^2)
        }
        else{
          B1[k]=0
        }
      }
      #update R1
      for (k in 1:p) {
        judge=mean((Xp[,p]-Xp[,-p]%*%B1-Z1[,-k]%*%R1[-k])*Z1[,k])
        if((abs(judge)-lambda3)>0){
          R1[k]=sign(judge)*(abs(judge)-lambda3)/mean(Z1[,k]^2)
        }
        else{
          R1[k]=0
        }
      }
      #update Beta2
      for (k in 1:(p-1)) {
        judge=mean((Xp[,p]-Xp[,-c(k,p)]%*%B2[-k]-Z2%*%N1)*Xp[,k]+
                     (Xp[,p]-Xp[,-c(k,p)]%*%(B1[-k]*B2[-k]*Delta[-k])-Z1Z2%*%R2N2)*
                     (Xp[,k]*B1[k]*Delta[k]))
        if((abs(judge)-lambda1)>0){
          B2[k]=sign(judge)*(abs(judge)-lambda1)/mean(Xp[,k]^2+
                                                      (Xp[,k]*B1[k]*Delta[k])^2)
        }
        else{
          B2[k]=0
        }
      }
      #update N1
      for (k in 1:p) {
        judge=mean((Xp[,p]-Xp[,-p]%*%B2-Z2[,-k]%*%N1[-k])*Z2[,k])
        if((abs(judge)-lambda3)>0){
          N1[k]=sign(judge)*(abs(judge)-lambda3)/mean(Z2[,k]^2)
        }
        else{
          N1[k]=0
        }
      }
      #end condition of iteration
      Q1=mean((Xp[,p]-Xp[,-p]%*%B1-Z1%*%R1)^2)+mean((Xp[,p]-Xp[,-p]%*%B2-Z2%*%N1)^2)+
        mean((Xp[,p]-Xp[,-p]%*%(B1*B2*Delta)-Z1Z2%*%R2N2)^2)
      Q2=lambda1*(sum(abs(B1))+sum(abs(B2)))+lambda2*sum(abs(Delta))+
        lambda3*(sum(abs(R1))+sum(abs(N1))+sum(abs(R2N2)))
      Q=Q1+Q2
      if(j==maxiter|abs((P-Q)/P)<0.01){
        res_1[i,-i]=B1
        res_2[i,-i]=B2
        res_3[i,-i]=B1*B2*Delta
        iter[i]=j
        b1=-ob*log(mean(Q1))
        b2=(sum(sign(abs(B1)))+sum(sign(abs(R1)))+sum(sign(abs(B2)))+
              sum(sign(abs(N1)))+sum(sign(abs(B1*B2*Delta)))+
              sum(sign(abs(R2N2))))*log(ob)
        Bic[i]=-b1+b2
        break
      }
      P=Q
    }
  }
  bic_value=sum(Bic)
  return(list(res_1,res_2,res_3,iter,bic_value))
}

