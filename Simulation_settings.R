############ Codes for generating simulated data 
library(MASS)
library(Matrix)
library(igraph)
##################################################
#Section1: setting design matrices and three kinds 
#of precision matrices of residual.

P_vec=c(20,60,100) # the dimension of variables
p=P_vec[1]
ob=200 # the number of observations
b=p/20 # the number of blocks
A=matrix(0,20,20) # A and B represent two design matrices
B=matrix(0,20,20)
M1=matrix(0,p,p)
M2=matrix(0,p,p)

for (k in 1:b) {
  A[6:12,6:12]=runif(n=49,min=0.4,max=1)
  #three kinds of overlapping regions for B:[6:12,6:12];[8:14,8:14];[10:16,10:16]
  #the case 1 of overlap
  B[6:12,6:12]=runif(n=49,min=-1,max=-0.4)
  M1[(20*k-19):(20*k),(20*k-19):(20*k)]=A
  M2[(20*k-19):(20*k),(20*k-19):(20*k)]=B
}

#three kinds of precison matrices of residual
mean=rep(0,p)
sigma=diag(p)
####setting 1: Power-law structure
pre_W_0=diag(p)
p2=7# region of overlap with dimensions 3/5/7
for (k in 1:b) {
  g=ba.game(p2,m=2,directed=F)
  Eg=as.data.frame(get.edgelist(g))
  subi=matrix(0,p2,p2)
  umin=0.5 # the signal strength 
  umax=1
  for (q in 1:dim(Eg)[1]) {
    i=Eg[q,1];j=Eg[q,2]
    ij=sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
    subi[i,j]=ij;subi[j,i]=ij
  }
  for (i in 1:p2) {
    subi[i,i] = sum(abs(subi[i,setdiff(1:p2,i)]))+0.1
  }
  #three cases of overlap
  pre_W_0[(20*k-14):(20*k-8),(20*k-14):(20*k-8)]=subi
  #pre_W_0[(20*k-12):(20*k-8),(20*k-12):(20*k-8)]=subi
  #pre_W_0[(20*k-10):(20*k-8),(20*k-10):(20*k-8)]=subi
}
W_1=solve(pre_W_0)

#####setting 2: Nearest neighbor structure
umin=0.2 # the signal strength
umax=0.4
p2=7 # the region of overlap with dimensions 3/5/7
pre_W_0=diag(p)
near=3 #the number of neighbors  1/2/3
for (k in 1:b) {
  subi = diag(1,p2,p2)
  point.matrix = matrix(runif(2*p2,0,1),p2,2)
  num.nearest = matrix(0,p2,near)
  for (j in 1:p2) {
    distancej = apply((t(point.matrix[-j,]) - point.matrix[j,])^2,2,sum)
    corj = setdiff(1:p,j)[order(distancej)[1:near]]
    num.nearest[j,] = corj
  }
  for (oi in 1:dim(num.nearest)[1]) {
    for (oj in 1:dim(num.nearest)[2]) {
      i=oi;j=num.nearest[oi,oj]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
  }
  for (i in 1:p2) {
    subi[i,i] = sum(abs(subi[i,setdiff(1:p2,i)]))+0.1
  }
  #three cases of overlap
  pre_W_0[(20*k-14):(20*k-8),(20*k-14):(20*k-8)]=subi
  #pre_W_0[(20*k-12):(20*k-8),(20*k-12):(20*k-8)]=subi
  #pre_W_0[(20*k-10):(20*k-8),(20*k-10):(20*k-8)]=subi
}
W_2=solve(pre_W_0)

#####setting 3: Banded structure
W_0=matrix(0,p,p)
for (i in 1:p) {
  for (j in 1:p) {
    W_0[i,j]=(0.6)^(abs(i-j)) # the signal strength
  }
}
W_3=W_0*sign(abs(M1))*sign(abs(M2))
diag(W_3)=1
  
#########################
###Section 2: generating simulated data 
W=W_1

#generating data for methods to calculate
dataZ1=mvrnorm(n=ob,mean,sigma)
dataZ2=mvrnorm(n=ob,mean,sigma)#Z2
dataW=mvrnorm(n=ob,mean,W)
dataXp=dataZ1%*%M1+dataZ2%*%M2+dataW

#generating 500000 to approximate the real precision matrices
data_Z1=mvrnorm(n=500000,mean,sigma)
data_Z2=mvrnorm(n=500000,mean,sigma)
data_W=mvrnorm(n=500000,mean,W)
TrueData=data_Z1%*%M1+data_Z2%*%M2+data_W

int_XZ1=solve(cov(cbind(TrueData,data_Z1)))[1:p,1:p]
int_XZ2=solve(cov(cbind(TrueData,data_Z2)))[1:p,1:p]
int_XZ1Z2=solve(cov(cbind(TrueData,data_Z1,data_Z2)))[1:p,1:p]
real_XZ1=matrix(0,p,(p-1))
real_XZ2=matrix(0,p,(p-1))
real_XZ1Z2=matrix(0,p,(p-1))
for(i in 1:p){
  real_XZ1[i,]=int_XZ1[i,-i]
  real_XZ2[i,]=int_XZ2[i,-i]
  real_XZ1Z2[i,]=int_XZ1Z2[i,-i]
}
real_XZ1[abs(real_XZ1)<=0.01]=0
real_XZ2[abs(real_XZ2)<=0.01]=0
real_XZ1Z2[abs(real_XZ1Z2)<=0.05]=0
XZ1_sps=sign(abs(real_XZ1))
XZ2_sps=sign(abs(real_XZ2))
XZ1Z2_sps=sign(abs(real_XZ1Z2))
rm(data_Z1,data_Z2,data_W,TrueData)
gc()
