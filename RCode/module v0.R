
setwd("C:/Github/Study_CrossoverDesign/RCode")
#=============================================
#simulation parameter
#=============================================
sim_time=3000;cor_par=1/6
seq_size=c(25,50,100,200,300)
cros_type=c('ABBA','ABBBAA','AABABABAA','ABCBCACAB','BACACBBCA','BBAACBCAC')
#=============================================
#true value of params with treat-seq-time
#=============================================
param_222=c(1.0,0.7,0.3,0.2);param_223=c(1.0,0.7,0.3,0.3,0.2);param_233=c(1.0,0.7,0.3,0.3,0.2,0.2);param_333=c(1.0,0.7,0.7,0.3,0.3,0.2,0.2)
#=================================================
#Link of Yist:Y11,Y12,Y21,Y22
#row is vector of param,column is x.mat of Yist
#=================================================
#ABBA
xmat_222=matrix(c(1,1,1,1, 0,1,1,0, 0,1,0,1, 0,0,1,1), nrow = 4, ncol = 4,byrow = TRUE)
#ABBBAA
xmat_223=matrix(c(1,1,1,1,1,1, 0,1,1,1,0,0, 0,1,0,0,1,0, 0,0,1,0,0,1, 0,0,0,1,1,1), nrow = 5, ncol = 6,byrow = TRUE)
#AABABABAA
xmat_233=matrix(c(1,1,1,1,1,1,1,1,1, 0,0,1,0,1,0,1,0,0, 0,1,0,0,1,0,0,1,0, 0,0,1,0,0,1,0,0,1, 0,0,0,1,1,1,0,0,0, 0,0,0,0,0,0,1,1,1), nrow = 6, ncol = 9,byrow = TRUE)
#ABCBCACAB
xmat_333.1=matrix(c(1,1,1,1,1,1,1,1,1, 0,1,0,1,0,0,0,0,1, 0,0,1,0,1,0,1,0,0, 0,1,0,0,1,0,0,1,0, 0,0,1,0,0,1,0,0,1, 0,0,0,1,1,1,0,0,0, 0,0,0,0,0,0,1,1,1), nrow = 7, ncol = 9,byrow = TRUE)
#BACACBBCA
xmat_333.2=matrix(c(1,1,1,1,1,1,1,1,1, 1,0,0,0,0,1,1,0,0, 0,0,1,0,1,0,0,1,0, 0,1,0,0,1,0,0,1,0, 0,0,1,0,0,1,0,0,1, 0,0,0,1,1,1,0,0,0, 0,0,0,0,0,0,1,1,1), nrow = 7, ncol = 9,byrow = TRUE)
#BBAACBCAC
xmat_333.3=matrix(c(1,1,1,1,1,1,1,1,1, 1,1,0,0,0,1,0,0,0, 0,0,0,0,1,0,1,0,1, 0,1,0,0,1,0,0,1,0, 0,0,1,0,0,1,0,0,1, 0,0,0,1,1,1,0,0,0, 0,0,0,0,0,0,1,1,1), nrow = 7, ncol = 9,byrow = TRUE)
#========================================
#Function of output object
#========================================
#Par.values==TrueMean
Mean.True<-function(Par.values,x.mat){  return(exp(Par.values%*%x.mat))}

#I.true
Matrix.I<-function(cors.type,params,x.mat){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  mat.I<- matrix(0, nrow = length(params), ncol = length(params))
  for (i in 1:length(params) ){ mat.I[i,]<-Mean%*%(t(x.mat)*x.mat[i,])/num.seq}
  return(mat.I)
}
#======================================================================
#Generate data
#======================================================================

Data.ind<-function(cros.type,mean.true,seq.size){
    num.seq<-if (nchar(cros.type)==9) 3 else 2
    data <-  matrix(0, nrow = seq.size, ncol = nchar(cros.type))
    for (j in 1:length(mean.true)) {data[,j] =rpois(seq.size, lambda = mean.true[1,j])}
    return(data)
  }
#--------------------------------------------------------------------
Data.cor<-function(cros.type,mean.true,seq.size,cor.par){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  nui<-replicate(2,rgamma(n=seq.size,shape=1/cor.par,scale=cor.par))
  mean.seq<-matrix(mean.true,ncol=nchar(cros.type)/num.seq,nrow=num.seq)
  mean.cor<-list();data <-matrix(0, nrow = seq.size, ncol = nchar(cros.type))
  for (i in 1:num.seq){
    m<-split(c(outer(nui[,i],mean.seq[,i],  function(x, y) x * y)), ceiling(seq_along(c(outer(nui[,i],mean.seq[,i],  function(x, y) x * y)))/seq.size))
    mean.cor<-append(mean.cor,m)
    }
  mean.cor<- t(matrix(unlist(mean.cor), ncol = seq.size, byrow = TRUE))
  for (i in 1:nchar(cros.type)){
    list_poisson <- unlist(lapply(mean.cor[,i], FUN = function(x, y) rpois(y, x), y = 1))
    data[,i]<-list_poisson
  }
  return(data)
}
#--------------------------------------------------------------------
#get correlation matrix
data.cor<-Data.cor(cros.type = cros_type[1],mean.true =Mean.True(param_222,xmat_222),seq.size =100000,cor.par = cor_par )
cor(data.cor)
#======================================================================
#MLE of different type of crossover design
#======================================================================
##param<-c(tao,eta,gamma,delta)
MLE.ABBA<-function(data,seq.size){
  y.sum=colSums(data)
  tao.hat=log(y.sum[1]/seq.size)
  eta.hat=0.5*(log(y.sum[2])+log(y.sum[3])-log(y.sum[1])-log(y.sum[4]))
  gamma.hat=0.5*(log(y.sum[2])+log(y.sum[4])-log(y.sum[1])-log(y.sum[3]))
  delta.hat=0.5*(log(y.sum[3])+log(y.sum[4])-log(y.sum[1])-log(y.sum[2]))
  MLE.i = matrix(c(tao.hat, eta.hat, gamma.hat, delta.hat), nrow=1,ncol=4)
  return(MLE.i)
}
##param<-c(tao,eta,gamma1,gamma2,delta)
MLE.ABBBAA<-function(data,seq.size){
  y.sum=colSums(data)
  tao.hat=log(y.sum[1]/seq.size)
  eta_hat =0.25*(log(y.sum[2])+log(y.sum[3])+2*log(y.sum[4])-log(y.sum[5])-log(y.sum[6])-2*log(y.sum[1]))
  gamma1_hat = 0.5*(log(y.sum[2])+log(y.sum[5])-log(y.sum[4])-log(y.sum[1]))
  gamma2_hat = 0.5*(log(y.sum[3])+log(y.sum[6])-log(y.sum[4])-log(y.sum[1]))
  delta_hat = 0.5*(log(y.sum[4])-log(y.sum[1]))-0.25*(log(y.sum[2])+log(y.sum[3])-log(y.sum[4])-log(y.sum[6]))
  MLE.i = matrix(c(tao.hat, eta.hat, gamma1.hat, gamma2.hat, delta.hat), nrow=1,ncol=5)
  return(MLE.i)
}

##param<-c(tao,eta,gamma1,gamma2,delta1,delta2)
NegLL.AABABABAA <- function(param) { sum(factorial(data))
  -sum(param[1]*data[,1]-exp(param[1])+(param[1]+param[3])*data[,2]-exp(param[1]+param[3])+(param[1]+param[2]+param[4])*data[,3]-exp(param[1]+param[2]+param[4])
       +(param[1]+param[5])*data[,4]-exp(param[1]+param[5])+(sum(param[1:3])+param[5])*data[,5]-exp(sum(param[1:3])+param[5])+(param[1]+param[4]+param[5])*data[,6]-exp(param[1]+param[4]+param[5])
       +(param[1]+param[2]+param[6])*data[,7]-exp(param[1]+param[2]+param[6])+(param[1]+param[3]+param[5])*data[,8]-exp(param[1]+param[3]+param[5])+(param[1]+param[4]+param[6])*data[,9]-exp(param[1]+param[4]+param[6])
  )}
MLE.AABABABAA<-function(data,seq.size){
  AABABABAA<-optim(param <- c(0.95, 0.65, 0.25, 0.25, 0.15, 0.15), NegLL.AABABABAA, hessian=TRUE)
  return(AABABABAA$par)
}

##param<-c(tao,eta1,eta2,gamma1,gamma2,delta1,delta2)
NegLL.ABCBCACAB <- function(param) { sum(factorial(data))
  -sum(param[1]*data[,1]-exp(param[1])+(sum(param[1:3]))*data[,2]-exp(sum(param[1:3]))+(param[1]+param[3]+param[5])*data[,3]-exp(param[1]+param[3]+param[5])
       +(param[1]+param[2]+param[6])*data[,4]-exp(param[1]+param[2]+param[6])+(param[1]+parma[3]+param[4]+param[6])*data[,5]-exp(param[1]+parma[3]+param[4]+param[6])+(param[1]+param[5]+param[6])*data[,6]-exp(param[1]+param[5]+param[6])
       +(param[1]+param[3]+param[7])*data[,7]-exp(param[1]+param[3]+param[7])+(param[1]+param[4]+param[7])*data[,8]-exp(param[1]+param[4]+param[7])+(param[1]+param[2]+param[5]+param[7])*data[,9]-exp(param[1]+param[2]+param[5]+param[7])
  )}
MLE.ABCBCACAB<-function(data,seq.size){
  ABCBCACAB<-optim(param <- c(0.95, 0.65, 0.65, 0.25, 0.25, 0.15, 0.15), NegLL.ABCBCACAB, hessian=TRUE)
  return(ABCBCACAB$par)
}
##param<-c(tao,eta1,eta2,gamma1,gamma2,delta1,delta2)
NegLL.BACACBBCA <- function(param) { sum(factorial(data))
  -sum( sum(param[1:2])*data[,1]-exp(sum(param[1:2]))+(param[1]+param[4])*data[,2]-exp(param[1]+param[4])+(param[1]+param[3]+param[5])*data[,3]-exp(param[1]+param[3]+param[5])
       +(param[1]+param[6])*data[,4]-exp(param[1]+param[6])+(param[1]+parma[3]+param[4]+param[6])*data[,5]-exp(param[1]+parma[3]+param[4]+param[6])+(param[1]+param[2]+param[5]+param[6])*data[,6]-exp(param[1]+param[2]+param[5]+param[6])
       +(param[1]+param[2]+param[7])*data[,7]-exp(param[1]+param[2]+param[7])+(param[1]+param[3]+param[4]+param[7])*data[,8]-exp(param[1]+param[3]+param[4]+param[7])+(param[1]+param[5]+param[7])*data[,9]-exp(param[1]+param[5]+param[7])
  )}
MLE.BACACBBCA<-function(data,seq.size){
  BACACBBCA<-optim(param <- c(0.95, 0.65, 0.65, 0.25, 0.25, 0.15, 0.15), NegLL.BACACBBCA, hessian=TRUE)
  return(BACACBBCA$par)
}
##param<-c(tao,eta1,eta2,gamma1,gamma2,delta1,delta2)
NegLL.BBAACBCAC <- function(param) { sum(factorial(data))
  -sum( sum(param[1:2])*data[,1]-exp(sum(param[1:2]))+(param[1]+param[2]+param[4])*data[,2]-exp(param[1]+param[2]+param[4])+(param[1]+param[5])*data[,3]-exp(param[1]+param[5])
        +(param[1]+param[6])*data[,4]-exp(param[1]+param[6])+(param[1]+parma[3]+param[4]+param[6])*data[,5]-exp(param[1]+parma[3]+param[4]+param[6])+(param[1]+param[2]+param[5]+param[6])*data[,6]-exp(param[1]+param[2]+param[5]+param[6])
        +(param[1]+param[3]+param[7])*data[,7]-exp(param[1]+param[3]+param[7])+(param[1]+param[4]+param[7])*data[,8]-exp(param[1]+param[4]+param[7])+(param[1]+param[3]+param[5]+param[7])*data[,9]-exp(param[1]+param[3]+param[5]+param[7])
  )}
MLE.BBAACBCAC<-function(data,seq.size){
  BBAACBCAC<-optim(param <- c(0.95, 0.65, 0.65, 0.25, 0.25, 0.15, 0.15), NegLL.BBAACBCAC, hessian=TRUE)
  return(BBAACBCAC$par)
}
#--------------------------------------------------------------------
#I.hat & V.hat
#--------------------------------------------------------------------
Matrix.IV<-function(cros.type,mle.values,x.mat,seq.size,data){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  mat.I<- matrix(0, nrow = length(mle.values), ncol = length(mle.values))
  mat.V<- matrix(0, nrow = length(mle.values), ncol = length(mle.values))
  mat.score<- matrix(0, nrow = seq.size, ncol = length(mle.values))
  mean.est<-exp(mle.values%*%x.mat)
  for (i in 1:length(mle.values) ){ 
    mat.I[i,]<-mean.est%*%(t(x.mat)*x.mat[i,])/num.seq
    mat.score[,i]<-sweep(data, 2, mean.est[1,])%*%x.mat[i,]
    mat.V[i,]<- colSums(mat.score*mat.score[,i])/(seq.size*num.seq)
  }
  mat.V<-Matrix::forceSymmetric(mat.V,uplo="L")
  list.IV <- list("I.hat" = mat.I, "V.hat" =mat.V)
  return(list.IV) 
}
#========================================
#format that output to excel
#========================================
#obj=obj.ind;num.param=length(param_222);MATS=MATS.ind
Output.Format<-function(obj,num.param,MATS){
  space<-matrix(c('','',''), nrow=1,ncol = length(c('','','')))
  obj.names<-matrix(0,ncol = num.param,length(obj))
  for (i in 1:length(obj)){obj.names[i,]<-cbind(obj[i],space)}
  #format result that we want to show in excel
  seq.result<-matrix(c('tao_hat','eta_hat','gamma_hat','delta_hat'),ncol = num.param,nrow = 1)
  for (i in 1:length(obj)){
    temp<-rbind(obj.names[i,],MATS[[i]])
    seq.result<-rbind(seq.result,temp)
    }
  
  return(seq.result)
}
#--------------------------------------------------------------------
#========================================
#main
#result:to store result of each seq_size
#========================================
result.ind <- list()
result.cor <- list()
#simulation for ABBA

for (seq in seq_size){
  
  MLE.ind<-matrix(0, nrow = sim_time, ncol = length(param_222))
  MLE.cor<-matrix(0, nrow = sim_time, ncol = length(param_222))
  I.ind<- 0 ; I.cor<- 0 ; V.ind<-0 ;V.cor<-0;invI.ind<-0;invI.cor<-0
  set.seed(7354)
  for (i in 1:sim_time){
    
    data.ind<-Data.ind(cros.type = cros_type[1],mean.true =Mean.True(param_222,xmat_222),seq.size = seq )
    MLE.ind.i<-MLE.ABBA(data.ind,seq.size=seq)
    mean.est<-Mean.True(MLE.ind.i,xmat_222)
    IV.ind.i<-Matrix.IV(cros.type=cros_type[1], mle.values=MLE.ind.i, x.mat=xmat_222, seq.size=seq, data=data.ind)
    
    data.cor<-Data.cor(cros.type = cros_type[1],mean.true =Mean.True(param_222,xmat_222),seq.size =seq,cor.par = cor_par )
    MLE.cor.i<-MLE.ABBA(data.cor,seq.size=seq)
    mean.est<-Mean.True(MLE.cor.i,xmat_222)
    IV.cor.i<-Matrix.IV(cros.type=cros_type[1], mle.values=MLE.cor.i, x.mat=xmat_222, seq.size=seq, data=data.cor)
    
    
    MLE.ind[i,]<-MLE.ind.i
    I.ind.i<-IV.ind.i$I.hat
    V.ind.i<-Matrix::forceSymmetric(IV.ind.i$V.hat,uplo="L")
    invI.ind.i<-solve(I.ind.i)
    I.ind<-I.ind+I.ind.i
    V.ind<-V.ind+V.ind.i
    invI.ind<-invI.ind+invI.ind.i
    
    MLE.cor[i,]<-MLE.cor.i
    I.cor.i<-IV.cor.i$I.hat
    V.cor.i<-Matrix::forceSymmetric(IV.cor.i$V.hat,uplo="L")
    invI.cor.i<-solve(I.cor.i)
    I.cor<-I.cor+I.cor.i
    V.cor<-V.cor+V.cor.i
    invI.cor<-invI.ind+invI.cor.i
    
    
    }
  
  #IVI.cor<-solve(I.cor/sim_time)%*%(as.matrix(V.cor)/sim_time)%*%solve(I.cor/sim_time)
  IVI.cor<-(invI.cor/sim_time)%*%(as.matrix(V.cor)/sim_time)%*%(invI.cor/sim_time)
  #store result in MATS for seq
  MATS.ind <- list(signif(t(as.matrix(colMeans(MLE.ind))),5),signif(I.ind/sim_time,5),signif(as.matrix(V.ind)/sim_time,5),signif(cov(MLE.ind)*seq*2,5),signif(invI.ind/sim_time,5))
  MATS.cor <- list(signif(t(as.matrix(colMeans(MLE.cor))),5),signif(I.cor/sim_time,5),signif(as.matrix(V.cor)/sim_time,5),signif(cov(MLE.cor)*seq*2,5),signif(IVI.cor,5))
  
  obj.ind<-c('MLE.closeform','I.hat.ind','V.hat.ind','NS.MLE.closeform.ind','inv(I.hat.ind)')
  obj.cor<-c('MLE.closeform','I.hat.cor','V.hat.cor','NS.MLE.closeform.cor','inv(I.hat.cor).V.hat.cor.inv(I.hat.cor)')
  
  #format result that we want to show in excel
  seq.res.ind<-Output.Format(obj=obj.ind,num.param=length(param_222),MATS=MATS.ind)
  seq.res.cor<-Output.Format(obj=obj.cor,num.param=length(param_222),MATS=MATS.cor)
  #store result of all seq_size
  num=seq_size[which(seq_size == seq)]
  result.ind <- append(result.ind, list( num = seq.res.ind))
  result.cor <- append(result.cor, list( num = seq.res.cor))
  
}

library(xlsx)
write.xlsx(result.ind, file = paste0(Sys.Date(),'ABBAind.xlsx'))
write.xlsx(result.cor, file = paste0(Sys.Date(),'ABBAcor.xlsx'))
