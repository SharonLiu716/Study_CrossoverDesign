rm(list=ls(all=TRUE)) 
setwd("C:/Users/User/Documents/Study_CrossoverDesign/RCode")
#=============================================
#simulation parameter
#=============================================
sim_time=5000;cor_par=1
seq_size=c(25,50,100,200,300)
cros_type=c('ABBA','ABBBAA','AABABABAA','ABCBCACAB','BACACBBCA','BBAACBCAC')
#=============================================
#true value of params with treat-seq-time
#=============================================
param_222=c(0.5,1.2,-0.7,0.2)
#=================================================
#Link of Yist:Y11,Y12,Y21,Y22
#row is vector of param,column is x.mat of Yist
#=================================================
#ABBAlog(0.5)
xmat_222=matrix(c(1,1,1,1, 0,1,1,0, 0,1,0,1, 0,0,1,1), nrow = 4, ncol = 4,byrow = TRUE)
#========================================
#Function of output object
#========================================
#Par.values==TrueMean
Mean.True<-function(Par.values,x.mat){  return(exp(Par.values%*%x.mat))}
#Mean.True(Par.values=param_222,x.mat=xmat_222)
#I.true
Matrix.I<-function(cros.type,params,x.mat){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  mat.I<- matrix(0, nrow = length(params), ncol = length(params))
  Mean<-Mean.True(Par.values=params,x.mat=x.mat)
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

#======================================================================
#MLE of different type of crossover design
#======================================================================
##param<-c(tao,eta,gamma,delta)
MLE.ABBA<-function(data){
  y.sum=colMeans(data)
  tao.hat=log(y.sum[1])
  eta.hat=0.5*(log(y.sum[2])+log(y.sum[3])-log(y.sum[1])-log(y.sum[4]))
  gamma.hat=0.5*(log(y.sum[2])+log(y.sum[4])-log(y.sum[1])-log(y.sum[3]))
  delta.hat=0.5*(log(y.sum[3])+log(y.sum[4])-log(y.sum[1])-log(y.sum[2]))
  MLE.i = matrix(c(tao.hat, eta.hat, gamma.hat, delta.hat), nrow=1,ncol=4)
  return(MLE.i)
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
#True I、correlation、true V
#========================================
set.seed(7353)
ABBA.I<-Matrix.I(cros.type=cros_type[1],params=param_222,x.mat=xmat_222)
#get correlation matrix
df.ind<-Data.ind(cros.type = cros_type[1],mean.true =Mean.True(param_222,xmat_222),seq.size =100000 )
df.cor<-Data.cor(cros.type = cros_type[1],mean.true =Mean.True(param_222,xmat_222),seq.size =100000,cor.par = cor_par )
cor(df.ind)
cor(df.cor)
Matrix.IV(cros.type=cros_type[1],mle.values=param_222,x.mat=xmat_222,seq.size=100000,data=df.ind)
Matrix.IV(cros.type=cros_type[1],mle.values=param_222,x.mat=xmat_222,seq.size=100000,data=df.cor)


#===================================================================
#Loglikelihood of ABBA,MLE.ABBA under H0,Matrix AB
#===================================================================
loglik.ABBA<-function(param,data){ 
  
   sum( param[1]*data[,1]-exp(param[1]))
   +sum(sum(param[1:3])*data[,2]-exp(sum(param[1:3])))
   +sum((param[1]+param[2]+param[4])*data[,3]-exp(param[1]+param[2]+param[4]))
  +sum((param[1]+param[3]+param[4])*data[,4]-exp(param[1]+param[3]+param[4]))
   }


MLE.ABBAnull<-function(data,seq.size,eta.null){
  y.sum=colMeans(data)
  tao.hatnull<-log(y.sum[1])
  gamma.hatnull=log(y.sum[2])-log(y.sum[1])-eta.null
  delta.hatnull=log(y.sum[4])-log(y.sum[2])+eta.null
  MLE.i.null = matrix(c(tao.hatnull,  eta.null, gamma.hatnull, delta.hatnull), nrow=1,ncol=4)
  return(MLE.i.null)
}
Matrix.AB<-function(Mat.I,Mat.V){
  Mat.I<-I.ind.i
  diag.I<-as.matrix(diag(Mat.I),byrow = TRUE)
  diag.V<-as.matrix(diag(Mat.V),byrow = TRUE)
  #I.theta.theta
  I.tt<-diag.I[2];V.tt<-diag.V[2];I.pp<-Mat.I[-c(2),-c(2)];V.pp<-Mat.V[-c(2),-c(2)]
  I.tp<-as.matrix(Mat.I[2,]);I.tp<-as.matrix(I.tp[-c(2),]);V.tp<-as.matrix(Mat.V[2,]);V.tp<-as.matrix(V.tp[-c(2),])
  #mat.A
  Mat.A<-I.tt-t(I.tp)%*%solve(I.pp)%*%I.tp
  Mat.B<-V.tt-2*t(I.tp)%*%solve(I.pp)%*%V.tp+t(I.tp)%*%solve(I.pp)%*%V.pp%*%solve(I.pp)%*%I.tp
  Mat.AB<- list("Mat.A" = Mat.A, "Mat.B" =Mat.B)
  return(Mat.AB)
}

#========================================
#main
#result:to store result of each seq_size
#========================================
result.ind <- list();result.cor <- list();pvalue.ind<-list();pvalue.cor<-list()
#simulation for ABBA

for (seq in seq_size){
  
  MLE.ind<-matrix(0, nrow = sim_time, ncol = length(param_222))
  MLE.cor<-matrix(0, nrow = sim_time, ncol = length(param_222))
  I.ind<- 0 ; I.cor<- 0 ; V.ind<-0 ;V.cor<-0;invI.ind<-0;invI.cor<-0
  # Create df for statistics
  statistics.columns <- c("Wald.na","Wald.rb","LR.na","LR.rb") 
  df.statistics.ind <- data.frame(matrix(nrow = sim_time, ncol = length(statistics.columns))) 
  colnames(df.statistics.ind) <- statistics.columns
  df.statistics.cor = data.frame(matrix(nrow = sim_time, ncol = length(statistics.columns))) 
  colnames(df.statistics.cor) <- statistics.columns
  seq=100
  set.seed(7353)
  for (i in 1:sim_time){
    #====================================================
    #I,V,IVI
    #====================================================
    #independent
    data.ind<-Data.ind(cros.type = cros_type[1],mean.true =Mean.True(param_222,xmat_222),seq.size = seq )
    MLE.ind.i<-MLE.ABBA(data.ind)
    mean.est<-Mean.True(MLE.ind.i,xmat_222)
    IV.ind.i<-Matrix.IV(cros.type=cros_type[1], mle.values=MLE.ind.i, x.mat=xmat_222, seq.size=seq, data=data.ind)
    
    #correlated
    data.cor<-Data.cor(cros.type = cros_type[1],mean.true =Mean.True(param_222,xmat_222),seq.size =seq,cor.par = cor_par )
    MLE.cor.i<-MLE.ABBA(data.cor)
    mean.est<-Mean.True(MLE.cor.i,xmat_222)
    IV.cor.i<-Matrix.IV(cros.type=cros_type[1], mle.values=MLE.cor.i, x.mat=xmat_222, seq.size=seq, data=data.cor)
    
    #store result of MLE,I,V,inv.I
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
    #=========================================================================
    #output:matrix AB, Wald statistics,LR statistics, Score statistics
    #=========================================================================
    matA.ind.i<-Matrix.AB(I.ind.i,as.matrix(V.ind.i))$Mat.A
    matB.ind.i<-Matrix.AB(I.ind.i,as.matrix(V.ind.i))$Mat.B
    matB.ind.i/(matA.ind.i^2)
    theta.null=0
    df.statistics.ind[i, "Wald.na"] <- seq*2*matA.ind.i*(MLE.ind.i[,2]-theta.null)^2
    df.statistics.ind[i, "Wald.rb"] <- seq*2*(matA.ind.i^2)/matB.ind.i*((MLE.ind.i[,2]-theta.null)^2)
    MLE.null<-MLE.ABBAnull(data = data.ind, seq.size = seq, eta.null=theta.null)
    df.statistics.ind[i, "LR.na"] <-2*(loglik.ABBA(param = MLE.ind.i,data = data.ind)-loglik.ABBA(param = as.vector(MLE.null),data = data.ind))
    df.statistics.ind[i, "LR.rb"] <-2*matA.ind.i/matB.ind.i*(loglik.ABBA(param = MLE.ind.i,data = data.ind)-loglik.ABBA(param = as.vector(MLE.null),data = data.ind))
    #wald.test(Sigma = cov(MLE.ind), b = t(MLE.ind.i), Terms = 2)
    matA.cor.i<-Matrix.AB(I.cor.i,as.matrix(V.cor.i))$Mat.A
    matB.cor.i<-Matrix.AB(I.cor.i,as.matrix(V.cor.i))$Mat.B
    df.statistics.cor[i, "Wald.na"] <- matA.cor.i*(MLE.cor.i[,2]-theta.null)^2
    df.statistics.cor[i, "Wald.rb"] <- (matA.cor.i^2)/matB.cor.i*((MLE.cor.i[,2]-theta.null)^2)
    MLE.null<-MLE.ABBAnull(data = data.cor, seq.size = seq, eta.null=theta.null)
    df.statistics.cor[i, "LR.na"] <-2*(loglik.ABBA(param = MLE.cor.i,data = data.cor)-loglik.ABBA(param = as.vector(MLE.null),data = data.cor))
    df.statistics.cor[i, "LR.rb"] <-2*matA.cor.i*(loglik.ABBA(param = MLE.cor.i,data = data.cor)-loglik.ABBA(param = as.vector(MLE.null),data = data.cor))/matB.cor.i
    
    
    # #Wald.naive
    # Wald.na.ind<-seq.size*2*Mat.A*(theta.hat-theta.null)^2
    # #Wald.robust
    # Wald.rb<-seq.size*2*(Mat.A^2)*((theta.hat-theta.null)^2)/Mat.B
    # #LR.naive
    # MLE.null<-MLE.ABBAnull(data = data.ind, seq.size = seq,eta.null=0)
    # LR.na<-2*(loglik.ABBA(param = MLE.ind.i,data = data.ind)-loglik.ABBA(param = as.vector(MLE.null),data = data.ind))
    # #LR.robust
    # LR.rb<-2*Mat.A*(loglik.ABBA(param = MLE.ind.i,data = data.ind)-loglik.ABBA(param = as.vector(MLE.null),data = data.ind))/Mat.B
    # 
    
    }
  
  #=========================================
  #output for each seq.size:I, V, IVI
  #=========================================
  IVI.cor<-solve(I.cor/sim_time)%*%(as.matrix(V.cor)/sim_time)%*%solve(I.cor/sim_time)
  #IVI.cor<-(invI.cor/sim_time)%*%(as.matrix(V.cor)/sim_time)%*%(invI.cor/sim_time)
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
  
  #testing
  
  below <- function(x) {return(length(x[x>qchisq(0.95,1)])/sim_time)}
  print(paste('ind',seq,sapply(df.statistics.ind, below)))
  pvalue.ind<-append(pvalue.ind,sapply(df.statistics.ind, below))
  print(paste('cor',seq,sapply(df.statistics.cor, below)))
  pvalue.cor<-append(pvalue.cor,sapply(df.statistics.cor, below))
}

library(xlsx)
write.xlsx(result.ind, file = paste0(Sys.Date(),'ABBAind.xlsx'))
write.xlsx(result.cor, file = paste0(Sys.Date(),'ABBAcor.xlsx'))


