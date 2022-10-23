set.seed(7353)
setwd("C:/Github/Study_CrossoverDesign/RCode")
sim_time=2000
seq_size<-c(25,50,100,200)
cor_param=0.1
cros_type='ABBA'
#tao,eta,gamma,delta
params = matrix(c(1.0, 0.7, 0.3, 0.2), nrow=1,ncol=4)
#link of mean.i11,mean.mean.i12,mean.i21,mean.22
x.mat<-matrix(c(1,1,1,1,0,1,1,0,0,1,0,1,0,0,1,1), nrow = 4, ncol = 4,byrow=TRUE)
mean_true=exp(params%*%x.mat)
#(seq_size)/(seq_size*2)
i.tt<-sum(mean_true)/2
i.ee<-(mean_true[2]+mean_true[3])/2
i.gg<-(mean_true[2]+mean_true[4])/2
i.dd<-(mean_true[3]+mean_true[4])/2
i.eg<-mean_true[2]/2
i.ed<-mean_true[3]/2
i.gd<-mean_true[4]/2
II<-matrix(c(i.tt,i.ee,i.gg,i.dd,i.ee,i.ee,i.eg,i.ed,
             i.gg,i.eg,i.gg,i.gd,i.dd,i.ed,i.gd,i.dd), nrow = 4, ncol = 4,byrow=TRUE)
I.true<- matrix(0, nrow = 4, ncol = 4)
for (i in 1:length(params) ){ I.true[i,]<-0.5*mean_true%*%(t(x.mat)*x.mat[i,])}
#the form of matrix that export to excel
#output object
obj<-c('MLE.optim','MLE.closeform','I.hat.hessian','I.hat','V.hat.matix','NS.MLE.optim','NS.MLE.closeform','inv(I.hat)')
#name of output object
space<-matrix(c('','',''), nrow=1,ncol = 3)
obj.names<-matrix(0,ncol = 4,length(obj))
for (i in 1:length(obj)){obj.names[i,]<-cbind(obj[i],space)}
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

#================================================
#Matrix V
#================================================


#=========================
#outlier detect
#=========================
library(data.table)
library(dplyr)
del.range<-c(.01, .99)
remove.outlier<-function(yi,mean){
  data.new<-yi[do.call(between, c(list(yi), quantile(yi, del.range, names=F)))]
  while (length(data.new)!=length(yi)){data.new<-append(data.new, rpois(length(yi)-length(data.new), lambda = mean))}
  return(data.new) }
#==========================
#main
#==========================
set.seed(7353)
result <- list()#to store result of each seq_size
#use optim to obtain MLE、I_hat、V_hat
for (seq in seq_size){
  MLE.optim<-matrix(0, nrow = sim_time, ncol = 4)
  I.optim=0
  MLE.closeform<-matrix(0, nrow = sim_time, ncol = 4)
  I.hat<- 0;V.hat<- 0;invI.hat<-0
  #seq=50
  for (i in 1:sim_time){
    pi=0.5
    data <-  matrix(0, nrow = seq, ncol = 4)
    for (j in 1:length(mean_true)) { data[,j] =rpois(seq, lambda = mean_true[1,j])#}
       data[,j]<-remove.outlier(yi =data[,j], mean = mean_true[1,j])}
    #for (i in 1:4){  data[,i]<-remove.outlier(yi =data[,i], mean = mean_true[,i])}
    negll <- function(param) { sum(factorial(data))
       -sum( param[1]*data[,1]-exp(param[1])+sum(param[1:3])*data[,2]-exp(sum(param[1:3]))+(param[1]+param[2]+param[4])*data[,3]-exp(param[1]+param[2]+param[4])+(param[1]+param[3]+param[4])*data[,4]-exp(param[1]+param[3]+param[4]) )}
     ABBA<-optim(param <- c(0.95, 0.65, 0.25, 0.15), negll, hessian=TRUE)
     MLE.optim[i,]<-ABBA$par
     I.optim <-I.optim+ABBA$hessian
    y.sum=colMeans(data)
    tao.hat=log(y.sum[1])
    eta.hat=(log(y.sum[2])+log(y.sum[3])-log(y.sum[1])-log(y.sum[4]))/2
    gamma.hat=(log(y.sum[2])-log(y.sum[1])+log(y.sum[4])-log(y.sum[3]))/2
    delta.hat=(log(y.sum[3])+log(y.sum[4])-log(y.sum[1])-log(y.sum[2]))/2
    MLE.i = matrix(c(tao.hat, eta.hat, gamma.hat, delta.hat), nrow=1,ncol=4)
    mean_est<-exp(MLE.i%*%x.mat)
    I.hat.i<- matrix(0, nrow = 4, ncol = 4)
    V.hat.i<- matrix(0, nrow = 4, ncol = 4)
    score<- matrix(0, nrow = seq, ncol = 4)
    for (k in 1:length(params) ){ 
      I.hat.i[k,]<-pi*mean_est%*%(t(x.mat)*x.mat[k,])
      score[,k]<-(sweep(data, 2, mean_est[1,]))%*%x.mat[k,]
      V.hat.i[k,]<- colSums(score*score[,k])/(seq*2)
     }
     
    V.hat.i<-Matrix::forceSymmetric(V.hat.i,uplo="L")
    
    # i.tt<-sum(mean_est)/2
    # i.ee<-(mean_est[2]+mean_est[3])/2
    # i.gg<-(mean_est[2]+mean_est[4])/2
    # i.dd<-(mean_est[3]+mean_est[4])/2
    # i.eg<-mean_est[2]/2
    # i.ed<-mean_est[3]/2
    # i.gd<-mean_est[4]/2
    # I.hat.i<-matrix(c(i.tt,i.ee,i.gg,i.dd,i.ee,i.ee,i.eg,i.ed,
    #              i.gg,i.eg,i.gg,i.gd,i.dd,i.ed,i.gd,i.dd), nrow = 4, ncol = 4,byrow=TRUE)
    # var<-diag(cov(data))
    # cov<-c(mean(data[,1]*data[,2])-mean_est[1]*mean_est[2],mean(data[,3]*data[,4])-mean_est[3]*mean_est[4])
    # v.tt = i.tt + sum(cov)
    # v.ee = i.ee
    # v.gg = i.gg
    # v.dd = i.dd + cov[2]
    # v.te = i.ee + sum(cov)/2
    # v.tg = i.gg + sum(cov)/2
    # v.td = i.dd + cov[2]
    # v.eg = i.eg + cov[2]/2
    # v.ed = i.ed + cov[2]/2
    # v.gd = i.gd + cov[2]/2
    # V.hat.i = matrix( c(v.tt, v.te, v.tg, v.td, v.te, v.ee, v.eg, v.ed, 
    #               v.tg, v.eg, v.gg, v.gd, v.td, v.ed, v.gd, v.dd),
    #             nrow=4, ncol=4, byrow = TRUE)
    MLE.closeform[i,]<- MLE.i
    I.hat<-I.hat+I.hat.i
    V.hat<-V.hat+V.hat.i
    invI.hat<-invI.hat+solve(I.hat.i)
    
    
  }
  print(paste('seq=',seq))
  print(invI.hat/sim_time)
  print(cov(MLE.closeform)*seq*2)
  
  
  #store result in MATS for seq
  MATS <- list(signif(t(as.matrix(colMeans(MLE.optim))),5),signif(t(as.matrix(colMeans(MLE.closeform))),5),signif(I.optim/(sim_time*seq*2),5),signif(I.hat/sim_time,5),signif(as.matrix(V.hat/sim_time),5),signif(cov(MLE.optim)*seq*2,5),signif(cov(MLE.closeform)*seq*2,5),signif(invI.hat/sim_time,5))
  #format result that we want to show in excel
   seq.result<-matrix(c('tao_hat','eta_hat','gamma_hat','delta_hat'),ncol = 4,nrow = 1)
   for (i in 1:length(obj)){
     temp<-rbind(obj.names[i,],MATS[[i]])
     seq.result<-rbind(seq.result,temp)}
   #store result of all seq_size
   num=seq_size[which(seq_size == seq)]
   result <- append(result, list( num = seq.result))
  
}
# MATS
library(xlsx)
write.xlsx(result, file = paste0(Sys.Date(),del.range[2]-del.range[1],'.xlsx'))


summary(MLE.closeform)
summary(MLE.optim)

