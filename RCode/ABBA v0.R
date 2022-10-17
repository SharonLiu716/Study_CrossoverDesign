set.seed(110225021)
setwd("C:/Github/Study_CrossoverDesign/RCode")
sim_time=2000
seq_size<-c(25,50,100,200)
cor_param=0.1
cros_type='ABBA'
#tao,eta,gamma,delta
params = matrix(c(1.0, 0.7, 0.3, 0.2), nrow=1,ncol=4)
#link of mean.i11,mean.mean.i12,mean.i21,mean.22
x.mat  <-  matrix(c(1,0,0,0,1,1,1,0,1,1,0,1,1,0,1,1), nrow = 4, ncol = 4)
mean_true=exp(params%*%x.mat)
pi=0.5#(seq_size)/(seq_size*2)
I.true<- matrix(0, nrow = 4, ncol = 4)
for (i in 1:length(params) ){ I.true[i,]<-0.5*mean_true%*%(t(x.mat)*x.mat[i,])}
#the form of matrix that export to excel
#output object
obj<-c('MLE.optim','MLE.closeform','I.hat.hessian','I.hat','V.hat.matix','NS.MLE.optim','NS.MLE.closeform','inv(I.hat)')
#name of output object
space<-matrix(c('','',''), nrow=1,ncol = 3)
obj.names<-matrix(0,ncol = 4,length(obj))
for (i in 1:length(obj)){obj.names[i,]<-cbind(obj[i],space)}


result <- list()#to store result of each seq_size
#use optim to obtain MLE、I_hat、V_hat
for (seq in seq_size){
  MLE.optim<-matrix(0, nrow = sim_time, ncol = 4)
  I.optim=0
  MLE.closeform<-matrix(0, nrow = sim_time, ncol = 4)
  I.hat<- 0
  V.hat<- 0
  for (i in 1:sim_time){
    pi=(seq)/(seq*2)
    data <-  matrix(0, nrow = seq, ncol = 4)
    for (j in 1:length(mean_true)) {data[,j] =rpois(seq, lambda = mean_true[1,j])}
    negll <- function(param) { sum(factorial(data))
      -sum( param[1]*data[,1]-exp(param[1])+sum(param[1:3])*data[,2]-exp(sum(param[1:3]))+(param[1]+param[2]+param[4])*data[,3]-exp(param[1]+param[2]+param[4])+(param[1]+param[3]+param[4])*data[,4]-exp(param[1]+param[3]+param[4]) )}
    ABBA<-optim(param <- c(0.95, 0.65, 0.25, 0.15), negll, hessian=TRUE)
    MLE.optim[i,]<-ABBA$par
    I.optim <-I.optim+ABBA$hessian
  
    y.sum=colSums(data)
    tao.hat=log(y.sum[1]/seq)
    eta.hat=0.5*(log(y.sum[2])+log(y.sum[3])-log(y.sum[1])-log(y.sum[4]))
    gamma.hat=0.5*(log(y.sum[2])+log(y.sum[4])-log(y.sum[1])-log(y.sum[3]))
    delta.hat=0.5*(log(y.sum[3])+log(y.sum[4])-log(y.sum[1])-log(y.sum[2]))
    MLE.i = matrix(c(tao.hat, eta.hat, gamma.hat, delta.hat), nrow=1,ncol=4)
    mean_est<-exp(MLE.i%*%x.mat)
    I.hat.i<- matrix(0, nrow = 4, ncol = 4)
    V.hat.i<- matrix(0, nrow = 4, ncol = 4)
    score<- matrix(0, nrow = seq, ncol = 4)
    for (k in 1:length(params) ){ 
       I.hat.i[k,]<-pi*mean_est%*%(t(x.mat)*x.mat[k,])
       score[,k]<-sweep(data, 2, mean_est[1,])%*%x.mat[k,]
       V.hat.i[k,]<- colSums(score*score[,k])/(seq*2)
    }
    
    V.hat.i<-Matrix::forceSymmetric(V.hat.i,uplo="L")
    MLE.closeform[i,]<- MLE.i
    I.hat<-I.hat+I.hat.i
    V.hat<-V.hat+V.hat.i
    
    
  }
  
  
  #remove outlier of MLE.closeform
  #MLE.closeform<-apply(MLE.closeform, 2, sort)
  #MLE.closeform <- MLE.closeform[11:(nrow(MLE.closeform) - 10),]
  
  #store result in MATS for seq
  MATS <- list(signif(t(as.matrix(colMeans(MLE.optim))),5),signif(t(as.matrix(colMeans(MLE.closeform))),5),signif(I.optim/(sim_time*seq*2),5),signif(I.hat/sim_time,5),signif(as.matrix(V.hat/sim_time),5),signif(cov(MLE.optim)*seq*2,5),signif(cov(MLE.closeform)*seq*2,5),signif(solve(I.hat/sim_time),5))
  #format result that we want to show in excel
  seq.result<-matrix(c('tao_hat','eta_hat','gamma_hat','delta_hat'),ncol = 4,nrow = 1)
  for (i in 1:length(obj)){
    temp<-rbind(obj.names[i,],MATS[[i]])
    seq.result<-rbind(seq.result,temp)}
  #store result of all seq_size
  num=seq_size[which(seq_size == seq)]
  result <- append(result, list( num = seq.result))
  
}

write.xlsx(result, file = paste0(Sys.Date(),'.xlsx'))


summary(MLE.closeform)
summary(MLE.optim)

