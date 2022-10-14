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

printmat <- function(name,mat) {
  out <- capture.output(mat)
  out[1] <- paste0(name, out[1])
  cat(paste(out, collapse = "\n"))
}



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
    V.hat.i[k,]<-colSums(score*score[,k])/(seq*2)
  }
  MLE.closeform[i,]<- MLE.i
  I.hat<-I.hat+I.hat.i
  V.hat<-V.hat+V.hat.i
  
}
  
  print(paste('seq_size',seq))
  printmat('MLE optim',colMeans(MLE.optim))
  print(' ')
  printmat('MLE closeform',colMeans(MLE.closeform))
  print(' ')
  printmat('I_optim',I.optim/(sim_time*seq*2))
  print(' ')
  printmat('I_hat',I.hat/sim_time)
  print(' ')
  printmat('V_hat',V.hat/sim_time)
  print(' ')
  printmat('N*S_optim',cov(MLE.optim)*seq*2)
  print(' ')
  printmat('N*S_closeform',cov(MLE.closeform)*seq*2)
  print(' ')
  
  printmat('inv(I_hat)',solve(I.hat/sim_time))
  print(' ')
  
  }

seq=50
for (i in 1:sim_time){
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
    
    
    I.hat.i[k,]<-mean_est%*%((t(x.mat)*x.mat[k,]))/2
    score[,k]<-sweep(data, 2, mean_est[1,])%*%x.mat[k,]
    V.hat.i[k,]<-colSums(score*score[,k])/(seq*2)
  }
  
  MLE.closeform[i,]<- MLE.i
  I.hat<-I.hat+I.hat.i
  V.hat<-V.hat+V.hat.i
}



colMeans(MLE.optim)
colMeans(MLE.closeform)
I.hat/sim_time
V.hat/sim_time
I.optim/(sim_time*seq*2)



#check hessian==I?

I
ABBA$hessian/200
solve(ABBA$hessian)*200
solve(I)


#start to do simulation(write as class)
#use close form to simulate
y.sum=colSums(data)
tao.hat=log(y.sum[1]/100)
eta.hat=0.5*(log(y.sum[2])+log(y.sum[3])-log(y.sum[1])-log(y.sum[4]))
gamma.hat=0.5*(log(y.sum[2])+log(y.sum[4])-log(y.sum[1])-log(y.sum[3]))
delta.hat=0.5*(log(y.sum[3])+log(y.sum[4])-log(y.sum[1])-log(y.sum[2]))
mle = matrix(c(tao.hat, eta.hat, gamma.hat, delta.hat), nrow=1,ncol=4)
mean_est<-exp(mle%*%x.mat)
I.hat<- matrix(0, nrow = 4, ncol = 4)
for (i in 1:nchar(cros_type) ){ I.hat[i,]<-pi*mean_est%*%(t(x.mat)*x.mat[i,])}
score<- matrix(0, nrow = 100, ncol = 4)
for (i in 1:length(params) ){ score[,i]<-sweep(data, 2, mean_est[1,])%*%x.mat[i,]} 
V.hat<- matrix(0, nrow = 4, ncol = 4)
for (i in 1:length(params) ){ V.hat[i,]<-colSums(score*score[,i])/200}

ABBA$hessian/200
solve(I.hat)

