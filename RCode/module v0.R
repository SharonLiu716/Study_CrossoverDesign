#simulation parameter
sim_time=2000
seq_size=c(25,50,100,200)
cros_type=c('ABBA','ABBBAA','AABABABAA','ABCBCACAB','BACACBBCA','BBAACBCAC')
#true value of params with treat-seq-time 
param_222=c(1.0,0.7,0.3,0.2);param_223=c(1.0,0.7,0.3,0.3,0.2);param_233=c(1.0,0.7,0.3,0.3,0.2,0.2);param_333=c(1.0,0.7,0.7,0.3,0.3,0.2,0.2)
#Link of Yist:Y11,Y12,Y21,Y22
#row is vector of param,column is x.mat of Yist 
xmat_222=matrix(c(1,0,0,0, 1,1,1,0, 1,1,0,1, 1,0,1,1), nrow = 4, ncol = 4)
xmat_223=matrix(c(1,0,0,0,0, 1,1,1,0,0, 1,1,0,1,0, 1,1,0,0,1, 1,0,1,0,1, 1,0,0,1,1), nrow = 5, ncol = 6)
xmat_233=matrix(c(1,0,0,0,0,0, 1,0,1,0,0,0, 1,1,0,1,0,0, 1,0,0,0,1,0, 1,1,1,0,0,0, 1,0,0,1,1,0, 1,1,0,0,0,1, 1,0,1,0,0,1, 1,0,0,1,0,1), nrow = 6, ncol = 9)
#ABCBCACAB
xmat_333.1=matrix(c(1,0,0,0,0,0,0, 1,1,0,1,0,0,0, 1,0,1,0,1,0,0, 1,1,0,0,0,1,0, 1,0,1,1,0,1,0, 1,0,0,0,1,1,0, 1,0,1,0,0,0,1, 1,0,0,1,0,0,1, 1,1,0,0,1,0,1), nrow = 7, ncol = 9)
#BACACBBCA
xmat_333.2=matrix(c(1,1,0,0,0,0,0, 1,0,0,1,0,0,0, 1,0,1,0,1,0,0, 1,0,0,0,0,1,0, 1,0,1,1,0,1,0, 1,1,0,0,1,1,0, 1,1,0,0,0,0,1, 1,0,1,1,0,0,1, 1,0,0,0,1,0,1), nrow = 7, ncol = 9)
#BBAACBCAC
xmat_333.3=matrix(c(1,1,0,0,0,0,0, 1,1,0,1,0,0,0, 1,0,0,0,1,0,0, 1,0,0,0,0,1,0, 1,0,1,1,0,1,0, 1,1,0,0,1,1,0, 1,0,1,0,0,0,1, 1,0,0,1,0,0,1, 1,0,1,0,1,0,1), nrow = 7, ncol = 9)



#Function of output object

#Par.values==TrueMean
Mean.True<-function(Par.values,x.mat){  return(exp(Par.values%*%x.mat))}

#I.true
Matrix.I<-function(cors.type,params,x.mat){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  mat.I<- matrix(0, nrow = length(params), ncol = length(params))
  for (i in 1:length(params) ){ mat.I[i,]<-Mean%*%(t(x.mat)*x.mat[i,])/num.seq}
  return(mat.I)
}

#Generate data
cros.type="ABBA"
mean.true=Mean.True(param_222,xmat_222)
seq.size=10
num.seq=2
cor.par=3
Data_Generate<-function(cros.type,mean.true,seq.size,data.type,cor.par,eta0){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  data <-  matrix(0, nrow = seq.size, ncol = nchar(cros.type))
  if (data.type == 'ind')
    for (j in 1:length(mean.true)) {data[,j] =rpois(seq.size, lambda = mean.true[1,j])}
  else{
    #生成一組nui乘上兩個mean待補mean[3:4]
    nui<-rgamma(n=seq.size,shape=1/cor.par,scale=cor.par)
    d<-split(c(outer(nui,mean.true[1:2],  function(x, y) x * y)), ceiling(seq_along(c(outer(nui,mean.true[1:2],  function(x, y) x * y)))/10))
    data[,1]<-unlist(d[1],use.names = FALSE)
    data[,2]<-unlist(d[2],use.names = FALSE)
  }
  
}
  






#MLE of different type of crossover design
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

#I.hat & V.hat
Matrix.IV<-function(cros.type,mle.values,x.mat,seq.size,data){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  mat.I<- matrix(0, nrow = length(mle.values), ncol = length(mle.values))
  mat.V<- matrix(0, nrow = length(mle.values), ncol = length(mle.values))
  mat.score<- matrix(0, nrow = seq.size, ncol = length(mle.values))
  mean.est<-exp(mle.values%*%x.mat)
  for (i in 1:length(mle.values) ){ 
    mat.I[i,]<-mean.est%*%(t(x.mat)*x.mat[i,])/num.seq
    mat.score[,i]<-sweep(data, 2, mean.est[1,])%*%x.mat[i,]
    mat.V[i,]<- colSums(mat.score*mat.score[,i])/(seq_size*num.seq)
  }
  mat.V<-Matrix::forceSymmetric(mat.V,uplo="L")
  return(mat.I,mat.V)
}

mean_est<-exp(MLE.i%*%x.mat)