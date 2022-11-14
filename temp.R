rm(list=ls(all=TRUE)) 
setwd("C:/Users/User/Documents/Study_CrossoverDesign/RCode")
sim_time=1000;cor_par=2
param=c(0.2,0,-0.7,-0.2)
xmat=matrix(c(1,1,1,1, 0,1,1,0, 0,1,0,1, 0,0,1,1), nrow = 4, ncol = 4,byrow = TRUE)
mean_true=exp(param%*%xmat)
Mean.True<-function(Par.values,x.mat){  return(exp(Par.values%*%x.mat))}
Data.ind<-function(cros.type,mean.true,seq.size){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  data <-  matrix(0, nrow = seq.size, ncol = nchar(cros.type))
  for (j in 1:length(mean.true)) {data[,j] =rpois(seq.size, lambda = mean.true[1,j])}
  return(data)
}
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
MLE.ABBA<-function(data){
  y.sum=colMeans(data)
  tao.hat=log(y.sum[1])
  eta.hat=0.5*(log(y.sum[2])+log(y.sum[3])-log(y.sum[1])-log(y.sum[4]))
  gamma.hat=0.5*(log(y.sum[2])+log(y.sum[4])-log(y.sum[1])-log(y.sum[3]))
  delta.hat=0.5*(log(y.sum[3])+log(y.sum[4])-log(y.sum[1])-log(y.sum[2]))
  MLE.i = matrix(c(tao.hat, eta.hat, gamma.hat, delta.hat), nrow=1,ncol=4)
  return(MLE.i)
}
MLE.ABBAnull<-function(data,seq.size,eta.null){
  y.sum=colMeans(data)
  tao.hatnull<-log(y.sum[1])
  gamma.hatnull=log(y.sum[2])-log(y.sum[1])-eta.null
  delta.hatnull=log(y.sum[4])-log(y.sum[2])+eta.null
  MLE.i.null = matrix(c(tao.hatnull,  eta.null, gamma.hatnull, delta.hatnull), nrow=1,ncol=4)
  return(MLE.i.null)
}
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
loglik.ABBA<-function(param,data){ 
  y.sum=colSums(data)
  n=nrow(data)
  ll= param[1]*y.sum[1]-n*exp(param[1])
      +sum(param[1:3])*y.sum[2]-n*exp(sum(param[1:3]))
      +(param[1]+param[2]+param[4])*y.sum[3]-n*exp(param[1]+param[2]+param[4])
      +(param[1]+param[3]+param[4])*y.sum[4]-n*exp(param[1]+param[3]+param[4])   
  return(ll)
}
set.seed(7353)
Matrix.AB<-function(Mat.I,Mat.V,loc){
  # loc=2
  # Mat.I<-I.cor/sim_time
  # Mat.V<-as.matrix(V.cor/sim_time)
  diag.I<-as.matrix(diag(Mat.I),byrow = TRUE)
  diag.V<-as.matrix(diag(Mat.V),byrow = TRUE)
  #I.theta.theta
  I.tt<-diag.I[loc];V.tt<-diag.V[loc];I.pp<-Mat.I[-c(loc),-c(loc)];V.pp<-Mat.V[-c(loc),-c(loc)]
  I.tp<-t(as.matrix(Mat.I[loc,]));I.tp<-I.tp[,-c(loc)];V.tp<-t(as.matrix(Mat.V[loc,]));V.tp<-V.tp[,-c(loc)]
  #mat.A
  Mat.A<-I.tt-I.tp%*%solve(I.pp)%*%as.matrix(I.tp)
  Mat.B<-V.tt - 2*I.tp %*% solve(I.pp) %*% V.tp + I.tp %*% solve(I.pp) %*% V.pp %*% solve(I.pp) %*% as.matrix(I.tp)
  Mat.AB<- list("Mat.A" = Mat.A, "Mat.B" =Mat.B)
  return(Mat.AB)
}
I.cor/sim_time+2*seq*cov(MLE.cor)
V.cor/sim_time
2*seq*cov(MLE.cor)
Mat.I[1,1]
cov(data.cor[,3],data.cor[,4])
(mean(data.cor[,3]*data.cor[,4])-mean(data.cor[,3])*mean(data.cor[,4]))
sum(diag(var(data.cor)))+2*1.0102041+2*-0.2346939
diag(var(data.ind))[1]+2*(cov(data.ind[,1],data.ind[,2]) +cov(data.ind[,3],data.ind[,4]))
colMeans(data.cor)
V.ind/sim_time
V.cor/sim_time
v.aa =  I[1,1]+ cov(data.cor[,1],data.cor[,2]) +cov(data.cor[,3],data.cor[,4])
v.ee = Mat.I[2,2]
v.gg = Mat.I[3,3]
v.dd = Mat.I[4,4] + cov(data.cor[,3],data.cor[,4])
  
v.ae = Mat.I[2,2] +  cov(data.cor[,1],data.cor[,2])/2 +cov(data.cor[,3],data.cor[,4])/2

v.ag = Mat.I[3,3]+  cov(data.cor[,1],data.cor[,2])/2 +cov(data.cor[,3],data.cor[,4])/2

v.ad = Mat.I[4,4] + cov(data.cor[,3],data.cor[,4])/2

v.eg = Mat.I[2,3] + cov(data.cor[3,],data.cor[4,])/2

v.ed = Mat.I[2,4]+ cov(data.cor[3,],data.cor[4,])/2
v.gd = Mat.I[3,4] + cov(data.cor[3,],data.cor[4,])/2
V = matrix( c(v.aa, v.ae, v.ag, v.ad, v.ae, v.ee, v.eg, v.ed, 
              v.ag, v.eg, v.gg, v.gd, v.ad, v.ed, v.gd, v.dd),
            nrow=4, ncol=4, byrow = TRUE)
I.tt - I.tp %*% solve(I.pp) %*% as.matrix(I.tp)
V.tt - 2 * I.tp %*% solve(I.pp) %*% V.tp + 
  I.tp %*% solve(I.pp) %*% V.pp %*% solve(I.pp) %*% as.matrix(I.tp)
set.seed(11018022)



#========================================
#main
#result:to store result of each seq_size
#========================================
result.ind <- list();result.cor <- list();pvalue.ind<-list();pvalue.cor<-list()
#simulation for ABBA
seq=50
MLE.ind<-matrix(0, nrow = sim_time, ncol = length(param))
MLE.cor<-matrix(0, nrow = sim_time, ncol = length(param))
MLE.null.ind<-matrix(0, nrow = sim_time, ncol = length(param))
MLE.null.cor<-matrix(0, nrow = sim_time, ncol = length(param))
I.ind<- 0 ; I.cor<- 0 ; V.ind<-0 ;V.cor<-0;invI.ind<-0;invI.cor<-0
# Create df for statistics
statistics.columns <- c("Wald.na","Wald.rb",'Score.na','Score.rb',"LR.na") 
df.statistics.ind <- data.frame(matrix(nrow = sim_time, ncol = length(statistics.columns))) 
colnames(df.statistics.ind) <- statistics.columns
df.statistics.cor = data.frame(matrix(nrow = sim_time, ncol = length(statistics.columns))) 
colnames(df.statistics.cor) <- statistics.columns
matA.ind<-c();matB.ind<-c();matA.cor<-c();matB.cor<-c()
indt.A<-c();indt.B<-c();indg.A<-c();indg.B<-c();indd.A<-c();indd.B<-c()
cort.A<-c();cort.B<-c();corg.A<-c();corg.B<-c();cord.A<-c();cord.B<-c()
set.seed(11018022)
for (i in 1:sim_time){
  #====================================================
  #I,V,IVI
  #====================================================
  theta.null=0
  X = c(rep(0,seq), rep(1,2*seq), rep(0,seq))
  Z = c(rep(0,seq), rep(1,seq), rep(0,seq), rep(1,seq))
  G = c(rep(0,2*seq), rep(1,2*seq))
  #---------------------------
  #independent
  #---------------------------
  data.ind<-Data.ind(cros.type = 'ABBA',mean.true = mean_true,seq.size = seq )
  #---------------------------
  #use package to do testing
  #---------------------------
  Y <- c(data.ind[,1],data.ind[,2],data.ind[,3],data.ind[,4])
  df.ind = data.frame(Y,X,Z,G)
  mod.1 <- glm(Y ~ X + Z + G, family = poisson(link = "log"), df.ind)
  mod.0 <- glm(Y ~ Z + G, family = poisson(link = "log"), df.ind)
  MLE.ind[i,]<-c( mod.1$coefficients[1],mod.1$coefficients[2],mod.1$coefficients[3],mod.1$coefficients[4])
  MLE.null.ind[i,]<-c( mod.0$coefficients[1],theta.null,mod.0$coefficients[2],mod.0$coefficients[3])
  #---------------------------
  #use closeform to do testing
  #---------------------------
  # MLE.ind[i,]<- MLE.ABBA(data.ind)
  # MLE.null.ind[i,]<-MLE.ABBAnull(data = data.ind, seq.size = seq, eta.null=theta.null)
  mean.est<-Mean.True(MLE.ind[i,],xmat)
  IV.ind.i<-Matrix.IV(cros.type='ABBA', mle.values=MLE.ind[i,], x.mat=xmat, seq.size=seq, data=data.ind)
  #---------------------------
  #correlated
  #---------------------------
  data.cor<-Data.cor(cros.type = 'ABBA',mean.true =mean_true,seq.size =seq,cor.par = cor_par )
  #use package to do testing
  Y <- c(data.cor[,1],data.cor[,2],data.cor[,3],data.cor[,4])
  df.cor = data.frame(Y,X,Z,G)
  mod.1 <- glm(Y ~ X + Z + G, family = poisson(link = "log"), df.cor)
  mod.0 <- glm(Y ~ Z + G, family = poisson(link = "log"), df.cor)
  MLE.cor[i,]<-c( mod.1$coefficients[1],mod.1$coefficients[2],mod.1$coefficients[3],mod.1$coefficients[4])
  MLE.null.cor[i,]<-c( mod.0$coefficients[1],theta.null,mod.0$coefficients[2],mod.0$coefficients[3])
  #---------------------------
  #use closeform to do testing
  #---------------------------
  # MLE.cor[i,]<- MLE.ABBA(data.cor)
  # MLE.null.cor[i,]<- MLE.ABBAnull(data = data.cor, seq.size = seq, eta.null=theta.null)
  mean.est<-Mean.True(MLE.cor[i,],xmat)
  IV.cor.i<-Matrix.IV(cros.type='ABBA', mle.values=MLE.cor[i,], x.mat=xmat, seq.size=seq, data=data.cor)
  
  #store result of MLE,I,V,inv.I
  
  I.ind.i<-IV.ind.i$I.hat
  V.ind.i<-Matrix::forceSymmetric(IV.ind.i$V.hat,uplo="L")
  invI.ind.i<-solve(I.ind.i)
  I.ind<-I.ind+I.ind.i
  V.ind<-V.ind+V.ind.i
  invI.ind<-invI.ind+invI.ind.i
  
  I.cor.i<-IV.cor.i$I.hat
  V.cor.i<-Matrix::forceSymmetric(IV.cor.i$V.hat,uplo="L")
  invI.cor.i<-solve(I.cor.i)
  I.cor<-I.cor+I.cor.i
  V.cor<-V.cor+V.cor.i
  invI.cor<-invI.cor+invI.cor.i
  #=========================================================================
  #output:matrix AB, Wald statistics,LR statistics, Score statistics
  #=========================================================================
  matA.ind[i]<-Matrix.AB(I.ind.i,as.matrix(V.ind.i),2)$Mat.A
  matB.ind[i]<-Matrix.AB(I.ind.i,as.matrix(V.ind.i),2)$Mat.B
  indt.A[i]<-Matrix.AB(I.ind.i,as.matrix(V.ind.i),1)$Mat.A
  indt.B[i]<-Matrix.AB(I.ind.i,as.matrix(V.ind.i),1)$Mat.B
  indg.A[i]<-Matrix.AB(I.ind.i,as.matrix(V.ind.i),3)$Mat.A
  indg.B[i]<-Matrix.AB(I.ind.i,as.matrix(V.ind.i),3)$Mat.B
  indd.A[i]<-Matrix.AB(I.ind.i,as.matrix(V.ind.i),4)$Mat.A
  indd.B[i]<-Matrix.AB(I.ind.i,as.matrix(V.ind.i),4)$Mat.B
  
  
  df.statistics.ind[i, "Wald.na"] <- seq*2*matA.ind[i]*(MLE.ind[i,2]-theta.null)^2
  df.statistics.ind[i, "Wald.rb"] <- seq*2*(matA.ind[i]^2)/matB.ind[i]*((MLE.ind[i,2]-theta.null)^2)
  df.statistics.ind[i, "LR.na"] <- 2*(loglik.ABBA(param = MLE.ind[i,],data = data.ind)-loglik.ABBA(param = as.vector(MLE.null.ind[i,]),data = data.ind))
  # df.statistics.ind[i, "LR.rb"] <- 2*matA.ind.i[i]/matB.ind.i[i]*(loglik.ABBA(param = MLE.ind.i,data = data.ind)-loglik.ABBA(param = as.vector(MLE.null),data = data.ind))
  mean.null<-Mean.True(MLE.null.ind[i,],xmat)
  df.statistics.ind[i, "Score.na"] <- (sum(data.ind[,2]-mean.null[,2])+sum(data.ind[,3]-mean.null[,3]))^2 / (matA.ind[i]*seq*2)
  df.statistics.ind[i, "Score.rb"] <- (sum(data.ind[,2]-mean.null[,2])+sum(data.ind[,3]-mean.null[,3]))^2 / (matB.ind[i]*seq*2)
  
  
  # #wald.test(Sigma = cov(MLE.ind), b = t(MLE.ind.i), Terms = 2)
  
  matA.cor[i]<-Matrix.AB(I.cor.i,as.matrix(V.cor.i),2)$Mat.A
  matB.cor[i]<-Matrix.AB(I.cor.i,as.matrix(V.cor.i),2)$Mat.B
  cort.A[i]<-Matrix.AB(I.cor.i,as.matrix(V.cor.i),1)$Mat.A
  cort.B[i]<-Matrix.AB(I.cor.i,as.matrix(V.cor.i),1)$Mat.B
  corg.A[i]<-Matrix.AB(I.cor.i,as.matrix(V.cor.i),3)$Mat.A
  corg.B[i]<-Matrix.AB(I.cor.i,as.matrix(V.cor.i),3)$Mat.B
  cord.A[i]<-Matrix.AB(I.cor.i,as.matrix(V.cor.i),4)$Mat.A
  cord.B[i]<-Matrix.AB(I.cor.i,as.matrix(V.cor.i),4)$Mat.B
  
  
  
  
  df.statistics.cor[i, "Wald.na"] <- seq*2*matA.cor[i]*(MLE.cor[i,2]-theta.null)^2
  df.statistics.cor[i, "Wald.rb"] <- seq*2*(matA.cor[i]^2)/matB.cor[i]*((MLE.cor[i,2]-theta.null)^2)
  df.statistics.cor[i, "LR.na"] <- 2*(loglik.ABBA(param = MLE.cor[i,],data = data.cor)-loglik.ABBA(param = as.vector(MLE.null.cor[i,]),data = data.cor))
  # df.statistics.cor[i, "LR.rb"] <- 2*matA.cor[i]/matB.cor[i]*(loglik.ABBA(param = MLE.cor.i,data = data.cor)-loglik.ABBA(param = as.vector(MLE.null),data = data.cor))
  mean.null<-Mean.True(MLE.null.cor[i,],xmat)
  df.statistics.cor[i, "Score.na"] <- (sum(data.cor[,2]-mean.null[,2])+sum(data.cor[,3]-mean.null[,3]))^2 / (matA.cor[i]*seq*2)
  df.statistics.cor[i, "Score.rb"] <- (sum(data.cor[,2]-mean.null[,2])+sum(data.cor[,3]-mean.null[,3]))^2 / (matB.cor[i]*seq*2)
}


#=========================================
#output for each seq.size:I, V, IVI
#=========================================

#store result in MATS for seq
MATS.ind <- list(signif(t(as.matrix(colMeans(MLE.ind))),5),signif(diag(cov(MLE.ind)),5),signif(mean(1/(2*seq*matA.ind)) ,5),signif(mean(matB.ind/(2*seq*matA.ind^2)) ,5))
MATS.cor <- list(signif(t(as.matrix(colMeans(MLE.cor))),5),signif(diag(cov(MLE.cor)),5),signif(mean(1/(2*seq*matA.cor)) ,5),signif(mean(matB.cor/(2*seq*matA.cor^2)) ,5))
#testing
below <- function(x) {return(length(x[x>qchisq(0.95,1)])/sim_time)}
print(paste('ind',seq,sapply(df.statistics.ind, below)))
print(paste('cor',seq,sapply(df.statistics.cor, below)))

mean(matB.cor)
sum(MLE.ind[1,])
sum(MLE.null.ind[1,])
2*(loglik.ABBA(param = MLE.ind[702,],data = data.ind)-loglik.ABBA(param = as.vector(MLE.null.ind[702,]),data = data.ind))

mean(1/(2*seq*indt.A))
mean(indt.B/(2*seq*indt.A^2))
mean(1/(2*seq*indg.A))
mean(indg.B/(2*seq*indg.A^2))
mean(1/(2*seq*indd.A))
mean(indd.B/(2*seq*indd.A^2))

mean(1/(2*seq*cort.A))
mean(cort.B/(2*seq*cort.A^2))
mean(1/(2*seq*corg.A))
mean(corg.B/(2*seq*corg.A^2))
mean(1/(2*seq*cord.A))
mean(cord.B/(2*seq*cord.A^2))


