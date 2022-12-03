rm(list=ls(all=TRUE)) 
require(faraway)
require(numDeriv)
require(MASS)

setwd("C:/Users/User/Documents/Study_CrossoverDesign/RCode")
sim_time=10000;cor_par=1;seq=50
param=c(1.2,0,1.0,0.2)#c(1.2,0,1.0,0.2)#c(0.3,0,0.4,-1.0)
xmat=matrix(c(1,1,1,1, 0,1,1,0, 0,1,0,1, 0,0,1,1), nrow = 4, ncol = 4,byrow = TRUE)
mean.true=exp(param%*%xmat)
X = c(rep(0,seq), rep(1,2*seq), rep(0,seq))
Z = c(rep(0,seq), rep(1,seq), rep(0,seq), rep(1,seq))
G = c(rep(0,2*seq), rep(1,2*seq))

MLE.optim<-matrix(0, nrow = sim_time, ncol = 4)
MLE.glm<-matrix(0, nrow = sim_time, ncol = 4)
MLE.closeform<-matrix(0, nrow = sim_time, ncol = 4)
invI.optim=0;I.glm=0;I.optim=0;invI.glm=0;
I.cf=0;V.cf=0;invI.closeform=0;I0.cf=0;V0.cf=0;invI0.closeform=0
tao<-c();eta<-c();gam<-c();del<-c()
tao.0<-c();eta.0<-c();gam.0<-c();del.0<-c()
W.na1<-0;W.rb1<-0;LR.na1<-0;LR.rb1<-0;S.na1<-0;S.rb1<-0
eta.var<-c()
#========================================================
#independent data
#========================================================

set.seed(110225021)
for (i in 1:sim_time){
  #mean=r*(1-p)/p
  # y11<-rnbinom(n = seq, size = 2, prob = 0.6)
  # y12<-rnbinom(n = seq, size = 3, prob = 0.6)
  # y21<-rnbinom(n = seq, size = 2, prob = 0.8)
  # y22<-rnbinom(n = 100000, size = 3, prob = 0.8)
  # Y <- c(y11,y12,y21,y22)
  
  y11<-rpois(seq, lambda = mean.true[1])
  y12<-rpois(seq, lambda = mean.true[2])
  y21<-rpois(seq, lambda = mean.true[3])
  y22<-rpois(seq, lambda = mean.true[4])
  
  
 
  # Y <- c(y11,y12,y21,y22)
  # df.ind = data.frame(Y,X,Z,G)
  # #---------------------------------
  # #loglikelihood own define
  # #---------------------------------
  # negll <- function(param) { 
  #   -sum( param[1]*y11-exp(param[1])+sum(param[1:3])*y12-exp(sum(param[1:3]))+(param[1]+param[2]+param[4])*y21-exp(param[1]+param[2]+param[4])+(param[1]+param[3]+param[4])*y22-exp(param[1]+param[3]+param[4]) )}
  # ABBA<-optim(param <- c(0.45, -0.05, 0.35, -0.25), negll, hessian=TRUE)
  # MLE.optim[i,]<-ABBA$par
  # I.optim<-I.optim+(ABBA$hessian)/(2*seq)
  # invI.optim <-invI.optim+solve((ABBA$hessian)/(2*seq))
  # #---------------------------------
  # #GLM
  # #---------------------------------
  # mod.1 <- glm(Y ~ X + Z + G, family = poisson(link = "log"), df.ind)
  # MLE.glm[i,]<-coef(mod.1)
  # I.glm<-I.glm+solve(vcov(mod.1))/(2*seq)
  # invI.glm <-invI.glm+vcov(mod.1)*(2*seq)
  #---------------------------------
  #closeformL:exp(.)
  #---------------------------------
  
  tao[i]<-mean(y11)
  eta[i]<-sqrt( sum(y12)*sum(y21)/(sum(y11)*sum(y22)) )
  gam[i]<-sqrt( sum(y12)*sum(y22)/(sum(y11)*sum(y21)) )
  del[i]<-sqrt( sum(y21)*sum(y22)/(sum(y11)*sum(y12)) )
  
  # cov1<-mean( (y11-tao[i])*(y12-tao[i]*eta[i]*gam[i]) )
  # cov2<-mean( (y21-tao[i]*eta[i]*del[i])*(y22-tao[i]*gam[i]*del[i]) )
  #Matrix I & Matrix V
  
  i.tt<-( tao[i]+tao[i]*eta[i]*gam[i] + tao[i]*eta[i]*del[i]+tao[i]*gam[i]*del[i] )/2
  i.gg<-( tao[i]*eta[i]*gam[i] + tao[i]*gam[i]*del[i] )/2 
    
  i.ee<-( tao[i]*eta[i]*gam[i] + tao[i]*eta[i]*del[i] )/2 
  
  i.dd<- ( tao[i]*eta[i]*del[i]+tao[i]*gam[i]*del[i] )/2
  i.eg<- ( tao[i]*eta[i]*gam[i] )/2
  i.ed<- ( tao[i]*eta[i]*del[i] )/2
  i.gd<- ( tao[i]*gam[i]*del[i] )/2
  I<-matrix(c(i.tt,i.ee,i.gg,i.dd,
              i.ee,i.ee,i.eg,i.ed,
              i.gg,i.eg,i.gg,i.gd,
              i.dd,i.ed,i.gd,i.dd),nrow=4, ncol=4, byrow = TRUE)
  I.cf<-I.cf+I
  
  cov1<-mean( (y11-tao[i])*(y12-tao[i]*eta[i]*gam[i]) )
  cov2<-mean( (y21-tao[i]*eta[i]*del[i])*(y22-tao[i]*gam[i]*del[i]) )
  
  invI.closeform <- invI.closeform + solve(I)
  v.tt = i.tt + cov1+cov2#cov(y11,y12)+cov(y21,y22)
  v.ee = i.ee
  v.gg = i.gg 
  v.dd = i.dd + cov2#cov(y21,y22)
  
  v.te = i.ee + cov1/2+cov2/2#cov(y11,y12)/2+cov(y21,y22)/2
  v.tg = i.gg + cov1/2+cov2/2#cov(y11,y12)/2+cov(y21,y22)/2
  v.td = i.dd + cov1#cov(y21,y22)
    
  
  v.eg = i.eg + cov2/2#cov(y21,y22)/2
  v.ed = i.ed + cov2/2#cov(y21,y22)/2
  v.gd = i.gd + cov2/2#cov(y21,y22)/2

  V.cf<-V.cf+matrix( c(v.tt, v.te, v.tg, v.td,
                       v.te, v.ee, v.eg, v.ed, 
                       v.tg, v.eg, v.gg, v.gd,
                       v.td, v.ed, v.gd, v.td),nrow=4, ncol=4, byrow = TRUE)
  
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.tt, i.gg, i.dd, i.gg, i.gg, i.gd, 
                  i.dd, i.gd, i.dd), nrow = 3, ncol = 3)

  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  v.pe = matrix(c(v.te, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.tt, v.tg, v.td, v.tg, v.gg, v.gd, v.td, 
                  v.gd, v.td), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  
  eta.var[i]<-0.25*(1/sum(y11)+1/sum(y12)+1/sum(y21)+1/sum(y22))
  
  #wald test
  wna = log(eta[i]) * A * log(eta[i]) * 2*seq
  wrb = log(eta[i]) * A^2 / B * log(eta[i]) * 2*seq
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  
  #log LR test :par=exp(.)
  lik<-function(par){ 
    ll=sum(log(par[1])*y11-par[1]+log(par[1]*par[2]*par[3])*y12-par[1]*par[2]*par[3])+
       sum(log(par[1]*par[2]*par[4])*y21-par[1]*par[2]*par[4]+ log(par[1]*par[3]*par[4])*y22-par[1]*par[3]*par[4])
    
    # ll= par[1]*sum(y11)-seq*exp(param[1])
    #   +sum(par[1:3])*sum(y12)-seq*exp(sum(par[1:3]))
    #   +(par[1]+par[2]+par[4])*sum(y21)-seq*exp(par[1]+par[2]+param[4])
    #   +(par[1]+par[3]+par[4])*sum(y22)-seq*exp(param[1]+param[3]+param[4])   
    return(ll)
  }
  #null MLE
  Y <- c(y11,y12,y21,y22)
  df.ind = data.frame(Y,X,Z,G)
  mod.0 <- glm(Y ~ Z + G, family = poisson(link = "log"), df.ind)
  tao.0[i]<-exp(mod.0$coefficients[1])
  gam.0[i]<-exp(mod.0$coefficients[2])
  del.0[i]<-exp(mod.0$coefficients[3])
  #null cov
  cov1<-mean( (y11-tao.0[i])*(y12-tao.0[i]*gam.0[i]) )
  cov2<-mean( (y21-tao.0[i]*del.0[i])*(y22-tao.0[i]*gam.0[i]*del.0[i]) )
  
  #LR Test
  l1 = lik(c(tao[i], eta[i], gam[i], del[i]) )
  l0 = lik(c(tao.0[i], 1, gam.0[i], del.0[i]))
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2*A/B*(l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  
  #score test
  i0.tt<-( tao.0[i]+tao.0[i]*gam.0[i] + tao.0[i]*del[i]+tao[i]*gam[i]*del[i] )/2
  i0.gg<-( tao.0[i]*gam.0[i] + tao.0[i]*gam.0[i]*del.0[i] )/2 
  
  i0.ee<-( tao.0[i]*gam.0[i] + tao.0[i]*del.0[i] )/2 
  
  i0.dd<- ( tao.0[i]*del.0[i]+tao.0[i]*gam.0[i]*del.0[i] )/2
  i0.eg<- ( tao.0[i]*gam.0[i] )/2
  i0.ed<- ( tao.0[i]*del.0[i] )/2
  i0.gd<- ( tao.0[i]*gam.0[i]*del.0[i] )/2
  I0<-matrix(c(i0.tt,i0.ee,i0.gg,i0.dd,
               i0.ee,i0.ee,i0.eg,i0.ed,
               i0.gg,i0.eg,i0.gg,i0.gd,
               i0.dd,i0.ed,i0.gd,i0.dd),nrow=4, ncol=4, byrow = TRUE)
  I0.cf<-I0.cf+I0
  
  invI0.closeform <- invI0.closeform + solve(I0)
  v0.tt = i0.tt + cov1+cov2
  v0.ee = i0.ee
  v0.gg = i0.gg 
  v0.dd.0 = i0.dd + cov2
  
  v0.te = i0.ee + cov1/2+cov2/2
  v0.tg = i0.gg + cov1/2+cov2/2
  v0.td = i0.dd + cov2
  
  
  v0.eg = i0.eg + cov2/2
  v0.ed = i0.ed + cov2/2
  v0.gd = i0.gd + cov2/2
  
  V0.cf<-V.cf+matrix( c(v0.tt, v0.te, v0.tg, v0.td,
                       v0.te, v0.ee, v.eg, v0.ed, 
                       v0.tg, v0.eg, v0.gg, v0.gd,
                       v0.td, v0.ed, v0.gd, v0.td),nrow=4, ncol=4, byrow = TRUE)
  
  i0.ep = matrix(c(i0.ee, i0.eg, i0.ed), ncol = 3)
  i0.pp = matrix(c(i0.tt, i0.gg, i0.dd, i0.gg, i0.gg, i0.gd, 
                  i0.dd, i0.gd, i0.dd), nrow = 3, ncol = 3)
  
  A0 = i0.ee - i0.ep %*% solve(i0.pp) %*% t(i0.ep)
  v0.pe = matrix(c(v0.te, v0.eg, v0.ed), nrow = 3)
  v0.pp = matrix(c(v0.tt, v0.tg, v0.td, v0.tg, v0.gg, v0.gd, v0.td, 
                  v0.gd, v0.td), nrow = 3, ncol = 3)
  B0 = v0.ee - 2 * i0.ep %*% solve(i0.pp) %*% v0.pe + 
    i0.ep %*% solve(i0.pp) %*% v0.pp %*% solve(i0.pp) %*% t(i0.ep)
  s0 = seq * ( mean(y12) - tao.0[i]*gam.0[i] + mean(y21)- tao.0[i]*del.0[i] )
  sna = s0 / A0 * s0 / (2*seq)
  srb = s0 / B0 * s0 / (2*seq)
  if( sna <= qchisq(0.95, 1) )  S.na1 = S.na1+1
  if( srb <= qchisq(0.95, 1) )  S.rb1 = S.rb1+1
  
}

colMeans(MLE.optim)
I.optim/sim_time
invI.optim/sim_time

colMeans(MLE.glm)
I.glm/sim_time
invI.glm/sim_time

log(c(mean(tao),mean(eta),mean(gam),mean(del)))
cov.m<-2*seq*matrix(c(var(log(tao)),cov(log(tao),log(eta)),cov(log(tao),log(gam)),cov(log(tao),log(del)), 
                cov(log(tao),log(eta)),var(log(eta)),cov(log(gam),log(eta)),cov(log(del),log(eta)), 
                cov(log(tao),log(gam)),cov(log(gam),log(eta)),var(log(gam)),cov(log(gam),log(del)), 
                cov(log(tao),log(del)),cov(log(del),log(eta)),cov(log(gam),log(del)),var(log(del))), nrow = 4, ncol = 4,byrow = TRUE)

I.cf/sim_time
V.cf/sim_time
(invI.closeform/sim_time)%*%(V.cf/sim_time)%*%(invI.closeform/sim_time)
mean(eta.var)*2*seq
var(log(eta))*2*seq
1/A
A/B/B
1/Matrix.AB(I.cf/sim_time,V.cf/sim_time,2)$Mat.A
1/Matrix.AB(I.cf/sim_time,V.cf/sim_time,2)$Mat.B
wna/sim_time
wrb/sim_time
wna/sim_time
wrb/sim_time
#========================================================
#correlated data
#========================================================
MLE.optim<-matrix(0, nrow = sim_time, ncol = 4)
MLE.glm<-matrix(0, nrow = sim_time, ncol = 4)
MLE.closeform<-matrix(0, nrow = sim_time, ncol = 4)
invI.optim=0;invI.glm=0;invI.closeform=0;I.cf=0;V.cf=0;I.glm=0;I.optim=0
tao<-c();eta<-c();gam<-c();del<-c()
mean.cor<-matrix(0, nrow = seq, ncol = nchar('ABBA'))
data.cor<-matrix(0, nrow = seq, ncol = nchar('ABBA'))
set.seed(110225021)

for (i in 1:sim_time){
   nui<-replicate(2,rgamma(n=seq,shape=1/cor_par,scale=cor_par))
   mean.cor[,1]<-mean.true[,1]*nui[,1]
   mean.cor[,2]<-mean.true[,2]*nui[,1]
   mean.cor[,3]<-mean.true[,3]*nui[,2]
   mean.cor[,4]<-mean.true[,4]*nui[,2]
   
   for (i in 1:nchar('ABBA')){
    list_poisson <- unlist(lapply(mean.cor[,i], FUN = function(x) rpois(1, x)))
    data.cor[,i] <- list_poisson
  }
   y11<-data.cor[,1];y12<-data.cor[,2];y21<-data.cor[,3];y22<-data.cor[,4]
   
  # y1 <- rbvpois(100000, mean.true[1],  mean.true[2], 2.0)
  # y2 <- rbvpois(100000, mean.true[3],  mean.true[4], 1.5)
  # y11<-y1[,1];y12<-y1[,2];y21<-y2[,1];y22<-y2[,2]
  
  # l1 <- c(1.33, 2); l2 <- c(0.5, 0.75) # lambda for each new variable
  # y1 <- genCorGen(seq, nvars = 2, params1 = l1, dist = "poisson", rho = .4, corstr = "cs", wide = TRUE,cnames='y11,y12')
  # y1 <-as.matrix(y1[,c('y11','y12')])
  # y2 <- genCorGen(seq, nvars = 2, params1 = l2, dist = "poisson", rho = .4, corstr = "cs", wide = TRUE,cnames='y21,y22')
  # y2 <-as.matrix(y2[,c('y21','y22')])
  # y11<-y1[,1];y12<-y1[,2];y21<-y2[,1];y22<-y2[,2]
  # Y <- c(y11,y12,y21,y22)
  
  #---------------------------------
  #closeform MLE:exp(.)
  #---------------------------------
  
  tao[i]<-mean(y11)
  eta[i]<-sqrt( sum(y12)*sum(y21)/(sum(y11)*sum(y22)) )
  gam[i]<-sqrt( sum(y12)*sum(y22)/(sum(y11)*sum(y21)) )
  del[i]<-sqrt( sum(y21)*sum(y22)/(sum(y11)*sum(y12)) )
  
  # cov1<-mean( (y11-tao[i])*(y12-tao[i]*eta[i]*gam[i]) )
  # cov2<-mean( (y21-tao[i]*eta[i]*del[i])*(y22-tao[i]*gam[i]*del[i]) )
  #Matrix I & Matrix V
  
  i.tt<-( tao[i]+tao[i]*eta[i]*gam[i] + tao[i]*eta[i]*del[i]+tao[i]*gam[i]*del[i] )/2
  i.gg<-( tao[i]*eta[i]*gam[i] + tao[i]*gam[i]*del[i] )/2 
  
  i.ee<-( tao[i]*eta[i]*gam[i] + tao[i]*eta[i]*del[i] )/2 
  
  i.dd<- ( tao[i]*eta[i]*del[i]+tao[i]*gam[i]*del[i] )/2
  i.eg<- ( tao[i]*eta[i]*gam[i] )/2
  i.ed<- ( tao[i]*eta[i]*del[i] )/2
  i.gd<- ( tao[i]*gam[i]*del[i] )/2
  I<-matrix(c(i.tt,i.ee,i.gg,i.dd,
              i.ee,i.ee,i.eg,i.ed,
              i.gg,i.eg,i.gg,i.gd,
              i.dd,i.ed,i.gd,i.dd),nrow=4, ncol=4, byrow = TRUE)
  I.cf<-I.cf+I
  
  cov1<-mean( (y11-tao[i])*(y12-tao[i]*eta[i]*gam[i]) )
  cov2<-mean( (y21-tao[i]*eta[i]*del[i])*(y22-tao[i]*gam[i]*del[i]) )
  
  invI.closeform <- invI.closeform + solve(I)
  v.tt = i.tt + cov1+cov2#cov(y11,y12)+cov(y21,y22)
  v.ee = i.ee
  v.gg = i.gg 
  v.dd = i.dd + cov2#cov(y21,y22)
  
  v.te = i.ee + cov1/2+cov2/2#cov(y11,y12)/2+cov(y21,y22)/2
  v.tg = i.gg + cov1/2+cov2/2#cov(y11,y12)/2+cov(y21,y22)/2
  v.td = i.dd + cov1#cov(y21,y22)
  
  
  v.eg = i.eg + cov2/2#cov(y21,y22)/2
  v.ed = i.ed + cov2/2#cov(y21,y22)/2
  v.gd = i.gd + cov2/2#cov(y21,y22)/2
  
  V.cf<-V.cf+matrix( c(v.tt, v.te, v.tg, v.td,
                       v.te, v.ee, v.eg, v.ed, 
                       v.tg, v.eg, v.gg, v.gd,
                       v.td, v.ed, v.gd, v.td),nrow=4, ncol=4, byrow = TRUE)
  
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.tt, i.gg, i.dd, i.gg, i.gg, i.gd, 
                  i.dd, i.gd, i.dd), nrow = 3, ncol = 3)
  
  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  v.pe = matrix(c(v.te, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.tt, v.tg, v.td, v.tg, v.gg, v.gd, v.td, 
                  v.gd, v.td), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  
  
  #wald test
  wna = log(eta[i]) * A * log(eta[i]) * 2*seq
  wrb = log(eta[i]) * A^2 / B * log(eta[i]) * 2*seq
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  
  #log LR test :par=exp(.)
  lik<-function(par){ 
    ll=sum(log(par[1])*y11-par[1]+log(par[1]*par[2]*par[3])*y12-par[1]*par[2]*par[3])+
      sum(log(par[1]*par[2]*par[4])*y21-par[1]*par[2]*par[4]+ log(par[1]*par[3]*par[4])*y22-par[1]*par[3]*par[4])
    
    # ll= par[1]*sum(y11)-seq*exp(param[1])
    #   +sum(par[1:3])*sum(y12)-seq*exp(sum(par[1:3]))
    #   +(par[1]+par[2]+par[4])*sum(y21)-seq*exp(par[1]+par[2]+param[4])
    #   +(par[1]+par[3]+par[4])*sum(y22)-seq*exp(param[1]+param[3]+param[4])   
    return(ll)
  }
  #null MLE
  Y<-c(y11,y12,y21,y22)
  df.cor = data.frame(Y,X,Z,G)
  mod.0 <- glm(Y ~ Z + G, family = poisson(link = "log"), df.cor)
  tao.0[i]<-exp(mod.0$coefficients[1])
  gam.0[i]<-exp(mod.0$coefficients[2])
  del.0[i]<-exp(mod.0$coefficients[3])
  #null cov
  cov1<-mean( (y11-tao.0[i])*(y12-tao.0[i]*gam.0[i]) )
  cov2<-mean( (y21-tao.0[i]*del.0[i])*(y22-tao.0[i]*gam.0[i]*del.0[i]) )
  
  #LR Test
  l1 = lik(c(tao[i], eta[i], gam[i], del[i]) )
  l0 = lik(c(tao.0[i], 1, gam.0[i], del.0[i]))
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2*A/B*(l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  
  #score test
  i0.tt<-( tao.0[i]+tao.0[i]*gam.0[i] + tao.0[i]*del[i]+tao[i]*gam[i]*del[i] )/2
  i0.gg<-( tao.0[i]*gam.0[i] + tao.0[i]*gam.0[i]*del.0[i] )/2 
  
  i0.ee<-( tao.0[i]*gam.0[i] + tao.0[i]*del.0[i] )/2 
  
  i0.dd<- ( tao.0[i]*del.0[i]+tao.0[i]*gam.0[i]*del.0[i] )/2
  i0.eg<- ( tao.0[i]*gam.0[i] )/2
  i0.ed<- ( tao.0[i]*del.0[i] )/2
  i0.gd<- ( tao.0[i]*gam.0[i]*del.0[i] )/2
  I0<-matrix(c(i0.tt,i0.ee,i0.gg,i0.dd,
               i0.ee,i0.ee,i0.eg,i0.ed,
               i0.gg,i0.eg,i0.gg,i0.gd,
               i0.dd,i0.ed,i0.gd,i0.dd),nrow=4, ncol=4, byrow = TRUE)
  I0.cf<-I0.cf+I0
  
  invI0.closeform <- invI0.closeform + solve(I0)
  v0.tt = i0.tt + cov1+cov2
  v0.ee = i0.ee
  v0.gg = i0.gg 
  v0.dd.0 = i0.dd + cov2
  
  v0.te = i0.ee + cov1/2+cov2/2
  v0.tg = i0.gg + cov1/2+cov2/2
  v0.td = i0.dd + cov2
  
  
  v0.eg = i0.eg + cov2/2
  v0.ed = i0.ed + cov2/2
  v0.gd = i0.gd + cov2/2
  
  V0.cf<-V.cf+matrix( c(v0.tt, v0.te, v0.tg, v0.td,
                        v0.te, v0.ee, v.eg, v0.ed, 
                        v0.tg, v0.eg, v0.gg, v0.gd,
                        v0.td, v0.ed, v0.gd, v0.td),nrow=4, ncol=4, byrow = TRUE)
  
  i0.ep = matrix(c(i0.ee, i0.eg, i0.ed), ncol = 3)
  i0.pp = matrix(c(i0.tt, i0.gg, i0.dd, i0.gg, i0.gg, i0.gd, 
                   i0.dd, i0.gd, i0.dd), nrow = 3, ncol = 3)
  
  A0 = i0.ee - i0.ep %*% solve(i0.pp) %*% t(i0.ep)
  v0.pe = matrix(c(v0.te, v0.eg, v0.ed), nrow = 3)
  v0.pp = matrix(c(v0.tt, v0.tg, v0.td, v0.tg, v0.gg, v0.gd, v0.td, 
                   v0.gd, v0.td), nrow = 3, ncol = 3)
  B0 = v0.ee - 2 * i0.ep %*% solve(i0.pp) %*% v0.pe + 
    i0.ep %*% solve(i0.pp) %*% v0.pp %*% solve(i0.pp) %*% t(i0.ep)
  s0 = seq * ( mean(y12) - tao.0[i]*gam.0[i] + mean(y21)- tao.0[i]*del.0[i] )
  sna = s0 / A0 * s0 / (2*seq)
  srb = s0 / B0 * s0 / (2*seq)
  if( sna <= qchisq(0.95, 1) )  S.na1 = S.na1+1
  if( srb <= qchisq(0.95, 1) )  S.rb1 = S.rb1+1
  
}


for (i in 1:sim_time){
   
  #  nui<-replicate(2,rgamma(n=seq,shape=1/cor_par,scale=cor_par))
  #  mean.cor[,1]<-mean.true[,1]*nui[,1]
  #  mean.cor[,2]<-mean.true[,2]*nui[,1]
  #  mean.cor[,3]<-mean.true[,3]*nui[,2]
  #  mean.cor[,4]<-mean.true[,4]*nui[,2]
  # 
  #  for (i in 1:nchar('ABBA')){
  #    list_poisson <- unlist(lapply(mean.cor[,i], FUN = function(x) rpois(1, x)))
  #    data.cor[,i]<-list_poisson
  #  }
  # y11<-data.cor[,1];y12<-data.cor[,2];y21<-data.cor[,3];y22<-data.cor[,4]
  # Y<-c(y11,y12,y21,y22)
  
  # y1 <- rbvpois(100000, mean.true[1],  mean.true[2], 2.0)
  # y2 <- rbvpois(100000, mean.true[3],  mean.true[4], 1.5)
  # y11<-y1[,1];y12<-y1[,2];y21<-y2[,1];y22<-y2[,2]
  
  # l1 <- c(1.33, 2); l2 <- c(0.5, 0.75) # lambda for each new variable
  # y1 <- genCorGen(seq, nvars = 2, params1 = l1, dist = "poisson", rho = .4, corstr = "cs", wide = TRUE,cnames='y11,y12')
  # y1 <-as.matrix(y1[,c('y11','y12')])
  # y2 <- genCorGen(seq, nvars = 2, params1 = l2, dist = "poisson", rho = .4, corstr = "cs", wide = TRUE,cnames='y21,y22')
  # y2 <-as.matrix(y2[,c('y21','y22')])
  # y11<-y1[,1];y12<-y1[,2];y21<-y2[,1];y22<-y2[,2]
  # Y <- c(y11,y12,y21,y22)
  df.cor = data.frame(Y,X,Z,G)
  #---------------------------------
  #loglikelihood own define
  #---------------------------------
  negll <- function(param) { 
    -sum( param[1]*y11-exp(param[1])+sum(param[1:3])*y12-exp(sum(param[1:3]))+(param[1]+param[2]+param[4])*y21-exp(param[1]+param[2]+param[4])+(param[1]+param[3]+param[4])*y22-exp(param[1]+param[3]+param[4]) )}
  ABBA<-optim(param <- c(0.45, -0.05, 0.35, -0.25), negll, hessian=TRUE)
  MLE.optim[i,]<-ABBA$par
  I.optim<-I.optim+(ABBA$hessian)/(2*seq)
  invI.optim <-invI.optim+solve((ABBA$hessian)/(2*seq))
  #---------------------------------
  #GLM
  #---------------------------------
  mod.1 <- glm(Y ~ X + Z + G, family = poisson(link = "log"), df.cor)
  MLE.glm[i,]<-coef(mod.1)
  I.glm<-I.glm+solve(vcov(mod.1))/(2*seq)
  invI.glm <-invI.glm+vcov(mod.1)*(2*seq)
  #---------------------------------
  #closeform
  #---------------------------------
  tao[i]<-sum(y11)
  eta[i]<-sqrt( sum(y12)*sum(y21)/(sum(y11)*sum(y22)))
  gam[i]<-sqrt( sum(y12)*sum(y22)/(sum(y11)*sum(y21)))
  del[i]<-sqrt( sum(y21)*sum(y22)/(sum(y11)*sum(y12)))
  
  i.tt<-( tao[i] + tao[i]*eta[i]*gam[i] + tao[i]*eta[i]*del[i]+tao[i]*gam[i]*del[i] )/2
  i.te<-( tao[i]*eta[i]*gam[i] + tao[i]*eta[i]*del[i] )/2
  i.tg<-( tao[i]*eta[i]*gam[i] + tao[i]*gam[i]*del[i] )/2
  i.td<-( tao[i]*eta[i]*del[i] + tao[i]*gam[i]*del[i] )/2
  i.eg<-tao[i]*eta[i]*gam[i]/2
  i.ed<-tao[i]*eta[i]*del[i]/2
  i.gd<-tao[i]*gam[i]*del[i]/2
  I<-matrix( c(i.tt, i.te, i.tg, i.td, i.te, i.te, i.eg, i.ed, 
               i.tg, i.eg, i.tg, i.gd, i.td, i.ed, i.gd, i.td),nrow=4, ncol=4, byrow = TRUE)
  I.cf<-I.cf+I
  
  invI.closeform <- invI.closeform + solve(I)
  v.tt = i.tt + cov(y11,y12)+cov(y21,y22)
  v.te = i.te + cov(y11,y12)/2+cov(y21,y22)/2
  v.tg = i.tg + cov(y11,y12)/2+cov(y21,y22)/2
  v.td = i.td + cov(y21,y22)
  
  v.ee = i.te
  v.eg = i.eg + cov(y21,y22)/2
  v.ed = i.ed + cov(y21,y22)/2
  
  v.gg = i.tg 
  v.gd = i.gd + cov(y21,y22)/2
  
  V.cf<-V.cf+matrix( c(v.tt, v.te, v.tg, v.td, v.te, v.ee, v.eg, v.ed, 
                       v.tg, v.eg, v.gg, v.gd, v.td, v.ed, v.gd, v.td),nrow=4, ncol=4, byrow = TRUE)
  
}

colMeans(MLE.optim)
I.optim/sim_time
invI.optim/sim_time

colMeans(MLE.glm)
I.glm/sim_time
invI.glm/sim_time

log(c(mean(tao),mean(eta),mean(gam),mean(del)))

I.cf/sim_time
invI.closeform/sim_time
V.cf/sim_time

(invI.optim/sim_time)%*%(V.cf/sim_time)%*%(invI.optim/sim_time)
(invI.glm/sim_time)%*%(V.cf/sim_time)%*%(invI.glm/sim_time)
(invI.closeform/sim_time)%*%(V.cf/sim_time)%*%(invI.closeform/sim_time)


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
  mle.values[2]<-0
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


result.ind <- list();result.cor <- list();pvalue.ind<-list();pvalue.cor<-list()
#simulation for ABBA
seq=50
MLE.ind<-matrix(0, nrow = sim_time, ncol = length(param_222))
MLE.cor<-matrix(0, nrow = sim_time, ncol = length(param_222))
MLE.null.ind<-matrix(0, nrow = sim_time, ncol = length(param_222))
MLE.null.cor<-matrix(0, nrow = sim_time, ncol = length(param_222))
I.ind<- 0 ; I.cor<- 0 ; V.ind<-0 ;V.cor<-0;invI.ind<-0;invI.cor<-0
p1<-c();p2<-c()
alp = c();eta = c();ga = c();del = c();alp.0 = c();ga.0 = c();del.0 = c()
seq=100
II = 0;VV = 0;I.inv = 0;I.inv_V_I.inv = 0
set.seed(7353)
for (i in 1:sim_time) {
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
  y11 = rpois(100, lambda = mean_true[1])
  y12 = rpois(100, lambda = mean_true[2])
  y21 = rpois(100, lambda = mean_true[3])
  y22 = rpois(100, lambda = mean_true[4])
  
  p1[i] = 0
  p2[i] = 0
  for (j in 1:seq) {
    if(data.ind[j,1]==1 && data.ind[j,2]==1)  p1[i] = p1[i]+1
    if(data.ind[j,3]==1 && data.ind[j,4]==1)  p2[i] = p2[i]+1
  }
  p1[i]=p1[i]/seq
  p2[i]=p2[i]/seq
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
  mean.est<-Mean.True(MLE.ind[i,],xmat_222)
  IV.ind.i<-Matrix.IV(cros.type='ABBA', mle.values=MLE.ind[i,], x.mat=xmat_222, seq.size=seq, data=data.ind)
  #store result of MLE,I,V,inv.I
  
  I.ind.i<-IV.ind.i$I.hat
  V.ind.i<-Matrix::forceSymmetric(IV.ind.i$V.hat,uplo="L")
  invI.ind.i<-solve(I.ind.i)
  I.ind<-I.ind+I.ind.i
  V.ind<-V.ind+V.ind.i
  invI.ind<-invI.ind+invI.ind.i
  ##another way
 
  alp[i] = mod.1$coefficients[1]
  eta[i] = mod.1$coefficients[2]
  ga[i] = mod.1$coefficients[3]
  del[i] = mod.1$coefficients[4]
  
  alp.0[i] = mod.0$coefficients[1]
  ga.0[i] = mod.0$coefficients[2]
  del.0[i] = mod.0$coefficients[3]
  
  i.aa = ( exp(alp[i]) + exp(alp[i]+ga[i]) +
             exp(alp[i]+del[i]) + exp(alp[i]+ga[i]+del[i]) )/2
  i.ee = ( exp(alp[i]+ga[i]) + exp(alp[i]+del[i]) )/2
  i.gg = ( exp(alp[i]+ga[i]) + exp(alp[i]+ga[i]+del[i]) )/2
  i.dd = ( exp(alp[i]+del[i]) + exp(alp[i]+ga[i]+del[i]) )/2
  i.eg = exp(alp[i]+ga[i])/2
  i.ed = exp(alp[i]+del[i])/2
  i.gd = exp(alp[i]+ga[i]+del[i])/2
  I = matrix( c(i.aa, i.ee, i.gg, i.dd, i.ee, i.ee, i.eg, i.ed, 
                i.gg, i.eg, i.gg, i.gd, i.dd, i.ed, i.gd, i.dd),
              nrow=4, ncol=4, byrow = TRUE)
  v.aa = i.aa + (p1[i]-exp(alp[i])*exp(alp[i]+ga[i])) +
    (p2[i]-exp(alp[i]+del[i])*exp(alp[i]+ga[i]+del[i]))
  v.ee = i.ee
  v.gg = i.gg
  v.dd = i.dd + 
    (p2[i]-exp(alp[i]+del[i])*exp(alp[i]+ga[i]+del[i]))
  v.ae = i.ee + 
    (p1[i]-exp(alp[i])*exp(alp[i]+ga[i]))/2 +
    (p2[i]-exp(alp[i]+del[i])*exp(alp[i]+ga[i]+del[i]))/2
  v.ag = i.gg + (p1[i]-exp(alp[i])*exp(alp[i]+ga[i]))/2 +
    (p2[i]-exp(alp[i]+del[i])*exp(alp[i]+ga[i]+del[i]))/2
  v.ad = i.dd + 
    (p2[i]-exp(alp[i]+del[i])*exp(alp[i]+ga[i]+del[i]))
  v.eg = i.eg + 
    (p2[i]-exp(alp[i]+del[i])*exp(alp[i]+ga[i]+del[i]))/2
  v.ed = i.ed + 
    (p2[i]-exp(alp[i]+del[i])*exp(alp[i]+ga[i]+del[i]))/2
  v.gd = i.gd + 
    (p2[i]-exp(alp[i]+del[i])*exp(alp[i]+ga[i]+del[i]))/2
  V = matrix( c(v.aa, v.ae, v.ag, v.ad, v.ae, v.ee, v.eg, v.ed, 
                v.ag, v.eg, v.gg, v.gd, v.ad, v.ed, v.gd, v.dd),
              nrow=4, ncol=4, byrow = TRUE)
  
  
  II = II + I
  VV = VV + V
  I.inv = I.inv + solve(I)
  I.inv_V_I.inv = I.inv_V_I.inv + solve(I) %*% V %*% solve(I)
  
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                  i.gd, i.dd), nrow = 3, ncol = 3)
  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  v.pe = matrix(c(v.ae, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.aa, v.ag, v.ad, v.ag, v.gg, v.gd, v.ad, 
                  v.gd, v.dd), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  }

II/sim_time
VV/sim_time
I.inv/sim_time
I.inv_V_I.inv/sim_time

I.ind/sim_time
V.ind/sim_time
(invI.ind/sim_time)%*%(V.ind/sim_time)%*%(invI.ind/sim_time)
set.seed(11018022)

cov(MLE.ind)*2*seq

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


