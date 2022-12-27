rm(list=ls(all=TRUE)) 
require(faraway)
require(numDeriv)
require(MASS)
require(extraDistr)
require(simstudy)
setwd("C:/Users/User/Documents/Study_CrossoverDesign/ABBA")
#===========================================================
#參數說明
# - sim_time:模擬次數
# - cor_par:用random effect生成相關性資料的參數Gamma(alpha,beta)的beta(alpha*beta=1)
# - seq:每組人數，AB一組BA一組，共兩組。
# - param:參數真值(tao,eta(treatment effect),gamma(period effect),delta(sequence effect))
# - xmat:自變量的值
# - mean.true:平均數.真值
# - X、Z、G:自變量，用來fit glm
#===========================================================
sim_time=10000;seq=150
param=c(1.2,0,1.0,0.2)#c(1.2,0,1.0,0.2)#c(0.3,0,0.4,-1.0)
xmat=matrix(c(1,1,1,1, 0,1,1,0, 0,1,0,1, 0,0,1,1), nrow = 4, ncol = 4,byrow = TRUE)
mean.true=exp(param%*%xmat)

MLE.optim<-matrix(0, nrow = sim_time, ncol = 4)
MLE.glm<-matrix(0, nrow = sim_time, ncol = 4)
MLE.closeform<-matrix(0, nrow = sim_time, ncol = 4)
invI.optim=0;I.glm=0;I.optim=0;invI.glm=0

I.cf=0;V.cf=0;invI.closeform=0;I0.cf=0;V0.cf=0;invI0.closeform=0
tao<-c();eta<-c();gam<-c();del<-c()
tao.0<-c();eta.0<-c();gam.0<-c();del.0<-c()
W.na1<-0;W.rb1<-0;LR.na1<-0;LR.rb1<-0;S.na1<-0;S.rb1<-0
eta.var<-c();var.na<-c();var.rb<-c()
#========================================================
#independent data
#========================================================
set.seed(110225021)
for (i in 1:sim_time){
  #mean=r*(1-p)/p用負二項分配生成個數型資料
  # y11<-rnbinom(n = seq, size = 2, prob = 0.6)
  # y12<-rnbinom(n = seq, size = 3, prob = 0.6)
  # y21<-rnbinom(n = seq, size = 2, prob = 0.8)
  # y22<-rnbinom(n = 100000, size = 3, prob = 0.8)
  # Y <- c(y11,y12,y21,y22)
  
  #用卜瓦松分配生成個數型資料
  y11<-rpois(seq, lambda = mean.true[1])
  y12<-rpois(seq, lambda = mean.true[2])
  y21<-rpois(seq, lambda = mean.true[3])
  y22<-rpois(seq, lambda = mean.true[4])
  Y <- c(y11,y12,y21,y22)
  X = c(rep(0,seq), rep(1,2*seq), rep(0,seq))
  Z = c(rep(0,seq), rep(1,seq), rep(0,seq), rep(1,seq))
  G = c(rep(0,2*seq), rep(1,2*seq))
  df.ind = data.frame(Y,X,Z,G)
  #---------------------------------
  #MLE
  #---------------------------------
  mod.1 <- glm(Y ~ X + Z + G, family = poisson(link = "log"), df.ind)
  tao[i]<-mod.1$coefficients[1]
  eta[i]<-mod.1$coefficients[2]
  gam[i]<-mod.1$coefficients[3]
  del[i]<-mod.1$coefficients[4]
  #---------------------------------
  #null MLE
  #---------------------------------
  mod.0 <- glm(Y ~ Z + G, family = poisson(link = "log"), df.ind)
  tao.0[i]<-mod.0$coefficients[1]
  gam.0[i]<-mod.0$coefficients[2]
  del.0[i]<-mod.0$coefficients[3]
  
  #---------------------------------
  #Matrix I & Matrix V
  #---------------------------------
  i.tt<-( exp(tao[i]) + exp(tao[i]+eta[i]+gam[i]) + exp(tao[i]+eta[i]+del[i]) + exp(tao[i]+gam[i]+del[i]) )/2
  i.ee<-( exp(tao[i]+eta[i]+gam[i]) + exp(tao[i]+eta[i]+del[i]) )/2 
  i.gg<-( exp(tao[i]+eta[i]+gam[i]) + exp(tao[i]+gam[i]+del[i]) )/2 
  i.dd<- ( exp(tao[i]+eta[i]+del[i])+exp(tao[i]+gam[i]+del[i]) )/2
  i.eg<- ( exp(tao[i]+eta[i]+gam[i]) )/2
  i.ed<- ( exp(tao[i]+eta[i]+del[i]) )/2
  i.gd<- ( exp(tao[i]+gam[i]+del[i]) )/2
  I<-matrix(c(i.tt,i.ee,i.gg,i.dd,
              i.ee,i.ee,i.eg,i.ed,
              i.gg,i.eg,i.gg,i.gd,
              i.dd,i.ed,i.gd,i.dd),nrow=4, ncol=4, byrow = TRUE)
  I.cf<-I.cf+I
  invI.closeform <- invI.closeform + solve(I)
  #---------------------------------
  #Matrix A 
  #---------------------------------
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.tt, i.gg, i.dd, i.gg, i.gg, i.gd, 
                  i.dd, i.gd, i.dd), nrow = 3, ncol = 3)
  
  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  #---------------------------------
  #Matrix V
  #---------------------------------
  v.tt = i.tt + cov(y11,y12)+cov(y21,y22)
  v.ee = i.ee
  v.gg = i.gg 
  v.dd = i.dd + cov(y21,y22)
  v.te = i.ee + cov(y11,y12)/2 + cov(y21,y22)/2
  v.tg = i.gg + cov(y11,y12)/2 + cov(y21,y22)/2
  v.td = i.dd + cov(y21,y22)
  v.eg = i.eg + cov(y21,y22)/2
  v.ed = i.ed + cov(y21,y22)/2
  v.gd = i.gd + cov(y21,y22)/2
  
  V.cf<-V.cf+matrix( c(v.tt, v.te, v.tg, v.td,
                       v.te, v.ee, v.eg, v.ed, 
                       v.tg, v.eg, v.gg, v.gd,
                       v.td, v.ed, v.gd, v.td),nrow=4, ncol=4, byrow = TRUE)
  
  #---------------------------------
  #Matrix B
  #---------------------------------
  v.pe = matrix(c(v.te, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.tt, v.tg, v.td, v.tg, v.gg, v.gd, v.td, 
                  v.gd, v.td), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  
  var.na[i]<-1/A
  var.rb[i]<-B/A/A
  #---------------------------------
  #log LR test 
  #---------------------------------
  lik<-function(par){ 
    ll=sum(par[1]*y11-exp(par[1])+(par[1]+par[2]+par[3])*y12-exp(par[1]+par[2]+par[3]))+
      sum((par[1]+par[2]+par[4])*y21-exp(par[1]+par[2]+par[4])+ (par[1]+par[3]+par[4])*y22-exp(par[1]+par[3]+par[4]))
    return(ll)
  }
  
  #---------------------------------
  #LR Test
  #---------------------------------
  l1 = lik(c(tao[i], eta[i], gam[i], del[i]) )
  l0 = lik(c(tao.0[i], 0, gam.0[i], del.0[i]) )
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2 * A/B * (l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  #---------------------------------
  #null matrix I & matrix V 
  #---------------------------------
  i0.tt<-( exp(tao.0[i]) + exp(tao.0[i]+gam.0[i]) + exp(tao.0[i]+del[i]) + exp(tao[i]+gam[i]+del[i]) )/2
  i0.gg<-( exp(tao.0[i]+gam.0[i]) + exp(tao.0[i]+gam.0[i]+del.0[i]) )/2 
  i0.ee<-( exp(tao.0[i]+gam.0[i]) + exp(tao.0[i]+del.0[i]) )/2 
  i0.dd<- ( exp(tao.0[i]+del.0[i]) + exp(tao.0[i]+gam.0[i]+del.0[i]) )/2
  i0.eg<- ( exp(tao.0[i]+gam.0[i]) )/2
  i0.ed<- ( exp(tao.0[i]+del.0[i]) )/2
  i0.gd<- ( exp(tao.0[i]+gam.0[i]+del.0[i]) )/2
  I0<-matrix(c(i0.tt,i0.ee,i0.gg,i0.dd,
               i0.ee,i0.ee,i0.eg,i0.ed,
               i0.gg,i0.eg,i0.gg,i0.gd,
               i0.dd,i0.ed,i0.gd,i0.dd),nrow=4, ncol=4, byrow = TRUE)
  I0.cf<-I0.cf+I0
  invI0.closeform <- invI0.closeform + solve(I0)
  
  v0.tt = i0.tt + cov(y11,y12)+cov(y21,y22)
  v0.ee = i0.ee
  v0.gg = i0.gg 
  v0.dd = i0.dd + cov(y21,y22)
  v0.te = i0.ee + cov(y11,y12)/2 + cov(y21,y22)/2
  v0.tg = i0.gg + cov(y11,y12)/2 + cov(y21,y22)/2
  v0.td = i0.dd + cov(y21,y22)
  v0.eg = i0.eg + cov(y21,y22)/2
  v0.ed = i0.ed + cov(y21,y22)/2
  v0.gd = i0.gd + cov(y21,y22)/2
  
  V0.cf<-V0.cf+matrix( c(v0.tt, v0.te, v0.tg, v0.td,
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
  #---------------------------------
  #score test
  #---------------------------------
  s0 = seq * ( mean(y12) - exp(tao.0[i]+gam.0[i]) + mean(y21)- exp(tao.0[i]+del.0[i]) )
  sna = s0 / A0 * s0 / (2*seq)
  srb = s0 / B0 * s0 / (2*seq)
  if( sna <= qchisq(0.95, 1) )  S.na1 = S.na1+1
  if( srb <= qchisq(0.95, 1) )  S.rb1 = S.rb1+1
  #---------------------------------
  #wald test
  #---------------------------------
  wna = eta[i] * A0 * eta[i] * 2*seq
  wrb = eta[i] * A0^2 / B0 * eta[i] * 2*seq
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  
}

#========================================================
#correlated data
#========================================================
set.seed(110225021)
for (i in 1:sim_time){
  
  #用copula生成相關性資料,rho=0.7
  l1 <- c(mean.true[1], mean.true[2]); l2 <- c(mean.true[3], mean.true[4]) # lambda for each new variable
  y1 <- genCorGen(seq, nvars = 2, params1 = l1, dist = "poisson", rho = .7, corstr = "cs", wide = TRUE,cnames='y11,y12')
  y1 <-as.matrix(y1[,c('y11','y12')])
  y2 <- genCorGen(seq, nvars = 2, params1 = l2, dist = "poisson", rho = .7, corstr = "cs", wide = TRUE,cnames='y21,y22')
  y2 <-as.matrix(y2[,c('y21','y22')])
  y11<-y1[,1];y12<-y1[,2];y21<-y2[,1];y22<-y2[,2]
  Y <- c(y11,y12,y21,y22)
  X = c(rep(0,seq), rep(1,2*seq), rep(0,seq))
  Z = c(rep(0,seq), rep(1,seq), rep(0,seq), rep(1,seq))
  G = c(rep(0,2*seq), rep(1,2*seq))
  df.cor = data.frame(Y,X,Z,G)
  #---------------------------------
  #MLE
  #---------------------------------
  mod.1 <- glm(Y ~ X + Z + G, family = poisson(link = "log"), df.cor)
  tao[i]<-mod.1$coefficients[1]
  eta[i]<-mod.1$coefficients[2]
  gam[i]<-mod.1$coefficients[3]
  del[i]<-mod.1$coefficients[4]
  #---------------------------------
  #null MLE
  #---------------------------------
  mod.0 <- glm(Y ~ Z + G, family = poisson(link = "log"), df.cor)
  tao.0[i]<-mod.0$coefficients[1]
  gam.0[i]<-mod.0$coefficients[2]
  del.0[i]<-mod.0$coefficients[3]
  
  #---------------------------------
  #Matrix I & Matrix V
  #---------------------------------
  i.tt<-( exp(tao[i]) + exp(tao[i]+eta[i]+gam[i]) + exp(tao[i]+eta[i]+del[i]) + exp(tao[i]+gam[i]+del[i]) )/2
  i.ee<-( exp(tao[i]+eta[i]+gam[i]) + exp(tao[i]+eta[i]+del[i]) )/2 
  i.gg<-( exp(tao[i]+eta[i]+gam[i]) + exp(tao[i]+gam[i]+del[i]) )/2 
  i.dd<- ( exp(tao[i]+eta[i]+del[i])+exp(tao[i]+gam[i]+del[i]) )/2
  i.eg<- ( exp(tao[i]+eta[i]+gam[i]) )/2
  i.ed<- ( exp(tao[i]+eta[i]+del[i]) )/2
  i.gd<- ( exp(tao[i]+gam[i]+del[i]) )/2
  I<-matrix(c(i.tt,i.ee,i.gg,i.dd,
              i.ee,i.ee,i.eg,i.ed,
              i.gg,i.eg,i.gg,i.gd,
              i.dd,i.ed,i.gd,i.dd),nrow=4, ncol=4, byrow = TRUE)
  I.cf<-I.cf+I
  invI.closeform <- invI.closeform + solve(I)
  #---------------------------------
  #Matrix A 
  #---------------------------------
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.tt, i.gg, i.dd, i.gg, i.gg, i.gd, 
                  i.dd, i.gd, i.dd), nrow = 3, ncol = 3)
  
  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  #---------------------------------
  #Matrix V
  #---------------------------------
  v.tt = i.tt + cov(y11,y12)+cov(y21,y22)
  v.ee = i.ee
  v.gg = i.gg 
  v.dd = i.dd + cov(y21,y22)
  v.te = i.ee + cov(y11,y12)/2 + cov(y21,y22)/2
  v.tg = i.gg + cov(y11,y12)/2 + cov(y21,y22)/2
  v.td = i.dd + cov(y21,y22)
  v.eg = i.eg + cov(y21,y22)/2
  v.ed = i.ed + cov(y21,y22)/2
  v.gd = i.gd + cov(y21,y22)/2
  
  V.cf<-V.cf+matrix( c(v.tt, v.te, v.tg, v.td,
                       v.te, v.ee, v.eg, v.ed, 
                       v.tg, v.eg, v.gg, v.gd,
                       v.td, v.ed, v.gd, v.td),nrow=4, ncol=4, byrow = TRUE)
  
  #---------------------------------
  #Matrix B
  #---------------------------------
  v.pe = matrix(c(v.te, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.tt, v.tg, v.td, v.tg, v.gg, v.gd, v.td, 
                  v.gd, v.td), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  
  var.na[i]<-1/A
  var.rb[i]<-B/A/A
  #---------------------------------
  #log LR test 
  #---------------------------------
  lik<-function(par){ 
    ll=sum(par[1]*y11-exp(par[1])+(par[1]+par[2]+par[3])*y12-exp(par[1]+par[2]+par[3]))+
       sum((par[1]+par[2]+par[4])*y21-exp(par[1]+par[2]+par[4])+ (par[1]+par[3]+par[4])*y22-exp(par[1]+par[3]+par[4]))
    return(ll)
  }
  
  #---------------------------------
  #LR Test
  #---------------------------------
  l1 = lik(c(tao[i], eta[i], gam[i], del[i]) )
  l0 = lik(c(tao.0[i], 0, gam.0[i], del.0[i]) )
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2 * A/B * (l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  #---------------------------------
  #null matrix I & matrix V 
  #---------------------------------
  i0.tt<-( exp(tao.0[i]) + exp(tao.0[i]+gam.0[i]) + exp(tao.0[i]+del[i]) + exp(tao[i]+gam[i]+del[i]) )/2
  i0.gg<-( exp(tao.0[i]+gam.0[i]) + exp(tao.0[i]+gam.0[i]+del.0[i]) )/2 
  i0.ee<-( exp(tao.0[i]+gam.0[i]) + exp(tao.0[i]+del.0[i]) )/2 
  i0.dd<- ( exp(tao.0[i]+del.0[i]) + exp(tao.0[i]+gam.0[i]+del.0[i]) )/2
  i0.eg<- ( exp(tao.0[i]+gam.0[i]) )/2
  i0.ed<- ( exp(tao.0[i]+del.0[i]) )/2
  i0.gd<- ( exp(tao.0[i]+gam.0[i]+del.0[i]) )/2
  I0<-matrix(c(i0.tt,i0.ee,i0.gg,i0.dd,
               i0.ee,i0.ee,i0.eg,i0.ed,
               i0.gg,i0.eg,i0.gg,i0.gd,
               i0.dd,i0.ed,i0.gd,i0.dd),nrow=4, ncol=4, byrow = TRUE)
  I0.cf<-I0.cf+I0
  invI0.closeform <- invI0.closeform + solve(I0)
  
  v0.tt = i0.tt + cov(y11,y12)+cov(y21,y22)
  v0.ee = i0.ee
  v0.gg = i0.gg 
  v0.dd = i0.dd + cov(y21,y22)
  v0.te = i0.ee + cov(y11,y12)/2 + cov(y21,y22)/2
  v0.tg = i0.gg + cov(y11,y12)/2 + cov(y21,y22)/2
  v0.td = i0.dd + cov(y21,y22)
  v0.eg = i0.eg + cov(y21,y22)/2
  v0.ed = i0.ed + cov(y21,y22)/2
  v0.gd = i0.gd + cov(y21,y22)/2
  
  V0.cf<-V0.cf+matrix( c(v0.tt, v0.te, v0.tg, v0.td,
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
  #---------------------------------
  #score test
  #---------------------------------
  s0 = seq * ( mean(y12) - exp(tao.0[i]+gam.0[i]) + mean(y21)- exp(tao.0[i]+del.0[i]) )
  sna = s0 / A0 * s0 / (2*seq)
  srb = s0 / B0 * s0 / (2*seq)
  if( sna <= qchisq(0.95, 1) )  S.na1 = S.na1+1
  if( srb <= qchisq(0.95, 1) )  S.rb1 = S.rb1+1
  #---------------------------------
  #wald test
  #---------------------------------
  wna = eta[i] * A0 * eta[i] * 2*seq
  wrb = eta[i] * A0^2 / B0 * eta[i] * 2*seq
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  
}



#sample variance of MLE
cov.m<-2*seq*matrix(c( var(tao),cov(tao,eta),cov(tao,gam),cov(tao,del), 
                      cov(tao,eta),var(eta),cov(gam,eta),cov(del,eta), 
                      cov(tao,gam),cov(gam,eta),var(gam),cov(gam,del), 
                      cov(tao,del),cov(del,eta),cov(gam,del),var(del) ), nrow = 4, ncol = 4,byrow = TRUE)
#matrix I hat、V hat、invI*V*invI
I.cf/sim_time
V.cf/sim_time
(invI.closeform/sim_time) %*% (V.cf/sim_time) %*% (invI.closeform/sim_time)

#MLE
c(mean(tao),mean(eta),mean(gam),mean(del))

2*seq*var(eta)
mean(var.na)
mean(var.rb)

#p-value
1-W.na1/sim_time
1-W.rb1/sim_time
1-LR.na1/sim_time
1-LR.rb1/sim_time
1-S.na1/sim_time
1-S.rb1/sim_time
