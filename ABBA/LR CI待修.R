rm(list=ls(all=TRUE)) 
require(faraway)
require(numDeriv)
require(MASS)
require(extraDistr)
require(simstudy)
library(emplik)
setwd("C:/Users/User/Documents/Study_CrossoverDesign/ABBA")
#===========================================================
#參數說明
# - sim_time:模擬次數
# - cor_par:rho
# - seq:每組人數，AB一組BA一組，共兩組。
# - param:參數真值(tao,eta(treatment effect),gamma(period effect),delta(sequence effect))
# - xmat:自變量的值
# - mean.true:平均數.真值
# - X、Z、G:自變量，用來fit glm
#===========================================================
#N=25,50的LRCI for rho 0.7
sim_time=5000;cor_par=0.7;seq=50
param=c(1.2,0,1.0,0.2)
xmat=matrix(c(1,1,1,1, 0,1,1,0, 0,1,0,1, 0,0,1,1), nrow = 4, ncol = 4,byrow = TRUE)
mean.true=exp(param%*%xmat)

# ABBA.sim<-function(seq){
#   print(seq)
#MLE, null MLE, Matrix I, Matrix V, invI
tao<-c();eta<-c();gam<-c();del<-c()
tao.0<-c();eta.0<-c();gam.0<-c();del.0<-c()
I.cf=0;V.cf=0;invI=0;I0.cf=0;V0.cf=0;invI0=0

#for testing
W.na1<-0;W.rb1<-0;LR.na1<-0;LR.rb1<-0;S.na1<-0;S.rb1<-0
#for C.I., AL, CP
W.naup<-c();W.nalw<-c();Wna.len<-c();Wna.cp<-c()
W.rbup<-c();W.rblw<-c();Wrb.len<-c();Wrb.cp<-c()
LR.naup<-c();LR.nalw<-c();LRna.len<-c();LRna.cp<-c()
LR.rbup<-c();LR.rblw<-c();LRrb.len<-c();LRrb.cp<-c()
S.naup<-c();S.nalw<-c();Sna.len<-c();Sna.cp<-c()
S.rbup<-c();S.rblw<-c();Srb.len<-c();Srb.cp<-c()
#variance of eta
var.na<-c();var.rb<-c()

#========================================================
#correlated data
#========================================================
set.seed(110225021)
for (i in 1:sim_time){
  #用copula生成相關性資料,rho=0.4
  m1 <- c(mean.true[1], mean.true[2]); m2 <- c(mean.true[3], mean.true[4]) # lambda for each new variable
  y1 <- genCorGen(seq, nvars = 2, params1 = m1, dist = "poisson", rho = cor_par, corstr = "cs", wide = TRUE,cnames='y11,y12')
  y1 <-as.matrix(y1[,c('y11','y12')])
  y2 <- genCorGen(seq, nvars = 2, params1 = m2, dist = "poisson", rho = cor_par, corstr = "cs", wide = TRUE,cnames='y21,y22')
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
  invI <- invI + solve(I)
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
  
  #=====================================
  #LR Test
  #=====================================
  l1 = lik(c(tao[i], eta[i], gam[i], del[i]) )
  l0 = lik(c(tao.0[i], 0, gam.0[i], del.0[i]) )
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2 * A/B * (l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  
  #=====================================
  #LR test C.I. lower and upper
  #=====================================
  #Naive
  #-------------------------------------
  LR.na <- function(e){
    Mod.0 <- glm(formula = Y~Z+G+offset(e*X),family = poisson(link = "log"), data = df.cor)
    a = Mod.0$coefficients[1]
    g = Mod.0$coefficients[2]
    d = Mod.0$coefficients[3]
    l.0 = lik(c(a, e, g, d))
    u = 2*(l1-l.0) - qchisq(0.95, 1)
    return(u)
  }
  
  LR.naup[i] = uniroot(LR.na, c(eta[i],eta[i]+1.5))$root
  LR.nalw[i] = uniroot(LR.na, c(eta[i]-1.5,eta[i]))$root
  LRna.len[i] = LR.naup[i] - LR.nalw[i]
  LRna.cp[i] = ifelse((LR.nalw[i] < 0 & LR.naup[i] > 0), 1, 0)
  #-------------------------------------
  #Robust
  #-------------------------------------
  LR.rb <- function(e){
    Mod.0 <- glm(formula = Y~Z+G+offset(e*X),family = poisson(link = "log"), data = df.cor)
    a = Mod.0$coefficients[1]
    g = Mod.0$coefficients[2]
    d = Mod.0$coefficients[3]
    l.0 = lik(c(a, e, g, d))
    u = 2*A/B*(l1-l.0) - qchisq(0.95, 1)
    return(u)
  }
  
  LR.rbup[i] = uniroot(LR.rb, c(0,eta[i]+2))$root
  LR.rblw[i] = uniroot(LR.rb, c(eta[i]-1.5,eta[i]))$root
  LRrb.len[i] = LR.rbup[i] - LR.rblw[i]
  LRrb.cp[i] = ifelse((LR.rblw[i] < 0 & LR.rbup[i] > 0), 1, 0)}

# Profile likelihood of sigma
R.sigma = Vectorize( function(sigma) return( lik(c(tao[i], sigma, gam[i], del[i]) ) ))
# Visualising the profile likelihood of sigma
curve(R.sigma,-5,2, n = 1000, lwd =3, col = "blue", cex.axis = 2, cex.lab = 1.5, main =  expression(paste("Profile likelihood of ", sigma), ylim = c(0,1)),
      xlab = ~sigma, ylab = "Profile")
#abline(h = 0.147, lwd = 2, col = "red")
abline(v = eta[i], lwd = 2, col = "purple", lty = 2)

#MLE  
mle<-c(tao.hat = mean(tao),eta.hat = mean(eta),gam.hat = mean(gam),del.hat =mean(del))

#var
Variance<-c(sv = 2*seq*var(eta), na = mean(var.na), rb = mean(var.rb))

#p-value
Alpha<-c( Wna = 1-W.na1/sim_time, Wrb=1-W.rb1/sim_time, LRna=1-LR.na1/sim_time, LRrb=1-LR.rb1/sim_time, Sna=1-S.na1/sim_time, Srb=1-S.rb1/sim_time)

#Wald C.I.
Wna.CI<-c( lb = mean(W.nalw),ub = mean(W.naup), AL = mean(Wna.len), cp = sum(Wna.cp)/sim_time)
Wrb.CI<-c( lb = mean(W.rblw), ub = mean(W.rbup), AL = mean(Wrb.len), cp = sum(Wrb.cp)/sim_time)

#LR C.I.
LRna.CI<-c( lb = mean(LR.nalw),ub = mean(LR.naup), AL = mean(LRna.len), cp = sum(LRna.cp)/sim_time)
LRrb.CI<-c( lb = mean(LR.rblw),ub = mean(LR.rbup), AL = mean(LRrb.len), cp = sum(LRrb.cp)/sim_time)

#Score C.I.
Sna.CI<-c( lb = mean(S.nalw), ub = mean(S.naup), AL = mean(Sna.len), cp = sum(Sna.cp)/sim_time)
Srb.CI<-c( lb = mean(S.rblw), ub = mean(S.rbup), AL = mean(Srb.len), cp = sum(Srb.cp)/sim_time)

df=c( N = seq, MLE = mle, Var = Variance, a = Alpha, Wna = Wna.CI, Wrb = Wrb.CI ,LRna =LRna.CI, LRrb=LRrb.CI, Sna =Sna.CI, LRrb=Srb.CI )

# return(df)
# }



res<-lapply(c(75,100,125,150), ABBA.sim)
library(readr)
write_csv(data.frame(df), "C:\\Users\\User\\Documents\\Study_CrossoverDesign\\ABBA\\rho0.7.csv")