rm(list=ls(all=TRUE)) 
require(faraway)
require(numDeriv)
require(MASS)
require(extraDistr)
require(simstudy)
setwd("C:/Users/User/Documents/Study_CrossoverDesign/RCode")
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
sim_time=10000;seq=100
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
W.naup<-c();W.nalw<-c();Wna.len<-c();Wna.cp<-c()
W.rbup<-c();W.rblw<-c();Wrb.len<-c();Wrb.cp<-c()
LR.naup<-c();LR.nalw<-c();LRna.len<-c();LRna.cp<-c()
LR.rbup<-c();LR.rblw<-c();LRrb.len<-c();LRrb.cp<-c()
S.naup<-c();S.nalw<-c();Sna.len<-c();Sna.cp<-c()
S.rbup<-c();S.rblw<-c();Srb.len<-c();Srb.cp<-c()
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
  X = c(rep(0,seq), rep(1,2*seq), rep(0,seq))
  Z = c(rep(0,seq), rep(1,seq), rep(0,seq), rep(1,seq))
  G = c(rep(0,2*seq), rep(1,2*seq))
  Y <- c(y11,y12,y21,y22)
  df.ind = data.frame(Y,X,Z,G)
  #-------------------------------------------------------
  #closeform MLE:exp(.)估計量為指數形式e.g. eta[i]=exp(eta)
  #-------------------------------------------------------
  tao[i]<-mean(y11)
  eta[i]<-sqrt( sum(y12)*sum(y21)/(sum(y11)*sum(y22)) )
  gam[i]<-sqrt( sum(y12)*sum(y22)/(sum(y11)*sum(y21)) )
  del[i]<-sqrt( sum(y21)*sum(y22)/(sum(y11)*sum(y12)) )
  #-------------------------------------------------------
  #Matrix I & Matrix V
  #-------------------------------------------------------
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
  invI.closeform <- invI.closeform + solve(I)
  
  cov1<-mean( (y11-tao[i])*(y12-tao[i]*eta[i]*gam[i]) )
  cov2<-mean( (y21-tao[i]*eta[i]*del[i])*(y22-tao[i]*gam[i]*del[i]) )
  
  v.tt = i.tt + cov1+cov2
  v.ee = i.ee
  v.gg = i.gg 
  v.dd = i.dd + cov2
  v.te = i.ee + cov1/2+cov2/2
  v.tg = i.gg + cov1/2+cov2/2
  v.td = i.dd + cov1
  v.eg = i.eg + cov2/2
  v.ed = i.ed + cov2/2
  v.gd = i.gd + cov2/2
  
  V.cf<-V.cf+matrix( c(v.tt, v.te, v.tg, v.td,
                       v.te, v.ee, v.eg, v.ed, 
                       v.tg, v.eg, v.gg, v.gd,
                       v.td, v.ed, v.gd, v.td),nrow=4, ncol=4, byrow = TRUE)
  #-------------------------------------------------------
  #計算AB矩陣
  #-------------------------------------------------------
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.tt, i.gg, i.dd, i.gg, i.gg, i.gd, 
                  i.dd, i.gd, i.dd), nrow = 3, ncol = 3)
  
  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  v.pe = matrix(c(v.te, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.tt, v.tg, v.td, v.tg, v.gg, v.gd, v.td, 
                  v.gd, v.td), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  var.na[i]<-1/A
  var.rb[i]<-B/A/A
  
  #===================================
  #likelihood, par:[tao,eta,gam,del]
  #===================================
  l<-function(par){ 
    ll=sum(par[1]*y11-exp(par[1])+(par[1]+par[2]+par[3])*y12-exp(par[1]+par[2]+par[3]))+
      sum((par[1]+par[2]+par[4])*y21-exp(par[1]+par[2]+par[4])+ (par[1]+par[3]+par[4])*y22-exp(par[1]+par[3]+par[4]))
    return(ll)
  }
  #===================================
  #null MLE:用GLE算null之下的MLE
  #===================================
  mod.0 <- glm(Y ~ Z + G, family = poisson(link = "log"), df.ind)
  tao.0[i]<-exp(mod.0$coefficients[1])
  gam.0[i]<-exp(mod.0$coefficients[2])
  del.0[i]<-exp(mod.0$coefficients[3])
  
  #----------------------------------------
  #Matrix I0 & Matrix V0
  #----------------------------------------
  i0.tt<-( tao.0[i]+tao.0[i]*gam.0[i] + tao.0[i]*del[i]+tao[i]*gam[i]*del[i] )/2
  i0.gg<-( tao.0[i]*gam.0[i] + tao.0[i]*gam.0[i]*del.0[i] )/2 
  i0.ee<-( tao.0[i]*gam.0[i] + tao.0[i]*del.0[i] )/2 
  i0.dd<- ( tao.0[i]*del.0[i]+tao.0[i]*gam.0[i]*del.0[i] )/2
  i0.eg<- ( tao.0[i]*gam.0[i] )/2
  i0.ed<- ( tao.0[i]*del.0[i] )/2
  i0.gd<- ( tao.0[i]*gam.0[i]*del.0[i] )/2
  #----------------------------------------
  #計算null cov
  #----------------------------------------
  cov1<-mean( (y11-tao.0[i])*(y12-tao.0[i]*gam.0[i]) )
  cov2<-mean( (y21-tao.0[i]*del.0[i])*(y22-tao.0[i]*gam.0[i]*del.0[i]) )
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
  #----------------------------------------
  #計算A0B0
  #----------------------------------------
  i0.ep = matrix(c(i0.ee, i0.eg, i0.ed), ncol = 3)
  i0.pp = matrix(c(i0.tt, i0.gg, i0.dd, i0.gg, i0.gg, i0.gd, 
                   i0.dd, i0.gd, i0.dd), nrow = 3, ncol = 3)
  
  A0 = i0.ee - i0.ep %*% solve(i0.pp) %*% t(i0.ep)
  v0.pe = matrix(c(v0.te, v0.eg, v0.ed), nrow = 3)
  v0.pp = matrix(c(v0.tt, v0.tg, v0.td, v0.tg, v0.gg, v0.gd, v0.td, 
                   v0.gd, v0.td), nrow = 3, ncol = 3)
  B0 = v0.ee - 2 * i0.ep %*% solve(i0.pp) %*% v0.pe + 
    i0.ep %*% solve(i0.pp) %*% v0.pp %*% solve(i0.pp) %*% t(i0.ep)
  Eta<-log(eta[i])
  #=================================================
  #wald test
  #=================================================
  wna = log(eta[i]) * A * log(eta[i]) * 2*seq
  wrb = log(eta[i]) * A^2 / B * log(eta[i]) * 2*seq
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  
  #--------------------------------------------
  #wald confidence interval
  #--------------------------------------------
  #W.eta = Vectorize( function(e) return( (e-log(eta[i]))^2 ))
  
  W.na <- function(e){
    mod <- glm(formula = Y~Z+G+offset(e*X),family = poisson, df.ind)
    t = exp(mod$coefficients[1])
    g = exp(mod$coefficients[2])
    d = exp(mod$coefficients[3])
    ii.tt<-( t+t*exp(e)*g + t*exp(e)*d+t*g*d )/2
    ii.gg<-( t*exp(e)*g + t*g*d )/2 
    ii.ee<-( t*exp(e)*g + t*exp(e)*d )/2 
    ii.dd<- ( t*exp(e)*d+t*g*d )/2
    ii.eg<- ( t*exp(e)*g )/2
    ii.ed<- ( t*exp(e)*d )/2
    ii.gd<- ( t*g*d )/2
    ii.ep = matrix(c(ii.ee, ii.eg, ii.ed), ncol = 3)
    ii.pp = matrix(c(ii.tt, ii.gg, ii.dd, ii.gg, ii.gg, ii.gd, 
                    ii.dd, ii.gd, ii.dd), nrow = 3, ncol = 3)
    
    AA = ii.ee - ii.ep %*% solve(ii.pp) %*% t(ii.ep)
    
    w.na = 2*seq * (log(eta[i])-e)^2 * AA - qchisq(0.95, 1)
    return(w.na)
  }
  #W.na = Vectorize( function(e) return(2*seq * (log(eta[i])-e)^2 * A0 - qchisq(0.95, 1)))
  
  #lower and upper
  W.naup[i]<-uniroot(W.na,c(log(eta[i]),log(eta[i])+0.3))$root
  W.nalw[i]<-uniroot(W.na,c(log(eta[i])-0.3,log(eta[i])))$root
  Wna.len[i] = W.naup[i] -  W.nalw[i]
  Wna.cp[i] = ifelse((W.nalw[i] < 0 & W.naup[i] > 0), 1, 0)
  
  W.rb <- function(e){
    mod <- glm(formula = Y~Z+G+offset(e*X),family = poisson, df.ind)
    t = exp(mod$coefficients[1])
    g = exp(mod$coefficients[2])
    d = exp(mod$coefficients[3])
    ii.tt<-( t+t*exp(e)*g + t*exp(e)*d+t*g*d )/2
    ii.gg<-( t*exp(e)*g + t*g*d )/2 
    ii.ee<-( t*exp(e)*g + t*exp(e)*d )/2 
    ii.dd<- ( t*exp(e)*d+t*g*d )/2
    ii.eg<- ( t*exp(e)*g )/2
    ii.ed<- ( t*exp(e)*d )/2
    ii.gd<- ( t*g*d )/2
    
    ii.ep = matrix(c(ii.ee, ii.eg, ii.ed), ncol = 3)
    ii.pp = matrix(c(ii.tt, ii.gg, ii.dd, ii.gg, ii.gg, ii.gd, 
                     ii.dd, ii.gd, ii.dd), nrow = 3, ncol = 3)
    AA = ii.ee - ii.ep %*% solve(ii.pp) %*% t(ii.ep)
    
    cov.1<-mean( (y11-t)*(y12-t*exp(e)*g) )
    cov.2<-mean( (y21-t*exp(e)*d)*(y22-t*g*d) )
    
    vv.tt = i.tt + cov.1+cov.2
    vv.ee = i.ee
    vv.gg = i.gg 
    vv.dd = i.dd + cov.2
    vv.te = i.ee + cov.1/2+cov.2/2
    vv.tg = i.gg + cov.1/2+cov.2/2
    vv.td = i.dd + cov.1
    vv.eg = i.eg + cov.2/2
    vv.ed = i.ed + cov.2/2
    vv.gd = i.gd + cov.2/2
    
    vv.pe = matrix(c(vv.te, vv.eg, vv.ed), nrow = 3)
    vv.pp = matrix(c(vv.tt, vv.tg, vv.td, v.tg, v.gg, v.gd, v.td, 
                    vv.gd, vv.td), nrow = 3, ncol = 3)
    BB = vv.ee - 2 * ii.ep %*% solve(ii.pp) %*% vv.pe + 
      ii.ep %*% solve(ii.pp) %*% vv.pp %*% solve(ii.pp) %*% t(ii.ep)
    
    w.rb = 2*seq * (log(eta[i])-e)^2 * AA^2/BB - qchisq(0.95, 1)
    return(w.rb)
  }
  #W.rb = Vectorize( function(e) return(2*seq * (log(eta[i])-e)^2 * A0^2/B0 - qchisq(0.95, 1)))
  W.rbup[i]<-uniroot(W.rb,c(log(eta[i]),log(eta[i])+0.3))$root
  W.rblw[i]<-uniroot(W.rb,c(log(eta[i])-0.3,log(eta[i])))$root
  Wrb.len[i] = W.rbup[i] -  W.rblw[i]
  Wrb.cp[i] = ifelse((W.rblw[i] < 0 & W.rbup[i] > 0), 1, 0)
  
  
  # 
  #=================================================
  #LR Test
  #=================================================
  l1 = l(log(c(tao[i], eta[i], gam[i], del[i])) )
  l0 = l(log(c(tao.0[i], 1, gam.0[i], del.0[i])) )
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2*A/B*(l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  #--------------------------------------------
  #LR confidence interval
  #--------------------------------------------
  
  LR.na <- function(e){
    mod <- glm(Y ~ Z + G + offset(e*X), family = poisson(link = "log"), df.ind)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    l.0 = l(c(a, e, g, d))
    u = 2 * (l1-l.0) - qchisq(0.95, 1)
    return(u)
  }
  
  # LR.na<-Vectorize( function(e) return(2*(l1-l(log(c(tao.0[i], exp(e), gam.0[i], del.0[i])))) - qchisq(0.95, 1)))

  LR.naup[i] = uniroot(LR.na, c(log(eta[i]),log(eta[i])+0.2))$root
  LR.nalw[i] = uniroot(LR.na, c(log(eta[i])-0.2,log(eta[i])))$root
  LRna.len[i] = LR.naup[i] - LR.nalw[i]
  LRna.cp[i] = ifelse((LR.nalw[i] < 0 & LR.naup[i] > 0), 1, 0)
  
  LR.rb <- function(e){
    mod <- glm(Y ~ Z + G + offset(e*X), family = poisson(link = "log"), df.ind)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    l.0 = l(c(a, e, g, d))
    u = 2*A/B*(l1-l.0) - qchisq(0.95, 1)
    return(u)
  }
  #LR.rb<-Vectorize( function(e) return(2* A/B *(l1-l(log(c(tao.0[i], exp(e), gam.0[i], del.0[i])))) - qchisq(0.95, 1)))
  LR.rbup[i] = uniroot(LR.rb, c(log(eta[i]),log(eta[i])+0.15))$root
  LR.rblw[i] = uniroot(LR.rb, c(log(eta[i])-0.5,log(eta[i])))$root
  LRrb.len[i] = LR.rbup[i] - LR.rblw [i]
  LRrb.cp[i] = ifelse((LR.rblw[i] < 0 & LR.rbup[i] > 0), 1, 0)
  #==========================================
  #Score Test
  #==========================================
  s0 = seq * ( mean(y12) - tao.0[i]*gam.0[i] + mean(y21)- tao.0[i]*del.0[i] )
  sna = s0 / A0 * s0 / (2*seq)
  srb = s0 / B0 * s0 / (2*seq)
  if( sna <= qchisq(0.95, 1) )  S.na1 = S.na1+1
  if( srb <= qchisq(0.95, 1) )  S.rb1 = S.rb1+1
  # S.na <- function(e){
  #   
  #   mod <- glm(formula = Y~Z+G+offset(e*X),family = poisson, df.ind)
  #   t = exp(mod$coefficients[1])
  #   g = exp(mod$coefficients[2])
  #   d = exp(mod$coefficients[3])
  #   ii.tt<-( t+t*exp(e)*g + t*exp(e)*d+t*g*d )/2
  #   ii.gg<-( t*exp(e)*g + t*g*d )/2 
  #   ii.ee<-( t*exp(e)*g + t*exp(e)*d )/2 
  #   ii.dd<- ( t*exp(e)*d+t*g*d )/2
  #   ii.eg<- ( t*exp(e)*g )/2
  #   ii.ed<- ( t*exp(e)*d )/2
  #   ii.gd<- ( t*g*d )/2
  #   
  #   ii.ep = matrix(c(ii.ee, ii.eg, ii.ed), ncol = 3)
  #   ii.pp = matrix(c(ii.tt, ii.gg, ii.dd, ii.gg, ii.gg, ii.gd, 
  #                    ii.dd, ii.gd, ii.dd), nrow = 3, ncol = 3)
  #   AA = ii.ee - ii.ep %*% solve(ii.pp) %*% t(ii.ep)
  #   
  #   ss0 = seq * ( mean(y12) - t*g + mean(y21)- t*d )
  #   s.na = ss0 / AA * ss0 / (2*seq)- qchisq(0.95, 1)
  #   return(s.na)
  # }
  # ss0 = Vectorize( function(e) return( seq * ( mean(y12) - tao.0[i]*exp(e)*gam.0[i] + mean(y21)- tao.0[i]*exp(e)*del.0[i] )))
  # s.na = Vectorize( function(e) return(  ss0(e) / A0 * ss0(e) / (2*seq))- qchisq(0.95, 1))
  # 
  # S.na<-Vectorize( function(e) return(2*(l1-l(log(c(tao.0[i], exp(e), gam.0[i], del.0[i])))) - qchisq(0.95, 1)))
  # 
  # 
  # S.naup[i] = uniroot(S.na, c(log(eta[i]),log(eta[i])+0.1))$root
  # S.nalw[i] = uniroot(S.na, c(log(eta[i])-0.8,log(eta[i])))$root
  # Sna.len[i] = S.naup[i] - S.nalw[i]
  # 
  # S.rb <- function(e){
  #   mod <- glm(formula = Y~Z+G+offset(e*X),family = poisson, df.ind)
  #   t = exp(mod$coefficients[1])
  #   g = exp(mod$coefficients[2])
  #   d = exp(mod$coefficients[3])
  #   ii.tt<-( t+t*exp(e)*g + t*exp(e)*d+t*g*d )/2
  #   ii.gg<-( t*exp(e)*g + t*g*d )/2 
  #   ii.ee<-( t*exp(e)*g + t*exp(e)*d )/2 
  #   ii.dd<- ( t*exp(e)*d+t*g*d )/2
  #   ii.eg<- ( t*exp(e)*g )/2
  #   ii.ed<- ( t*exp(e)*d )/2
  #   ii.gd<- ( t*g*d )/2
  #   
  #   ii.ep = matrix(c(ii.ee, ii.eg, ii.ed), ncol = 3)
  #   ii.pp = matrix(c(ii.tt, ii.gg, ii.dd, ii.gg, ii.gg, ii.gd, 
  #                    ii.dd, ii.gd, ii.dd), nrow = 3, ncol = 3)
  #   AA = ii.ee - ii.ep %*% solve(ii.pp) %*% t(ii.ep)
  #   
  #   cov.1<-mean( (y11-t)*(y12-t*exp(e)*g) )
  #   cov.2<-mean( (y21-t*exp(e)*d)*(y22-t*g*d) )
  #   
  #   vv.tt = i.tt + cov.1+cov.2
  #   vv.ee = i.ee
  #   vv.gg = i.gg 
  #   vv.dd = i.dd + cov.2
  #   vv.te = i.ee + cov.1/2+cov.2/2
  #   vv.tg = i.gg + cov.1/2+cov.2/2
  #   vv.td = i.dd + cov.1
  #   vv.eg = i.eg + cov.2/2
  #   vv.ed = i.ed + cov.2/2
  #   vv.gd = i.gd + cov.2/2
  #   
  #   vv.pe = matrix(c(vv.te, vv.eg, vv.ed), nrow = 3)
  #   vv.pp = matrix(c(vv.tt, vv.tg, vv.td, v.tg, v.gg, v.gd, v.td, 
  #                    vv.gd, vv.td), nrow = 3, ncol = 3)
  #   BB = vv.ee - 2 * ii.ep %*% solve(ii.pp) %*% vv.pe + 
  #     ii.ep %*% solve(ii.pp) %*% vv.pp %*% solve(ii.pp) %*% t(ii.ep)
  #   
  #   
  #  
  #   ss0 = seq * ( mean(y12) - t*g + mean(y21)- t*d )
  #   s.rb = ss0 / BB * ss0 / (2*seq)- qchisq(0.95, 1)
  #   
  #   return(s.rb)
  # }
  # S.rbup[i] = uniroot(S.rb, c(log(eta[i]),log(eta[i])+0.2))$root
  # S.rblw[i] = uniroot(S.rb, c(log(eta[i])-0.2,log(eta[i])))$root
  # Srb.len[i] = S.rbup[i] - S.rblw[i]
  # curve(s.na,-0.1,0.1, n = 1000, lwd =3, col = "blue", cex.axis = 2, cex.lab = 1.5, main =  expression(paste("Wald function of ", mu), ylim = c(0,1)),
  #       xlab = ~mu, ylab = "LR")
  # abline(h = qchisq(0.95, df = 1)*A, lwd = 2, col = "red")
  # abline(v = log(eta[i]), lwd = 2, col = "purple", lty = 2)
  
}

#MLE
log(c(mean(tao),mean(log(etap[i])),mean(gam),mean(del)))
#sample variance of MLE
cov.m<-2*seq*matrix(c(var(log(tao)),cov(log(tao),log(eta)),cov(log(tao),log(gam)),cov(log(tao),log(del)), 
                      cov(log(tao),log(eta)),var(log(eta)),cov(log(gam),log(eta)),cov(log(del),log(eta)), 
                      cov(log(tao),log(gam)),cov(log(gam),log(eta)),var(log(gam)),cov(log(gam),log(del)), 
                      cov(log(tao),log(del)),cov(log(del),log(eta)),cov(log(gam),log(del)),var(log(del))), nrow = 4, ncol = 4,byrow = TRUE)
#matrix I hat、V hat、invI*V*invI
I.cf/sim_time
V.cf/sim_time
(invI.closeform/sim_time) %*% (V.cf/sim_time) %*% (invI.closeform/sim_time)

2*seq*var(log(eta))
mean(var.na)
mean(var.rb)
mean(eta.var)*2*seq
#p-value
1-W.na1/sim_time
1-W.rb1/sim_time
1-LR.na1/sim_time
1-LR.rb1/sim_time
1-S.na1/sim_time
1-S.rb1/sim_time

#Wald Naive C.I.
mean(W.nalw)
mean(W.naup)
mean(Wna.len)
sum(Wna.cp)/sim_time
#Wald Robust C.I.
mean(W.rblw)
mean(W.rbup)
mean(Wrb.len)
sum(Wrb.cp)/sim_time
#LR Naive C.I.
mean(LR.nalw)
mean(LR.naup)
mean(LRna.len)
sum(LRna.cp)/sim_time
#LR Robust C.I.
mean(LR.rblw)
mean(LR.rbup)
mean(LRrb.len)
sum(LRrb.cp)/sim_time
#========================================================
#correlated data
#========================================================

set.seed(110225021)
for (i in 1:sim_time){
  
  #用copula生成相關性資料,rho=0.4
  l1 <- c(mean.true[1], mean.true[2]); l2 <- c(mean.true[3], mean.true[4]) # lambda for each new variable
  y1 <- genCorGen(seq, nvars = 2, params1 = l1, dist = "poisson", rho = .4, corstr = "cs", wide = TRUE,cnames='y11,y12')
  y1 <-as.matrix(y1[,c('y11','y12')])
  y2 <- genCorGen(seq, nvars = 2, params1 = l2, dist = "poisson", rho = .4, corstr = "cs", wide = TRUE,cnames='y21,y22')
  y2 <-as.matrix(y2[,c('y21','y22')])
  y11<-y1[,1];y12<-y1[,2];y21<-y2[,1];y22<-y2[,2]
  Y <- c(y11,y12,y21,y22)
  
  #-------------------------------------------------------
  #closeform MLE:exp(.)估計量為指數形式e.g. eta[i]=exp(eta)
  #-------------------------------------------------------
  
  tao[i]<-mean(y11)
  eta[i]<-sqrt( sum(y12)*sum(y21)/(sum(y11)*sum(y22)) )
  gam[i]<-sqrt( sum(y12)*sum(y22)/(sum(y11)*sum(y21)) )
  del[i]<-sqrt( sum(y21)*sum(y22)/(sum(y11)*sum(y12)) )
  
  
  
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
  invI.closeform <- invI.closeform + solve(I)
  #確認covariance和R計算的covariance是否相同，下方的V用自己算的covariance而非R內建的
  cov1<-mean( (y11-tao[i])*(y12-tao[i]*eta[i]*gam[i]) )
  cov2<-mean( (y21-tao[i]*eta[i]*del[i])*(y22-tao[i]*gam[i]*del[i]) )
  
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
  
  #計算AB矩陣
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.tt, i.gg, i.dd, i.gg, i.gg, i.gd, 
                  i.dd, i.gd, i.dd), nrow = 3, ncol = 3)
  
  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  v.pe = matrix(c(v.te, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.tt, v.tg, v.td, v.tg, v.gg, v.gd, v.td, 
                  v.gd, v.td), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  var.na[i]<-1/A
  var.rb[i]<-B/A/A
  
  #=================================================
  #wald test
  #=================================================
  wna = log(eta[i]) * A * log(eta[i]) * 2*seq
  wrb = log(eta[i]) * A^2 / B * log(eta[i]) * 2*seq
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  #--------------------------------------------
  #wald confidence interval
  #--------------------------------------------
  lik<-function(par){ 
    ll=sum(log(par[1])*y11-par[1]+log(par[1]*par[2]*par[3])*y12-par[1]*par[2]*par[3])+
      sum(log(par[1]*par[2]*par[4])*y21-par[1]*par[2]*par[4]+ log(par[1]*par[3]*par[4])*y22-par[1]*par[3]*par[4])
    return(ll)
  }
  #null MLE:用GLE算null之下的MLE
  Y <- c(y11,y12,y21,y22)
  df.ind = data.frame(Y,X,Z,G)
  mod.0 <- glm(Y ~ Z + G, family = poisson(link = "log"), df.ind)
  tao.0[i]<-exp(mod.0$coefficients[1])
  gam.0[i]<-exp(mod.0$coefficients[2])
  del.0[i]<-exp(mod.0$coefficients[3])
  #計算null cov
  cov1<-mean( (y11-tao.0[i])*(y12-tao.0[i]*gam.0[i]) )
  cov2<-mean( (y21-tao.0[i]*del.0[i])*(y22-tao.0[i]*gam.0[i]*del.0[i]) )
  
  l1 = lik(c(tao[i], log(eta[i])[i], gam[i], del[i]) )
  l0 = lik(c(tao.0[i], 1, gam.0[i], del.0[i]))
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2*A/B*(l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  
  l<-function(par){ 
    ll=sum(par[1]*y11-exp(par[1])+(par[1]+par[2]+par[3])*y12-exp(par[1]+par[2]+par[3]))+
      sum((par[1]+par[2]+par[4])*y21-exp(par[1]+par[2]+par[4])+ (par[1]+par[3]+par[4])*y22-exp(par[1]+par[3]+par[4]))
    return(ll)
    
  }
  
  # Profile likelihood of exp(log(eta[i]))
  # R.e = Vectorize( function(e) return(  ( l(log(c(tao[i], log(eta[i])[i], gam[i], del[i])) )/l(log(c(tao.0[i], e, gam.0[i], del.0[i]))) ) ))
  # Lambda function
  # log(eta[i]) <-  Vectorize( function(e) 2*R.e(e) )
  #----------------------------------------
  # Wald approximation function for mu
  #-----------------------------------------
  # hess.e <-  hessian(log(eta[i]) , x = Eta)
  # 95% Confidence interval
  # lev <- 2*qchisq(0.95, df = 1)/hess.e
  W.eta = Vectorize( function(e) return( (e-Eta)^2 ))
  # W.etaC = Vectorize( function(e) return(W.eta(e)-lev))
  W.na = Vectorize( function(e) return(2*seq * W.eta(e) * A - qchisq(0.95, 1)))
  #lower and upper
  W.naup[i]<-uniroot(W.na,c(Eta,Eta+4))$root
  W.nalw[i]<-uniroot(W.na,c(Eta-4,Eta))$root
  Wna.len[i] = W.naup[i] -  W.nalw[i]
  Wna.cp[i] = ifelse((W.nalw[i] < 0 & W.naup[i] > 0), 1, 0)
  
  W.rb = Vectorize( function(e) return(2*seq * W.eta(e) * A^2/B - qchisq(0.95, 1)))
  W.rbup[i]<-uniroot(W.rb,c(Eta,Eta+4))$root
  W.rblw[i]<-uniroot(W.rb,c(Eta-4,Eta))$root
  # Wrb.len[i] = W.rbup[i] -  W.rblw[i]
  # Wrb.cp[i] = ifelse((W.rblw[i] < 0 & W.rbup[i] > 0), 1, 0)
  # 
  # # Profile likelihood of exp(eta)
  # R.e = Vectorize( function(e) return(  ( l(log(c(tao[i], eta[i], gam[i], del[i])) )/l(log(c(tao.0[i], e, gam.0[i], del.0[i]))) ) ))
  # # Lambda function
  # Eta <-  Vectorize( function(e) 2*R.e(e) )
  # #----------------------------------------
  # # Wald approximation function for mu
  # #-----------------------------------------
  # hess.e <-  hessian(Eta , x = eta[i])
  # # Threshold
  lev <- qchisq(0.95, df = 1)*B/A
  # W.mu = Vectorize( function(e) return( (e-Eta)^2 ))
  # 
  curve(W.rb,-0.2,0.2, n = 1000, lwd =3, col = "blue", cex.axis = 2, cex.lab = 1.5, main =  expression(paste("Wald function of ", e)),
      xlab = ~e, ylab = "Wald")
  abline(h = lev, lwd = 2, col = "red")
  abline(v = Eta, lwd = 2, col = "purple", lty = 2)
}

set.seed(110225021)
for (i in 1:sim_time){
  
  
  #用copula生成相關性資料,rho=0.4
  l1 <- c(mean.true[1], mean.true[2]); l2 <- c(mean.true[3], mean.true[4]) # lambda for each new variable
  y1 <- genCorGen(seq, nvars = 2, params1 = l1, dist = "poisson", rho = .4, corstr = "cs", wide = TRUE,cnames='y11,y12')
  y1 <-as.matrix(y1[,c('y11','y12')])
  y2 <- genCorGen(seq, nvars = 2, params1 = l2, dist = "poisson", rho = .4, corstr = "cs", wide = TRUE,cnames='y21,y22')
  y2 <-as.matrix(y2[,c('y21','y22')])
  y11<-y1[,1];y12<-y1[,2];y21<-y2[,1];y22<-y2[,2]
  Y <- c(y11,y12,y21,y22)
  
  #---------------------------------
  #closeform MLE:exp(.)
  #---------------------------------
  # mod.1 <- glm(Y ~ X + Z + G, family = poisson(link = "log"), df.cor)
  # tao[i]<-exp(mod.1$coefficients[1])
  # eta[i]<-exp(mod.1$coefficients[2])
  # gam[i]<-exp(mod.1$coefficients[3])
  # del[i]<-exp(mod.1$coefficients[4])
  
  tao[i]<-mean(y11)
  eta[i]<-sqrt( sum(y12)*sum(y21)/(sum(y11)*sum(y22)) )
  gam[i]<-sqrt( sum(y12)*sum(y22)/(sum(y11)*sum(y21)) )
  del[i]<-sqrt( sum(y21)*sum(y22)/(sum(y11)*sum(y12)) )
  
  eta.var[i]<-0.25*(1/sum(y11)+1/sum(y12)+1/sum(y21)+1/sum(y22))
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
  invI.closeform <- invI.closeform + solve(I)
  
  cov1<-mean( (y11-tao[i])*(y12-tao[i]*eta[i]*gam[i]) )
  cov2<-mean( (y21-tao[i]*eta[i]*del[i])*(y22-tao[i]*gam[i]*del[i]) )
  
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
  
  var.na[i]<-1/A
  var.rb[i]<-B/A/A
  
  #wald test
  wna = Eta * A * Eta * 2*seq
  wrb = Eta * A^2 / B * Eta * 2*seq
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  
  #log LR test :par=exp(.)
  lik<-function(par){ 
    ll=sum(log(par[1])*y11-par[1]+log(par[1]*par[2]*par[3])*y12-par[1]*par[2]*par[3])+
      sum(log(par[1]*par[2]*par[4])*y21-par[1]*par[2]*par[4]+ log(par[1]*par[3]*par[4])*y22-par[1]*par[3]*par[4])
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
  
  
  W.na <- function(e){
    
    mod <- glm(Y ~ Z + G + offset(e*X), family = poisson(link = "log"), df.cor)
    t<-exp(mod$coefficients[1])
    g<-exp(mod$coefficients[2])
    d<-exp(mod$coefficients[3])
    i.tt<-( t+t*e*g + t*e*d+t*g*d )/2
    i.gg<-( t*e*g + t*g*d )/2 
    
    i.ee<-( t*e*g + t*e*d )/2 
    
    i.dd<- ( t*e*d+t*g*d )/2
    i.eg<- ( t*e*g )/2
    i.ed<- ( t*e*d )/2
    i.gd<- ( t*g*d )/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.tt, i.gg, i.dd, i.gg, i.gg, i.gd, 
                    i.dd, i.gd, i.dd), nrow = 3, ncol = 3)
    
    Aa = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
    
    w.na = 2*seq * (Eta-e)^2 * Aa - qchisq(0.95, 1)
    return(w.na)
  }
  
  #lower and upper
  W.naup[i]<-uniroot(W.na,c(Eta,Eta+2))$root
  W.nalw[i]<-uniroot(W.na,c(Eta-2,Eta))$root
  Wna.len[i] = W.naup[i] -  W.nalw[i]
  Wna.cp[i] = ifelse((W.nalw[i] < 0 & W.naup[i] > 0), 1, 0)
  
  W.rb <- function(e){
    mod <- glm(Y ~ Z + G + offset(e*X), family = poisson(link = "log"), df.cor)
    t<-exp(mod$coefficients[1])
    g<-exp(mod$coefficients[2])
    d<-exp(mod$coefficients[3])
    i.tt<-( t+t*e*g + t*e*d+t*g*d )/2
    i.gg<-( t*e*g + t*g*d )/2 
    
    i.ee<-( t*e*g + t*e*d )/2 
    
    i.dd<- ( t*e*d+t*g*d )/2
    i.eg<- ( t*e*g )/2
    i.ed<- ( t*e*d )/2
    i.gd<- ( t*g*d )/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.tt, i.gg, i.dd, i.gg, i.gg, i.gd, 
                    i.dd, i.gd, i.dd), nrow = 3, ncol = 3)
    
    Aa = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
    
    c1<-mean( (y11-t)*(y12-t*e*g) )
    c2<-mean( (y21-t*e*d)*(y22-t*g*d) )
    
    v.tt = i.tt + c1+c2#cov(y11,y12)+cov(y21,y22)
    v.ee = i.ee
    v.gg = i.gg 
    v.dd = i.dd + c2#cov(y21,y22)
    
    v.te = i.ee + c1/2+c2/2#cov(y11,y12)/2+cov(y21,y22)/2
    v.tg = i.gg + c1/2+c2/2#cov(y11,y12)/2+cov(y21,y22)/2
    v.td = i.dd + c1#cov(y21,y22)
    
    
    v.eg = i.eg + c2/2#cov(y21,y22)/2
    v.ed = i.ed + c2/2#cov(y21,y22)/2
    v.gd = i.gd + c2/2#cov(y21,y22)/2
    
    v.pe = matrix(c(v.te, v.eg, v.ed), nrow = 3)
    v.pp = matrix(c(v.tt, v.tg, v.td, v.tg, v.gg, v.gd, v.td, 
                    v.gd, v.td), nrow = 3, ncol = 3)
    Ba = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
      i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
    
    
    w.rb = 2*seq * (eta[i]-e)^2 * Aa^2/Ba - qchisq(0.95, 1)
    return(w.rb)
  }
  #lower and upper
  W.rbup[i]<-uniroot(W.rb,c(Eta,Eta+1))$root
  W.rblw[i]<-uniroot(W.rb,c(Eta-1,Eta))$root
  Wrb.len[i] = W.rbup[i] -  W.rblw[i]
  Wrb.cp[i] = ifelse((W.rblw[i] < 0 & W.rbup[i] > 0), 1, 0)
}


#MLE
log(c(mean(tao),mean(eta),mean(g),mean(del)))
#sample variance of MLE
cov.m<-2*seq*matrix(c(var(log(tao)),cov(log(tao),log(eta)),cov(log(tao),log(g)),cov(log(tao),log(del)), 
                      cov(log(tao),log(eta)),var(log(eta)),cov(log(g),log(eta)),cov(log(del),log(eta)), 
                      cov(log(tao),log(g)),cov(log(g),log(eta)),var(log(g)),cov(log(g),log(del)), 
                      cov(log(tao),log(del)),cov(log(del),log(eta)),cov(log(g),log(del)),var(log(del))), nrow = 4, ncol = 4,byrow = TRUE)
#matrix I hat、V hat、invI*V*invI
I.cf/sim_time
V.cf/sim_time
(invI.closeform/sim_time) %*% (V.cf/sim_time) %*% (invI.closeform/sim_time)

2*seq*var(log(eta))
mean(var.na)
mean(var.rb)
mean(eta.var)*2*seq
#p-value
1-W.na1/sim_time
1-W.rb1/sim_time
1-LR.na1/sim_time
1-LR.rb1/sim_time
1-S.na1/sim_time
1-S.rb1/sim_time