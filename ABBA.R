f <- function(x){
  y = log(x/(1-x))
  return(y)
}
f.inv <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}
m <- function(x){
  y = exp(x)/(1+exp(x))^2
  return(y)
}

#without correlation
n = 500
r = 1000
q = 1.5
X = c(rep(0,n), rep(1,2*n), rep(0,n))
Z = c(rep(0,n), rep(1,n), rep(0,n), rep(1,n))
G = c(rep(0,2*n), rep(1,2*n))
alp = c()
eta = c()
ga = c()
del = c()
alp.0 = c()
ga.0 = c()
del.0 = c()
p1 <- c()
p2 <- c()
II = 0
VV = 0
I.inv = 0
I.inv_V_I.inv = 0
LR.na1 = 0
LR.rb1 = 0
LR.na2 = 0
LR.rb2 = 0
LRna.lb <- c()
LRna.ub <- c()
LRrb.lb <- c()
LRrb.ub <- c()
LRna.len <- c()
LRrb.len <- c()
S.na1 = 0
S.rb1 = 0
S.na2 = 0
S.rb2 = 0
Sna.lb <- c()
Sna.ub <- c()
Srb.lb <- c()
Srb.ub <- c()
Sna.len <- c()
Srb.len <- c()
W.na1 = 0
W.rb1 = 0
W.na2 = 0
W.rb2 = 0
Wna.lb <- c()
Wna.ub <- c()
Wrb.lb <- c()
Wrb.ub <- c()
Wna.len <- c()
Wrb.len <- c()
set.seed(110225007)
for (i in 1:r) {
  y11 = rbinom(n,1,f.inv(0.5))
  y12 = rbinom(n,1,f.inv(0.7))
  y21 = rbinom(n,1,f.inv(1.3))
  y22 = rbinom(n,1,f.inv(1.5))
  
  p1[i] = 0
  p2[i] = 0
  for (j in 1:n) {
    if(y11[j]==1 && y12[j]==1)  p1[i] = p1[i]+1
    if(y21[j]==1 && y22[j]==1)  p2[i] = p2[i]+1
  }
  p1[i] = p1[i]/n
  p2[i] = p2[i]/n
  
  Y <- c(y11,y12,y21,y22)
  D = data.frame(Y,X,Z,G)
  lmod <- glm(cbind(Y,1-Y)~X+Z+G, family = binomial, D)
  alp[i] = lmod$coefficients[1]
  eta[i] = lmod$coefficients[2]
  ga[i] = lmod$coefficients[3]
  del[i] = lmod$coefficients[4]
  lmod.0 <- glm(cbind(Y,1-Y)~Z+G, family = binomial, D)
  alp.0[i] = lmod.0$coefficients[1]
  ga.0[i] = lmod.0$coefficients[2]
  del.0[i] = lmod.0$coefficients[3]
  
  i.aa = ( m(alp[i]) + m(alp[i]+ga[i]) +
             m(alp[i]+del[i]) + m(alp[i]+ga[i]+del[i]) )/2
  i.ee = ( m(alp[i]+ga[i]) + m(alp[i]+del[i]) )/2
  i.gg = ( m(alp[i]+ga[i]) + m(alp[i]+ga[i]+del[i]) )/2
  i.dd = ( m(alp[i]+del[i]) + m(alp[i]+ga[i]+del[i]) )/2
  i.eg = m(alp[i]+ga[i])/2
  i.ed = m(alp[i]+del[i])/2
  i.gd = m(alp[i]+ga[i]+del[i])/2
  I = matrix( c(i.aa, i.ee, i.gg, i.dd, i.ee, i.ee, i.eg, i.ed, 
                i.gg, i.eg, i.gg, i.gd, i.dd, i.ed, i.gd, i.dd),
              nrow=4, ncol=4, byrow = TRUE)
  v.aa = i.aa + (p1[i]-f.inv(alp[i])*f.inv(alp[i]+ga[i])) +
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))
  v.ee = i.ee
  v.gg = i.gg
  v.dd = i.dd + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))
  v.ae = i.ee + 
    (p1[i]-f.inv(alp[i])*f.inv(alp[i]+ga[i]))/2 +
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
  v.ag = i.gg + (p1[i]-f.inv(alp[i])*f.inv(alp[i]+ga[i]))/2 +
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
  v.ad = i.dd + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))
  v.eg = i.eg + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
  v.ed = i.ed + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
  v.gd = i.gd + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
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
  
  #log LR test
  lik <- function(t){
    u = t[1]*sum(y11) - n*log(1+exp(t[1])) + (t[1]+t[2]+t[3])*sum(y12) - 
      n*log(1+exp(t[1]+t[2]+t[3])) + (t[1]+t[2]+t[4])*sum(y21) - 
      n*log(1+exp(t[1]+t[2]+t[4])) + (t[1]+t[3]+t[4])*sum(y22) - 
      n*log(1+exp(t[1]+t[3]+t[4]))
    return(u)
  }
  l1 = lik(c(alp[i], eta[i], ga[i], del[i]))
  l0 = lik(c(alp.0[i], 0, ga.0[i], del.0[i]))
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2*A/B*(l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  
  lr.na <- function(e){
    mod.0 <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
                 family = binomial, data = D)
    a = mod.0$coefficients[1]
    g = mod.0$coefficients[2]
    d = mod.0$coefficients[3]
    l.0 = lik(c(a, e, g, d))
    u = 2*(l1-l.0) - qchisq(0.95, 1)
    return(u)
  }
  lrna1 = uniroot(lr.na, c(eta[i]-q,eta[i]))
  lrna2 = uniroot(lr.na, c(eta[i],eta[i]+q))
  LRna.lb[i] = lrna1$root
  LRna.ub[i] = lrna2$root
  LRna.len[i] = LRna.ub[i] - LRna.lb[i]
  
  lr.rb <- function(e){
    mod.0 <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
                 family = binomial, data = D)
    a = mod.0$coefficients[1]
    g = mod.0$coefficients[2]
    d = mod.0$coefficients[3]
    l0 = lik(c(a, e, g, d))
    u = 2*A/B*(l1-l0) - qchisq(0.95, 1)
    return(u)
  }
  lrrb1 = uniroot(lr.rb, c(eta[i]-1.5,eta[i]))
  lrrb2 = uniroot(lr.rb, c(eta[i],eta[i]+1.5))
  LRrb.lb[i] = lrrb1$root
  LRrb.ub[i] = lrrb2$root
  LRrb.len[i] = LRrb.ub[i] - LRrb.lb[i]
  
  if(LRna.lb[i]<0 && LRna.ub[i]>0)  LR.na2 = LR.na2+1
  if(LRrb.lb[i]<0 && LRrb.ub[i]>0)  LR.rb2 = LR.rb2+1
  
  #score test
  i.aa = ( m(alp.0[i]) + m(alp.0[i]+ga.0[i]) +
             m(alp.0[i]+del.0[i]) + m(alp.0[i]+ga.0[i]+del.0[i]) )/2
  i.ee = ( m(alp.0[i]+ga.0[i]) + m(alp.0[i]+del.0[i]) )/2
  i.gg = ( m(alp.0[i]+ga.0[i]) + m(alp.0[i]+ga.0[i]+del.0[i]) )/2
  i.dd = ( m(alp.0[i]+del.0[i]) + m(alp.0[i]+ga.0[i]+del.0[i]) )/2
  i.eg = m(alp.0[i]+ga.0[i])/2
  i.ed = m(alp.0[i]+del.0[i])/2
  i.gd = m(alp.0[i]+ga.0[i]+del.0[i])/2
  
  v.aa = i.aa + (p1[i]-f.inv(alp.0[i])*f.inv(alp.0[i]+ga.0[i])) +
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))
  v.ee = i.ee
  v.gg = i.gg
  v.dd = i.dd + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))
  v.ae = i.ee + 
    (p1[i]-f.inv(alp.0[i])*f.inv(alp.0[i]+ga.0[i]))/2 +
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  v.ag = i.gg + (p1[i]-f.inv(alp.0[i])*f.inv(alp.0[i]+ga.0[i]))/2 +
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  v.ad = i.dd + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))
  v.eg = i.eg + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  v.ed = i.ed + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  v.gd = i.gd + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                  i.gd, i.dd), nrow = 3, ncol = 3)
  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  v.pe = matrix(c(v.ae, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.aa, v.ag, v.ad, v.ag, v.gg, v.gd, v.ad, 
                  v.gd, v.dd), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  
  s0 = n * ( mean(y12) - f.inv(alp.0[i]+ga.0[i]) +
               mean(y21) - f.inv(alp.0[i]+del.0[i]) )
  sna = s0 / A * s0 / (2*n)
  srb = s0 / B * s0 / (2*n)
  if( sna <= qchisq(0.95, 1) )  S.na1 = S.na1+1
  if( srb <= qchisq(0.95, 1) )  S.rb1 = S.rb1+1
  
  st.na <- function(e){
    mod <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
               family = binomial, D)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    i.aa = ( m(a) + m(a+e+g) + m(a+e+d) + m(a+g+d) )/2
    i.ee = ( m(a+e+g) + m(a+e+d) )/2
    i.gg = ( m(a+e+g) + m(a+g+d) )/2
    i.dd = ( m(a+e+d) + m(a+g+d) )/2
    i.eg = m(a+e+g)/2
    i.ed = m(a+e+d)/2
    i.gd = m(a+g+d)/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                    i.gd, i.dd), nrow = 3, ncol = 3)
    A0 = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
    
    s0 = n * ( mean(y12) - f.inv(a+e+g) +
                 mean(y21) - f.inv(a+e+d) )
    s.na = s0 / A0 * s0 / (2*n) - qchisq(0.95, 1)
    return(s.na)
  }
  scna1 = uniroot(st.na, c(eta[i]-q,eta[i]))
  scna2 = uniroot(st.na, c(eta[i],eta[i]+q))
  Sna.lb[i] = scna1$root
  Sna.ub[i] = scna2$root
  Sna.len[i] = Sna.ub[i] - Sna.lb[i]
  
  st.rb <- function(e){
    mod <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
               family = binomial, D)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    i.aa = ( m(a) + m(a+e+g) + m(a+e+d) + m(a+g+d) )/2
    i.ee = ( m(a+e+g) + m(a+e+d) )/2
    i.gg = ( m(a+e+g) + m(a+g+d) )/2
    i.dd = ( m(a+e+d) + m(a+g+d) )/2
    i.eg = m(a+e+g)/2
    i.ed = m(a+e+d)/2
    i.gd = m(a+g+d)/2
    v.aa = i.aa + (p1[i]-f.inv(a)*f.inv(a+e+g)) +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.ee = i.ee
    v.gg = i.gg
    v.dd = i.dd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.ae = i.ee + (p1[i]-f.inv(a)*f.inv(a+e+g))/2 +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ag = i.gg + (p1[i]-f.inv(a)*f.inv(a+e+g))/2 +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ad = i.dd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.eg = i.eg + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ed = i.ed + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.gd = i.gd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                    i.gd, i.dd), nrow = 3, ncol = 3)
    v.pe = matrix(c(v.ae, v.eg, v.ed), nrow = 3)
    v.pp = matrix(c(v.aa, v.ag, v.ad, v.ag, v.gg, v.gd, v.ad, 
                    v.gd, v.dd), nrow = 3, ncol = 3)
    B0 = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
      i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
    
    s0 = n * ( mean(y12) - f.inv(a+e+g) +
                 mean(y21) - f.inv(a+e+d) )
    s.rb = s0 / B0 * s0 / (2*n) - qchisq(0.95, 1)
    return(s.rb)
  }
  scrb1 = uniroot(st.rb, c(eta[i]-q,eta[i]))
  scrb2 = uniroot(st.rb, c(eta[i],eta[i]+q))
  Srb.lb[i] = scrb1$root
  Srb.ub[i] = scrb2$root
  Srb.len[i] = Srb.ub[i] - Srb.lb[i]
  
  if(Sna.lb[i]<0 && Sna.ub[i]>0)  S.na2 = S.na2+1
  if(Srb.lb[i]<0 && Srb.ub[i]>0)  S.rb2 = S.rb2+1
  
  #wald test
  wna = eta[i] * A * eta[i] * 2*n
  wrb = eta[i] * A^2 / B * eta[i] * 2*n
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  
  wt.na <- function(e){
    mod <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
               family = binomial, D)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    i.aa = ( m(a) + m(a+e+g) + m(a+e+d) + m(a+g+d) )/2
    i.ee = ( m(a+e+g) + m(a+e+d) )/2
    i.gg = ( m(a+e+g) + m(a+g+d) )/2
    i.dd = ( m(a+e+d) + m(a+g+d) )/2
    i.eg = m(a+e+g)/2
    i.ed = m(a+e+d)/2
    i.gd = m(a+g+d)/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                    i.gd, i.dd), nrow = 3, ncol = 3)
    A0 = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
    
    w.na = 2*n * (eta[i]-e)^2 * A0 - qchisq(0.95, 1)
    return(w.na)
  }
  wdna1 = uniroot(wt.na, c(eta[i]-q,eta[i]))
  wdna2 = uniroot(wt.na, c(eta[i],eta[i]+q))
  Wna.lb[i] = wdna1$root
  Wna.ub[i] = wdna2$root
  Wna.len[i] = Wna.ub[i] - Wna.lb[i]
  
  wt.rb <- function(e){
    mod <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
               family = binomial, D)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    i.aa = ( m(a) + m(a+e+g) + m(a+e+d) + m(a+g+d) )/2
    i.ee = ( m(a+e+g) + m(a+e+d) )/2
    i.gg = ( m(a+e+g) + m(a+g+d) )/2
    i.dd = ( m(a+e+d) + m(a+g+d) )/2
    i.eg = m(a+e+g)/2
    i.ed = m(a+e+d)/2
    i.gd = m(a+g+d)/2
    
    v.aa = i.aa + (p1[i]-f.inv(a)*f.inv(a+e+g)) +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.ee = i.ee
    v.gg = i.gg
    v.dd = i.dd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.ae = i.ee + (p1[i]-f.inv(a)*f.inv(a+e+g))/2 +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ag = i.gg + (p1[i]-f.inv(a)*f.inv(a+e+g))/2 +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ad = i.dd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.eg = i.eg + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ed = i.ed + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.gd = i.gd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                    i.gd, i.dd), nrow = 3, ncol = 3)
    A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
    v.pe = matrix(c(v.ae, v.eg, v.ed), nrow = 3)
    v.pp = matrix(c(v.aa, v.ag, v.ad, v.ag, v.gg, v.gd, v.ad, 
                    v.gd, v.dd), nrow = 3, ncol = 3)
    B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
      i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
    
    s0 = n * ( mean(y21) - f.inv(a+e+g) +
                 mean(y12) - f.inv(a+e+d) )
    w.rb = 2*n * (eta[i]-e)^2 * A^2/B - qchisq(0.95, 1)
    return(w.rb)
  }
  wdrb1 = uniroot(wt.rb, c(eta[i]-q,eta[i]))
  wdrb2 = uniroot(wt.rb, c(eta[i],eta[i]+q))
  Wrb.lb[i] = wdrb1$root
  Wrb.ub[i] = wdrb2$root
  Wrb.len[i] = Wrb.ub[i] - Wrb.lb[i]
  
  if(Wna.lb[i]<0 && Wna.ub[i]>0)  W.na2 = W.na2+1
  if(Wrb.lb[i]<0 && Wrb.ub[i]>0)  W.rb2 = W.rb2+1
}
c(mean(alp), mean(eta), mean(ga), mean(del))
c(mean(alp.0), mean(ga.0), mean(del.0))
c(mean(p1), mean(p2))
II/r
VV/r
I.inv/r
I.inv_V_I.inv/r
matrix(c(var(alp), cov(alp,eta), cov(alp,ga), cov(alp,del),
         cov(alp,eta), var(eta), cov(eta,ga), cov(eta,del),
         cov(alp,ga), cov(eta,ga), var(ga), cov(ga,del),
         cov(alp,del), cov(eta,del), cov(ga,del), var(del)),
       nrow = 4, ncol = 4) * (2*n)
c(LR.na1, LR.rb1, LR.na2, LR.rb2)/r
c(mean(LRna.lb), mean(LRna.ub))
c(mean(LRrb.lb), mean(LRrb.ub))
c(mean(LRna.len), mean(LRrb.len))
c(S.na1, S.rb1, S.na2, S.rb2)/r
c(mean(Sna.lb), mean(Sna.ub))
c(mean(Srb.lb), mean(Srb.ub))
c(mean(Sna.len), mean(Srb.len))
c(W.na1, W.rb1, W.na2, W.rb2)/r
c(mean(Wna.lb), mean(Wna.ub))
c(mean(Wrb.lb), mean(Wrb.ub))
c(mean(Wna.len), mean(Wrb.len))





#with correlation
library(MASS)
n = 500
r = 1000
q = 1.5
X = c(rep(0,n), rep(1,2*n), rep(0,n))
Z = c(rep(0,n), rep(1,n), rep(0,n), rep(1,n))
G = c(rep(0,2*n), rep(1,2*n))
alp = c()
eta = c()
ga = c()
del = c()
alp.0 = c()
ga.0 = c()
del.0 = c()
p1 <- c()
p2 <- c()
II = 0
VV = 0
I.inv = 0
I.inv_V_I.inv = 0
LR.na1 = 0
LR.rb1 = 0
LR.na2 = 0
LR.rb2 = 0
LRna.lb <- c()
LRna.ub <- c()
LRrb.lb <- c()
LRrb.ub <- c()
LRna.len <- c()
LRrb.len <- c()
S.na1 = 0
S.rb1 = 0
S.na2 = 0
S.rb2 = 0
Sna.lb <- c()
Sna.ub <- c()
Srb.lb <- c()
Srb.ub <- c()
Sna.len <- c()
Srb.len <- c()
W.na1 = 0
W.rb1 = 0
W.na2 = 0
W.rb2 = 0
Wna.lb <- c()
Wna.ub <- c()
Wrb.lb <- c()
Wrb.ub <- c()
Wna.len <- c()
Wrb.len <- c()
cm1 <- matrix(c(1,0.7,0.7,1), nrow = 2, ncol = 2)
cm2 <- matrix(c(1,0.9,0.9,1), nrow = 2, ncol = 2)
set.seed(110225007)
for (i in 1:r) {
  y1 <- mvrnorm(n, mu=c(0,0), Sigma=cm1)
  y2 <- mvrnorm(n, mu=c(0,0), Sigma=cm2)
  y11 <- y1[1:n]
  y12 <- y1[(n+1):(2*n)]
  y21 <- y2[1:n]
  y22 <- y2[(n+1):(2*n)]
  for (j in 1:n) {
    if( (1-pnorm(y11[j])) < f.inv(0.5) )  y11[j] = 1
    else  y11[j] = 0
    if( (1-pnorm(y12[j])) < f.inv(0.7) )  y12[j] = 1
    else  y12[j] = 0
    if( (1-pnorm(y21[j])) < f.inv(1.3) )  y21[j] = 1
    else  y21[j] = 0
    if( (1-pnorm(y22[j])) < f.inv(1.5) )  y22[j] = 1
    else  y22[j] = 0
  }
  
  p1[i] = 0
  p2[i] = 0
  for (j in 1:n) {
    if(y11[j]==1 && y12[j]==1)  p1[i] = p1[i]+1
    if(y21[j]==1 && y22[j]==1)  p2[i] = p2[i]+1
  }
  p1[i] = p1[i]/n
  p2[i] = p2[i]/n
  
  Y <- c(y11,y12,y21,y22)
  D = data.frame(Y,X,Z,G)
  lmod <- glm(cbind(Y,1-Y)~X+Z+G, family = binomial, D)
  alp[i] = lmod$coefficients[1]
  eta[i] = lmod$coefficients[2]
  ga[i] = lmod$coefficients[3]
  del[i] = lmod$coefficients[4]
  lmod.0 <- glm(cbind(Y,1-Y)~Z+G, family = binomial, D)
  alp.0[i] = lmod.0$coefficients[1]
  ga.0[i] = lmod.0$coefficients[2]
  del.0[i] = lmod.0$coefficients[3]
  
  i.aa = ( m(alp[i]) + m(alp[i]+ga[i]) +
             m(alp[i]+del[i]) + m(alp[i]+ga[i]+del[i]) )/2
  i.ee = ( m(alp[i]+ga[i]) + m(alp[i]+del[i]) )/2
  i.gg = ( m(alp[i]+ga[i]) + m(alp[i]+ga[i]+del[i]) )/2
  i.dd = ( m(alp[i]+del[i]) + m(alp[i]+ga[i]+del[i]) )/2
  i.eg = m(alp[i]+ga[i])/2
  i.ed = m(alp[i]+del[i])/2
  i.gd = m(alp[i]+ga[i]+del[i])/2
  I = matrix( c(i.aa, i.ee, i.gg, i.dd, i.ee, i.ee, i.eg, i.ed, 
                i.gg, i.eg, i.gg, i.gd, i.dd, i.ed, i.gd, i.dd),
              nrow=4, ncol=4, byrow = TRUE)
  v.aa = i.aa + (p1[i]-f.inv(alp[i])*f.inv(alp[i]+ga[i])) +
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))
  v.ee = i.ee
  v.gg = i.gg
  v.dd = i.dd + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))
  v.ae = i.ee + 
    (p1[i]-f.inv(alp[i])*f.inv(alp[i]+ga[i]))/2 +
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
  v.ag = i.gg + (p1[i]-f.inv(alp[i])*f.inv(alp[i]+ga[i]))/2 +
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
  v.ad = i.dd + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))
  v.eg = i.eg + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
  v.ed = i.ed + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
  v.gd = i.gd + 
    (p2[i]-f.inv(alp[i]+del[i])*f.inv(alp[i]+ga[i]+del[i]))/2
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
  
  #log LR test
  lik <- function(t){
    u = t[1]*sum(y11) - n*log(1+exp(t[1])) + (t[1]+t[2]+t[3])*sum(y12) - 
      n*log(1+exp(t[1]+t[2]+t[3])) + (t[1]+t[2]+t[4])*sum(y21) - 
      n*log(1+exp(t[1]+t[2]+t[4])) + (t[1]+t[3]+t[4])*sum(y22) - 
      n*log(1+exp(t[1]+t[3]+t[4]))
    return(u)
  }
  l1 = lik(c(alp[i], eta[i], ga[i], del[i]))
  l0 = lik(c(alp.0[i], 0, ga.0[i], del.0[i]))
  if( (2*(l1-l0))<=qchisq(0.95, 1) )  LR.na1 = LR.na1+1
  if( (2*A/B*(l1-l0))<=qchisq(0.95, 1) )  LR.rb1 = LR.rb1+1
  
  lr.na <- function(e){
    mod.0 <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
                 family = binomial, data = D)
    a = mod.0$coefficients[1]
    g = mod.0$coefficients[2]
    d = mod.0$coefficients[3]
    l.0 = lik(c(a, e, g, d))
    u = 2*(l1-l.0) - qchisq(0.95, 1)
    return(u)
  }
  lrna1 = uniroot(lr.na, c(eta[i]-q,eta[i]))
  lrna2 = uniroot(lr.na, c(eta[i],eta[i]+q))
  LRna.lb[i] = lrna1$root
  LRna.ub[i] = lrna2$root
  LRna.len[i] = LRna.ub[i] - LRna.lb[i]
  
  lr.rb <- function(e){
    mod.0 <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
                 family = binomial, data = D)
    a = mod.0$coefficients[1]
    g = mod.0$coefficients[2]
    d = mod.0$coefficients[3]
    l0 = lik(c(a, e, g, d))
    u = 2*A/B*(l1-l0) - qchisq(0.95, 1)
    return(u)
  }
  lrrb1 = uniroot(lr.rb, c(eta[i]-1.5,eta[i]))
  lrrb2 = uniroot(lr.rb, c(eta[i],eta[i]+1.5))
  LRrb.lb[i] = lrrb1$root
  LRrb.ub[i] = lrrb2$root
  LRrb.len[i] = LRrb.ub[i] - LRrb.lb[i]
  
  if(LRna.lb[i]<0 && LRna.ub[i]>0)  LR.na2 = LR.na2+1
  if(LRrb.lb[i]<0 && LRrb.ub[i]>0)  LR.rb2 = LR.rb2+1
  
  #score test
  i.aa = ( m(alp.0[i]) + m(alp.0[i]+ga.0[i]) +
             m(alp.0[i]+del.0[i]) + m(alp.0[i]+ga.0[i]+del.0[i]) )/2
  i.ee = ( m(alp.0[i]+ga.0[i]) + m(alp.0[i]+del.0[i]) )/2
  i.gg = ( m(alp.0[i]+ga.0[i]) + m(alp.0[i]+ga.0[i]+del.0[i]) )/2
  i.dd = ( m(alp.0[i]+del.0[i]) + m(alp.0[i]+ga.0[i]+del.0[i]) )/2
  i.eg = m(alp.0[i]+ga.0[i])/2
  i.ed = m(alp.0[i]+del.0[i])/2
  i.gd = m(alp.0[i]+ga.0[i]+del.0[i])/2
  
  v.aa = i.aa + (p1[i]-f.inv(alp.0[i])*f.inv(alp.0[i]+ga.0[i])) +
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))
  v.ee = i.ee
  v.gg = i.gg
  v.dd = i.dd + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))
  v.ae = i.ee + 
    (p1[i]-f.inv(alp.0[i])*f.inv(alp.0[i]+ga.0[i]))/2 +
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  v.ag = i.gg + (p1[i]-f.inv(alp.0[i])*f.inv(alp.0[i]+ga.0[i]))/2 +
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  v.ad = i.dd + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))
  v.eg = i.eg + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  v.ed = i.ed + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  v.gd = i.gd + 
    (p2[i]-f.inv(alp.0[i]+del.0[i])*f.inv(alp.0[i]+ga.0[i]+del.0[i]))/2
  
  i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
  i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                  i.gd, i.dd), nrow = 3, ncol = 3)
  A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
  v.pe = matrix(c(v.ae, v.eg, v.ed), nrow = 3)
  v.pp = matrix(c(v.aa, v.ag, v.ad, v.ag, v.gg, v.gd, v.ad, 
                  v.gd, v.dd), nrow = 3, ncol = 3)
  B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
    i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
  
  s0 = n * ( mean(y12) - f.inv(alp.0[i]+ga.0[i]) +
               mean(y21) - f.inv(alp.0[i]+del.0[i]) )
  sna = s0 / A * s0 / (2*n)
  srb = s0 / B * s0 / (2*n)
  if( sna <= qchisq(0.95, 1) )  S.na1 = S.na1+1
  if( srb <= qchisq(0.95, 1) )  S.rb1 = S.rb1+1
  
  st.na <- function(e){
    mod <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
               family = binomial, D)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    i.aa = ( m(a) + m(a+e+g) + m(a+e+d) + m(a+g+d) )/2
    i.ee = ( m(a+e+g) + m(a+e+d) )/2
    i.gg = ( m(a+e+g) + m(a+g+d) )/2
    i.dd = ( m(a+e+d) + m(a+g+d) )/2
    i.eg = m(a+e+g)/2
    i.ed = m(a+e+d)/2
    i.gd = m(a+g+d)/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                    i.gd, i.dd), nrow = 3, ncol = 3)
    A0 = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
    
    s0 = n * ( mean(y12) - f.inv(a+e+g) +
                 mean(y21) - f.inv(a+e+d) )
    s.na = s0 / A0 * s0 / (2*n) - qchisq(0.95, 1)
    return(s.na)
  }
  scna1 = uniroot(st.na, c(eta[i]-q,eta[i]))
  scna2 = uniroot(st.na, c(eta[i],eta[i]+q))
  Sna.lb[i] = scna1$root
  Sna.ub[i] = scna2$root
  Sna.len[i] = Sna.ub[i] - Sna.lb[i]
  
  st.rb <- function(e){
    mod <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
               family = binomial, D)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    i.aa = ( m(a) + m(a+e+g) + m(a+e+d) + m(a+g+d) )/2
    i.ee = ( m(a+e+g) + m(a+e+d) )/2
    i.gg = ( m(a+e+g) + m(a+g+d) )/2
    i.dd = ( m(a+e+d) + m(a+g+d) )/2
    i.eg = m(a+e+g)/2
    i.ed = m(a+e+d)/2
    i.gd = m(a+g+d)/2
    v.aa = i.aa + (p1[i]-f.inv(a)*f.inv(a+e+g)) +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.ee = i.ee
    v.gg = i.gg
    v.dd = i.dd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.ae = i.ee + (p1[i]-f.inv(a)*f.inv(a+e+g))/2 +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ag = i.gg + (p1[i]-f.inv(a)*f.inv(a+e+g))/2 +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ad = i.dd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.eg = i.eg + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ed = i.ed + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.gd = i.gd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                    i.gd, i.dd), nrow = 3, ncol = 3)
    v.pe = matrix(c(v.ae, v.eg, v.ed), nrow = 3)
    v.pp = matrix(c(v.aa, v.ag, v.ad, v.ag, v.gg, v.gd, v.ad, 
                    v.gd, v.dd), nrow = 3, ncol = 3)
    B0 = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
      i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
    
    s0 = n * ( mean(y12) - f.inv(a+e+g) +
                 mean(y21) - f.inv(a+e+d) )
    s.rb = s0 / B0 * s0 / (2*n) - qchisq(0.95, 1)
    return(s.rb)
  }
  scrb1 = uniroot(st.rb, c(eta[i]-q,eta[i]))
  scrb2 = uniroot(st.rb, c(eta[i],eta[i]+q))
  Srb.lb[i] = scrb1$root
  Srb.ub[i] = scrb2$root
  Srb.len[i] = Srb.ub[i] - Srb.lb[i]
  
  if(Sna.lb[i]<0 && Sna.ub[i]>0)  S.na2 = S.na2+1
  if(Srb.lb[i]<0 && Srb.ub[i]>0)  S.rb2 = S.rb2+1
  
  #wald test
  wna = eta[i] * A * eta[i] * 2*n
  wrb = eta[i] * A^2 / B * eta[i] * 2*n
  if( wna<=qchisq(0.95, 1) )  W.na1 = W.na1+1
  if( wrb<=qchisq(0.95, 1) )  W.rb1 = W.rb1+1
  
  wt.na <- function(e){
    mod <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
               family = binomial, D)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    i.aa = ( m(a) + m(a+e+g) + m(a+e+d) + m(a+g+d) )/2
    i.ee = ( m(a+e+g) + m(a+e+d) )/2
    i.gg = ( m(a+e+g) + m(a+g+d) )/2
    i.dd = ( m(a+e+d) + m(a+g+d) )/2
    i.eg = m(a+e+g)/2
    i.ed = m(a+e+d)/2
    i.gd = m(a+g+d)/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                    i.gd, i.dd), nrow = 3, ncol = 3)
    A0 = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
    
    w.na = 2*n * (eta[i]-e)^2 * A0 - qchisq(0.95, 1)
    return(w.na)
  }
  wdna1 = uniroot(wt.na, c(eta[i]-q,eta[i]))
  wdna2 = uniroot(wt.na, c(eta[i],eta[i]+q))
  Wna.lb[i] = wdna1$root
  Wna.ub[i] = wdna2$root
  Wna.len[i] = Wna.ub[i] - Wna.lb[i]
  
  wt.rb <- function(e){
    mod <- glm(formula = cbind(Y,1-Y)~Z+G+offset(e*X), 
               family = binomial, D)
    a = mod$coefficients[1]
    g = mod$coefficients[2]
    d = mod$coefficients[3]
    i.aa = ( m(a) + m(a+e+g) + m(a+e+d) + m(a+g+d) )/2
    i.ee = ( m(a+e+g) + m(a+e+d) )/2
    i.gg = ( m(a+e+g) + m(a+g+d) )/2
    i.dd = ( m(a+e+d) + m(a+g+d) )/2
    i.eg = m(a+e+g)/2
    i.ed = m(a+e+d)/2
    i.gd = m(a+g+d)/2
    
    v.aa = i.aa + (p1[i]-f.inv(a)*f.inv(a+e+g)) +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.ee = i.ee
    v.gg = i.gg
    v.dd = i.dd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.ae = i.ee + (p1[i]-f.inv(a)*f.inv(a+e+g))/2 +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ag = i.gg + (p1[i]-f.inv(a)*f.inv(a+e+g))/2 +
      (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ad = i.dd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))
    v.eg = i.eg + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.ed = i.ed + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    v.gd = i.gd + (p2[i]-f.inv(a+e+d)*f.inv(a+g+d))/2
    
    i.ep = matrix(c(i.ee, i.eg, i.ed), ncol = 3)
    i.pp = matrix(c(i.aa, i.gg, i.dd, i.gg, i.gg, i.gd, i.dd, 
                    i.gd, i.dd), nrow = 3, ncol = 3)
    A = i.ee - i.ep %*% solve(i.pp) %*% t(i.ep)
    v.pe = matrix(c(v.ae, v.eg, v.ed), nrow = 3)
    v.pp = matrix(c(v.aa, v.ag, v.ad, v.ag, v.gg, v.gd, v.ad, 
                    v.gd, v.dd), nrow = 3, ncol = 3)
    B = v.ee - 2 * i.ep %*% solve(i.pp) %*% v.pe + 
      i.ep %*% solve(i.pp) %*% v.pp %*% solve(i.pp) %*% t(i.ep)
    
    s0 = n * ( mean(y21) - f.inv(a+e+g) +
                 mean(y12) - f.inv(a+e+d) )
    w.rb = 2*n * (eta[i]-e)^2 * A^2/B - qchisq(0.95, 1)
    return(w.rb)
  }
  wdrb1 = uniroot(wt.rb, c(eta[i]-q,eta[i]))
  wdrb2 = uniroot(wt.rb, c(eta[i],eta[i]+q))
  Wrb.lb[i] = wdrb1$root
  Wrb.ub[i] = wdrb2$root
  Wrb.len[i] = Wrb.ub[i] - Wrb.lb[i]
  
  if(Wna.lb[i]<0 && Wna.ub[i]>0)  W.na2 = W.na2+1
  if(Wrb.lb[i]<0 && Wrb.ub[i]>0)  W.rb2 = W.rb2+1
}
c(mean(alp), mean(eta), mean(ga), mean(del))
c(mean(alp.0), mean(ga.0), mean(del.0))
c(mean(p1), mean(p2))
II/r
VV/r
I.inv/r
I.inv_V_I.inv/r
matrix(c(var(alp), cov(alp,eta), cov(alp,ga), cov(alp,del),
         cov(alp,eta), var(eta), cov(eta,ga), cov(eta,del),
         cov(alp,ga), cov(eta,ga), var(ga), cov(ga,del),
         cov(alp,del), cov(eta,del), cov(ga,del), var(del)),
       nrow = 4, ncol = 4) * (2*n)
c(LR.na1, LR.rb1, LR.na2, LR.rb2)/r
c(mean(LRna.lb), mean(LRna.ub))
c(mean(LRrb.lb), mean(LRrb.ub))
c(mean(LRna.len), mean(LRrb.len))
c(S.na1, S.rb1, S.na2, S.rb2)/r
c(mean(Sna.lb), mean(Sna.ub))
c(mean(Srb.lb), mean(Srb.ub))
c(mean(Sna.len), mean(Srb.len))
c(W.na1, W.rb1, W.na2, W.rb2)/r
c(mean(Wna.lb), mean(Wna.ub))
c(mean(Wrb.lb), mean(Wrb.ub))
c(mean(Wna.len), mean(Wrb.len))




