Matrix.I<-function(cros.type,params,x.mat){
  num.seq<-if (nchar(cros.type)==9) 3 else 2
  mat.I<- matrix(0, nrow = length(params), ncol = length(params))
  Mean<-exp(params%*%x.mat)
  for (i in 1:length(params) ){ mat.I[i,]<-Mean%*%(t(x.mat)*x.mat[i,])/num.seq}
  return(mat.I)
}

Matrix.AB<-function(Mat.I,Mat.V){
  Mat.I<-I.ind.i
  diag.I<-as.matrix(diag(Mat.I),byrow = TRUE)
  diag.V<-as.matrix(diag(Mat.V),byrow = TRUE)
  #I.theta.theta
  I.tt<-diag.I[2];V.tt<-diag.V[2];I.pp<-Mat.I[-c(2),-c(2)];V.pp<-Mat.V[-c(2),-c(2)]
  I.tp<-as.matrix(Mat.I[2,]);I.tp<-as.matrix(I.tp[-c(2),]);V.tp<-as.matrix(Mat.V[2,]);V.tp<-as.matrix(V.tp[-c(2),])
  #mat.A
  Mat.A<-I.tt-t(I.tp)%*%solve(I.pp)%*%I.tp
  Mat.B<-V.tt-2*t(I.tp)%*%solve(I.pp)%*%V.tp+t(I.tp)%*%solve(I.pp)%*%V.pp%*%solve(I.pp)%*%I.tp
  Mat.AB<- list("Mat.A" = Mat.A, "Mat.B" =Mat.B)
  return(Mat.AB)
}

t=1.2;e=0.3;g=1.0;d=0.2
mle<-c(t,e,g,d)
Matrix.I('ABBA',mle,xmat)
i.tt<-( exp(t) + exp(t+e+g) + exp(t+e+d) + exp(t+g+d) )/2
i.ee<-( exp(t+e+g) + exp(t+e+d) )/2 
i.gg<-( exp(t+e+g) + exp(t+g+d) )/2 
i.dd<- ( exp(t+e+d)+exp(t+g+d) )/2
i.eg<- ( exp(t+e+g) )/2
i.ed<- ( exp(t+e+d) )/2
i.gd<- ( exp(t+g+d) )/2
I<-matrix(c(i.tt,i.ee,i.gg,i.dd,
            i.ee,i.ee,i.eg,i.ed,
            i.gg,i.eg,i.gg,i.gd,
            i.dd,i.ed,i.gd,i.dd),nrow=4, ncol=4, byrow = TRUE)
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


