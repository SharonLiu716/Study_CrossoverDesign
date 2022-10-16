#Python

class FinanceTerm():
  def __int__(self,name):
    self.name=name
  def pv_transform(FV,coupon_ratem,C0,i,t):
    PMT = -C0*coupon_rate
    C0+=PMT*(1/i - 1/(i*(1+i)**t))
    C0+=FV/((1+i)**t)
    return C0

FinanceTerm.pv_transform(100000,0.08,-100000,0.05,3)

#R
FinanceTerm<-setRefClass("FinanceTerm",fields = list(name="character",bundle="numeric"),
                         methods = list(
                           pv_transform=function(){
                             FV = bundle[1]
                             coupon_rate=bundle[2]
                             C0=bundle[3]
                             i=bundle[4]
                             t=bundle[5]
                             PMT= -C0*coupon_rate
                             C0=C0 + PMT*(1/i-1/(i*(1+i)**t))
                             C0=C0+FV/((1+i)**t)
                             return(C0)
                           }
                         )
                         )
bond_c1033<-FinanceTerm(name=c("FV","coupon_rate","C0","i","t"),
                         bundle=c(100000,0.08,-100000,0.05,3))
bond_c1033$pv_transform()