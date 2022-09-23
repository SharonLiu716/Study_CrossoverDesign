# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 20:33:37 2022

@author: 懿萱
"""
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import scipy
from scipy import optimize
from math import gamma

def f2(x):
    return gamma(1+ (2/x))-3*(gamma(1+ (1/x))**2)
a2 = scipy.optimize.brentq(f2, 0.1, 10000) 
l2=1/gamma(1+ (1/a2))

def f3(x):
    return gamma(1+ (2/x))-4*(gamma(1+ (1/x))**2)
a3 = scipy.optimize.brentq(f3, 0.1, 10000) 
l3=1/gamma(1+ (1/a3))

B=[0.5,1.0,1.0]


def Model(Beta,Size,params):
    M = pd.DataFrame(columns = ['Mu1', 'Mu2', 'Mu3'])
    for i in ['Mu1', 'Mu2', 'Mu3']:
        X1=np.random.uniform(low=params[0],high=params[1],size=Size)
        X2=np.random.uniform(low=params[0],high=params[1],size=Size)
        M[i]=Beta[0]+Beta[1]*X1+Beta[2]*X2
    return M
Mu=Model(B,10000,[1,3])
R=[1,3,5]

def Yij(V,Mu,R):
    df=pd.DataFrame(columns = ['Y1', 'Y2', 'Y3'])
    for (idx,mu) in enumerate(Mu.columns):    
        Mean=V[idx]*Mu[mu] #v*mu
        P=R[idx]/(R[idx]+Mean) #p
        df['Y'+str(idx+1)]=[ np.random.negative_binomial(R[idx], p) for p in P]
    Cor=[ pearsonr(df['Y1'], df['Y2'])[1],pearsonr(df['Y2'], df['Y3'])[1],pearsonr(df['Y1'], df['Y3'])[1]]
    
    return Cor#df,


V=pd.DataFrame()
V1=pd.DataFrame([np.random.uniform(low=0.5,high=1.5,size=10000) for i in range(3)]).transpose()
V2=pd.DataFrame([np.random.uniform(low=0,high=2,size=10000) for i in range(3)]).transpose()
V3=pd.DataFrame([np.random.gamma(shape=1, scale=1,size=10000) for i in range(3)]).transpose()
V4=pd.DataFrame([np.random.gamma(shape=1/2, scale=2,size=10000) for i in range(3)]).transpose()
V5=pd.DataFrame([np.random.gamma(shape=1/3, scale=3,size=10000) for i in range(3)]).transpose()
V6=pd.DataFrame([np.random.weibull(1,size=10000)*1 for i in range(3)]).transpose()
V7=pd.DataFrame([np.random.weibull(a2,size=10000)*l2 for i in range(3)]).transpose()
V8=pd.DataFrame([np.random.weibull(a3,size=10000)*l3 for i in range(3)]).transpose()
from statistics import variance 
variance(V1[0])
var,cor=[],[]
V=[V1,V2,V3,V4,V5,V6,V7,V8]
for (idx,vi) in enumerate(V):
    var.append([variance(vi[0]),variance(vi[1]),variance(vi[2])])
    cor.append(Yij(vi,Mu,R))

Var=pd.DataFrame(var,index=[*range(1,9)],columns = ['1','2','3'])
RHO=pd.DataFrame(cor,index=[*range(1,9)],columns = ['12','23','13'])

