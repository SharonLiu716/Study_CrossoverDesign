# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.optimize import fsolve

import scipy
from scipy import optimize
from math import gamma
alpha, eta1, eta2, gamma1, gamma2, delta1, delta2 = 1.0, 0.32, 0.37, 1.21, 1.27, 0.22, 0.28
params=[np.exp(alpha), np.exp(alpha + eta1+ gamma1), np.exp(alpha + eta2 + gamma2),
        np.exp(alpha + eta1 + delta1), np.exp(alpha + eta2 + gamma1 + delta1), np.exp(alpha + gamma2 + delta1),
        np.exp(alpha + eta2 +delta2), np.exp(alpha + gamma1 + delta2), np.exp(alpha + eta1 + gamma2 + delta2)]

mean_true=[ params[i:i + 3] for i in range(0, len(params), 3) ]
sample_size=100000
Mu_cor = pd.DataFrame(columns = ['Mu11', 'Mu21', 'Mu31','Mu12', 'Mu22', 'Mu32','Mu13', 'Mu23', 'Mu33'])

for i in range(0,len(params),3):
    nu=np.random.gamma(1,1,sample_size)
    for (ix,param) in enumerate(params[i:i + 3]):
        Mu_cor[Mu_cor.columns[ix+i]]=pd.DataFrame(np.multiply(np.array(param).repeat(sample_size), nu).T)

data=pd.DataFrame(columns = ['Y11', 'Y21', 'Y31','Y12', 'Y22', 'Y32','Y13', 'Y23', 'Y33'])
for (idx,mu) in enumerate(Mu_cor.columns):
    data[data.columns[idx]]=pd.DataFrame([np.random.poisson(p) for p in Mu_cor[mu]])
Mu_cor=np.array([np.multiply(np.array(mean_true[ix][i]).repeat(sample_size),np.random.gamma(1,1,sample_size)) for i in range(len(mean_true[0])) for ix in range(len(mean_true))]).reshape(9,sample_size)
Mu_cor=pd.DataFrame(Mu_cor.transpose(),columns = ['Mu11', 'Mu21', 'Mu31','Mu12', 'Mu22', 'Mu32','Mu13', 'Mu23', 'Mu33'])
data=pd.DataFrame(columns = ['Y11', 'Y21', 'Y31','Y12', 'Y22', 'Y32','Y13', 'Y23', 'Y33'])
for (idx,mu) in enumerate(Mu_cor.columns):
    data[data.columns[idx]]=pd.DataFrame([np.random.poisson(p) for p in Mu_cor[mu]])
data.mean()
data.corr() 
data.cov()
Mu = pd.DataFrame(columns = ['Mu11', 'Mu21', 'Mu31','Mu12', 'Mu22', 'Mu32','Mu13', 'Mu23', 'Mu33'])
df=pd.DataFrame(columns = ['Y11', 'Y21', 'Y31','Y12', 'Y22', 'Y32','Y13', 'Y23', 'Y33'])
for (idx,param) in enumerate(params): 
    #temp=param+np.random.gamma(shape=1, scale=1,size=25).reshape(25,1) 
    Mu[Mu.columns[idx]]=pd.DataFrame([param*np.random.gamma(shape=1, scale=1) ]).transpose()
    df[df.columns[idx]]=pd.DataFrame([np.random.poisson(Mu[Mu.columns[idx]],size=30)]).transpose()

df.corr() 
df.cov() 
df.mean()

alpha,eta,gamma,delta=1.0, 0.3, 1.2, 0.2
tm=[np.exp(alpha),np.exp(alpha+eta+gamma),np.exp(alpha+eta+delta),np.exp(alpha+gamma+delta)]
  
rho=0.75

def solve_randomeffect(params):
    rho,theta1,theta2,beta=0.7,1.6487,4.4817,params
    return (1+beta/theta1)*(1+beta/theta2)-rho**(-2)
solved=fsolve(solve_randomeffect,[0.1])

lam,lam_=np.exp(alpha)+solved,np.exp(alpha+eta+gamma)+solved
y11,y21=np.random.poisson(lam,100000),np.random.poisson(lam_,100000)
np.corr(y11,y21)
lam*lam_
100000*np.corrcoef(y11,y21)
