# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 12:10:21 2022

@author: a7086
"""

import scipy
import math
from math import factorial
import numpy as np
import pandas as pd
from scipy.optimize import minimize 
import numpy.linalg as lin
from scipy.stats import poisson

'''ABBA'''
alpha,eta,gamma,delta=1.0, 0.5, 0.2, 0.2
params = np.array([ [alpha], [eta], [gamma], [delta]])
covariate=[[1,0,0,0],[1,1,1,0],[1,1,0,1],[1,0,1,1]]
#true mean
tm=np.exp(np.array([ xijk  for xijk in covariate]).dot(params)).reshape(1,4).tolist()[0]

seq_size=100000
pi=(seq_size)/(seq_size*4)
I=pi*np.array([
               [sum(tm), tm[1]+tm[2], tm[1]+tm[3], sum(tm[2:4])],
               [tm[1]+tm[2],sum(tm[1:3]), tm[1], tm[2]],
               [tm[1]+tm[3],tm[1], tm[1]+tm[3],tm[3]],
               [sum(tm[2:4]),tm[2],tm[3],sum(tm[2:4])]
               ])

np.random.seed(980716)
mean_true=tm
cor_par=[]
def Corr_Poisson(mean_true,sample_size,cor_par):
    '''
    Algorithm(for a sequence data):
        -given true mean 1 、true mean 2、correlated_mean、ratio k
        -generated u1,u2,u3~iid U(0,1)
        -X1=ppf(u1,tm1-correlated_mean)、X2=ppf(u2,tm2-k*correlated_mean)、X12=ppf(u3,correlated_mean)、X12_=ppf(u3,k*correlated_mean)
        -Y1=X1+X12、Y2=X2+X12_
    mean_true:mean of Yi
    sample_size:seq_size(number of x1+x2)
    cor_par:correlated mean    
    '''    
    dict_parmas={'K':[mean_true[1]/mean_true[0],mean_true[3]/mean_true[2]],
                 'Lam12':[cor_par[0],cor_par[1]],
                 'Lam1':[mean_true[0]-cor_par[0],mean_true[2]-cor_par[1]],
                 'Lam2':[mean_true[1]-cor_par[0]*mean_true[1]/mean_true[0],mean_true[3]-cor_par[1]*mean_true[3]/mean_true[2]],
                 'Lam12_':[cor_par[0]*mean_true[1]/mean_true[0],cor_par[1]*mean_true[3]/mean_true[2]]}
    df=pd.DataFrame()
    for i in range(0,int(len(tm)/2)):
        u1,u2,u3=np.random.uniform(0,1,size=sample_size),np.random.uniform(0,1,size=sample_size),np.random.uniform(0,1,size=sample_size)
        Y1,Y2,Y12,Y12_=scipy.stats.poisson.ppf(u1, dict_parmas['Lam1'][i]),poisson.ppf(u2, dict_parmas['Lam2'][i]),poisson.ppf(u3,dict_parmas['Lam12'][i]),poisson.ppf(u3, dict_parmas['Lam12_'][i])
        X1,X2=pd.Series(Y1+Y12),pd.Series(Y2+Y12_)
        df=pd.concat([df, X1,X2], axis=1)
    
    df.columns = ['Yi11', 'Yi21', 'Yi12', 'Yi22']
    return df

temp=Corr_Poisson(mean_true=tm,sample_size=seq_size, cor_par=[0.6,1.2])
temp.corr()
sample_size=200000
sample_size=int(sample_size/2)

u1,u2,u3=np.random.uniform(0,1,size=sample_size),np.random.uniform(0,1,size=sample_size),np.random.uniform(0,1,size=sample_size)
np.random.seed(980716)
data=pd.DataFrame()
for i in range(0,int(len(tm)/2)):
    u1,u2,u3=np.random.uniform(0,1,size=sample_size),np.random.uniform(0,1,size=sample_size),np.random.uniform(0,1,size=sample_size)
    Y1,Y2,Y12,Y12_=scipy.stats.poisson.ppf(u1, temp['Lam1'][i]),poisson.ppf(u2, temp['Lam2'][i]),poisson.ppf(u3,temp['Lam12'][i]),poisson.ppf(u3, temp['Lam12_'][i])
    X1,X2=pd.Series(Y1+Y12),pd.Series(Y2+Y12_)
    data=pd.concat([data, X1,X2], axis=1)

data.corr()
data.mean()


k=tm[1]/tm[0]
lam12=0.6
#correlation
(np.mean(poisson.ppf(u3,lam12)*poisson.ppf(u3, k*lam12))-k*lam12**(2))/(tm[0]*np.sqrt(k))
#lam12=np.float(solve(x/(((tm[0]+x)*(tm[1]+x))**(0.5))-0.1, x)[0])
mu1,mu2,mu12,mu12_=tm[0]-lam12,tm[1]-k*lam12,lam12,k*lam12
Y1,Y2,Y12,Y12_=scipy.stats.poisson.ppf(u1, mu1),poisson.ppf(u2, mu2),poisson.ppf(u3,mu12),poisson.ppf(u3, mu12_)
X1,X2=Y1+Y12,Y2+Y12_


lam12=1.2
#lam12=np.float(solve(x/(((tm[2]+x)*(tm[3]+x))**(0.5))-0.1, x)[0])
k=tm[3]/tm[2]
#correlation
(np.mean(poisson.ppf(u3,lam12)*poisson.ppf(u3, k*lam12))-k*lam12**(2))/(tm[2]*np.sqrt(k))
mu1,mu2,mu12,mu12_=tm[2]-lam12,tm[3]-k*lam12,lam12,k*lam12
Y1,Y2,Y12,Y12_=scipy.stats.poisson.ppf(u1, mu1),poisson.ppf(u2, mu2),poisson.ppf(u3,mu12),poisson.ppf(u3, mu12_)
X1,X2=Y1+Y12,Y2+Y12_
np.corrcoef(X1,X2)


'''ABBBAA'''
import scipy
import math
from math import factorial
import numpy as np
import pandas as pd
from scipy.optimize import minimize 
from scipy.stats import poisson
import numpy.linalg as lin


alpha, eta, gamma1, gamma2, delta = 1.0, 0.5, 0.2, 0.2, 0.2
params = np.array([ [alpha], [eta], [gamma1], [gamma2], [delta]])
covariate=[[1,0,0,0,0],[1,1,1,0,0],[1,1,0,1,0],[1,1,0,0,1],[1,0,1,0,1],[1,0,0,1,1]]
#true mean
tm=np.exp(np.array([ xijk  for xijk in covariate]).dot(params)).reshape(1,6).tolist()[0]
gamma_param=[0.6,1.2]
sample_size=100000

np.random.seed(980716)
u1,u2,u3=np.random.uniform(0,1,size=sample_size),np.random.uniform(0,1,size=sample_size),np.random.uniform(0,1,size=sample_size)

k1121,k1131,k2131=tm[1]/tm[0],tm[2]/tm[0],tm[2]/tm[1]
lam11,lam12,lam13=0.6,0.6,1.09
#correlation
np.random.seed(980716)
print('rho for Y11 & Y21:',(np.mean(poisson.ppf(u3,lam11)*poisson.ppf(u3, k1121*lam11))-k1121*lam11**(2))/(tm[0]*np.sqrt(k1121)))
print('rho for Y11 & Y31:',(np.mean(poisson.ppf(u3,lam12)*poisson.ppf(u3, k1131*lam12))-k1131*lam12**(2))/(tm[0]*np.sqrt(k1131)))
print('rho for Y21 & Y31:',(np.mean(poisson.ppf(u3,lam13)*poisson.ppf(u3, k2131*lam13))-k2131*lam13**(2))/(tm[1]*np.sqrt(k2131)))

mu1_1121,mu2,mu12,mu12_=tm[0]-lam12,tm[1]-k*lam12,lam12,k*lam12
Y1,Y2,Y12,Y12_=scipy.stats.poisson.ppf(u1, mu1),poisson.ppf(u2, mu2),poisson.ppf(u3,mu12),poisson.ppf(u3, mu12_)
X1,X2=Y1+Y12,Y2+Y12_




