# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 09:58:46 2022

@author: cherl
"""
import sys
sys.executable
sys.path
sys.prefix



import scipy
import math
from math import factorial
import numpy as np
import pandas as pd
import scipy.stats
from scipy.optimize import minimize, rosen, rosen_der
from scipy.stats import poisson
import numpy.linalg as lin
from scipy.optimize import fsolve
from sympy import *
from torch.autograd import Variable, grad
import autograd.numpy as ag
from autograd import grad, jacobian, hessian
sim_time,seq_size,cor_param,cros_type=1000,100,0.1,'ABBA'
tao,eta,gamma,delta=1.0, 0.67, 0.23, 0.12
params = np.array([tao,eta,gamma,delta])
covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
mean_true=np.exp(np.dot(params.transpose(),covariate))
pi=(seq_size)/(seq_size*2)
I=pi*np.array([[np.dot(mean_true,covariate[i]*covariate[j]) for j in range(len(cros_type))] for i in range(len(cros_type))])
df=pd.DataFrame(np.array([poisson.rvs(p, size=seq_size) for p in mean_true]).T.tolist(),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
df=pd.DataFrame(np.random.poisson(lam=mean_true, size=(seq_size, len(cros_type))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
df.mean()
df.corr()
mean_true
#MLE
def logL(params,data):
    tao_, eta_, gamma_,delta_=params
    factorize=pd.DataFrame(np.array([[math.factorial(data.at[i,col]) for i in range(len(data))] for col in data.columns]).T.tolist(),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    logL=np.sum(tao_*data['Yi11']-np.exp(tao_)+(tao_+eta_+gamma_)*data['Yi12']-np.exp((tao_+eta_+gamma_))+(tao_+eta_+delta_)*data['Yi21']-np.exp((tao_+eta_+delta_))+(tao_+gamma_+delta_)*data['Yi22'] -np.exp((tao_+gamma_+delta_)))-np.log(factorize).sum().sum()       
    return -logL

def logL_mle(data):
    res = minimize(fun=lambda par, data: logL(par, data),x0=np.array([0.95, 0.62, 0.18, 0.07]).astype(float), args=(data,),method='BFGS')
    #tao_mle, eta_mle, gamma_mle, delta_mle = 
    return res.x,res.hess_inv

mle=[]
sigma=0
for i in range(sim_time):
    df=pd.DataFrame(np.random.poisson(lam=mean_true, size=(seq_size, len(cros_type))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    mle_i,sigma_i=logL_mle(df)
    mle.append(list(mle_i))
    sigma+=sigma_i


mle_df=pd.DataFrame (mle, columns = ['tao_hat', 'eta_hat','gamma_hat','delta_hat'])
mle_df.mean()
mle_df.cov()*(seq_size*2)
sigma
lin.inv(I)
lin.inv(sigma)
I

'''
用來驗證logL有沒有算錯 但optim的值不正確!
def logL(params,data):
    tao_, eta_, gamma_,delta_=params
    params = np.array([tao_, eta_, gamma_,delta_])
    covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
    logmean=np.dot(params.transpose(),covariate)
    factorize=pd.DataFrame(np.array([[math.factorial(data.at[i,col]) for i in range(len(data))] for col in data.columns]).T.tolist(),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    logL=np.sum((np.log(mean_true)*data-mean_true-np.log(factorize))).sum()
    
    return -logL'''


score=[np.dot((data-mean_true).to_numpy(),covariate[i].transpose()) for i in range(len(covariate))]
V=np.array([[sum(score[i]*score[j]) for j in range(len(covariate))] for i in range(4)])/200