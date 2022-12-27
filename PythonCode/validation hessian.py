# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 18:30:09 2022

@author: cherl
"""

import scipy
import math
from math import factorial
import numpy as np
import pandas as pd
import scipy.stats
from scipy.optimize import minimize 
from scipy.stats import poisson
import numpy.linalg as lin
from scipy.optimize import fsolve
from sympy import *
import autograd.numpy as ag
from autograd import grad, jacobian, hessian
sim_time=1000
alpha,eta,gamma,delta=1.0, 0.67, 0.23, 0.12
params = np.array([alpha,eta,gamma,delta])
covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
tm=np.exp(np.dot(params.transpose(),covariate))
#true mean
gamma_param=0.1#beta
seq_size=100
pi=(seq_size)/(seq_size*2)
I=pi*np.array([[np.dot(tm,covariate[i]*covariate[j]) for j in range(4)] for i in range(4)])
data=pd.DataFrame(np.array([poisson.rvs(p, size=seq_size) for p in tm]).T.tolist(),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
score=[np.dot((data-tm).to_numpy(),covariate[i].transpose()) for i in range(len(covariate))]
V=np.array([[sum(score[i]*score[j]) for j in range(len(covariate))] for i in range(4)])/200

true_value=lin.inv(I).dot(V_hat).dot(lin.inv(I))

data.corr()
def L(params):
    alpha,eta,gamma,delta=params
    f1,f2=0,0
    for i in range(len(data)):
        f1=f1+alpha*data.at[i,'Yi11']-np.exp(alpha)+(alpha+eta+gamma)*data.at[i,'Yi12']-np.exp(alpha+eta+gamma)# - np.log(factorial(data.at[i,'Yi11']))-np.log(factorial(data.at[i,'Yi21']))
        f2=f2+(alpha+eta+delta)*data.at[i,'Yi21']-np.exp(alpha+eta+delta)+(alpha+gamma+delta)*data.at[i,'Yi22']-np.exp(alpha+gamma+delta)#-np.log(factorial(data.at[i,'Yi12']))-np.log(factorial(data.at[i,'Yi22']))
    
    return f1+f2
initial_guess = np.array([0.93,0.60,0.2,0.1])
solve = scipy.optimize.minimize (lambda params: -L(params), initial_guess)
lin.inv(solve.hess_inv)
MLE_opt=solve.x
hessian_ = hessian(L)
I_ = hessian_(np.array([alpha,eta,gamma,delta]))
