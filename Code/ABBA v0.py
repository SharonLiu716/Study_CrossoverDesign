# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 09:58:46 2022

@author: cherl
"""
import sys
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
from autograd import grad
import autograd.numpy as np
from autograd import grad, jacobian, hessian
sim_time,seq_size,cor_param,cros_type=1000,100,0.1,'ABBA'
tao,eta,gamma,delta=1.0, 0.7, 0.3, 0.2
params = np.array([tao,eta,gamma,delta])
covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
mean_true=np.exp(np.dot(params.transpose(),covariate))
pi=(seq_size)/(seq_size*2)
np.array([[covariate[i]*covariate[j] for j in range(len(cros_type))] for i in range(len(cros_type))])
I=pi*np.array([[np.dot(mean_true,covariate[i]*covariate[j]) for j in range(len(cros_type))] for i in range(len(cros_type))])

data=pd.DataFrame(np.random.poisson(lam=mean_true, size=(seq_size, len(cros_type))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
data=np.random.poisson(lam=mean_true, size=(seq_size, len(cros_type))).astype(float)
len(data[:,0])
mean_true
#MLE
def logL(params):
    tao_, eta_, gamma_,delta_=params.astype(float)
    factorize=sum(sum(np.log(np.array([[math.factorial(data.at[i,col]) for i in range(len(data))] for col in data.columns]))))
    logLL=np.sum(tao_*data['Yi11']-np.exp(tao_)+(tao_+eta_+gamma_)*data['Yi12']-np.exp((tao_+eta_+gamma_))+(tao_+eta_+delta_)*data['Yi21']-np.exp((tao_+eta_+delta_))+(tao_+gamma_+delta_)*data['Yi22'] -np.exp((tao_+gamma_+delta_)))-factorize       
    return logLL

def logL(params):
    tao_, eta_, gamma_,delta_=params.astype(float)
    factorize=sum(sum(np.log(np.array([[math.factorial(data[i,j]) for i in range(len(data))] for j in range(len(cros_type))]))))
    #logLL=np.sum(tao_*data[:,0]-np.exp(tao_)+(tao_+eta_+gamma_)*data[:,1]-np.exp((tao_+eta_+gamma_))+(tao_+eta_+delta_)*data[:,2]-np.exp((tao_+eta_+delta_))+(tao_+gamma_+delta_)*data[:,3] -np.exp((tao_+gamma_+delta_)))-factorize           
    #return logLL
    logLL=np.sum(tao_*data[:,0].astype(float)-np.exp(tao_)+(tao_+eta_+gamma_)*data[:,1].astype(float)-np.exp((tao_+eta_+gamma_))+(tao_+eta_+delta_)*data[:,2].astype(float)-np.exp((tao_+eta_+delta_))+(tao_+gamma_+delta_)*data[:,3].astype(float) -np.exp((tao_+gamma_+delta_)))-factorize       
    return logLL
def logL(params,data):
    tao_, eta_, gamma_,delta_=params
    factorize=pd.DataFrame(np.array([[math.factorial(data.at[i,col]) for i in range(len(data))] for col in data.columns]).T.tolist(),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    logL=np.sum(tao_*data['Yi11']-np.exp(tao_)+(tao_+eta_+gamma_)*data['Yi12']-np.exp((tao_+eta_+gamma_))+(tao_+eta_+delta_)*data['Yi21']-np.exp((tao_+eta_+delta_))+(tao_+gamma_+delta_)*data['Yi22'] -np.exp((tao_+gamma_+delta_)))-np.log(factorize).sum().sum()       
    return -logL
res = minimize(fun=lambda par: -logL(par),x0=np.array([0.95, 0.65, 0.25, 0.15], dtype=float))
mle=list(res.x)
hessian_ = hessian(logL)
VCM = lin.inv(-hessian_(res.x))
VCM*200
lin.inv(I)
-hessian_(res.x)/200
I
score=[np.dot((data-mean_true),covariate[i].transpose()) for i in range(len(covariate))]
V=np.array([[sum(score[i]*score[j]) for j in range(len(covariate))] for i in range(4)])/200
def Estimate(cors_type,mean_true,sample_size,data_type,cor_par,eta0):

    if data_type == 'ind': 
        data=pd.DataFrame(np.random.poisson(lam=mean_true, size=(seq_size, len('ABBA'))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    else:   
        
        Mu_cor = pd.DataFrame(columns = ['Mu11', 'Mu12', 'Mu21', 'Mu22'])
        for i in range(0,len(mean_true),2):            
            nu=np.random.gamma(1/cor_par,cor_par,sample_size)
            for (ix,param) in enumerate(mean_true[i:i + 2]):
                Mu_cor[Mu_cor.columns[ix+i]]=pd.DataFrame(np.multiply(np.array(param).repeat(sample_size), nu).T)
                

        data=pd.DataFrame(columns = ['Yi11', 'Yi12', 'Yi21', 'Yi22'])
        for (idx,mu) in enumerate(Mu_cor.columns):
            data[data.columns[idx]]=pd.DataFrame([np.random.poisson(p) for p in Mu_cor[mu]])
   
    '''MLE''' 
    alpha_hat=np.log(data['Yi11'].mean())
    eta_hat=0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi21'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi22'].sum()))
    gamma_hat=0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi22'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi21'].sum()))
    delta_hat=0.5*(np.log(data['Yi21'].sum())+np.log(data['Yi22'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi12'].sum()))
    # MLE
    estimate= pd.DataFrame({'alpha_hat': alpha_hat, 
                            'eta_hat': eta_hat, 
                            'gamma_hat':gamma_hat, 
                            'delta_hat': delta_hat},index=[0])
    
    '''Mean'''
    mle = np.array([alpha_hat,eta_hat,gamma_hat,delta_hat])
    covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
    etimate_of_mean=np.exp(np.dot(mle.transpose(),covariate))
    
    '''Estimate matrix I''' 
    #sm=list(etimate_of_mean.loc[0,:])#data.mean()
    pi=sample_size/(sample_size*2)
    I_hat=pi*np.array([[np.dot(etimate_of_mean,covariate[i]*covariate[j]) for j in range(4)] for i in range(4)])
    score=[np.dot((data-etimate_of_mean).to_numpy(),covariate[i].transpose()) for i in range(len(covariate))]
    V_hat=np.array([[sum(score[i]*score[j]) for j in range(len(covariate))] for i in range(4)])/(sample_size*2)
    return estimate, I_hat, V_hat

np.random.seed(980716)
mle_ind,mle_cor=pd.DataFrame(),pd.DataFrame()
I_ind, V_ind,I_cor, V_cor = 0, 0, 0, 0
for i in range(sim_time):    
    #independent
    mle_i, I_i, V_i=Estimate(cors_type='ABBA',mean_true=tm,sample_size=seq_size,data_type='ind',cor_par=gamma_param,eta0=0)
    #mle_i, I_i, V_i, LR_na, LR_rb, Wald_na, Wald_rb=Estimate(mean_true=tm, sample_size=seq_size,data_type='ind',cor_par=gamma_param,eta0=0)
    mle_ind = mle_ind.append(mle_i,ignore_index=True)
    #mu_ind = mu_ind.append(mu_i,ignore_index=True)
    I_ind+=I_i
    V_ind+=V_i
    #correlated
    mle_i, I_i, V_i=Estimate(cors_type='ABBA',mean_true=tm,sample_size=seq_size,data_type='cor',cor_par=gamma_param,eta0=0)
    mle_cor = mle_cor.append(mle_i,ignore_index=True)
    I_cor+=I_i
    V_cor+=V_i
    

I_ind=I_ind/sim_time
V_ind=V_ind/sim_time
I_cor=I_cor/sim_time
V_cor=V_cor/sim_time
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

for seq in seq_size: 
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