# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 18:31:21 2022

@author: 懿萱
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
sim_time=1000
alpha,eta,gamma,delta=1.0, 0.7, 0.3, 0.2
params = np.array([alpha,eta,gamma,delta])
covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
tm=np.exp(np.dot(params.transpose(),covariate))
#true mean
gamma_param=0.1#beta
seq_size=100
pi=(seq_size)/(seq_size*2)
I=pi*np.array([[np.dot(tm,covariate[i]*covariate[j]) for j in range(4)] for i in range(4)])
#data=pd.DataFrame(np.random.poisson(lam=tm, size=(seq_size, len('ABBA'))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])

score=[np.dot((data-tm).to_numpy(),covariate[i].transpose()) for i in range(len(covariate))]
V=np.array([[sum(score[i]*score[j]) for j in range(len(covariate))] for i in range(4)])/200

true_value=lin.inv(I).dot(V_hat).dot(lin.inv(I))

def logL(params,data):
    tao_, eta_, gamma_,delta_=params
    factorize=pd.DataFrame(np.array([[math.factorial(data.at[i,col]) for i in range(len(data))] for col in data.columns]).T.tolist(),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    logL=np.sum(tao_*data['Yi11']-np.exp(tao_)+(tao_+eta_+gamma_)*data['Yi12']-np.exp((tao_+eta_+gamma_))+(tao_+eta_+delta_)*data['Yi21']-np.exp((tao_+eta_+delta_))+(tao_+gamma_+delta_)*data['Yi22'] -np.exp((tao_+gamma_+delta_)))-np.log(factorize).sum().sum()       
    return -logL

def logL_mle(data):
    res = minimize(fun=lambda par, data: logL(par, data),x0=np.array([0.95, 0.62, 0.18, 0.07]).astype(float), args=(data,),method='BFGS')
    #tao_mle, eta_mle, gamma_mle, delta_mle = 
    return res.x,res.hess_inv

def matrix_AB(I_,V_):
    I_eta=I_[1, 1]
    I_etapsy=np.array([I_[1,0],I_[1, 2],I_[1, 3]])
    I_psy=np.array([[I_[0,0],I_[0, 2],I_[0, 3]],[I_[2,0],I_[2, 2],I_[2, 3]],[I_[3,0],I_[3, 2],I_[3, 3]]])
    V_eta=V_[1, 1]
    V_etapsy=np.array([V_[1,0],V_[1, 2],V_[1, 3]])
    V_psy=np.array([[V_[0,0],V_[0, 2],V_[0, 3]],[V_[2,0],V_[2, 2],V_[2, 3]],[V_[3,0],V_[3, 2],V_[3, 3]]])
    A=I_eta-I_etapsy.dot(lin.inv(I_psy)).dot(I_etapsy.reshape(3,1))
    B=V_eta-2*I_etapsy.dot(lin.inv(I_psy)).dot(V_etapsy.reshape(3,1))+I_etapsy.dot(lin.inv(I_psy)).dot(V_psy).dot(lin.inv(I_psy)).dot(I_etapsy.reshape(3,1))
    return A,B


def Estimate(mean_true,sample_size,data_type,cor_par,eta0):
    #mean_true=tm
    #sample_size=100000
    #cor_par=gamma_param
    #sample_size=int(sample_size/2)
    ''' 假設樣本之間為獨立生成的Poisson分配，求其MLE、I matrix、V matrix
    Input：
        -mean_true:用參數真值算出來的樣本平均數，用以生成獨立卜瓦松樣本y11、y21、y12、y22
        -sample_size：yij的樣本數(總樣本數為sample_sizs*4)
    Output：
        -estimate(dataframe):MLE
        -etimate_of_mean(dataframe)：mean of yij
        -I_hat：Fisher information matrix
        -V_hat：variance matrix 
    Validate：
        1. V=I?
        2.sample variance of MLE = diag(I**(-1))
    '''
    if data_type == 'ind': 
        data=pd.DataFrame(np.random.poisson(lam=mean_true, size=(seq_size, len('ABBA'))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    else:
        '''
        #np.random.seed(980716)
        data=Corr_Poisson(mean_true=mean_true,sample_size=sample_size, cor_par=cor_par)
        #data.corr()
        '''
        
        #np.random.seed(980716)
        Mu_cor = pd.DataFrame(columns = ['Mu11', 'Mu12', 'Mu21', 'Mu22'])
        for i in range(0,len(mean_true),2):            
            nu=np.random.gamma(1/cor_par,cor_par,sample_size)
            #print(i,nu)
            for (ix,param) in enumerate(mean_true[i:i + 2]):
                Mu_cor[Mu_cor.columns[ix+i]]=pd.DataFrame(np.multiply(np.array(param).repeat(sample_size), nu).T)
                #print(i,ix,param)

        data=pd.DataFrame(columns = ['Yi11', 'Yi12', 'Yi21', 'Yi22'])
        for (idx,mu) in enumerate(Mu_cor.columns):
            data[data.columns[idx]]=pd.DataFrame([np.random.poisson(p) for p in Mu_cor[mu]])
        
    #data.corr()
    #data.cov()
   
    '''MLE''' 
    #exp(alpha)
    #alpha_hat=np.log(np.mean(data['Yi11']))
    alpha_hat=np.log(data['Yi11'].sum())
    #exp(Eta)
    #eta_hat=0.5*(np.log(np.mean(data['Yi12']))+np.log(np.mean(data['Yi21']))-np.log(np.mean(data['Yi11']))-np.log(np.mean(data['Yi22'])))
    eta_hat=0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi21'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi22'].sum()))
    #exp(gamma)
    #gamma_hat=0.5*(np.log(np.mean(data['Yi12']))+np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi11']))-np.log(np.mean(data['Yi21'])))
    gamma_hat=0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi22'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi21'].sum()))
    #exp(delta)
    #delta_hat=0.5*(np.log(np.mean(data['Yi21']))+np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi11']))-np.log(np.mean(data['Yi12'])))
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
    '''
    I_hat=pi*np.array([[sum(sm), sm[1]+sm[2], sm[1]+sm[3], sum(sm[2:4])],
                   [sm[1]+sm[2],sum(sm[1:3]), sm[1], sm[2]],
                   [sm[1]+sm[3],sm[1], sm[1]+sm[3],sm[3]],
                   [sum(sm[2:4]),sm[2],sm[3],sum(sm[2:4])]
                   ])
   
    #Diagonal variance matrix
    Var=data.var(ddof=1)
    #Var=[var(data[res],sm[idx]) for (idx,res) in enumerate(data.columns)]
    Cov=[ data.Yi11.cov(data.Yi12), data.Yi21.cov(data.Yi22)]
    #Cov=[cov(data['Yi11'],data['Yi21'],sm[0],sm[1]),cov(data['Yi12'],data['Yi22'],sm[2],sm[3])]
    V_hat=pi*np.array([
                    [sum(Var)+2*sum(Cov), Var[1]+Cov[0]+Var[2]+Cov[1] , Var[1]+Cov[0]+Var[3]+Cov[1], sum(Var[2:4])+2*Cov[1]],
                    [Var[1]+Cov[0]+Var[2]+Cov[1], Var[1]+Var[2], Var[1]+Cov[1], Var[2]+Cov[1]],
                    [Var[1]+Cov[0]+Var[3]+Cov[1], Var[1]+Cov[1], Var[1]+Var[3], Var[3]+Cov[1]],
                    [sum(Var[2:4])+2*Cov[1], Var[2]+Cov[1],  Var[3]+Cov[1], sum(Var[2:4])+2*Cov[1]]
                    ])
    A_hat,B_hat=matrix_AB(I_hat,V_hat)
    
    eta_null=eta0
    gamma_null=np.log(np.mean(data['Yi21']))-np.log(np.mean(data['Yi11']))-eta_null
    delta_null=eta_null+np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi21']))
    LR_naive=2*(L(data=data,params=[alpha_hat,eta_hat,gamma_hat,delta_hat])-L(data=data,params=[alpha_hat,eta_null,gamma_null,delta_null])) 
    LR_robust=2*(A_hat/B_hat)*(L(data=data,params=[alpha_hat,eta_hat,gamma_hat,delta_hat])-L(data=data,params=[alpha_hat,eta_null,gamma_null,delta_null]))
    Wald_naive=(sample_size*4)*(eta_hat-eta_null)*A_hat*(eta_hat-eta_null)
    Wald_robust=(sample_size*4)*(eta_hat-eta_null)*(A_hat/B_hat)*(eta_hat-eta_null)'''
    
    return estimate, I_hat, V_hat#, LR_naive, LR_robust,Wald_naive,Wald_robust

np.random.seed(980716)
mle_ind,mle_cor=pd.DataFrame(),pd.DataFrame()
I_ind, V_ind,I_cor, V_cor = 0, 0, 0, 0
LR_ind_na,LR_ind_rb,LR_cor_na,LR_cor_rb=np.empty((0,1), float),np.empty((0,1), float),np.empty((0,1), float),np.empty((0,1), float)
Wald_ind_na,Wald_ind_rb,Wald_cor_na,Wald_cor_rb=np.empty((0,1), float),np.empty((0,1), float),np.empty((0,1), float),np.empty((0,1), float)
for i in range(sim_time):    
    #independent
    mle_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size,data_type='ind',cor_par=gamma_param,eta0=0)
    #mle_i, I_i, V_i, LR_na, LR_rb, Wald_na, Wald_rb=Estimate(mean_true=tm, sample_size=seq_size,data_type='ind',cor_par=gamma_param,eta0=0)
    mle_ind = mle_ind.append(mle_i,ignore_index=True)
    #mu_ind = mu_ind.append(mu_i,ignore_index=True)
    I_ind+=I_i
    V_ind+=V_i
    #LR_ind_na=np.append(LR_ind_na, LR_na)
    #LR_ind_rb=np.append(LR_ind_rb, LR_rb)
    #Wald_ind_na=np.append(Wald_ind_na, Wald_na)
    #Wald_ind_rb=np.append(Wald_ind_rb, Wald_rb)
    #correlated
    mle_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size,data_type='cor',cor_par=gamma_param,eta0=0)
    #mle_i, I_i, V_i, LR_na, LR_rb, Wald_na, Wald_rb=Estimate(mean_true=tm, sample_size=seq_size,data_type='cor',cor_par=gamma_param,eta0=0)
    mle_cor = mle_cor.append(mle_i,ignore_index=True)
    I_cor+=I_i
    V_cor+=V_i
    #LR_cor_na=np.append(LR_cor_na, LR_na)
    #LR_cor_rb=np.append(LR_cor_rb, LR_rb)
    #Wald_cor_na=np.append(Wald_cor_na, Wald_na)
    #Wald_cor_rb=np.append(Wald_cor_rb, Wald_rb)

I_ind=I_ind/sim_time
V_ind=V_ind/sim_time
I_cor=I_cor/sim_time
V_cor=V_cor/sim_time

print('Independent Data, seq_size = ',seq_size)
print('True Value of parameters\n alpha：%f, eta：%f, gamma：%f, delta：%f' %(1.0, 0.5, 0.2, 0.2))
print('MLE\n',mle_ind.mean())
print('Sample variance of estimates\n',2*seq_size*mle_ind.var(ddof=1))

print('Inverse of matrix I\n',lin.inv(I_ind))
print('Sample variance matrix\n',2*seq_size*mle_ind.cov(ddof=1))
print('inv(I)*V*inv(I)\n',lin.inv(I_ind).dot(V_ind).dot(lin.inv(I_ind)))

print('True Value of matrix I\n',I)
print('Estimate of matrix I\n',I_ind)
print('Estimate of matrix V\n',V_ind)

print('LR naive',np.mean(LR_ind_na))
print('LR robust',np.mean(LR_ind_rb))
print('Wald naive',np.mean(Wald_ind_na))
print('Wald robust',np.mean(Wald_ind_rb))

print('LR naive p-value',scipy.stats.t.sf(abs(np.mean(LR_ind_na)), df=1)*2)
print('LR robust p-value',scipy.stats.t.sf(abs(np.mean(LR_ind_rb)), df=1)*2)
print('Wald naive p-value',scipy.stats.t.sf(abs(np.mean(Wald_ind_na)), df=1)*2)
print('Wald robust p-value',scipy.stats.t.sf(abs(np.mean(Wald_ind_rb)), df=1)*2)




  
# find p-value for two-tailed test
#scipy.stats.t.sf(abs(np.mean(LR_ind_na)), df=1)*2

print('Correlated Data (alpha =10 ,beta = 0.1), seq_size = ',seq_size*2)
print('True Value of parameters\n alpha：%f, eta：%f, gamma：%f, delta：%f' %(1.0, 0.5, 0.2, 0.2))#1.0, 1.0, 1.0, 0.5
print('MLE\n',mle_cor.mean())
print('Sample variance of estimates\n',2*seq_size*mle_cor.var(ddof=1))

print('Inverse of matrix I\n',lin.inv(I_cor))
print('inv(I)*V*inv(I)\n',lin.inv(I_cor).dot(V_cor).dot(lin.inv(I_cor)))
print('Sample variance matrix\n',2*seq_size*mle_cor.cov(ddof=1))

print('Estimate of matrix I\n',I_cor)
print('Estimate of matrix V\n',V_cor)

print('LR naive',np.mean(LR_cor_na))
print('LR robust',np.mean(LR_cor_rb))
print('Wald naive',np.mean(Wald_cor_na))
print('Wald robust',np.mean(Wald_cor_rb))

print('LR naive p-value',scipy.stats.t.sf(abs(np.mean(LR_cor_na)), df=1)*2)
print('LR robust p-value',scipy.stats.t.sf(abs(np.mean(LR_cor_rb)), df=1)*2)
print('Wald naive p-value',scipy.stats.t.sf(abs(np.mean(Wald_cor_na)), df=1)*2)
print('Wald robust p-value',scipy.stats.t.sf(abs(np.mean(Wald_cor_rb)), df=1)*2)



def L(data,params):
    alpha, eta, gamma1, gamma2, delta1, delta2=params
    f1,f2=0,0
    for i in range(len(data)):
        f1=f1+alpha*data.at[i,'Yi11']+(alpha+eta+gamma)*data.at[i,'Yi21']
        f2=f2+(alpha+eta+delta)*data.at[i,'Yi12']+(alpha + gamma + delta)*data.at[i,'Yi22']
    
    return f1+f2
def Test(I_,V_):
    I_eta=I_[1, 1]
    I_etapsy=np.array([I_[1,0],I_[1, 2],I_[1, 3]])
    I_psy=np.array([[I_[0,0],I_[0, 2],I_[0, 3]],[I_[2,0],I_[2, 2],I_[2, 3]],[I_[3,0],I_[3, 2],I_[3, 3]]])
    V_eta=V_[1, 1]
    V_etapsy=np.array([V_[1,0],V_[1, 2],V_[1, 3]])
    V_psy=np.array([[V_[0,0],V_[0, 2],V_[0, 3]],[V_[2,0],V_[2, 2],V_[2, 3]],[V_[3,0],V_[3, 2],V_[3, 3]]])
    A=I_eta-I_etapsy.dot(lin.inv(I_psy)).dot(I_etapsy.reshape(3,1))
    B=V_eta-2*I_etapsy.dot(lin.inv(I_psy)).dot(V_etapsy.reshape(3,1))+I_etapsy.dot(lin.inv(I_psy)).dot(V_psy).dot(lin.inv(I_psy)).dot(I_etapsy.reshape(3,1))
    return A,B

A_ind,B_ind=Test(I_ind,V_ind)
A_cor,B_cor=Test(I_cor,V_cor)

'''
use random effect to generate data
--Yi~Poisson(mu_i*eta), eta~Gamma(alpha,beta)
--rho=np.sqrt(( (1+beta/mu_1)*(1+beta/mu_2) ))**(-1)
--fixed rho beta to find values of parameter 


rho,beta=0.7,3
def param_value(params):    
    alpha_,eta_,gamma_,delta_=params[0],params[1],params[2],params[3]
    mu_11,mu_21,mu_12,mu_22=np.exp(alpha_),np.exp(alpha_+eta_+gamma_),np.exp(alpha_+eta_+delta_),np.exp(alpha_+gamma_+delta_)
    return [np.sqrt(( (1+beta/mu_11)*(1+beta/mu_21) ))**(-1)-rho,np.sqrt(( (1+beta/mu_12)*(1+beta/mu_22) ))**(-1)-rho]
def param_value(params):    
    alpha_,eta_,gamma_,delta_=params
    #mu_11,mu_21,mu_12,mu_22=np.exp(alpha_),np.exp(alpha_+eta_+gamma_),np.exp(alpha_+eta_+delta_),np.exp(alpha_+gamma_+delta_)
    return [np.sqrt(( (1+beta/np.exp(alpha_))*(1+beta/np.exp(alpha_+eta_+gamma_)) ))**(-1)-rho,np.sqrt(( (1+beta/np.exp(alpha_+eta_+delta_))*(1+beta/np.exp(alpha_+gamma_+delta_)) ))**(-1)-rho]

solved=fsolve(param_value,[0.1,0.5,0.5,0.1])

# 二元二次方程組
alpha_ = Symbol('alpha')
eta_=  Symbol('eta')
gamma_=  Symbol('gamma')
delta_=  Symbol('delta')
solved_value=solve([sqrt(( (1+beta/exp(alpha_))*(1+beta/exp(alpha_+eta_+gamma_)) ))**(-1)-rho,sqrt(( (1+beta/exp(alpha_+eta_+delta_))*(1+beta/exp(alpha_+gamma_+delta_)) ))**(-1)-rho], [alpha_,eta_,gamma_,delta_])
print(solved_value)'''







