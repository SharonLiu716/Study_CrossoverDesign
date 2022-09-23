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
from scipy.optimize import minimize 
import numpy.linalg as lin

alpha,eta,gamma,delta=1.0, 0.3, 1.2, 0.2
true_mean=[np.exp(alpha),np.exp(alpha+eta+gamma),np.exp(alpha+eta+delta),np.exp(alpha+gamma+delta)]
seq_size=200
p1,p2=0.5,0.5
I=np.array([[p1*(true_mean[0]+true_mean[1])+p2*(true_mean[2]+true_mean[3]), p1*true_mean[1]+p2*true_mean[2] ,p1*(true_mean[1])+p2*true_mean[3] , p2*(true_mean[2]+true_mean[3])],
            [p1*true_mean[1]+p2*true_mean[2],p1*true_mean[1]+p2*true_mean[2], p1*true_mean[1], p2*true_mean[2] ],
            [p1*(true_mean[1])+p2*true_mean[3],p1*true_mean[1], p1*true_mean[1]+p2*true_mean[3],p2*true_mean[3]],
            [p2*(true_mean[2]+true_mean[3]),p2*true_mean[2],p2*true_mean[3],p2*(true_mean[2]+true_mean[3])]
            ])

def var(x,mean_x): # * removed
    return sum((x - mean_x)**2)/len(x)
def cov(x, y,mean_x,mean_y): 
    return sum((x - mean_x) * (y - mean_y)) / len(x) 

def Estimate(mean_true,sample_size):
    
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
    #generate data
    data = {'Yi11': np.random.poisson(lam=mean_true[0], size=sample_size),
            'Yi21': np.random.poisson(lam=mean_true[1], size=sample_size),
            'Yi12': np.random.poisson(lam=mean_true[2], size=sample_size),  
            'Yi22': np.random.poisson(lam=mean_true[3], size=sample_size)}
    data = pd.DataFrame(data)
    '''MLE''' 
    #exp(alpha)
    alpha_hat=np.mean(data['Yi11'])
    #exp(Eta)
    eta_hat=np.sqrt(np.mean(data['Yi21'])*np.mean(data['Yi12'])/(np.mean(data['Yi11'])*np.mean(data['Yi22'])))
    #exp(gamma)
    gamma_hat=np.mean(data['Yi21'])/np.mean(data['Yi11'])/eta_hat
    #exp(delta)
    delta_hat=np.mean(data['Yi12'])/np.mean(data['Yi11'])/eta_hat
    
    # MLE
    estimate= pd.DataFrame({'alpha_hat': np.log(alpha_hat), 
                            'eta_hat': np.log(eta_hat), 
                            'gamma_hat':np.log(gamma_hat), 
                            'delta_hat': np.log(delta_hat)},index=[0])
    
    '''Mean'''
    etimate_of_mean = pd.DataFrame({'Mean_11': alpha_hat, 
                       'Mean_21': alpha_hat*eta_hat*gamma_hat, 
                       'Mean_12': alpha_hat*eta_hat*delta_hat, 
                       'Mean_22': alpha_hat*gamma_hat*delta_hat},index=[0])
    
    '''Estimate matrix I''' 
    p1,p2=0.5,0.5
    I_hat=np.matrix([[p1*sum(data.mean()), p1*data['Yi21'].mean()+p2*data['Yi12'].mean() , p1*data['Yi21'].mean()+p2*data['Yi22'].mean() , p2*(data['Yi12'].mean()+data['Yi22'].mean())],
                   [p1*data['Yi21'].mean()+p2*data['Yi12'].mean(),p1*data['Yi21'].mean()+p2*data['Yi12'].mean(), p1*data['Yi21'].mean(), p2*data['Yi12'].mean()],
                   [p1*data['Yi21'].mean()+p2*data['Yi22'].mean(), p1*data['Yi21'].mean(), p1*data['Yi21'].mean()+p2*data['Yi22'].mean(), p2*data['Yi22'].mean()],
                   [p2*(data['Yi12'].mean()+data['Yi22'].mean()), p2*data['Yi12'].mean(), p2*data['Yi22'].mean(), p2*(data['Yi12'].mean()+data['Yi22'].mean())]])
   
    '''Diagonal variance matrix'''
    mean=list(etimate_of_mean.loc[0,:])
    '''
    Var=[var(data[res],mean[idx]) for (idx,res) in enumerate(data.columns)]
    Cov=[cov(data['Yi11'],data['Yi21'],mean[0],mean[1]),cov(data['Yi12'],data['Yi22'],mean[2],mean[3])]
    V_hat=np.matrix([[p1*sum(Var), p1*( Var[1]+Cov[0] ) + p2*( Var[2]+Cov[1] ) , p1*( Var[1]+Cov[0] )+p2*( Var[3]+Cov[1] ) , p2*(Var[2]+Var[3])],
                    [p1*( Var[1]+Cov[0] ) + p2*( Var[2]+Cov[1] ), p1*Var[1]+p2*Var[2], p1*Var[1] + p2*Cov[1], p2*(Var[2]+Cov[1])],
                    [p1*( Var[1]+Cov[0] )+p2*( Var[3]+Cov[1] ), p1*Var[1] + p2*Cov[1], p1*Var[1] + p2*Var[3], p2*(Var[3]+Cov[1])],
                    [p2*(Var[2]+Var[3]), p2*(Var[2]+Cov[1]), p2*(Var[3]+Cov[1]), p2*(Var[2]+Var[3])]])'''
    '''Diagonal variance matrix'''
    Var=data.var(ddof=1)
    #Var=[var(data[res],sm[idx]) for (idx,res) in enumerate(data.columns)]
    Cov=[ data.Yi11.cov(data.Yi21), data.Yi12.cov(data.Yi22)]
    #Cov=[cov(data['Yi11'],data['Yi21'],sm[0],sm[1]),cov(data['Yi12'],data['Yi22'],sm[2],sm[3])]
    V_hat=p1*np.array([
                    [sum(Var)+2*sum(Cov), Var[1]+Cov[0]+Var[2]+Cov[1] , Var[1]+Cov[0]+Var[3]+Cov[1], sum(Var[2:4])+2*Cov[1]],
                    [Var[1]+Cov[0]+Var[2]+Cov[1], Var[1]+Var[2], Var[1]+Cov[1], Var[2]+Cov[1]],
                    [Var[1]+Cov[0]+Var[3]+Cov[1], Var[1]+Cov[1], Var[1]+Var[3], Var[3]+Cov[1]],
                    [sum(Var[2:4])+2*Cov[1], Var[2]+Cov[1],  Var[3]+Cov[1], sum(Var[2:4])+2*Cov[1]]
                    ])

    return estimate, etimate_of_mean, I_hat,V_hat

mle, mu=pd.DataFrame(),pd.DataFrame()
I_, V_ = 0, 0
for i in range(2000):
    mle_i, mu_i, I_i, V_i=Estimate(mean_true=true_mean, sample_size=seq_size)
    mle = mle.append(mle_i,ignore_index=True)
    mu = mu.append(mu_i,ignore_index=True)
    I_+=I_i
    V_+=V_i

I_=I_/2000
V_=V_/2000

print('True Value of parameters\n alpha：%f, eta：%f, gamma：%f, delta：%f' %(1.0, 0.3, 1.2, 0.2))
print('MLE\n',mle.mean())
print('Sample variance of parameters\n',2*seq_size*mle.cov(ddof=1))
print('Inverse of matrix I\n',lin.inv(I_))
print('True Value of matrix I\n',I)
print('Estimate of matrix I\n',I_)
print('Estimate of matrix V\n',V_)


'''
def L(params):
    alpha,eta,gamma,delta=params
    f1,f2=0,0
    for i in range(len(data)):
        f1=f1+alpha*data.at[i,'Yi11']-np.exp(alpha)+(alpha+eta+gamma)*data.at[i,'Yi21']-np.exp(alpha+eta+gamma) - np.log(factorial(data.at[i,'Yi11']))-np.log(factorial(data.at[i,'Yi21']))
        f2=f2+(alpha+eta+delta)*data.at[i,'Yi12']-np.exp(alpha+eta+delta)+(alpha+gamma+delta)*data.at[i,'Yi22']-np.exp(alpha+gamma+delta)-np.log(factorial(data.at[i,'Yi12']))-np.log(factorial(data.at[i,'Yi22']))
    
    return f1+f2
initial_guess = np.array([0.013,0.02,1.03,1.032])
solve = scipy.optimize.minimize (lambda params: -L(params), initial_guess)
MLE_opt=np.exp(solve.x)
hessian_ = hessian(L)
I_ = hessian_(np.array([alpha,eta,gamma,delta]))/-240
'''