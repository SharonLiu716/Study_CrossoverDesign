# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 14:54:52 2022

@author: a7086
"""

import scipy
import math
from math import factorial
import numpy as np
import pandas as pd
from scipy.stats import poisson
from scipy.optimize import minimize 
import numpy.linalg as lin
alpha, eta, gamma2, gamma3, delta = 1.0, 0.5, 0.2, 0.2, 0.2
#true mean
tm=[np.exp(alpha), np.exp(alpha + eta + gamma2), np.exp(alpha + eta + gamma3),
             np.exp(alpha + eta + delta), np.exp(alpha + gamma2 + delta),np.exp(alpha + gamma3 + delta)]

pi=60/120
I = pi*np.array([[(tm[0] + tm[1]+ tm[2]) +(tm[3] + tm[4]+ tm[5]),(tm[1]+ tm[2]) +  tm[3], tm[1] +tm[4], tm[2] + tm[5], (tm[3] + tm[4]+ tm[5])],
              [(tm[1]+ tm[2]) + tm[3],  (tm[1]+ tm[2]) + tm[3]  , tm[1], tm[2], tm[3]],
              [tm[1] + tm[4],  tm[1], tm[1]+tm[4] , 0, tm[4]], 
              [tm[2] + tm[5],  tm[2], 0, tm[2]+  tm[5],tm[5] ],
              [(tm[3] + tm[4]+ tm[5]),tm[3], tm[4],tm[5],(tm[3] + tm[4]+ tm[5])]
              ])

def var(x,mean_x): # * removed
    return sum((x - mean_x)**2)/len(x)
def cov(x, y,mean_x,mean_y): 
    return sum((x - mean_x) * (y - mean_y)) / len(x) 



def Estimate(mean_true,sample_size,data_type):
    
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
    Mu = pd.DataFrame(columns = ['Mu11', 'Mu21','Mu31', 'Mu12', 'Mu22','Mu32'])
    data=pd.DataFrame(columns = ['Yi11', 'Yi21','Yi31', 'Yi12', 'Yi22','Yi32'])
    for (idx,param) in enumerate(mean_true):
        if data_type == 'ind':
            #for ind data        
            Mu[Mu.columns[idx]]=pd.DataFrame(np.array(param).repeat(sample_size))
            data[data.columns[idx]]=[ np.random.poisson(p) for p in Mu[Mu.columns[idx]]]
        else:
            #for correlated data
            Mu[Mu.columns[idx]]=pd.DataFrame([param*np.random.gamma(shape=1, scale=1,size=sample_size) ]).transpose()
            data[data.columns[idx]]=[ np.random.poisson(p) for p in Mu[Mu.columns[idx]]]
    
    '''MLE''' 
    
    # alpha
    alpha_hat = np.mean(data['Yi11'])
    # Eta
    eta_hat = np.sqrt(np.mean(data['Yi31']) * np.mean(data['Yi12']) / (np.mean(data['Yi11']) * np.mean(data['Yi32'])))
    # gamma2
    gamma2_hat = np.sqrt(np.mean(data['Yi22']) * np.mean(data['Yi21']) / (np.mean(data['Yi11']) * np.mean(data['Yi12'])))
    # gamma3
    gamma3_hat = np.sqrt(np.mean(data['Yi32']) * np.mean(data['Yi31']) / (np.mean(data['Yi11']) * np.mean(data['Yi12'])))
    # delta
    delta_hat =  np.sqrt(np.mean(data['Yi32']) * np.mean(data['Yi12']) / (np.mean(data['Yi11']) * np.mean(data['Yi31'])))
    
    # MLE
    estimate= pd.DataFrame({'alpha_hat': np.log(alpha_hat), 
                            'eta_hat': np.log(eta_hat), 
                            'gamma2_hat':np.log(gamma2_hat), 
                            'gamma3_hat':np.log(gamma3_hat),
                            'delta_hat': np.log(delta_hat)},index=[0])
    
    '''Mean'''
    etimate_of_mean = pd.DataFrame({'Mean_11': alpha_hat, 
                                    'Mean_21': alpha_hat*eta_hat*gamma2_hat, 
                                    'Mean_31': alpha_hat*eta_hat*gamma3_hat,
                                    'Mean_12': alpha_hat*eta_hat*delta_hat, 
                                    'Mean_22': alpha_hat*gamma2_hat*delta_hat,
                                    'Mean_32': alpha_hat*gamma3_hat*delta_hat},index=[0])
    
    '''Estimate matrix I''' 
    sm=data.mean()
    pi=60/120
    I_hat = pi*np.array([[sum(sm),sum(sm[1:4]), sm[1] +sm[4], sm[2] + sm[5], sum(sm[3:6])],
                  [sum(sm[1:4]), sum(sm[1:4]), sm[1], sm[2], sm[3]],
                  [sm[1] +sm[4],  sm[1], sm[1]+sm[4] , 0, sm[4]], 
                  [sm[2] + sm[5],  sm[2], 0, sm[2]+tm[5],sm[5] ],
                  [sum(sm[3:6]),sm[3], sm[4],sm[5],sum(sm[3:6])]
                  ])
    '''Diagonal variance matrix'''
    
    Var=[var(data[res],sm[idx]) for (idx,res) in enumerate(data.columns)]
    Cov=[cov(data['Yi11'],data['Yi21'],sm[0],sm[1]),cov(data['Yi11'],data['Yi31'],sm[0],sm[2]),cov(data['Yi21'],data['Yi31'],sm[1],sm[2]),
         cov(data['Yi12'],data['Yi22'],sm[3],sm[4]),cov(data['Yi12'],data['Yi32'],sm[3],sm[5]),cov(data['Yi22'],data['Yi32'],sm[4],sm[5])]

    V_hat=pi*np.array([[sum(Var)+2*sum(Cov),sum(Var[1:4])+sum(Cov)-Cov[2]-Cov[5], Var[1]+Var[4]+sum(Cov)-Cov[2]-Cov[5], Var[2]+Var[5]+sum(Cov)-Cov[0]-Cov[3], sum(Var[3:6])+sum(2*Cov[3:6])],
                    [sum(Var[1:4])+sum(Cov)-Cov[2]-Cov[5], sum(Var[1:4])+2*Cov[2], Var[1]+Cov[2]+Cov[3], Var[2]+Cov[2]+Cov[4], Var[3]+Cov[3]+Cov[4]],
                    [Var[1]+Var[4]+sum(Cov)-Cov[2]-Cov[5], Var[1]+Cov[2]+Cov[3], Var[1]+Var[4], Cov[2]+Cov[5],  (Var[4]+Cov[3]+Cov[5])],
                    [Var[2]+Var[5]+sum(Cov)-Cov[0]-Cov[3], Var[2]+Cov[2]+Cov[4], Cov[2]+Cov[5], Var[2]+Var[5], (Var[5]+Cov[4]+Cov[5])],
                    [sum(Var[3:6])+sum(2*Cov[3:6]), Var[3]+Cov[3]+Cov[4], (Var[4]+Cov[3]+Cov[5]),Var[5]+Cov[4]+Cov[5],sum(Var[3:6])+sum(2*Cov[3:6])]
                    ])
    return estimate, etimate_of_mean, I_hat,V_hat

np.random.seed(980716)
mle_ind, mu_ind,mle_cor, mu_cor=pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
I_ind, V_ind,I_cor, V_cor = 0, 0, 0, 0
for i in range(1000):
    #independent
    mle_i, mu_i, I_i, V_i=Estimate(mean_true=tm, sample_size=30,data_type='ind')
    mle_ind = mle_ind.append(mle_i,ignore_index=True)
    mu_ind = mu_ind.append(mu_i,ignore_index=True)
    I_ind+=I_i
    V_ind+=V_i
    #correlated
    mle_i, mu_i, I_i, V_i=Estimate(mean_true=tm, sample_size=30,data_type='cor')
    mle_cor = mle_cor.append(mle_i,ignore_index=True)
    mu_cor = mu_cor.append(mu_i,ignore_index=True)
    I_cor+=I_i
    V_cor+=V_i

I_ind=I_ind/1000
V_ind=V_ind/1000
I_cor=I_cor/1000
V_cor=V_cor/1000

print('Independent Data')
print('True Value of parameters\n alpha：%f, eta：%f, gamma2：%f, gamma3：%f, delta：%f' %(1.0, 0.3, 1.21, 1.27, 0.2))
print('MLE\n',mle_ind.mean())
print('Sample variance of parameters\n',60*mle_ind.var(ddof=1))
print('Inverse of matrix I\n',lin.inv(I_ind))
print('True Value of matrix I\n',I)
print('Sample variance matrix\n',60*mle_ind.cov(ddof=1))
print('inv(I)*V*inv(I)\n',lin.inv(I_ind).dot(V_ind).dot(lin.inv(I_ind)))
print('Estimate of matrix I\n',I_ind)
print('Estimate of matrix V\n',V_ind)

print('Correlated Data')
print('True Value of parameters\n alpha：%f, eta：%f, gamma2：%f, gamma3：%f, delta：%f' %(1.0, 0.3, 1.21, 1.27, 0.2))
print('MLE\n',mle_cor.mean())
print('Sample variance of parameters\n',60*mle_cor.var(ddof=1))
print('Inverse of matrix I\n',lin.inv(I_cor))
print('Estimate of matrix I\n',I_cor)
print('Estimate of matrix V\n',V_cor)
print('inv(I)*V*inv(I)\n',lin.inv(I_cor).dot(V_cor).dot(lin.inv(I_cor)))
print('Sample variance matrix\n',60*mle_cor.cov(ddof=1))