# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:21:32 2022

@author: 懿萱
"""

import scipy
import math
from math import factorial
import numpy as np
import pandas as pd
from scipy.optimize import minimize 
import numpy.linalg as lin
alpha, eta, gamma2, gamma3, delta = 1.0, 0.3, 1.21, 1.27, 0.2
#true mean
tm=[np.exp(alpha), np.exp(alpha + eta + gamma2), np.exp(alpha + eta + gamma3),
             np.exp(alpha + eta + delta), np.exp(alpha + gamma2 + delta),np.exp(alpha + gamma3 + delta)]
seq_size=200
p1,p2=0.5,0.5
I = np.array([[p1 * (tm[0] + tm[1]+ tm[2]) + p2 * (tm[3] + tm[4]+ tm[5]), p1 * (tm[1]+ tm[2]) + p2 * tm[3], p1 * tm[1] + p2 * tm[4], p1 * tm[2] + p2 * tm[5], p2 * (tm[3] + tm[4]+ tm[5])],
              [p1 * (tm[1]+ tm[2]) + p2 * tm[3],  p1 * (tm[1]+ tm[2]) + p2 * tm[3]  , p1 * tm[1], p1 * tm[2], p2* tm[3]],
              [p1 * tm[1] + p2 * tm[4],  p1 * tm[1], p1 * tm[1]+p2 * tm[4] , 0, p2*tm[4]], 
              [p1 * tm[2] + p2 * tm[5], p1 * tm[2], 0, p1 * tm[2]+ p2 * tm[5],p2*tm[5] ],
              [p2 * (tm[3] + tm[4]+ tm[5]),p2* tm[3], p2*tm[4],p2*tm[5],p2*(tm[3] + tm[4]+ tm[5])]
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
            'Yi31': np.random.poisson(lam=mean_true[2], size=sample_size),
            'Yi12': np.random.poisson(lam=mean_true[3], size=sample_size),
            'Yi22': np.random.poisson(lam=mean_true[4], size=sample_size),
            'Yi32': np.random.poisson(lam=mean_true[5], size=sample_size)}
    data = pd.DataFrame(data)
    
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
    p1,p2=0.5,0.5
    I_hat = np.array([[p1 * sum(data.mean()), p1 * (np.mean(data['Yi21'])+np.mean(data['Yi31']))+ p2 *np.mean(data['Yi12']), p1 * (np.mean(data['Yi21']))+ p2 *np.mean(data['Yi22']), p1 * (np.mean(data['Yi31']))+ p2 *np.mean(data['Yi32']), p2 *(np.mean(data['Yi12'])+np.mean(data['Yi22'])+np.mean(data['Yi32']))],
                      [p1 * (np.mean(data['Yi21'])+np.mean(data['Yi31']))+ p2 *np.mean(data['Yi12']) ,p1 * (np.mean(data['Yi21'])+np.mean(data['Yi31']))+ p2 *np.mean(data['Yi12']),p1 * (np.mean(data['Yi21'])),p1 * (np.mean(data['Yi31'])),p2 * (np.mean(data['Yi12']))],
                      [p1 * (np.mean(data['Yi21']))+ p2 *np.mean(data['Yi22']),p1 * (np.mean(data['Yi21'])),p1 * (np.mean(data['Yi21']))+ p2 *np.mean(data['Yi22']),0,p2 * (np.mean(data['Yi22']))],
                      [p1 * (np.mean(data['Yi31']))+ p2 *np.mean(data['Yi32']),p1 * (np.mean(data['Yi31'])),0,p1 * (np.mean(data['Yi31']))+p2 * (np.mean(data['Yi32'])),p2 * (np.mean(data['Yi32']))],
                      [p2 *(np.mean(data['Yi12'])+np.mean(data['Yi22'])+np.mean(data['Yi32'])),p2 * (np.mean(data['Yi12'])),p2 * (np.mean(data['Yi22'])),p2 * (np.mean(data['Yi32'])),p2 *(np.mean(data['Yi12'])+np.mean(data['Yi22'])+np.mean(data['Yi32']))]
                      ])
    '''Diagonal variance matrix
    mean=list(etimate_of_mean.loc[0,:])
    Var=[var(data[res],mean[idx]) for (idx,res) in enumerate(data.columns)]
    Cov=[cov(data['Yi11'],data['Yi21'],mean[0],mean[1]),cov(data['Yi12'],data['Yi22'],mean[3],mean[4]),
         cov(data['Yi21'],data['Yi31'],mean[1],mean[2]),cov(data['Yi22'],data['Yi32'],mean[4],mean[5]),
         cov(data['Yi11'],data['Yi31'],mean[0],mean[2]),cov(data['Yi12'],data['Yi32'],mean[3],mean[5])]

    V_hat=np.array([[p1 * sum(Var), p1 *( Var[1]+Var[2]+Cov[0]+2*Cov[2]+Cov[4] ) + p2*( Var[3]+Cov[1]+Cov[5] ), p1*( Var[1]+Cov[0]+Cov[2] ) + p2*(Var[4]+ Cov[1] +Cov[3]) , p1*(Var[2]++Cov[2]+Cov[4])+p2*(Var[5]+Cov[3]+Cov[5]) , p2*(Var[3]+Var[4]+Var[5]+2*( Cov[1]+Cov[3]+Cov[5]))],
                    [p1 *( Var[1]+Var[2]+Cov[0]+2*Cov[2]+Cov[4] ) + p2*( Var[3]+Cov[1]+Cov[5] ), p1*( Var[1]+Var[2]+2*Cov[2]) + p2* Var[3],  p1*( Var[1]+Cov[2] ) + p2*Cov[1], p1*(Var[2]+Cov[2])+p2*Cov[5], p2*(Var[3]+Cov[1]+Cov[5])],
                    [p1*( Var[1]+Cov[0]+Cov[2] ) + p2*(Var[4]+ Cov[1] +Cov[3]), p1*( Var[1]+Cov[2] ) + p2*Cov[1],  p1*Var[1]+p2*Var[4], p1*Cov[2]+p2*Cov[3],  p2*(Var[4]+Cov[1]+Cov[3])],
                    [p1*(Var[2]++Cov[2]+Cov[4])+p2*(Var[5]+Cov[3]+Cov[5]), p1*(Var[2]+Cov[2])+p2*Cov[5], p1*Cov[2]+p2*Cov[3], p2*(Var[2]+Var[5]),p2*(Var[4]+Cov[3]+Cov[5])],
                    [p2*(Var[3]+Var[4]+Var[5]+2*( Cov[1]+Cov[3]+Cov[5])), p2*(Var[3]+Cov[1]+Cov[5]), p2*(Var[4]+Cov[1]+Cov[3]), p2*(Var[4]+Cov[3]+Cov[5]),p2*(Var[3]+Var[4]+Var[5])]
                    ])'''
    
    '''Diagonal variance matrix'''
    Var=data.var(ddof=1)
    Cov=[ data.Yi11.cov(data.Yi21), data.Yi11.cov(data.Yi31),data.Yi21.cov(data.Yi31), data.Yi12.cov(data.Yi22), data.Yi12.cov(data.Yi32),data.Yi22.cov(data.Yi32)]
    
    V_hat=p1*np.array([[sum(Var)+2*sum(Cov),sum(Var[1:4])+sum(Cov)+Cov[2]-Cov[5], Var[1]+Var[4]+sum(Cov)-Cov[1]-Cov[4], Var[2]+Var[5]+sum(Cov)-Cov[0]-Cov[3], sum(Var[3:6])+2*sum(Cov[3:6])],
                    [sum(Var[1:4])+sum(Cov)+Cov[2]-Cov[5], sum(Var[1:4])+2*Cov[2], Var[1]+Cov[2]+Cov[3], Var[2]+Cov[2]+Cov[4], Var[3]+Cov[3]+Cov[4]],
                    [Var[1]+Var[4]+sum(Cov)-Cov[1]-Cov[4], Var[1]+Cov[2]+Cov[3], Var[1]+Var[4], Cov[2]+Cov[5],  (Var[4]+Cov[3]+Cov[5])],
                    [Var[2]+Var[5]+sum(Cov)-Cov[0]-Cov[3], Var[2]+Cov[2]+Cov[4], Cov[2]+Cov[5], Var[2]+Var[5], (Var[5]+Cov[4]+Cov[5])],
                    [sum(Var[3:6])+2*sum(Cov[3:6]), Var[3]+Cov[3]+Cov[4], (Var[4]+Cov[3]+Cov[5]),Var[5]+Cov[4]+Cov[5],sum(Var[3:6])+2*sum(Cov[3:6])]
                    ])
    return estimate, etimate_of_mean, I_hat,V_hat

mle, mu=pd.DataFrame(),pd.DataFrame()
I_, V_ = 0, 0
for i in range(1000):
    mle_i, mu_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size)
    mle = mle.append(mle_i,ignore_index=True)
    mu = mu.append(mu_i,ignore_index=True)
    I_+=I_i
    V_+=V_i

I_=I_/1000
V_=V_/1000
np.set_printoptions(suppress=True,precision=5)
print('True Value of parameters\n alpha：%f, eta：%f, gamma2：%f, gamma3：%f, delta：%f' %(1.0, 0.3, 1.21, 1.27, 0.2))
print('MLE\n',mle.mean())
print('Sample variance of parameters\n',2*seq_size*mle.cov())
print('Inverse of matrix I\n',lin.inv(I))
print('True Value of matrix I\n',I)
print('Estimate of matrix I\n',I_)
print('Estimate of matrix V\n',V_)
