# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 16:12:52 2022

@author: user1
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:56:49 2022

@author: a7086
"""

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
alpha, eta, gamma1, gamma2, delta1, delta2 = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
#true mean
tm=[np.exp(alpha), np.exp(alpha + gamma1), np.exp(alpha + eta + gamma2),
    np.exp(alpha +  delta1), np.exp(alpha + eta + gamma1 + delta1),np.exp(alpha + gamma2 + delta1),
    np.exp(alpha + eta +delta2), np.exp(alpha + gamma1 + delta2), np.exp(alpha + gamma2 + delta2)]

p1,p2,p3=50/150,50/150,50/150
I = np.array([  [p1 * (tm[0] + tm[1]+ tm[2]) + p2 * (tm[3] + tm[4]+ tm[5])+ p3 * (tm[6] + tm[7]+ tm[8]),
                 p1 * tm[2] + p2 * tm[4] + p3 * tm[6],
                 p1 * tm[1] + p2 * tm[4] + p3 * tm[7],
                 p1 * tm[2] + p2 * tm[5] + p3 * tm[8],
                 p2 * (tm[3] + tm[4]+ tm[5]),
                 p3 * (tm[6] + tm[7]+ tm[8])],
                [p1 * tm[2] + p2 * tm[4] + p3 * tm[6],
                 p1 * tm[2] + p2 * tm[4]+ p3 * tm[6], 
                 p2 * tm[4], 
                 p1 * tm[2], 
                 p2 * tm[4], 
                 p3 * tm[6]],
                [p1 * tm[1] + p2 * tm[4] + p3 * tm[7],
                 p2 * tm[4],
                 p1 * tm[1] + p2 * tm[4] + p3 * tm[7],
                 0, 
                 p2 * tm[4], 
                 p3 * tm[7]],
                [p1 * tm[2] + p2 * tm[5] + p3 * tm[8],
                 p1 * tm[2],
                 0,
                 p1 * tm[2] + p2 * tm[5]+ p3 * tm[8], 
                 p2 * tm[5], 
                 p3 * tm[8]],
                [p2 * (tm[3] + tm[4]+ tm[5]),
                 p2 * tm[4],
                 p2 * tm[4],
                  p2 * tm[5],
                 p2 * (tm[3] + tm[4]+ tm[5]), 
                 0 ],
                [p3 * (tm[6] + tm[7]+ tm[8]),
                 p3 * tm[6],
                 p3 * tm[7],
                 p3 * tm[8],
                 0,
                 p3 * (tm[6] + tm[7]+ tm[8])]
                ])


def var(x,mean_x): # * removed
    return sum((x - mean_x)**2)/len(x)
def cov(x, y,mean_x,mean_y): 
    return sum((x - mean_x) * (y - mean_y)) / len(x) 

def Estimate(mean_true,sample_size):
    #mean_true=tm
    #sample_size=25
    
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
            'Yi32': np.random.poisson(lam=mean_true[5], size=sample_size),
            'Yi13': np.random.poisson(lam=mean_true[6], size=sample_size),
            'Yi23': np.random.poisson(lam=mean_true[7], size=sample_size),
            'Yi33': np.random.poisson(lam=mean_true[8], size=sample_size)}
    data = pd.DataFrame(data)
    
    '''MLE''' 
    mean=np.log(data.mean())
    
    # alpha
    alpha_hat = mean[0]
    # Eta1
    eta_hat = mean[4]+mean[0]-mean[1]-mean[3]
    # gamma1
    gamma1_hat = mean[1]-mean[0]
    # gamma2
    gamma2_hat = mean[5]-mean[3]
    # delta1
    delta1_hat = mean[1]-mean[0]
    # delta
    delta2_hat = mean[7]-mean[1]
    
    # MLE
    estimate= pd.DataFrame({'alpha_hat': alpha_hat, 
                            'eta_hat': eta_hat, 
                            'gamma1_hat':gamma1_hat, 
                            'gamma2_hat':gamma2_hat,
                            'delta1_hat': delta1_hat,
                            'delta2_hat': delta2_hat},index=[0])
     
    '''Mean'''
    etimate_of_mean = pd.DataFrame({'Mean_11': np.exp(alpha_hat), 
                                    'Mean_21': np.exp(alpha_hat + gamma1_hat), 
                                    'Mean_31': np.exp(alpha_hat + eta_hat + gamma2_hat),
                                    'Mean_12': np.exp(alpha_hat +  delta1_hat), 
                                    'Mean_22': np.exp(alpha_hat + eta_hat + gamma1_hat + delta1_hat),
                                    'Mean_32': np.exp(alpha_hat + gamma2_hat + delta1_hat),
                                    'Mean_13': np.exp(alpha_hat + eta_hat +delta2_hat), 
                                    'Mean_23': np.exp(alpha_hat + gamma1_hat + delta2_hat),
                                    'Mean_33': np.exp(alpha_hat + gamma2_hat + delta2_hat)},index=[0])
    sm=list(etimate_of_mean.loc[0,:])
    return estimate, etimate_of_mean


mle, mu=pd.DataFrame(),pd.DataFrame()
I_, V_ = 0, 0
for i in range(1000):
    mle_i, mu_i=Estimate(mean_true=tm, sample_size=25)
    mle = mle.append(mle_i,ignore_index=True)
    mu = mu.append(mu_i,ignore_index=True)

I_=I_/1000
V_=V_/1000
np.set_printoptions(suppress=True,precision=5)
print('True Value of parameters\n alpha：%f, eta：%f, gamma1：%f, gamma2：%f, delta1：%f, delta1：%f' %(1.0, 0.3, 1.21, 1.27, 0.22, 0.28))
print('MLE\n',mle.mean())
print('Sample variance of parameters\n',50*mle.var())


mle, mu=pd.DataFrame(),pd.DataFrame()
I_, V_ = 0, 0
for i in range(1000):
    mle_i, mu_i, I_i, V_i=Estimate(mean_true=tm, sample_size=25)
    mle = mle.append(mle_i,ignore_index=True)
    mu = mu.append(mu_i,ignore_index=True)
    I_+=I_i
    V_+=V_i

I_=I_/1000
V_=V_/1000
np.set_printoptions(suppress=True,precision=5)
print('True Value of parameters\n alpha：%f, eta：%f, gamma1：%f, gamma2：%f, delta1：%f, delta1：%f' %(1.0, 0.3, 1.21, 1.27, 0.22, 0.28))
print('MLE\n',mle.mean())
print('Sample variance of parameters\n',50*mle.var())
print('Inverse of matrix I\n',lin.inv(I))
print('True Value of matrix I\n',I)
print('Estimate of matrix I\n',I_)
print('Estimate of matrix V\n',V_)

TEMP=V_-lin.inv(I_).dot(V_).dot(lin.inv(I_))
lin.inv(TEMP)
