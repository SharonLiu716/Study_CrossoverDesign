# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 20:53:24 2022

@author: a7086
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
import autograd.numpy as np
from autograd import grad, jacobian, hessian

alpha, eta, gamma1, gamma2, delta1, delta2 = 1.0, 0.3, 1.21, 1.27, 0.22, 0.28
#true mean
tm=[np.exp(alpha), np.exp(alpha + gamma1), np.exp(alpha + eta + gamma2),
    np.exp(alpha +  delta1), np.exp(alpha + eta + gamma1 + delta1),np.exp(alpha + gamma2 + delta1),
    np.exp(alpha + eta +delta2), np.exp(alpha + gamma1 + delta2), np.exp(alpha + gamma2 + delta2)]
seq_size=50
pi=seq_size/(seq_size*3)
I = pi*np.array([  [sum(tm), tm[2]+tm[4]+tm[6], tm[1]+tm[4]+tm[7], tm[2]+tm[5]+tm[8], sum(tm[3:6]), sum(tm[6:9])],
                    [tm[2]+tm[4]+tm[6], tm[2]+tm[4]+tm[6], tm[4], tm[2], tm[4], tm[6]],
                    [tm[1]+tm[4]+tm[7],tm[4],tm[1]+tm[4]+tm[7],0, tm[4],tm[7]],
                    [tm[2]+tm[5]+tm[8],tm[2],0,tm[2]+tm[5]+tm[8], tm[5], tm[8]],
                    [sum(tm[3:6]),tm[4],tm[4],tm[5],sum(tm[3:6]),0],
                    [sum(tm[6:9]),tm[6],tm[7],tm[8],0,sum(tm[6:9])]
                ])



'''
def L(params):
    alpha, eta, gamma1, gamma2, delta1, delta2=params
    f1,f2,f3=0,0,0
    for i in range(len(data)):
        f1=f1+alpha*data.at[i,'Yi11']+(alpha+gamma1)*data.at[i,'Yi21']+(alpha+eta+gamma2)*data.at[i,'Yi31']-np.exp(alpha)-np.exp(alpha+gamma1)-np.exp(alpha+eta+gamma2)
        f2=f2+(alpha+delta1)*data.at[i,'Yi12']+(alpha + eta + gamma1 + delta1)*data.at[i,'Yi22']+(alpha + gamma2 + delta1)*data.at[i,'Yi32']-np.exp(alpha+delta1)-np.exp(alpha+eta+gamma1+delta1)-np.exp(alpha+gamma2+delta1)
        f3=f3+(alpha + eta +delta2)*data.at[i,'Yi13']+(alpha + gamma1 + delta2)*data.at[i,'Yi23']+(alpha + gamma2 + delta2)*data.at[i,'Yi33']-np.exp(alpha+eta+delta2)-np.exp(alpha+gamma1+delta2)-np.exp(alpha+gamma2+delta2)

    return f1+f2+f3
initial_guess = np.array([0.841, 0.218, 1.11, 1.223, 0.202, 0.262])
solve = scipy.optimize.minimize (lambda params: -L(params), initial_guess)
MLE_opt=np.exp(solve.x)
hessian_ = hessian(L)
I_ = hessian_(np.array([alpha,eta,gamma,delta]))/-240
'''

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
    
    def L(params):
        alpha, eta, gamma1, gamma2, delta1, delta2=params
        f1,f2,f3=0,0,0
        for i in range(len(data)):
            f1=f1+alpha*data.at[i,'Yi11']+(alpha+gamma1)*data.at[i,'Yi21']+(alpha+eta+gamma2)*data.at[i,'Yi31']-np.exp(alpha)-np.exp(alpha+gamma1)-np.exp(alpha+eta+gamma2)
            f2=f2+(alpha+delta1)*data.at[i,'Yi12']+(alpha + eta + gamma1 + delta1)*data.at[i,'Yi22']+(alpha + gamma2 + delta1)*data.at[i,'Yi32']-np.exp(alpha+delta1)-np.exp(alpha+eta+gamma1+delta1)-np.exp(alpha+gamma2+delta1)
            f3=f3+(alpha + eta +delta2)*data.at[i,'Yi13']+(alpha + gamma1 + delta2)*data.at[i,'Yi23']+(alpha + gamma2 + delta2)*data.at[i,'Yi33']-np.exp(alpha+eta+delta2)-np.exp(alpha+gamma1+delta2)-np.exp(alpha+gamma2+delta2)

        return f1+f2+f3
    initial_guess = np.array([0.841, 0.218, 1.11, 1.223, 0.202, 0.262])
    solve = scipy.optimize.minimize (lambda params: -L(params), initial_guess)
    MLE_opt=solve.x
    # MLE
    estimate= pd.DataFrame({'alpha_hat': MLE_opt[0], 
                            'eta_hat': MLE_opt[1], 
                            'gamma1_hat':MLE_opt[2], 
                            'gamma2_hat':MLE_opt[3],
                            'delta1_hat': MLE_opt[4],
                            'delta2_hat': MLE_opt[5]},index=[0])  
     
    '''Mean'''
    etimate_of_mean = pd.DataFrame({'Mean_11': np.exp(MLE_opt[0]), 
                                    'Mean_21': np.exp(MLE_opt[0] + MLE_opt[2]), 
                                    'Mean_31': np.exp(MLE_opt[0] + MLE_opt[1] + MLE_opt[3]),
                                    'Mean_12': np.exp(MLE_opt[0] +  MLE_opt[4]), 
                                    'Mean_22': np.exp(MLE_opt[0] + MLE_opt[1] + MLE_opt[2] + MLE_opt[4]),
                                    'Mean_32': np.exp(MLE_opt[0] + MLE_opt[3] + MLE_opt[4]),
                                    'Mean_13': np.exp(MLE_opt[0] + MLE_opt[1] + MLE_opt[5]), 
                                    'Mean_23': np.exp(MLE_opt[0] + MLE_opt[2] + MLE_opt[5]),
                                    'Mean_33': np.exp(MLE_opt[0] + MLE_opt[3] + MLE_opt[5])},index=[0])
    
        
    '''Estimate matrix I'''
    sm=data.mean()##有問題!!MLE!=真實的mean
    p1=seq_size/(seq_size*3)
    I_hat = p1*np.array([  [sum(sm), sm[2]+sm[4]+sm[6], sm[1]+sm[4]+sm[7], sm[2]+sm[5]+sm[8], sum(sm[3:6]), sum(sm[6:9])],
                        [sm[2]+sm[4]+sm[6], sm[2]+sm[4]+sm[6], sm[4], sm[2], sm[4], sm[6]],
                        [sm[1]+sm[4]+sm[7],sm[4],sm[1]+sm[4]+sm[7],0, sm[4],sm[7]],
                        [sm[2]+sm[5]+sm[8],sm[2],0,sm[2]+sm[5]+sm[8], sm[5], sm[8]],
                        [sum(sm[3:6]),sm[4],sm[4],sm[5],sum(sm[3:6]),0],
                        [sum(sm[6:9]),sm[6],sm[7],sm[8],0,sum(sm[6:9])]
                    ])
    
    
    '''Diagonal variance matrix'''
    Var=[var(data[res],sm[idx]) for (idx,res) in enumerate(data.columns)]
    Cov=[cov(data['Yi11'],data['Yi21'],sm[0],sm[1]), cov(data['Yi11'],data['Yi31'],sm[0],sm[2]), cov(data['Yi21'],data['Yi31'],sm[1],sm[2]),
        cov(data['Yi12'],data['Yi22'],sm[3],sm[4]), cov(data['Yi12'],data['Yi32'],sm[3],sm[5]), cov(data['Yi22'],data['Yi32'],sm[4],sm[5]),
        cov(data['Yi13'],data['Yi23'],sm[6],sm[7]), cov(data['Yi13'],data['Yi33'],sm[6],sm[8]), cov(data['Yi23'],data['Yi33'],sm[7],sm[8])]

    V_hat=pi*np.array([
        [sum(Var), Var[2]+Var[4]+Var[6]+sum(Cov[1:4])+sum(Cov[5:8]),Var[1]+Var[4]+Var[7]+Cov[0]+Cov[2]+Cov[3]+Cov[5]+Cov[6]+Cov[8],Var[2]+Var[5]+Var[8]+Cov[1]+Cov[2]+Cov[4]+Cov[5]+Cov[7]+Cov[8],sum(Var[3:6])+2*sum(Cov[3:6]),sum(Var[6:9])+2*sum(Cov[6:9])],
        [Var[2]+Var[4]+Var[6]+sum(Cov[1:4])+sum(Cov[5:8]),Var[2]+Var[4]+Var[6],Var[4]+Cov[2]+Cov[6],Var[2]+Cov[5]+Cov[7],Var[4]+Cov[3]+Cov[5],Var[6]+Cov[6]+Cov[7]],
        [Var[1]+Var[4]+Var[7]+Cov[0]+Cov[2]+Cov[3]+Cov[5]+Cov[6]+Cov[8],Var[4]+Cov[2]+Cov[6],Var[1]+Var[4]+Var[7],Cov[2]+Cov[5]+Cov[8],Var[4]+Cov[3]+Cov[5],Var[5]+Cov[6]+Cov[8]],
        [Var[2]+Var[5]+Var[8]+Cov[1]+Cov[2]+Cov[4]+Cov[5]+Cov[7]+Cov[8],Var[2]+Cov[5]+Cov[7],Cov[2]+Cov[5]+Cov[8],Var[2]+Var[5]+Var[8],Var[5]+sum(Cov[4:6]),Var[8]+sum(Cov[7:9])],
        [sum(Var[3:6])+2*sum(Cov[3:6]),Var[4]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[5]+sum(Cov[4:6]),sum(Var[3:6])+2*sum(Cov[3:6]),0],
        [sum(Var[6:9])+2*sum(Cov[6:9]),Var[6]+Cov[6]+Cov[7],Var[5]+Cov[6]+Cov[8],Var[8]+sum(Cov[7:9]),0,sum(Var[6:9])+2*sum(Cov[6:9])]
                                    ])

    return estimate, etimate_of_mean, I_hat,V_hat


np.random.seed(980716)
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

print('True Value of parameters\n alpha：%f, eta：%f, gamma1：%f, gamma2：%f, delta1：%f, delta1：%f' %(1.0, 0.3, 1.21, 1.27, 0.22, 0.28))
print('MLE\n',mle.mean())
print('Sample variance of parameters\n',3*seq_size*mle.cov())
print('Inverse of matrix I\n',lin.inv(I))
print('Inverse of matrix I_hat\n',lin.inv(I_))
print('inv(I)*V*inv(I)\n',lin.inv(I_).dot(V_).dot(lin.inv(I_)))
print('True Value of matrix I\n',I)
print('Estimate of matrix I\n',I_)
print('Estimate of matrix V\n',V_)

TEMP=V_-lin.inv(I_).dot(V_).dot(lin.inv(I_))
lin.inv(TEMP)
