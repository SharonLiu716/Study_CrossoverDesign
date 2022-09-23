# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 16:45:54 2022

@author: a7086
"""

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
from scipy.stats import poisson
import numpy.linalg as lin
import autograd.numpy as np
from autograd import grad, jacobian, hessian
sim_time=1000
alpha, eta, gamma1, gamma2, delta1, delta2 = 1.0, 0.5, 0.2, 0.2, 0.2, 0.2
tm=[np.exp(alpha), np.exp(alpha + gamma1), np.exp(alpha + eta + gamma2),
    np.exp(alpha +  delta1), np.exp(alpha + eta + gamma1 + delta1),np.exp(alpha + gamma2 + delta1),
    np.exp(alpha + eta +delta2), np.exp(alpha + gamma1 + delta2), np.exp(alpha + gamma2 + delta2)]
params = np.array([ [alpha], [eta], [gamma1], [gamma2], [delta1], [delta2]])
covariate=[[1,0,0,0,0,0],[1,0,1,0,0,0],[1,1,0,1,0,0],[1,0,0,0,1,0],[1,1,1,0,1,0],[1,0,0,1,1,0],[1,1,0,0,0,1],[1,0,1,0,0,1],[1,0,0,1,0,1]]
#true mean
tm=np.exp(np.array([ xijk   for xijk in covariate]).dot(params))
gamma_param=0.1
seq_size=75#150,225,300
pi=(seq_size)/(seq_size*9)
I = pi*np.array([  [sum(tm), tm[2]+tm[4]+tm[6], tm[1]+tm[4]+tm[7], tm[2]+tm[5]+tm[8], sum(tm[3:6]), sum(tm[6:9])],
                    [tm[2]+tm[4]+tm[6], tm[2]+tm[4]+tm[6], tm[4], tm[2], tm[4], tm[6]],
                    [tm[1]+tm[4]+tm[7],tm[4],tm[1]+tm[4]+tm[7],0, tm[4],tm[7]],
                    [tm[2]+tm[5]+tm[8],tm[2],0,tm[2]+tm[5]+tm[8], tm[5], tm[8]],
                    [sum(tm[3:6]),tm[4],tm[4],tm[5],sum(tm[3:6]),0],
                    [sum(tm[6:9]),tm[6],tm[7],tm[8],0,sum(tm[6:9])]
                ])




def Estimate(mean_true,sample_size,data_type,cor_par,eta0):
    #mean_true=tm
    #sample_size=100000
    #cor_par=gamma_param
    sample_size=int(sample_size/3)    
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
        data=pd.DataFrame(np.array([np.random.poisson(lam=p, size=sample_size) for p in mean_true]).T.tolist(),columns =['Yi11', 'Yi21','Yi31', 'Yi12', 'Yi22','Yi32','Yi13', 'Yi23','Yi33'])
    else:     
        #np.random.seed(980716)
        Mu_cor = pd.DataFrame(columns = ['Mu11', 'Mu21','Mu31', 'Mu12', 'Mu22','Mu32','Mu13', 'Mu23','Mu33'])
        for i in range(0,len(mean_true),3):            
            nu=np.random.gamma(1/cor_par,cor_par,sample_size)
            #print(i,nu)
            for (ix,param) in enumerate(mean_true[i:i + 3]):
                Mu_cor[Mu_cor.columns[ix+i]]=pd.DataFrame(np.multiply(np.array(param).repeat(sample_size), nu).T)
                #print(i,ix,param)

        data=pd.DataFrame(columns = ['Yi11', 'Yi21','Yi31', 'Yi12', 'Yi22','Yi32','Yi13', 'Yi23','Yi33'])
        for (idx,mu) in enumerate(Mu_cor.columns):
            data[data.columns[idx]]=pd.DataFrame([np.random.poisson(p) for p in Mu_cor[mu]])
        
    #temp=data.corr()
    #data.cov()
    
    '''MLE'''   
    def L(params):
        alpha, eta, gamma1, gamma2, delta1, delta2=params
        f1,f2,f3=0,0,0
        for i in range(len(data)):
            f1=f1+alpha*data.at[i,'Yi11']+(alpha+gamma1)*data.at[i,'Yi21']+(alpha+eta+gamma2)*data.at[i,'Yi31']-np.exp(alpha)-np.exp(alpha+gamma1)-np.exp(alpha+eta+gamma2)
            f2=f2+(alpha+delta1)*data.at[i,'Yi12']+(alpha + eta + gamma1 + delta1)*data.at[i,'Yi22']+(alpha + gamma2 + delta1)*data.at[i,'Yi32']-np.exp(alpha+delta1)-np.exp(alpha+eta+gamma1+delta1)-np.exp(alpha+gamma2+delta1)
            f3=f3+(alpha + eta +delta2)*data.at[i,'Yi13']+(alpha + gamma1 + delta2)*data.at[i,'Yi23']+(alpha + gamma2 + delta2)*data.at[i,'Yi33']-np.exp(alpha+eta+delta2)-np.exp(alpha+gamma1+delta2)-np.exp(alpha+gamma2+delta2)

        return f1+f2+f3
    initial_guess = np.array([0.841, 0.318, 0.11, 0.223, 0.202, 0.262])
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
    p1=(seq_size)/(seq_size*9)
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
        [sum(Var)+2*sum(Cov), Var[2]+Var[4]+Var[6]+sum(Cov[1:4])+sum(Cov[5:8]),Var[1]+Var[4]+Var[7]+Cov[0]+Cov[2]+Cov[3]+Cov[5]+Cov[6]+Cov[8],Var[2]+Var[5]+Var[8]+Cov[1]+Cov[2]+Cov[4]+Cov[5]+Cov[7]+Cov[8],sum(Var[3:6])+2*sum(Cov[3:6]),sum(Var[6:9])+2*sum(Cov[6:9])],
        [Var[2]+Var[4]+Var[6]+sum(Cov[1:4])+sum(Cov[5:8]),Var[2]+Var[4]+Var[6],Var[4]+Cov[2]+Cov[6],Var[2]+Cov[5]+Cov[7],Var[4]+Cov[3]+Cov[5],Var[6]+Cov[6]+Cov[7]],
        [Var[1]+Var[4]+Var[7]+Cov[0]+Cov[2]+Cov[3]+Cov[5]+Cov[6]+Cov[8],Var[4]+Cov[2]+Cov[6],Var[1]+Var[4]+Var[7],Cov[2]+Cov[5]+Cov[8],Var[4]+Cov[3]+Cov[5],Var[5]+Cov[6]+Cov[8]],
        [Var[2]+Var[5]+Var[8]+Cov[1]+Cov[2]+Cov[4]+Cov[5]+Cov[7]+Cov[8],Var[2]+Cov[5]+Cov[7],Cov[2]+Cov[5]+Cov[8],Var[2]+Var[5]+Var[8],Var[5]+sum(Cov[4:6]),Var[8]+sum(Cov[7:9])],
        [sum(Var[3:6])+2*sum(Cov[3:6]),Var[4]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[5]+sum(Cov[4:6]),sum(Var[3:6])+2*sum(Cov[3:6]),0],
        [sum(Var[6:9])+2*sum(Cov[6:9]),Var[6]+Cov[6]+Cov[7],Var[5]+Cov[6]+Cov[8],Var[8]+sum(Cov[7:9]),0,sum(Var[6:9])+2*sum(Cov[6:9])]
                                    ])

    return estimate, etimate_of_mean, I_hat,V_hat


np.random.seed(980716)
mle_ind, mu_ind,mle_cor, mu_cor=pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
I_ind, V_ind,I_cor, V_cor = 0, 0, 0, 0
for i in range(sim_time):
    #independent
    mle_i, mu_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size,data_type='ind',cor_par=gamma_param,eta0=0)
    mle_ind = mle_ind.append(mle_i,ignore_index=True)
    mu_ind = mu_ind.append(mu_i,ignore_index=True)
    I_ind+=I_i
    V_ind+=V_i
    #correlated
    mle_i, mu_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size,data_type='cor',cor_par=gamma_param,eta0=0)
    mle_cor = mle_cor.append(mle_i,ignore_index=True)
    mu_cor = mu_cor.append(mu_i,ignore_index=True)
    I_cor+=I_i
    V_cor+=V_i

I_ind=I_ind/sim_time
V_ind=V_ind/sim_time
I_cor=I_cor/sim_time
V_cor=V_cor/sim_time


np.set_printoptions(suppress=True,precision=5)
print('Independent Data, seq_size = ',seq_size)
print('True Value of parameters\n alpha：%f, eta：%f, gamma1：%f, gamma2：%f, delta1：%f, delta1：%f' %(1.0, 0.5, 0.2, 0.2, 0.2, 0.2))
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
scipy.stats.t.sf(abs(np.mean(LR_ind_na)), df=1)*2

print('Correlated Data (alpha =10 ,beta = 0.1), seq_size = ',seq_size)
print('True Value of parameters\n alpha：%f, eta：%f, gamma1：%f, gamma2：%f, delta1：%f, delta1：%f' %(1.0, 0.5, 0.2, 0.2, 0.2, 0.2))
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