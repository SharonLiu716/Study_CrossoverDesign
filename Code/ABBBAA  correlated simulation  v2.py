# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 11:20:12 2022

@author: a7086
"""

import scipy
import math
from math import factorial
import numpy as np
import pandas as pd
from scipy.optimize import minimize 
from scipy.stats import poisson
import numpy.linalg as lin

sim_time=3000
alpha, eta, gamma1, gamma2, delta = 1.0, 0.7, 0.12, 0.33, 0.21
params = np.array([ [alpha], [eta], [gamma1], [gamma2], [delta]])
covariate=[[1,0,0,0,0],[1,1,1,0,0],[1,1,0,1,0],[1,1,0,0,1],[1,0,1,0,1],[1,0,0,1,1]]
#true mean
tm=np.exp(np.array([ xijk  for xijk in covariate]).dot(params)).reshape(1,6).tolist()[0]
gamma_param=0.1
seq_size=100
pi=(seq_size)/(seq_size*2)
I = pi*np.array([[sum(tm),sum(tm[1:4]), tm[1] +tm[4], tm[2] + tm[5], sum(tm[3:6])],
              [sum(tm[1:4]), sum(tm[1:4]),  tm[1], tm[2], tm[3]],
              [tm[1] + tm[4],  tm[1], tm[1] + tm[4] , 0, tm[4]], 
              [tm[2] + tm[5],  tm[2], 0, tm[2]+ tm[5],tm[5] ],
              [sum(tm[3:6]),tm[3], tm[4],tm[5],sum(tm[3:6])]
              ])

def L(data,params):
    alpha, eta, gamma1, gamma2, delta,=params
    f1,f2=0,0
    for i in range(len(data)):
        f1=f1+alpha*data.at[i,'Yi11']+(alpha+eta+gamma1)*data.at[i,'Yi21']+(alpha+eta+gamma2)*data.at[i,'Yi31']
        f2=f2+(alpha+eta+delta)*data.at[i,'Yi12']+(alpha + gamma1 + delta)*data.at[i,'Yi22']+(alpha + gamma2 + delta)*data.at[i,'Yi32']    
    return f1+f2
def matrix_AB(I_,V_):
    loc=np.arange(0,len(I_)).tolist()
    loc.pop(1)
    I_eta=I_[1, 1]
    I_etapsy=np.array([I_[1,i] for i in loc])#np.array([I_[1,0],I_[1, 2],I_[1, 3],I_[1, 4]])
    I_psy=np.array([[I_[i,0],I_[i, 2],I_[i, 3],I_[i, 4]] for i in loc])
    #np.array([[I_[0,0],I_[0, 2],I_[0, 3],I_[0, 4]],[I_[2,0],I_[2, 2],I_[2, 3],I_[2, 4]],[I_[3,0],I_[3, 2],I_[3, 3],I_[3, 4]],[I_[4,0],I_[4, 2],I_[4, 3],I_[4, 4]]])
    V_eta=V_[1, 1]
    V_etapsy=np.array([V_[1,i] for i in loc])#np.array([V_[1,0],V_[1, 2],V_[1, 3],V_[1, 4]])
    V_psy=np.array([[V_[i,0],V_[i, 2],V_[i, 3],V_[i, 4]] for i in loc])
    #np.array([[V_[0,0],V_[0, 2],V_[0, 3],V_[0, 4]],[V_[2,0],V_[2, 2],V_[2, 3],V_[2, 4]],[V_[3,0],V_[3, 2],V_[3, 3],V_[3, 4]],[V_[4,0],V_[4, 2],V_[4, 3],V_[4, 4]]])
    A=I_eta-I_etapsy.dot(lin.inv(I_psy)).dot(I_etapsy.reshape(len(I_)-1,1))
    B=V_eta-2*I_etapsy.dot(lin.inv(I_psy)).dot(V_etapsy.reshape(len(I_)-1,1))+I_etapsy.dot(lin.inv(I_psy)).dot(V_psy).dot(lin.inv(I_psy)).dot(I_etapsy.reshape(len(I_)-1,1))
    return A,B


def Estimate(mean_true,sample_size,data_type,cor_par,eta0):
    #mean_true=tm
    #sample_size=100000
    #cor_par=gamma_param
    #sample_size=int(sample_size/3)
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
        data=pd.DataFrame(np.array([np.random.poisson(lam=p, size=sample_size) for p in mean_true]).T.tolist(),columns = ['Yi11', 'Yi21','Yi31', 'Yi12', 'Yi22','Yi32'])
    else:     
        #np.random.seed(980716)
        Mu_cor = pd.DataFrame(columns = ['Mu11', 'Mu21','Mu31', 'Mu12', 'Mu22','Mu32'])
        for i in range(0,len(mean_true),3):            
            nu=np.random.gamma(1/cor_par,cor_par,sample_size)
            #print(i,nu)
            for (ix,param) in enumerate(mean_true[i:i + 3]):
                Mu_cor[Mu_cor.columns[ix+i]]=pd.DataFrame(np.multiply(np.array(param).repeat(sample_size), nu).T)
                #print(i,ix,param)

        data=pd.DataFrame(columns = ['Yi11', 'Yi21','Yi31', 'Yi12', 'Yi22','Yi32'])
        for (idx,mu) in enumerate(Mu_cor.columns):
            data[data.columns[idx]]=pd.DataFrame([np.random.poisson(p) for p in Mu_cor[mu]])
        
    #data.corr()
    #data.cov()
   
    '''MLE''' 
    
    # alpha
    alpha_hat = np.log(np.mean(data['Yi11']))
    # Eta
    eta_hat =0.25*(np.log(np.mean(data['Yi21']))+np.log(np.mean(data['Yi31']))+2*np.log(np.mean(data['Yi12']))-np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi32']))-2*np.log(np.mean(data['Yi11'])))
    # gamma2
    gamma2_hat = 0.5*(np.log(np.mean(data['Yi21']))+np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi12']))-np.log(np.mean(data['Yi11'])))
    # gamma3
    gamma3_hat = 0.5*(np.log(np.mean(data['Yi31']))+np.log(np.mean(data['Yi32']))-np.log(np.mean(data['Yi12']))-np.log(np.mean(data['Yi11'])))
    # delta
    delta_hat = 0.5*(np.log(np.mean(data['Yi12']))-np.log(np.mean(data['Yi11'])))-0.25*(np.log(np.mean(data['Yi21']))+np.log(np.mean(data['Yi31']))-np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi32'])))
    
    # MLE
    estimate= pd.DataFrame({'alpha_hat': alpha_hat, 
                            'eta_hat': eta_hat, 
                            'gamma2_hat':gamma2_hat, 
                            'gamma3_hat':gamma3_hat,
                            'delta_hat': delta_hat},index=[0])
    
    '''Mean'''
    etimate_of_mean = pd.DataFrame({'Mean_11': np.exp(alpha_hat), 
                                    'Mean_21': np.exp(alpha_hat+eta_hat+gamma2_hat), 
                                    'Mean_31': np.exp(alpha_hat+eta_hat+gamma3_hat),
                                    'Mean_12': np.exp(alpha_hat+eta_hat+delta_hat), 
                                    'Mean_22': np.exp(alpha_hat+gamma2_hat+delta_hat),
                                    'Mean_32': np.exp(alpha_hat+gamma3_hat+delta_hat)},index=[0])
    
    '''Estimate matrix I''' 
    sm=list(etimate_of_mean.loc[0,:])#data.mean()
    pi=sample_size/(sample_size*2)
    I_hat = pi*np.array([[sum(sm),sum(sm[1:4]), sm[1] +sm[4], sm[2] + sm[5], sum(sm[3:6])],
                  [sum(sm[1:4]), sum(sm[1:4]), sm[1], sm[2], sm[3]],
                  [sm[1] +sm[4],  sm[1], sm[1]+sm[4] , 0, sm[4]], 
                  [sm[2] + sm[5],  sm[2], 0, sm[2]+tm[5],sm[5] ],
                  [sum(sm[3:6]),sm[3], sm[4],sm[5],sum(sm[3:6])]
                  ])
    '''Diagonal variance matrix'''
    Var=data.var(ddof=1)
    Cov=[ data.Yi11.cov(data.Yi21), data.Yi11.cov(data.Yi31),data.Yi21.cov(data.Yi31), data.Yi12.cov(data.Yi22), data.Yi12.cov(data.Yi32),data.Yi22.cov(data.Yi32)]
    
    V_hat=pi*np.array([[sum(Var)+2*sum(Cov),sum(Var[1:4])+sum(Cov)+Cov[2]-Cov[5], Var[1]+Var[4]+sum(Cov)-Cov[1]-Cov[4], Var[2]+Var[5]+sum(Cov)-Cov[0]-Cov[3], sum(Var[3:6])+2*sum(Cov[3:6])],
                    [sum(Var[1:4])+sum(Cov)+Cov[2]-Cov[5], sum(Var[1:4])+2*Cov[2], Var[1]+Cov[2]+Cov[3], Var[2]+Cov[2]+Cov[4], Var[3]+Cov[3]+Cov[4]],
                    [Var[1]+Var[4]+sum(Cov)-Cov[1]-Cov[4], Var[1]+Cov[2]+Cov[3], Var[1]+Var[4], Cov[2]+Cov[5],  (Var[4]+Cov[3]+Cov[5])],
                    [Var[2]+Var[5]+sum(Cov)-Cov[0]-Cov[3], Var[2]+Cov[2]+Cov[4], Cov[2]+Cov[5], Var[2]+Var[5], (Var[5]+Cov[4]+Cov[5])],
                    [sum(Var[3:6])+2*sum(Cov[3:6]), Var[3]+Cov[3]+Cov[4], (Var[4]+Cov[3]+Cov[5]),Var[5]+Cov[4]+Cov[5],sum(Var[3:6])+2*sum(Cov[3:6])]
                    ])
    '''
    A_hat,B_hat=matrix_AB(I_hat,V_hat)
    ##推導完待修目前是ABBA的
    eta_null=eta0
    gamma1_null=
    gamma2_null=
    delta_null=
    LR_naive=2*(L(data=data,params=[alpha_hat,eta_hat,gamma2_hat,gamma3_hat,delta_hat])-L(data=data,params=[alpha_hat,eta_null,gamma2_null,gamma3_null,delta_null])) 
    LR_robust=2*(A_hat/B_hat)*(L(data=data,params=[alpha_hat,eta_hat,gamma2_hat,gamma3_hat,delta_hat])-L(data=data,params=[alpha_hat,eta_null,gamma2_null,gamma3_null,delta_null]))
    Wald_naive=(sample_size*4)*(eta_hat-eta_null)*A_hat*(eta_hat-eta_null)
    Wald_robust=(sample_size*4)*(eta_hat-eta_null)*(A_hat/B_hat)*(eta_hat-eta_null)
    '''
    return estimate, I_hat, V_hat#, LR_naive, LR_robust,Wald_naive,Wald_robust



np.random.seed(980716)
mle_ind, mu_ind,mle_cor, mu_cor=pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
I_ind, V_ind,I_cor, V_cor = 0, 0, 0, 0
for i in range(sim_time):
    
    #independent
    mle_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size,data_type='ind',cor_par=gamma_param,eta0=0)
    mle_ind = mle_ind.append(mle_i,ignore_index=True)
    I_ind+=I_i
    V_ind+=V_i
    #correlated
    mle_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size,data_type='cor',cor_par=gamma_param,eta0=0)
    mle_cor = mle_cor.append(mle_i,ignore_index=True)
    I_cor+=I_i
    V_cor+=V_i
    '''
np.random.seed(980716)
mle_ind,mle_cor=pd.DataFrame(),pd.DataFrame()
I_ind, V_ind,I_cor, V_cor = 0, 0, 0, 0
LR_ind_na,LR_ind_rb,LR_cor_na,LR_cor_rb=np.empty((0,1), float),np.empty((0,1), float),np.empty((0,1), float),np.empty((0,1), float)
Wald_ind_na,Wald_ind_rb,Wald_cor_na,Wald_cor_rb=np.empty((0,1), float),np.empty((0,1), float),np.empty((0,1), float),np.empty((0,1), float)
for i in range(sim_time):    
    #independent
    #mle_i, mu_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size,data_type='ind',cor_par=gamma_param)
    mle_i, I_i, V_i, LR_na, LR_rb, Wald_na, Wald_rb=Estimate(mean_true=tm, sample_size=seq_size,data_type='ind',cor_par=gamma_param,eta0=0)
    mle_ind = mle_ind.append(mle_i,ignore_index=True)
    #mu_ind = mu_ind.append(mu_i,ignore_index=True)
    I_ind+=I_i
    V_ind+=V_i
    LR_ind_na=np.append(LR_ind_na, LR_na)
    LR_ind_rb=np.append(LR_ind_rb, LR_rb)
    Wald_ind_na=np.append(Wald_ind_na, Wald_na)
    Wald_ind_rb=np.append(Wald_ind_rb, Wald_rb)
    #correlated
    #mle_i, mu_i, I_i, V_i=Estimate(mean_true=tm, sample_size=seq_size,data_type='cor',cor_par=gamma_param,eta0=0)
    mle_i, I_i, V_i, LR_na, LR_rb, Wald_na, Wald_rb=Estimate(mean_true=tm, sample_size=seq_size,data_type='cor',cor_par=gamma_param,eta0=0)
    mle_cor = mle_cor.append(mle_i,ignore_index=True)
    I_cor+=I_i
    V_cor+=V_i
    LR_cor_na=np.append(LR_cor_na, LR_na)
    LR_cor_rb=np.append(LR_cor_rb, LR_rb)
    Wald_cor_na=np.append(Wald_cor_na, Wald_na)
    Wald_cor_rb=np.append(Wald_cor_rb, Wald_rb)'''
    
I_ind=I_ind/sim_time
V_ind=V_ind/sim_time
I_cor=I_cor/sim_time
V_cor=V_cor/sim_time

np.set_printoptions(suppress=True,precision=5)
print('Independent Data, seq_size = ',seq_size)
print('True Value of parameters\n alpha：%f, eta：%f, gamma2：%f, gamma3：%f, delta：%f' %(1.0, 0.7, 0.12, 0.33, 0.21))
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
print('True Value of parameters\n alpha：%f, eta：%f, gamma2：%f, gamma3：%f, delta：%f' %(1.0, 0.5, 0.2, 0.2, 0.2))
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