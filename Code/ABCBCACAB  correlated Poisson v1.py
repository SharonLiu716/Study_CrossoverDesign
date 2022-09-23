# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 17:20:39 2022

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
sim_time=1000
alpha, eta1, eta2, gamma1, gamma2, delta1, delta2 = 1.0, 0.55, 0.45, 0.21, 0.26, 0.13, 0.17
params = np.array([ [alpha], [eta1], [eta2], [gamma1], [gamma2], [delta1], [delta2]])
covariate=[[1,0,0,0,0,0,0],[1,1,0,1,0,0,0],[1,0,1,0,1,0,0],[1,1,0,0,0,1,0],[1,0,1,1,0,1,0],[1,0,0,0,1,1,0],[1,0,1,0,0,0,1],[1,0,0,1,0,0,1],[1,1,0,0,1,0,1]]
#true mean
tm=np.exp(np.array([ xijk   for xijk in covariate]).dot(params)).reshape(1,9).tolist()[0]
gamma_param=0.1
seq_size=200
pi=(seq_size)/(seq_size*3)
I = pi*np.array([  [sum(tm), tm[1]+tm[3]+tm[8], tm[2]+tm[4]+tm[6], tm[1]+tm[4]+tm[7], tm[2]+tm[5]+tm[8], tm[3]+tm[4]+tm[5], tm[6]+tm[7]+tm[8]],
                [tm[1]+tm[3]+tm[8],tm[1]+tm[3]+tm[8], 0, tm[1], tm[8], tm[3], tm[8] ],
                [tm[2]+tm[4]+tm[6],0, tm[2]+tm[4]+tm[6], tm[4], tm[2], tm[4], tm[6]],
                [tm[1]+tm[4]+tm[7],tm[1],tm[4],tm[1]+tm[4]+tm[7],0, tm[4], tm[7]],
                [tm[2]+tm[5]+tm[8],tm[8], tm[2],0,tm[2]+tm[5]+tm[8], tm[5], tm[8]],
                [tm[3]+tm[4]+tm[5],tm[3], tm[4], tm[4],tm[5],tm[3]+tm[4]+tm[5],0],
                [tm[6]+tm[7]+tm[8],tm[8],tm[6], tm[7],tm[8],0,tm[6]+tm[7]+tm[8]]
                ])


def matrix_AB(I_,V_):
    loc=np.arange(0,len(I_)).tolist()
    loc.pop(1)
    I_eta=I_[1, 1]
    I_etapsy=np.array([I_[1,i] for i in loc])#np.array([I_[1,0],I_[1, 2],I_[1, 3],I_[1, 4]])
    I_psy=np.array([[I_[i,0],I_[i, 2],I_[i, 3],I_[i, 4],I_[i, 5],I_[i, 6]] for i in loc])
    #np.array([[I_[0,0],I_[0, 2],I_[0, 3],I_[0, 4]],[I_[2,0],I_[2, 2],I_[2, 3],I_[2, 4]],[I_[3,0],I_[3, 2],I_[3, 3],I_[3, 4]],[I_[4,0],I_[4, 2],I_[4, 3],I_[4, 4]]])
    V_eta=V_[1, 1]
    V_etapsy=np.array([V_[1,i] for i in loc])#np.array([V_[1,0],V_[1, 2],V_[1, 3],V_[1, 4]])
    V_psy=np.array([[V_[i,0],V_[i, 2],V_[i, 3],V_[i, 4],V_[i, 5],V_[i, 6]] for i in loc])
    #np.array([[V_[0,0],V_[0, 2],V_[0, 3],V_[0, 4]],[V_[2,0],V_[2, 2],V_[2, 3],V_[2, 4]],[V_[3,0],V_[3, 2],V_[3, 3],V_[3, 4]],[V_[4,0],V_[4, 2],V_[4, 3],V_[4, 4]]])
    A=I_eta-I_etapsy.dot(lin.inv(I_psy)).dot(I_etapsy.reshape(len(I_)-1,1))
    B=V_eta-2*I_etapsy.dot(lin.inv(I_psy)).dot(V_etapsy.reshape(len(I_)-1,1))+I_etapsy.dot(lin.inv(I_psy)).dot(V_psy).dot(lin.inv(I_psy)).dot(I_etapsy.reshape(len(I_)-1,1))
    return A,B

def Estimate(mean_true,sample_size,data_type,cor_par,eta0):
    #mean_true=tm
    #sample_size=75
    #cor_par=gamma_param
    
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
        alpha, eta1, eta2, gamma1, gamma2, delta1, delta2=params
        f1,f2,f3=0,0,0
        for i in range(len(data)):
            f1=f1+alpha*data.at[i,'Yi11']+(alpha+eta1+gamma1)*data.at[i,'Yi21']+(alpha+eta2+gamma2)*data.at[i,'Yi31']-np.exp(alpha)-np.exp(alpha+eta1+gamma1)-np.exp(alpha+eta2+gamma2)
            f2=f2+(alpha+eta1+delta1)*data.at[i,'Yi12']+(alpha+eta2+gamma1+delta1)*data.at[i,'Yi22']+(alpha+gamma2+delta1)*data.at[i,'Yi32']-np.exp(alpha+eta1+delta1)-np.exp(alpha+eta2+gamma1+delta1)-np.exp(alpha+gamma2+delta1)
            f3=f3+(alpha+eta2+delta2)*data.at[i,'Yi13']+(alpha+gamma1+delta2)*data.at[i,'Yi23']+(alpha+eta1+gamma2+delta2)*data.at[i,'Yi33']-np.exp(alpha+eta2+delta2)-np.exp(alpha+gamma1+delta2)-np.exp(alpha+eta1+gamma2+delta2)
        return f1+f2+f3
    initial_guess = np.array([0.95, 0.372, 0.354, 0.236, 0.226, 0.206, 0.264])
    solve = scipy.optimize.minimize (lambda params: -L(params), initial_guess)
    MLE_opt=solve.x
    
    # MLE
    estimate= pd.DataFrame({'alpha_hat': MLE_opt[0], 
                            'eta1_hat': MLE_opt[1], 
                            'eta2_hat': MLE_opt[2], 
                            'gamma1_hat': MLE_opt[3], 
                            'gamma2_hat': MLE_opt[4],
                            'delta1_hat': MLE_opt[5],
                            'delta2_hat': MLE_opt[6]},index=[0])
    
    '''Mean'''   
    etimate_of_mean = pd.DataFrame({'Mean_11': np.exp(MLE_opt[0]), 
                                    'Mean_21': np.exp(MLE_opt[0] + MLE_opt[1] + MLE_opt[3]), 
                                    'Mean_31': np.exp(MLE_opt[0] + MLE_opt[2] + MLE_opt[4]),
                                    'Mean_12': np.exp(MLE_opt[0] + MLE_opt[1] + MLE_opt[5]), 
                                    'Mean_22': np.exp(MLE_opt[0] + MLE_opt[2] + MLE_opt[3] + MLE_opt[5]),
                                    'Mean_32': np.exp(MLE_opt[0] + MLE_opt[4] + MLE_opt[5]),
                                    'Mean_13': np.exp(MLE_opt[0] + MLE_opt[2] + MLE_opt[6]), 
                                    'Mean_23': np.exp(MLE_opt[0] + MLE_opt[3] + MLE_opt[6]),
                                    'Mean_33':  np.exp(MLE_opt[0] + MLE_opt[1] + MLE_opt[4] + MLE_opt[6])},index=[0])
    
    sm=list(etimate_of_mean.loc[0,:])    
    '''Estimate matrix I'''
    p1=(seq_size)/(seq_size*3)
    I_hat =p1*np.array([
        [sum(sm), sm[1]+sm[3]+sm[8], sm[2]+sm[4]+sm[6], sm[1]+sm[4]+sm[7], sm[2]+sm[5]+sm[8], sm[3]+sm[4]+sm[5], sm[6]+sm[7]+sm[8]],
        [sm[1]+sm[3]+sm[8],sm[1]+sm[3]+sm[8], 0, sm[1], sm[8], sm[3], sm[8] ],
        [sm[2]+sm[4]+sm[6],0, sm[2]+sm[4]+sm[6], sm[4], sm[2], sm[4], sm[6]],
        [sm[1]+sm[4]+sm[7],sm[1], sm[4],sm[1]+sm[4]+sm[7],0, sm[4], sm[7]],
        [sm[2]+sm[5]+sm[8],sm[8], sm[2],0,sm[2]+sm[5]+sm[8], sm[5], sm[8]],
        [sm[3]+sm[4]+sm[5],sm[3], sm[4], sm[4],sm[5],sm[3]+sm[4]+sm[5],0],
        [sm[6]+sm[7]+sm[8],sm[8], sm[6], sm[7],sm[8],0,sm[6]+sm[7]+sm[8]]
                      ])
    
    #Diagonal variance matrix
    Var=data.var(ddof=1)
    Cov=[ data.Yi11.cov(data.Yi21), data.Yi11.cov(data.Yi31),data.Yi21.cov(data.Yi31), 
          data.Yi12.cov(data.Yi22), data.Yi12.cov(data.Yi32),data.Yi22.cov(data.Yi32),
          data.Yi13.cov(data.Yi23), data.Yi13.cov(data.Yi33),data.Yi23.cov(data.Yi33)]

    V_hat=p1*np.array([
        [sum(Var)+2*sum(Cov), Var[1]+Var[3]+Var[8]+sum(Cov)-Cov[1]-Cov[5]-Cov[6], Var[2]+Var[4]+Var[6]+sum(Cov)-Cov[0]-Cov[4]-Cov[8], Var[1]+Var[4]+Var[7]+sum(Cov)-Cov[1]-Cov[4]-Cov[7], Var[2]+Var[5]+Var[8]+sum(Cov)-Cov[0]-Cov[3]-Cov[6], sum(Var[3:6])+2*sum(Cov[3:6]), sum(Var[6:9])+2*sum(Cov[6:9])],
        [Var[1]+Var[3]+Var[8]+sum(Cov)-Cov[1]-Cov[5]-Cov[6],Var[1]+Var[3]+Var[8], Cov[2]+Cov[3]+Cov[7], Var[1]+Cov[3]+Cov[8], Var[8]+Cov[2]+Cov[4], Var[3]+Cov[3]+Cov[5], Var[8]+Cov[7]+Cov[8]],
        [Var[2]+Var[4]+Var[6]+sum(Cov)-Cov[0]-Cov[4]-Cov[8],Cov[2]+Cov[3]+Cov[7],Var[2]+Var[4]+Var[6], Cov[2]+Var[4]+Cov[6], Var[2]+Cov[5]+Cov[7], Var[4]+Cov[3]+Cov[5], Var[6]+Cov[6]+Cov[8]],
        [Var[1]+Var[4]+Var[7]+sum(Cov)-Cov[1]-Cov[4]-Cov[7],Var[1]+Cov[3]+Cov[8],Cov[2]+Var[4]+Cov[6], Var[1]+Var[4]+Var[7], Cov[2]+Cov[5]+Cov[8], Var[4]+Cov[3]+Cov[5], Var[7]+Cov[6]+Cov[8]],
        [Var[2]+Var[5]+Var[8]+sum(Cov)-Cov[0]-Cov[3]-Cov[6],Var[8]+Cov[2]+Cov[4],Var[2]+Cov[5]+Cov[7],Cov[2]+Cov[5]+Cov[8],Var[2]+Var[5]+Var[8],Var[5]+Cov[4]+Cov[5], Var[8]+Cov[7]+Cov[8]],
        [sum(Var[3:6])+2*sum(Cov[3:6]),Var[3]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[5]+Cov[4]+Cov[5],sum(Var[3:6])+2*sum(Cov[3:6]),0],
        [sum(Var[6:9])+2*sum(Cov[6:9]),Var[8]+Cov[7]+Cov[8],Var[6]+Cov[6]+Cov[8],Var[7]+Cov[6]+Cov[8],Var[8]+Cov[7]+Cov[8],0,sum(Var[6:9])+2*sum(Cov[6:9])]
                    ])
    '''
    V_hat=p1*np.array([
        [sum(Var), Var[1]+Cov[0]+Cov[2]+Var[3]+Cov[3]+Cov[4]+Var[8]+Cov[7]+Cov[8], Var[2]+Cov[2]+Cov[3]+Var[4]+Cov[3]+Cov[4]+Var[6]+Cov[7]+Cov[8], Var[1]+Cov[0]+Cov[2]+Var[4]+Cov[3]+Cov[5]+Var[7]+Cov[6]+Cov[8], Var[2]+Cov[1]+Cov[2]+Var[5]+Cov[1]+Cov[5]+Var[8]+Cov[7]+Cov[8], sum(Var[3:6])+2*sum(Cov[3:6]), sum(Var[6:9])+2*sum(Cov[6:9])],
        [Var[1]+Cov[0]+Cov[2]+Var[3]+Cov[3]+Cov[4]+Var[8]+Cov[7]+Cov[8],Var[1]+Var[3]+Var[8], Cov[2]+Cov[3]+Cov[7], Var[1]+Cov[3]+Cov[8], Var[8]+Cov[2]+Cov[4], Var[3]+Cov[3]+Cov[5], Var[8]+Cov[7]+Cov[8]],
        [Var[2]+Cov[2]+Cov[3]+Var[4]+Cov[3]+Cov[4]+Var[6]+Cov[7]+Cov[8],Cov[2]+Cov[3]+Cov[7],Var[2]+Var[4]+Var[6], Cov[2]+Var[4]+Cov[6], Var[2]+Cov[5]+Cov[7], Var[4]+Cov[3]+Cov[5], Var[6]+Cov[6]+Cov[8]],
        [Var[1]+Cov[0]+Cov[2]+Var[4]+Cov[3]+Cov[5]+Var[7]+Cov[6]+Cov[8],Var[1]+Cov[3]+Cov[8],Cov[2]+Var[4]+Cov[6], Var[1]+Var[4]+Var[7], Cov[2]+Cov[5]+Cov[8], Var[4]+Cov[3]+Cov[5], Var[7]+Cov[6]+Cov[8]],
        [Var[2]+Cov[1]+Cov[2]+Var[5]+Cov[1]+Cov[5]+Var[8]+Cov[7]+Cov[8],Var[8]+Cov[2]+Cov[4],Var[2]+Cov[5]+Cov[7],Cov[2]+Cov[5]+Cov[8],Var[2]+Var[5]+Var[8],Var[5]+Cov[5]+Cov[8], Var[8]+Cov[7]+Cov[8]],
        [sum(Var[3:6])+2*sum(Cov[3:6]),Var[3]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[5]+Cov[5]+Cov[8],sum(Var[3:6])+2*sum(Cov[3:6]),0],
        [sum(Var[6:9])+2*sum(Cov[6:9]),Var[8]+Cov[7]+Cov[8],Var[6]+Cov[6]+Cov[8],Var[7]+Cov[6]+Cov[8],Var[8]+Cov[7]+Cov[8],0,sum(Var[6:9])+2*sum(Cov[6:9])]
                    ]) '''   
                    
    return estimate, I_hat,V_hat



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

I_ind=I_ind/sim_time
V_ind=V_ind/sim_time
I_cor=I_cor/sim_time
V_cor=V_cor/sim_time

np.set_printoptions(suppress=True,precision=5)
print('Independent Data, seq_size = ',seq_size)
print('True Value of parameters\n alpha：%f, eta1：%f, eta1：%f,gamma1：%f, gamma2：%f, delta1：%f, delta1：%f' %(1.0, 0.55, 0.45, 0.21, 0.26, 0.13, 0.17))
print('MLE\n',mle_ind.mean())
print('Sample variance of estimates\n',3*seq_size*mle_ind.var(ddof=1))

print('Inverse of matrix I\n',lin.inv(I_ind))
print('S\n',mle_ind.cov(ddof=1))
print('N*S\n',3*seq_size*mle_ind.cov(ddof=1))
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

ns=3*seq_size*mle_ind.cov(ddof=1)
# find p-value for two-tailed test
scipy.stats.t.sf(abs(np.mean(LR_ind_na)), df=1)*2

print('Correlated Data (alpha =10 ,beta = 0.1), seq_size = ',seq_size)
print('True Value of parameters\n alpha：%f, eta1：%f, eta1：%f,gamma1：%f, gamma2：%f, delta1：%f, delta1：%f' %(1.0, 0.55, 0.45, 0.21, 0.26, 0.13, 0.17))
print('MLE\n',mle_cor.mean())
print('Sample variance of estimates\n',3*seq_size*mle_cor.var(ddof=1))

print('Inverse of matrix I\n',lin.inv(I_cor))
print('S\n',mle_cor.cov(ddof=1))
print('N*S\n',3*seq_size*mle_cor.cov(ddof=1))
print('inv(I)*V*inv(I)\n',lin.inv(I_cor).dot(V_cor).dot(lin.inv(I_cor)))


print('Estimate of matrix I\n',I_cor)
print('Estimate of matrix V\n',V_cor)
ns_=3*seq_size*mle_cor.cov(ddof=1)
print('LR naive',np.mean(LR_cor_na))
print('LR robust',np.mean(LR_cor_rb))
print('Wald naive',np.mean(Wald_cor_na))
print('Wald robust',np.mean(Wald_cor_rb))

print('LR naive p-value',scipy.stats.t.sf(abs(np.mean(LR_cor_na)), df=1)*2)
print('LR robust p-value',scipy.stats.t.sf(abs(np.mean(LR_cor_rb)), df=1)*2)
print('Wald naive p-value',scipy.stats.t.sf(abs(np.mean(Wald_cor_na)), df=1)*2)
print('Wald robust p-value',scipy.stats.t.sf(abs(np.mean(Wald_cor_rb)), df=1)*2)

dict_ind={'MLE':mle_ind.mean(),
          'N*S':3*seq_size*mle_ind.var(ddof=1),
          'S':mle_ind.var(ddof=1),
          'inv I':lin.inv(I_ind),
          'N*Cov':3*seq_size*mle_ind.cov(ddof=1),
          'invI.V.invI':lin.inv(I_ind).dot(V_ind).dot(lin.inv(I_ind)),
          'I':I_ind,
          'V':V_ind}
dict_cor={'MLE':mle_cor.mean(),
          'N*S':3*seq_size*mle_cor.var(ddof=1),
          'S':mle_cor.var(ddof=1),
          'inv I':lin.inv(I_cor),
          'N*Cov':3*seq_size*mle_cor.cov(ddof=1),
          'invI.V.invI':lin.inv(I_cor).dot(V_cor).dot(lin.inv(I_cor)),
          'I':I_cor,
          'V':V_cor}