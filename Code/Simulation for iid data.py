# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 11:30:51 2022

@author: a7086
"""

import scipy
import math
from math import factorial
import numpy as np
import pandas as pd
from scipy.optimize import minimize 
import numpy.linalg as lin


''' 假設樣本之間為獨立生成的Poisson分配，求其最大概似估計量、費雪資訊矩陣、變異數矩陣
Input：
    -tm:用參數真值算出來的樣本平均數，用以生成獨立卜瓦松樣本 Yijk
    -sample_size：Yijk的樣本數(總樣本數為sample_sizs*4，每個序列的樣本數為sample_size*2)
Output：
    -estimate(dataframe):MLE
    -etimate_of_mean(dataframe)：mean of yij
    -I_hat：Fisher information matrix
    -V_hat：variance matrix 
Validate：
    1. V=I?
    2.sample variance of MLE = diag(I**(-1))
'''

#AB|BA、ABB|BAA、AAB|ABA|BAA、ABC|BCA|CAB、BAC|ACB|BCA、BBA|ACB|CAC
design='ABC|BCA|CAB'
design.split('|')

 

def true_value(design):
    #判斷交叉設計的序列數量、期間個數、藥劑個數
    paras_count=[len(design.split('|')), len(design.split('|')[0]), len([letter for letter, count in Counter(design).items()])-1]
    if all(element == 3 for element in paras_count):
        tv = list(map(int,input("\nEnter the value of alpha、eta1、eta2、gamma1、gamma2、delta1、delta2 by space : ").strip().split()))[:n]
  
print("\nList is - ", a)
    #序列個數
    if len(design.split('|')) == 2:
        delta=np.float(input('delta :'))
    else:
        delta1=np.float(input('delta1:'))
        delta2=np.float(input('delta2:'))
    #期間個數
    if len(design.split('|')[0]) == 2:
        gamma=np.float(input('gamma :'))
    else:
        gamma1=np.float(input('gamma1:'))
        gamma2=np.float(input('gamma2:'))
    #藥劑個數
    from collections import Counter
    treatment=[letter for letter, count in Counter(design).items()]
    if len([letter for letter, count in Counter(design).items()])-1 == 2:
        eta=np.float(input('delta :'))
    else:
        eta1=np.float(input('eta1:'))
        eta2=np.float(input('eta2:'))



class Independent_Simulation():
    '''def __init__
    Parameters
    ----------
    data_gen_type:independent/correlated
    design : type of cross over desing
        -type:ABBA、ABBBAA、AABABABAA、ABCBCACAB、BACACBBCA、BBAACBCAC
    value_of_tm: 根據design的參數個數給值
        用input&if-else給值
    
    Returns
    -------
    MLE、sample_mean、I_hat、V_hat
    
    '''
    
    def __init__(self, data_gen_type='independent', design='ABBA'):
        # data_version: original: monthly raw data from RQ, updated: 經過欄位縮短&變換
        self.data_gen_type = data_gen_type
        self.design = design
        self.data=self.data_generate(self.data_gen_type, self.mean_true, self.sample_size)
        self.mean_true=self.true_mean()            
        self.mapping_dict = pd.read_excel(mapping_dict, sheet_name=None) #sheetname=none可以匯入所有sheet，但以dict呈現
        self.merge_columns = merge_columns
        self.action = action
        self.version = version
        self.key1,self.key2 = self.PartKey(self.mapping_dict['Part Type'])        
        self.updated=self.Merge(self.DF,self.mapping_dict['Child'],self.mapping_dict['Part Type'],self.mapping_dict['Defect Type'],self.mapping_dict['Product Type'],self.merge_columns)

    def true_value(self,self.design):
        #判斷交叉設計的序列數量、藥劑個數、期間個數
        if len(self.design)/2 == 2:
            
            tv=np.float(input('parameter value for ABBA:'))
        elif self.design=='ABBBAA':
            tv=np.float(input('parameter value for ABBA:'))
        elif self.design=='ABBBAA':
            tv=np.float(input('parameter value for ABBA:'))
        
            
    def data_generate(self,data_type,mean_true,sample_size):
        if data_type=='independent':
            #generate data(以下以ABBA為例，需修改成通用格式)
            data = {'Yi11': np.random.poisson(lam=mean_true[0], size=sample_size),
                    'Yi21': np.random.poisson(lam=mean_true[1], size=sample_size),
                    'Yi12': np.random.poisson(lam=mean_true[2], size=sample_size),  
                    'Yi22': np.random.poisson(lam=mean_true[3], size=sample_size)}
            data = pd.DataFrame(data)
        else:
            #generate correlated data
            Mu = pd.DataFrame(columns = ['Mu11', 'Mu21', 'Mu12', 'Mu22'])
            data=pd.DataFrame(columns = ['Yi11', 'Yi21', 'Yi12', 'Yi22'])
            for (idx,param) in enumerate(mean_true): 
                
                Mu[Mu.columns[idx]]=pd.DataFrame([param*np.random.gamma(shape=1, scale=1,size=sample_size) ]).transpose()
                data[data.columns[idx]]=[ np.random.poisson(p) for p in Mu[Mu.columns[idx]]]


#參數真值
alpha, eta1, eta2, gamma1, gamma2, delta1, delta2 = 1.0, 0.32, 0.37, 1.21, 1.27, 0.22, 0.28
#true mean
tm=[np.exp(alpha), np.exp(alpha + eta1+ gamma1), np.exp(alpha + eta2 + gamma2),
    np.exp(alpha + eta1 + delta1), np.exp(alpha + eta2 + gamma1 + delta1),np.exp(alpha + gamma2 + delta1),
    np.exp(alpha + eta2 +delta2), np.exp(alpha + gamma1 + delta2), np.exp(alpha + eta1 + gamma2 + delta2)]
#費雪資訊矩陣(真實值)
pi=50/150
I = pi*np.array([  [sum(tm), tm[1]+tm[3]+tm[8], tm[2]+tm[4]+tm[6], tm[1]+tm[4]+tm[7], tm[2]+tm[5]+tm[8], tm[3]+tm[4]+tm[5], tm[6]+tm[7]+tm[8]],
                [tm[1]+tm[3]+tm[8],tm[1]+tm[3]+tm[8], 0, tm[1], tm[8], tm[3], tm[8] ],
                [tm[2]+tm[4]+tm[6],0, tm[2]+tm[4]+tm[6], tm[4], tm[2], tm[4], tm[6]],
                [tm[1]+tm[4]+tm[7],tm[1],tm[4],tm[1]+tm[4]+tm[7],0, tm[4], tm[7]],
                [tm[2]+tm[5]+tm[8],tm[8], tm[2],0,tm[2]+tm[5]+tm[8], tm[5], tm[8]],
                [tm[3]+tm[4]+tm[5],tm[3], tm[4], tm[4],tm[5],tm[3]+tm[4]+tm[5],0],
                [tm[6]+tm[7]+tm[8],tm[8],tm[6], tm[7],tm[8],0,tm[6]+tm[7]+tm[8]]
                ])

def var(x,mean_x): # * removed
    return sum((x - mean_x)**2)/len(x)
def cov(x, y,mean_x,mean_y): 
    return sum((x - mean_x) * (y - mean_y)) / len(x) 

def Estimate(mean_true,sample_size):
    #mean_true=tm
    #sample_size=25
    
    
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
        alpha, eta1, eta2, gamma1, gamma2, delta1, delta2=params
        f1,f2,f3=0,0,0
        for i in range(len(data)):
            f1=f1+alpha*data.at[i,'Yi11']+(alpha+eta1+gamma1)*data.at[i,'Yi21']+(alpha+eta2+gamma2)*data.at[i,'Yi31']-np.exp(alpha)-np.exp(alpha+eta1+gamma1)-np.exp(alpha+eta2+gamma2)
            f2=f2+(alpha+eta1+delta1)*data.at[i,'Yi12']+(alpha+eta2+gamma1+delta1)*data.at[i,'Yi22']+(alpha+gamma2+delta1)*data.at[i,'Yi32']-np.exp(alpha+eta1+delta1)-np.exp(alpha+eta2+gamma1+delta1)-np.exp(alpha+gamma2+delta1)
            f3=f3+(alpha+eta2+delta2)*data.at[i,'Yi13']+(alpha+gamma1+delta2)*data.at[i,'Yi23']+(alpha+eta1+gamma2+delta2)*data.at[i,'Yi33']-np.exp(alpha+eta2+delta2)-np.exp(alpha+gamma1+delta2)-np.exp(alpha+eta1+gamma2+delta2)
        return f1+f2+f3
    initial_guess = np.array([0.95, 0.302, 0.354, 1.198, 1.258, 0.206, 0.264])
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
    p1=50/150
    I_hat =p1*np.array([
        [sum(sm), sm[1]+sm[3]+sm[8], sm[2]+sm[4]+sm[6], sm[1]+sm[4]+sm[7], sm[2]+sm[5]+sm[8], sm[3]+sm[4]+sm[5], sm[6]+sm[7]+sm[8]],
        [sm[1]+sm[3]+sm[8],sm[1]+sm[3]+sm[8], 0, sm[1], sm[8], sm[3], sm[8] ],
        [sm[2]+sm[4]+sm[6],0, sm[2]+sm[4]+sm[6], sm[4], sm[2], sm[4], sm[6]],
        [sm[1]+sm[4]+sm[7],sm[1], sm[4],sm[1]+sm[4]+sm[7],0, sm[4], sm[7]],
        [sm[2]+sm[5]+sm[8],sm[8], sm[2],0,sm[2]+sm[5]+sm[8], sm[5], sm[8]],
        [sm[3]+sm[4]+sm[5],sm[3], sm[4], sm[4],sm[5],sm[3]+sm[4]+sm[5],0],
        [sm[6]+sm[7]+sm[8],sm[8], sm[6], sm[7],sm[8],0,sm[6]+sm[7]+sm[8]]
                      ])
    
    '''
    [sum(sm), sm[1]+sm[3]+sm[8], sm[2]+sm[4]+sm[6], sm[1]+sm[4]+sm[7], sm[2]+sm[5]+sm[8], sm[3]+sm[4]+sm[5], sm[6]+sm[7]+sm[8]],
    [sm[1]+sm[3]+sm[8],sm[1]+sm[3]+sm[8], 0, sm[1], sm[8], sm[3], sm[8] ],
    [sm[2]+sm[4]+sm[6],0, sm[2]+sm[4]+sm[6], sm[4], sm[2], sm[4], sm[6]],
    [sm[1]+sm[4]+sm[7],sm[1], sm[4],sm[1]+sm[4]+sm[7],0, sm[4], sm[7]],
    [sm[2]+sm[5]+sm[8],sm[8], sm[2],0,sm[2]+sm[5]+sm[8], sm[5], sm[8]],
    [sm[3]+sm[4]+sm[5],sm[3], sm[4], sm[4],sm[5],sm[3]+sm[4]+sm[5],0],
    [sm[6]+sm[7]+sm[8],sm[8], sm[6], sm[7],sm[8],0,sm[6]+sm[7]+sm[8]]
    '''
    #Diagonal variance matrix
    
    Var=[var(data[res],sm[idx]) for (idx,res) in enumerate(data.columns)]
    Cov=[cov(data['Yi11'],data['Yi21'],sm[0],sm[1]), cov(data['Yi11'],data['Yi31'],sm[0],sm[2]), cov(data['Yi21'],data['Yi31'],sm[1],sm[2]),
         cov(data['Yi12'],data['Yi22'],sm[3],sm[4]), cov(data['Yi12'],data['Yi32'],sm[3],sm[5]), cov(data['Yi22'],data['Yi32'],sm[4],sm[5]),
         cov(data['Yi13'],data['Yi23'],sm[6],sm[7]), cov(data['Yi13'],data['Yi33'],sm[6],sm[8]), cov(data['Yi23'],data['Yi33'],sm[7],sm[8])]

    V_hat=p1*np.array([
        [sum(Var), Var[1]+Cov[0]+Cov[2]+Var[3]+Cov[3]+Cov[4]+Var[8]+Cov[7]+Cov[8], Var[2]+Cov[2]+Cov[3]+Var[4]+Cov[3]+Cov[4]+Var[6]+Cov[7]+Cov[8], Var[1]+Cov[0]+Cov[2]+Var[4]+Cov[3]+Cov[5]+Var[7]+Cov[6]+Cov[8], Var[2]+Cov[1]+Cov[2]+Var[5]+Cov[1]+Cov[5]+Var[8]+Cov[7]+Cov[8], sum(Var[3:6])+2*sum(Cov[3:6]), sum(Var[6:9])+2*sum(Cov[6:9])],
        [Var[1]+Cov[0]+Cov[2]+Var[3]+Cov[3]+Cov[4]+Var[8]+Cov[7]+Cov[8],Var[1]+Var[3]+Var[8], Cov[2]+Cov[3]+Cov[7], Var[1]+Cov[3]+Cov[8], Var[8]+Cov[2]+Cov[4], Var[3]+Cov[3]+Cov[5], Var[8]+Cov[7]+Cov[8]],
        [Var[2]+Cov[2]+Cov[3]+Var[4]+Cov[3]+Cov[4]+Var[6]+Cov[7]+Cov[8],Cov[2]+Cov[3]+Cov[7],Var[2]+Var[4]+Var[6], Cov[2]+Var[4]+Cov[6], Var[2]+Cov[5]+Cov[7], Var[4]+Cov[3]+Cov[5], Var[6]+Cov[6]+Cov[8]],
        [Var[1]+Cov[0]+Cov[2]+Var[4]+Cov[3]+Cov[5]+Var[7]+Cov[6]+Cov[8],Var[1]+Cov[3]+Cov[8],Cov[2]+Var[4]+Cov[6], Var[1]+Var[4]+Var[7], Cov[2]+Cov[5]+Cov[8], Var[4]+Cov[3]+Cov[5], Var[7]+Cov[6]+Cov[8]],
        [Var[2]+Cov[1]+Cov[2]+Var[5]+Cov[1]+Cov[5]+Var[8]+Cov[7]+Cov[8],Var[8]+Cov[2]+Cov[4],Var[2]+Cov[5]+Cov[7],Cov[2]+Cov[5]+Cov[8],Var[2]+Var[5]+Var[8],Var[5]+Cov[5]+Cov[8], Var[8]+Cov[7]+Cov[8]],
        [sum(Var[3:6])+2*sum(Cov[3:6]),Var[3]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[5]+Cov[5]+Cov[8],sum(Var[3:6])+2*sum(Cov[3:6]),0],
        [sum(Var[6:9])+2*sum(Cov[6:9]),Var[8]+Cov[7]+Cov[8],Var[6]+Cov[6]+Cov[8],Var[7]+Cov[6]+Cov[8],Var[8]+Cov[7]+Cov[8],0,sum(Var[6:9])+2*sum(Cov[6:9])]
                    ])
        
                    
    return estimate, etimate_of_mean, I_hat,V_hat

'''
Var[0]=y11
Var[1]=y21
Var[2]=y31
Var[3]=y12
Var[4]=y22
Var[5]=y32
Var[6]=y13
Var[7]=y23
Var[8]=y33



[sum(Var), Var[1]+Cov[0]+Cov[2]+Var[3]+Cov[3]+Cov[4]+Var[8]+Cov[7]+Cov[8], Var[2]+Cov[2]+Cov[3]+Var[4]+Cov[3]+Cov[4]+Var[6]+Cov[7]+Cov[8], Var[1]+Cov[0]+Cov[2]+Var[4]+Cov[3]+Cov[5]+Var[7]+Cov[6]+Cov[8], Var[2]+Cov[1]+Cov[2]+Var[5]+Cov[1]+Cov[5]+Var[8]+Cov[7]+Cov[8], sum(Var[3:6])+2*sum(Cov[3:6]), sum(Var[6:9])+2*sum(Cov[6:9])],
[Var[1]+Cov[0]+Cov[2]+Var[3]+Cov[3]+Cov[4]+Var[8]+Cov[7]+Cov[8],Var[1]+Var[3]+Var[8], Cov[2]+Cov[3]+Cov[7], Var[1]+Cov[3]+Cov[8], Var[8]+Cov[2]+Cov[4], Var[3]+Cov[3]+Cov[5], Var[8]+Cov[7]+Cov[8]],
[Var[2]+Cov[2]+Cov[3]+Var[4]+Cov[3]+Cov[4]+Var[6]+Cov[7]+Cov[8],Cov[2]+Cov[3]+Cov[7],Var[2]+Var[4]+Var[6], Cov[2]+Var[4]+Cov[6], Var[2]+Cov[5]+Cov[7], Var[4]+Cov[3]+Cov[5], Var[6]+Cov[6]+Cov[8]],
[Var[1]+Cov[0]+Cov[2]+Var[4]+Cov[3]+Cov[5]+Var[7]+Cov[6]+Cov[8],Var[1]+Cov[3]+Cov[8],Cov[2]+Var[4]+Cov[6], Var[1]+Var[4]+Var[7], Cov[2]+Cov[5]+Cov[8], Var[4]+Cov[3]+Cov[5], Var[7]+Cov[6]+Cov[8]],
[Var[2]+Cov[1]+Cov[2]+Var[5]+Cov[1]+Cov[5]+Var[8]+Cov[7]+Cov[8],Var[8]+Cov[2]+Cov[4],Var[2]+Cov[5]+Cov[7],Cov[2]+Cov[5]+Cov[8],Var[2]+Var[5]+Var[8],Var[5]+Cov[5]+Cov[8], Var[8]+Cov[7]+Cov[8]],
[sum(Var[3:6])+2*sum(Cov[3:6]),Var[3]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[4]+Cov[3]+Cov[5],Var[5]+Cov[5]+Cov[8],sum(Var[3:6])+2*sum(Cov[3:6]),0],
[sum(Var[6:9])+2*sum(Cov[6:9]),Var[8]+Cov[7]+Cov[8],Var[6]+Cov[6]+Cov[8],Var[7]+Cov[6]+Cov[8],Var[8]+Cov[7]+Cov[8],0,sum(Var[6:9])+2*sum(Cov[6:9])]
                ])

'''
np.random.seed(980716)
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
np.set_printoptions(suppress=True,precision=4)
print('True Value of parameters\n alpha：%f, eta1：%f, eta1：%f,gamma1：%f, gamma2：%f, delta1：%f, delta1：%f' %(1.0, 0.32, 0.37, 1.21, 1.27, 0.22, 0.28))
print('MLE\n',mle.mean())
print('Sample variance of parameters\n',75*mle.cov())
print('Inverse of matrix I\n',lin.inv(I))
print('Inverse of matrix I_hat\n',lin.inv(I_))
print('True Value of matrix I\n',I)
print('Estimate of matrix I\n',I_)
print('Estimate of matrix V\n',V_)