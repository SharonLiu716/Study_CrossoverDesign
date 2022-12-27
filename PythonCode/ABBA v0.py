# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 09:58:46 2022

@author: cherl
"""
import sys
import scipy
import math
import datetime
from math import factorial
import numpy as np
import pandas as pd
import scipy.stats
from scipy.optimize import minimize, rosen, rosen_der
from scipy.stats import poisson
import numpy.linalg as lin
from scipy.optimize import fsolve
from sympy import *
from autograd import grad
import autograd.numpy as np
from autograd import grad, jacobian, hessian
sim_time,seq_size,cor_param,cros_type=1000,[25,50,100,200],0.1,'ABBA'
tao,eta,gamma,delta=1.0, 0.7, 0.3, 0.2
params = np.array([tao,eta,gamma,delta])
covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
mean_true=np.exp(np.dot(params.transpose(),covariate))
pi=(seq_size)/(seq_size*2)
np.array([[covariate[i]*covariate[j] for j in range(len(cros_type))] for i in range(len(cros_type))])
I=pi*np.array([[np.dot(mean_true,covariate[i]*covariate[j]) for j in range(len(cros_type))] for i in range(len(cros_type))])

#data=pd.DataFrame(np.random.poisson(lam=mean_true, size=(seq_size, len(cros_type))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
#data=np.random.poisson(lam=mean_true, size=(seq_size, len(cros_type))).astype(float)

#MLE
np.random.seed(55043)
MLE_optim=pd.DataFrame()
I_hessian=0
for i in range(sim_time):
    data_=np.random.poisson(lam=mean_true, size=(200, len(cros_type))).astype(float)
    def logL(params):
        tao_, eta_, gamma_,delta_=params.astype(float)
        factorize=sum(sum(np.log(np.array([[math.factorial(data_[i,j]) for i in range(len(data_))] for j in range(len(cros_type))]).astype(float))))
        logLL=np.sum(tao_*data_[:,0].astype(float)-np.exp(tao_)+(tao_+eta_+gamma_)*data_[:,1].astype(float)-np.exp((tao_+eta_+gamma_))+(tao_+eta_+delta_)*data_[:,2].astype(float)-np.exp((tao_+eta_+delta_))+(tao_+gamma_+delta_)*data_[:,3].astype(float) -np.exp((tao_+gamma_+delta_)))-factorize       
        return logLL
    res = minimize(fun=lambda par: -logL(par),x0=np.array([0.95, 0.65, 0.25, 0.15], dtype=float))
    mle= pd.DataFrame({'tao_hat': res.x[0],'eta_hat': res.x[1],'gamma_hat':res.x[2],'delta_hat': res.x[3]},index=[0])   
    MLE_optim = MLE_optim.append(mle,ignore_index=True)
    hessian_ = hessian(logL)
    I_h = -hessian_(res.x)
    I_hessian+=I_h

MLE_optim.mean()
MLE_optim.cov()*400
I_hessian/(1000*400)

def Estimate(cors_type,mean_true,sample_size,data_type,cor_par,eta0):

    if data_type == 'ind': 
        data=pd.DataFrame(np.random.poisson(lam=mean_true, size=(sample_size, len('ABBA'))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    else:   
        
        Mu_cor = pd.DataFrame(columns = ['Mu11', 'Mu12', 'Mu21', 'Mu22'])
        for i in range(0,len(mean_true),2):            
            nu=np.random.gamma(1/cor_par,cor_par,sample_size)
            for (ix,param) in enumerate(mean_true[i:i + 2]):
                Mu_cor[Mu_cor.columns[ix+i]]=pd.DataFrame(np.multiply(np.array(param).repeat(sample_size), nu).T)
                

        data=pd.DataFrame(columns = ['Yi11', 'Yi12', 'Yi21', 'Yi22'])
        for (idx,mu) in enumerate(Mu_cor.columns):
            data[data.columns[idx]]=pd.DataFrame([np.random.poisson(p) for p in Mu_cor[mu]])
   
    #MLE
    tao_hat=np.log(data['Yi11'].mean())
    eta_hat=0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi21'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi22'].sum()))
    gamma_hat=0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi22'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi21'].sum()))
    delta_hat=0.5*(np.log(data['Yi21'].sum())+np.log(data['Yi22'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi12'].sum()))
    # MLE
    estimate= pd.DataFrame({'alpha_hat': tao_hat, 
                            'eta_hat': eta_hat, 
                            'gamma_hat':gamma_hat, 
                            'delta_hat': delta_hat},index=[0])    
    
    '''Mean'''
    mle = np.array([tao_hat,eta_hat,gamma_hat,delta_hat])
    covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
    etimate_of_mean=np.exp(np.dot(mle.transpose(),covariate))
    
    '''Estimate matrix I''' 
    #sm=list(etimate_of_mean.loc[0,:])#data.mean()
    pi=sample_size/(sample_size*2)
    I_hat=pi*np.array([[np.dot(etimate_of_mean,covariate[i]*covariate[j]) for j in range(4)] for i in range(4)])
    score=[np.dot((data-etimate_of_mean).to_numpy(),covariate[i].transpose()) for i in range(len(covariate))]
    V_hat=np.array([[sum(score[i]*score[j]) for j in range(len(covariate))] for i in range(4)])/(sample_size*2)
    Var=data.var(ddof=1)
    Cov=[ data.Yi11.cov(data.Yi12), data.Yi21.cov(data.Yi22)]
    V = pi*np.array([
                    [sum(Var)+2*sum(Cov), Var[1]+Var[2]+sum(Cov), Var[1]+Var[3]+sum(Cov), sum(Var[2:4])+2*Cov[1]],
                    [Var[1]+Var[2]+sum(Cov),Var[1]+Var[2],Var[1]+Cov[1],Var[2]+Cov[1]],
                    [Var[1]+Var[3]+sum(Cov),Var[1]+Cov[1],Var[1]+Var[3],Var[3]+Cov[1]],
                    [sum(Var[2:4])+2*Cov[1],Var[2]+Cov[1],Var[3]+Cov[1],sum(Var[2:4])+2*Cov[1]]
                     ])
    return estimate, I_hat, V_hat, V



def Simulation(runtime,design_type,pv,seqsize):
    np.set_printoptions(suppress=True,precision=4)
    mle_ind,mle_cor=pd.DataFrame(),pd.DataFrame()
    I_ind, V_indm,I_cor, V_corm = 0, 0, 0, 0
    V_indc, V_corc = 0, 0
    for i in range(runtime):    
        #independent
        mle_i, I_i, V_i_matrix, V_i_close=Estimate(cors_type=design_type,mean_true=pv,sample_size=seqsize,data_type='ind',cor_par=0.1,eta0=0)
        #mle_i, I_i, V_i, LR_na, LR_rb, Wald_na, Wald_rb=Estimate(mean_true=tm, sample_size=seq_size,data_type='ind',cor_par=gamma_param,eta0=0)
        mle_ind = mle_ind.append(mle_i,ignore_index=True)
        #mu_ind = mu_ind.append(mu_i,ignore_index=True)
        I_ind+=I_i
        V_indm+=V_i_matrix    
        V_indc+=V_i_close
        '''
        #correlated
        mle_i, I_i,  V_i_matrix, V_i_close, V_i_hessian=Estimate(cors_type='ABBA',mean_true=mean_true,sample_size=seq_size,data_type='cor',cor_par=cor_param,eta0=0)
        mle_cor = mle_cor.append(mle_i,ignore_index=True)
        I_cor+=I_i
        V_corm+=V_i_matrix
        V_corc+=V_i_close
        V_corh+=V_i_hessian'''
        
    I_ind=I_ind/runtime
    V_indm=V_indm/runtime
    V_indc=V_indc/runtime
    
    '''
    I_cor=I_cor/sim_time
    V_corm=V_corm/sim_time
    V_corc=V_corc/sim_time
    V_corh=V_corh/sim_time'''
    
    
    return mle_ind.mean(),I_ind,V_indm,V_indc,mle_ind.cov(ddof=1)#, mle_cor.mean(),I_cor,V_corm,V_corc,V_corh,mle_cor.cov(ddof=1)
#mle_ind.mean(),lin.inv(I_ind).dot(V_ind).dot(lin.inv(I_ind)),mle_ind.cov(ddof=1),mle_cor.mean(),lin.inv(I_cor).dot(V_cor).dot(lin.inv(I_cor)),mle_cor.cov(ddof=1)

np.random.seed(980716)
output_path="C:/Github/Study_CrossoverDesign/SimOutput"
today=str(datetime.date.today().month)+'0'+str(datetime.date.today().day)
writer=pd.ExcelWriter(output_path+'/{}_{}.xlsx'.format("ABBA",today), engine="openpyxl")
num_seq= 2
    
for seq in seq_size:        
    MLE_ind,I_ind_,V_ind_m,V_ind_c,cov_ind=Simulation( runtime=sim_time, design_type="ABBA", pv = mean_true, seqsize=seq)#,MLE_cor,I_cor_,V_cor_m,V_cor_c,V_cor_h,cov_cor
    sys.stdout.write('\rReading : '+"ABBA"+str(seq))
    np.set_printoptions(suppress=True,precision=5)
    MLE_ind.to_excel(writer, sheet_name="MLE_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')        
    pd.DataFrame(I_ind_).to_excel(writer, sheet_name="I_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
    pd.DataFrame(V_ind_m).to_excel(writer, sheet_name="V_ind_matrix"+str(seq), engine='openpyxl', encoding='utf_8_sig')
    pd.DataFrame(V_ind_c).to_excel(writer, sheet_name="V_ind_closeform"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        #pd.DataFrame(V_ind_h).to_excel(writer, sheet_name="V_ind_hessian"+str(seq), engine='openpyxl', encoding='utf_8_sig')
    pd.DataFrame(lin.inv(I_ind_)).to_excel(writer, sheet_name="inv(I)_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        #pd.DataFrame(lin.inv(I_ind_).dot(V_ind_).dot(lin.inv(I_ind_))).to_excel(writer, sheet_name="invI_V_inv_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
    pd.DataFrame((num_seq*seq*cov_ind)).to_excel(writer, sheet_name="NS_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
    writer.save()
writer.close()
    
        '''
        MLE_cor.to_excel(writer, sheet_name="MLE_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')        
        pd.DataFrame(I_cor_).to_excel(writer, sheet_name="I_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        #pd.DataFrame(V_cor_).to_excel(writer, sheet_name="V_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        #pd.DataFrame(lin.inv(I_cor_)).to_excel(writer, sheet_name="inv(I)_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame(lin.inv(I_cor_).dot(V_cor_).dot(lin.inv(I_cor_))).to_excel(writer, sheet_name="invI_V_inv_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame((num_seq*seq*cov_cor)).to_excel(writer, sheet_name="NS_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')'''
    






'''
用來驗證logL有沒有算錯 但optim的值不正確!
def logL(params,data):
    tao_, eta_, gamma_,delta_=params
    params = np.array([tao_, eta_, gamma_,delta_])
    covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
    logmean=np.dot(params.transpose(),covariate)
    factorize=pd.DataFrame(np.array([[math.factorial(data.at[i,col]) for i in range(len(data))] for col in data.columns]).T.tolist(),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
    logL=np.sum((np.log(mean_true)*data-mean_true-np.log(factorize))).sum()
    
    return -logL'''