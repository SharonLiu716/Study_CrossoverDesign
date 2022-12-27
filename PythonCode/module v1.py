# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:43:59 2022

@author: cherl
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 14:49:36 2022

@author: cherl
"""

'''
class CrossoverDesign()
    def __init__
        Parameters
        ----------
        cros_type : 交叉設計方式 e.g 'ABBA'
        data_type : 'ind'/'cor'
        cor_par：生成相關性資料的Gamma參數值(beta)
       mean_true:參數真實值，以list形式輸入，根據參數數量決定list長度
       sample_size；每個序列的樣本數nj
       eta0 = 欲檢定的的參數值
    
    Returns
    -------
    更新完的dataframe
    

'''
import sys
import scipy
import math
import openpyxl
import numpy as np
import pandas as pd
import datetime
import scipy.stats
from scipy.optimize import minimize 
import numpy.linalg as lin
from scipy.optimize import fsolve

class CrossoverDesign():
    def __init__(self, cros_type,data_type,cor_par,params_value,sample_size,eta0):
        # 基本input變數
        self.cros_type = cros_type                
        self.data_type = data_type
        self.cor_par = cor_par     
        self.params_value = params_value
        self.sample_size=sample_size
        self.eta0=eta0
        #由以下函數生成的變數
        #self.pi=self.PI(self.cros_type,self.sample_size)
        self.params,self.covariate=self.Link(self.cros_type,self.params_value)
        self.mean_true=self.Mean(self.params,self.covariate)
        self.data= self.Data_Generate(self.cros_type,self.mean_true,self.sample_size,self.data_type,self.cor_par,self.eta0)
        self.estimate = self.MLE(self.cros_type, self.data)
        self.mean_estimate=self.Mean(self.estimate.to_numpy()[0],self.covariate)
        self.I=self.MatrixI(self.cros_type, self.params, self.covariate, self.mean_true, self.sample_size)
        self.V_hat=self.MatrixV(self.cros_type, self.estimate.to_numpy()[0], self.covariate, self.mean_estimate, self.data, self.sample_size)
        self.I_hat=self.MatrixI(self.cros_type, self.estimate.to_numpy()[0], self.covariate, self.mean_estimate, self.sample_size)
        
    #資料生成
    def Data_Generate(self,cros_type,mean_true,sample_size,data_type,cor_par,eta0):
    
        #num_seq:序列數量決定我們要生成幾組nui
        num_seq= 3 if len(cros_type)==9 else 2   
        if data_type == 'ind': 
            data=pd.DataFrame(np.random.poisson(lam=mean_true, size=(sample_size, len(cros_type))))
        else:
            Mu_cor = pd.DataFrame(columns = [*range(0,len(cros_type))])
            for i in range(0,len(mean_true),num_seq):            
                nu=np.random.gamma(1/cor_par,cor_par,sample_size)
                for (ix,param) in enumerate(mean_true[i:i + num_seq]):
                    Mu_cor[Mu_cor.columns[ix+i]]=pd.DataFrame(np.multiply(np.array(param).repeat(sample_size), nu).T)
            data=pd.DataFrame(columns = [*range(0,len(cros_type))])
            for (idx,mu) in enumerate(Mu_cor.columns):
                data[data.columns[idx]]=pd.DataFrame([np.random.poisson(p) for p in Mu_cor[mu]])
        #rename data column, mu_column,pi,covariance        
        if len(cros_type)==4:
            data.columns=['Yi11', 'Yi12','Yi21', 'Yi22']
            
        elif len(cros_type)==6:
            data.columns=['Yi11', 'Yi12','Yi13','Yi21','Yi22', 'Yi23']
            
        elif len(cros_type)==9:
            data.columns=['Yi11', 'Yi12','Yi13','Yi21','Yi22', 'Yi23','Yi31','Yi32', 'Yi33']
            
        return data
    
    def Link(self,cros_type,params_value):
        '''
        用途:計算平均值
        1.cors_type:交叉設計
        2.params_value:若帶入的是真值則為真實平均值，若為MLE則為估計量
        
        '''
        
        #找出交叉設計中B藥與C藥的組別與時間點
        loc_B=np.argwhere(np.array(list(cros_type)) =='B')
        loc_C=np.argwhere(np.array(list(cros_type)) =='C')
        #只有B藥        
        if len(loc_C) == 0:
            #2X2 e.g. ABBA
            if len(cros_type)==4:
                tao,eta,gamma,delta=params_value
                params = np.array([tao,eta,gamma,delta])
                covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
                covariate[1,tuple(loc_B.T)] = 1
            #2X3 e.g. ABBBAA
            if len(cros_type)==6:
                tao,eta,gamma1,gamma2,delta=params_value
                params = np.array([tao,eta,gamma1,gamma2,delta])
                covariate=np.array([[1,1,1,1,1,1],[0,0,0,0,0,0],[0,1,0,0,1,0],[0,0,1,0,0,1],[0,0,0,1,1,1]])   
                covariate[1,tuple(loc_B.T)] = 1
            #3X3 e.g.AABABABAA
            if len(cros_type)==9:
                tao,eta,gamma1,gamma2,delta1,delta2=params_value
                params = np.array([ tao,eta,gamma1,gamma2,delta1,delta2])
                covariate=np.array([[1,1,1]*3,[0,0,0]*3,[0,1,0]*3,[0,0,1]*3,[0,0,0,1,1,1,0,0,0],[0,0,0,0,0,0,1,1,1]])
                covariate[1,tuple(loc_B.T)] = 1
        #3X3 有B藥與C藥    
        elif (len(loc_C) != 0) and (len(cros_type)==9):
            tao,eta1,eta2,gamma1,gamma2,delta1,delta2=params_value
            params = np.array([tao,eta1,eta2,gamma1,gamma2,delta1,delta2])
            covariate=np.array([[1,1,1]*3,[0,0,0]*3,[0,0,0]*3,[0,1,0]*3,[0,0,1]*3,[0,0,0,1,1,1,0,0,0],[0,0,0,0,0,0,1,1,1]])
            covariate[1,tuple(loc_B.T)] = 1
            covariate[2,tuple(loc_C.T)] = 1
        
        return params,covariate
    
    def Mean(self,params,covariate):
        return np.exp(np.dot(params,covariate))
    
    def MatrixI(self,cros_type,params,covariate,mean_,sample_size):
        num_seq= 3 if len(cros_type)==9 else 2
        pi=sample_size/(sample_size*num_seq)        
        matrixI=pi*np.array([[np.dot(mean_,covariate[i]*covariate[j]) for j in range(len(params))] for i in range(len(params))])
        return matrixI
    
    def MatrixV(self,cros_type,params,covariate,mean_,data,sample_size):
        num_seq= 3 if len(cros_type)==9 else 2
        score=[np.dot((data-mean_).to_numpy(),covariate[i].transpose()) for i in range(len(params))]
        matrixV=np.array([[sum(score[i]*score[j]) for j in range(len(params))] for i in range(len(params))])/(sample_size*num_seq)
        return matrixV
          
      
    def MLE(self,cros_type,data):
        '''
        用途：計算矩陣V(Poisson模型推導結果)
        變數解釋：
        1.cros_type:序列設計方式 e.g.'ABBA'
        2.data:生成的資料
        3.initial_guess
        initial_guess = np.array([0.841, 0.318, 0.11, 0.223, 0.202, 0.262])
        initial_guess = np.array([0.95, 0.372, 0.354, 0.236, 0.226, 0.206, 0.264])
        '''
        #3X3 pi=(seq_size)/(seq_size*3)
        #3X2、2X2 pi=(seq_size)/(seq_size*2)
        initial_guess = np.array(self.params_value)-0.05
        
        if cros_type=='ABBA':
            tao_hat=np.log(data['Yi11'].mean())
            eta_hat=0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi21'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi22'].sum()))
            gamma_hat=0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi22'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi21'].sum()))
            delta_hat=0.5*(np.log(data['Yi21'].sum())+np.log(data['Yi22'].sum())-np.log(data['Yi11'].sum())-np.log(data['Yi12'].sum()))
            
            # MLE
            estimate= pd.DataFrame({'tao_hat': tao_hat, 
                                    'eta_hat': eta_hat, 
                                    'gamma_hat':gamma_hat, 
                                    'delta_hat': delta_hat},index=[0])
        elif cros_type=='ABBBAA':       
            tao_hat = np.log(np.mean(data['Yi11']))
            eta_hat =0.25*(np.log(data['Yi12'].sum())+np.log(data['Yi13'].sum())+2*np.log(data['Yi21'].sum())-np.log(data['Yi22'].sum())-np.log(data['Yi23'].sum())-2*np.log(data['Yi11'].sum()))
            gamma1_hat = 0.5*(np.log(data['Yi12'].sum())+np.log(data['Yi22'].sum())-np.log(data['Yi21'].sum())-np.log(data['Yi11'].sum()))
            gamma2_hat = 0.5*(np.log(data['Yi13'].sum())+np.log(data['Yi23'].sum())-np.log(data['Yi21'].sum())-np.log(data['Yi11'].sum()))
            delta_hat = 0.5*(np.log(data['Yi21'].sum())-np.log(data['Yi11'].sum()))-0.25*(np.log(data['Yi12'].sum())+np.log(data['Yi13'].sum())-np.log(data['Yi22'].sum())-np.log(data['Yi23'].sum()))
            '''
            eta_hat =0.25*(np.log(np.mean(data['Yi12']))+np.log(np.mean(data['Yi13']))+2*np.log(np.mean(data['Yi21']))-np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi23']))-2*np.log(np.mean(data['Yi11'])))
            gamma1_hat = 0.5*(np.log(np.mean(data['Yi12']))+np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi21']))-np.log(np.mean(data['Yi11'])))
            gamma2_hat = 0.5*(np.log(np.mean(data['Yi13']))+np.log(np.mean(data['Yi23']))-np.log(np.mean(data['Yi21']))-np.log(np.mean(data['Yi11'])))
            delta_hat = 0.5*(np.log(np.mean(data['Yi21']))-np.log(np.mean(data['Yi11'])))-0.25*(np.log(np.mean(data['Yi12']))+np.log(np.mean(data['Yi13']))-np.log(np.mean(data['Yi22']))-np.log(np.mean(data['Yi23'])))
            '''
            # MLE
            estimate= pd.DataFrame({'tao': tao_hat, 
                                    'eta_hat': eta_hat, 
                                    'gamma1_hat':gamma1_hat, 
                                    'gamma2_hat':gamma2_hat,
                                    'delta_hat': delta_hat},index=[0])
        elif cros_type=='AABABABAA':
            def L_332(params):
                tao, eta, gamma1, gamma2, delta1, delta2=params
                f1,f2,f3=0,0,0
                for i in range(len(data)):
                    f1=f1+tao*data.at[i,'Yi11']+(tao+gamma1)*data.at[i,'Yi21']+(tao+eta+gamma2)*data.at[i,'Yi31']-np.exp(tao)-np.exp(tao+gamma1)-np.exp(tao+eta+gamma2)
                    f2=f2+(tao+delta1)*data.at[i,'Yi12']+(tao+eta+gamma1+delta1)*data.at[i,'Yi22']+(tao+gamma2+delta1)*data.at[i,'Yi32']-np.exp(tao+delta1)-np.exp(tao+eta+gamma1+delta1)-np.exp(tao+gamma2+delta1)
                    f3=f3+(tao+eta+delta2)*data.at[i,'Yi13']+(tao+gamma1+delta2)*data.at[i,'Yi23']+(tao+gamma2+delta2)*data.at[i,'Yi33']-np.exp(tao+eta+delta2)-np.exp(tao+gamma1+delta2)-np.exp(tao+gamma2+delta2)
                return f1+f2+f3
            solve = scipy.optimize.minimize (lambda params: -L_332(params), initial_guess)
            MLE_opt=solve.x
            estimate= pd.DataFrame({'tao_hat': MLE_opt[0], 
                                    'eta_hat': MLE_opt[1], 
                                    'gamma1_hat': MLE_opt[2], 
                                    'gamma2_hat': MLE_opt[3], 
                                    'delta1_hat': MLE_opt[4],
                                    'delta2_hat': MLE_opt[5],
                                    },index=[0])
        elif cros_type=='ABCBCACAB':
            def L_3331(params):
                tao, eta1, eta2, gamma1, gamma2, delta1, delta2=params
                f1,f2,f3=0,0,0
                for i in range(len(data)):
                    f1=f1+tao*data.at[i,'Yi11']+(tao+eta1+gamma1)*data.at[i,'Yi21']+(tao+eta2+gamma2)*data.at[i,'Yi31']-np.exp(tao)-np.exp(tao+eta1+gamma1)-np.exp(tao+eta2+gamma2)
                    f2=f2+(tao+eta1+delta1)*data.at[i,'Yi12']+(tao+eta2+gamma1+delta1)*data.at[i,'Yi22']+(tao+gamma2+delta1)*data.at[i,'Yi32']-np.exp(tao+eta1+delta1)-np.exp(tao+eta2+gamma1+delta1)-np.exp(tao+gamma2+delta1)
                    f3=f3+(tao+eta2+delta2)*data.at[i,'Yi13']+(tao+gamma1+delta2)*data.at[i,'Yi23']+(tao+eta1+gamma2+delta2)*data.at[i,'Yi33']-np.exp(tao+eta2+delta2)-np.exp(tao+gamma1+delta2)-np.exp(tao+eta1+gamma2+delta2)
                return f1+f2+f3
            solve = scipy.optimize.minimize (lambda params: -L_3331(params), initial_guess)
            MLE_opt=solve.x
            estimate= pd.DataFrame({'tao_hat': MLE_opt[0], 
                                    'eta1_hat': MLE_opt[1], 
                                    'eta2_hat': MLE_opt[2], 
                                    'gamma1_hat': MLE_opt[3], 
                                    'gamma2_hat': MLE_opt[4],
                                    'delta1_hat': MLE_opt[5],
                                    'delta2_hat': MLE_opt[6]},index=[0])
        elif cros_type=='BACACBBCA':        
            def L_3332(params):
                tao, eta1, eta2, gamma1, gamma2, delta1, delta2=params
                f1,f2,f3=0,0,0
                for i in range(len(data)):
                    f1=f1+(tao+eta1)*data.at[i,'Yi11']+(tao+gamma1)*data.at[i,'Yi21']+(tao+eta2+gamma2)*data.at[i,'Yi31']-np.exp(tao+eta1)-np.exp(tao+gamma1)-np.exp(tao+eta2+gamma2)
                    f2=f2+(tao+delta1)*data.at[i,'Yi12']+(tao+eta2+gamma1+delta1)*data.at[i,'Yi22']+(tao+eta1+gamma2+delta1)*data.at[i,'Yi32']-np.exp(tao+delta1)-np.exp(tao+eta2+gamma1+delta1)-np.exp(tao+eta1+gamma2+delta1)
                    f3=f3+(tao+eta1+delta2)*data.at[i,'Yi13']+(tao+eta2+gamma1+delta2)*data.at[i,'Yi23']+(tao+gamma2+delta2)*data.at[i,'Yi33']-np.exp(tao+eta1+delta2)-np.exp(tao+eta2+gamma1+delta2)-np.exp(tao+gamma2+delta2)
                return f1+f2+f3
            solve = scipy.optimize.minimize (lambda params: -L_3332(params), initial_guess)
            MLE_opt=solve.x
            estimate= pd.DataFrame({'tao_hat': MLE_opt[0], 
                                    'eta1_hat': MLE_opt[1], 
                                    'eta2_hat': MLE_opt[2], 
                                    'gamma1_hat': MLE_opt[3], 
                                    'gamma2_hat': MLE_opt[4],
                                    'delta1_hat': MLE_opt[5],
                                    'delta2_hat': MLE_opt[6]},index=[0])
        elif cros_type=='BBAACBCAC':
            def L_3333(params):
                tao, eta1, eta2, gamma1, gamma2, delta1, delta2=params
                f1,f2,f3=0,0,0
                for i in range(len(data)):
                    f1=f1+(tao+eta1)*data.at[i,'Yi11']+(tao+eta1+gamma1)*data.at[i,'Yi21']+(tao+gamma2)*data.at[i,'Yi31']-np.exp(tao+eta1)-np.exp(tao+eta1+gamma1)-np.exp(tao+gamma2)
                    f2=f2+(tao+delta1)*data.at[i,'Yi12']+(tao+eta2+gamma1+delta1)*data.at[i,'Yi22']+(tao+eta1+gamma2+delta1)*data.at[i,'Yi32']-np.exp(tao+delta1)-np.exp(tao+eta2+gamma1+delta1)-np.exp(tao+eta1+gamma2+delta1)
                    f3=f3+(tao+eta2+delta2)*data.at[i,'Yi13']+(tao+gamma1+delta2)*data.at[i,'Yi23']+(tao+eta2+gamma2+delta2)*data.at[i,'Yi33']-np.exp(tao+eta2+delta2)-np.exp(tao+gamma1+delta2)-np.exp(tao+eta2+gamma2+delta2)
                return f1+f2+f3
            solve = scipy.optimize.minimize (lambda params: -L_3333(params), initial_guess)
            MLE_opt=solve.x
            estimate= pd.DataFrame({'tao_hat': MLE_opt[0], 
                                    'eta1_hat': MLE_opt[1], 
                                    'eta2_hat': MLE_opt[2], 
                                    'gamma1_hat': MLE_opt[3], 
                                    'gamma2_hat': MLE_opt[4],
                                    'delta1_hat': MLE_opt[5],
                                    'delta2_hat': MLE_opt[6]},index=[0])
        return estimate

        
 

sim_time=2000
design_type=['ABBA','ABBBAA']#,'AABABABAA','ABCBCACAB','BACACBBCA','BBAACBCAC']   
pv_by_design=np.array([[1.0,0.7,0.3,0.2],[1.0,0.7,0.3,0.3,0.2],[1.0,0.7,0.3,0.3,0.2,0.2],[1.0,0.7,0.7,0.3,0.3,0.2,0.2],[1.0,0.7,0.7,0.3,0.3,0.2,0.2],[1.0,0.7,0.7,0.3,0.3,0.2,0.2]])#根據不同交叉設計產生的參數須給定真值
seq_size=np.array([25,50,100,200])
#next:simulation and export to excel

def Simulation(runtime,design_type,pv,seqsize):
    np.set_printoptions(suppress=True,precision=4)
    mle_ind,mle_cor=pd.DataFrame(),pd.DataFrame()
    I_ind, V_ind,I_cor, V_cor = 0, 0, 0, 0
    for i in range(runtime):    
        #independent
        CROS_ind=CrossoverDesign(cros_type=design_type,data_type='ind',cor_par=0.1,params_value=pv,sample_size=seqsize,eta0=0)
        mle_i=CROS_ind.estimate
        mle_ind = mle_ind.append(mle_i,ignore_index=True)
        I_i=CROS_ind.I_hat
        V_i=CROS_ind.V_hat
        I_ind+=I_i
        V_ind+=V_i
        CROS_cor=CrossoverDesign(cros_type=design_type,data_type='cor',cor_par=0.1,params_value=pv,sample_size=seqsize,eta0=0)
        mle_i=CROS_cor.estimate
        mle_cor = mle_cor.append(mle_i,ignore_index=True)
        I_i=CROS_cor.I_hat
        V_i=CROS_cor.V_hat
        I_cor+=I_i
        V_cor+=V_i
    

    I_ind=I_ind/sim_time
    V_ind=V_ind/sim_time
    I_cor=I_cor/sim_time
    V_cor=V_cor/sim_time
    
    return mle_ind.mean(),I_ind,V_ind,mle_ind.cov(ddof=1),mle_cor.mean(),I_cor,V_cor,mle_cor.cov(ddof=1)
#mle_ind.mean(),lin.inv(I_ind).dot(V_ind).dot(lin.inv(I_ind)),mle_ind.cov(ddof=1),mle_cor.mean(),lin.inv(I_cor).dot(V_cor).dot(lin.inv(I_cor)),mle_cor.cov(ddof=1)

np.random.seed(980716)
output_path="C:/Github/Study_CrossoverDesign/SimOutput"
today=str(datetime.date.today().month)+'0'+str(datetime.date.today().day)

for (ix,design) in enumerate(design_type):
    #writer = openpyxl.Workbook()
    
    writer=pd.ExcelWriter(output_path+'/{}_{}.xlsx'.format(design,today), engine="openpyxl")
    num_seq= 2 if len(design)!=9 else 3
    
    for seq in seq_size:        
        MLE_ind,I_ind_,V_ind_,cov_ind,MLE_cor,I_cor_,V_cor_,cov_cor=Simulation( runtime=sim_time, design_type=design, pv = pv_by_design[ix], seqsize=seq)
        sys.stdout.write('\rReading : '+design+str(seq))
        np.set_printoptions(suppress=True,precision=5)
        MLE_ind.to_excel(writer, sheet_name="MLE_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')        
        pd.DataFrame(I_ind_).to_excel(writer, sheet_name="I_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame(V_ind_).to_excel(writer, sheet_name="V_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame(lin.inv(I_ind_)).to_excel(writer, sheet_name="inv(I)_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        #pd.DataFrame(lin.inv(I_ind_).dot(V_ind_).dot(lin.inv(I_ind_))).to_excel(writer, sheet_name="invI_V_inv_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame((num_seq*seq*cov_ind)).to_excel(writer, sheet_name="NS_ind"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        MLE_cor.to_excel(writer, sheet_name="MLE_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')        
        pd.DataFrame(I_cor_).to_excel(writer, sheet_name="I_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        #pd.DataFrame(V_cor_).to_excel(writer, sheet_name="V_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        #pd.DataFrame(lin.inv(I_cor_)).to_excel(writer, sheet_name="inv(I)_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame(lin.inv(I_cor_).dot(V_cor_).dot(lin.inv(I_cor_))).to_excel(writer, sheet_name="invI_V_inv_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame((num_seq*seq*cov_cor)).to_excel(writer, sheet_name="NS_cor"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        writer.save()
    
    writer.close()
    

 
        '''MLE_ind,mb_ind,cov_ind,MLE_cor,mb_cor,cov_cor=Simulation(runtime=2000, design_type=design, pv=pv_by_design[ix], seqsize=seq)
        sys.stdout.write('\rReading : '+design+str(seq))
        np.set_printoptions(suppress=True,precision=5)
        MLE_ind.to_excel(writer, sheet_name="MLE"+str(seq), engine='openpyxl', encoding='utf_8_sig')        
        pd.DataFrame(mb_ind).to_excel(writer, sheet_name="invI_V_inv"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame(cov_ind).to_excel(writer, sheet_name="S"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame((num_seq*seq*cov_ind)).to_excel(writer, sheet_name="NS"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        MLE_cor.to_excel(writer, sheet_name="MLE"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame(mb_cor).to_excel(writer, sheet_name="invI_V_inv"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame(cov_cor).to_excel(writer, sheet_name="S"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        pd.DataFrame((num_seq*seq*cov_cor)).to_excel(writer, sheet_name="NS"+str(seq), engine='openpyxl', encoding='utf_8_sig')
        writer.save()
    writer.close()'''

