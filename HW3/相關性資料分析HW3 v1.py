# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 21:30:09 2022

@author: 懿萱
"""

import numpy as np
import pandas as pd
import scipy
import statsmodels.api as sm
import statsmodels.formula.api as smf

path="C:\\Users\\a7086\\PycharmProjects\\pythonProject\\2022HW-3cross.xls"  #'C://Users//cherl//GitHub//NCU//Study//HW3//2022HW-3cross.xls'
data = pd.read_excel(path)
data[['treatment','time']]=data[['treatment','time']].replace([ 1, 2, 3], [ 0, 1, 2])
data['treatment']=data['treatment'].astype('category')
'''Model1 Y~alpha+beta1*treatment'''
#Logit OLS
x = sm.add_constant(data['treatment'])
logit1 = sm.Logit(data['Y'], x)
res_logit1 = logit1.fit(method='newton')
print(res_logit1.summary())


#GLM
bin_model1 = sm.GLM(data['Y'], x, family=sm.families.Binomial())
glm_results1 = bin_model1.fit()
print(glm_results1.summary())
from statsmodels.genmod.generalized_estimating_equations import GEE
#GEE-independet
b = sm.families.Binomial()
ind = sm.cov_struct.Independence()
gee_ind1 =smf.gee("Y~treatment", "ID", data, cov_struct=ind, family=b)
#GEE(data['Y'], x,data['ID'],None,b,ind)

res_gee_ind1 = gee_ind1.fit()
print(res_gee_ind1.summary())
res_gee_ind1.bse.values
res_gee_ind1.cov_params_default
res_gee_ind1.cov_naive
res_gee_ind1.cov_robust

#GEE-Exchangeable
exg = sm.cov_struct.Exchangeable()
gee_exg1 = smf.gee("Y~treatment", "ID", data, cov_struct=exg, family=b)
res_gee_exg1 = gee_exg1.fit()
print(res_gee_exg1.summary())

res_gee_exg1.bse.values
res_gee_exg1.cov_params_default
res_gee_exg1.cov_naive
res_gee_exg1.cov_robust
res_gee_exg1.resid

'''Model2 Y~alpha+beta1*treatment+beta2*time'''
#Logit OLS
x = sm.add_constant(data[['treatment','time']])
logit2 = sm.Logit(data['Y'], x)
res_logit2 = logit2.fit(method='newton')
print(res_logit2.summary())


#GLM
bin_model2 = sm.GLM(data['Y'], x, family=sm.families.Binomial())
glm_results2 = bin_model2.fit()
print(glm_results2.summary())

#GEE-independet
b = sm.families.Binomial()
ind = sm.cov_struct.Independence()
gee_ind2 = smf.gee("Y~treatment+time", "ID", data, cov_struct=ind, family=b)
res_gee_ind2 = gee_ind2.fit()
print(res_gee_ind2.summary())
#GEE-Exchangeable
exg = sm.cov_struct.Exchangeable()
gee_exg2 = smf.gee("Y~treatment+time", "ID", data, cov_struct=exg, family=b)
res_gee_exg2 = gee_exg2.fit()
print(res_gee_exg2.summary())

'''Model1 v.s. Model2'''
print('Model 1 v.s. Model 2')
print('LR_statistic(ind)',-2*(gee_ind1.fit().llf-gee_ind2.fit().llf))
print('p-value(ind)',scipy.stats.chi2.sf(-2*(gee_ind1.fit().llf-gee_ind2.fit().llf), 2))
print('LR_statistic(exchangeable)',-2*(gee_exg1.fit().llf-gee_exg2.fit().llf))
print('p-value(exchangeable)',scipy.stats.chi2.sf(-2*(gee_exg1.fit().llf-gee_exg2.fit().llf), 2))


'''Model3 Y~alpha+beta1*treatment+beta2*time+beta3*time*treatment'''
#Logit OLS
data['time:treatment']=data['treatment']*data['time']
x = sm.add_constant(data[['treatment','time','time:treatment']])
logit3 = sm.Logit(data['Y'], x)
res_logit3 = logit3.fit(method='newton')
print(res_logit3.summary())


#GLM
bin_model3 = sm.GLM(data['Y'], x, family=sm.families.Binomial())
glm_results3 = bin_model3.fit()
print(glm_results3.summary())

#GEE-independet
b = sm.families.Binomial()
ind = sm.cov_struct.Independence()
gee_ind3 = smf.gee("Y~treatment+time+time*treatment", "ID", data, cov_struct=ind, family=b)
res_gee_ind3 = gee_ind3.fit()
print(res_gee_ind3.summary())
#GEE-Exchangeable
exg = sm.cov_struct.Exchangeable()
gee_exg3 = smf.gee("Y~treatment+time+time*treatment", "ID", data, cov_struct=exg, family=b)
res_gee_exg3 = gee_exg3.fit()
print(res_gee_exg3.summary())


'''Model2 v.s. Model3'''
print('Model 2 v.s. Model 3')
print('LR_statistic(ind)',2*(gee_ind3.fit().llf-gee_ind2.fit().llf))
print('p-value(ind)',scipy.stats.chi2.sf(2*(gee_ind3.fit().llf-gee_ind2.fit().llf), 2))
print('LR_statistic(exchangeable)',2*(gee_exg3.fit().llf-gee_exg2.fit().llf))
print('p-value(exchangeable)',scipy.stats.chi2.sf(2*(gee_exg3.fit().llf-gee_exg2.fit().llf), 2))



