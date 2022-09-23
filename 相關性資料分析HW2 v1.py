# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 21:30:09 2022

@author: 懿萱
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

data = pd.read_excel('C://Users//cherl//GitHub//NCU//data.xlsx')

g = sm.families.Gaussian().identity()

ind = sm.cov_struct.Independence()

mod_g = smf.gee("bmi~index1+age+aog+sex", "id", data,cov_struct=ind, family=g)

res_g = mod_g.fit()

print(res_g.summary())

ga = sm.families.Gamma()

ind = sm.cov_struct.Independence()

mod_ga = smf.gee("bmi~index1+age+aog+sex", "id", data,cov_struct=ind, family=ga)

res_ga = mod_ga.fit()

print(res_ga.summary())

p = sm.families.Gaussian()

ind = sm.cov_struct.Independence()

mod_p = smf.gee("bmi~index1+age+aog+sex", "id", data,cov_struct=ind, family=p)

res_p = mod_p.fit()

print(res_p.summary())