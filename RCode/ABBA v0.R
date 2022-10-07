set.seed(716)
setwd("C:/Github/Study_CrossoverDesign/RCode")
sim_time,seq_size,cor_param,cros_type=1000,100,0.1,'ABBA'
params = c(1.0, 0.67, 0.23, 0.12)#tao,eta,gamma,delta
covariate=np.array([[1,1,1,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]])
mean_true=np.exp(np.dot(params.transpose(),covariate))
pi=(seq_size)/(seq_size*2)
I=pi*np.array([[np.dot(mean_true,covariate[i]*covariate[j]) for j in range(len(cros_type))] for i in range(len(cros_type))])
df=pd.DataFrame(np.array([poisson.rvs(p, size=seq_size) for p in mean_true]).T.tolist(),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
df=pd.DataFrame(np.random.poisson(lam=mean_true, size=(seq_size, len(cros_type))),columns = ['Yi11', 'Yi12','Yi21', 'Yi22'])
df.mean()
df.corr()
data<-data.frame(rpois(100, lambda = 3))