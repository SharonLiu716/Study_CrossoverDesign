soliq<-read.table('C:/Users/cherl/GitHub/NCU/data.txt', col.names=c("id", "aog","age", "timeill", "sex","index1","bmi"))
write.csv(soliq, 'C:/Users/cherl/GitHub/NCU/data.csv')
soliq$id <- factor(soliq$id)
soliq$sex <- factor(soliq$sex)


glm1<-glm(bmi~index1+age+aog+factor(sex)+timeill,data=soliq,family="gaussian")
summary(glm1)

glm2<-glm(bmi~index1+as.factor(age)+aog+factor(sex)+factor(timeill),data=soliq,family="Gamma")
summary(glm2)

glm3<-glm(bmi~index1+age+aog+factor(sex)+factor(timeill),data=soliq,family="poisson")
summary(glm3)

library(geepack)
fit0<-geeglm(bmi~factor(index1)+factor(age)+factor(aog)+factor(sex)+factor(timeill),data=soliq,family=gaussian,id=id,corstr = 'exchangeable')
summary(fit0)
'''
proc sort data=skin; by id year;
run;
proc genmod data=skin;
class id yearcat;
model y=year trt*year / dist=poisson link=log type3 wald waldci;
repeated subject=id / withinsubject=yearcat type=un;
run;
-----------------------------------------------------------------------------------------------
Analysis Of GEE Parameter Estimates
Empirical Standard Error Estimates
Standard 95% Confidence
Parameter Estimate Error       Limits       Z    Pr > |Z|
Intercept -1.3341  0.0815 -1.4938 -1.1743 -16.37 <.0001
year      -0.0090  0.0271 -0.0622  0.0441 -0.33  0.7392
year*trt   0.0429  0.0319 -0.0195  0.1053  1.35  0.1781'''
