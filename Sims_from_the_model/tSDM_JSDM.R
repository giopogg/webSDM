###### Script to test JSDM vs tSDM

rm(list=ls())

library(dplyr)
library(BayesComm)
set.seed(17121993)

####################################
## Gaussian data


# Simulate data

# 100 sites, one environmental variable x, two species y1 and y2.
# y1 only depends linearly on x, while y2 depends linearly on both y1 and x.

n=1000
x= runif(n)
# y1 has intercept=1 and "niche parameter"=3
y1=rnorm(n=n,mean=1+3*x,1)
# y2 has intercept=1, "niche parameter"=2 and "biotic parameter"=1
y2=rnorm(n=n,mean=1+2*x+y1,1)


#SDM
SDM1=lm(y1~x)
SDM2=lm(y2~x)
summary(SDM2)
# We infer a wrong "niche" parameter=5 for y2, instead of the true one that is equal to 2.

# gaussian JSDM (https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Multivariate-Linear-Models.pdf)
JSDM= lm(cbind(y1,y2)~x)
summary(JSDM)

# As we know, at least for the gaussian case, the ML estimators of the niche parameter coincide with
# indipendent linear regressions (SDM). Therefore, just as above, we find the wrong niche parameter (=5) for species 2


## trophic SDM
tSDM1=lm(y1~x)
summary(tSDM1)
tSDM2=lm(y2~x+y1)
summary(tSDM2)
# We correctly infer the right "niche" parameter=2 (well, easy because the simulated data are sampled from this model)


###### With presence-absence

# We keep the same framework of above but we simulate presence-absence data using the probit model
n=2000
x= runif(n)

# y1 has intercept=1 and "niche parameter"=-1
y1_lat=rnorm(n=n,mean=1-x,1)
y1=ifelse(y1_lat>0,1,0)
# y2 has intercept=2, "niche parameter"=-2 and "biotic parameter"=-2
y2_lat=rnorm(n=n,mean=2-2*x-2*y1,1)
y2=ifelse(y2_lat>0,1,0)

#JSDM (BayesComm)
JSDM=BC(Y=cbind(y1,y2),X=as.matrix(x),model="full",its=1000)

#compute the posterior mean of niche coefficients from the MCMC chains
B_jsdm =cbind(JSDM$trace$B[[1]],JSDM$trace$B[[2]])
meanB_s = apply(B_jsdm,2,mean)
meanB_s=matrix(meanB_s,ncol=2)
colnames(meanB_s)=c("y1","y2")
rownames(meanB_s)=c("Intercept","x")
meanB_s

#again, we do not infer the right niche parameters

## trophic SDM (just two SDMs but including y1 to predict y2)
tSDM1=glm(y1~x,family=binomial(link = "probit"))
summary(tSDM1)
tSDM2=glm(y2~x+y1,family=binomial(link = "probit"))
summary(tSDM2)
#again we infer the right niche parameters (again, easy because the simulated data are sampled from this model)

