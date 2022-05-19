##########################################################################################################
############ TEST tSDM
##########################################################################################################
rm(list=ls())
library(igraph)
library(GGally)
library(intergraph)
library(gridExtra)
library(dismo)
library(igraph)
library(Matrix)
library(cheddar)
library(cowplot)
library(GGally)
library(intergraph)
library(gridExtra)
library(ggpubr)
library(dismo)
library(coda)
library(transport)
library(dplyr)
library(rstan)
library(rstanarm)
library(arm)
library(glmnet)
library(formula.tools)
library(abind)
library(brms)
library(stringr)
library(bayesplot)
library(parallel)

#library(pcalg) not working with R 4.0

wd="~/Documents/GitHub/trophicSDM/Sims_from_the_model/"
setwd(wd)
source('../tSDM_functions.R')


n=1000 #number of sites
n_oos= 500#number of sites out of sample
S=6 #number of species
#This is to chose if we want to include a quadratic term for X_2 or not
linear=T
# if yes, chose whether we want a raw polynomial for X_2 (poly=T) or an orthogonal one (poly=F)
if(!linear) poly=F 

# This is to choose how we create data:
# "true logit" : the correct way of drawing from a logit regression. The linear term (V), that has no noise, is transformed with the inverse logit function to obtain the probabilities (prob), then Y is sample randomly from a bernoulli with probability prob : https://stats.stackexchange.com/questions/46523/how-to-simulate-artificial-data-for-logistic-regression/46525
# the problem of this approach is that it is hard (when the regression coefficients are randomly chose, as here, to obtain probabilities of presence for presences that are well separated from absences. Indeed we have very low TSS and AUC
simul_type="true_logit"
# "logit_threshold": simulates the linear term (V) with a gaussian noise, then transformes it with an inverse logit function to obtain the probabilities (prob), and Y is 1 if prob>0.5, 0 otherwise. This is not the correct way to simulate, but we have well separated probabilities of presences/absences (indeed AUC and TSS are 1)
#simul_type="logit_threshold"
# "logit_prob_covariate" : like logit, but instead of using binary presence/absences as predictors we use the probabilities of presence. BUT this does not work, since when we fit the model we use presence/absences. Therefore, the model does not even retrieve the parameters
#simul_type="logit_prob_covariate"

#number of covariates (including polynomial terms)
if(linear){ K=2 }else{ K=3 }#number of covariates (including polynomial terms) }

### Build the newtork !

A=matrix(0,nrow=S,ncol=S)
#or sample random DAG with pcalg (not working for me)
while(!is.connected(graph_from_adjacency_matrix(A))){
  A=A_w=matrix(0,nrow=S,ncol=S)
  for(j in 1:(S-1)){A[j,c((j+1):S)]=sample(x=c(0,1),size=S-j,replace=T)}
}

# transpose A since we want high numbers to be top predators (just for a matter of visualization)
A=t(A) 

# Do a loop to create simulated data:
# we create data until we have that each species has a prevalence between 0.1 and 0.9 (if this is the case, cond=T, last line of the while loop)
cond=F
k=0 #just a counter of loops
while(!cond){
k=k+1 #update counter

#create arrays for data
Y=V=prob=matrix(0,nrow=n+n_oos,ncol=S)

#create weights of the network (i.e. the regression coefficients)
nonzeros=which(A!=0)
for(i in nonzeros) A_w[i]=runif(1,min=-2,max=2) #min and max can be tuned (or even sampled)

# sample environmental covariates
X_1= runif(n+n_oos,max=1)
X_2= runif(n+n_oos,max=1)

# build the matrix to create the data, depending on whether we want to include polynomial terms (and how)
if(linear){
  X=cbind(1,X_1,X_2)
}else{
  if(!poly){
    X=cbind(1,X_1,X_2,X_2^2) 
  }else{
    X=cbind(1,X_1,poly(X_2,degree=2))
  }
}

#create regression coefficients
B=matrix(runif(S*(K+1),min=-2,max=2),ncol=K+1) #+1 bc of the intercept

if(simul_type=="logit_threshold"){
  #because otherwise it's not random
 sigma_2 = runif(S,min = 0.5,max=1)
}

#sample data for each species
for(j in 1:S){
  
if(simul_type=="true_logit"){
  V[,j]=X%*%B[j,]+Y%*%A_w[j,]
  prob[,j]=1/(1+exp(-V[,j]))
  Y[,j]=rbinom(n=n+n_oos,size=1,prob=prob[,j])
}
if(simul_type=="logit_threshold"){
    V[,j]=X%*%B[j,]+Y%*%A_w[j,]+rnorm(n+n_oos,0,sqrt(sigma_2[j]))
    prob[,j]=1/(1+exp(-V[,j]))
    for(i in 1:nrow(V)){ Y[i,j] = ifelse(prob[i,j]>0.5,1,0)}
}
  
if(simul_type=="logit_prob_covariate"){
    V[,j]=X%*%B[j,]+prob%*%A_w[j,]
    prob[,j]=1/(1+exp(-V[,j]))
    Y[,j]=rbinom(n=n+n_oos,size=1,prob=prob[,j])
}

}

# compute AUC and TSS for simulated data. If very low, we can't hope to do any better!
AUC_simul=TSS_simul=vector(length=ncol(Y))

for(j in 1:ncol(Y)){

  e_joint = evaluate(p=prob[which(Y[,j]==1),j], a=prob[which(Y[,j]==0),j])
  AUC_simul[j]=e_joint@auc
  TSS_simul[j]=max(e_joint@TPR+e_joint@TNR-1)
  
}


# conditions, we can add as much as we want. For now we keep only the condition that each species prevalence has to be between 0.1 and 0.9
# the commented condition is to guarantee high AUC for each species (but takes a lot of time to run, actually maybe impossible with AUC>0.8, should be lowered)
#if(length(which(AUC_simul>0.8))==S){
if(length(which(colSums(Y)/(n+n_oos) <0.9 & colSums(Y)/(n+n_oos)>0.1)) ==S) { cond = T}
}
#}

summary(Y)
summary(prob)

# Important: Y names should not be just raw numbers
colnames(Y)=colnames(A)=rownames(A)=paste0("Y_",as.character(1:S))

#build graph from A
G=graph_from_adjacency_matrix(A)

# Split dataset in train-test
prob_test=as.data.frame(prob[(n+1):(n+n_oos),])
prob=as.data.frame(prob[1:n,])

X_test=data.frame(X_1=X_1[(n+1):(n+n_oos)],X_2=X_2[(n+1):(n+n_oos)])
X_raw=data.frame(X_1=X_1[1:n],X_2=X_2[1:n]) 

Y_test=as.data.frame(Y[(n+1):(n+n_oos),])
Y=as.data.frame(Y[1:n,])



# Create formula
if(linear){
  env.form=as.list(rep("~X_1+X_2",S))
}else{
if(!poly){
env.form=as.list(rep("~X_1+X_2+I(X_2^2)",S))
}else{
  env.form=as.list(rep("~X_1+poly(X_2,2)",S))
}
}

env.form=lapply(env.form,FUN=as.formula)
names(env.form)=colnames(Y)

#notice that if one doesn't want to use I(X_2^2) but an ortogonal poly instead, he can just 
#create the X matrix in such a way, and use formula = NULL (or any kind of formula including
# the poly term as an additional variable)

# fit the models
m_glm = trophicSDM(Y=Y,X=X_raw,G=G,env.formula=env.form,penal=NULL,method="glm",family=binomial(link = "logit"),iter=1000)
#m_glm_net = trophicSDM(Y=Y,X=X_raw,G=G,env.formula=env.form,penal="elasticnet",method="glm",family=binomial(link = "logit"),iter=2000)
m_stan = trophicSDM(Y=Y,X=X_raw,G=G,env.formula=env.form,penal=NULL,method="stan_glm",family=binomial(link = "logit"),iter=1000,run.parallel=F)
#m_stan_hs =trophicSDM(Y=Y,X=X_raw,G=G,env.formula=env.form,penal="horshoe",method="stan_glm",family=binomial(link = "logit"),iter=2000)
m_bayes = trophicSDM(Y=Y,X=X_raw,G=G,env.formula=env.form,penal=NULL,method="bayesglm",family=binomial(link = "logit"),iter=1000)


# show traceplots for m_stan (also to check if there are some identification issues)
for(j in 1:S){
  traceplot(mcmc(as.matrix(m_stan$model[[paste0("Y_",j)]])))
}

###############
# Analyse results: obtain posterior means and 95% credible regions for each parameter of each species, and compare against the true ones

estimates_stan=estimates_glm= data.frame(p.est=double(),
                              est.97=double(),
                              est.02=double(),
                              true=double(),
                              type=factor(levels = c("biotic","abiotic")),
                              sp.id=factor())

for(j in 1:S){
  estimates_stan=rbind(estimates_stan,data.frame(p.est=coef(m_stan$model[[paste0("Y_",j)]]),
                                                 est.02=posterior_interval(m_stan$model[[paste0("Y_",j)]],prob = 0.95)[,1],
                                                 est.97=posterior_interval(m_stan$model[[paste0("Y_",j)]],prob = 0.95)[,2],
                                                 true=c(B[j,],A_w[j,which(A_w[j,]!=0)]),
                                                 type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(which(A_w[j,]!=0))))),
                                                 sp.id=as.factor(rep(j,length(coef(m_stan$model[[paste0("Y_",j)]]))))
                                    ))
  
  estimates_glm=rbind(estimates_glm,data.frame(p.est=coef(m_glm$model[[paste0("Y_",j)]]),
                                               est.02=confint(m_glm$model[[paste0("Y_",j)]],prob = 0.95)[,1],
                                               est.97=confint(m_glm$model[[paste0("Y_",j)]],prob = 0.95)[,2],
                                               true=c(B[j,],A_w[j,which(A_w[j,]!=0)]),
                                               type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(which(A_w[j,]!=0))))),
                                               sp.id=as.factor(rep(j,length(coef(m_glm$model[[paste0("Y_",j)]]))))
                                    ))
                                            
}



# plot pred vs true for stan_glm
p1=ggplot(data=estimates_stan,aes(x=true,y=p.est,colour=type,shape=sp.id))+ geom_point()+
  scale_shape_manual(values = c(1:2,4:(S+1))) +
  geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type))+
  geom_abline(intercept = 0, slope = 1, colour = "#66B2FF" , linetype= "dashed")+
  #coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+ 
  ggtitle("stan_glm()")+
xlab("True") + ylab("Inferred")+theme_minimal()

# plot pred vs true for glm
p2=ggplot(data=estimates_glm,aes(x=true,y=p.est,colour=type,shape=sp.id))+ geom_point()+
  scale_shape_manual(values = c(1:2,4:(S+1))) +
  geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type))+
  geom_abline(intercept = 0, slope = 1, colour = "#66B2FF" , linetype= "dashed")+
  #coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+ 
  ggtitle("glm()")+
  xlab("True") + ylab("Inferred")+theme_minimal()

#plot together with the trophic graph
#pdf(paste0(wd,"PLOTS/Simulations/Binary_simul_firsttry.pdf"))
grid.arrange(arrangeGrob(p1,p2,top="True vs inferred regression coefficients"),
             ggnet2(intergraph::asNetwork(G),arrow.gap = 0.05,arrow.size = 10,label=TRUE)+ggtitle("Known graph"),
             ncol=2)
dev.off()


####################################################################################################################
# Predictions
####################################################################################################################

##### In-sample, no err_prop
Xnew=X_raw #Xnew are the environmental variables that we used to fit the model (actually if Xnew=NULL, this is the default choice) 
error_prop_sample=1 #error_prop_sample is an important parameter. The greater it is, the largest the error propagation (and therefore the predictive credible regions). see tSDM_functions for more details

# This was just to check that it worked
#pred_stan=trophicSDM_predict(m=m_stan,Xnew=Xnew,binary.resp = T,prob.cov=F,pred_samples=1000,error_prop_sample=error_prop_sample)

# Just to check if we could get error out of some combinations of binary.resp and prob.cov (we don't!)
#pred_stan=trophicSDM_predict(m=m_stan,Xnew=Xnew,binary.resp = F,prob.cov=F,pred_samples=100,error_prop_sample=10)
#pred_stan=trophicSDM_predict(m=m_stan,Xnew=Xnew,binary.resp = T,prob.cov=T,pred_samples=100,error_prop_sample=10)
#pred_stan=trophicSDM_predict(m=m_stan,Xnew=Xnew,binary.resp = F,prob.cov=T,pred_samples=100,error_prop_sample=10)

#pred_bayes=trophicSDM_predict(m=m_bayes,Xnew,binary.resp = T,prob.cov=F,pred_samples=1000,error_prop_sample=error_prop_sample)
#pred_glm=trophicSDM_predict(m=m_glm,Xnew,binary.resp = T,prob.cov=F,pred_samples=1,error_prop_sample=error_prop_sample)

###########################################################################################################
# Check probabilities of presence pred vs known (because we simulated :) )

#run predictions
pred_stan_bin=trophicSDM_predict(m=m_stan,Xnew=Xnew,binary.resp = T,prob.cov=F,pred_samples=100,mode="out",run.parallel=T,verbose=T)

colnames(prob)=colnames(Y)
#reorder the columns of prob to be sure that we are comparing the same things
prob=prob[,names(topo_sort(G,mode="in"))]
Y=Y[,names(topo_sort(G,mode="in"))]

p.mean.stan.temp=lapply(pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,mean, MARGIN = 1))
p.mean.stan_bin=do.call(cbind, p.mean.stan.temp)
colnames(p.mean.stan_bin)= names(topo_sort(G,mode="in"))

p.q975.stan.temp=lapply(pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.975))
p.q975.stan_bin=do.call(cbind, p.q975.stan.temp)
colnames(p.q975.stan_bin)= names(topo_sort(G,mode="in"))

p.q025.stan.temp=lapply(pred_stan_bin$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.025))
p.q025.stan_bin=do.call(cbind, p.q025.stan.temp)
colnames(p.q025.stan_bin)= names(topo_sort(G,mode="in"))



table_no_err_pr_prob_bin=data.frame(obs=as.vector(as.matrix(prob)),pred=as.vector(as.matrix(p.mean.stan_bin)),
                 est.02=as.vector(as.matrix(p.q025.stan_bin)),est.97=as.vector(as.matrix(p.q975.stan_bin)),sp.name=rep(colnames(prob),each=nrow(prob)))

#plot predicted vs true probabilities of presence 
ggplot(data=table_no_err_pr_prob_bin,mapping=aes(x=obs,y=pred,col=sp.name))+geom_point(alpha=0.3)+
  #geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=sp.name),alpha=0.3)+
  geom_abline(slope=1,intercept = 0)+ggtitle(paste0("In-sample predictions stan_glm() \n error_prop_sample = ",error_prop_sample))+
  xlab("Observed") + ylab("predicted")+theme_minimal()+coord_fixed()

# notice that when simul_type="true_logit" we get very very weird results. This is because in the true data there is not neat separation of the probabilities between presencences and absences (indeed very low AUC and TSS), therefore,
# we have these weird "stripes" in the plots that correspond to different values of presences/absences of the preys. For example, if Y_2 eats Y_1, we will have two stripes, one where Y_1=1 and one where Y_1=0.
# This because prob_hat[,2] ~ B%*% X + Y_hat%*%A_w (or prob%*%A_w if prob.cov=T).
# To make this work, we need to have prob and Y that are well separated (with high AUC and TSS)


# Compute R2 and proportion of confidence interval that match
R2_no_err_pr_bin=sapply(names(topo_sort(G,mode="in")),function(x){cor(p.mean.stan_bin[,x],prob[,x])^2})
hist(R2_no_err_pr_bin)


CR_match_no_err_pr_bin=sapply(names(topo_sort(G,mode="in")), function(x) mean(sapply(1:n,FUN=function(y) 
  ifelse(prob[y,x]>p.q025.stan_bin[y,x] & prob[y,x]<p.q975.stan_bin[y,x],1,0))))
hist(CR_match_no_err_pr_bin)

##### JUST DEBUGS
# aa=prob[,2]
# bb=p.mean.stan_bin[,2]
# plot(p.mean.stan_bin[which(aa<bb),2],prob[which(aa<bb),2])
# plot(p.mean.stan_bin[which(aa>bb),2],prob[which(aa>bb),2])
# plot(p.mean.stan_bin[which(aa<bb),1],prob[which(aa<bb),1])
# plot(p.mean.stan_bin[which(aa>bb),1],prob[which(aa>bb),1])
# 
# summary(Y[which(aa<bb),1])
# summary(Y[which(aa>bb),1])
# 
# summary(p.mean.stan_bin[which(aa<bb),1])
# summary(p.mean.stan_bin[which(aa>bb),1])


###################################################################
# What happens when we use probabilities of presence for predictions? -> set prob.cov = T

pred_stan_prob_cov=trophicSDM_predict(m=m_stan,Xnew=Xnew,binary.resp = T,prob.cov=T,pred_samples=100,error_prop_sample=error_prop_sample)

p.mean.stan.temp=lapply(pred_stan_prob_cov$sp.prediction,FUN=function(x) apply(x$predictions.prob,mean, MARGIN = 1))
p.mean.stan_prob_cov=do.call(cbind, p.mean.stan.temp)
colnames(p.mean.stan_prob_cov)= names(topo_sort(G,mode="in"))

p.q975.stan.temp=lapply(pred_stan_prob_cov$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.975))
p.q975.stan_prob_cov=do.call(cbind, p.q975.stan.temp)
colnames(p.q975.stan_prob_cov)= names(topo_sort(G,mode="in"))

p.q025.stan.temp=lapply(pred_stan_prob_cov$sp.prediction,FUN=function(x) apply(x$predictions.prob,quantile, MARGIN = 1,probs=0.025))
p.q025.stan_prob_cov=do.call(cbind, p.q025.stan.temp)
colnames(p.q025.stan_prob_cov)= names(topo_sort(G,mode="in"))

table_no_err_pr_prob_cov=data.frame(obs=as.vector(as.matrix(prob)),pred=as.vector(as.matrix(p.mean.stan_prob_cov)),
                 est.02=as.vector(as.matrix(p.q025.stan_prob_cov)),est.97=as.vector(as.matrix(p.q975.stan_prob_cov)),sp.name=rep(colnames(prob),each=nrow(prob)))

#plot predicted vs observed
ggplot(data=table_no_err_pr_prob_cov,mapping=aes(x=obs,y=pred,col=sp.name))+geom_point(alpha=0.3)+
  geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=sp.name),alpha=0.3)+
  geom_abline(slope=1,intercept = 0)+ggtitle(paste0("In-sample predictions stan_glm() \n error_prop_sample = ",error_prop_sample))+
  xlab("Observed") + ylab("predicted")+theme_minimal()+coord_fixed()


# Compute R2 and proportion of confidence interval that match
R2_no_err_pr_prob_cov=sapply(names(topo_sort(G,mode="in")),function(x){cor(p.mean.stan_prob_cov[,x],prob[,x])^2})
hist(R2_no_err_pr_prob_cov)


CR_match_no_err_pr_prob_cov=sapply(names(topo_sort(G,mode="in")), function(x) mean(sapply(1:n,FUN=function(y) 
  ifelse(prob[y,x]>p.q025.stan_prob_cov[y,x] & prob[y,x]<p.q975.stan_prob_cov[y,x],1,0))))
hist(CR_match_no_err_pr_prob_cov)


####### Compare when prob_cov and when not

table_vs_prob_cov_no_err_pr = data.frame("bin_pred"=table_no_err_pr_prob_bin$pred,"prob.cov_pred"=table_no_err_pr_prob_cov$pred,"sp.name"=table_no_err_prob_bin$sp.name)

ggplot(data=table_vs_prob_cov_no_err_pr,aes(x=bin_pred,y=prob.cov_pred,col=sp.name))+geom_point()

cor(table_vs_prob_cov_no_err_pr$bin_pred,table_vs_prob_cov_no_err_pr$prob.cov_pred)+ggtitle("predictive means")
#0.997

#very similar correlation for predictive means


#but what about quantiles?

table_vs_prob_cov_no_err_pr_q.025 = data.frame("bin_pred"=table_no_err_pr_prob_bin$est.02,"prob.cov_pred"=table_no_err_pr_prob_cov$est.02,"sp.name"=table_no_err_prob_bin$sp.name)

ggplot(data=table_vs_prob_cov_no_err_pr_q.025,aes(x=bin_pred,y=prob.cov_pred,col=sp.name))+geom_point()+ggtitle("predictive 0.025 quantile")

cor(table_vs_prob_cov_no_err_pr$bin_pred,table_vs_prob_cov_no_err_pr$prob.cov_pred)
#0.997

table_vs_prob_cov_no_err_pr_q.975 = data.frame("bin_pred"=table_no_err_pr_prob_bin$est.97,"prob.cov_pred"=table_no_err_pr_prob_cov$est.97,"sp.name"=table_no_err_prob_bin$sp.name)

ggplot(data=table_vs_prob_cov_no_err_pr_q.975,aes(x=bin_pred,y=prob.cov_pred,col=sp.name))+geom_point()+ggtitle("predictive 0.975 quantile")

cor(table_vs_prob_cov_no_err_pr$bin_pred,table_vs_prob_cov_no_err_pr$prob.cov_pred)

#Strongly correlated, but when we use binary data we have larger credible intervals 

###########################################################################################################
# Check if we predict 0/1 well (AUC, TSS)

#We compute AUC and TSS in three cases: when using binary predictions as newdata,
#when using the predicted probabilities of presence, and the AUC and TSS of the simulated data (are probabilities of presence well separated from presence/absences)
AUC_prob_cov=TSS_prob_cov=AUC_bin=TSS_bin=AUC_simul=TSS_simul=vector(length=ncol(Y))
names(AUC_prob_cov)=names(TSS_prob_cov)=names(AUC_bin)=names(TSS_bin)=colnames(Y)

for(j in 1:ncol(Y)){
  e_joint = evaluate(p=p.mean.stan_prob_cov[which(Y[,j]==1),j], a=p.mean.stan_prob_cov[which(Y[,j]==0),j])
  AUC_prob_cov[j]=e_joint@auc
  TSS_prob_cov[j]=max(e_joint@TPR+e_joint@TNR-1)
  
  e_joint = evaluate(p=p.mean.stan_bin[which(Y[,j]==1),j], a=p.mean.stan_bin[which(Y[,j]==0),j])
  AUC_bin[j]=e_joint@auc
  TSS_bin[j]=max(e_joint@TPR+e_joint@TNR-1)
  
  e_joint = evaluate(p=prob[which(Y[,j]==1),j], a=prob[which(Y[,j]==0),j])
  AUC_simul[j]=e_joint@auc
  TSS_simul[j]=max(e_joint@TPR+e_joint@TNR-1)
  
}

AUC_bin
AUC_prob_cov
TSS_bin
TSS_prob_cov
#when using the threshold, both AUC and TSS are 1 (of course)
AUC_simul
TSS_simul

########################################################################################################
####### Test vs SDM


