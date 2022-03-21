##########################################################################################################
############ TEST tSDM
##########################################################################################################
rm(list=ls())
library(igraph)
library(ggplot2)
#library(pcalg) not working with R 4.0

wd="~/Documents/Phd/Futureweb/Code/"
setwd(wd)
source('/Users/poggiatg/Documents/Phd/Futureweb/Code/MY_SCRIPTS/tSDM_functions.R')


n=1000 #number of sites
n_oos= 500#number of sites out of sample
S=6 #number of species
#This is to chose if we want to include a quadratic term for X_2 or not
linear=T
# if yes, chose whether we want a raw polynomial for X_2 (poly=T) or an orthogonal one (poly=F)
if(!linear) poly=F 

#number of covariates (including polynomial terms) 
if(linear){ K=2 }else{ K=3 }


### Build the newtork !
A=matrix(0,nrow=S,ncol=S)
#or sample random DAG with pcalg (not working for me)
while(!is.connected(graph_from_adjacency_matrix(A))){
A=A_w=matrix(0,nrow=S,ncol=S)
for(j in 1:(S-1)){A[j,c((j+1):S)]=sample(x=c(0,1),size=S-j,replace=T)}
}
# transpose A since we want high numbers to be top predators (just for a matter of visualization)
A=t(A) 


#create weights (i.e. regression coefficients, -2 and 2 can be tuned)
nonzeros=which(A!=0)
for(i in nonzeros) A_w[i]=runif(1,min=-2,max=2) 
# sample regression coefficients (-5 and 5 can be tuned)
B=matrix(runif(S*(K+1),min=-5,max=5),ncol=K+1) #+1 bc of the intercept

# sample environmental covariates
X_1= runif(n+n_oos,max=10)
X_2= runif(n+n_oos,max=10)

#assamble them depending on what we want
if(linear){
  X=cbind(1,X_1,X_2)
}else{
if(!poly){
  X=cbind(1,X_1,X_2,X_2^2) 
}else{
  X=cbind(1,X_1,poly(X_2,degree=2))
}
}

# random noise
sigma_2 = runif(S,min = 0.5,max=1)

# create empy arrays for data 
Y=V=matrix(0,nrow=n+n_oos,ncol=S)

#sample data for each species
for(j in 1:S){
    Y[,j]=X%*%B[j,]+Y%*%A_w[j,]+rnorm(n+n_oos,0,sqrt(sigma_2[j]))
}

head(Y)

# Important: Y names should not be just raw numbers
colnames(Y)=colnames(A)=rownames(A)=paste0("Y_",as.character(1:S))

#build graph from A
G=graph_from_adjacency_matrix(A)

# Split dataset in train-test
X_test=data.frame(X_1=X_1[(n+1):(n+n_oos)],X_2=X_2[(n+1):(n+n_oos)])
X_raw=data.frame(X_1=X_1[1:n],X_2=X_2[1:n]) 

Y_test=as.data.frame(Y[(n+1):(n+n_oos),])
Y=as.data.frame(Y[1:n,])


# Create formulas for the environmental part of each species
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


###############
# Run models
m_stan = trophicSDM(Y=Y,X=X_raw,G=G,formulas=env.form,penal=NULL,method="stan_glm",sp.formula = NULL,sp.partition = NULL,family=gaussian(),iter=2000)
m_bayes = trophicSDM(Y=Y,X=X_raw,G=G,formulas=env.form,penal=NULL,method="bayesglm",sp.formula = NULL,sp.partition = NULL,family=gaussian(),iter=2000)
m_glm = trophicSDM(Y=Y,X=X_raw,G=G,formulas=env.form,penal=NULL,method="glm",sp.formula = NULL,sp.partition = NULL,family=gaussian(),iter=2000)

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
                                                 est.02=posterior_interval(m_stan$model[[paste0("Y_",j)]],prob = 0.95)[-nrow(posterior_interval(m_stan$model[[paste0("Y_",j)]])),1],
                                                 est.97=posterior_interval(m_stan$model[[paste0("Y_",j)]],prob = 0.95)[-nrow(posterior_interval(m_stan$model[[paste0("Y_",j)]])),2],
                                                 true=c(B[j,],A_w[j,which(A_w[j,]!=0)]),
                                                 type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(which(A_w[j,]!=0))))),
                                                 sp.id=as.factor(rep(j,length(coef(m_stan$model[[paste0("Y_",j)]]))))
  ))
  
  estimates_glm=rbind(estimates_glm,data.frame(p.est=coef(m_glm$model[[paste0("Y_",j)]]),
                                               est.02=confint(m_glm$model[[paste0("Y_",j)]],prob = 0.95)[-nrow(posterior_interval(m_stan$model[[paste0("Y_",j)]])),1],
                                               est.97=confint(m_glm$model[[paste0("Y_",j)]],prob = 0.95)[-nrow(posterior_interval(m_stan$model[[paste0("Y_",j)]])),2],
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

# plot together with the trophic network
#pdf(paste0(wd,"PLOTS/Simulations/Gaussian_simul_firsttry.pdf"))
grid.arrange(arrangeGrob(p1,p2,top="True vs inferred regression coefficients"),
             ggnet2(intergraph::asNetwork(G),arrow.gap = 0.05,arrow.size = 10,label=TRUE)+ggtitle("Known graph"),
             ncol=2)
dev.off()

#############################
# Predictions

# In-sample predictions
Xnew=X_raw  #Xnew are the environmental variables that we used to fit the model (actually if Xnew=NULL, this is the default choice) 
pred_samples=100 # the number of predictive samples we want for each species
error_prop_sample=10 #error_prop_sample is an important parameter. The greater it is, the largest the error propagation (and therefore the predictive credible regions). see tSDM_functions for more details

#run predictions
pred_stan=trophicSDM_predict(m=m_stan,Xnew=Xnew,pred_samples=pred_samples,error_prop_sample=error_prop_sample)
pred_bayes=trophicSDM_predict(m=m_bayes,Xnew,pred_samples=pred_samples,error_prop_sample=pred_samples)
pred_glm=trophicSDM_predict(m=m_glm,Xnew,pred_samples=1,error_prop_sample=1)

p.mean.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,mean, MARGIN = 1))
p.mean.stan=do.call(cbind, p.mean.stan.temp)
colnames(p.mean.stan)= names(topo_sort(G,mode="in"))

p.q975.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,quantile, MARGIN = 1,probs=0.975))
p.q975.stan=do.call(cbind, p.q975.stan.temp)
colnames(p.q975.stan)= names(topo_sort(G,mode="in"))

p.q025.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,quantile, MARGIN = 1,probs=0.025))
p.q025.stan=do.call(cbind, p.q025.stan.temp)
colnames(p.q025.stan)= names(topo_sort(G,mode="in"))


#reorder the columns of Y
Y=Y[,names(topo_sort(G,mode="in"))]

table=data.frame(obs=as.vector(as.matrix(Y)),pred=as.vector(as.matrix(p.mean.stan)),
                 est.02=as.vector(as.matrix(p.q025.stan)),est.97=as.vector(as.matrix(p.q975.stan)),sp.name=rep(colnames(Y),each=nrow(Y)))

ggplot(data=table,mapping=aes(x=obs,y=pred,col=sp.name))+geom_point(alpha=0.3)+
  geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=sp.name),alpha=0.5)+
  geom_abline(slope=1,intercept = 0)+ggtitle(paste0("In-sample predictions stan_glm() \n error_prop_sample = ",error_prop_sample))+
  xlab("Observed") + ylab("predicted")+theme_minimal()+coord_fixed()

# Compute R2 and proportion of confidence interval that match
R2_err_pr=sapply(names(topo_sort(G,mode="in")),function(x){cor(p.mean.stan[,x],Y[,x])^2})
hist(R2_err_pr)


CR_match_err_pr=sapply(names(topo_sort(G,mode="in")), function(x) mean(sapply(1:n,FUN=function(y) 
  ifelse(Y[y,x]>p.q025.stan[y,x] & Y[y,x]<p.q975.stan[y,x],1,0))))
hist(CR_match_err_pr)

# How do change in n_sample/error_pred sample affect the results?
# In-sample predictions
pred_samples=100
error_prop_sample=1
pred_stan=trophicSDM_predict(m=m_stan,Xnew=Xnew,penal=NULL,pred_samples=pred_samples,error_prop_sample=error_prop_sample)


p.mean.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,mean, MARGIN = 1))
p.mean.stan=do.call(cbind, p.mean.stan.temp)
colnames(p.mean.stan)= names(topo_sort(G,mode="in"))

p.q975.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,quantile, MARGIN = 1,probs=0.975))
p.q975.stan=do.call(cbind, p.q975.stan.temp)
colnames(p.q975.stan)= names(topo_sort(G,mode="in"))

p.q025.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,quantile, MARGIN = 1,probs=0.025))
p.q025.stan=do.call(cbind, p.q025.stan.temp)
colnames(p.q025.stan)= names(topo_sort(G,mode="in"))

table=data.frame(obs=as.vector(as.matrix(Y)),pred=as.vector(as.matrix(p.mean.stan)),
                 est.02=as.vector(as.matrix(p.q025.stan)),est.97=as.vector(as.matrix(p.q975.stan)),sp.name=rep(colnames(Y),each=nrow(Y)))

ggplot(data=table,mapping=aes(x=obs,y=pred,col=sp.name))+geom_point(alpha=0.3)+
  geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=sp.name),alpha=0.3)+
  geom_abline(slope=1,intercept = 0)+ggtitle(paste0("In-sample predictions stan_glm() \n error_prop_sample = ",error_prop_sample))+
  xlab("Observed") + ylab("predicted")+theme_minimal()+coord_fixed()

# Compute R2 and proportion of confidence interval that match
R2_no_err_pr=sapply(names(topo_sort(G,mode="in")),function(x){cor(p.mean.stan[,x],Y[,x])^2})
hist(R2_no_err_pr)


CR_match_no_err_pr=sapply(names(topo_sort(G,mode="in")), function(x) mean(sapply(1:n,FUN=function(y) 
  ifelse(Y[y,x]>p.q025.stan[y,x] & Y[y,x]<p.q975.stan[y,x],1,0))))
hist(CR_match_no_err_pr)


##################################################################################
# Out-of-sample predictions
Y_test=Y=Y_test[,names(topo_sort(G,mode="in"))]

pred_samples=100
error_prop_sample=1

pred_stan=trophicSDM_predict(m=m_stan,Xnew=X_test,pred_samples=pred_samples,error_prop_sample=error_prop_sample)



p.mean.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,mean, MARGIN = 1))
p.mean.stan=do.call(cbind, p.mean.stan.temp)
colnames(p.mean.stan)= names(topo_sort(G,mode="in"))

p.q975.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,quantile, MARGIN = 1,probs=0.975))
p.q975.stan=do.call(cbind, p.q975.stan.temp)
colnames(p.q975.stan)= names(topo_sort(G,mode="in"))

p.q025.stan.temp=lapply(pred_stan$sp.prediction,FUN=function(x) apply(x,quantile, MARGIN = 1,probs=0.025))
p.q025.stan=do.call(cbind, p.q025.stan.temp)
colnames(p.q025.stan)= names(topo_sort(G,mode="in"))

table=data.frame(obs=as.vector(as.matrix(Y_test)),pred=as.vector(as.matrix(p.mean.stan)),
                 est.02=as.vector(as.matrix(p.q025.stan)),est.97=as.vector(as.matrix(p.q975.stan)),sp.name=rep(colnames(Y),each=nrow(Y)))

ggplot(data=table,mapping=aes(x=obs,y=pred,col=sp.name))+geom_point(alpha=0.3)+
  geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=sp.name),alpha=0.3)+
  geom_abline(slope=1,intercept = 0)+ggtitle(paste0("Out-of-sample predictions stan_glm() \n error_prop_sample = ",error_prop_sample))+
  xlab("Observed") + ylab("predicted")+theme_minimal()+coord_fixed()


# Compute R2 and proportion of confidence interval that match
R2_no_err_pr=sapply(names(topo_sort(G,mode="in")),function(x){cor(p.mean.stan[,x],Y[,x])^2})
hist(R2_no_err_pr)


CR_match_no_err_pr=sapply(names(topo_sort(G,mode="in")), function(x) mean(sapply(1:n_oos,FUN=function(y) 
  ifelse(Y_test[y,x]>p.q025.stan[y,x] & Y_test[y,x]<p.q975.stan[y,x],1,0))))
hist(CR_match_no_err_pr)


