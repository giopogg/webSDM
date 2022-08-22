##########################################################################################################
############ Test tSDM by Giovanni Poggiato

## This scripts simply simulates from our trophic model and tests whether our functions can retrieve the model parameters and predictions
##########################################################################################################


# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with metanetwork.  If not, see <http://www.gnu.org/licenses/>



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


wd="~/Documents/GitHub/trophicSDM/Sims_from_the_model/"
setwd(wd)

source("../tSDM_functions.R") # This should be updated to call the function from the package 
n=1000 #number of sites
S=6 #number of species
#This is to chose if we want to include a quadratic term for X_2 or not
linear=T
# if yes, chose whether we want a raw polynomial for X_2 (poly=T) or an orthogonal one (poly=F)
if(!linear) poly=F 

# This is to choose how we create data:
# "true logit" : the correct way of drawing from a logit regression. The linear term (V), that has no noise, is transformed with the inverse logit function to obtain the probabilities (prob), then Y is sample randomly from a bernoulli with probability prob : https://stats.stackexchange.com/questions/46523/how-to-simulate-artificial-data-for-logistic-regression/46525
# the problem of this approach is that it is hard (when the regression coefficients are randomly chose, as here, to obtain probabilities of presence for presences that are well separated from absences. Indeed we have very low TSS and AUC
simul_type="true_logit"


#number of covariates (including polynomial terms)
if(linear){ K=2 }else{ K=3 }#number of covariates (including polynomial terms) }


# The following lines simulate from the model with random coefficients and random graph. Then, they compare it to the true coefficients. We repeate the following 100 times and plot all the results together.

estimates_stan_all = data.frame()
for(s in 1:100){
  
  print(s)
  ### Build the DAG !

  A=matrix(0,nrow=S,ncol=S)

  while(!is.connected(graph_from_adjacency_matrix(A))){
    A=A_w=matrix(0,nrow=S,ncol=S)
    for(j in 1:(S-1)){A[j,c((j+1):S)]=sample(x=c(0,1),size=S-j,replace=T)}
  }

  # transpose A since we want high numbers to be top predators (just for a matter of visualization)
  A=t(A) 

  # Do a loop to create simulated data:
  # we create data until we have that each species has a prevalence between 0.1 and 0.9 (if this is the case, cond   =T, last line of the while loop)
  cond=F
  k=0 #just a counter of loops
  while(!cond){
  k=k+1 #update counter

  #create arrays for data
  Y=V=prob=matrix(0,nrow=n,ncol=S)

  #create weights of the network (i.e. the regression coefficients)
  nonzeros=which(A!=0)
  for(i in nonzeros) A_w[i]=runif(1,min=-2,max=2) #min and max can be tuned (or even sampled)

  # sample environmental covariates
  X_1= runif(n,max=1)
  X_2= runif(n,max=1)

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


  #sample data from the model for each species
  for(j in 1:S){
  
    # Y directly influences the 'link' layer, from which we compute the probabiliy and we sample the realised occurrences.
    V[,j]=X%*%B[j,]+Y%*%A_w[j,] 
    prob[,j]=1/(1+exp(-V[,j]))
    Y[,j]=rbinom(n=n,size=1,prob=prob[,j])

  }


#We check that each species prevalence has to be between 0.1 and 0.9

  if(length(which(colSums(Y)/(n) <0.9 & colSums(Y)/(n)>0.1)) ==S) { cond = T}
  }


  summary(Y)
  summary(prob)

  colnames(Y)=colnames(A)=rownames(A)=paste0("Y",as.character(1:S))

  #build graph from A
  G=graph_from_adjacency_matrix(A)

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


  # fit the model
  m_stan = trophicSDM(Y = Y, X = X, G = G, env.formula = env.form, 
                    penal = NULL, method = "stan_glm", family = binomial(link = "logit"),
                    iter = 1000, run.parallel = F)

  # show traceplots for m_stan (also to check if there are some identification issues)
  # for(j in 1:S){
  #   traceplot(mcmc(as.matrix(m_stan$model[[paste0("Y",j)]])))
  # }

###############
# Analyse results: obtain posterior means and 95% credible regions for each parameter of each species, and compare against the true ones

  estimates_stan = data.frame(p.est=double(),
                              est.97=double(),
                              est.02=double(),
                              true=double(),
                              type=factor(levels = c("biotic","abiotic")),
                              sp.id=factor(),
                              rep = double())

  for(j in 1:S){
    estimates_stan=
      rbind(estimates_stan,data.frame(p.est=coef(m_stan$model[[paste0("Y",j)]]),
                                    est.02=posterior_interval(m_stan$model[[paste0("Y",j)]],prob = 0.95)[,1],
                                    est.97=posterior_interval(m_stan$model[[paste0("Y",j)]],prob = 0.95)[,2],
                                    true=c(B[j,],A_w[j,which(A_w[j,]!=0)]),
                                    type=as.factor(c(rep("abiotic",K+1),rep("biotic",length(which(A_w[j,]!=0))))),
                                    sp.id=as.factor(rep(j,length(coef(m_stan$model[[paste0("Y",j)]])))),
                                    rep = rep(s,length(coef(m_stan$model[[paste0("Y",j)]])))
                                    ))
  
                                            
  }
  estimates_stan_all = rbind(estimates_stan_all,  estimates_stan)

# plot inferred vs true coefficients
# p1=ggplot(data=estimates_stan,aes(x=true,y=p.est,colour=type,shape=sp.id))+ geom_point()+
#   scale_shape_manual(values = c(1:2,4:(S+1))) +
#   geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type))+
#   geom_abline(intercept = 0, slope = 1, colour = "#66B2FF" , linetype= "dashed")+
#   #coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+ 
#   ggtitle("stan_glm()")+
# xlab("True") + ylab("Inferred")+theme_minimal()


# grid.arrange(p1,
#              ggnet2(intergraph::asNetwork(G),arrow.gap = 0.05,arrow.size = 10,label=TRUE)+ggtitle("Known graph"),
#              ncol=2)
# ggsave("InferredVsSim.pdf")

}

# Analyse the results

p1 = ggplot(data=estimates_stan_all,aes(x=true,y=p.est,col=type))+ geom_point(alpha=0.5)+
     geom_linerange(mapping=aes(ymin=est.02, ymax=est.97,colour=type),alpha=0.5)+
     geom_abline(intercept = 0, slope = 1, colour = "#66B2FF" , linetype= "dashed")+
     #coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+ 
     ggtitle("Inferred versus true coefficients")+
     xlab("True") + ylab("Inferred")+theme_minimal()
  
ggsave("Sim_from_model_100.pdf")
length(which(estimates_stan_all$true < estimates_stan_all$est.97 & estimates_stan_all$true > estimates_stan_all$est.02))/ nrow(estimates_stan_all)

