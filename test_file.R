library(igraph)
library(parallel)
library(rstanarm)
library(dplyr)
library(rstan)
library(rstanarm)
library(arm)
library(glmnet)
library(formula.tools)
library(abind)
library(brms)
library(piecewiseSEM)
library(stringr)
library(broom)
library(jtools)
library(ggstance)
library(ggplot2)
library(gridExtra)
#library(devtools)
#install_github("MarcOhlmann/metanetwork")
library(metanetwork)
# test functions


load("~/Desktop/SimulData.RData")
for(j in list.files("~/Documents/GitHub/trophicSDM/trophicSDM/R")) {
  source(paste0("~/Documents/GitHub/trophicSDM/trophicSDM/R/",j))
}

env.formula = env.form
penal = NULL
method = "stan_glm"
family = binomial(link = "logit")
iter = 100
chains = 2
run.parallel = F
sp.formula = NULL
mode = "prey"
sp.partition = NULL
#sp.partition = list(c("Y1","Y2","Y3"),"Y4",c("Y5","Y6"))
verbose=F
# sp.formula = list()
# for(j in 1:ncol(Y)){
#   if(length(names(neighbors(G,colnames(Y)[j],mode="out")))>0){
#   sp.formula[[j]] = paste(paste0("X_1*",names(neighbors(G,colnames(Y)[j],mode="out"))),collapse = "+")
#   }
# }
# names(sp.formula) = colnames(Y)

colnames(Y) = gsub("Y","S", colnames(Y))
G <- set.vertex.attribute(G, "name", value=colnames(Y))
names(env.form) = colnames(Y)
env.formula = env.form

tSDM = trophicSDM(Y = Y, X = X, G = G, env.formula = env.form, sp.formula = NULL,
                    sp.partition = NULL,
                    penal = "elasticnet", method = "glm", family = binomial(link = "logit"),
                    iter = 100, run.parallel = F, verbose=F)

evaluateModelFit(tSDM)
plotG(tSDM)
aa = predictFundamental(tSDM, fullPost =F)

class(tSDM)

# Check for every combinations. i.e. method = "glm"& penal = NULL/horshoe or
# method = "stan_glm", penal = NULL/horshoe/coeff.signs()
SDM = tSDM$model$S6
class(SDM)
SDM$coef
coef(SDM)
plot(SDM)
coef(SDM, standardise = T)
summary(tSDM)
SDM

plot(tSDM)

m = m_stan
Xnew = NULL
binary.resp = T
prob.cov = F
pred_samples = 1000
run.parallel = T
verbose = F
filter.table = NULL
# filter.table = list()
# for(j in 1:ncol(Y)) filter.table[[j]] = rbinom(nrow(X),1,0.5)
# names(filter.table) = colnames(Y)
fullPost = T

aa = predict(m_stan, prob.cov = F, fullPost = F)


