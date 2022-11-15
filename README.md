R package ‘webSDM’
================
Giovanni Poggiato
09/11/22

# <img src="man/figures/logo_webSDM.png" align="right" width="300px"/>

# webSDM - including known trophic interactions in SDM

webSDM implements the trophic SDM model described in Poggiato et al.,
“Integrating trophic food webs in species distribution models improves
ecological niche estimation and prediction”. Trophic SDM integrate a
known trophic interaction network in SDM, thus providing more realistic
and ecologically sound predictions. We here after present a quick
introduction to some of the feature of the package, but we refer to the
vignettes and the help functions for a more complete documentation of
our package.

## Installing the R package

``` r
library(devtools)
## Le chargement a nécessité le package : usethis
# Run to install webSDM
install_github("giopogg/webSDM")
## Downloading GitHub repo giopogg/webSDM@HEAD
## Skipping 2 packages ahead of CRAN: StanHeaders, rstan
## Installing 1 packages: plyr
## Downloading GitHub repo paul-buerkner/brms@HEAD
## Installing 1 packages: V8
library(webSDM)
## Registered S3 method overwritten by 'GGally':
##   method from   
##   +.gg   ggplot2
library(ggplot2)
library(rstanarm)
## Le chargement a nécessité le package : Rcpp
## This is rstanarm version 2.21.3
## - See https://mc-stan.org/rstanarm/articles/priors for changes to default priors!
## - Default priors may change, so it's safest to specify priors, even if equivalent to the defaults.
## - For execution on a local, multicore CPU with excess RAM we recommend calling
##   options(mc.cores = parallel::detectCores())
```

## Fit a trophic SDM to data

Fitting a trophic SDM require the species distribution data, the
environmental covariates and a known trophic interaction network.  
We load a simulated datasets, where Y contains the species distribution
data (a site x species matrix), X the environmental covariates (a sites
x covariates matrix) and G is an `igraph` object describing the
interaction network (with links going from predators to prey). We
specify, for every species, a quadratic relationship with respect to the
environment, and we fit a model assuming bottom-up control (i.e.,
predators are modeled as a function of preys). We fit the trophicSDM in
the Bayesian framework.

``` r
data(X, Y, G)

m = trophicSDM(Y = Y, X = X, G = G, env.formula = "~ X_1 + X_2",
               family = binomial(link = "logit"),
               mode = "prey", method = "stan_glm")
```

### Inferred coefficients

We can see the formula of each glm in the argument `$form.all`

``` r
m$form.all
## $Y1
## [1] "y ~ X_1 + X_2"
## 
## $Y2
## [1] "y ~ X_1 + X_2"
## 
## $Y3
## [1] "y ~ X_1 + X_2"
## 
## $Y5
## [1] "y ~ X_1+X_2+Y1+Y2+Y3"
## 
## $Y4
## [1] "y ~ X_1+X_2+Y3"
## 
## $Y6
## [1] "y ~ X_1+X_2+Y3+Y5"
```

We can have a first look to regression coefficients using the function
`plot`.

``` r
plot(m)
```

![](man/figures/unnamed-chunk-4-1.png)<!-- -->

We can access to the regression coefficients (eventually standardised)
with the function `coef`.

``` r
coef(m, standardise = T, level = 0.9)
## $Y1
##               estimate         5%        95%
## (Intercept)  0.3579785  0.3579785  0.3579785
## X_1         -0.2046204 -0.2581289 -0.1507982
## X_2          0.2550685  0.1927555  0.3167255
## 
## $Y2
##                estimate         5%         95%
## (Intercept) -0.27418710 -0.2741871 -0.27418710
## X_1          0.38494658  0.3230859  0.45113422
## X_2         -0.07451045 -0.1262927 -0.01844076
## 
## $Y3
##               estimate         5%        95%
## (Intercept) -0.5879114 -0.5879114 -0.5879114
## X_1          0.2645129  0.2007260  0.3262837
## X_2          0.1954173  0.1306161  0.2625365
## 
## $Y5
##                estimate          5%         95%
## (Intercept) -0.62455617 -0.62455617 -0.62455617
## X_1          0.20118853  0.14522254  0.26381736
## X_2         -0.30476466 -0.36915856 -0.24523493
## Y1           0.44149251  0.37271295  0.51996169
## Y2          -0.02182216 -0.07798091  0.03246396
## Y3          -0.20867017 -0.26749779 -0.15107126
## 
## $Y4
##                estimate          5%         95%
## (Intercept)  0.81087991  0.81087991  0.81087991
## X_1         -0.17537046 -0.23262792 -0.11708815
## X_2         -0.01573358 -0.07228362  0.04171549
## Y3          -0.07727005 -0.12803793 -0.01765965
## 
## $Y6
##                estimate          5%         95%
## (Intercept)  1.32794077  1.32794077  1.32794077
## X_1         -0.03391749 -0.09782827  0.02242149
## X_2          0.09007290  0.03310870  0.15240385
## Y3          -0.24091000 -0.30688487 -0.17868017
## Y5          -0.09929429 -0.16268869 -0.03896562
```

We can visualise the biotic coefficients with the function
`plotG_inferred`

``` r
plotG_inferred(m)
```

![](man/figures/unnamed-chunk-6-1.png)<!-- -->

We can also visualise the importance of each (set of) variable with the
function \`computeVariableImportance’

``` r
VarImpo = computeVariableImportance(m, 
                                    groups = list("Abiotic" = c("X_1","X_2"),
                                                  "Biotic" = c("Y1","Y2", "Y3", "Y4", "Y5", "Y6")))
VarImpo = apply(VarImpo, 2, function(x) x/(x[1]+x[2]))
tab = reshape2::melt(VarImpo)
tab$Var2 = factor(tab$Var2, levels = colnames(Y))
ggplot(tab, aes(x = Var2, y = value, fill = Var1)) + geom_bar(stat="identity") +
  theme_classic()
```

![](man/figures/varImpo-1.png)<!-- -->

### Local objects

We can access each local model in the field `$model` and then select a
given local model, that can be analysed using some the implemented
methods

``` r
m$model$Y5
## ==================================================================
## Local SDMfit for species Y5 with no penalty , fitted using stan_glm 
## ==================================================================
## * Useful S3 methods
##     print(), coef(), plot(), predict()
##     $model gives the stanreg class object 
## ==================================================================
```

### Predictions

We can predict with a fitted trophic SDM in multiple ways. The most
straightforward way is to use the function `predict`. We can for example
obtain 50 draws from the posterior predictive distribution of species
(pred_samples = 50) using predicted presence-absences of species to
predict their predators (prob.cov = TRUE). When we don’t specify Xnew,
the function sets Xnew = X by default. We can thus obtain the fitted
values of the model and compute goodness of fit metrics. Notice that
other metrics like AIC or loo are also available.

``` r
Ypred = predict(m, fullPost = FALSE, pred_samples = 50, prob.cov = FALSE)
# predict returns a list contaning for each species the predictive samples at each site
# But since we specified fullPost = FALSE it only give back the predictive mean and quantiles 
Ypred = do.call(cbind,
                lapply(Ypred, function(x) x$predictions.mean))

Ypred = Ypred[,colnames(Y)]
evaluateModelFit(m, Ynew = Y, Ypredicted = Ypred)
##         auc        tss species
## 1 0.6733744 0.28064516      Y1
## 2 0.7049049 0.31464180      Y2
## 3 0.6689794 0.27546580      Y3
## 4 0.5970780 0.16303148      Y4
## 5 0.6196847 0.20630517      Y5
## 6 0.5456740 0.08833226      Y6

m$log.lik
## [1] -3638.561
m$AIC
## [1] 7337.122
loo(m)
## [1] -3650.704
```

We can also evaluate the quality of model predictions using K-fold
cross-validation:

``` r
CV = trophicSDM_CV(m, K = 3, prob.cov = T, run.parallel = FALSE)
# Transfom in a table
Ypred = CV$meanPred
# Re order columns
Ypred = Ypred[,colnames(Y)]
evaluateModelFit(m, Ynew = Y, Ypredicted = Ypred)
```

We can now evaluate species probabilities of presence for the
environmental conditions X_1 = 0.5 and X_2 = 0.5.

``` r
Pred = predict(m, Xnew = data.frame(X_1 = 0.5, X_2 = 0.5), fullPost = F)
t(do.call(cbind, Pred))
##    predictions.mean predictions.q025 predictions.q975
## Y1 0.6347592        0.6021979        0.6771479       
## Y2 0.6858699        0.6443526        0.7173587       
## Y3 0.7195125        0.6943108        0.7541284       
## Y5 0.3715789        0.08845326       0.7382257       
## Y4 0.488158         0.4306166        0.5685854       
## Y6 0.6415191        0.5085278        0.7978324
```

We can also obtain an estimation of the fundamental niche, that
corresponds, in the bottom-up approach, to the probability of presence
of a species given that its preys are present. We can for example
compute the probability of presence of species for the environmental
conditions X_1 = 0.5 and X_2 = 0.5 assuming all their preys to be
present.

``` r
Ypred = predictFundamental(m, fullPost = FALSE, pred_samples = 100, Xnew = data.frame(X_1 = 0.5, X_2 = 0.5))
```

### Frequentist approach

Notice that we can also fit a trophic SDM in the frequentist approach.

``` r
m = trophicSDM(Y = Y, X = X, G = G, env.formula = "~ X_1 + X_2",
               family = binomial(link = "logit"),
               mode = "prey", method = "glm")
## [1] "--- Species Y1 ---"
## [1] "--- Species Y2 ---"
## [1] "--- Species Y3 ---"
## [1] "--- Species Y5 ---"
## [1] "--- Species Y4 ---"
## [1] "--- Species Y6 ---"
```

All the above-mentioned functions are also available in the frequentist
framework, with adaptations when necessary (e.g. coefficients are
significant or not depending on the p-values instead of the credible
interval). However, error propagation is not available and we can only
obtain one prediction for each species and site, instead of multiple
samples in the Bayesian case. ### Penalisation We can infer a sparse
model by specifying `penal = "horshoe"` if we set `method = "stan_glm"`
(i.e. in the Bayesian framework), or `penal = "elasticnet"` if we set
`method = "glm"` (i.e. in the frequentist framework).

``` r
m = trophicSDM(Y = Y, X = X, G = G, env.formula = "~ X_1 + X_2",
               family = binomial(link = "logit"),
               mode = "prey", method = "glm", penal = "elasticnet")
## [1] "--- Species Y1 ---"
## [1] "--- Species Y2 ---"
## [1] "--- Species Y3 ---"
## [1] "--- Species Y5 ---"
## [1] "--- Species Y4 ---"
## [1] "--- Species Y6 ---"

m = trophicSDM(Y = Y, X = X, G = G, env.formula = "~ X_1 + X_2",
               family = binomial(link = "logit"),
               mode = "prey", method = "stan_glm", penal = "horshoe")
## [1] "--- Species Y1 ---"
## [1] "--- Species Y2 ---"
## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#bulk-ess
## [1] "--- Species Y3 ---"
## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#bulk-ess
## [1] "--- Species Y5 ---"
## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#bulk-ess
## [1] "--- Species Y4 ---"
## [1] "--- Species Y6 ---"
## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#bulk-ess
## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#tail-ess
```

### Composite variables

We can include composite variables using the arguments `sp.formula` and
`sp.partition`. We refer to the vignette ‘Composite variable’ for an
exhaustive description of these arguments and how to use them to obtain
any kind of model specification. For example, we can model species as a
function of of their prey richness.

``` r
m = trophicSDM(Y = Y, X = X, G = G, env.formula = "~ X_1 + X_2",
               sp.formula = "richness",
               family = binomial(link = "logit"),
               mode = "prey", method = "glm")
## [1] "--- Species Y1 ---"
## [1] "--- Species Y2 ---"
## [1] "--- Species Y3 ---"
## [1] "--- Species Y5 ---"
## [1] "--- Species Y4 ---"
## [1] "--- Species Y6 ---"
m$form.all
## $Y1
## [1] "y ~ X_1 + X_2"
## 
## $Y2
## [1] "y ~ X_1 + X_2"
## 
## $Y3
## [1] "y ~ X_1 + X_2"
## 
## $Y5
## [1] "y ~ X_1+X_2+I(Y1 + Y2 + Y3)"
## 
## $Y4
## [1] "y ~ X_1+X_2+I(Y3)"
## 
## $Y6
## [1] "y ~ X_1+X_2+I(Y3 + Y5)"
```

## Author

This package is currently developed by Giovanni Poggiato from
Laboratoire d’Ecologie Alpine. It is supported by the ANR GAMBAS. The
framework implemented in this package is described in: “Integrating
trophic food webs in species distribution models improves ecological
niche estimation and prediction” Poggiato Giovanni, Jérémy Andréoletti,
Laura J. Pollock and Wilfried Thuiller. In preparation.
