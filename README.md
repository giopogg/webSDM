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
# Run to install webSDM
install.packages('webSDM', repos = "http://cran.us.r-project.org")
library(webSDM)
library(ggplot2)
library(rstanarm)
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
## (Intercept)  0.3475526  0.3475526  0.3475526
## X_1         -0.2009638 -0.2632488 -0.1446234
## X_2          0.2585566  0.1934955  0.3192463
## 
## $Y2
##                estimate         5%         95%
## (Intercept) -0.25840889 -0.2584089 -0.25840889
## X_1          0.38338523  0.3217044  0.44840152
## X_2         -0.07590404 -0.1331136 -0.02143486
## 
## $Y3
##               estimate         5%        95%
## (Intercept) -0.6131422 -0.6131422 -0.6131422
## X_1          0.2678163  0.2049716  0.3316497
## X_2          0.1959962  0.1304534  0.2667688
## 
## $Y5
##                estimate          5%         95%
## (Intercept) -0.60887640 -0.60887640 -0.60887640
## X_1          0.20226902  0.13688919  0.26445533
## X_2         -0.30719639 -0.36633457 -0.24113583
## Y1           0.44107544  0.37239328  0.50833208
## Y2          -0.02521752 -0.08470379  0.03257408
## Y3          -0.20786330 -0.26715917 -0.14850745
## 
## $Y4
##                estimate          5%         95%
## (Intercept)  0.78354528  0.78354528  0.78354528
## X_1         -0.17449487 -0.22787550 -0.11480076
## X_2         -0.01111425 -0.07165551  0.04516293
## Y3          -0.07563451 -0.13088922 -0.01919380
## 
## $Y6
##                estimate          5%         95%
## (Intercept)  1.33375272  1.33375272  1.33375272
## X_1         -0.03631624 -0.09686909  0.02993860
## X_2          0.08569472  0.02941836  0.14677667
## Y3          -0.24180669 -0.29984779 -0.17718661
## Y5          -0.10172193 -0.16300415 -0.04541945
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
##         auc       tss species
## 1 0.6732428 0.2776740      Y1
## 2 0.7047877 0.3131087      Y2
## 3 0.6690612 0.2793124      Y3
## 4 0.5985382 0.1766574      Y4
## 5 0.6121942 0.1979403      Y5
## 6 0.5592781 0.1074632      Y6

m$log.lik
## [1] -3638.963
m$AIC
## [1] 7337.927
loo(m)
## [1] -3651.522
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
## Y1 0.6349063        0.6027953        0.664655        
## Y2 0.68883          0.6706062        0.7136599       
## Y3 0.7196283        0.6963811        0.7456787       
## Y5 0.3966708        0.09034863       0.7058777       
## Y4 0.4887261        0.4434229        0.5765817       
## Y6 0.6243699        0.5044867        0.8205877
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
```

All the above-mentioned functions are also available in the frequentist
framework, with adaptations when necessary (e.g. coefficients are
significant or not depending on the p-values instead of the credible
interval). However, error propagation is not available and we can only
obtain one prediction for each species and site, instead of multiple
samples in the Bayesian case.

### Penalisation

We can infer a sparse model by specifying `penal = "horshoe"` if we set
`method = "stan_glm"` (i.e. in the Bayesian framework), or
`penal = "elasticnet"` if we set `method = "glm"` (i.e. in the
frequentist framework).

``` r
m = trophicSDM(Y = Y, X = X, G = G, env.formula = "~ X_1 + X_2",
               family = binomial(link = "logit"),
               mode = "prey", method = "glm", penal = "elasticnet")

m = trophicSDM(Y = Y, X = X, G = G, env.formula = "~ X_1 + X_2",
               family = binomial(link = "logit"),
               mode = "prey", method = "stan_glm", penal = "horshoe")
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
