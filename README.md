# trophicSDM
Combines a given trophic (acyclic!) network with species distribution and environmental data by including preys occurrences as predictors for their predators. Predictions are obtaining starting from basal species which are then used to predict their predators and so on up to top predators.

tSDM_functions.R : the core script with the necessary functions
tSDM_test_gauss.R : simulates gaussian data from the model, fits tSDMs and predicts. Check if parameters are retrived and if predictions coincide.
tSDM_test_logit.R : simulates binary data from the model, fits tSDMs and predicts. Check if parameters are retrived and if predictions coincide.
