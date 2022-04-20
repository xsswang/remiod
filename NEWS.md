# remiod 1.0.1


## New features
* `Ordinal Probit Model`: new model option under `algorithm='jags'`. The option
  can be requested with argument `models = 'opm'` (or named vectors, for example,
  `models = c(y1='opm',y2='opm')`). 
* `algorithm='tang_seq'`: computation is speed up.
* `get_subset`: function is exported for extracting specified set of parameters
  from MCMC samples. The extracted samples have a class of 'mcmc', and can be 
  directly used by functions/packages which take a 'mcmc' object.


## Minor improvements and bug fixes
* Bug fix: when using a J2R imputation of continuous outcome through GLM
  with `family = gaussian()` using JAGS.



