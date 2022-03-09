
library("remiod")

Sys.setenv(IS_CHECK = "true")

run_model <- function(){
  data(schizow)

  test_j = remiod(formula = y6 ~ tx + y0 + y1 + y3, data=schizow, trtvar = 'tx',
                  model_order = NULL,  method="MAR", ord_cov_dummy = FALSE,
                  algorithm = 'jags', n.adapt = 10, n.chains = 1, include = FALSE,
                  n.iter = 50, mess=FALSE, warn=FALSE, seed=35103)
  return(test_j)
}

model <- run_model()

test_that("models run", {expect_s3_class(model, "remiod")})

test_that("MCMC is mcmc.list", {expect_s3_class(model$mc.mar$MCMC, "mcmc.list")})

test_that("No Variable has missing value", {
  res <- apply(model$mi.data, 2, function(x) sum(is.na(x))==0)
  testthat::expect_true(all(res), label = "None contains NA")
})

Sys.setenv(IS_CHECK = "")
