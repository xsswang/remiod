
library("remiod")

Sys.setenv(IS_CHECK = "true")

run_model <- function(){
  data(schizow)
  test_j = remiod(formula = y6 ~ tx + y0 + y1 + y3, data=schizow, trtvar = 'tx',
                  model_order = NULL,  method="MAR", ord_cov_dummy = FALSE,
                  algorithm = 'jags', n.adapt = 50, n.chains = 2, include = FALSE,
                  n.iter = 50, mess=FALSE, warn=FALSE, seed=35103)
  return(summary(object = test_j, outcome = c("y6","y1")))
}

model <- run_model()

test_that("summary run", {expect_s3_class(model, "summary.remiod")})

expect_equal(length(model$res), length(model$outcome))

expect_equal(names(model$res), model$outcome)

expect_equal(sum(colnames(model$res[[1]]$regcoef) %in% c("GR-crit","MCE/SD")), 2)

Sys.setenv(IS_CHECK = "")
