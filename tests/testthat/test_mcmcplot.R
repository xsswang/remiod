
library("remiod")

Sys.setenv(IS_CHECK = "true")

run_model <- function(){
  data(schizow)
  test_j = remiod(formula = y6 ~ tx + y0 + y1 + y3, data=schizow, trtvar = 'tx',
                  model_order = NULL,  method="MAR", ord_cov_dummy = FALSE,
                  algorithm = 'jags', n.adapt = 50, n.chains = 2, include = FALSE,
                  n.iter = 50, mess=FALSE, warn=FALSE, seed=35103)

  return(mcmcplot(object = test_j, what="trace"))
}

model <- run_model()

test_that("plotting run", {expect_s3_class(model, "ggplot")})

expect_equal(length(unique(model$data$chain)), 2)

Sys.setenv(IS_CHECK = "")
