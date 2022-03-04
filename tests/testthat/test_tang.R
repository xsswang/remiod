
library("remiod")

Sys.setenv(IS_CHECK = "true")

run_model <- function(){
  data(schizo)
  schizow = reshape2::dcast(schizo, id + tx ~ week, value.var = "imps79o")
  colnames(schizow) = c(colnames(schizow)[1:2], paste0("y",colnames(schizow)[-c(1:2)]))
  
  schizow[,colnames(schizow)[-c(1:2)]] = lapply(schizow[,colnames(schizow)[-c(1:2)]], 
                                                function(x) factor(x, levels = c("1", "2", "3", "4"), ordered=TRUE))
  schizow1 = subset(schizow, select=c("id","tx","y0","y1","y3","y6"))
  
  test_j = remiod(formula = y6 ~ tx + y0 + y1 + y3, data=schizow1, trtvar = 'tx', 
                  model_order = NULL,  method="MAR", ord_cov_dummy = FALSE,
                  algorithm = 'tang_seq', n.adapt = 10, n.chains = 1, include = FALSE,
                  n.iter = 50, mess=FALSE, warn=FALSE, seed=35103)
  return(test_j)
}

model <- run_model()

test_that("models run", {expect_s3_class(model, "remiod")})

test_that("MCMC is mcmc.list", {expect_s3_class(model$mc.mar$MCMC, "mcmc.list")})

test_that("No Variable has missing value", {
  res <- lapply(model$mi.data, function(x) sum(is.na(x))==0)
  testthat::expect_true(all(unlist(res)), label = "None contains NA")
})

Sys.setenv(IS_CHECK = "")