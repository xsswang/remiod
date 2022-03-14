#' wide format of ordinal response of Schizophrenia data.
#'
#' A dataset containing the treatment and ordinal responses measured at baseline and 3 post-baseline visits
#'
#' @format A data frame with 437 rows and 5 variables:
#' \describe{
#'   \item{tx}{treatment, 1 for treated and 0 for placebo}
#'   \item{y0}{ordinal response at the baseline}
#'   \item{y1, y3, y6}{ordinal response at the post-baseline week 1, 3, and 6.}
#' }
#' @source long-to-wise tranformation of schizo data, i.e. schizow = data.table::dcast(schizo, id + tx ~ week, value.var = "imps79o")
"schizow"



#' wide format of binary response of Schizophrenia data.
#'
#' A dataset containing the treatment and binary responses measured at baseline and 3 post-baseline visits
#'
#' @format A data frame with 437 rows and 5 variables:
#' \describe{
#'   \item{tx}{treatment, 1 for treated and 0 for placebo}
#'   \item{y0}{binary response at the baseline}
#'   \item{y1, y3, y6}{binary response at the post-baseline week 1, 3, and 6.}
#' }
#' @source long-to-wise tranformation of schizo data, i.e. schizob = data.table::dcast(schizo, id + tx ~ week, value.var = "imps79b")
"schizob"
