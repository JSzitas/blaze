## usethis namespace: start
#' @useDynLib blaze, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

## ensure module gets loaded
Rcpp::loadModule("blaze_arima", TRUE)
Rcpp::loadModule("blaze_ar", TRUE)
Rcpp::loadModule("blaze_auto_ar", TRUE)
Rcpp::loadModule("blaze_inoise", TRUE)

