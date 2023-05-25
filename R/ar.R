#' BlazeAR
#'
#' @description Blaze Autoregressive Model Interface
#'
#' @param y The time series to model; an ordered numeric vector or a time-series object.
#' @param p The lag order of the model - a non-negative integer. 
#' @param include_mean Whether to include a mean in the estimation, defaults to TRUE.
#' @param include_drift Whether to include drift in the estimation, defaults to FALSE.
#' @param xreg A list of exogenous regressors, by default an empty list, please do not set to NULL.
#' @return A fitted model.
#'
#' @details TODO.
#' @export
#' @importFrom methods new
#' @importFrom stats is.ts
ar <- function(y, p = 1, include_mean = TRUE, include_drift = FALSE,
               xreg = list()) {
  ar_model <- new( blaze_ar, c(y), c(p), xreg, include_mean, include_drift)
  ar_model$fit()
  return(structure(list(model = ar_model), class = "BlazeAR"))
}
#' @importFrom generics forecast
#' @exportS3Method
forecast.BlazeAR <- function(object, h = 10, newxreg = list(), ...) {
  return(object$model$forecast(h,newxreg))
}
#' @importFrom stats fitted
#' @exportS3Method
fitted.BlazeAR <- function(object, ...) {
  return(object$model$fitted())
}
#' @importFrom stats residuals
#' @exportS3Method
residuals.BlazeAR <- function(object, ...) {
  return(object$model$residuals())
}
