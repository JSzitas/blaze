#' BlazeAutoAR
#'
#' @description Blaze Automatic Autoregressive Model Interface
#'
#' @param y The time series to model; an ordered numeric vector or a time-series object.
#' @param min_p The minimum lag order of the model - a non-negative integer.
#' @param max_p The maximum lag order of the model - a non-negative integer. 
#' @param include_mean Whether to include a mean in the estimation, defaults to TRUE.
#' @param include_drift Whether to include drift in the estimation, defaults to FALSE.
#' @param xreg A list of exogenous regressors, by default an empty list, please do not set to NULL.
#' @param fitting_method The method to use; currently one of "AIC", "AICc", "BIC".
#' @return A fitted model.
#'
#' @details TODO.
#' @export
#' @importFrom methods new
#' @importFrom stats is.ts
autoar <- function(y, min_p = 1, max_p = min(length(y)-1, floor(10 * log10(length(y)))),
                   include_mean = TRUE, include_drift = FALSE,
                   xreg = list(), fitting_method = "AIC") {
  ar_model <- new(blaze_auto_ar, c(y), c(min_p), c(max_p), xreg, include_mean,
                  include_drift, fitting_method)
  ar_model$fit()
  return(structure(list(model = ar_model), class = "BlazeAutoAR"))
}
#' @importFrom generics forecast
#' @exportS3Method
forecast.BlazeAutoAR <- function(object, h = 10, newxreg = list(), ...) {
  return(object$model$forecast(h,newxreg))
}
#' @importFrom stats fitted
#' @exportS3Method
fitted.BlazeAutoAR <- function(object, ...) {
  return(object$model$fitted())
}
#' @importFrom stats residuals
#' @exportS3Method
residuals.BlazeAutoAR <- function(object, ...) {
  return(object$model$residuals())
}
