#' BlazeArima
#'
#' @description Blaze Arima Model Interface
#'
#' @param y The time series to model; an ordered numeric vector or a time-series object.
#' @param pdq The p,d,q order of the model, defaults to c(1,0,0).
#' @param PDQ The PDQ (seasonal) order of the model, defaults to c(0,0,0).
#' @param season The seasonal period of the time series, defaults to 1
#' @param include_mean Whether to include a mean in the estimation, defaults to TRUE.
#' @param include_drift Whether to include drift in the estimation, defaults to FALSE.
#' @param fitting_method The fitting method to use, one of c("CSS-ML", "CSS", "ML", "LASSO").
#' See details.
#' @param xreg A list of exogenous regressors, by default an empty list, please do not set to NULL.
#' @param ss_init The steady state initialization approach, one of c("Gardner", "Rossignol").
#' @param transform_parameters Whether to transform parameters when estimating the State Space model,
#' only impactful for solvers 'CSS-ML' and 'ML'. Defaults to TRUE.
#' @param kappa the prior variance (as a multiple of the innovations variance) for the past observations in a differenced model.
#' Do not reduce this - set to a high value for a reason. Defaults to 1000000.
#' @return A fitted model.
#'
#' @details TODO.
#' @export
#' @importFrom methods new
#' @importFrom stats is.ts
arima <- function(y, pdq = c(1,0,0), PDQ = c(0,0,0),
                  season = 1, include_mean = TRUE, include_drift = FALSE,
                  fitting_method = c("CSS-ML", "CSS", "ML", "LASSO"),
                  xreg = list(), ss_init = c("Gardner", "Rossignol"),
                  transform_parameters = TRUE, kappa = 1000000 ) {
  # cast time series to raw vectors
  # if(stats::is.ts(y)) { y <- c(y) }
  # ss_init <- match.arg(ss_init)
  # fitting_method <- match.arg(fitting_method)
  arima_model <- new(
    blaze_arima,
    c(y), c(pdq, PDQ, season), xreg, ss_init[1], fitting_method[1],
    c(include_mean, include_drift, transform_parameters), kappa)
  arima_model$fit()
  return(structure(list(model = arima_model), class = "BlazeArima"))
}
#' @importFrom generics forecast
#' @exportS3Method
forecast.BlazeArima <- function(object, h = 10, newxreg = list(), ...) {
  return(object$model$forecast(h,newxreg))
}
#' @importFrom stats fitted
#' @exportS3Method
fitted.BlazeArima <- function(object, ...) {
  return(object$model$fitted())
}
#' @importFrom stats residuals
#' @exportS3Method
residuals.BlazeArima <- function(object, ...) {
  return(object$model$residuals())
}
