#' BlazeINoise
#'
#' @description Blaze Integrated Noise Model Interface
#'
#' @param y The time series to model; an ordered numeric vector or a time-series object.
#' @param differences The regular differencing order, by default 1.
#' @param seasonal_difference The seasonal difference, by default 10.
#' @param stabilize Whether to apply a variance stabilizing transform, defaults to TRUE.
#' @return A fitted model.
#'
#' @details TODO.
#' @export
#' @importFrom methods new
inoise <- function(y, differences = 1, seasonal_difference = 10, stabilize = TRUE) {
  # cast time series to raw vectors
  inoise_model <- new(blaze_inoise, c(y), differences, seasonal_difference,
                     stabilize)
  inoise_model$fit()
  return(structure(list(model = inoise_model), class = "BlazeINoise"))
}
#' @importFrom generics forecast
#' @exportS3Method
forecast.BlazeINoise <- function(object, h = 10, produce_se = TRUE,
                                 samples = 10000) {
  return(object$model$forecast(h, produce_se, samples))
}
