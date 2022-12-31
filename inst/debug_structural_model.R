pkgload::load_all(compile = TRUE)


arima_obj <- new( BlazeArima, c(lynx)[1:100], c(3,1,1,0,0,0,10), list(), "Gardner", c(TRUE,TRUE), 1000000  )
arima_obj$fit()
cpp_mod <- arima_obj$get_structural_model()

source("inst/arima_debug.R")
arima2(c(lynx)[1:100], order = c(3,1,1), seasonal = list( order = c(0,0,0), period = 10),
       include.mean = TRUE,
       transform.pars = TRUE, kappa = 1000000, method = "CSS")

