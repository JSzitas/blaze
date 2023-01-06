remove(list=ls())
pkgload::load_all(compile = TRUE)

p <- 2
d <- 1
q <- 1
P <- 0
D <- 0
Q <- 0
season <- 10


arima_obj <- new( BlazeArima, c(lynx)[1:100], c(p,d,q,P,D,Q,season), list(), "Gardner", c(TRUE,TRUE), 1000000  )
arima_obj$fit()
arima_obj$get_coef()

arima_obj$forecast(14, list()) -> cpp_fcst

arima_mod <- arima(c(lynx)[1:100], order = c(p,d,q), seasonal = list( order = c(P,D,Q), period = season),
                   include.mean = TRUE,
                   transform.pars = TRUE, kappa = 1000000, method = "CSS")

predict(arima_mod, 14) -> r_fcst

plot.ts(lynx[101:114])
lines(cpp_fcst$forecast, col = "red")
lines(c(r_fcst$pred), col = "green")



