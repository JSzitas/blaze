remove(list = ls())
pkgload::load_all(compile = TRUE)

p <- 3
d <- 1
q <- 1
P <- 1
D <- 0
Q <- 0
season <- 10
use_mean <- TRUE
use_drift <- FALSE
transform <- FALSE

train <- c(lynx)[1:100]
test <- lynx[101:114]

# train[ c(4,7,26)] <- NA


arima_obj <- new(BlazeArima, train, c(p, d, q, P, D, Q, season), list(), "Gardner", "CSS", c(use_mean, use_drift, transform), 1000000)
arima_obj$fit()

arima_obj$forecast(14, list()) -> cpp_fcst

arima_mod <- forecast::Arima(ts(train,frequency = season),
  order = c(p, d, q), seasonal = list(order = c(P, D, Q)),
  include.constant = use_mean,
  include.drift = use_drift,
  transform.pars = transform, kappa = 1000000, method = "CSS"
)
forecast::forecast(arima_mod, 14) -> r_fcst

plot.ts(c(lynx)[91:114])
lines(c( rep(NA, 10), cpp_fcst$forecast), col = "red")
lines(c( rep(NA, 10), r_fcst$mean), col = "green")
