remove(list = ls())
pkgload::load_all(compile = TRUE)

p <- 3
d <- 1
q <- 2
P <- 1
D <- 1
Q <- 1
season <- 10
use_mean <- TRUE

train <- c(lynx)[1:100]
test <- lynx[101:114]

# train[ c(4,7,26)] <- NA


arima_obj <- new(BlazeArima, train, c(p, d, q, P, D, Q, season), list(), "Gardner", "CSS-ML", c(use_mean, TRUE), 1000000)
arima_obj$fit()

arima_obj$forecast(14, list()) -> cpp_fcst

arima_mod <- arima(train,
  order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = season),
  include.mean = use_mean,
  transform.pars = TRUE, kappa = 1000000, method = "CSS-ML"
)
predict(arima_mod, 14) -> r_fcst

plot.ts(test)
lines(cpp_fcst$forecast, col = "red")
lines(c(r_fcst$pred), col = "green")
