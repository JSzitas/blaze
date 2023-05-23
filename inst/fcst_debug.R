remove(list = ls())
pkgload::load_all(compile = TRUE)

p <- 2
d <- 1
q <- 1
P <- 1
D <- 1
Q <- 1
season <- 10
use_mean <- FALSE
use_drift <- FALSE
transform <- TRUE

train <- c(lynx)[1:100]
test <- lynx[101:114]

# blaze::ar(train, p = 2) -> cpp
# forecast(cpp, h = 14) -> fcst





# train[ c(4,7,26)] <- NA


arima_cpp <- blaze::arima(train, c(p, d, q), c(P, D, Q), season,
                          use_mean, use_drift, fitting_method = "ML",
                          xreg = list(),
                          ss_init = "Gardner",
                          transform_parameters = transform, 1000000)
forecast(arima_cpp, 14, list()) -> cpp_fcst

arima_r <- forecast::Arima(ts(train,frequency = season),
                           order = c(p, d, q), seasonal = list(order = c(P, D, Q)),
                           include.mean = use_mean, include.drift = use_drift,
                           transform.pars = transform, kappa = 1000000, method = "ML"
)
forecast::forecast(arima_r, 14) -> r_fcst

plot.ts(c(lynx)[91:114])
lines(c( rep(NA, 10), cpp_fcst$forecast), col = "red")
lines(c( rep(NA, 10), r_fcst$mean), col = "green")
