remove(list=ls())
pkgload::load_all(compile = TRUE)

p <- 2
d <- 0
q <- 1
P <- 0
D <- 0
Q <- 0
season <- 1

drift_xreg <- seq_len(length(lynx))
x_reg <- list( xreg = cos(sqrt(lynx ^ 2)), drift = drift_xreg )
train_x_reg <- list( x_reg[[1]][1:100], x_reg[[2]][1:100] )
test_x_reg <- list( x_reg[[1]][101:114], x_reg[[2]][101:114] )

arima_obj <- new( BlazeArima, c(lynx)[1:100], c(p,d,q,P,D,Q,season), train_x_reg, "Gardner", c(TRUE,TRUE), 1000000  )
arima_obj$fit()
# arima_obj$get_coef()
arima_obj$forecast(14, test_x_reg) -> cpp_fcst

arima_mod <- arima(c(lynx)[1:100], order = c(p,d,q), seasonal = list( order = c(P,D,Q), period = season),
                   include.mean = TRUE,
                   xreg = as.matrix( as.data.frame(train_x_reg)),
                   transform.pars = TRUE, kappa = 1000000, method = "CSS")
predict(arima_mod, 14, newxreg =  as.matrix( as.data.frame(test_x_reg))) -> r_fcst

plot.ts(lynx[101:114])
lines(cpp_fcst$forecast, col = "red")
lines(c(r_fcst$pred), col = "green")
