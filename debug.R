source("R/arima_debug.R")
source("R/debug_make_arima.R")

arima2(lynx, order = c(1,4,2))

Rcpp::sourceCpp("src/test.cpp")
do.call( try_make_arima, debug_arima[[2]][1:3])

debug_arima <- do.call(debug_make_arima, debug_arima[[2]][1:3])


# checked:
# print(debug_arima[[1]]$V)
# print(debug_arima[[1]]$Z)
# print(debug_arima[[1]]$a)
# print(debug_arima[[1]]$P)
# TODO:
# print(debug_arima$Pn)
