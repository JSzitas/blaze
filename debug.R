source("R/arima_debug.R")

arima2(lynx, order = c(2,0,2))

Rcpp::sourceCpp("src/new_arima.cpp")

.Call( stats:::C_ARIMA_transPars,
       debug_arima_trans_pars[[1]],
       debug_arima_trans_pars[[2]],
       debug_arima_trans_pars[[3]])

arima_transform_parameters( debug_arima_trans_pars[[1]],
                            debug_arima_trans_pars[[2]],
                            debug_arima_trans_pars[[3]])

