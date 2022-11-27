source("R/arima_debug.R")

arima2(lynx, order = c(1,4,2))

Rcpp::sourceCpp("src/arima.cpp")
# do.call( try_make_arima, debug_arima[[2]][1:3])

seasonal_period <- 10
normal_diff <- 3
seasonal_order_diff <- 2

# r_ts_conv <- function(a,b) {
#   result <- rep(0, length(a) + length(b)-1)
#   for (i in seq_len(length(a))) {
#     for (j in seq_len(length(b))) {
#       result[i + j-1] = result[i+j-1] + a[i] * b[j];
#     }
#   }
#   return(result)
# }
#
# r_ts_conv_n <- function( n_diff = 4) {
#
#   a <- 1
#   b <- c(1,-1)
#   for( k in seq_len(n_diff) ) {
#     temp = rep(0, n_diff+1)
#     for (i in seq_len(k)) {
#       for (j in seq_len(2)) {
#         temp[i + j-1] = temp[i+j-1] + a[i] * b[j];
#       }
#     }
#     a = temp
#   }
#   return(-a[-1L])
# }



# Delta <- 1
# for (i in seq_len(normal_diff)) {
#   Delta <- ts_conv(Delta, c(1, -1))
#   print(Delta)
# }
# for (i in seq_len(seasonal_order_diff)) {
#   Delta <- ts_conv(Delta, c(1, rep.int(0, seasonal_period - 1), -1))
#   print(Delta)
# }
# Delta <- -Delta[-1L]


# debug_arima <- do.call(debug_make_arima, debug_arima[[2]][1:3])


# checked:
# print(debug_arima[[1]]$V)
# print(debug_arima[[1]]$Z)
# print(debug_arima[[1]]$a)
# print(debug_arima[[1]]$P)
# TODO:
# print(debug_arima$Pn)
