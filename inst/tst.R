

# bench::press( .grid = expand.grid(diffs = c(1,3,5,8, 10),
#                                   lags = c(1,3,5,8,10,12)),
#               {
#                 x <- rnorm(200)
#                 bench::mark(r = diff(x, lags, diffs),
#                             cpp = diff_(x, lags, diffs),
#                             cpp2 = diff2_(x, lags, diffs),
#                             iterations = 1000
#                             )
#                 }
#               ) -> rs

lynx2 <- log(lynx + rnorm(length(lynx), sd = 50)^2)
plot.ts(lynx2)

rs <- arima4( lynx, c(4,3,2), xreg = lynx2 )






