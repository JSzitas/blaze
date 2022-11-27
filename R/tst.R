

bench::press( .grid = expand.grid(diffs = c(1,3,5,8, 10),
                                  lags = c(1,3,5,8,10,12)),
              {
                x <- rnorm(200)
                bench::mark(r = diff(x, lags, diffs),
                            cpp = diff_(x, lags, diffs),
                            cpp2 = diff2_(x, lags, diffs),
                            iterations = 1000
                            )
                }
              ) -> rs
