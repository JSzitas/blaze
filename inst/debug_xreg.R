# debugging file for handling xreg
remove(list=ls())
pkgload::load_all(compile=TRUE)

d <- 0
D <- 0
season <- 1
intercept = TRUE

drift_xreg <- seq_len(length(lynx))
x_reg <- list( xreg = cos(sqrt(lynx ^ 2))[1:100], drift = drift_xreg[1:100] )

res_coef <- test_xreg(c(lynx)[1:100], x_reg, d, D, season, intercept)
print(res_coef)

