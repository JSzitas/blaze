remove(list=ls())
pkgload::load_all(compile = TRUE)
source("tst.R")

# lynx_mdl <- blaze::autoar(lynx, 1, 12)
lynx_mdl <- blaze::ar(lynx, 8)
forecast(lynx_mdl, h = 10)
# lynx_mdl$model$get_coef()
predict(ar_ols(lynx, order.max = 8), n.ahead = 10)
# lynx_seasons <- find_seasons(c(lynx))
# print(find_period(c(lynx)))
# print(lynx_seasons)

# elec_seasons <- find_seasons(c(tsibbledata::vic_elec$Demand))
# print(elec_seasons)
# print(find_period(c(tsibbledata::vic_elec$Demand)))



