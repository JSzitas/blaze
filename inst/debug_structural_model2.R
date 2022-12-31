# compare structural models
pkgload::load_all(compile = TRUE)

p <- 3
d <- 1
q <- 1
P <- 3
D <- 1
Q <- 2
season <- 10

arima_obj <- new( BlazeArima, c(lynx)[1:100], c(p,d,q,P,D,Q,season), list(), "Gardner", c(TRUE,TRUE), 1000000  )
arima_obj$fit()

source("inst/arima_debug.R")

arima_mod <- arima2(c(lynx)[1:100],
                    order = c(p,d,q),
                    seasonal = list( order = c(P,D,Q), period = season),
                   include.mean = TRUE,
                   transform.pars = TRUE, kappa = 1000000, method = "CSS")

compare_structural_models <- function( r_mod, cpp_mod, tol = 0.01) {
  assertthat::assert_that(length(cpp_mod) == length(r_mod), msg = "model length unequal")

  names(r_mod) <- names(cpp_mod)

  setNames(
    purrr::map_lgl(names(r_mod),
                   ~ any( abs( c(r_mod[[.x]]) - cpp_mod[[.x]] ) > tol)),
    names(r_mod)
    )
}





