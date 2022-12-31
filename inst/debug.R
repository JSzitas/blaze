# debug initialization and fitting
Rcpp::sourceCpp("src/export.cpp")

arima_obj <- new( BlazeArima, c(lynx), c(2,0,1,0,0,0,1), list(), "Gardner", c(TRUE,TRUE), 1000000  )
print("Cpp coef: ")
arima_obj$fit()
arima_obj$get_coef()

print("Intended (R) coef: ")
coef(arima(c(lynx), order = c(2,0,1), seasonal = list( order = c(0,0,0), period = 1),
            include.mean = TRUE,
            transform.pars = TRUE, kappa = 1000000, method = "CSS"))

arima_obj2 <- new( BlazeArima, c(lynx), c(4,1,2,0,0,0,10), list(), "Gardner", c(TRUE,TRUE), 1000000  )
print("Cpp coef: ")
arima_obj2$fit()
arima_obj2$get_coef()

print("Intended (R) coef: ")
coef(arima(c(lynx), order = c(4,1,2), seasonal = list( order = c(0,0,0), period = 10),
            include.mean = TRUE,
            transform.pars = TRUE, kappa = 1000000, method = "CSS"))

arima_obj3 <- new( BlazeArima, c(lynx), c(3,0,1,2,0,0,10), list(), "Gardner", c(TRUE,TRUE), 1000000  )
print("Cpp coef: ")
arima_obj3$fit()
arima_obj3$get_coef()

print("Intended (R) coef: ")
coef(arima(c(lynx), order = c(3,0,1), seasonal = list( order = c(2,0,0), period = 10),
            include.mean = TRUE,
            transform.pars = TRUE, kappa = 1000000, method = "CSS"))

arima_obj4 <- new( BlazeArima, c(lynx), c(3,1,1,2,1,2,10), list(), "Gardner", c(TRUE,TRUE), 1000000  )
print("Cpp coef: ")
arima_obj4$fit()
arima_obj4$get_coef()

print("Intended (R) coef: ")
coef(arima(c(lynx), order = c(3,1,1), seasonal = list( order = c(2,1,2), period = 10),
           include.mean = TRUE,
           transform.pars = TRUE, kappa = 1000000, method = "CSS"))


arima_obj5 <- new( BlazeArima, c(lynx), c(3,2,1,3,0,1,12), list(), "Gardner", c(TRUE,TRUE), 1000000  )
print("Cpp coef: ")
arima_obj5$fit()
print(arima_obj5$get_coef())

print("Intended (R) coef: ")
coef(arima(c(lynx), order = c(3,2,1), seasonal = list( order = c(3,0,1), period = 12),
           include.mean = TRUE,
           transform.pars = TRUE, kappa = 1000000, method = "CSS"))


arima_obj6 <- new( BlazeArima, c(lynx), c(2,3,0,0,0,0,1), list(), "Gardner", c(TRUE,TRUE), 1000000  )
print("Cpp coef: ")
arima_obj6$fit()
print(arima_obj6$get_coef())


print("Intended (R) coef: ")
coef(arima(c(lynx), order = c(2,3,0), seasonal = list( order = c(0,0,0), period = 1),
           include.mean = TRUE,
           transform.pars = TRUE, kappa = 1000000, method = "CSS"))

# arima_obj <- new( BlazeArima, c(lynx), c(2,0,1,0,0,0,1), list(), "Gardner", c(TRUE,TRUE), 1000000  )
# # print("Cpp coef: ")
# arima_obj$fit()
# arima_obj$get_coef()
#
# arima_obj$get_structural_model()
#
# res = arima_obj$forecast(12, list())
# res






