#ifndef AUTO_ARIMA
#define AUTO_ARIMA

#include "arima/structures/fitting_method.h"
#include "arima/structures/structural_model.h"
#include "arima/structures/ss_init.h"

#include <cmath>

struct arima_restrictions{
  arima_restrictions(const size_t d = -1,
                     const size_t D = -1,
                     const size_t max_p = 5,
                     const size_t max_q = 5,
                     const size_t max_P = 2,
                     const size_t max_Q = 2,
                     const size_t max_order = 5,
                     const size_t max_d = 2,
                     const size_t max_D = 1,
                     const size_t start_p = 2,
                     const size_t start_q = 2,
                     const size_t start_P = 1,
                     const size_t start_Q = 1,
                     const bool stationary = false,
                     const bool seasonal = true) :
  d(d), D(D), max_p(max_p), max_q(max_q), max_P(max_P), max_Q(max_Q),
  max_order(max_order), max_d(max_d), max_D(max_D), start_p(start_p),
  start_q(start_q), start_P(start_P), start_Q(start_Q), stationary(stationary),
  seasonal(seasonal) {}
private:
  const size_t d, D, max_p, max_q, max_P, max_Q,
  max_order, max_d, max_D, start_p, start_q, start_P, start_Q;
  const bool stationary, seasonal;
};

enum IC{
  AIC,
  BIC,
  AICC
};

template <typename scalar_t = double> class AutoArima{
public:
  AutoArima(std::vector<double> y,
            const arima_restrictions restrictions,
            const IC criterium = IC::AIC,
            const bool stepwise = true,
            const std::vector<std::vector<scalar_t>> xreg = {{}},
            const size_t n_models = 94,
            const fitting_method method = ML,
            const SSinit ss_init = SSinit::Gardner,
            const bool approximate = false,
            const bool try_drift = true,
            const bool try_mean = true,
            const bool try_box_cox = true,
            const bool standardize = true,
            const scalar_t kappa = 1000000) {
    // if(is_constant(y)) {
    //   if(all_is_nan(y)) {
    //     Arima();
    //   }
    // }
  }
  void fit() {

  }



};


{
  if (is.constant(x)) {
    if (all(is.na(x)))
      stop("All data are missing")
      if (allowmean) {
        fit <- Arima(x, order = c(0, 0, 0), fixed = mean(x,
                                  na.rm = TRUE), ...)
      }
      else {
        fit <- Arima(x, order = c(0, 0, 0), include.mean = FALSE,
                                                   ...)
      }
      fit$x <- orig.x
        fit$series <- series
        fit$call <- match.call()
        fit$call$x <- data.frame(x = x)
        fit$constant <- TRUE
        return(fit)
  }

  series <- deparse(substitute(y))
    x <- as.ts(x)

    orig.x <- x
      missing <- is.na(x)
      firstnonmiss <- head(which(!missing), 1)
      lastnonmiss <- tail(which(!missing), 1)
      serieslength <- sum(!missing[firstnonmiss:lastnonmiss])
      x <- subset(x, start = firstnonmiss)

      ic <- match.arg(ic)
        test <- match.arg(test)
        seasonal.test <- match.arg(seasonal.test)
        if (seasonal) {
          m <- frequency(x)
        }
        else {
          m <- 1
        }
        if (m < 1) {
          m <- 1
        }
        else {
          m <- round(m)
        }
        max.p <- min(max.p, floor(serieslength/3))
          max.q <- min(max.q, floor(serieslength/3))
          max.P <- min(max.P, floor(serieslength/3/m))
          max.Q <- min(max.Q, floor(serieslength/3/m))
          if (serieslength <= 3L) {
            ic <- "aic"
          }
          if (!is.null(lambda)) {
            x <- BoxCox(x, lambda)
            lambda <- attr(x, "lambda")
            attr(lambda, "biasadj") <- biasadj
          }
          if (stationary) {
            d <- D <- 0
          }
          if (m == 1) {
            D <- max.P <- max.Q <- 0
          }
          else if (is.na(D)) {
            D <- do.call("nsdiffs", c(list(xx, test = seasonal.test,
                                    max.D = max.D), seasonal.test.args))
              if (D > 0 && !is.null(xregg)) {
                diffxreg <- diff(xregg, differences = D, lag = m)
                if (any(apply(diffxreg, 2, is.constant))) {
                  D <- D - 1
                }
              }
              if (D > 0) {
                dx <- diff(xx, differences = D, lag = m)
                if (all(is.na(dx)))
                  D <- D - 1
              }
          }
          if (D > 0) {
            dx <- diff(xx, differences = D, lag = m)
          }
          else {
            dx <- xx
          }
          if (!is.null(xregg)) {
            if (D > 0) {
              diffxreg <- diff(xregg, differences = D, lag = m)
            }
            else {
              diffxreg <- xregg
            }
          }
          if (is.na(d)) {
            d <- do.call("ndiffs", c(list(dx, test = test, max.d = max.d),
                               test.args))
              if (d > 0 && !is.null(xregg)) {
                diffxreg <- diff(diffxreg, differences = d, lag = 1)
                if (any(apply(diffxreg, 2, is.constant))) {
                  d <- d - 1
                }
              }
              if (d > 0) {
                diffdx <- diff(dx, differences = d, lag = 1)
                if (all(is.na(diffdx)))
                  d <- d - 1
              }
          }
          if (d > 0) {
            dx <- diff(dx, differences = d, lag = 1)
          }
          if (length(dx) == 0L)
            stop("Not enough data to proceed")
            else if (is.constant(dx)) {
              if (is.null(xreg)) {
                if (D > 0 && d == 0) {
                  fit <- Arima(x, order = c(0, d, 0), seasonal = list(order = c(0,
                                                                      D, 0), period = m), include.constant = TRUE,
                                                                      fixed = mean(dx/m, na.rm = TRUE), method = method,
                                                                                           ...)
                }
                else if (D > 0 && d > 0) {
                  fit <- Arima(x, order = c(0, d, 0), seasonal = list(order = c(0,
                                                                      D, 0), period = m), method = method, ...)
                }
                else if (d == 2) {
                  fit <- Arima(x, order = c(0, d, 0), method = method,
                               ...)
                }
                else if (d < 2) {
                  fit <- Arima(x, order = c(0, d, 0), include.constant = TRUE,
                               fixed = mean(dx, na.rm = TRUE), method = method,
                                                  ...)
                }
                else {
                  stop("Data follow a simple polynomial and are not suitable for ARIMA modelling.")
                }
              }
              else {
                if (D > 0) {
                  fit <- Arima(x, order = c(0, d, 0), seasonal = list(order = c(0,
                                                                      D, 0), period = m), xreg = xreg, method = method,
                                                                      ...)
                }
                else {
                  fit <- Arima(x, order = c(0, d, 0), xreg = xreg,
                               method = method, ...)
                }
              }
              fit$x <- orig.x
                fit$series <- series
                fit$call <- match.call()
                fit$call$x <- data.frame(x = x)
                return(fit)
            }
            if (m > 1) {
              if (max.P > 0) {
                max.p <- min(max.p, m - 1)
              }
              if (max.Q > 0) {
                max.q <- min(max.q, m - 1)
              }
            }
            if (approximation) {
              if (!is.null(truncate)) {
                tspx <- tsp(x)
                if (length(x) > truncate) {
                  x <- ts(tail(x, truncate), end = tspx[2], frequency = tspx[3])
                }
              }
              if (D == 0) {
                fit <- try(stats::arima(x, order = c(0, d, 0), xreg = xreg,
                                        ...), silent = TRUE)
              }
              else {
                fit <- try(stats::arima(x, order = c(0, d, 0), seasonal = list(order = c(0,
                                                                               D, 0), period = m), xreg = xreg, ...), silent = TRUE)
              }
              if (!is.element("try-error", class(fit))) {
                offset <- -2 * fit$loglik - serieslength * log(fit$sigma2)
              }
              else {
                offset <- 0
              }
            }
            else {
              offset <- 0
            }
            allowdrift <- allowdrift & (d + D) == 1
            allowmean <- allowmean & (d + D) == 0
            constant <- allowdrift | allowmean

              if (!stepwise) {
                bestfit <- search.arima(x, d, D, max.p, max.q, max.P,
                                        max.Q, max.order, stationary, ic, trace, approximation,
                                        method = method, xreg = xreg, offset = offset, allowdrift = allowdrift,
                                        allowmean = allowmean, parallel = parallel, num.cores = num.cores,
                                                                                       ...)
                                        bestfit$call <- match.call()
                bestfit$call$x <- data.frame(x = x)
                bestfit$lambda <- lambda
                bestfit$x <- orig.x
                bestfit$series <- series
                bestfit$fitted <- fitted.Arima(bestfit)
                return(bestfit)
              }
              if (length(x) < 10L) {
                start.p <- min(start.p, 1L)
                start.q <- min(start.q, 1L)
                start.P <- 0L
                start.Q <- 0L
              }
              p <- start.p <- min(start.p, max.p)
                q <- start.q <- min(start.q, max.q)
                P <- start.P <- min(start.P, max.P)
                Q <- start.Q <- min(start.Q, max.Q)
                results <- matrix(NA, nrow = nmodels, ncol = 8)
                bestfit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                                                D, Q), constant = constant, ic, trace, approximation,
                                                method = method, offset = offset, xreg = xreg, ...)
                results[1, ] <- c(p, d, q, P, D, Q, constant, bestfit$ic)
                fit <- myarima(x, order = c(0, d, 0), seasonal = c(0, D,
                                            0), constant = constant, ic, trace, approximation, method = method,
                                            offset = offset, xreg = xreg, ...)
                results[2, ] <- c(0, d, 0, 0, D, 0, constant, fit$ic)
                if (fit$ic < bestfit$ic) {
                  bestfit <- fit
                  p <- q <- P <- Q <- 0
                }
                k <- 2
                if (max.p > 0 || max.P > 0) {
                  fit <- myarima(x, order = c(max.p > 0, d, 0), seasonal = c((m >
                                                                                1) & (max.P > 0), D, 0), constant = constant, ic,
                                                                                  trace, approximation, method = method, offset = offset,
                                                                                  xreg = xreg, ...)
                  results[k + 1, ] <- c(max.p > 0, d, 0, (m > 1) & (max.P >
                                                                      0), D, 0, constant, fit$ic)
                  if (fit$ic < bestfit$ic) {
                    bestfit <- fit
                    p <- (max.p > 0)
                    P <- (m > 1) & (max.P > 0)
                    q <- Q <- 0
                  }
                  k <- k + 1
                }
                if (max.q > 0 || max.Q > 0) {
                  fit <- myarima(x, order = c(0, d, max.q > 0), seasonal = c(0,
                                              D, (m > 1) & (max.Q > 0)), constant = constant, ic,
                                              trace, approximation, method = method, offset = offset,
                                              xreg = xreg, ...)
                  results[k + 1, ] <- c(0, d, max.q > 0, 0, D, (m > 1) &
                    (max.Q > 0), constant, fit$ic)
                  if (fit$ic < bestfit$ic) {
                    bestfit <- fit
                    p <- P <- 0
                    Q <- (m > 1) & (max.Q > 0)
                    q <- (max.q > 0)
                  }
                  k <- k + 1
                }
                if (constant) {
                  fit <- myarima(x, order = c(0, d, 0), seasonal = c(0,
                                              D, 0), constant = FALSE, ic, trace, approximation,
                                              method = method, offset = offset, xreg = xreg, ...)
                  results[k + 1, ] <- c(0, d, 0, 0, D, 0, 0, fit$ic)
                  if (fit$ic < bestfit$ic) {
                    bestfit <- fit
                    p <- q <- P <- Q <- 0
                  }
                  k <- k + 1
                }
                startk <- 0
                while (startk < k && k < nmodels) {
                  startk <- k
                  if (P > 0 && newmodel(p, d, q, P - 1, D, Q, constant,
                                        results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q), seasonal = c(P -
                        1, D, Q), constant = constant, ic, trace, approximation,
                          method = method, offset = offset, xreg = xreg,
                          ...)
                          results[k, ] <- c(p, d, q, P - 1, D, Q, constant,
                                            fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        P <- (P - 1)
                        next
                      }
                  }
                  if (Q > 0 && newmodel(p, d, q, P, D, Q - 1, constant,
                                        results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                                                  D, Q - 1), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p, d, q, P, D, Q - 1, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        Q <- (Q - 1)
                        next
                      }
                  }
                  if (P < max.P && newmodel(p, d, q, P + 1, D, Q, constant,
                                            results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q), seasonal = c(P +
                        1, D, Q), constant = constant, ic, trace, approximation,
                          method = method, offset = offset, xreg = xreg,
                          ...)
                          results[k, ] <- c(p, d, q, P + 1, D, Q, constant,
                                            fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        P <- (P + 1)
                        next
                      }
                  }
                  if (Q < max.Q && newmodel(p, d, q, P, D, Q + 1, constant,
                                            results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                                                  D, Q + 1), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p, d, q, P, D, Q + 1, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        Q <- (Q + 1)
                        next
                      }
                  }
                  if (Q > 0 && P > 0 && newmodel(p, d, q, P - 1, D, Q -
                      1, constant, results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q), seasonal = c(P -
                        1, D, Q - 1), constant = constant, ic, trace,
                          approximation, method = method, offset = offset,
                          xreg = xreg, ...)
                      results[k, ] <- c(p, d, q, P - 1, D, Q - 1, constant,
                                        fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        Q <- (Q - 1)
                        P <- (P - 1)
                        next
                      }
                  }
                  if (Q < max.Q && P > 0 && newmodel(p, d, q, P - 1, D,
                                                     Q + 1, constant, results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q), seasonal = c(P -
                        1, D, Q + 1), constant = constant, ic, trace,
                          approximation, method = method, offset = offset,
                          xreg = xreg, ...)
                      results[k, ] <- c(p, d, q, P - 1, D, Q + 1, constant,
                                        fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        Q <- (Q + 1)
                        P <- (P - 1)
                        next
                      }
                  }
                  if (Q > 0 && P < max.P && newmodel(p, d, q, P + 1, D,
                                                     Q - 1, constant, results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q), seasonal = c(P +
                        1, D, Q - 1), constant = constant, ic, trace,
                          approximation, method = method, offset = offset,
                          xreg = xreg, ...)
                      results[k, ] <- c(p, d, q, P + 1, D, Q - 1, constant,
                                        fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        Q <- (Q - 1)
                        P <- (P + 1)
                        next
                      }
                  }
                  if (Q < max.Q && P < max.P && newmodel(p, d, q, P + 1,
                                                         D, Q + 1, constant, results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q), seasonal = c(P +
                        1, D, Q + 1), constant = constant, ic, trace,
                          approximation, method = method, offset = offset,
                          xreg = xreg, ...)
                      results[k, ] <- c(p, d, q, P + 1, D, Q + 1, constant,
                                        fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        Q <- (Q + 1)
                        P <- (P + 1)
                        next
                      }
                  }
                  if (p > 0 && newmodel(p - 1, d, q, P, D, Q, constant,
                                        results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p - 1, d, q), seasonal = c(P,
                                                  D, Q), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p - 1, d, q, P, D, Q, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        p <- (p - 1)
                        next
                      }
                  }
                  if (q > 0 && newmodel(p, d, q - 1, P, D, Q, constant,
                                        results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q - 1), seasonal = c(P,
                                                  D, Q), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p, d, q - 1, P, D, Q, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        q <- (q - 1)
                        next
                      }
                  }
                  if (p < max.p && newmodel(p + 1, d, q, P, D, Q, constant,
                                            results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p + 1, d, q), seasonal = c(P,
                                                  D, Q), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p + 1, d, q, P, D, Q, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        p <- (p + 1)
                        next
                      }
                  }
                  if (q < max.q && newmodel(p, d, q + 1, P, D, Q, constant,
                                            results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p, d, q + 1), seasonal = c(P,
                                                  D, Q), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p, d, q + 1, P, D, Q, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        q <- (q + 1)
                        next
                      }
                  }
                  if (q > 0 && p > 0 && newmodel(p - 1, d, q - 1, P, D,
                                                 Q, constant, results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p - 1, d, q - 1), seasonal = c(P,
                                                  D, Q), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p - 1, d, q - 1, P, D, Q, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        q <- (q - 1)
                        p <- (p - 1)
                        next
                      }
                  }
                  if (q < max.q && p > 0 && newmodel(p - 1, d, q + 1, P,
                                                     D, Q, constant, results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p - 1, d, q + 1), seasonal = c(P,
                                                  D, Q), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p - 1, d, q + 1, P, D, Q, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        q <- (q + 1)
                        p <- (p - 1)
                        next
                      }
                  }
                  if (q > 0 && p < max.p && newmodel(p + 1, d, q - 1, P,
                                                     D, Q, constant, results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p + 1, d, q - 1), seasonal = c(P,
                                                  D, Q), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p + 1, d, q - 1, P, D, Q, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        q <- (q - 1)
                        p <- (p + 1)
                        next
                      }
                  }
                  if (q < max.q && p < max.p && newmodel(p + 1, d, q +
                      1, P, D, Q, constant, results[1:k, ])) {
                    k <- k + 1
                    if (k > nmodels)
                      next
                      fit <- myarima(x, order = c(p + 1, d, q + 1), seasonal = c(P,
                                                  D, Q), constant = constant, ic, trace, approximation,
                                                  method = method, offset = offset, xreg = xreg,
                                                  ...)
                                                  results[k, ] <- c(p + 1, d, q + 1, P, D, Q, constant,
                                                                    fit$ic)
                      if (fit$ic < bestfit$ic) {
                        bestfit <- fit
                        q <- (q + 1)
                        p <- (p + 1)
                        next
                      }
                  }
                  if (allowdrift || allowmean) {
                    if (newmodel(p, d, q, P, D, Q, !constant, results[1:k,
                    ])) {
                      k <- k + 1
                      if (k > nmodels)
                        next
                        fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                                                    D, Q), constant = !constant, ic, trace, approximation,
                                                    method = method, offset = offset, xreg = xreg,
                                                    ...)
                                                    results[k, ] <- c(p, d, q, P, D, Q, !constant,
                                                                      fit$ic)
                        if (fit$ic < bestfit$ic) {
                          bestfit <- fit
                          constant <- !constant
                        }
                    }
                  }
                }
                if (approximation && !is.null(bestfit$arma)) {

                  icorder <- order(results[, 8])
                    nmodels <- sum(!is.na(results[, 8]))
                    for (i in seq(nmodels)) {
                      k <- icorder[i]
                      fit <- myarima(x, order = c(results[k, 1], d, results[k,
                                                  3]), seasonal = c(results[k, 4], D, results[k,
                                                  6]), constant = results[k, 7] == 1, ic, trace,
                                                  approximation = FALSE, method = method, xreg = xreg,
                                                  ...)
                                                  if (fit$ic < Inf) {
                                                    bestfit <- fit
                                                    break
                                                  }
                    }
                }
}


#endif
