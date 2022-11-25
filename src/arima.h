#ifndef ARIMA_WRAPPER
#define ARIMA_WRAPPER

#include "initializers.h"
#include "utils.h"
#include "structural_model.h"
#include "xreg.h"

// defines the arima structure
struct arima_kind{
   arima_kind( int p, int d, int q, int P, int D, int Q, int s_period ) :
    p(p), d(d), q(q), P(P), D(D), Q(Q), s_period(s_period){}
  int p, d, q, P, D, Q, s_period;
};

template <typename U=double> class Arima {
  Arima<U>(std::vector<U> & y,
           arima_kind kind,
           std::vector<std::vector<U>> xreg = {{}},
           bool intercept = true,
           bool transform_parameters = true,
           SSinit ss_init = Gardner,
           U kappa = 1000000 ){
    this->y = y;
    this->xreg = xreg;
    this->intercept = intercept;
    this->transform_parameters = transform_parameters;
    this->ss_init = ss_init;
    this->kind = kind;
    // this->arma_structure = std::vector<int>{ order[0], order[2], seas_order[0],
    //                                          seas_order[2], seasonal_period,
    //                                          order[1], seas_order[1]};
    // for(int i=0; i < 4; i++) {
    //   this->n_arma_coef += this->arma_structure[i];
    // }

                  // seasonal$period <- frequency(x)
                  // arma <- as.integer(c(order[-2L], seasonal$order[-2L], seasonal$period,
                  //                      order[2L], seasonal$order[2L]))
                  // narma <- sum(arma[1L:4L])


                // Delta <- 1
                // for (i in seq_len(order[2L])) Delta <- Delta %+% c(1, -1)
                //   for (i in seq_len(seasonal$order[2L])) Delta <- Delta %+%
                //     c(1, rep.int(0, seasonal$period - 1), -1)
                //     Delta <- -Delta[-1L]
                //   nd <- order[2L] + seasonal$order[2L]
                //   n.used <- sum(!is.na(x)) - length(Delta)
                //     if (is.null(xreg)) {
                //       ncxreg <- 0L
                //     }
                //     else {
                //       nmxreg <- deparse1(substitute(xreg))
                //       if (NROW(xreg) != n)
                //         stop("lengths of 'x' and 'xreg' do not match")
                //         ncxreg <- NCOL(xreg)
                //         xreg <- as.matrix(xreg)
                //         storage.mode(xreg) <- "double"
                //     }
                //     class(xreg) <- NULL
                //     if (ncxreg > 0L && is.null(colnames(xreg)))
                //       colnames(xreg) <- if (ncxreg == 1L)
                //         nmxreg
                //         else paste0(nmxreg, 1L:ncxreg)
                //           if (include.mean && (nd == 0L)) {
                //             xreg <- cbind(intercept = rep(1, n), xreg = xreg)
                //             ncxreg <- ncxreg + 1L
                //           }
                //           if (method == "CSS-ML") {
                //             anyna <- anyNA(x)
                //             if (ncxreg)
                //               anyna <- anyna || anyNA(xreg)
                //               if (anyna)
                //                 method <- "ML"
                //           }
                //           if (method == "CSS" || method == "CSS-ML") {
                //             ncond <- order[2L] + seasonal$order[2L] * seasonal$period
                //             ncond1 <- order[1L] + seasonal$period * seasonal$order[1L]
                //             ncond <- ncond + if (!missing(n.cond))
                //               max(n.cond, ncond1)
                //               else ncond1
                //           }
                //           else ncond <- 0
                //           if (is.null(fixed))
                //             fixed <- rep(NA_real_, narma + ncxreg)
                //             else if (length(fixed) != narma + ncxreg)
                //               stop("wrong length for 'fixed'")
                //               mask <- is.na(fixed)
                //               no.optim <- !any(mask)
                //               if (no.optim)
                //                 transform.pars <- FALSE
                //                 if (transform.pars) {
                //                   ind <- arma[1L] + arma[2L] + seq_len(arma[3L])
                //                   if (any(!mask[seq_len(arma[1L])]) || any(!mask[ind])) {
                //                     warning("some AR parameters were fixed: setting transform.pars = FALSE")
                //                     transform.pars <- FALSE
                //                   }
                //                 }
                //                 init0 <- rep.int(0, narma)
                //                   parscale <- rep(1, narma)
                //                   if (ncxreg) {
                //                     cn <- colnames(xreg)
                //                     orig.xreg <- (ncxreg == 1L) || any(!mask[narma + 1L:ncxreg])
                //                     if (!orig.xreg) {
                //                       S <- svd(na.omit(xreg))
                //                       xreg <- xreg %*% S$v
                //                     }
                //                     dx <- x
                //                       dxreg <- xreg
                //                       if (order[2L] > 0L) {
                //                         dx <- diff(dx, 1L, order[2L])
                //                         dxreg <- diff(dxreg, 1L, order[2L])
                //                       }
                //                       if (seasonal$period > 1L && seasonal$order[2L] > 0) {
                //                         dx <- diff(dx, seasonal$period, seasonal$order[2L])
                //                         dxreg <- diff(dxreg, seasonal$period, seasonal$order[2L])
                //                       }
                //                       fit <- if (length(dx) > ncol(dxreg))
                //                         lm(dx ~ dxreg - 1, na.action = na.omit)
                //                         else list(rank = 0L)
                //                           if (fit$rank == 0L) {
                //                             fit <- lm(x ~ xreg - 1, na.action = na.omit)
                //                           }
                //                           isna <- is.na(x) | apply(xreg, 1L, anyNA)
                //                             n.used <- sum(!isna) - length(Delta)
                //                             init0 <- c(init0, coef(fit))
                //                             ses <- summary(fit)$coefficients[, 2L]
                //                           parscale <- c(parscale, 10 * ses)
                //                   }


  };
  void fit(){
    // this should just proceed with fitting, not do things which can be done in
    // the constructor


  };
  forecast_result<U> predict(int h = 10, std::vector<std::vector<U>> newxreg = {{}}){
     // validate xreg length
     if( newxreg.size() != this->xreg.size() ) {
       return forecast_result<U>(0) ;
     }
     auto res = kalman_forecast(h, this-> structural_arma_model);
     auto xreg_adjusted = std::vector<U>(h);
     if( this->reg_coef.size() > 0 ) {
       // get the result of xreg regression
       xreg_adjusted = predict(this->reg_coef, newxreg);
     }

     res.forecast = xreg_adjusted + res.forecast;
     for( int i = 0; i < res.se.size(); i++ ) {
       res.se[i] = res.se[i] * sigma2;
     }
    return res;
  };
  std::vector<U> y;
  arima_kind kind;
  std::vector<U> residuals;
  std::vector<std::vector<U>> xreg;
  lm_coef<U> reg_coef;
  bool intercept;
  bool transform_parameters;
  SSinit ss_init;
  U sigma2;
};

#endif
