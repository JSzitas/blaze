#ifndef STATIONARITY_TESTS
#define STATIONARITY_TESTS

#include "utils/utils.h"

enum kpss_type {
  mu,
  tau
};

enum lag_type {
  short_lag,
  long_lag,
  nil_lag
};

template <typename scalar_t> std::vector<scalar_t> trend_residuals(
    const std::vector<scalar_t>& x) {
  const scalar_t n = x.size();
  scalar_t scale = 0;
  for( size_t i = 0; i < n; i++) {
    scale += i * i;
  }
  scale = 1/scale;
  // compute coefficient
  scalar_t coef = 0;
  for(size_t i = 0; i < n; i++) {
    coef += scale * i * x[i];
  }
  // compute residuals
  std::vector<scalar_t> residuals(n);
  for(size_t i = 0; i < n; i++) {
    residuals[i] = x[i] - (coef * i);
  }
  return residuals;
}

template <typename scalar_t> std::vector<scalar_t> mean_residuals(
    const std::vector<scalar_t>& x) {
  const size_t n = x.size();
  auto mean_val = mean(x);
  // compute residuals
  std::vector<scalar_t> residuals(n);
  for(size_t i = 0; i < n; i++) {
    residuals[i] = x[i] - mean_val;
  }
  return residuals;
}

template <typename scalar_t> std::vector<scalar_t> cumsum(
    const std::vector<scalar_t> & x) {
  std::vector<scalar_t> result(x.size());
  scalar_t sum = 0;
  for(size_t i = 0; i < x.size(); i++) {
    sum += x[i];
    result[i] = sum;
  }
  return result;
}

template <typename scalar_t> scalar_t mse(const std::vector<scalar_t> &x) {
  scalar_t result = 0;
  size_t n = x.size();
  for( auto &val:x ) {
    if(!std::isnan(val)) {
      result += val*val;
    }
    else {
      n--;
    }
  }
  return result/n;
}

template <typename scalar_t> size_t alpha_to_level(const scalar_t alpha) {
  if(alpha >= 0.1) return 0;
  else if(alpha >= 0.05) return 1;
  else if(alpha >= 0.025) return 2;
  else return 3;
}

template <typename scalar_t> bool kpss_test(
    const std::vector<scalar_t> &y,
    const scalar_t alpha = 0.05,
    const kpss_type type = kpss_type::mu,
    const lag_type lags = lag_type::short_lag,
    const bool use_lag = true)
{
  const size_t n = size_omit_nan(y);
  size_t lag_max = 0;
  if(use_lag) {
    lag_max = trunc( 4 * std::pow(n/100, 0.25));
  }
  else if (lags == lag_type::short_lag) {
    lag_max = trunc(4 * std::pow((n/100), 0.25));
  }
  else if (lags == lag_type::long_lag) {
    lag_max = trunc(12 * std::pow((n/100), 0.25));
  }
  else {
    lag_max = 0;
  }
  std::vector<scalar_t> residuals(n);
  scalar_t crit_value = 0;
  if (type == kpss_type::mu) {
    residuals = mean_residuals(y);
    constexpr std::array<double, 4> crit_vals = {0.347, 0.463, 0.574, 0.739};
    crit_value = crit_vals[alpha_to_level(alpha)];
  }
  else if (type == kpss_type::tau) {
    residuals = trend_residuals(y);
    constexpr std::array<double, 4> crit_vals = {0.119, 0.146, 0.176, 0.216};
    crit_value = crit_vals[alpha_to_level(alpha)];
  }
  auto S = cumsum(residuals);
  auto nominator = mse(S)/n;
  auto s2 = mse(residuals);
  scalar_t denominator = s2;
  if(lag_max != 0) {
    std::vector<scalar_t> covs(lag_max);
    for( size_t i = 0; i < lag_max; i++ ) {
      for(size_t j = i+1; j < n; j++) {
        covs[i] += residuals[j] * residuals[j-i-1];
      }
    }
    for(size_t i = 0; i < lag_max; i++) {
      denominator += 2/n * (1 - (i + 1)/(lag_max+1)) * covs[i];
    }
  }
  // alternative hypothesis is stationarity - thus
  // if we are over the critical value we have rejected stationarity
  return (nominator/denominator) > crit_value;
}

enum adf_test_type {
  none,
  drift,
  trend
};

enum select_lag {
  fixed,
  aic,
  bic
};

// template <typename scalar_t> bool augmented_dickey_fuller_test(
//   const std::vector<scalar_t> & y,
//   const adf_test_type type = adf_test_type::none,
//   const size_t lags = 1,
//   const select_lag lag_select = select_lag::fixed) {
//
//     // selectlags <- match.arg(selectlags)
//     // type <- match.arg(type)
//     // if (ncol(as.matrix(y)) > 1)
//       // stop("\ny is not a vector or univariate time series.\n")
//       // if (any(is.na(y)))
//         // stop("\nNAs in y.\n")
//         // y <- as.vector(y)
//         // lag <- as.integer(lags)
//         // if (lag < 0)
//           // stop("\nLags must be set to an non negative integer value.\n")
//           // CALL <- match.call()
//           // DNAME <- deparse(substitute(y))
//           // x.name <- deparse(substitute(y))
//           // lags <- lags + 1
//
//     // lags += 1;
//     // auto lag = lags;
//
//   auto z = diff(y, 1, 1);
//   const size_t n = z.size();
//
//
//         // z <- diff(y)
//           // n <- length(z)
//           x <- embed(z, lags)
//           z.diff <- x[, 1]
//         z.lag.1 <- y[lags:n]
//         tt <- lags:n
//           if (lags > 1) {
//           if (selectlags != "Fixed") {
//             critRes <- rep(NA, lags)
//             for (i in 2:(lags)) {
//               z.diff.lag = x[, 2:i]
//               if (type == "none")
//                 result <- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
//                 if (type == "drift")
//                   result <- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
//                   if (type == "trend")
//                     result <- lm(z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
//                     critRes[i] <- AIC(result, k = switch(selectlags,
//                                                          AIC = 2, BIC = log(length(z.diff))))
//             }
//             lags <- which.min(critRes)
//           }
//           z.diff.lag = x[, 2:lags]
//           if (type == "none") {
//             result <- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
//             tau <- coef(summary(result))[1, 3]
//             teststat <- as.matrix(tau)
//             colnames(teststat) <- "tau1"
//           }
//           if (type == "drift") {
//             result <- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
//             tau <- coef(summary(result))[2, 3]
//             phi1.reg <- lm(z.diff ~ -1 + z.diff.lag)
//             phi1 <- anova(phi1.reg, result)$F[2]
//             teststat <- as.matrix(t(c(tau, phi1)))
//             colnames(teststat) <- c("tau2", "phi1")
//           }
//           if (type == "trend") {
//             result <- lm(z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
//             tau <- coef(summary(result))[2, 3]
//             phi2.reg <- lm(z.diff ~ -1 + z.diff.lag)
//             phi3.reg <- lm(z.diff ~ z.diff.lag)
//             phi2 <- anova(phi2.reg, result)$F[2]
//             phi3 <- anova(phi3.reg, result)$F[2]
//             teststat <- as.matrix(t(c(tau, phi2, phi3)))
//             colnames(teststat) <- c("tau3", "phi2", "phi3")
//           }
//         }
//           else {
//             if (type == "none") {
//               result <- lm(z.diff ~ z.lag.1 - 1)
//               tau <- coef(summary(result))[1, 3]
//               teststat <- as.matrix(tau)
//               colnames(teststat) <- "tau1"
//             }
//             if (type == "drift") {
//               result <- lm(z.diff ~ z.lag.1 + 1)
//               phi1.reg <- lm(z.diff ~ -1)
//               phi1 <- anova(phi1.reg, result)$F[2]
//               tau <- coef(summary(result))[2, 3]
//               teststat <- as.matrix(t(c(tau, phi1)))
//               colnames(teststat) <- c("tau2", "phi1")
//             }
//             if (type == "trend") {
//               result <- lm(z.diff ~ z.lag.1 + 1 + tt)
//               phi2.reg <- lm(z.diff ~ -1)
//               phi3.reg <- lm(z.diff ~ 1)
//               phi2 <- anova(phi2.reg, result)$F[2]
//               phi3 <- anova(phi3.reg, result)$F[2]
//               tau <- coef(summary(result))[2, 3]
//               teststat <- as.matrix(t(c(tau, phi2, phi3)))
//               colnames(teststat) <- c("tau3", "phi2", "phi3")
//             }
//           }
//           rownames(teststat) <- "statistic"
//           testreg <- summary(result)
//             res <- residuals(testreg)
//             if (n < 25)
//               rowselec <- 1
//             if (25 <= n & n < 50)
//               rowselec <- 2
//             if (50 <= n & n < 100)
//               rowselec <- 3
//             if (100 <= n & n < 250)
//               rowselec <- 4
//             if (250 <= n & n < 500)
//               rowselec <- 5
//             if (n >= 500)
//               rowselec <- 6
//             if (type == "none") {
//               cval.tau1 <- rbind(c(-2.66, -1.95, -1.6), c(-2.62, -1.95,
//                                    -1.61), c(-2.6, -1.95, -1.61), c(-2.58, -1.95, -1.62),
//                                    c(-2.58, -1.95, -1.62), c(-2.58, -1.95, -1.62))
//               cvals <- t(cval.tau1[rowselec, ])
//               testnames <- "tau1"
//             }
//             if (type == "drift") {
//               cval.tau2 <- rbind(c(-3.75, -3, -2.63), c(-3.58, -2.93,
//                                    -2.6), c(-3.51, -2.89, -2.58), c(-3.46, -2.88, -2.57),
//                                    c(-3.44, -2.87, -2.57), c(-3.43, -2.86, -2.57))
//               cval.phi1 <- rbind(c(7.88, 5.18, 4.12), c(7.06, 4.86,
//                                    3.94), c(6.7, 4.71, 3.86), c(6.52, 4.63, 3.81), c(6.47,
//                                    4.61, 3.79), c(6.43, 4.59, 3.78))
//               cvals <- rbind(cval.tau2[rowselec, ], cval.phi1[rowselec,
//               ])
//               testnames <- c("tau2", "phi1")
//             }
//             if (type == "trend") {
//               cval.tau3 <- rbind(c(-4.38, -3.6, -3.24), c(-4.15, -3.5,
//                                    -3.18), c(-4.04, -3.45, -3.15), c(-3.99, -3.43, -3.13),
//                                    c(-3.98, -3.42, -3.13), c(-3.96, -3.41, -3.12))
//               cval.phi2 <- rbind(c(8.21, 5.68, 4.67), c(7.02, 5.13,
//                                    4.31), c(6.5, 4.88, 4.16), c(6.22, 4.75, 4.07), c(6.15,
//                                    4.71, 4.05), c(6.09, 4.68, 4.03))
//               cval.phi3 <- rbind(c(10.61, 7.24, 5.91), c(9.31, 6.73,
//                                    5.61), c(8.73, 6.49, 5.47), c(8.43, 6.49, 5.47),
//                                    c(8.34, 6.3, 5.36), c(8.27, 6.25, 5.34))
//               cvals <- rbind(cval.tau3[rowselec, ], cval.phi2[rowselec,
//               ], cval.phi3[rowselec, ])
//               testnames <- c("tau3", "phi2", "phi3")
//             }
//             colnames(cvals) <- c("1pct", "5pct", "10pct")
//               rownames(cvals) <- testnames
//               new("ur.df", y = y, model = type, cval = cvals, lags = lag,
//                   teststat = teststat, testreg = testreg, res = res, test.name = "Augmented Dickey-Fuller Test")
//   }

// ur.pp
//   function (x, type = c("Z-alpha", "Z-tau"), model = c("constant",
//                         "trend"), lags = c("short", "long"), use.lag = NULL)
//   {
//     x <- na.omit(as.vector(x))
//     n <- length(x)
//     y <- x[-1]
//     y.l1 <- x[-n]
//     n <- n - 1
//     lags <- match.arg(lags)
//     model <- match.arg(model)
//     type <- match.arg(type)
//     if (!(is.null(use.lag))) {
//       lmax <- as.integer(use.lag)
//       if (lmax < 0) {
//         warning("\nuse.lag has to be positive and integer; lags='short' used.")
//         lmax <- trunc(4 * (n/100)^0.25)
//       }
//     }
//     else if (lags == "short") {
//       lmax <- trunc(4 * (n/100)^0.25)
//     }
//     else if (lags == "long") {
//       lmax <- trunc(12 * (n/100)^0.25)
//     }
//     if (model == "trend") {
//       cval <- as.matrix(t(c(-3.9638 - 8.353/n - 47.44/(n^2),
//                             -3.4126 - 4.039/n - 17.83/(n^2), -3.1279 - 2.418/n -
//                               7.58/(n^2))))
//       colnames(cval) <- c("1pct", "5pct", "10pct")
//       rownames(cval) <- "critical values"
//       model <- "with intercept and trend"
//       trend <- (1:n) - n/2
//       test.reg <- summary(lm(y ~ y.l1 + trend))
//       res <- residuals(test.reg)
//       my.tstat <- coef(test.reg)[1, 3]
//       beta.tstat <- coef(test.reg)[3, 3]
//       res <- residuals(test.reg)
//       s <- 1/n * (sum(res^2))
//       myybar <- (1/n^2) * sum((y - mean(y))^2)
//       myy <- (1/n^2) * sum(y^2)
//       mty <- (n^(-5/2)) * (t(1:n) %*% y)
//       my <- (n^(-3/2)) * sum(y)
//       idx <- 1:lmax
//       coprods <- sapply(idx, function(l) t(res[-c(1:l)]) %*%
//         (res[-c((n - l + 1):n)]))
//       weights <- 1 - idx/(lmax + 1)
//       sig <- s + (2/n) * (t(weights) %*% coprods)
//       lambda <- 0.5 * (sig - s)
//       lambda.prime <- lambda/sig
//       M <- (1 - n^(-2)) * myy - 12 * mty^2 + 12 * (1 + 1/n) *
//         mty * my - (4 + 6/n + 2/n^2) * my^2
//       my.stat <- sqrt(s/sig) * my.tstat - lambda.prime * sqrt(sig) *
//         my/(sqrt(M) * sqrt((M + my^2)))
//       beta.stat <- sqrt(s/sig) * beta.tstat - lambda.prime *
//         sqrt(sig) * (0.5 * my - mty)/(sqrt(M/12) * sqrt(myybar))
//       aux.stat <- as.matrix(c(round(my.stat, 4), round(beta.stat,
//                                     4)))
//       rownames(aux.stat) <- c("Z-tau-mu", "Z-tau-beta")
//       colnames(aux.stat) <- "aux. Z statistics"
//       if (type == "Z-tau") {
//         tstat <- (coef(test.reg)[2, 1] - 1)/coef(test.reg)[2,
//                                  2]
//         teststat <- sqrt(s/sig) * tstat - lambda.prime *
//           sqrt(sig)/sqrt(M)
//       }
//       else if (type == "Z-alpha") {
//         alpha <- coef(test.reg)[2, 1]
//         teststat <- n * (alpha - 1) - lambda/M
//         cval <- as.matrix(t(c(NA, NA, NA)))
//       }
//     }
//     else if (model == "constant") {
//       cval <- as.matrix(t(c(-3.4335 - 5.999/n - 29.25/(n^2),
//                             -2.8621 - 2.738/n - 8.36/(n^2), -2.5671 - 1.438/n -
//                               4.48/(n^2))))
//       colnames(cval) <- c("1pct", "5pct", "10pct")
//       rownames(cval) <- "critical values"
//       model <- "with intercept"
//       test.reg <- summary(lm(y ~ y.l1))
//       my.tstat <- coef(test.reg)[1, 3]
//       res <- residuals(test.reg)
//       s <- 1/n * (sum(res^2))
//       myybar <- (1/n^2) * sum((y - mean(y))^2)
//       myy <- (1/n^2) * sum(y^2)
//       my <- (n^(-3/2)) * sum(y)
//       idx <- 1:lmax
//       coprods <- sapply(idx, function(l) t(res[-c(1:l)]) %*%
//         (res[-c((n - l + 1):n)]))
//       weights <- 1 - idx/(lmax + 1)
//       sig <- s + (2/n) * (t(weights) %*% coprods)
//       lambda <- 0.5 * (sig - s)
//       lambda.prime <- lambda/sig
//       my.stat <- sqrt(s/sig) * my.tstat + lambda.prime * sqrt(sig) *
//         my/(sqrt(myy) * sqrt(myybar))
//       aux.stat <- as.matrix(round(my.stat, 4))
//       rownames(aux.stat) <- "Z-tau-mu"
//       colnames(aux.stat) <- "aux. Z statistics"
//       if (type == "Z-tau") {
//         tstat <- (coef(test.reg)[2, 1] - 1)/coef(test.reg)[2,
//                                  2]
//         teststat <- sqrt(s/sig) * tstat - lambda.prime *
//           sqrt(sig)/sqrt(myybar)
//       }
//       else if (type == "Z-alpha") {
//         alpha <- coef(test.reg)[2, 1]
//         teststat <- n * (alpha - 1) - lambda/myybar
//         cval <- as.matrix(t(c(NA, NA, NA)))
//       }
//     }
//     new("ur.pp", y = y, type = type, model = model, lag = as.integer(lmax),
//         cval = cval, teststat = as.numeric(teststat), testreg = test.reg,
//         auxstat = aux.stat, res = res, test.name = "Phillips-Perron")
//   }


#endif
