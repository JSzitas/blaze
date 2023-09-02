#ifndef BLAZE_HOLT_WINTERS
#define BLAZE_HOLT_WINTERS

#include <stdexcept>

#include "typedefs/seasonality.h"

#include "utils/seasonality_utils.h"
#include "utils/boxcox.h"

template <typename scalar_t,
          const SeasonalityType seas_type = SeasonalityType::A,
          const bool exponential = false> class HoltWinters {
public:
  HoltWinters<scalar_t, seas_type, exponential>(){};
  HoltWinters<scalar_t, seas_type, exponential>(
    const std::vector<scalar_t> &y,
    size_t seasonality = 0, 
    scalar_t alpha = 0.5, 
    scalar_t beta = std::nan("0.0"), 
    scalar_t gamma = std::nan("0.0"), 
    scalar_t phi = std::nan("0.0"), 
    scalar_t lambda = std::nan("0.0"), 
    const bool use_box_cox = false,
    const bool bias_adjust = false) {
    
    if(use_box_cox) {
        this->transformer = BoxCoxTransformer<scalar_t>(y, lambda);
    }
    if(std::isnan(phi))  this->phi = 1.0;
    if(alpha < 0 || alpha > 1) {
      std::cout << "Invalid value of alpha supplied - initializing to 0.5." << std::endl;
      this->alpha = 0.5;
    }
    if(beta < 0 || beta > 1) {
      std::cout << "Invalid value of beta supplied - initializing to 0.3." << std::endl;
      this->beta = 0.3;
    }
    if(gamma < 0 || gamma > 1) {
      std::cout << "Invalid value of gamma supplied - initializing to 0.3." << std::endl;
      this->beta = 0.3;
    }
    if(std::isnan(gamma) || gamma < 0.00001) {
      if(seas_type == SeasonalityType::M) {
        for( auto & val:y) {
          if(val <= 0) {
            throw std::invalid_argument("All data must be positive if using multiplicative seasonality in HoltWinters.");
          }; 
        }
      }
    }
    if (m <= 1) {
      gamma <- FALSE
    }
    if (!is.null(gamma) && is.logical(gamma) && !gamma) {
      seasonal <- "none"
      l.start <- x[1L]
      s.start <- 0
      if (is.null(beta) || !is.logical(beta) || beta) {
        if (!exponential) {
          b.start <- x[2L] - x[1L]
        }
        else {
          b.start <- x[2L]/x[1L]
        }
      }
    }
    else {
      l.start <- mean(x[1:m])
      b.start <- (mean(x[m + (1:m)]) - l.start)/m
      if (seasonal == "additive") {
        s.start <- x[1:m] - l.start
      }
      else {
        s.start <- x[1:m]/l.start
      }
    }
    lower <- c(0, 0, 0, 0)
      upper <- c(1, 1, 1, 1)
      if (!is.null(beta) && is.logical(beta) && !beta) {
        trendtype <- "N"
      }
      else if (exponential) {
        trendtype <- "M"
      }
      else {
        trendtype <- "A"
      }
      if (seasonal == "none") {
        seasontype <- "N"
      }
      else if (seasonal == "multiplicative") {
        seasontype <- "M"
      }
      else {
        seasontype <- "A"
      }
      optim.start <- initparam(alpha = alpha, beta = beta, gamma = gamma, 
                               phi = 1, trendtype = trendtype, seasontype = seasontype, 
                               damped = FALSE, lower = lower, upper = upper, m = m)
        error <- function(p, select) {
          if (select[1] > 0) {
            alpha <- p[1L]
          }
          if (select[2] > 0) {
            beta <- p[1L + select[1]]
          }
          if (select[3] > 0) {
            gamma <- p[1L + select[1] + select[2]]
          }
          zzhw(x, lenx = lenx, alpha = alpha, beta = beta, gamma = gamma, 
               seasonal = seasonal, m = m, dotrend = (!is.logical(beta) || 
                 beta), doseasonal = (!is.logical(gamma) || gamma), 
                   exponential = exponential, phi = phi, l.start = l.start, 
                   b.start = b.start, s.start = s.start)$SSE
        }
      select <- as.numeric(c(is.null(alpha), is.null(beta), is.null(gamma)))
        if (sum(select) > 0) {
          sol <- optim(optim.start, error, method = "L-BFGS-B", 
                       lower = lower[select], upper = upper[select], select = select)
          if (sol$convergence || any(sol$par < 0 | sol$par > 1)) {
            if (sol$convergence > 50) {
              if (warnings) {
                warning(gettextf("optimization difficulties: %s", 
                                 sol$message), domain = NA)
              }
            }
            else {
              stop("optimization failure")
            }
          }
          if (select[1] > 0) {
            alpha <- sol$par[1L]
          }
          if (select[2] > 0) {
            beta <- sol$par[1L + select[1]]
          }
          if (select[3] > 0) {
            gamma <- sol$par[1L + select[1] + select[2]]
          }
        }
        final.fit <- zzhw(x, lenx = lenx, alpha = alpha, beta = beta, 
                          gamma = gamma, seasonal = seasonal, m = m, dotrend = (!is.logical(beta) || 
                            beta), doseasonal = (!is.logical(gamma) || gamma), 
                              exponential = exponential, phi = phi, l.start = l.start, 
                              b.start = b.start, s.start = s.start)
          tspx <- tsp(x)
          fitted <- ts(final.fit$fitted, frequency = m, start = tspx[1])
          res <- ts(final.fit$residuals, frequency = m, start = tspx[1])
          if (!is.null(lambda)) {
            fitted <- InvBoxCox(fitted, lambda, biasadj, var(final.fit$residuals))
            attr(lambda, "biasadj") <- biasadj
          }
          states <- matrix(final.fit$level, ncol = 1)
            colnames(states) <- "l"
          if (trendtype != "N") {
            states <- cbind(states, b = final.fit$trend)
          }
          if (seasontype != "N") {
            nr <- nrow(states)
            nc <- ncol(states)
            for (i in 1:m) states <- cbind(states, final.fit$season[(m - 
                 i) + (1:nr)])
              colnames(states)[nc + (1:m)] <- paste("s", 1:m, sep = "")
          }
          states <- ts(states, frequency = m, start = tspx[1] - 1/m)
            damped <- (phi < 1)
            // if (seasonal == "additive") {
            //   components <- c("A", trendtype, seasontype, damped)
            // }
            // else if (seasonal == "multiplicative") {
            //   components <- c("M", trendtype, seasontype, damped)
            // }
            // else if (seasonal == "none" && exponential) {
            //   components <- c("M", trendtype, seasontype, damped)
            // }
            // else {
            //   components <- c("A", trendtype, seasontype, damped)
            // }
            // initstate <- states[1, ]
            // param <- alpha
            //   names(param) <- "alpha"
            // if (trendtype != "N") {
            //   param <- c(param, beta = beta)
            //   names(param)[length(param)] <- "beta"
            // }
            // if (seasontype != "N") {
            //   param <- c(param, gamma = gamma)
            //   names(param)[length(param)] <- "gamma"
            // }
            // if (damped) {
            //   param <- c(param, phi = phi)
            //   names(param)[length(param)] <- "phi"
            // }
            if constexpr(seas_type == SeasonalityType::A) {
              this->sigma2 = mean(this->residuals);
            }
            else if constexpr(seas_type == SeasonalityType::M) {
              // otherwise a bit more involved 
              scalar_t temp = 0.0;
              for(size_t i = 0; i < this->y.size(); i++) {
                temp += std::pow(this->residuals[i]/this->fitted[i], 2)
              }
              temp /= this->y.size();
              sigma2 = temp;
            }
            // structure(list(fitted = fitted, residuals = res, components = components, 
            //                x = x, par = c(param, initstate), initstate = initstate, 
            //                states = states, SSE = final.fit$SSE, sigma2 = sigma2, 
            //                call = match.call(), m = m), class = "ets")
  };
  void fit(){
    
  }
private:
  template <const bool trend = false,
            const bool seasonal = false, 
            const bool update_all=false> scalar_t loss(
              const std::vector<scalar_t> &y,
                const size_t m, 
                std::vector<scalar_t> &s_start, 
                std::vector<scalar_t> &level_temp, 
                std::vector<scalar_t> &trend_temp, 
                std::vector<scalar_t> &season_temp, 
                std::vector<scalar_t> &residuals, 
                std::vector<scalar_t> &fitted, 
                scalar_t alpha = 0.0,
                scalar_t beta = 0.0,
                scalar_t gamma = 0.0, 
                scalar_t l_start = 0, // initial level? 
                scalar_t phi = 1.0,
                scalar_t b_start = 0) {
                if constexpr(!trend) {
                  beta = 0;
                  b_start = 0.0;
                }
                if constexpr(!seasonal) {
                  gamma = 0.0;
                  for(auto & val:s_start) {
                    if constexpr(seas_type == SeasonalityType::A) {
                      val = 0.0;
                    } else {
                      val = 1.0;
                    }
                  }
                }
                scalar_t last_level = l_start, last_season = s_start[0],
                                                                    last_trend = b_start, y_hat = 0.0,
                                                                    current_level = 0.0, current_trend = 0.0,
                                                                    result = 0.0;
                
                for(size_t i = 0; i < y.size(); i++) {
                  if constexpr(seasonal) {
                    if (i > m) {
                      last_season = season_temp[i - m];
                    }
                    else {
                      last_season = s_start[i];
                    }
                  } else {
                    last_season = seas_type != SeasonalityType::A;
                  }
                  if constexpr(seas_type == SeasonalityType::A) {
                    if constexpr(!exponential) {
                      y_hat = last_level + phi * last_trend + last_season;
                      current_level = alpha * (y[i] - last_season) + (1 - 
                        alpha) * (last_level + phi * last_trend);
                      current_trend = beta * (current_level - last_level) + (1 - 
                        beta) * phi * last_trend;
                      season_temp[i] = gamma * (y[i] - last_level - phi * 
                        last_trend) + (1 - gamma) * last_season;
                      
                    }
                    else {
                      y_hat = last_level * last_trend^phi + last_season;
                      current_level = alpha * (y[i] - last_season) + (1 - 
                        alpha) * (last_level * last_trend^phi);
                      current_trend = beta * (current_level/last_level) + (1 - beta) * 
                        last_trend^phi;
                      season_temp[i] = gamma * (y[i] - last_level * last_trend^phi) + 
                        (1 - gamma) * last_season;
                    }
                  }
                  // multiplicative seasonality
                  else if constexpr(seas_type == SeasonalityType::M){
                    if constexpr(!exponential) {
                      y_hat = (last_level + phi * last_trend) * last_season;
                      current_level = alpha * (y[i]/last_season) + (1 - 
                        alpha) * (last_level + phi * last_trend);
                      current_trend = beta * (current_level - last_level) + (1 - 
                        beta) * phi * last_trend;
                      season_temp[i] = gamma * (y[i]/(last_level + phi * 
                        last_trend)) + (1 - gamma) * last_season;
                    }
                    else {
                      y_hat = last_level * last_trend^phi * last_season; 
                      current_level = alpha * (y[i]/last_season) + (1 - 
                        alpha) * (last_level * last_trend^phi);
                      current_trend = beta * (current_level/last_level) + (1 - beta) * 
                        last_trend^phi;
                      season_temp[i] =- gamma * (y[i]/(last_level * last_trend^phi)) + 
                        (1 - gamma) * last_season;
                    }
                  }
                  last_level = current_level;
                  last_trend = current_trend;
                  // update residuals, fitted values, sum of squares and temporaries 
                  // as needed
                  scalar_t residual = y[i] - y_hat[i];
                  if constexpr(update_all) {
                    fitted[i] = y_hat;
                    residuals[i] = y[i] - y_hat[i];
                    trend_temp[i] = current_trend;
                    level_temp[i] = current_level;
                  }
                  result += std::pow(residual, 2);
                }
                return result;
              }
  
  
              
};


// HoltWintersZZ = function (x, alpha = NULL, beta = NULL, gamma = NULL, seasonal = c("additive", 
//                                                                      "multiplicative"), 
//                                                                      exponential = FALSE, phi = NULL, lambda = NULL, 
//                                                                      biasadj = FALSE, warnings = TRUE) 
//   {
//     // x <- as.ts(x)
//     // seasonal <- match.arg(seasonal)
//     // m <- frequency(x)
//     // lenx <- length(x)
//     if (!is.null(lambda)) {
//       x <- BoxCox(x, lambda)
//       lambda <- attr(x, "lambda")
//     }
//     if (is.null(phi) || !is.numeric(phi)) {
//       phi <- 1
//     }
//     if (!is.null(alpha) && !is.numeric(alpha)) {
//       stop("cannot fit models without level ('alpha' must not be 0 or FALSE).")
//     }
//     if (!all(is.null(c(alpha, beta, gamma))) && any(c(alpha, 
//                      beta, gamma) < 0 | c(alpha, beta, gamma) > 1)) {
//       stop("'alpha', 'beta' and 'gamma' must be within the unit interval.")
//     }
//     if ((is.null(gamma) || gamma > 0)) {
//       if (seasonal == "multiplicative" && any(x <= 0)) {
//         stop("data must be positive for multiplicative Holt-Winters.")
//       }
//     }
//     if (m <= 1) {
//       gamma <- FALSE
//     }
//     if (!is.null(gamma) && is.logical(gamma) && !gamma) {
//       seasonal <- "none"
//       l.start <- x[1L]
//       s.start <- 0
//       if (is.null(beta) || !is.logical(beta) || beta) {
//         if (!exponential) {
//           b.start <- x[2L] - x[1L]
//         }
//         else {
//           b.start <- x[2L]/x[1L]
//         }
//       }
//     }
//     else {
//       l.start <- mean(x[1:m])
//       b.start <- (mean(x[m + (1:m)]) - l.start)/m
//       if (seasonal == "additive") {
//         s.start <- x[1:m] - l.start
//       }
//       else {
//         s.start <- x[1:m]/l.start
//       }
//     }
//     lower <- c(0, 0, 0, 0)
//       upper <- c(1, 1, 1, 1)
//       if (!is.null(beta) && is.logical(beta) && !beta) {
//         trendtype <- "N"
//       }
//       else if (exponential) {
//         trendtype <- "M"
//       }
//       else {
//         trendtype <- "A"
//       }
//       if (seasonal == "none") {
//         seasontype <- "N"
//       }
//       else if (seasonal == "multiplicative") {
//         seasontype <- "M"
//       }
//       else {
//         seasontype <- "A"
//       }
//       optim.start <- initparam(alpha = alpha, beta = beta, gamma = gamma, 
//                                phi = 1, trendtype = trendtype, seasontype = seasontype, 
//                                damped = FALSE, lower = lower, upper = upper, m = m)
//         error <- function(p, select) {
//           if (select[1] > 0) {
//             alpha <- p[1L]
//           }
//           if (select[2] > 0) {
//             beta <- p[1L + select[1]]
//           }
//           if (select[3] > 0) {
//             gamma <- p[1L + select[1] + select[2]]
//           }
//           zzhw(x, lenx = lenx, alpha = alpha, beta = beta, gamma = gamma, 
//                seasonal = seasonal, m = m, dotrend = (!is.logical(beta) || 
//                  beta), doseasonal = (!is.logical(gamma) || gamma), 
//                    exponential = exponential, phi = phi, l.start = l.start, 
//                    b.start = b.start, s.start = s.start)$SSE
//         }
//       select <- as.numeric(c(is.null(alpha), is.null(beta), is.null(gamma)))
//         if (sum(select) > 0) {
//           sol <- optim(optim.start, error, method = "L-BFGS-B", 
//                        lower = lower[select], upper = upper[select], select = select)
//           if (sol$convergence || any(sol$par < 0 | sol$par > 1)) {
//             if (sol$convergence > 50) {
//               if (warnings) {
//                 warning(gettextf("optimization difficulties: %s", 
//                                  sol$message), domain = NA)
//               }
//             }
//             else {
//               stop("optimization failure")
//             }
//           }
//           if (select[1] > 0) {
//             alpha <- sol$par[1L]
//           }
//           if (select[2] > 0) {
//             beta <- sol$par[1L + select[1]]
//           }
//           if (select[3] > 0) {
//             gamma <- sol$par[1L + select[1] + select[2]]
//           }
//         }
//         final.fit <- zzhw(x, lenx = lenx, alpha = alpha, beta = beta, 
//                           gamma = gamma, seasonal = seasonal, m = m, dotrend = (!is.logical(beta) || 
//                             beta), doseasonal = (!is.logical(gamma) || gamma), 
//                               exponential = exponential, phi = phi, l.start = l.start, 
//                               b.start = b.start, s.start = s.start)
//           tspx <- tsp(x)
//           fitted <- ts(final.fit$fitted, frequency = m, start = tspx[1])
//           res <- ts(final.fit$residuals, frequency = m, start = tspx[1])
//           if (!is.null(lambda)) {
//             fitted <- InvBoxCox(fitted, lambda, biasadj, var(final.fit$residuals))
//             attr(lambda, "biasadj") <- biasadj
//           }
//           states <- matrix(final.fit$level, ncol = 1)
//             colnames(states) <- "l"
//           if (trendtype != "N") {
//             states <- cbind(states, b = final.fit$trend)
//           }
//           if (seasontype != "N") {
//             nr <- nrow(states)
//             nc <- ncol(states)
//             for (i in 1:m) states <- cbind(states, final.fit$season[(m - 
//                  i) + (1:nr)])
//               colnames(states)[nc + (1:m)] <- paste("s", 1:m, sep = "")
//           }
//           states <- ts(states, frequency = m, start = tspx[1] - 1/m)
//             damped <- (phi < 1)
//             if (seasonal == "additive") {
//               components <- c("A", trendtype, seasontype, damped)
//             }
//             else if (seasonal == "multiplicative") {
//               components <- c("M", trendtype, seasontype, damped)
//             }
//             else if (seasonal == "none" && exponential) {
//               components <- c("M", trendtype, seasontype, damped)
//             }
//             else {
//               components <- c("A", trendtype, seasontype, damped)
//             }
//             initstate <- states[1, ]
//             param <- alpha
//               names(param) <- "alpha"
//             if (trendtype != "N") {
//               param <- c(param, beta = beta)
//               names(param)[length(param)] <- "beta"
//             }
//             if (seasontype != "N") {
//               param <- c(param, gamma = gamma)
//               names(param)[length(param)] <- "gamma"
//             }
//             if (damped) {
//               param <- c(param, phi = phi)
//               names(param)[length(param)] <- "phi"
//             }
//             if (components[1] == "A") {
//               sigma2 <- mean(res^2)
//             }
//             else {
//               sigma2 <- mean((res/fitted)^2)
//             }
//             structure(list(fitted = fitted, residuals = res, components = components, 
//                            x = x, par = c(param, initstate), initstate = initstate, 
//                            states = states, SSE = final.fit$SSE, sigma2 = sigma2, 
//                            call = match.call(), m = m), class = "ets")
//   }

#endif
