upARIMA <- function(mod, phi, theta, r, SS.G) {
  p <- length(phi)
  q <- length(theta)
  mod$phi <- phi
  mod$theta <- theta
  r <- max(p, q + 1L)
  if (p > 0)
    mod$T[1L:p, 1L] <- phi
  if (r > 1L)
    mod$Pn[1L:r, 1L:r] <- if (SS.G)
      .Call(stats:::C_getQ0, phi, theta)
  else .Call(stats:::C_getQ0bis, phi, theta, tol = 0)
  else mod$Pn[1L, 1L] <- if (p > 0)
    1/(1 - phi^2)
  else 1
  mod$a[] <- 0
  mod
}
arimaSS <- function(y, mod) {
  ## next call changes Z components a, P, Pn so beware!
  .Call(stats:::C_ARIMA_Like, y, mod, 0L, TRUE)
}
arCheck <- function(ar) {
  p <- max(which(c(1, -ar) != 0)) - 1
  if (!p)
    return(TRUE)
  all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
}
maInvert <- function(ma) {
  q <- length(ma)
  q0 <- max(which(c(1, ma) != 0)) - 1L
  if (!q0)
    return(ma)
  roots <- polyroot(c(1, ma[1L:q0]))
  ind <- Mod(roots) < 1
  if (all(!ind))
    return(ma)
  if (q0 == 1)
    return(c(1/ma[1L], rep.int(0, q - q0)))
  roots[ind] <- 1/roots[ind]
  x <- 1
  for (r in roots) x <- c(x, 0) - c(0, x)/r
  c(Re(x[-1L]), rep.int(0, q - q0))
}

make_delta_r <- function( n_diff, seas_period, n_seas_diff ) {
  Delta <- 1
  for (i in seq_len(n_diff)) {
    Delta <- ts_conv(Delta, c(1, -1))
  }
  for (i in seq_len(n_seas_diff)) {
    Delta <- ts_conv(Delta, c(1, rep.int(0, seas_period - 1), -1))
  }
  return(-Delta[-1L])
}

arimaSS <- function(y, mod) {
  ## next call changes Z components a, P, Pn so beware!
  .Call(stats:::C_ARIMA_Like, y, mod, 0L, TRUE)
}

upARIMA <- function(mod, phi, theta, r, SS.G) {
  p <- length(phi)
  q <- length(theta)
  mod$phi <- phi
  mod$theta <- theta
  r <- max(p, q + 1L)
  if (p > 0)
    mod$T[1L:p, 1L] <- phi
  if (r > 1L)
    mod$Pn[1L:r, 1L:r] <- if (SS.G)
      .Call(stats:::C_getQ0, phi, theta)
  else .Call(stats:::C_getQ0bis, phi, theta, tol = 0)
  else mod$Pn[1L, 1L] <- if (p > 0)
    1/(1 - phi^2)
  else 1
  mod$a[] <- 0
  mod
}


arima3 <- function (x,
                    order = c(0L, 0L, 0L),
                    seasonal = list(order = c(0L,0L, 0L),
                                    period = 1),
                    xreg = NULL,
                    include.mean = TRUE,
                    transform.pars = TRUE,
                    fixed = NULL,
                    init = NULL,
                    method = c("CSS-ML",
                               "ML",
                               "CSS")[1],
                    n.cond,
                    SSinit = c("Gardner1980", "Rossignol2011")[1],
                    optim.method = "BFGS",
                    optim.control = list(),
                    kappa = 1e+06)
{
  "%+%" <- function(a, b) .Call(stats:::C_TSconv, a, b)
  SS.G <- SSinit == "Gardner1980"
  armafn <- function(p, trans) {
    par <- coef
    par[mask] <- p
    trarma <- .Call(stats:::C_ARIMA_transPars, par, arma, trans)
    if (is.null(Z <- tryCatch(upARIMA(mod, trarma[[1L]],
                                      trarma[[2L]], r, SS.G), error = function(e) NULL)))
      return(.Machine$double.xmax)
    if (ncxreg > 0)
      x <- x - xreg %*% par[narma + (1L:ncxreg)]
    res <- .Call(stats:::C_ARIMA_Like, x, Z, 0L, FALSE)
    s2 <- res[1L]/res[3L]
    0.5 * (log(s2) + res[2L]/res[3L])
  }
  armaCSS <- function(p) {# fixed, mask, arma, trarma, ncond, x, xreg, narma, ncxreg) {
    par <- as.double(fixed)
    par[mask] <- p
    trarma <- .Call(stats:::C_ARIMA_transPars, par, arma, FALSE)
    if (ncxreg > 0)
      x <- x - xreg %*% par[narma + (1L:ncxreg)]
    res <- .Call(stats:::C_ARIMA_CSS, x, arma, trarma[[1L]], trarma[[2L]],
                 as.integer(ncond), FALSE)
    0.5 * log(res)
  }
  arma <- as.integer(c(order[-2L], seasonal$order[-2L], seasonal$period,
                       order[2L], seasonal$order[2L]))
  narma <- sum(arma[1L:4L])
  xtsp <- tsp(x)
  tsp(x) <- NULL
  Delta <- make_delta_r( order[2L], seasonal$period, seasonal$order[2])
  nd <- order[2L] + seasonal$order[2L]
  n.used <- sum(!is.na(x)) - length(Delta)
  if (is.null(xreg)) {
    ncxreg <- 0L
  }
  else {
    nmxreg <- deparse1(substitute(xreg))
    if (NROW(xreg) != n)
      stop("lengths of 'x' and 'xreg' do not match")
    ncxreg <- NCOL(xreg)
    xreg <- as.matrix(xreg)
    storage.mode(xreg) <- "double"
  }
  class(xreg) <- NULL
  if (ncxreg > 0L && is.null(colnames(xreg)))
    colnames(xreg) <- if (ncxreg == 1L)
      nmxreg
  else paste0(nmxreg, 1L:ncxreg)
  if (include.mean && (nd == 0L)) {
    xreg <- cbind(intercept = rep(1, n), xreg = xreg)
    ncxreg <- ncxreg + 1L
  }
  if (method == "CSS-ML") {
    anyna <- anyNA(x)
    if (ncxreg)
      anyna <- anyna || anyNA(xreg)
    if (anyna)
      method <- "ML"
  }
  if (method == "CSS" || method == "CSS-ML") {
    ncond <- order[2L] + seasonal$order[2L] * seasonal$period
    ncond1 <- order[1L] + seasonal$period * seasonal$order[1L]
    ncond <- ncond + if (!missing(n.cond))
      max(n.cond, ncond1)
    else ncond1
  }
  else ncond <- 0
  if (transform.pars) {
    ind <- arma[1L] + arma[2L] + seq_len(arma[3L])
  }
  init0 <- rep.int(0, narma)
  parscale <- rep(1, narma)
  if (ncxreg) {
    cn <- colnames(xreg)
    orig.xreg <- (ncxreg == 1L) || any(!mask[narma + 1L:ncxreg])
    if (!orig.xreg) {
      S <- svd(na.omit(xreg))
      xreg <- xreg %*% S$v
    }
    dx <- x
    dxreg <- xreg
    if (order[2L] > 0L) {
      dx <- diff(dx, 1L, order[2L])
      dxreg <- diff(dxreg, 1L, order[2L])
    }
    if (seasonal$period > 1L && seasonal$order[2L] > 0) {
      dx <- diff(dx, seasonal$period, seasonal$order[2L])
      dxreg <- diff(dxreg, seasonal$period, seasonal$order[2L])
    }
    fit <- if (length(dx) > ncol(dxreg))
      lm(dx ~ dxreg - 1, na.action = na.omit)
    else list(rank = 0L)
    if (fit$rank == 0L) {
      fit <- lm(x ~ xreg - 1, na.action = na.omit)
    }
    isna <- is.na(x) | apply(xreg, 1L, anyNA)
    n.used <- sum(!isna) - length(Delta)
    init0 <- c(init0, coef(fit))
    ses <- summary(fit)$coefficients[, 2L]
    parscale <- c(parscale, 10 * ses)
  }
  if (n.used <= 0)
    stop("too few non-missing observations")
  init <- init0
  fixed <- rep(NA_real_, narma + ncxreg)
  coef <- as.double(fixed)
  mask <- is.na(fixed)
  if (!("parscale" %in% names(optim.control)))
    optim.control$parscale <- parscale[mask]
  if (method == "CSS") {
    res <- optim(init[mask], armaCSS, method = optim.method,
                 hessian = TRUE, control = optim.control)
    coef[mask] <- res$par
    trarma <- .Call(stats:::C_ARIMA_transPars, coef, arma, FALSE)
    mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa,
                     SSinit)
    if (ncxreg > 0) {
      x <- x - xreg %*% coef[narma + (1L:ncxreg)]
    }
    val <- .Call(stats:::C_ARIMA_CSS, x, arma, trarma[[1L]], trarma[[2L]],
                 as.integer(ncond), TRUE)
    sigma2 <- val[[1L]]
    var <- solve(res$hessian * n.used)
  }
  else {
    if (method == "CSS-ML") {
      res <- optim(init[mask], armaCSS, method = optim.method,
                 hessian = FALSE, control = optim.control)
      if (res$convergence == 0)
        init[mask] <- res$par
      if (arma[1L] > 0)
        if (!arCheck(init[1L:arma[1L]]))
          stop("non-stationary AR part from CSS")
      if (arma[3L] > 0)
        if (!arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
          stop("non-stationary seasonal AR part from CSS")
      ncond <- 0L
    }
    if (transform.pars) {
      init <- .Call(stats:::C_ARIMA_Invtrans, init, arma)
      if (arma[2L] > 0) {
        ind <- arma[1L] + 1L:arma[2L]
        init[ind] <- maInvert(init[ind])
      }
      if (arma[4L] > 0) {
        ind <- sum(arma[1L:3L]) + 1L:arma[4L]
        init[ind] <- maInvert(init[ind])
      }
    }
    trarma <- .Call(stats:::C_ARIMA_transPars, init, arma, transform.pars)
    mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa,
                     SSinit)
    res <- optim(init[mask], armafn, method = optim.method,
               hessian = TRUE, control = optim.control, trans = as.logical(transform.pars))
    if (res$convergence > 0)
      warning(gettextf("possible convergence problem: optim gave code = %d",
                       res$convergence), domain = NA)
    coef[mask] <- res$par
    if (transform.pars) {
      if (arma[2L] > 0L) {
        ind <- arma[1L] + 1L:arma[2L]
        if (all(mask[ind]))
          coef[ind] <- maInvert(coef[ind])
      }
      if (arma[4L] > 0L) {
        ind <- sum(arma[1L:3L]) + 1L:arma[4L]
        if (all(mask[ind]))
          coef[ind] <- maInvert(coef[ind])
      }
      if (any(coef[mask] != res$par)) {
        oldcode <- res$convergence
        res <- optim(coef[mask], armafn, method = optim.method,
                     hessian = TRUE, control = list(maxit = 0L,
                                                    parscale = optim.control$parscale), trans = TRUE)
        res$convergence <- oldcode
        coef[mask] <- res$par
      }
      A <- .Call(stats:::C_ARIMA_Gradtrans, as.double(coef), arma)
      A <- A[mask, mask]
      var <- crossprod(A, solve(res$hessian * n.used, A))
      coef <- .Call(stats:::C_ARIMA_undoPars, coef, arma)
    }
    else var <- solve(res$hessian * n.used)
    trarma <- .Call(stats:::C_ARIMA_transPars, coef, arma, FALSE)
    mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa,
                     SSinit)
    val <- if (ncxreg > 0L) {
      arimaSS(x - xreg %*% coef[narma + (1L:ncxreg)], mod)
    } else arimaSS(x, mod)
    sigma2 <- val[[1L]][1L]/n.used
  }
  value <- 2 * n.used * res$value + n.used + n.used * log(2 * pi)
  aic <- NA
  aic <- if (method != "CSS") value + 2 * sum(mask) + 2
  if (ncxreg > 0L) {
    nm <- c(nm, cn)
    if (!orig.xreg) {
      ind <- narma + 1L:ncxreg
      coef[ind] <- S$v %*% coef[ind]
      A <- diag(narma + ncxreg)
      A[ind, ind] <- S$v
      A <- A[mask, mask]
      var <- A %*% var %*% t(A)
    }
  }
  resid <- val[[2L]]
  structure(list(coef = coef, sigma2 = sigma2, var.coef = var,
                 mask = mask, loglik = -0.5 * value, aic = aic, arma = arma,
                 residuals = resid, call = match.call(),
                 code = res$convergence, n.cond = ncond, nobs = n.used,
                 model = mod), class = "Arima")
}

# predict.Arima <- function (object, n.ahead = 1L, newxreg = NULL, se.fit = TRUE,
#                            ...)
# {
#   myNCOL <- function(x) if (is.null(x))
#     0
#   else NCOL(x)
#   rsd <- object$residuals
#   xr <- object$call$xreg
#   xreg <- if (!is.null(xr))
#     eval.parent(xr)
#   else NULL
#   ncxreg <- myNCOL(xreg)
#   if (myNCOL(newxreg) != ncxreg)
#     stop("'xreg' and 'newxreg' have different numbers of columns")
#   xtsp <- tsp(rsd)
#   n <- length(rsd)
#   arma <- object$arma
#   coefs <- object$coef
#   narma <- sum(arma[1L:4L])
#   if (length(coefs) > narma) {
#     if (names(coefs)[narma + 1L] == "intercept") {
#       newxreg <- cbind(intercept = rep(1, n.ahead), newxreg)
#       ncxreg <- ncxreg + 1L
#     }
#     xm <- if (narma == 0)
#       drop(as.matrix(newxreg) %*% coefs)
#     else drop(as.matrix(newxreg) %*% coefs[-(1L:narma)])
#   }
#   else xm <- 0
#   if (arma[2L] > 0L) {
#     ma <- coefs[arma[1L] + 1L:arma[2L]]
#     if (any(Mod(polyroot(c(1, ma))) < 1))
#       warning("MA part of model is not invertible")
#   }
#   if (arma[4L] > 0L) {
#     ma <- coefs[sum(arma[1L:3L]) + 1L:arma[4L]]
#     if (any(Mod(polyroot(c(1, ma))) < 1))
#       warning("seasonal MA part of model is not invertible")
#   }
#   z <- KalmanForecast(n.ahead, object$model)
#   pred <- ts(z[[1L]] + xm, start = xtsp[2L] + deltat(rsd),
#              frequency = xtsp[3L])
#   if (se.fit) {
#     se <- ts(sqrt(z[[2L]] * object$sigma2), start = xtsp[2L] +
#                deltat(rsd), frequency = xtsp[3L])
#     list(pred = pred, se = se)
#   }
#   else pred
# }


