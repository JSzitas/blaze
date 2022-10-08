arma <- function (x, order = c(1, 1), lag = NULL, coef = NULL, include.intercept = TRUE,
          series = NULL, qr.tol = 1e-07, ...)
{
  seqN <- function(N) {
    if (0 == length(N))
      NULL
    else if (N <= 0)
      NULL
    else seq(N)
  }
  err <- function(coef) {
    u <- double(n)
    u[seqN(max.order)] <- 0
    u <- .C(tseries_arma, as.vector(x, mode = "double"),
            u = as.vector(u), as.vector(coef, mode = "double"),
            as.integer(lag$ar), as.integer(lag$ma), as.integer(ar.l),
            as.integer(ma.l), as.integer(max.order), as.integer(n),
            as.integer(include.intercept))$u
    return(sum(u^2))
  }
  resid <- function(coef) {
    u <- double(n)
    u[seqN(max.order)] <- 0
    u <- .C(tseries_arma, as.vector(x, mode = "double"),
            u = as.vector(u), as.vector(coef, mode = "double"),
            as.integer(lag$ar), as.integer(lag$ma), as.integer(ar.l),
            as.integer(ma.l), as.integer(max.order), as.integer(n),
            as.integer(include.intercept))$u
    return(u)
  }
  arma.init <- function() {
    k <- round(1.1 * log(n))
    e <- na.omit(drop(ar.ols(x, order.max = k, aic = FALSE,
                             demean = FALSE, intercept = include.intercept)$resid))
    ee <- embed(e, max.order + 1)
    xx <- embed(x[-(1:k)], max.order + 1)
    if (include.intercept == TRUE) {
      if (is.null(lag$ar))
        coef <- lm(xx[, 1] ~ ee[, lag$ma + 1])$coef
      else if (is.null(lag$ma))
        coef <- lm(xx[, 1] ~ xx[, lag$ar + 1])$coef
      else coef <- lm(xx[, 1] ~ xx[, lag$ar + 1] + ee[,
                                                      lag$ma + 1])$coef
      coef <- c(coef[-1], coef[1])
    }
    else {
      if (is.null(lag$ar))
        coef <- lm(xx[, 1] ~ ee[, lag$ma + 1] - 1)$coef
      else if (is.null(lag$ma))
        coef <- lm(xx[, 1] ~ xx[, lag$ar + 1] - 1)$coef
      else coef <- lm(xx[, 1] ~ xx[, lag$ar + 1] + ee[,
                                                      lag$ma + 1] - 1)$coef
    }
    return(coef)
  }
  if (!is.null(order) && !is.null(lag))
    warning("order is ignored")
  if (is.null(order) && is.null(lag))
    stop("order or lag must be given")
  if (is.null(lag) && !is.null(order))
    lag <- list(ar = seqN(order[1]), ma = seqN(order[2]))
  lag$ar <- unique(lag$ar)
  lag$ma <- unique(lag$ma)
  max.order <- max(unlist(lag), 0)
  ar.l <- length(lag$ar)
  ma.l <- length(lag$ma)
  if (NCOL(x) > 1)
    stop("x is not a vector or univariate time series")
  if (is.null(series))
    series <- deparse(substitute(x))
  ists <- is.ts(x)
  x <- as.ts(x)
  xfreq <- frequency(x)
  if (any(is.na(x)))
    stop("NAs in x")
  if (ists)
    xtsp <- tsp(x)
  n <- length(x)
  if (!is.null(unlist(lag)))
    if ((min(unlist(lag)) < 1) || (max(unlist(lag)) > (n -
                                                       1)))
      stop("invalid lag")
  ncoef <- length(unlist(lag)) + as.numeric(include.intercept)
  if (is.null(coef)) {
    if (!is.null(unlist(lag)))
      coef <- arma.init()
    else coef <- 0
  }
  if (length(coef) != ncoef)
    stop("invalid coef")
  md <- optim(coef, err, gr = NULL, hessian = TRUE, ...)
  coef <- md$par
  rank <- qr(md$hessian, qr.tol)$rank
  if (rank != ncoef) {
    vc <- matrix(NA, nrow = ncoef, ncol = ncoef)
    warning("singular Hessian")
  }
  else {
    vc <- 2 * md$value/n * solve(md$hessian)
    if (any(diag(vc) < 0))
      warning("Hessian negative-semidefinite")
  }
  e <- resid(coef)
  e[seqN(max.order)] <- NA
  f <- x - e
  if (ists) {
    attr(e, "tsp") <- xtsp
    attr(e, "class") <- "ts"
    attr(f, "tsp") <- xtsp
    attr(f, "class") <- "ts"
  }
  nam.ar <- if (!is.null(lag$ar))
    paste("ar", lag$ar, sep = "")
  else NULL
  nam.ma <- if (!is.null(lag$ma))
    paste("ma", lag$ma, sep = "")
  else NULL
  nam.int <- if (include.intercept)
    "intercept"
  else NULL
  nam.coef <- c(nam.ar, nam.ma, nam.int)
  names(coef) <- nam.coef
  colnames(vc) <- rownames(vc) <- nam.coef
  arma <- list(coef = coef, css = md$value, n.used = n, residuals = e,
               fitted.values = f, series = series, frequency = xfreq,
               call = match.call(), vcov = vc, lag = lag, convergence = md$convergence,
               include.intercept = include.intercept)
  class(arma) <- "arma"
  return(arma)
}
