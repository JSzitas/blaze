stats:::ar.yw.default
function (x, aic = TRUE, order.max = NULL, na.action = na.fail,
          demean = TRUE, series = NULL, ...)
{
  if (demean) {
    xm <- colMeans(x, na.rm = TRUE)
    x <- sweep(x, 2L, xm, check.margin = FALSE)
  }
  else xm <- rep.int(0, nser)
  n.used <- nrow(x)
  n.obs <- sum(!is.na(x[, 1]))
  order.max <- if (is.null(order.max))
    min(n.obs - 1L, floor(10 * log10(n.obs)))
  else floor(order.max)
  if (order.max < 1L)
    stop("'order.max' must be >= 1")
  else if (order.max >= n.obs)
    stop("'order.max' must be < 'n.obs'")
  xacf <- acf(x, type = "covariance", lag.max = order.max,
              plot = FALSE, demean = demean, na.action = na.pass)$acf
  r <- as.double(drop(xacf))
  # apparently solver toeplitz matrix equation
  z <- .Fortran(C_eureka, as.integer(order.max), r, r,
                coefs = double(order.max^2), vars = double(order.max),
                double(order.max))
  coefs <- matrix(z$coefs, order.max, order.max)
  partialacf <- array(diag(coefs), dim = c(order.max, 1L,
                                           1L))
  var.pred <- c(r[1L], z$vars)
  xaic <- n.obs * log(var.pred) + 2 * (0L:order.max) +
    2 * demean
  maic <- min(aic)
  xaic <- setNames(if (is.finite(maic))
    xaic - min(xaic)
    else ifelse(xaic == maic, 0, Inf), 0L:order.max)
  order <- if (aic)
    (0L:order.max)[xaic == 0L]
  else order.max
  ar <- if (order)
    coefs[order, seq_len(order)]
  else numeric()
  var.pred <- var.pred[order + 1L]
  var.pred <- var.pred * n.obs/(n.obs - (order + 1L))

  res <- list(order = order, ar = ar, var.pred = var.pred)
  res
}


spec.ar <- function (x, n.freq, order = NULL,
          method = "yule-walker", ...)
{
  series <- deparse1(substitute(x))
  x <- na.action(as.ts(x))
  xfreq <- frequency(x)
  nser <- NCOL(x)
  x <- ar(x, is.null(order), order, na.action = na.fail, method = method)
  order <- x$order
  if (missing(n.freq))
    n.freq <- 500
  freq <- seq.int(0, 0.5, length.out = n.freq)
  if (nser == 1) {
    coh <- phase <- NULL
    var.p <- as.vector(x$var.pred)
    spec <- if (order >= 1) {
      cs <- outer(freq, 1L:order, function(x, y) cos(2 *
                                                       pi * x * y)) %*% x$ar
      sn <- outer(freq, 1L:order, function(x, y) sin(2 *
                                                       pi * x * y)) %*% x$ar
      var.p/(xfreq * ((1 - cs)^2 + sn^2))
    }
    else rep.int(var.p/xfreq, length(freq))
  }
  spg.out <- list(freq = freq * xfreq, spec = spec)

  else spg.out
}
