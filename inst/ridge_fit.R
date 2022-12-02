

lag <- function( x, lags = 1 ) {
  X <- matrix( NA, ncol =  length(lags), nrow = length(x))

  current_lag <- 1
  for( lag in lags ) {
    id <- 1:(length(x)-lag)
    X[id+lag,current_lag] <- x[id]
    current_lag <- current_lag + 1
  }
  return(X)
}

get_lambda_path <- function( y, lags = 1:5, epsilon = .000001, K = 100 ) {

  x <- lag(y, lags)
  x <- na.omit(x)
  y <- y[(max(lags)+1):length(y)]

  ## Standardize variables: (need to use n instead of (n-1) as denominator)
  mysd <- function(z) sqrt(sum((z-mean(z))^2)/length(z))
  sx <- scale(x, scale = apply(x, 2, mysd))
  sx <- as.matrix(sx, ncol = ncol(x), nrow = nrow(x))

  ## Calculate lambda path (first get lambda_max):
  lambda_max <- max(abs(colSums(sx*y)))/length(y)
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                              length.out = K)), digits = 10)
  return(lambdapath)
}

rmse <- function( y_pred, y ) {
  sqrt(mean( (y-y_pred)^2))
}

auto_ridge_arima <- function( y, order, season = 1, valid = 10,#season,
                              n_lambda = 100,
                              scoring_rule = rmse,... ) {

  valid_y <- tail(y, valid)
  y <- head(y, length(y)-valid)

  lambda_path <- get_lambda_path(y, lags = max(c(1:season, order[-2L]) ), K = n_lambda)
  results <- list()
  p <- 1
  for( lambda in lambda_path ) {

    tryCatch({
      model <- suppressWarnings(arima_ridge(y, order, lambda = lambda,...))
      results[[p]] <- list( model = model,
                            valid_score = scoring_rule( c(predict(model, n.ahead = valid)$pred), valid_y),
                            lambda = lambda)

      p <- p+1
    }, error = function(e){})
  }

  return(results)
}

train_lynx <- lynx[1:100]
test_lynx <- lynx[101:114]

mdls <- auto_ridge_arima(train_lynx, order = c(5,2,3))

errs <- lapply( mdls, function(i) { rmse( predict(i[["model"]], 14)$pred, test_lynx) })
















