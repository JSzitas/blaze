

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

get_lambda_path <- function( y, lags = 1:2, epsilon = .000001, K = 100 ) {

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
  full_y <- y
  y <- head(y, length(y)-valid)

  lambda_path <- get_lambda_path(y, lags = max(c(1:season, order[-2L]) ), K = n_lambda)
  results <- list()
  p <- 1
  for( lambda in lambda_path ) {

    tryCatch({
      model <- suppressWarnings(arima_ridge(y, order, lambda = lambda,...))
      # this model needs to be refitted afterwards otherwise this is kihda meaningless
      score <- scoring_rule( c(predict(model, n.ahead = valid)$pred), valid_y)
      model <- suppressWarnings(arima_ridge(full_y, order, lambda = lambda,...))

      results[[p]] <- list( model = model,
                            valid_score = score,
                            lambda = lambda)

      p <- p+1
    }, error = function(e){})
  }

  return(results)
}

train_lynx <- lynx[1:100]
test_lynx <- lynx[101:114]

test_order = c(4,2,2)
test_s_order = list(order = c(1,0,0),12)


mdls <- auto_ridge_arima(train_lynx, order = test_order, seasonal = test_s_order )

preds <- lapply( mdls, function(i) { predict(i[["model"]], 14)$pred })
valid_scores <- lapply( mdls, function(i) { i[["valid_score"]] })
errs <- lapply( preds, function(i) { rmse( i, test_lynx) })

lambdas <- sapply( mdls, function(i) i[["lambda"]] )

df <- lapply( seq_along(mdls), function(i) {
  data.frame(
    predictions = preds[[i]],
    lambda = lambdas[[i]],
    index = seq_along(preds[[i]])
  )
}) %>%
  dplyr::bind_rows()

data.frame( predictions = test_lynx, index = seq_along(1:14), lambda = 0 ) %>%
  ggplot2::ggplot( ggplot2::aes( x = index, y = predictions, group = lambda ) ) +
  ggplot2::geom_line(size = 1.6) +
  ggplot2::geom_line( data = df, ggplot2::aes(color = lambda)) +
  ggplot2::geom_line( data = data.frame(predictions = predict(arima_ridge( train_lynx,
                                                                           lambda = 0,
                                                                           order = test_order,
                                                                           seasonal = test_s_order),
                                                              14)$pred,
                                        index = seq_along(1:14),
                                        lambda = 0 ),
                      color = "red", size = 1.2) +
  ggplot2::geom_line( data = df %>%
                        dplyr::filter(lambda == lambdas[[which.min(unlist(errs))]]),
                      color = "darkgreen", size = 1.2
                     ) +
  ggplot2::geom_line( data = df %>%
                        dplyr::filter(lambda == lambdas[[which.min(unlist(valid_scores))]]),
                      color = "purple", size = 1.2
  )

# errors on test set








