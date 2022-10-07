# optimise(lambda_coef_var, c(lower, upper),
#          x = x, .period = max(.period, 2))$minimum)


x <- forecast::wineind

lambda_coef_var <- function(lambda) {
  x <- split(x, (seq_along(x) - 1)%/%2)
  mu_h <- purrr::map_dbl(x, mean, na.rm = TRUE)
  sig_h <- purrr::map_dbl(x, sd, na.rm = TRUE)
  rat <- sig_h/mu_h^(1 - lambda)
  stats::sd(rat, na.rm = TRUE)/mean(rat, na.rm = TRUE)
}

# guerrero
guer <- feasts::guerrero(x)

bisection <- function( fun, lower, upper, tol = 0.01, max_iter = 1000 ) {
  iter <- 1

  # val_upper <- fun( upper )
  val_lower <- fun( lower )
  # find a midpoint between lower and upper
  # mid <- (upper+lower)/2
  print(val_lower)
  while( iter <= max_iter && ((upper - lower) > tol)) {
    mid <- (upper+lower)/2
    print(paste0( "Lower: ",lower, " Mid: ", mid, " Upper: ", upper))
    val_mid <- fun( mid )
    print( paste0(" f(mid)= ", val_mid, " f(lower)= ", val_lower))
    if( val_mid == 0 ) {
      break
    }
    else if( (val_lower * val_mid) > 0 ) {
      lower = mid
      val_lower = val_mid
    }
    else {
      upper = mid
    }
    iter <- iter + 1
  }
  return(mid)
}


ridders_method <- function(  fun, lower, upper, tol = 0.00001, max_iter = 1 ) {

  val_lower <- fun(lower)
  val_upper <- fun(upper)

  for( iter in seq_len(max_iter) ) {
    guess = (lower + upper)/2
    val_guess = fun(guess)
    if( val_guess^2 < (val_lower * val_upper) ) {
      stop("No real roots that can be found with Ridders method.")
    }
    new_guess = guess + (guess - lower) * sign( val_lower - val_upper)*val_guess/sqrt( val_guess^2 - val_lower * val_upper)
    # print(new_guess)
    if( min( c(abs(new_guess - upper), abs(new_guess - lower) )) < tol) {
      return(new_guess)
    }
    update_val = fun(new_guess)
    if( sign(val_guess) != sign(update_val) ) {
      # update guesses
      lower = guess
      val_lower = val_guess
      upper = new_guess
      val_upper = update_val
    }
    else if( sign(val_upper) != sign(update_val) ) {
      lower = new_guess
      val_lower = update_val
    }
    else {
      upper = new_guess
      val_upper = update_val
    }
  }
  return(new_guess)
}
