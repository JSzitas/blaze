

inv_1 <- function( roots ) {
  x <- 1
  for( r in roots ) {
    # resize x - add 1 to length
    # by padding with a zero - this is the same as initializing the zero
    # beforehand and using restricted length
    x <- c(x, 0) - c(0,x)/r
  }
  return(x)
}

# this says, zero out the first element and only do div over other than first
# element
div_non_firt_elem <- function( x, r ) {
  # result <- rep(NA, length(x)+1)
  # result[1] <- 0
  for( i in seq_len(length(x)) ) {
    result[i+1] = x[i]/r
  }
  return(result)
}

inv_2 <- function(roots) {
  result <- rep(0, length(roots)+1)
  result[1] <- 1

  for( i in seq_len(length(roots)) ) {
    r <- roots[i]
    result[seq_len(i+1)] <- result[seq_len(i+1)] - c( 0, result[seq_len(i)]/r )
  }
  return(result)
}

inv_2 <- function(roots) {
  result <- rep(0, length(roots)+1)
  result[1] <- 1

  for( i in seq_len(length(roots)) ) {
    r <- roots[i]
    result[seq_len(i+1)] <- result[seq_len(i+1)] - c( 0, result[seq_len(i)]/r )
  }
  return(result)
}

inv_3 <- function(roots) {
  result <- rep(0, length(roots)+1)
  result[1] <- 1

  for( i in seq_len(length(roots)) ) {
    r <- roots[i]
    temp <- result
    for( k in seq_len(i+1) ) {
      if( k != 1 ) {
        result[k] = temp[k] - temp[k-1]/r
      }
    }
    # result[seq_len(i+1)] <- result[seq_len(i+1)] - c( 0, result[seq_len(i)]/r )
    # result[seq_len(i+1)] -= c( 0, result[seq_len(i)]/r )
  }
  return(result)
}








