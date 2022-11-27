period <- function(x)
{
  n <- length(x)
  spec <- stats::spec.ar(c(x),plot=FALSE)
  if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
  {
    period <- round(1/spec$freq[which.max(spec$spec)])
    if(period==Inf) # Find next local maximum
    {
      j <- which(diff(spec$spec)>0)
      if(length(j)>0)
      {
        nextmax <- j[1] + which.max(spec$spec[j[1]:500])
        period <- round(1/spec$freq[nextmax])
      }
      else
        period <- 1
    }
  }
  else
    period <- 1
  return(period)
}

find_seasonalities <- function( y, max_iter = 5, aggregator = sum, upper_limit = 1500 ) {

  periods <- list()
  for( iter in seq_len(max_iter) ) {
    last_period <- period(y)
    if( last_period <= 1 ){
      break;
    }
    periods[[iter]] <- last_period
    y <- stats::aggregate(
      stats::ts(y, freq = last_period), # where last_period is last infered periodicity
      nfrequency = 1, # nfrequency always set to 1
      FUN = aggregator # ie mean
    )
  }
  x <- cumprod( unlist(periods))
  x[ x < upper_limit ]
}

x <- tsibbledata::vic_elec$Demand
# seas <- find_seasonalities(x)

make_seasonal_dummy <- function(seas) {
  lapply( seas,
          function(i){
            data.frame(diag(x = 1, nrow = i))
            }
          )
}
# seasonalities <- make_seasonal_dummy(seas)

# seas2 <- list(seasonalities[[1]])

seasonal_means <- function(x, upper_limit = 20000, max_iter = 10) {
  seas <- find_seasonalities(x, max_iter = max_iter, upper_limit = upper_limit)
  seas <- make_seasonal_dummy(seas)

  res <- purrr::map( seas, function(i) {
      purrr::map_dbl(i, function(j) {
        indices <- rep(j, ceiling(length(x)/length(j)))[seq_len(length(x))]
        unname(sum(x * indices)/sum(indices))
        })
  })
  return(res)
}


max_length <- floor(length(x)/max(purrr::map_dbl(rs,length))
) * max(purrr::map_dbl(rs,length))

x[seq_len(max_length)] %>%
  plot.ts

rep(unlist(rs[[1]]), max_length/length(rs[1]) )[
  seq_len( max_length)] %>%
  lines(col = "green")
rep(unlist(rs[[2]]), max_length/length(rs[[2]]) )[
  seq_len( max_length)] %>%
  lines(col = "red")
rep(unlist(rs[[3]]), max_length/length(rs[[3]]) )[
  seq_len( max_length)] %>%
  lines(col = "brown")

rep(unlist(rs[[4]]), max_length/length(rs[[4]]) )[
  seq_len( max_length)] %>%
  lines(col = "blue")

df <- data.frame( target = x[seq_len(max_length)],
                  seas1 = rep(unlist(rs[[1]]), max_length/length(rs[1]) )[
                    seq_len( max_length)],
                  seas2 = rep(unlist(rs[[2]]), max_length/length(rs[2]) )[
                    seq_len( max_length)],
                  seas3 = rep(unlist(rs[[3]]), max_length/length(rs[3]) )[
                    seq_len( max_length)],
                  seas4 = rep(unlist(rs[[4]]), max_length/length(rs[4]) )[
                    seq_len( max_length)]
                  )







