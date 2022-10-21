library(ggplot2)
library(gganimate)
nudge <- function( particles, scores ) {
  df <- data.frame( scores = scores, particles = particles )
  df <- df[order(df$scores),]
  best <- df[1,2]
  return( rowMeans( cbind(df[,2], best)))
}

swarm_drainer <- function( fun, lower, upper, n_particles = 100, tol = 0.001, max_iter = 100 ) {
  # create particles between lower and upper
  particles <- runif(n_particles, lower, upper)
  particle_mat <- matrix(NA, ncol =n_particles, nrow = max_iter)
  score_mat <- matrix(NA, ncol =n_particles, nrow = max_iter)
  iter <- 1
  while(iter <= max_iter) {
    particle_mat[iter,] <- particles
    # get particle scores
    scores <- sapply( particles, fun)
    score_mat[iter,] <- scores
    # update particles by "nudging" them towards better particles - they travel half the direction to
    # the nearest best particle
    particles <- c(nudge( particles, scores ))
    if( all((particles - particles[1] < tol)) ) {
      break;
    }
    iter <- iter + 1
  }
  return(list( parts = particle_mat, scores = score_mat))
}

good_fun <- function( x ){ -x^7-2.3*x^6-x^5-3*x^4+12/17*x^3-1/9*x^2+1/9*x * sign( x^3 ) }

# animate particles
swarm_drainer( lambda_coef_var, -20, 20, n_particles = 100, tol = 0.0001, max_iter = 250 ) -> rs

df <- rs[["parts"]]
df <- as.data.frame(na.omit(df))
df <- dplyr::mutate(df, index = seq_len(nrow(df)))
df <- tidyr::pivot_longer( df, cols = - tidyselect::last_col())
df <- dplyr::rename(df, "position" = "value")

df2 <- rs[["scores"]]
df2 <- as.data.frame(na.omit(df2))
df2 <- dplyr::mutate(df2, index = seq_len(nrow(df2)))
df2 <- tidyr::pivot_longer( df2, cols = - tidyselect::last_col())
df2 <- dplyr::rename(df2, "score" = "value")

df <- dplyr::left_join( df, df2, by = c("index","name") )

anim <- ggplot(df, aes(x = position, y = score, color = name)) +
  geom_point() +
  transition_states(index, transition_length = 2, state_length = 1) +
  enter_fade() +
  exit_fade()
animate(anim)
