# Curve reconstruction from noisy points
# 
# Algorithm from:
#   Curve reconstruction from unorganized points
#   In-Kwon Lee
#   Computer Aided Geometric Design (2000)
# 
# Author: Kieran Campbell, University of Oxford <kieranrcampbell@gmail.com>

#' Curve reconstruction from noisy points
reconstruct <- function(X, h = 10, niter = 5) {
  W <- weight_matrix(X, h)
  Y <- X

  for(i in 1:niter) {
      Y <- t(sapply(1:dim(Y)[1], point_transformation, X, W))
  }
  
  colnames(Y) <- colnames(X)
  return( Y )
}

point_transformation <- function(index, X, W) {

  ## subset data with non-zero weights
  weights <- W[index,]
  p <- X[index,]
  points <- X[weights > 0,]
  weights <- weights[weights > 0]
  
  ## initial fit and transformation to P*
  fit <- lm(points[,2] ~ points[,1], weights=weights)
  m <- coef(fit)[2]
  R <- rotation_from_gradient(m)
  
  x_star <- t(R %*% t(as.matrix(points)))
  p_star <- t(R %*% t(p))
  
  x_star <- t(apply(x_star, 1, '-', p_star))
  
  ## second non-linear weighted regression
  y <- x_star[,2]
  x <- x_star[,1]
  fit_nonlin <- lm(y ~ x + I(x^2), weights=weights)
  
  C <- coef(fit_nonlin)[1]
  names(C) <- NULL
  
  ## reverse M transformation
  z <- c(0, C) + p_star  
  R_inv <- rotation_from_gradient(-m)
  return( t(R_inv %*% t(z)))
}

weight_matrix <- function(X, h) {
  r <- as.matrix(dist(X))
  w <- 2 * r * r * r / h^3 - 3 * r * r / h^2 + 1  
  w <- w * 1 * (r < h)
}

percentile_r <- function(X) {
  r <- as.vector(dist(X))
  return( quantile(r, p = seq(from=0.01, to=0.1, by=0.01)))
}

#' Rotation matrix from a gradient
rotation_from_gradient <- function(m) {
  theta <- -atan(m)
  matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol=2, byrow=TRUE)
}

plot_transformation <- function(X, Y) {
  library(ggplot2)
  df <- data.frame(rbind(X,Y))
  df$Curve <- rep(c("Original","Transformed"), each=dim(X)[1])
  
  df_seg <- data.frame(cbind(X,Y))
  names(df_seg) <- c('x','y','xend','yend')
  ggplot(df) + geom_point( aes(x=x, y=y, color=Curve)) +
    theme_bw() + geom_segment(data=df_seg, aes(x=x, xend=xend, y=y, yend=yend), alpha=0.5, linetype=2)
}





