# Curve reconstruction from noisy points
# 
# Algorithm from:
#   Curve reconstruction from unorganized points
#   In-Kwon Lee
#   Computer Aided Geometric Design (2000)
# 
# Author: Kieran Campbell, University of Oxford <kieranrcampbell@gmail.com>

#' Curve reconstruction from noisy points
reconstruct <- function(X, h = 0.02, method=c('dist','mst','corr')) {
  names(X) <- c('x','y')
  method <- match.arg(method)
  
  W <- weight_matrix(X, h, method)
  Y <- data.frame(t(sapply(1:dim(Y)[1], point_transformation, X, W)))
  
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

weight_matrix <- function(X, h, method) {
  r <- as.matrix(dist(X))
  w <- 2 * r * r * r / h^3 - 3 * r * r / h^2 + 1  
  
  if(method == 'dist'){
    w <- w * 1 * (r < h)
  } else if(method == 'mst') {
    w <- w * call_collect(X, h)
  }
  return( w )
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

#' Improved moving least-squares using minimum spanning trees
#' 
#' 
call_collect <- function(X, h=0.02) {
  D <- as.matrix(dist(X))
  N_cells <- dim(X)[1]
  g <- graph.adjacency(D, mode='undirected', weighted=TRUE)
  mst <- minimum.spanning.tree(g)
  
  ## can plot with
  ## plot(mst, vertex.size=1, vertex.label=NA, layout=as.matrix(X))
  
  ## only use vertices less than h
  use_vertices <- D < h
  uv <- apply(use_vertices, 1, which)
  
  W <- matrix(0, ncol=N_cells, nrow=N_cells)
  for(i in 1:N_cells) {
    A <- collect(i, i, mst = mst, uv = uv)
    W[i,A] <- 1 
  }
  
  print('returning from call_collect')
  return( W )
}

collect <- function(P, P_star, mst, uv) {
  env <- new.env()
  env$A <- NULL
  
  collect1 <- function(P, P_star, mst, uv, env) {
    env$A <- c(P, env$A)
    edges  <- neighbors(mst, P)
    for(P_j in edges) {
      if( (!(P_j %in% env$A)) && (P_j %in% uv[[P_star]]) ) {
        ## P_j qualifies to be in A and its neighbours explored
        collect1(P_j, P_star, mst, uv, env)
      }
    }
  }
  
  collect1(P, P_star, mst, uv, env)
  return( env$A )
}



