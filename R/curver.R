# Curve reconstruction from noisy points
# 
# Algorithm from:
#   Curve reconstruction from unorganized points
#   In-Kwon Lee
#   Computer Aided Geometric Design (2000)
# 
# Author: Kieran Campbell, University of Oxford <kieranrcampbell@gmail.com>

#' Curve reconstruction from noisy points
reconstruct <- function(X, h = 0.02, method=c('dist','mst','corr'), 
                        niter = 1, use_weight_checker = TRUE) {
  names(X) <- c('x','y')
  method <- match.arg(method)
  
  W <- weight_matrix(X, h, method, use_weight_checker)
  Y <- X
  for(i in 1:niter) {
    Y <- data.frame(t(sapply(1:dim(Y)[1], point_transformation, Y, W)))
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

#' Construct the weight matrix using X, h and method
weight_matrix <- function(X, h, method, use_weight_checker) {
  r <- as.matrix(dist(X))
  w <- 2 * r * r * r / h^3 - 3 * r * r / h^2 + 1  
  
  w_all <- w
  if(method == 'dist'){
    w <- w * 1 * (r < h)
  } else if(method == 'mst') {
    w <- w * call_collect(X, h)
  }

#   if(use_weight_checker) w <- weight_matrix_checker(w, w_all)
  return( w )
}

#' Since regressions are performed using nearest neighbours suggested by W, 
#' we need to make sure there are enough so W has more than one point and something
#' relatively stable. As a result, we choose the value 5 to be the minimum number of
#' neighbours a point is allowed 
# weight_matrix_checker <- function(w, w_all, nn = 5) {
#   w_logical <- ceiling(w)
#   if(any(rowSums(w_logical) < nn)) {
#     warning('Some rows have less than 5 nearest neighbours, so increasing all to 5.\nConsider using a larger value for H')
#     problem_rows <- which(rowSums(w_logical) < nn)
#     
#   }
#   return( w )
# }

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
call_collect <- function(X, h=0.02, collect = c(1,2), ...) {
  collect <- match.arg(collect)
  
  D <- as.matrix(dist(X))
  N_cells <- dim(X)[1]
  g <- graph.adjacency(D, mode='undirected', weighted=TRUE)
  mst <- minimum.spanning.tree(g)
  
  ## can plot with
  ## plot(mst, vertex.size=1, vertex.label=NA, layout=as.matrix(X))
  
  W <- matrix(0, ncol=N_cells, nrow=N_cells)
  ## only use vertices less than h
  
  if(collect == 1) {
    use_vertices <- D < h
    uv <- apply(use_vertices, 1, which)
    
    for(i in 1:N_cells) {
      A <- collect(i, i, mst = mst, uv = uv)
      W[i,A] <- 1 
    }
  } else {
    for(i in 1:N_cells) {
      d <- D[i,]
      ## ... should include h_0, rho_0 and epsilon
      A <- collect2(X, d, i, i, mst, ...)
    }
  }
  
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

#' Implements the collect2 algorithm
#'
#' @param d Distance from P_star to all other points
collect2 <- function(X, d, P, P_star, mst, h_0, rho_0, epsilon) {
  h <- h_0
  A <- NULL
  repeat {
    uv <- which(d < h)
    A <- collect(P, P_star, mst, uv)
    if(cor(X[A,]) > rho_0) {
      break
    } else {
      h <- h + epsilon
    }
  }
  return( A )
}


