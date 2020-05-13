# Contents of file `nnm-functions.R`

#' Solve method for variance matrices.
#'
#' @param V Variance matrix.
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
#' @return Matrix solving system of equations and optionally the log-determinant.
#' @details This function is faster and more stable than `base::solve()` when `V` is known to be positive-definite.
#' @author Martin Lysy
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V) # cholesky decomposition
  if(missing(x)) x <- diag(nrow(V))
  # solve is O(ncol(C) * ncol(x)) with triangular matrices
  # using backward subsitution
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}


#' Count the number of observations within each cluster.
#'
#' @param z Cluster membership indicator vector of length `N`.
#' @param K Number of clusters (positive integer).
#' @return A named vector of length K consisting of number of observations (`N_k`) in each cluster.
freq_table <- function(z, K) {
  N_k <- sapply(1:K, function(k) sum(z == k))
  names(N_k) <- 1:K
  return(N_k)
}


#' Complete log-likelihood for the normal-normal-mixture model.
#'
#' @param y Observed data of dimension `N x q`.
#' @param Theta Unobserved data matrix of dimension `N x q`.
#' @param V Observed Fisher information array of dimension `q x q x N`.
#' @param Sigma Cluster variance array of dimension `q x q x K`.
#' @param mu Cluster means matrix of dimension `K x q`.
#' @param z Cluster membership indicator vector of length `N`.
#' @param rho Gaussian Mixture model weight of vector of length `K`.
#' @param N Number of observations (positive integer).
#' @param q Size of \eqn{\theta_i} (positive integer).
#' @param K Number of clusters (positive integer).
#' @param v_k Baseline degree of freedom for Sigma. Either a vector of length `K` or a scalar.
#' @param omega Baseline scale matrix for Sigma. Either a `q x q` matrix or a `q x q x K` array.
#' @return A list consisting of a vector of length N for \eqn{z}, and a N x K matrix for \eqn{lambda} (see details).
ll <- function(y, Theta, V, Sigma, mu, z, rho, N, q, K, v_k, omega) {
  if (length(v_k) == 1) {v_k <- rep(v_k, K)} # transform to a vector
  if (length(dim(omega)) == 2) {omega <- replicate(K, omega)} # make into 3D array
  N_k <- freq_table(z, K)
  prior <- obs_ll <- mean_ll <- 0
  for (k in 1:K) {
    S <- Sigma[,,k]
    log_det_S <- log(det(S))
    solveS <- solveV(S)    
    Theta_mu <- Theta[z==k,] - mu[rep(k, N_k[k]),]
    # prior
    prior <- prior + (v_k[k] + q + 1) * log_det_S + sum(diag(solveS %*% omega[,,k]))
    # mean contribution
    mean_ll <- mean_ll + sum(diag(tcrossprod(Theta_mu %*% solveS, Theta_mu)) + log_det_S)
  }
  for (i in 1:N) {
    # observation contribution
    obs_ll <- obs_ll + (y - Theta)[i,] %*% solveV(V[,,i]) %*% (y - Theta)[i,]
  }
  # membership contribution
  memb_ll <- sum(log(rho) * N_k)
  # assembled complete log-likelihood
  -0.5*(prior + obs_ll + mean_ll) + memb_ll
}


#' Sample from a categorical distribution.
#'
#' Performs one draw from a categorical distribution (i.e., a multinomial distribution with size `n = 1`) for each of multiple probability vectors.
#'
#' @param prob An `n_cat x n_prob` matrix of probability vectors, each of which is a column of `prob`.  The entries of `prob` must be nonnegative, but are internally normalized such that each column sums to one.
#' @return A vector of length `n_prob` of integers between 1 and `n_cat` indicating the category draw for each column of `prob`.
#' @author Martin Lysy
rcategorical <- function(prob) {
  if(any(prob < 0)) stop("prob must contain only nonnegative entries.")
  cprob <- apply(prob, 2, cumsum)
  u <- runif(ncol(prob)) * cprob[nrow(prob),]
  apply(sweep(cprob, 2, u) >= 0, 2, which.max)
}

#' Density of the Dirichlet distribution.
#'
#' @param x Observation vector of nonnegative entries which sum to one.
#' @param alpha Weight parameter: a vector of nonnegative entries the same length as `x`.
#' @param log Logical; whether to evaluate the density on the log scale.
#' @return The density or log-density evaluated at the inputs (scalar).
#' @author Martin Lysy
ddirichlet <- function(x, alpha, log = FALSE) {
  ld <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  ld <- ld + sum((alpha-1) * log(x))
  if(!log) ld <- exp(ld)
  ld
}

#' Random sampling from the Dirichlet distribution.
#'
#' @param n Number of random draws.
#' @param alpha Weight parameter: a vector of nonnegative entries.
#' @return A matrix of size `n x length(alpha)` of which each row is a random draw.
#' @author Martin Lysy
rdirichlet <- function(n, alpha) {
  K <- length(alpha) # number of categories
  X <- matrix(rgamma(n*K, shape = alpha), K, n)
  drop(t(sweep(X, 2, colSums(X), "/")))
}
