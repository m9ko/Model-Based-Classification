# Contents of file `gibbs-sampler.R`

#' Gibbs sampler for model based classification.
#'
#' @param z0 Initial estimate of the cluster membership indicator vector of length `N`.
#' @param y Observed data of dimension `N x q`.
#' @param V Observed Fisher information array of dimension `q x q x N`.
#' @param mu Cluster means matrix of dimension `K x q`.
#' @param Sigma Cluster variance array of dimension `q x q x K`.
#' @param rho Gaussian Mixture model weight of vector of length `K`.
#' @param v_k Baseline degree of freedom for Sigma. Either a vector of length `K` or a scalar.
#' @param omega Baseline scale matrix for Sigma. Either a `q x q` matrix or a `q x q x K` array.
#' @param N Number of observations (positive integer).
#' @param q Size of \eqn{\theta_i} (positive integer).
#' @param K Number of clusters (positive integer).
#' @param iter Number of Gibbs sampling iterations, default to 1500.
#' @param burn Number of warm-up iterations, default to 500.
#' @param returnOpt Boolean for whether \eqn{\Theta} and \eqn{z} should be returned. 
#' @param reord Boolean for whether the original order of clusters is recovered. Please see `relabel()`.
#' @param true_mu Ground truth means matrix of dimension `K x q`, if known, for reordering. Please see `relabel()`.
#'
#' @return A list consisting of:
#' \eqn{\rho} (a `K x (iter-burn)` matrix),
#' \eqn{\mu} (a `K x q x (iter-burn)` array),
#' \eqn{\Sigma} (a `q x q x K x (iter-burn)` array),
#' \eqn{\Lambda} (a `N x K` matrix),
#' \eqn{\Theta} (a `N x q x (iter-burn)` array, optional), and
#' \eqn{z} (a `N x (iter-burn)` matrix, optional).
#'
#' @examples
#' # This is a simulated example
#' N <- 100
#' q <- 5
#' K <- 4
#' true_z <- sample(1:K, N, TRUE)
#' true_mu <- matrix(rnorm(q*K, 0, 5), ncol = q)
#' true_Sigma <- replicate(K, expr = {mniw::rwish(1, diag(q)*runif(1, 0.01, 0.1), q+1)})
#' true_Theta <- mniw::rmNorm(N, true_mu[true_z,], true_Sigma[,,true_z])
#' V <- replicate(N, expr = {mniw::rwish(1, diag(q)*0.05, q+1)})
#' y <- mniw::rmNorm(N, true_Theta, V)
#' v_k <- q + 2
#' Var_y <- cov(y)
#' params <- init_params(y, K, type = "kmeans")
#' # usually, for the Gibbs sampler to converge, more iterations are needed
#' gs <- gibbs_sampler(params$z0, params$Theta, V, params$mu, params$Sigma,
#'                         params$rho, v_k, Var_y, N, q, K, iter=50, burn=10)
#'
#' @export
gibbs_sampler <- function(z0, y, V, mu, Sigma, rho, v_k, omega, N, q, K, iter = 1500, burn = 500, returnOpt = FALSE, reord = FALSE, true_mu = NULL) {
  actual_iter = iter - burn
  # placeholding arrays
  rho_matrix <- matrix(NA, nrow = K, ncol = actual_iter)
  mu_matrix <- array(matrix(NA, nrow = K, ncol = q), dim = c(K,q,actual_iter))
  Sigma_matrix <- array(matrix(NA, nrow = q, ncol = q), dim = c(q,q,K,actual_iter))
  Lambda_matrix <- matrix(0, nrow = N, ncol = K)
  if (returnOpt) {
    Theta_matrix <- array(matrix(NA, nrow = N, ncol = q), dim = c(N,q,actual_iter))
    z_matrix <- matrix(NA, nrow = N, ncol = actual_iter)
  }
  z <- z0 # initial z
  prev_mu <- mu # initial mu
  # run the iterations
  for (i in 1:iter) {
    # parameter updates
    Theta <- update_Theta(y, V, z, mu, Sigma, N)
    mu <- update_mu(z, prev_mu, Theta, Sigma, N, q, K)
    prev_mu <- mu
    Sigma <- update_Sigma(z, Theta, mu, q, K, v_k, omega)
    rho <- update_rho(z, K)
    zlist <- update_z(Sigma, Theta, mu, rho, N, q, K)
    z <- zlist$updated_z
    # save values after burn
    if (i > burn) {
      idx <- i - burn
      rho_matrix[,idx] <- rho
      mu_matrix[,,idx] <- mu
      Sigma_matrix[,,,idx] <- Sigma
      Lambda_matrix <- Lambda_matrix + zlist$lambda
      if (returnOpt) {
        Theta_matrix[,,idx] <- Theta
        z_matrix[,idx] <- z
      }
    }
  }
  # take average of lambda across iterations
  Lambda_matrix <- Lambda_matrix / matrix(rep(actual_iter, N*K), nrow = N, ncol = K)
  if (reord & !is.null(true_mu)) {
    new_order <- relabel(true_mu, mu_matrix, K)
    rho_matrix <- rho_matrix[new_order$order_param,]
    mu_matrix <- mu_matrix[new_order$order_param,,]
    Sigma_matrix <- Sigma_matrix[,,new_order$order_param,]
    Lambda_matrix <- Lambda_matrix[,new_order$order_param]
    if (returnOpt) {z_matrix <- apply(z_matrix, 2, function(x) {new_order$order_z[x]})}
  }
  # parameter list
  params_list <- list(rho = rho_matrix, mu = mu_matrix, Sigma = Sigma_matrix, Lambda = Lambda_matrix)
  # return Theta and z matrix as well if returnOpt = TRUE
  if (returnOpt) {
    params_list[["Theta"]] <- Theta_matrix
    params_list[["z"]] <- z_matrix
  }
  return(params_list)
}


#' Labels for recovering cluster order according to the ground truth mu.
#'
#' @param true_mu Ground truth mu of dimension `K x q`. Can be cluster means estimated with known membership of each observation.
#' @param sample_mu MCMC samples of mu of dimension `K x q x iteration` from [mbc::gibbs_sampler()].
#' @param K Number of clusters (positive integer).
#' @return A list consisting `order_param` (ordering for parameters) and `order_z` (ordering for z), both of length `K`.
#' @details The original order of the shuffled cluster labels is estimated with the ground truth mu and MCMC sample mu, by sequentially finding the minimum difference between the two.
relabel <- function(true_mu, sample_mu, K) {
  avg_mu <- apply(sample_mu, c(1,2), mean) # average of the MCMC sample mu
  order_param <- order_z <- numeric(K) # initialize order for parameter and z
  # matrix of squared difference between each true mu and each average mu
  dif <- sapply(1:K, function(k) colMeans((replicate(K, true_mu[k,]) - t(avg_mu))^2))
  # iterate for K times
  for (i in 1:K) {
    # get matrix coordinates of the minimum value in dif (minimum difference)
    scalar <- which.min(dif)
    x <- scalar %% K
    y <- scalar %/% K + 1
    if (x == 0) {
      x <- K
      y <- y - 1
    }
    order_param[y] <- x # assign x to index y
    order_z[x] <- y # assign y to index x
    # update dif to remove possibility for being chosen again
    dif[x,] <- Inf
    dif[,y] <- Inf
  }
  list(order_param = order_param, order_z = order_z)
}


#' Estimate the classification accuracy with MCMC estimate of membership probability and ground truth membership.
#'
#' @param prob MCMC samples of \eqn{\Lambda} of dimension `N x K` from [mbc::gibbs_sampler()].
#' @param true_z Ground truth \eqn{z} of dimension `N`.
#' @param order Optional true order of the membership indicator with length `K` obtained from [mbc::relabel()].
#' @return A list consisting `prob`, the probability of each observation being successfully classified with dimension `N x 2`,
#'         and `expectation` the overall probability of successfully classifying with length `2`.
#' @export 
estim_accuracy <- function(prob, true_z, order) {
  if (missing(order)) {order <- 1:ncol(prob)} # if prob is already in order
  reordered <- prob[,order]
  # use indicator matrix generated from true membership
  indicator <- as.matrix(sparseMatrix(i = true_z, j = 1:nrow(prob), x = 1))
  choice_prob <- diag(reordered %*% indicator)
  df <- data.frame('Success Rate' = choice_prob, 'Missed Rate' = 1 - choice_prob)
  return(list(prob = df, expectation = colMeans(df)))
}
