# Contents of file `update-functions.R`

#' Conditional update for the unobserved data \eqn{\Theta}.
#'
#' @param y Observed data of dimension `N x q`.
#' @param V Observed Fisher information array of dimension `q x q x N`.
#' @param z Cluster membership indicator vector of length `N`.
#' @param mu Cluster means matrix of dimension `K x q`.
#' @param Sigma Cluster variance array of dimension `q x q x K`.
#' @param N Number of observations (positive integer).
#' @param pars Return the means and variances.
#'
#' @return A `N x q` matrix for \eqn{\Theta} drawn from its conditional posterior distribution.
#' If pars is TRUE, a list of parameters of the conditional posterior distribution is returned.
#'
#' @details The conditional posterior distribution of \eqn{\Theta} is given by:
#' ```
#' \eqn{\theta_i | A - \{ \theta_i \} \sim N(G_i(y_i - \mu_{z_i}) + \mu_{z_i}, G_i V_i)},
#' \eqn{G_i = \Sigma_{z_i} (V_i + \Sigma_{z_i})^{-1}}.
#' ```
#' Hence `Theta` is sampled from the random-effects normal distribution given all other parameters.
update_Theta <- function(y, V, z, mu, Sigma, N, pars = FALSE) {
  # return parameters
  if (pars) {
    q <- ncol(mu)
    mean_Theta <- matrix(0, q, N)
    var_Theta <- array(0, dim = c(q,q,N))
    for (i in 1:N) {
      G <- Sigma[,,z[i]] %*% solveV(V[,,i] + Sigma[,,z[i]])
      mean_Theta[,i] <- t(G %*% (y[i,] - mu[z[i],])) + mu[z[i],]
      var_Theta[,,i] <- G %*% V[,,i]
    }
    return(list(mean_Theta = t(mean_Theta), var_Theta = var_Theta))
  }
  rRxNorm(N, y, V, mu[z,], Sigma[,,z])
}


#' Conditional update for \eqn{\mu}.
#'
#' @param z Cluster membership indicator vector of length `N`.
#' @param mu Original cluster means matrix of dimension `K x q`.
#' @param Theta Unobserved data matrix of dimension `N x q`.
#' @param Sigma Cluster variance array of dimension `q x q x K`.
#' @param N Number of observations (positive integer).
#' @param q Size of \eqn{\theta_i} (positive integer).
#' @param K Number of clusters (positive integer).
#' @param pars Return the means and variances.
#'
#' @return A `K x q` matrix for \eqn{\mu} drawn from its conditional posterior distribution.
#' If pars is TRUE, a list of parameters of the conditional posterior distribution is returned.
#'
#' @details The conditional posterior distribution of \eqn{\mu} is given by:
#' ```
#' \eqn{\mu_k | A - \{ \mu_1,...,\mu_K \} \sim N(\bar \theta_k, Sigma_k / N_k)},
#' \eqn{\sum_{i=1}^N 1[z_i = k]}.
#' ```
#' Hence `mu` is sampled from the multivariate normal distribution given all other parameters.
#' Moreover, updating `mu` is allowed for only nonempty clusters (`N_k` > 0).
update_mu <- function(z, mu, Theta, Sigma, N, q, K, pars = FALSE) {
  N_k <- freq_table(z, K)
  nonzero <- which(N_k != 0) # find nonempty clusters
  cluster_Matrix <- as.matrix(sparseMatrix(i = z, j = 1:N, x = 1))[nonzero,] # indicator matrix for observations
  # average of theta for each nonempty cluster
  sample_avg <- (cluster_Matrix / N_k[nonzero]) %*% Theta
  # variance of the average theta for each nonempty cluster
  dnom <- array(rep(N_k[nonzero], rep(q^2, length(nonzero))),
                dim = c(q,q,length(nonzero))) # make denominator (N_k) compatible with dimension of Sigma
  var_avg <- Sigma[,,nonzero] / dnom
  # only update mu for nonempty cluster
  if (pars) {
    mean_mu <- matrix(NA, K, q)
    var_mu <- array(NA, dim = c(q,q,K))
    mean_mu[nonzero, ] <- sample_avg
    var_mu[,,nonzero] <- var_avg
    return(list(mean_mu = mean_mu, var_mu = var_mu))
  }
  mu[nonzero,] <- rmNorm(length(nonzero), sample_avg, var_avg)
  mu
}


#' Conditional update for cluster variance \eqn{\Sigma}.
#'
#' @param z Cluster membership indicator vector of length `N`.
#' @param Theta Unobserved data matrix of dimension `N x q`.
#' @param mu Cluster means matrix of dimension `K x q`.
#' @param q Size of \eqn{\theta_i} (positive integer).
#' @param K Number of clusters (positive integer).
#' @param v_k Baseline degree of freedom for Sigma. Either a vector of length `K` or a scalar.
#' @param omega Baseline scale matrix for Sigma. Either a `q x q` matrix or a `q x q x K` array.
#' @param pars Return the scale parameter and degree of freedom.
#'
#' @return A `q x q x K` array of \eqn{\Sigma} drawn from its conditional posterior distribution.
#' If pars is TRUE, a list of parameters of the conditional posterior distribution is returned.
#'
#' @details The conditional posterior distribution of \eqn{\Sigma} is given by:
#' ```
#' \eqn{\Sigma_k | A - \{ \Sigma_1,...,\Sigma_K \} \sim InvWish(\Omega_k + 1[z_i=k] (\theta_i - \mu_k)(\theta_i - \mu_k)', N_k + v_k)}.
#' ```
#' Hence `Sigma` is sampled from the Inverse-Wishart distribution given all other parameters.
update_Sigma <- function(z, Theta, mu, q, K, v_k, omega, pars = FALSE) {
  N_k <- freq_table(z, K)
  if (length(dim(omega)) == 2) {omega <- replicate(K, omega)} # make into 3D array
  # list of scale matrix for the Inverse-Wishart distribution
  l <- lapply(1:K, function(k) {
    if (N_k[k] > 0) {
      Theta_mu <- Theta[z==k,] - mu[rep(k, N_k[k]),]
      # special case where N_k = 1 and must use `tcrossprod` instead of `crossprod`
      if (N_k[k] == 1) {cp <- tcrossprod(Theta_mu)}
      else {cp <- crossprod(Theta_mu)}
    }
    else {cp <- 0}
    return(cp + omega[,,k])
  })
  if (pars) {return(list(scale = l, df = N_k + v_k))}
  riwish(K, array(unlist(l), dim = c(q,q,K)), N_k + v_k)
}


#' Conditional update for Gaussian Mixture model weight \eqn{\rho}.
#'
#' @param z Cluster membership indicator vector of length `N`.
#' @param K Number of clusters (positive integer).
#' @param pars Return the vector length `K` of number of observations.
#'
#' @return A vector of length `K` of \eqn{\rho} drawn from its conditional posterior distribution.
#' If pars is TRUE, a list of parameters of the conditional posterior distribution is returned.
#'
#' @details The conditional posterior distribution of \eqn{\rho} is given by:
#' ```
#' \eqn{\rho | A - \{ \rho \} \sim Dirichlet(\alpha)},    \eqn{\alpha = (N_1 + 1 ,..., N_K + 1)}.
#' ```
#' Hence `rho` is sampled from the Dirichlet distribution given all other parameters.
update_rho <- function(z, K, pars = FALSE) {
  N_k <- freq_table(z, K)
  alpha <- N_k + 1
  if (pars) {return(list(alpha = alpha))}
  rdirichlet(1, alpha)
}


#' Conditional update for cluster membership indicator \eqn{z}.
#'
#' @param Sigma Cluster variance array of dimension `q x q x K`.
#' @param Theta Unobserved data matrix of dimension `N x q.`
#' @param mu Cluster means matrix of dimension `K x q`.
#' @param rho Gaussian Mixture model weight of vector of length `K`.
#' @param N Number of observations (positive integer).
#' @param q Size of \eqn{\theta_i} (positive integer).
#' @param K Number of clusters (positive integer).
#'
#' @return A list consisting of a vector of length `N` for \eqn{z}, and a `N x K` matrix for \eqn{lambda} (see details).
#'
#' @details The conditional posterior distribution of \eqn{z} is given by:
#' ```
#' \eqn{\kappa_{ik} = log \rho_k + log \phi(\theta_i | \mu_k, \Sigma_k) + CONST},
#' \eqn{z_i | A - \{ z \} \sim Multinomial(K, \lambda_i)},    \eqn{\lambda_{ik} = exp(\kappa_{ik}) / \sum^K_{m=1} exp(\kappa_{im}))}.
#' ```
#' Hence, `z` is sampled from the Multinomial distribution given all other parameters.
update_z <- function(Sigma, Theta, mu, rho, N, q, K) {
  exp_kappa <- sapply(1:K, function(k) {
    S <- Sigma[,,k]
    inS <- solveV(S)
    Theta_mu <- Theta - mu[rep(k, N),]
    # log likelihood of z
    rho[k] * exp(-0.5*(diag(Theta_mu %*% inS %*% t(Theta_mu)) + log(det(S))))
  })
  lambda <- exp_kappa / rowSums(exp_kappa) # normalize by row
  new_z <- rcategorical(t(lambda))
  list(updated_z = new_z, lambda = lambda)
}
