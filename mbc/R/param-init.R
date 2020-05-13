# Contents of file `param-init.R`

#' Initial cluster allocation.
#'
#' Initializes the clusters using the kmeans++ algorithm of Arthur & Vassilvitskii (2007).
#'
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param K Number of clusters (positive integer).
#' @param max_iter The maximum number of steps in the [stats::kmeans()] algorithm.
#' @param nstart The number of random starts in the [stats::kmeans()] algorithm.
#' @return A vector of length `N` consisting of integers between `1` and `K` specifying an initial clustering of the observations.
#' @author Martin Lysy
#' @references Arthur, D., Vassilvitskii, S. "k-means++: the advantages of careful seeding" *Proceedings of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms.* (2007): 1027–1035. <http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf>.
init_z <- function(y, K, max_iter = 10, nstart = 10) {
  # init++
  N <- nrow(y)
  p <- ncol(y)
  x <- t(y) # easier to use columns
  centers <- matrix(NA, p, K) # centers
  icenters <- rep(NA, K-1) # indices of centers
  minD <- rep(Inf, N) # current minimum distance vector
  # initialize
  icenters[1] <- sample(N,1)
  centers[,1] <- x[,icenters[1]]
  for(ii in 1:(K-1)) {
    D <- colSums((x - centers[,ii])^2) # new distance
    minD <- pmin(minD, D)
    icenters[ii+1] <- sample(N, 1, prob = minD)
    centers[,ii+1] <- x[,icenters[ii+1]]
  }
  centers <- t(centers)
  colnames(centers) <- rownames(x)
  # run kmeans with these centers
  km <- kmeans(x = y, centers = centers, nstart = nstart, iter.max = max_iter)
  km$cluster
}


#' Parameter initialization routine.
#'
#' @param y Observed data of dimension `N x q`.
#' @param K Number of clusters (positive integer).
#' @param type The type of initialization method. Can be either "kmeans" or "random".
#' @param max_iter The maximum number of steps in the [stats::kmeans()] algorithm.
#' @param nstart The number of random starts in the [stats::kmeans()] algorithm.
#'
#' @return A list consisting of:
#' z0: initial estimate of the cluster membership indicator vector of length `N`,
#' \eqn{\rho}: initial estimate of Gaussian Mixture model weight of vector of length `K`,
#' \eqn{\Theta}: inital estimate of unoberved data of dimension `N x q`,
#' \eqn{\mu}: initial estimate of cluster means matrix of dimension `K x q`.
#' \eqn{\Sigma}: initial estimate of cluster variance array of dimension `q x q x K`.
#'
#' @details This functions uitlizes \code{init_z} from [nnm-functions].
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
#' init_params(y, K, type = "kmeans")
#'
#' @export
init_params <- function(y, K, type = "kmeans", max_iter = 10, nstart = 10) {
  N <- nrow(y)
  q <- ncol(y)
  # k-means algorithm for initial z
  if (type == "kmeans") {
    z0 <- init_z(y, K, max_iter, nstart)}
  # random assignment for initial z
  if (type == "random") {
    z0 <- sample(1:K, nrow(y), TRUE)}
  N_k <- freq_table(z0, K)
  rho <- N_k / N
  Theta <- y
  mu <- replicate(K, colMeans(y)) # baseline mu as overall sample mean
  Sigma <- replicate(K, var(y)) # baseline Sigma as overall sample variance
  for (k in 1:K) {
    # update mu if N_k > 0
    if (N_k[k] > 0) {
      obs_k <- y[z0==k,]
      if (N_k[k] == 1) {
        mu[,k] <- obs_k}
      else { 
        mu[,k] <- colMeans(obs_k)
        # update Sigma if N_k >= q
        if (N_k[k] >= q) {
          Sigma[,,k] <- var(obs_k)}
      }}
  }
  return(list(z0 = z0, rho = rho, Theta = Theta, mu = t(mu), Sigma = Sigma))
}
