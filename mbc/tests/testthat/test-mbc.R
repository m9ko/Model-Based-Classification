# Contents of file `test-mbc.R`

test_that("The normalized conditional log-PDF and the unnormalized
          complete data log-posterior should differ by the the same
          constant for any value of the argument with all the conditions held fixed", {
            tol <- 1e-5 # tolerance level
            
            # randomly define parameter and test size
            N <- sample(100:500, 1)
            q <- sample(2:10, 1)
            K <- sample(2:10, 1)
            testsize <- sample(5:10, 1)
            
            # randomly generate data
            rho <- rdirichlet(1, sample(1:10, K))
            z <- sample(1:K, N, TRUE, rho)
            mu <- matrix(rnorm(q*K), ncol = q)
            Sigma <- replicate(K, expr = {rwish(1, diag(q), q+1)})
            V <- replicate(N, expr = {rwish(1, diag(q), q+1)})
            Theta <- rmNorm(N, mu[z,], Sigma[,,z])
            y <- rmNorm(N, Theta, V)
            Var_y <- cov(y)
            v_k <- q+2

            ## Test Theta ------------------------------------------------------------------

            # Theta test cases
            test_Theta <- replicate(testsize, update_Theta(y, V, z, mu, Sigma, N))
            # parameters of the conditional posterior distribution
            Theta_params <- update_Theta(y, V, z, mu, Sigma, N, pars = TRUE)
            
            # normalized conditional log-PDF
            test1 <- sapply(1:testsize, function(i) {
              sum(dmNorm(test_Theta[,,i], mu = Theta_params$mean_Theta,
                         Sigma = Theta_params$var_Theta, log = TRUE))
            })
            
            # complete conditional log-PDF
            test2 <- sapply(1:testsize, function(i) {
              Theta_i <- test_Theta[,,i]
              ll(y, Theta_i, V, Sigma, mu, z, rho, N, q, K, v_k, Var_y)
            })
            
            # the two calculations should be off by the same constant
            expect_equal(var(test1 - test2), 0, tolerance = tol)
            
            ## Test mu ------------------------------------------------------------------

            # mu test cases
            test_mu <- replicate(testsize, update_mu(z, mu, Theta, Sigma, N, q, K))
            # parameters of the conditional posterior distribution
            mu_params <- update_mu(z, mu, Theta, Sigma, N, q, K, pars = TRUE)

            # normalized conditional log-PDF
            test1 <- sapply(1:testsize, function(i) {
              N_k <- freq_table(z, K)
              nonzero <- which(N_k > 0)
              sum(dmNorm(test_mu[nonzero,,i],
                         mu = mu_params$mean_mu[nonzero,],
                         Sigma = mu_params$var_mu[,,nonzero], log = TRUE))
            })

            # complete conditional log-PDF
            test2 <- sapply(1:testsize, function(i) {
              mu_i <- test_mu[,,i]
              ll(y, Theta, V, Sigma, mu_i, z, rho, N, q, K, v_k, Var_y)
            })

            # the two calculations should be off by the same constant
            expect_equal(var(test1 - test2), 0, tolerance = tol)
            
            ## Test Sigma ------------------------------------------------------------------

            # Sigma test cases
            test_Sigma <- replicate(testsize, update_Sigma(z, Theta, mu, q, K, v_k, Var_y))
            # parameters of the conditional posterior distribution
            Sigma_params <- update_Sigma(z, Theta, mu, q, K, v_k, Var_y, pars = TRUE)

            # normalized conditional log-PDF
            test1 <- sapply(1:testsize, function(i) {
              sum(sapply(1:K, function(k) {
                diwish(test_Sigma[,,k,i], Sigma_params$scale[[k]], Sigma_params$df[k], log = TRUE)
              }))
            })

            # complete conditional log-PDF
            test2 <- sapply(1:testsize, function(i) {
              Sigma_i <- test_Sigma[,,,i]
              ll(y, Theta, V, Sigma_i, mu, z, rho, N, q, K, v_k, Var_y)
            })

            # the two calculations should be off by the same constant
            expect_equal(var(test1 - test2), 0, tolerance = tol)
            
            ## Update rho ------------------------------------------------------------------

            # rho test cases
            test_rho <- replicate(testsize, update_rho(z, K))
            # parameters of the conditional posterior distribution
            rho_params <- update_rho(z, K, pars = TRUE)

            # normalized conditional log-PDF
            test1 <- sapply(1:testsize, function(i) {
              sum(ddirichlet(test_rho[,i], rho_params$alpha, log = TRUE))
            })

            # complete conditional log-PDF
            test2 <- sapply(1:testsize, function(i) {
              rho_i <- test_rho[,i]
              ll(y, Theta, V, Sigma, mu, z, rho_i, N, q, K, v_k, Var_y)
            })

            # the two calculations should be off by the same constant
            expect_equal(var(test1 - test2), 0, tolerance = tol)
            
            ## Update z ------------------------------------------------------------------

            # z test cases
            test_z <- replicate(testsize, update_z(Sigma, Theta, mu, rho, N, q, K)$updated_z)
            # parameters of the conditional posterior distribution
            z_params <- update_z(Sigma, Theta, mu, rho, N, q, K)$lambda

            # normalized conditional log-PDF
            test1 <- sapply(1:testsize, function(i) {
              sum(sapply(1:N, function(j) {
                log(dmultinom(as.numeric(1:K == test_z[j,i]), prob = z_params[j,]))
              }))
            })

            # complete conditional log-PDF
            test2 <- sapply(1:testsize, function(i) {
              z_i <- test_z[,i]
              ll(y, Theta, V, Sigma, mu, z_i, rho, N, q, K, v_k, Var_y)
            })

            # the two calculations should be off by the same constant
            expect_equal(var(test1 - test2), 0, tolerance = tol)
            
})
