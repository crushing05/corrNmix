#' sim_data_static
#'
#' Simulation of single-season count data for the static correlated-detection N-mixture model
#' @param nRoutes Number of sample units
#' @param nStops Number of stops within each sample unit
#' @param alpha0 Intercept for detection model
#' @param alpha1 Linear effect of stop number on detection probability
#' @param alpha2 Quadratic effect of stop number on detection probability
#' @param beta0 Log expected number of individuals at each sampling unit
#' @param sigma Standard deviation around mean count
#' @param theta Vector containing the correlation terms
#' @param nSum Optional number of stops to aggregate
#' @export

sim_data_static <- function(nRoutes = 150, nStops = 50,
                             alpha0 = -0.5, alpha1 = -0.3, alpha2 = 0.2,
                             beta0 = 1, beta1 = 1.5, sigma = 0.5,
                             theta = c(0.3, 0.75),
                             nSum = NULL){


  ## Standardize stop number
  stop <- scale(seq(1:nStops))[,1]
  stop2 <- stop^2

  ## Stop detection probability
  lp <- alpha0 + alpha1 * stop + alpha2 * stop2
  p <- exp(lp)/(1 + exp(lp))

  ## Route-level abundance
  eta <- rnorm(nRoutes, 0, sigma)
  X <- rnorm(nRoutes)
  l.lam <- beta0 + beta1 * X + eta
  lam <- exp(l.lam)

  ## Equilibrium proportion of available sites
  theta1 <- theta[1] / (theta[1] + (1 - theta[2]))

  ## Stop level data
  c <- matrix(NA, nrow = nRoutes, ncol = nStops)  # True number of indvs at each stop
  y <- matrix(NA, nrow = nRoutes, ncol = nStops)  # Stop availability
  h <- matrix(NA, nrow = nRoutes, ncol = nStops)  # Observation number of indvs

  for(i in 1:nRoutes){
    c[i, 1] <- rpois(n = 1, lambda = lam[i] * 1/(theta1 * nStops)) * rbernoulli(n = 1, p = theta1)
    y[i, 1] <- ifelse(c[i, 1] == 0, 0, 1)
    h[i, 1] <- rbinom(n = 1, size = c[i, 1], prob = p[1])
    for(j in 2:nStops){
      c[i, j] <- rpois(n = 1, lambda = lam[i] * 1/(theta1 * nStops)) * rbernoulli(n = 1, p = theta[y[i, j - 1] + 1])
      y[i, j] <- ifelse(c[i, j] == 0, 0, 1)
      h[i, j] <- rbinom(n = 1, size = c[i, j], prob = p[j])
    }
  }

  ## Total detection probability
  N <- apply(c, 1, sum)

  ## Total route-level count
  n <- apply(h, 1, sum)

  if(!is.null(nSum)){
    ## Sum stop-level counts
    for(i in 1:nRoutes){
      ht <- unname(tapply(h[i,], (seq_along(h[i,])-1) %/% nSum, sum))
      if(i == 1){h2 <- ht}else{h2 <- rbind(h2, ht)}
    }

    ## New number of stops
    nStops <- dim(h2)[2]

    ## New scaled stop number
    stop <- scale(seq(1:nStops))[,1]
    stop2 <- stop^2

    ## New stop-level availability
    y <- h2
    y[y > 0] <- 1
    h <- h2
  }

  sim_data <- list(N = N, n = n, h = h, y = y, p = p, X = X, nRoutes = nRoutes, nStops = nStops,
                   beta0 = beta0, sigma = sigma, alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                   stop = stop, stop2 = stop2, theta = theta)

  return(sim_data)
}
