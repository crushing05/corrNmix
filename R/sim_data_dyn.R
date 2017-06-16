#' sim_data_dyn
#'
#' Simulation of multi-season count data for the dynamic correlated-detection N-mixture model
#' @param nRoutes Number of sample units
#' @param nStops Number of stops within each sample unit
#' @param nYears Number of years
#' @param alpha0 Intercept for detection model
#' @param alpha1 Linear effect of stop number on detection probability
#' @param alpha2 Quadratic effect of stop number on detection probability
#' @param beta0 Log expected number of individuals at each sampling unit in year 1
#' @param sigma1 Standard deviation around mean count in year 1
#' @param gamma Recruitment rate
#' @param sigma2 Standard deviation of random immigration from other sies
#' @param omega Survival rate
#' @param theta Vector containing the correlation terms
#' @param nSum Optional number of stops to aggregate
#' @export

sim_data_dyn<- function(nRoutes = 150, nStops = 50, nYears = 10,
                            alpha0 = -0.5, alpha1 = -0.3, alpha2 = 0.2,
                            beta0 = 2, sigma1 = 1,
                            gamma = 0.35, omega = 0.6, sigma2 = 0.1,
                            theta = c(0.3, 0.75),
                            nSum = NULL){


  ## Standardize stop number
  stop <- scale(seq(1:nStops))[,1]
  stop2 <- stop^2

  ## Stop detection probability
  lp <- alpha0 + alpha1 * stop + alpha2 * stop2
  p <- exp(lp)/(1 + exp(lp))

  ## Route-level abundance in year 1
  eta <- rnorm(nRoutes, 0, sigma1)
  lN1 <- beta0 + eta
  N1 <- round(exp(lN1))

  ## Abundance in years 2-nYears
  N <- matrix(NA, nrow = nRoutes, ncol = nYears)
  N[, 1] <- N1

  epsilon <- matrix(NA, nrow = nRoutes, ncol = nYears - 1)
  for(i in 1:nRoutes){
    for(t in 2:nYears){
      epsilon[i, t - 1] <- exp(rnorm(1, 0, sigma2))
      N[i, t] <- rpois(1, gamma * N[i, t - 1] + epsilon[i, t - 1]) + rbinom(1, N[i, t - 1], omega)
    }
  }


  ## Stop-level availability & total available stops each year
  for(t in 1:nYears){
    y1 <- matrix(NA, nrow = nRoutes, ncol = nStops)
    for(i in 1:nRoutes){
      y1[i, 1] <- rbinom(1, 1, prob = theta[1]/(theta[1] + (1 - theta[2])))

      for(j in 2:nStops){
        y1[i, j] <- rbinom(1, 1, prob = theta[y1[i, j - 1] + 1])
      }
    }
    y.sum1 <- apply(y1, 1, sum)
    if(t == 1){
      y <- y1
      y.sum <- y.sum1
    }else{
        y <- abind::abind(y, y1, along = 3)
        y.sum <- cbind(y.sum, y.sum1)}
  }


  ## Multinomial cell probabilities
  for(t in 1:nYears){
    pi <- matrix(NA, nrow = nRoutes, ncol = nStops)
    h1 <- matrix(NA, nrow = nRoutes, ncol = nStops)

    for(i in 1:nRoutes){
      pi[i,] <- (p/y.sum[i, t])*y[i,,t]
    }

    ## Total detection probability
    pdet <- apply(pi, 1, sum)

    ## Total route-level count
    n1 <- rbinom(nRoutes, N[, t], pdet)

    ## Conditional cell probabilities
    pic <- pi/pdet

    for(i in 1:nRoutes){
      h1[i,] <- rmultinom(1, size = n1[i], prob = pic[i,])
    }

    if(t == 1){
      n <- n1
      h <- h1
    }else{
      n <- cbind(n, n1)
      h <- abind::abind(h, h1, along = 3)
    }
  }


  if(!is.null(nSum)){
    ## Sum stop-level counts
    for(t in 1:nYears){
      for(i in 1:nRoutes){
        ht <- unname(tapply(h[i,,t], (seq_along(h[i,,t])-1) %/% nSum, sum))
        if(i == 1){h2 <- ht}else{h2 <- rbind(h2, ht)}
      }
      if(t == 1){h3 <- h2}else{h3 <- abind::abind(h3, h2, along = 3)}
    }

    ## New number of stops
    nStops <- dim(h2)[2]

    ## New scaled stop number
    stop <- scale(seq(1:5))[,1]
    stop2 <- stop^2

    ## New stop-level availability
    y <- h3
    y[y > 0] <- 1
    h <- h3
  }

  sim_data <- list(N = N, n = n, h = h, y = y, p = p, nRoutes = nRoutes, nStops = nStops,
                   beta0 = beta0, sigma1 = sigma1, alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                   stop = stop, stop2 = stop2, gamma = gamma, sigma2 = sigma2, omega = omega, theta = theta)

  return(sim_data)
}
