#' run_static
#'
#' Run static correlated-detection N-mixture model
#' @param sim_data Simulated data from sim_data_static()
#' @param nC Number of chains
#' @param nI Number of iterations
#' @param nB Number of iterations to discard burn-in
#' @param nT Thinning rate
#' @param Parallel Should chains be run in parallel (default = TRUE)
#' @export

run_static <- function(sim_data, nC = 3, nI = 30000, nB = 10000, nT = 20, Parallel = TRUE){

  jags.data <- list(h = sim_data$h, nRoutes = dim(sim_data$h)[1], nStops = dim(sim_data$h)[2],
                    n = sim_data$n, stop = sim_data$stop, stop2 = sim_data$stop2,
                    X = sim_data$X, C = 10^6)

  jags.param <- c("xpsi", "beta0", "beta1", "sigma", "alpha0", "alpha1", "alpha2", "y", "ns", "p", "lambda", "N")

  jags.inits <- function(){list(N = sim_data$n, y = sim_data$y,
                                alpha0 = sim_data$alpha0, alpha1 = sim_data$alpha1,
                                alpha2 = sim_data$alpha2,
                                beta0 = sim_data$beta0,
                                xpsi = sim_data$theta, tau = 1/(sim_data$sigma))}

  jags.fit <- jagsUI::jags(data = jags.data, parameters.to.save = jags.param, inits = jags.inits,
                          model.file = "inst/cor_N_static.jags", n.chains = nC, n.iter = nI,
                          n.burnin = nB, n.thin = nT, parallel = Parallel)

  return(jags.fit)
}
