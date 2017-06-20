sink(file="inst/cor_N_static.jags")
cat("
    model {

    #### Prior distributions
    beta0 ~ dnorm(0, 0.1)T(-10, 10)
    beta1 ~ dnorm(0, 0.1)T(-10, 10)
    tau ~ dgamma(0.1, 0.1)
    sigma <- sqrt(1/tau)

    alpha0 ~ dnorm(0, 0.1)T(-10, 10)
    alpha1 ~ dnorm(0, 0.1)T(-10, 10)
    alpha2 ~ dnorm(0, 0.1)T(-10, 10)

    xpsi[1] ~ dunif(0, 1)
    xpsi[2] ~ dunif(0, 1)


    ## Detection probability
    for(jj in 1:nStops){
      logit(p[jj]) <- alpha0 + alpha1 * stop[jj] + alpha2 * stop2[jj]
    }

    for (ii in 1:nRoutes) {

      ## Expected total count
      eta[ii] ~ dnorm(0, tau)
      log(lambda[ii]) <- beta0 + beta1 * X[ii] + eta[ii]


      ##  y: local availability at stop 1 -- 0 = locally unavailable,  1 = locally available
      y[ii, 1] ~ dbern(xpsi[1]/(xpsi[1] + (1 - xpsi[2])))


      ## Availability at stops 2-nStops
      for (jj in 2:nStops) {
        y[ii, jj] ~ dbern(xpsi[(y[ii, jj - 1] + 1)])
      } # jj


      ## Total number of available stops
      ns[ii] <- sum(y[ii,]) + C


      ## Multinomial cell probabilities
      pi[ii, 1] <- p[1] * (y[ii, 1]/ns[ii]) + C
      pic[ii, 1] <- pi[ii, 1] / pcap[ii]

      for(jj in 2:nStops){
        pi[ii, jj] <- p[jj] * (y[ii, jj]/ns[ii])
        pic[ii, jj] <- pi[ii, jj] / pcap[ii]
      }


      ## Total 'capture' probability
      pcap[ii] <- sum(pi[ii,])


      ## Multinomial likelihood
      h[ii, 1:nStops] ~ dmulti(pic[ii, 1:nStops], n[ii])
      n[ii] ~ dbin(pcap[ii], N[ii])
      N[ii] ~ dpois(lambda[ii])
    } # ii

    }
    ", fill=TRUE)
sink()
