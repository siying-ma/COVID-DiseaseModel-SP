# runs a MCMC on a data list or the simulated history from make_history()
build_nimble_model <- function(data, chains, iter, burnin) {
  
  covidCode <- nimbleCode({
    # uniform priors
    psi ~ dunif(0, 1)
    omega ~ dunif(0, 1)
    theta ~ dbeta(1, 1)
    pa ~ dbeta(1, 1)
    pb ~ dbeta(1, 1)
    
    # model
    for (i in 1:M) {
      z[i] ~ dbern(psi)
      z_b[i] ~ dbern(z[i] * omega)
      lag[i] ~ dnbinom(size = 1, prob = theta)
      lag_severe[i] <- min(lag[i]+1, k)
      b[i,1:k] <- c(z[i] * z_b[i] * c(0, rep(0, lag_severe[i]-1), rep(1, k-lag_severe[i])))
      
      #' j=1
      h_hosp[i,1] ~ dbern(0)
      h_lab[i,1] ~ dbern(z[i] * pa)
      
      #' j>=2
      for (j in 2:k) {
        h_hosp[i,j] ~ dbern(z[i] * z_b[i] * b[i,j] * prod(1-h_hosp[i,1:(j-1)]) * (1-pb)^(j-lag[i]-2) * pb)
        h_lab[i,j] ~ dbern(z[i]  * prod(1-h_lab[i,1:(j-1)]) * prod(1-h_hosp[i,1:j]) * (1-pa)^(j-1) * pa)
      }
    }
    
    N <- sum(z[1:M])
  })
  
  covidConsts <- list(M = data$M,
                      k = data$k)

  covidData <- list(h_lab  = data$h_lab_aug,
                    h_hosp = data$h_hosp_aug,
                    z      = data$z_init_1,
                    z_b    = data$z_b_init_1, 
                    lag    = data$lag_init_1) 

  covidInits <- list(psi     = runif(1, data$psi_init, 1),
                     omega  = runif(1, data$omega_init, 1),
                     theta   = data$theta_init,
                     pa      = data$pa_init,
                     pb      = data$pb_init,
                     z       = data$z_init_2,
                     z_b     = data$z_b_init_2,
                     lag     = data$lag_init_2)
  
  
  covidModel <- nimbleModel(code = covidCode,
                            constants = covidConsts,
                            data = covidData,
                            inits = covidInits)
  
  covidConf <- configureMCMC(model = covidModel,
                             monitors  = c("N", "psi", "omega", "theta", "pa", "pb"),
                             enableWAIC = TRUE)
  
  covidMCMC <- buildMCMC(conf = covidConf)
  
  compModel <- compileNimble(covidModel)
  
  compMCMC <- compileNimble(covidMCMC, project = compModel)
  
  mcmc_out <- runMCMC(mcmc = compMCMC,
                      niter = iter,
                      nburnin = burnin,
                      nchains = chains,
                      summary = TRUE,
                      WAIC = TRUE)
  
  return(mcmc_out)
}
