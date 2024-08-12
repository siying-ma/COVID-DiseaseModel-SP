#' Simulate histories for h_lab and h_hosp
#'
#' @param N population size
#' @param M augmented population size
#' @param k number of capture occasions
#' @param omaga_a probability of being symptomatic
#' @param pmega_b probability of showing severe symptom
#' @param theta probability of lag of severe symptom
#' @param pa capture probability of lab test
#' @param pb capture probability of hospitalization
#'
#' `b` latent indicator of first time that individual shows severe symptoms
#' `h_lab` capture history of lab test (patient can be captured via lab test any time after showing symptom)
#' `h_hosp` capture history of hospitalization (patient can be captured via hospitalization any time after showing severe symptom)
#' 
make_history <- function(N, M, k, omega, theta, pa, pb) {
  
  z_b <- rep(NA,N)
  lag <- rep(NA,N)
  lag_severe <-rep(NA,N)
  b <- matrix(NA, nrow = N, ncol = k)
  h_lab <- matrix(NA, nrow = N, ncol = k) 
  h_hosp <- matrix(NA, nrow = N, ncol = k) # to contain 2s for recaptures
  
  #' We use onset symptomatic time: all the individual are showing symptoms of COVID on Time 1
  #' First, simulate latent state `b` and use it to create `h_lab` and `h_hosp`
  for (i in 1:N) {
    #' `z_b` is the probability where a patient has severe symptom
    z_b[i] <- rbern(1, omega)
    
    #' The lag of severe symptom is following the geometric distribution starting from 0
    #' Since NIMBLE does not support `dgeom`, so we use negative binomial with size equal to 1 instead
    lag[i] <- rnbinom(1, size = 1, prob = theta)
    
    #' The actual lag of severe symptom starts from 1 (assumption) and then we restrict the lag within Time k
    lag_severe[i] <- min(lag[i] + 1, k)
    b[i,] <- c(z_b[i] * c(0, rep(0, lag_severe[i]-1), rep(1, k-lag_severe[i])))
    
    #' Generate observation matrix
    for (j in 1:k) {
      if (j == 1) {
        h_lab[i,j] <- rbern(1, pa)
        h_hosp[i,j] <- 0
      } else {
        #' must populate `h_hosp` first: 
        #' if available for both lab test and hospital in the same time period then only count `h_hosp` (Assumption)
        #' if a patient is captured from hospitalization, then he is no longer available for lab test (Assumption)
        h_hosp[i,j] <- rbern(1, b[i,j] * prod(1-h_hosp[i,1:(j-1)]) * (1-pb)^(j-lag[i]-2) * pb)
        h_lab[i,j] <- rbern(1, prod(1-h_lab[i,1:(j-1)]) * prod(1-h_hosp[i,1:j]) * (1-pa)^(j-1) * pa)
      }     
    }
  }

  #' augmented histories contain M-N 0s for h-lab and h_hosp
  h_lab_aug   <- rbind(h_lab, matrix(0, nrow = M - N, ncol = k))
  h_hosp_aug  <- rbind(h_hosp, matrix(0, nrow = M - N, ncol = k))
  
  #' init value of `z`, `z_a`, `z_b` and `lag` (partially observed)
  #' init value of `psi`, `theta`, `omega`, `pa` and `pb` (parameters)
  z_init <- rowSums(h_lab_aug + h_hosp_aug)
  z_init_1 <- z_init_2 <- z_init
  z_init_1[z_init == 0] <- NA
  z_init_1[z_init > 1]  <- 1
  z_init_2[z_init >= 1]  <- NA
  
  #' The init value of `z_a` is same as the init value as `z`
  #' Since if a patient is observed in hospital then he has to be symptomatic 
  #' (Assumption according to stratified Peterson estimator)
  
  z_b_init   <- rowSums(h_hosp_aug)
  z_b_init_1 <- z_b_init_2 <- z_b_init
  z_b_init_1[z_b_init==0] <- NA
  z_b_init_2[z_b_init==1] <- NA
  
  lag_obs <- rowSums(h_hosp_aug)
  lag_init_1 <- c(lag, rep(NA,M-N))
  lag_init_1[lag_obs==0] <- NA
  
  mean_lag <- mean(lag_init_1, na.rm = T)
  theta_init <- 1 / (mean_lag + 1)
  
  lag_init_2 <- rnbinom(M, size = 1, prob = theta_init)
  lag_init_2[lag_obs==1] <- NA
  
  n <- sum(z_init_1 == 1, na.rm = TRUE) # number captured
  print(paste("n =", n))
  n_2 <- sum(z_b_init_1 == 1, na.rm = T) # number recaptured
  print(paste("n_2 =", n_2))
  print(paste("h_lab =", sum(h_lab == 1)))
  print(paste("h_hosp =", sum(h_hosp == 1, na.rm = T)))
  
  psi_init <- n / M
  omega_init <- n_2 / M
  
  lag_init_lab <- apply(h_lab_aug, 1, function(x) ifelse(sum(x)==0, NA, which(x==1)))
  mean_lag_lab <- mean(lag_init_lab, na.rm = T)
  pa_init <- 1 / mean_lag_lab
  
  lag_init_hosp <- apply(h_hosp_aug, 1, function(x) ifelse(sum(x)==0, NA, which(x==1)))
  mean_lag_hosp <- mean((lag_init_hosp - lag_init_1 - 1), na.rm = T)
  pb_init <- 1 / mean_lag_hosp
  
  
  #' return list `result`
  result              <- NULL
  result$N            <- N
  result$n            <- n
  result$n_2          <- n_2
  result$M            <- M
  result$k            <- k
  result$omega.       <- omega
  result$pa           <- pa
  result$pb           <- pb
  result$theta        <- theta
  result$h_lab        <- h_lab
  result$h_hosp       <- h_hosp
  result$h_lab_aug    <- h_lab_aug
  result$h_hosp_aug   <- h_hosp_aug
  result$b            <- b
  result$z_b          <- z_b
  result$z_init_1     <- z_init_1
  result$z_init_2     <- z_init_2
  result$z_b_init_1   <- z_b_init_1
  result$z_b_init_2   <- z_b_init_2
  result$lag          <- lag
  result$lag_severe   <- lag_severe
  result$lag_init_1   <- lag_init_1
  result$lag_init_2   <- lag_init_2
  result$psi_init     <- psi_init
  result$omega_init   <- omega_init
  result$theta_init   <- theta_init
  result$pa_init      <- pa_init
  result$pb_init      <- pb_init

  return(result)
}

