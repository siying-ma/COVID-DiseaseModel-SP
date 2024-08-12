#' Simulate histories h_lab and h_hosp
#'
#' @param N population size
#' @param M augmented population size
#' @param k number of capture occasions
#' @param theta probability of lag of severe symptom
#' @param pa capture probability of lab test
#' @param pb capture probability of hospitalization
#'
#' `b` latent indicator of first time that individual shows severe symptoms
#' `h_lab` capture history of lab test (patient can be captured via lab test any time after showing symptom)
#' `h_hosp` capture history of hospitalization (patient can be captured via hospitalization any time after showing severe symptom)
#' 
make_history <- function(N, M, k, omega, pa, pb) {
  
  z_b <- rep(NA,N)
  #lag <- rep(NA,N)
  #lag_severe <-rep(NA,N)
  #b <- matrix(NA, nrow = N, ncol = k)
  h_lab <- matrix(NA, nrow = N, ncol = k) 
  h_hosp <- matrix(NA, nrow = N, ncol = k) # to contain 2s for recaptures
  
  #' We use onset symptomatic time: all the individual are showing symptoms of COVID on Time 1
  #' First, simulate latent state `b` and use it to create `h_lab` and `h_hosp`
  for (i in 1:N) {
    #' `z_b` is the probability where a patient has severe symptom
    z_b[i] <- rbern(1, omega)
    
    #' Generate observation matrix
    for (j in 1:k) {
      if (j == 1) {
        h_lab[i,j] <- rbern(1, pa)
        h_hosp[i,j] <- 0
      } else {
        #' must populate `h_hosp` first: 
        #' if available for both lab test and hospital in the same time period then only count `h_hosp` (Assumption)
        #' if a patient is captured from hospitalization, then he is no longer available for lab test (Assumption)
        h_hosp[i,j] <- rbern(1, z_b[i] * prod(1-h_hosp[i,1:(j-1)]) * (1-pb)^(j-2) * pb)
        h_lab[i,j] <- rbern(1, prod(1-h_lab[i,1:(j-1)]) * prod(1-h_hosp[i,1:j]) * (1-pa)^(j-1) * pa)
      }     
    }
  }
  
  case_wo_sym <- sample(1:nrow(h_lab), size = ceiling(0.3*nrow(h_lab)),replace = F)
  h_lab_1 <- h_lab[-case_wo_sym,]
  h_lab_2 <- h_lab[case_wo_sym,]
  h_hosp_1 <- h_hosp[-case_wo_sym,]
  h_hosp_2 <- h_hosp[case_wo_sym,]
  
  h_lab_asym_pos <- rowSums(h_lab_2)==1
  h_hosp_asym_pos <- rowSums(h_hosp_2) ==1
  
  h_asym_lab_hosp <- h_lab_asym_pos & h_hosp_asym_pos
  h_asym_hosp_only <- !h_lab_asym_pos & h_hosp_asym_pos
  h_asym_lab_only <- h_lab_asym_pos & !h_hosp_asym_pos
  h_asym_not_obs <- !h_lab_asym_pos & !h_hosp_asym_pos

  sum(h_asym_lab_hosp)
  sum(h_asym_lab_only)
  sum(h_asym_hosp_only)
  sum(h_asym_not_obs)
  
  h_lab_asym_lab_hosp <- h_lab_2[h_asym_lab_hosp,]
  h_hosp_asym_lab_hosp <- h_hosp_2[h_asym_lab_hosp,]
  h_gap_lab_hosp <- apply(h_hosp_asym_lab_hosp, 1, function(x) which(x==1)) -
    apply(h_lab_asym_lab_hosp, 1, function(x) which(x==1))
  h_lab_asym_lab_hosp <- matrix(NA, nrow = nrow(h_lab_asym_lab_hosp), ncol = k)
  h_hosp_asym_lab_hosp <- matrix(NA, nrow = nrow(h_lab_asym_lab_hosp), ncol = k)
  for (i in 1:nrow(h_hosp_asym_lab_hosp)) {
    h_hosp_asym_lab_hosp[i, 1:h_gap_lab_hosp[i]] <- rep(0, h_gap_lab_hosp[i])
  }
  
  h_lab_asym_hosp_only <- h_lab_2[h_asym_hosp_only,]
  h_hosp_asym_hosp_only <- h_hosp_2[h_asym_hosp_only,]
  h_lab_asym_hosp_only <- matrix(0, nrow = nrow(h_lab_asym_hosp_only), ncol = k)
  h_hosp_asym_hosp_only <- matrix(NA, nrow = nrow(h_hosp_asym_hosp_only), ncol = k)
  
  h_lab_asym_lab_only <- h_lab_2[h_asym_lab_only,]
  h_hosp_asym_lab_only <- h_hosp_2[h_asym_lab_only,]
  h_lab_asym_lab_only <- matrix(NA, nrow = nrow(h_lab_asym_lab_only), ncol = k)
  h_hosp_asym_lab_only <- matrix(0, nrow = nrow(h_hosp_asym_lab_only), ncol = k)
  
  h_lab_asym_not_obs<- h_lab_2[h_asym_not_obs,]
  h_hosp_asym_not_obs <- h_hosp_2[h_asym_not_obs,]
  
  h_lab_asym <- rbind(h_lab_asym_hosp_only,
                      h_lab_asym_lab_hosp,
                      h_lab_asym_lab_only)
  
  h_hosp_asym <- rbind(h_hosp_asym_hosp_only,
                       h_hosp_asym_lab_hosp,
                       h_hosp_asym_lab_only)
  
  #' augmented histories contain more 0s for h-lab and h_hosp
  h_lab_aug   <- rbind(h_lab_1, h_lab_asym_not_obs, matrix(0, nrow = M - N, ncol = k))
  h_hosp_aug  <- rbind(h_hosp_1, h_hosp_asym_not_obs, matrix(0, nrow = M - N, ncol = k))
  
  #' init value for `z`, `z_b` (partially observed) and `lag` (fully latent)
  z_b_init   <- rowSums(h_hosp_aug)
  z_b_init_1 <- z_b_init_2 <- z_b_init
  z_b_init_1[z_b_init==0] <- NA
  z_b_init_2[z_b_init==1] <- NA
  
  z_b_asym_init_1 <- c(rep(1, nrow(h_lab_asym_lab_hosp) + nrow(h_lab_asym_hosp_only)), rep(NA, nrow(h_lab_asym_lab_only)))
  z_b_asym_init_2 <- c(rep(NA, nrow(h_lab_asym_lab_hosp) + nrow(h_lab_asym_hosp_only)), rep(0, nrow(h_lab_asym_lab_only)))
  z_b_init_1 <- c(z_b_init_1, z_b_asym_init_1)
  z_b_init_2 <- c(z_b_init_2, z_b_asym_init_2)
  
  z_init <- rowSums(h_lab_aug + h_hosp_aug)
  z_init_1 <- z_init_2 <- z_init
  z_init_1[z_init == 0] <- NA
  z_init_1[z_init > 1]  <- 1
  z_init_2[z_init >= 1]  <- NA
  
  z_asym_init_1 <- c(rep(1, nrow(h_lab_asym)))
  z_asym_init_2 <- c(rep(NA, nrow(h_lab_asym)))
  z_init_1 <- c(z_init_1, z_asym_init_1)
  z_init_2 <- c(z_init_2, z_asym_init_2)
  
  n <- sum(z_init_1 == 1, na.rm = TRUE) # number captured
  print(paste("n =", n))
  n_2 <- sum(z_b_init_1 == 1, na.rm = T) # number recaptured
  print(paste("n_2 =", n_2))
  
  #' init value for `psi`, `theta`, `omega`, `pa` and `pb`
  psi_init <- n / M
  omega_init <- n_2 / M
  
  lag_init_hosp <- apply(h_hosp_aug, 1, function(x) ifelse(sum(x)==0, NA, which(x==1)))
  lag_init_hosp <- lag_init_hosp-1
  mean_lag_hosp <- mean(lag_init_hosp, na.rm = T)
  pb_init <- 1 / mean_lag_hosp
  
  lag_init_lab <- apply(h_lab_aug, 1, function(x) ifelse(sum(x)==0, NA, which(x==1)))
  mean_lag_lab <- mean(lag_init_lab, na.rm = T)
  pa_init <- 1 / mean_lag_lab
  
  #####Create inits for asymptomatic patients
  h_lab_asym_both_lab_hosp_init <- matrix(0, nrow = nrow(h_lab_asym_lab_hosp), ncol = k)
  h_hosp_asym_both_lab_hosp_init <- matrix(0, nrow = nrow(h_lab_asym_lab_hosp), ncol = k)
  for (i in 1:nrow(h_lab_asym_lab_hosp)) {
    h_hosp_asym_both_lab_hosp_init[i, 1:h_gap_lab_hosp[i]] <- rep(NA,h_gap_lab_hosp[i])
  }
  for (i in 1:nrow(h_lab_asym_both_lab_hosp_init)) {
    asym_lab_hosp <- rnbinom(1, size = 1, prob = pa_init)
    while (asym_lab_hosp+1 > k) {
      asym_lab_hosp <- rnbinom(1, size = 1, prob = pa_init)
    }
    h_lab_asym_both_lab_hosp_init[i, asym_lab_hosp+1] <- 1
    lab_hosp <- asym_lab_hosp+ 1 + h_gap_lab_hosp[i]
    if (lab_hosp <= k){
      h_hosp_asym_both_lab_hosp_init[i, lab_hosp] <- 1
    }
  }
  
  
  h_hosp_asym_hosp_init <- matrix(0, nrow = nrow(h_hosp_asym_hosp_only), ncol = k)
  for (i in 1:nrow(h_hosp_asym_hosp_only)) {
    asym_hosp <- rnbinom(1, size = 1, prob = pb_init)
    while (asym_hosp+2 > k) {
      asym_hosp <- rnbinom(1, size = 1, prob = pb_init)
    }
    h_hosp_asym_hosp_init[i, asym_hosp+2] <- 1
  }
  h_hosp_asym_hosp_init[,1] <- rep(NA, nrow(h_hosp_asym_hosp_init))
  
  
  h_lab_asym_only_lab_init <- matrix(0, nrow = nrow(h_lab_asym_lab_only), ncol = k)
  for (i in 1:nrow(h_lab_asym_lab_only)) {
    asym_lab <- rnbinom(1, size = 1, prob = pa_init)
    while(asym_lab+1 > k){
      asym_lab <- rnbinom(1, size = 1, prob = pa_init)
    }
    h_lab_asym_only_lab_init[i, asym_lab+1] <- 1
  }
  
  h_lab_asym_init <- rbind(matrix(NA, nrow = nrow(h_hosp_asym_hosp_init), ncol = k),
                           h_lab_asym_both_lab_hosp_init,
                           h_lab_asym_only_lab_init)
  h_hosp_asym_init <- rbind(h_hosp_asym_hosp_init,
                            h_hosp_asym_both_lab_hosp_init,
                            matrix(NA, nrow = nrow(h_lab_asym_only_lab_init), ncol = k))
  
  h_hosp_asym_start <- nrow(h_hosp_aug) + 1
  h_hosp_asym_end <- nrow(h_hosp_aug) + nrow(h_hosp_asym_both_lab_hosp_init) + nrow(h_hosp_asym_hosp_init)
  h_lab_asym_start <- nrow(h_lab_aug) + nrow(h_hosp_asym_hosp_init) + 1
  h_lab_asym_end <- M
  
  h_lab_aug_init <- matrix(NA, nrow = nrow(h_lab_aug), ncol = k)
  h_hosp_aug_init <- matrix(NA, nrow = nrow(h_hosp_aug), ncol = k)
  h_lab_aug_init <- rbind(h_lab_aug_init, h_lab_asym_init)
  h_hosp_aug_init <- rbind(h_hosp_aug_init, h_hosp_asym_init)
  dim(h_lab_aug_init)
  dim(h_hosp_aug_init)
  
  h_lab_aug <- rbind(h_lab_aug, h_lab_asym)
  h_hosp_aug <- rbind(h_hosp_aug, h_hosp_asym)
  dim(h_lab_aug)
  dim(h_hosp_aug)
  
  
  result <- NULL
  
  result$N            <- N
  result$n            <- n
  result$n_2          <- n_2
  result$M            <- M
  result$k            <- k
  result$omega        <- omega
  result$pa           <- pa
  result$pb           <- pb
  result$h_lab        <- h_lab
  result$h_lab_asym   <- h_lab_asym
  result$h_hosp       <- h_hosp
  result$h_hosp_asym  <- h_hosp_asym
  result$h_lab_aug    <- h_lab_aug
  result$h_lab_aug_init  <- h_lab_aug_init
  result$h_hosp_aug   <- h_hosp_aug
  result$h_hosp_aug_init <- h_hosp_aug_init
  result$h_lab_asym_start  <- h_lab_asym_start
  result$h_lab_asym_end    <- h_lab_asym_end
  result$h_hosp_asym_start <- h_hosp_asym_start
  result$h_hosp_asym_end   <- h_hosp_asym_end
  result$z_init_1     <- z_init_1
  result$z_init_2     <- z_init_2
  result$z_b          <- z_b
  result$z_b_init     <- z_b_init
  result$z_b_init_1   <- z_b_init_1
  result$z_b_init_2   <- z_b_init_2
  result$psi_init     <- psi_init
  result$omega_init   <- omega_init
  result$pa_init      <- pa_init
  result$pb_init      <- pb_init

  

  return(result)
}

