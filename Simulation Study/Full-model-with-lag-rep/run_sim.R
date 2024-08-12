library(tictoc)
library(nimble)
library(dplyr)
library(Rlab)

#set working directionary
#setwd("/Users/siyingma/Desktop/COVID/COVID-codes/Full-model-with-lag")
source("./make_history_with_lag_partially_obs.R")
source("./build_nimble_model.R")

# R script for running a simulation study. Uses a parameter matrix
# created by generate_parameter_matrix.R. Takes the task id as a command
# line argument,which corresponds to a row of the parameter matrix.
# Outputs model_statistics_* text files with the summary statistics for the coda
# samples.

tic()

print("finished loading libraries")

args <- commandArgs(trailingOnly = TRUE)
task_id_str <- args[1]
print(paste0("task id: ", task_id_str))

# number of repetitions for each simulation
reps <- 5

params <- c(
  N = 2967,
  M = 3760,
  k = 28,
  theta = 0.15,
  omega = 0.1057,
  pa = 0.1128,
  pb = 0.4
)

# generate data

# true population size
N <- as.integer(params["N"])

# augmented population size
M <- as.integer(params["M"])

# number of sampling occasions / strata
k <- as.integer(params["k"])

# true severe symptoms availability probability
omega <- as.double(params["omega"])

# true lag of severe symptoms probability
theta <- as.double(params["theta"])

# true lab test capture probability
pa <- as.double(params["pa"])

# true hospital capture probability
pb <- as.double(params["pb"])

result_table <- matrix(NA, nrow = reps, ncol = 26)
colnames(result_table) <- c("Seed",
                            "N_mean", "N_omega", "N_sd", "N_lower", "N_upper",
                            "omega_mean", "omega_omega", "omega_sd", "omega_lower", "omega_upper",
                            "pa_mean", "pa_omega", "pa_sd", "pa_lower", "pa_upper",
                            "pb_mean", "pb_omega", "pb_sd", "pb_lower", "pb_upper",
                            "theta_mean", "theta_omega", "theta_sd", "theta_lower", "theta_upper")

for (rep in 1:reps) { # uncomment this if want repetitions
  # setting seed
  seed <- as.integer(1e7 * runif(1))
  set.seed(seed)
  
  print(paste0("seed: ", seed))
  
  # simulate capture histories
  sim_out <- make_history(N, M, k, omega, theta, pa, pb)
  
  print(paste("sum(z_init_1) =", sum(sim_out$z_init_1 == 1, na.rm = T)))
  print(paste("n =", sim_out$n))
  
  # save simulation data for debugging
  #write_sim_data(task_id_str, sim_out$h_aug, sim_out$h_2, sim_out$obs_a_aug, sim_out$obs_b_aug)
  
  #-------------------- Call nimble from R -----------------------#
  #start_time <- Sys.time()
  mcmc_out <- build_nimble_model(sim_out, chains = 1, iter = 50000, burnin = 10000)
  #end_time <- Sys.time()
  #end_time - start_time
  
  if (is.null(mcmc_out)) {
    print("Error: model failed to work")
    next
  }
  
  result_table[rep, 1] <- seed
  result_table[rep, 2:6] <- mcmc_out$summary[1,]
  result_table[rep, 7:11] <- mcmc_out$summary[2,]
  result_table[rep, 12:16] <- mcmc_out$summary[3,]
  result_table[rep, 17:21] <- mcmc_out$summary[4,]
  result_table[rep, 22:26] <- mcmc_out$summary[6,]
  
  print(paste0("Done rep #", rep))
 }

saveRDS(result_table, paste0("./coverage_res_full",task_id_str,".rds"))
print(paste0("Done task #", task_id_str))
toc()
