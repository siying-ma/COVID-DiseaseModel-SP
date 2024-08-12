setwd("~/Desktop/COVID/ComputeCan/full-model-with-lag-rep/results")
coverage_res_full1 <- readRDS("./coverage_res_full1.rds")
coverage_res_full2 <- readRDS("./coverage_res_full2.rds")
coverage_res_full3 <- readRDS("./coverage_res_full3.rds")
coverage_res_full4 <- readRDS("./coverage_res_full4.rds")
coverage_res_full5 <- readRDS("./coverage_res_full5.rds")
coverage_res_full19 <- readRDS("./coverage_res_full19.rds")
coverage_res_full6 <- readRDS("./coverage_res_full6.rds")
coverage_res_full7 <- readRDS("./coverage_res_full7.rds")
coverage_res_full8 <- readRDS("./coverage_res_full8.rds")
coverage_res_full9 <- readRDS("./coverage_res_full9.rds")
coverage_res_full10 <- readRDS("./coverage_res_full10.rds")
coverage_res_full11 <- readRDS("./coverage_res_full11.rds")
coverage_res_full12 <- readRDS("./coverage_res_full12.rds")
coverage_res_full13 <- readRDS("./coverage_res_full13.rds")
coverage_res_full14 <- readRDS("./coverage_res_full14.rds")
coverage_res_full15 <- readRDS("./coverage_res_full15.rds")
coverage_res_full16 <- readRDS("./coverage_res_full16.rds")
coverage_res_full17 <- readRDS("./coverage_res_full17.rds")
coverage_res_full18 <- readRDS("./coverage_res_full18.rds")
coverage_res_full19 <- readRDS("./coverage_res_full19.rds")
coverage_res_full20 <- readRDS("./coverage_res_full20.rds")

coverage_res <- rbind(coverage_res_full1, coverage_res_full2, coverage_res_full3, coverage_res_full4,
                      coverage_res_full5, coverage_res_full6, coverage_res_full7, coverage_res_full8,
                      coverage_res_full9, coverage_res_full10, coverage_res_full11, coverage_res_full12,
                      coverage_res_full13, coverage_res_full14, coverage_res_full15, coverage_res_full16,
                      coverage_res_full17, coverage_res_full18, coverage_res_full19, coverage_res_full20)
coverage_res = data.frame(coverage_res)
#estimates
est_N <- mean(coverage_res$N_mean)
est_omega <- mean(coverage_res$omega_mean)
est_pa <- mean(coverage_res$pa_mean)
est_pb <- mean(coverage_res$pb_mean)
est_theta <- mean(coverage_res$theta_mean)

#coverage rate
cover_N <- sum(coverage_res$N_lower<=2967 & coverage_res$N_upper>=2967)
cover_N/100

cover_omega <- sum(coverage_res$omega_lower<=0.1057 & coverage_res$omega_upper>=0.1057)
cover_omega/100

cover_pa <- sum(coverage_res$pa_lower<=0.1128 & coverage_res$pa_upper>=0.1128)
cover_pa/100

cover_pb <- sum(coverage_res$pb_lower<=0.4 & coverage_res$pb_upper>=0.4)
cover_pb/100

cover_theta <- sum(coverage_res$theta_lower <= 0.15 & coverage_res$theta_upper >= 0.15)
cover_theta/100

#bias
bias_N <- sum(coverage_res$N_mean)/100 - 2967
bias_N

bias_omega <- sum(coverage_res$omega_mean)/100 - 0.1057
bias_omega

bias_pa  <- sum(coverage_res$pa_mean)/100 - 0.1128
bias_pa

bias_pb <- sum(coverage_res$pb_mean)/100 - 0.4
bias_pb

bias_theta <- sum(coverage_res$theta_mean)/100 - 0.15
bias_theta

#mean square error
mse_N <- sqrt(sum((coverage_res$N_mean-2967)^2)/100)
mse_N

mse_omega <- sqrt(sum((coverage_res$omega_mean-0.1057)^2)/100)
mse_omega

mse_pa <- sqrt(sum((coverage_res$pa_mean-0.1128)^2)/100)
mse_pa

mse_pb <- sqrt(sum((coverage_res$pb_mean-0.4)^2)/100)
mse_pb

mse_theta <- sum((coverage_res$theta_mean-0.15)^2)/100
mse_theta

library(xtable)
coverage_table <- data.frame(estimates = c(est_N, est_omega, est_pa, est_pb, est_theta),
                             coverage = c(cover_N, cover_omega, cover_pa, cover_pb, cover_theta),
                             bias = c(bias_N, bias_omega, bias_pa, bias_pb, bias_theta),
                             mse = c(mse_N, mse_omega, mse_pa, mse_pb, mse_theta))
round(coverage_table,5)
