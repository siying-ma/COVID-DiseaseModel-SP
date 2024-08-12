setwd("~/Desktop/COVID/ComputeCan/Simplified-model-without-lag-rep/Result")
coverage_res1 <- readRDS("./coverage_res1.rds")
coverage_res2 <- readRDS("./coverage_res2.rds")
coverage_res3 <- readRDS("./coverage_res3.rds")
coverage_res4 <- readRDS("./coverage_res4.rds")
coverage_res5 <- readRDS("./coverage_res5.rds")
coverage_res6 <- readRDS("./coverage_res6.rds")
coverage_res7 <- readRDS("./coverage_res7.rds")
coverage_res8 <- readRDS("./coverage_res8.rds")
coverage_res9 <- readRDS("./coverage_res9.rds")
coverage_res10 <- readRDS("./coverage_res10.rds")
coverage_res11 <- readRDS("./coverage_res11.rds")
coverage_res12 <- readRDS("./coverage_res12.rds")
coverage_res13 <- readRDS("./coverage_res13.rds")
coverage_res14 <- readRDS("./coverage_res14.rds")
coverage_res15 <- readRDS("./coverage_res15.rds")
coverage_res16 <- readRDS("./coverage_res16.rds")
coverage_res17 <- readRDS("./coverage_res17.rds")
coverage_res18 <- readRDS("./coverage_res18.rds")
coverage_res19 <- readRDS("./coverage_res19.rds")
coverage_res20 <- readRDS("./coverage_res20.rds")

coverage_res <- rbind(coverage_res1, coverage_res2, coverage_res3, coverage_res4,
                      coverage_res5, coverage_res6, coverage_res7, coverage_res8,
                      coverage_res9, coverage_res10, coverage_res11, coverage_res12,
                      coverage_res13, coverage_res14, coverage_res15, coverage_res16,
                      coverage_res17, coverage_res18, coverage_res19, coverage_res20)
coverage_res = data.frame(coverage_res)
#estimates
est_N <- mean(coverage_res$N_mean)
est_omega <- mean(coverage_res$omega_mean)
est_pa <- mean(coverage_res$pa_mean)
est_pb <- mean(coverage_res$pb_mean)

#coverage rate
cover_N <- sum(coverage_res$N_lower<=2967 & coverage_res$N_upper>=2967)
cover_N/100
#0.97

cover_omega <- sum(coverage_res$omega_lower<=0.1057 & coverage_res$omega_upper>=0.1057)
cover_omega/100
#0.95

cover_pa <- sum(coverage_res$pa_lower<=0.1128 & coverage_res$pa_upper>=0.1128)
cover_pa/100
#0.85

cover_pb <- sum(coverage_res$pb_lower<=0.0599 & coverage_res$pb_upper>=0.0599)
cover_pb/100
#0.85

#bias
bias_N <- sum(coverage_res$N_mean)/100 - 2967
bias_N
#2.449723

bias_omega <- sum(coverage_res$omega_mean)/100 - 0.1057
bias_omega
#0.004546708

bias_pa  <- sum(coverage_res$pa_mean)/100 - 0.1128
bias_pa
#0.000367445

bias_pb <- sum(coverage_res$pb_mean)/100 - 0.0599
bias_pb
#-0.004326595

#mean square error
mse_N <- sqrt(sum((coverage_res$N_mean-2967)^2)/100)
mse_N

mse_omega <- sqrt(sum((coverage_res$omega_mean-0.1057)^2)/100)
mse_omega

mse_pa <- sqrt(sum((coverage_res$pa_mean-0.1128)^2)/100)
mse_pa

mse_pb <- sqrt(sum((coverage_res$pb_mean-0.0599)^2)/100)
mse_pb


library(xtable)
coverage_table <- data.frame(estimates = c(est_N, est_omega, est_pa, est_pb),
                             coverage = c(cover_N, cover_omega, cover_pa, cover_pb),
                             bias = c(bias_N, bias_omega, bias_pa, bias_pb),
                             mse = c(mse_N, mse_omega, mse_pa, mse_pb))
round(coverage_table,5)
xtable(round(coverage_table,5))
