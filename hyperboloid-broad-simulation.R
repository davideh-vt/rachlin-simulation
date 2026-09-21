library(mvtnorm)
library(delaydiscount)
library(nlme)
library(multcomp) # TODO: Check if this was used
source("C:/Users/deh99/Documents/Classes/Work/hyperboloid-analysis-function.R")
source("C:/Users/deh99/Documents/Classes/Work/simulate-dataset-hyperboloid.R")
source("C:/Users/deh99/Documents/Classes/Work/hyperboloid-unconstrained-analysis.R")

 # Note: (X'X)^(-1) for the remedi dataset time points is:
# 2.1392701 -0.33509561
# -0.3350956  0.05624541

sim_results <- data.frame(matrix(nrow = 100, ncol = 58))
# need to add columns for true parameters
names(sim_results) <- c("EFT_ln_k", "EFT_s", "NCC_ln_k", "NCC_s", "sigma_sq", "Sigma_11", "Sigma_12", "Sigma_22", "n_EFT", "n_NCC",
                        "EFT_est_ln_k", "EFT_est_s", "NCC_est_ln_k", "NCC_est_s", "est_sigma_sq", "est_Sigma_11", "est_Sigma_12", "est_Sigma_22",
                        "est_con_sigma_sq", "est_con_g",
                        "red_est_ln_k", "red_est_s", "red_est_sigma_sq", "red_est_Sigma_11", "red_est_Sigma_12", "red_est_Sigma_22",
                        "red_est_con_sigma_sq", "red_est_con_g",
                        "mse_ln_k", "mse_s", "err_inprod",
                        "EFT_ln_k_mazur", "NCC_ln_k_mazur", "est_sigma_sq_mazur", "est_g_mazur",
                        "EFT_ln_k_nl", "EFT_s_nl", "NCC_ln_k_nl", "NCC_s_nl",
                        "var_ln_k_nl", "var_s_nl", "corr_ln_k_s_nl", "var_residual_nl",
                        "lik_ut_nl", "BIC_ut_nl", "deriv_offset", "p_value_nl", "p_value_nl_lrt",
                        "mse_ln_k_nl", "mse_s_nl", "err_inprod_nl",
                        "lik_uc", "lik_uc_red", "lik_con", "lik_con_red", "F_con", "lik_mazur", "F_mazur")

# Fill sim_results with parameters
sim_results$EFT_ln_k <- -5.86
sim_results$EFT_s <- 0.83
sim_results$NCC_ln_k <- -5.86
sim_results$NCC_s <- 0.83
sim_results$sigma_sq <- 1.359326
sim_results$Sigma_11 <- 9.604077
sim_results$Sigma_12 <- -1.1484446
sim_results$Sigma_22 <- 0.2002995
sim_results$n_EFT <- 35
sim_results$n_NCC <- 35

lmeControl(msMaxIter = 500)
time_points <- c(30, 90, 180, 365, 1095, 1825, 3650)
n_tp <- length(time_points)
x_mat <- matrix(nrow = n_tp, ncol = 2)
x_mat[,1] <- 1
x_mat[,2] <- log(time_points)

for(i in 1:100){
  # Simulate a dataset from the unconstrained hyperboloid model
  sim_data <- simulate_dataset_hyperboloid(groups = c("EFT", "NCC"),
                                           num_subj = c(sim_results$n_EFT[i], sim_results$n_NCC[i]),
                                           time_points = time_points,
                                           mean_param = matrix(c(sim_results$EFT_ln_k[i], sim_results$EFT_s[i], sim_results$NCC_ln_k[i], sim_results$NCC_s[i]), nrow = 2),
                                           sigma_sq = sim_results$sigma_sq[i],
                                           Sigma = matrix(c(sim_results$Sigma_11[i], sim_results$Sigma_12[i], sim_results$Sigma_12[i], sim_results$Sigma_22[i]), nrow = 2))
  
  # Make reduced data frame
  sim_data_red <- sim_data %>%
    mutate(group = "EFT")
  
  # Analyze under unconstrained model
  sim_huc <- hyperboloid_unconstrained_analysis(sim_data)
  sim_huc_red <- hyperboloid_unconstrained_analysis(sim_data_red)
  
  # Save unconstrained model parameters
  sim_results$EFT_est_ln_k[i] <- sim_huc$est_hyparms$ln_k[1]
  sim_results$NCC_est_ln_k[i] <- sim_huc$est_hyparms$ln_k[2]
  sim_results$EFT_est_s[i] <- sim_huc$est_hyparms$s[1]
  sim_results$NCC_est_s[i] <- sim_huc$est_hyparms$s[2]
  sim_results$est_sigma_sq[i] <- sim_huc$sigma_sq
  sim_results$est_Sigma_11[i] <- sim_huc$Sigma[1,1]
  sim_results$est_Sigma_12[i] <- sim_huc$Sigma[1,2]
  sim_results$est_Sigma_22[i] <- sim_huc$Sigma[2,2]
  sim_results$lik_uc[i] <- sim_huc$log_lik
  
  # Save unconstrained reduced model parameters
  sim_results$red_est_ln_k[i] <- sim_huc_red$est_hyparms$ln_k[1]
  sim_results$red_est_s[i] <- sim_huc_red$est_hyparms$s[1]
  sim_results$red_est_sigma_sq[i] <- sim_huc_red$sigma_sq
  sim_results$red_est_Sigma_11[i] <- sim_huc_red$Sigma[1,1]
  sim_results$red_est_Sigma_12[i] <- sim_huc_red$Sigma[1,2]
  sim_results$red_est_Sigma_22[i] <- sim_huc_red$Sigma[2,2]
  sim_results$lik_uc_red[i] <- sim_huc_red$log_lik

  # Fit constrained model
  sim_hc <- hyperboloid_analysis(sim_data)
  sim_results$est_con_sigma_sq[i] <- sim_hc$est_var["sigma_sq"]
  sim_results$est_con_g[i] <- sim_hc$est_var["g"]
  sim_results$lik_con[i] <- sim_hc$log_lik
  
  sim_hc_red <- hyperboloid_analysis(sim_data_red)
  sim_results$red_est_con_sigma_sq[i] <- sim_hc_red$est_var["sigma_sq"]
  sim_results$red_est_con_g[i] <- sim_hc_red$est_var["g"]
  sim_results$lik_con_red[i] <- sim_hc_red$log_lik
  
  # Next, need subject ln_k and s estimate differences
  true_subj_parms <- sim_data %>%
    select(subj, true_ln_k, true_s, group) %>%
    unique()
  # Estimates are in sim_huc$subj_parms
  ln_k_errs <- sim_huc$subj_parms$ln_k - true_subj_parms$true_ln_k
  s_errs <- sim_huc$subj_parms$s - true_subj_parms$true_s
  # Save ln(k) MSE, s MSE, and mean ln(k),s error inner product
  # (These components can be used to calculate Mahalanobis distance later)
  sim_results$mse_ln_k[i] <- mean(ln_k_errs^2)
  sim_results$mse_s[i] <- mean(s_errs^2)
  sim_results$err_inprod[i] <- mean(ln_k_errs*s_errs)
  
  # Linearized Mazur
  prep_sim_data <- delaydiscount::prepare_data_frame(sim_data)
  hyperbolic_fit <- delaydiscount::dd_hyperbolic_model(prep_sim_data)
  # Save parameters
  sim_results$EFT_ln_k_mazur[i] <- hyperbolic_fit$ln_k_mean$ln_k_mean[1]
  sim_results$NCC_ln_k_mazur[i] <- hyperbolic_fit$ln_k_mean$ln_k_mean[2]
  sim_results$est_sigma_sq_mazur[i] <- hyperbolic_fit$var[1]
  sim_results$est_g_mazur[i] <- hyperbolic_fit$var[2]
  
  # Save F-stat
  sim_results$F_mazur[i] <- hyperbolic_fit$model_test$F_stat
  
  # Nonlinear Estimated Rachlin (based on Young)
  # This "works" but is inelegant"
  best_lik <- -Inf
  for(test_s in c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
    try({
      nlRachFit_a<-nlme(indiff~1/(1+exp(logk)*delay^s), fixed=logk+s~as.factor(group), random=logk+s~1|as.factor(subj), data=sim_data,  start=c(-5.86,-5.86,test_s,test_s))
      if(nlRachFit_a$logLik > best_lik){
        nlRachFit <- nlRachFit_a
        best_lik <- nlRachFit$logLik
      }
    })
  }
  if(best_lik == -Inf){error("All iterations failed")}
  
  # Get hyperparameters
  sim_results$EFT_ln_k_nl[i] <- nlRachFit$coefficients$fixed[1]
  sim_results$EFT_s_nl[i] <- nlRachFit$coefficients$fixed[3]
  sim_results$NCC_ln_k_nl[i] <- nlRachFit$coefficients$fixed[1] + nlRachFit$coefficients$fixed[2]
  sim_results$NCC_s_nl[i] <- nlRachFit$coefficients$fixed[3] + nlRachFit$coefficients$fixed[4]
  # Get variance components
  vc <- VarCorr(nlRachFit)
  sim_results$var_ln_k_nl[i] <- as.numeric(vc[1,1])
  sim_results$var_s_nl[i] <- as.numeric(vc[2,1])
  sim_results$corr_ln_k_s_nl[i] <- as.numeric(vc[2,3])
  sim_results$var_residual_nl[i] <- as.numeric(vc[3,1])
  # Likelihood by NLME
  sim_results$lik_ut_nl[i] <- nlRachFit$logLik
  summNL <- summary(nlRachFit)
  sim_results$BIC_ut_nl[i] <- summNL$BIC  # NLME BIC treats N as total number of observations (num_subj * n_tp) and number of parameters as 8
  # 8 parameters likely to be:
  #  4 mean parameters, 2 for each group
  #  3 parameters for the random effect variance, which is a 2x2 symmetric matrix
  #  1 parameter for the independent random error
  # Save offset
  # You would add this offset to a likelihood based on transformed data to compare with untransformed
  # Or SUBTRACT from a likelihood based on untransformed data like the one returned by NLME
  sim_results$deriv_offset[i] <- sum(log(1/(sim_data$indiff - sim_data$indiff^2)))
  # Save p-value based on NLME
  # While NLME gives t-tests for effects, it does not give a result for the overall
  #  F-test that we care about, so we need to calculate an F-statistic ourselves
  #  based on what it does give us.
  # TODO: Check that this is reasonably accurate
  # First, get the estimated differences
  grp_diff <- summNL$coefficients$fixed[c(2,4)]
  # grp_diff <- summNL$tTable[c(2,4),1]
  # Then get the initial approximate var covar matrix
  prescale_covar <- summNL$varFix[c(2,4),c(2,4)]
  # Get scaling factor as the standard errors do not quite correspond to the var-covar matrix
  #  (but it is close, and the correlation does exactly correspond.)
  rescale_covar <- prescale_covar * (summNL$tTable[2,2]^2/prescale_covar[1,1])
  # Get denom degrees of freedom
  denom_df <- summNL$tTable[2,3]
  # Conduct F-test
  f_val <- t(grp_diff) %*% solve(rescale_covar) %*% grp_diff / 2
  sim_results$p_value_nl[i] <- pf(f_val, df1 = 2, df2 = denom_df, lower.tail = FALSE)
  
  # p-value based on reduced nonlinear model likelihood ratio test
  best_lik <- -Inf
  for(test_s in c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
    try({
      nlRachFit_red_a <- nlme(indiff~1/(1+exp(logk)*delay^s), fixed=logk+s~1, random=logk+s~1|as.factor(subj), data=sim_data,  start=c(-5.86,test_s))
      if(nlRachFit_red_a$logLik > best_lik){
        nlRachFit_red <- nlRachFit_red_a
        best_lik <- nlRachFit_red$logLik
      }
    })
  }
  if(best_lik == -Inf){error("All iterations failed")}
  sim_results$p_value_nl_lrt[i] <- pchisq(2*(nlRachFit$logLik - nlRachFit_red$logLik), df = 2, lower.tail = FALSE)
  

  # Get subject ln k estimates
  rand_effs <- nlRachFit$coefficients$random[[1]]
  nl_subj_ln_k_ests <- rand_effs[,1] + nlRachFit$coefficients$fixed[1] + ifelse(sim_data$group == "NCC", nlRachFit$coefficients$fixed[2], 0)
  nl_subj_s_ests <- rand_effs[,2] + nlRachFit$coefficients$fixed[3] + ifelse(sim_data$group == "NCC", nlRachFit$coefficients$fixed[4], 0)
  nl_subj_ln_k_errs <- sim_huc$subj_parms$ln_k-nl_subj_ln_k_ests
  nl_subj_s_errs <- sim_huc$subj_parms$s - nl_subj_s_ests
  sim_results$mse_ln_k_nl[i] <- mean(nl_subj_ln_k_errs^2)
  sim_results$mse_s_nl[i] <- mean(nl_subj_s_errs^2)
  sim_results$err_inprod_nl[i] <- mean(nl_subj_ln_k_errs*nl_subj_s_errs)
}

# Now implement formulas for important statistics

n <- sim_results$n_EFT + sim_results$n_NCC

# BIC
sim_results$BIC_lin_uc <- 8*log(n*n_tp) - 2*sim_results$lik_uc
# unconstrained linearized Rachlin model
  # 4 mean parameters
  # 3 random effect variance parameters
  # 1 independent variance parameter

sim_results$BIC_lin_con <- 6*log(n*n_tp) - 2*sim_results$lik_con
# constrained linearized Rachlin model
  # 4 mean parameters
  # 1 random effect variance parameter
  # 1 independent variance parameter

sim_results$BIC_mazur <- 4*log(n*n_tp) - 2*sim_results$lik_mazur
# linearized Mazur model
# 2 mean parameters
# 1 random effect variance parameter
# 1 independent variance parameter

sim_results$BIC_nl <- 8*log(n*n_tp) - 2*(sim_results$lik_ut_nl - sim_results$deriv_offset)

# TODO: Is it possible to formulate the constrained Rachlin F-test in terms of 
#  some quadratic form of the treatment difference similar to how the F-test in
#  Mazur could also be related to a t-test of the estimated subject ln(k)s?

# F-statistic for constrained model
sim_results$F_con <-
  (sim_results$red_est_con_g - sim_results$est_con_g)/(sim_results$est_con_g + 1)*(n-2)


# Likelihood from Linearized Mazur
sim_results$lik_mazur <- as.numeric(
  -n*n_tp/2*(log(2*pi)+1) - 
    n/2*(n_tp*log(sim_results$est_sigma_sq_mazur) + log(1+sim_results$est_g_mazur)))

# Mahalanobis distance

# Linearized estimates
sim_results$mh_dist_lin <- (
  sim_results$Sigma_22*sim_results$mse_ln_k -
    2*sim_results$Sigma_12*sim_results$err_inprod +
    sim_results$Sigma_11*sim_results$mse_s
)/(sim_results$Sigma_22*sim_results$Sigma_11 - sim_results$Sigma_12^2)  # divide by determinant

# Nonlinearized estimates
sim_results$mh_dist_nl <- (
  sim_results$Sigma_22*sim_results$mse_ln_k_nl -
    2*sim_results$Sigma_12*sim_results$err_inprod_nl +
    sim_results$Sigma_11*sim_results$mse_s_nl
)/(sim_results$Sigma_22*sim_results$Sigma_11 - sim_results$Sigma_12^2)

# Save the output
# TODO: Name the files by the simulation parameters
write.csv(sim_results, "out1.csv")


# p-values for nonlinearized based on LRT and on nlme estimate
#plot(sim_results$p_value_nl, sim_results$p_value_nl_lrt)
#hist(sim_results$p_value_nl)
#hist(sim_results$p_value_nl_lrt)

# Histograms of s hyperparameter estimates
#hist(c(sim_results$NCC_s_nl, sim_results$EFT_s_nl))  # nonlinear
#hist(c(sim_results$EFT_est_s, sim_results$NCC_est_s))  # linearized
#hist(c(sim_results$))

# Histogram of ln(k) hyperparameter estimates
#hist(c(sim_results$NCC_ln_k_nl, sim_results$EFT_ln_k_nl))
#hist(c(sim_results$EFT_est_ln_k, sim_results$NCC_est_ln_k))

# Scatterplot of ln(k), s hyperparameter estimate pairs (nonlinearized)
#plot(c(sim_results$NCC_ln_k_nl, sim_results$EFT_ln_k_nl),
#     c(sim_results$NCC_s_nl, sim_results$EFT_s_nl))


# Scatterplot of ln(k), s hyperparameter estimate pairs (linearized)
#plot(c(sim_results$NCC_est_ln_k, sim_results$EFT_est_ln_k),
#     c(sim_results$NCC_est_s, sim_results$EFT_est_s))

# Scatterplot of ln(k), s hyperparameter estimate difference pairs (nonlinearized)
#plot(c(sim_results$NCC_est_ln_k - sim_results$NCC_ln_k_nl, sim_results$EFT_est_ln_k - sim_results$EFT_ln_k_nl),
#     c(sim_results$NCC_est_s - sim_results$NCC_s_nl, sim_results$EFT_est_s - sim_results$EFT_s_nl))

# Histogram of p-values for likelihood ratio test using chi-square approximation for linearized unconstrained Rachlin hyperparameter difference
#hist(pchisq((sim_results$lik_uc - sim_results$lik_uc_red)*2, df = 2, lower.tail = FALSE))

# Histogram of p-values for F-test for constrained Rachlin
#hist(pf(sim_results$F_con, df1 = 2, df2 = 2*(n-2), lower.tail = FALSE))

# Histogram of p-values for F-test for Mazur
#hist(pf(sim_results$F_mazur, df1 = 1, df2 = n-2), lower.tail = FALSE)

# Plot of F-statistics for Mazur and constrained Rachlin
#plot(sim_results$F_con, sim_results$F_mazur)

# Plot of Rachlin F-stat and Rachlin likelihood difference
#plot(sim_results$F_con, sim_results$lik_con - sim_results$lik_con_red)
# sim_results$lik_con - sim_results$lik_con_red corresponds exactly to
# 70*(log(1+sim_results$red_est_con_g)-log(1+sim_results$est_con_g))
# (as it should, although the 70 will change for different sample sizes)

# Plot of Rachlin constrained and unconstrained likelihood difference
#plot(sim_results$lik_con - sim_results$lik_con_red, sim_results$lik_uc - sim_results$lik_uc_red)


# Histograms of BICs
#hist(sim_results$BIC_lin_uc)
#hist(sim_results$BIC_lin_con)
#hist(sim_results$BIC_mazur)
#hist(sim_results$BIC_nl)

# Comparisons
#hist(sim_results$BIC_lin_uc - sim_results$BIC_lin_con)
#hist(sim_results$BIC_lin_con - sim_results$BIC_mazur)
#hist(sim_results$BIC_lin_uc - sim_results$BIC_nl)     
#hist(sim_results$BIC_lin_con - sim_results$BIC_nl) 
