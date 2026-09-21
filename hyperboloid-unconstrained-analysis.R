library(dplyr)

hyperboloid_unconstrained_analysis <- function(dd_data){
  # Check that input passes preconditions
  check_input_preconditions(dd_data)
  
  # first, make sure that it does not have extra variables
  dd_data <- dd_data %>%
    dplyr::select(.data$subj, .data$group, .data$delay, .data$indiff) %>%
    dplyr::arrange(.data$group, .data$subj, .data$delay)
  
  # add transformations of variables
  dd_data$log_delay <- log(dd_data$delay)
  dd_data$indiff_transform <- log(1/dd_data$indiff - 1)
  
  # Fit hyperboloid model 
  #hyperboloid_group_fit <- lm(indiff_transform ~ group*log_delay, data = dd_data)
  # Fit subject model
  # TODO: note that this assumes all subjects are distinct
  #hyperboloid_subj_fit <- lm(indiff_transform ~ subj*log_delay, data = dd_data)
  
  # get xmat
  n_tp <- length(unique(dd_data$log_delay))
  x_mat <- matrix(nrow = n_tp, ncol = 2)
  x_mat[,1] <- 1
  x_mat[,2] <- unique(dd_data$log_delay)
  
  subj_effect_mat <- solve(t(x_mat) %*% x_mat) %*% t(x_mat)
  y_mat <- matrix(dd_data$indiff_transform, nrow = n_tp)
  subj_param <- subj_effect_mat %*% y_mat
  
  subj_param_df <- dd_data %>%
    dplyr::select(subj, group) %>%
    unique()
  subj_param_df$ln_k <- subj_param[1,]
  subj_param_df$s <- subj_param[2,]
  
  dd_data = merge(dd_data, subj_param_df) %>%
    dplyr::arrange(.data$group, .data$subj, .data$delay)
  
  # Calculate residual hyperboloid
  dd_data$fitted_indiff_transform <- dd_data$ln_k + dd_data$s*dd_data$log_delay
  dd_data$residual_hyperboloid <- dd_data$indiff_transform - dd_data$fitted_indiff_transform
  
  # Estimate hyperparameters
  est_hyparms <- subj_param_df %>%
    group_by(.data$group) %>%
    summarise(ln_k = mean(ln_k), s = mean(s))
  
  # Estimate variance parameters
  n <- length(unique(dd_data$subj))
  sse_z <- sum(dd_data$residual_hyperboloid^2)
  est_sigma_sq <- sse_z/(n*(n_tp-2))
  
  subjwise_hyparms <- dd_data %>%
    dplyr::select(group) %>%
    merge(est_hyparms)
  fitted_group_vals <- subjwise_hyparms$ln_k + subjwise_hyparms$s*dd_data$log_delay
  
  proj_mean_diff <- matrix(dd_data$fitted_indiff_transform - fitted_group_vals, 
                           ncol = n_tp, byrow = TRUE)
  
  # Do PCA
  pc_pmd <- princomp(proj_mean_diff)
  
  v1 <- pc_pmd$loadings[1:n_tp]
  v2 <- pc_pmd$loadings[(n_tp+1):(2*n_tp)]
  
  # These are also the first and second eigenvalues from pc_pmd
  # TODO: the t(v1) %*% v1 may just be 1
  u1 <- sum(diag(proj_mean_diff %*% v1 %*% solve(t(v1) %*% v1) %*% t(v1) %*% t(proj_mean_diff)))/n
  u2 <- sum(diag(proj_mean_diff %*% v2 %*% solve(t(v2) %*% v2) %*% t(v2) %*% t(proj_mean_diff)))/n
  
  # Get lambda values from u1 and u2
  lambda1 <- u1 - est_sigma_sq
  lambda2 <- u2 - est_sigma_sq
  
  XSXt <- lambda1*v1%*%t(v1) + lambda2*v2%*%t(v2)
  XtX_inv <- solve(t(x_mat) %*% x_mat)
  Sigma_hat <- XtX_inv %*% t(x_mat) %*% XSXt %*% x_mat %*% XtX_inv
  
  # TODO: Consider the scenario when a lambda is not positive
  
  return(list(
    dd_data = dd_data,
    subj_parms = subj_param_df,
    est_hyparms = est_hyparms,
    sigma_sq = est_sigma_sq,
    Sigma = Sigma_hat,
    log_lik = as.numeric(
      -(n*n_tp)*(0.5)*log(2*pi) - 
        n/2*((n_tp-2)*log(est_sigma_sq) + log(u1) + log(u2)) -
        n_tp/2*n)
  ))
}
