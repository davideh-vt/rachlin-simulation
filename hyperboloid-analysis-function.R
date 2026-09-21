library(dplyr)

hyperboloid_analysis <- function(dd_data){
  # Check that input passes preconditions
  check_input_preconditions(dd_data)
  
  # first, make sure that it does not have extra variables
  dd_data <- dd_data %>%
    dplyr::select(.data$subj, .data$group, .data$delay, .data$indiff) %>%
    dplyr::arrange(.data$group, .data$subj, .data$delay)
  
  # add transformations of variables
  dd_data$log_delay <- log(dd_data$delay)
  dd_data$indiff_transform <- log(1/dd_data$indiff - 1)
  
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
  sse_z <- sum(dd_data$residual_hyperboloid^2)
  subjwise_hyparms <- dd_data %>%
    dplyr::select(group) %>%
    merge(est_hyparms)
  fitted_group_vals <- subjwise_hyparms$ln_k + subjwise_hyparms$s*dd_data$log_delay
  sse_x <- sum((dd_data$indiff_transform-fitted_group_vals)^2)
  ssr_z_x <- sse_x - sse_z
  g_hat <- max(0, ssr_z_x/sse_z * (n_tp-2)/2 - 1)
  sigma_sq_hat <- ifelse(g_hat == 0,
                         sse_x/length(dd_data$subj),
                         sse_z/(length(subj_param_df$subj)*(n_tp-2)))
  
  
  est_var <- c(sigma_sq_hat, g_hat)
  names(est_var) <- c("sigma_sq", "g")
  
  # Likelihood
  n <- length(subj_param_df$subj)
  loglik = as.numeric(
    -(n*n_tp)*(0.5)*log(2*pi) - 
      n/2*(n_tp*log(sigma_sq_hat) + 2*log(1+g_hat)) - 
      (n_tp/2*n))
  
  result <- list(dd_data, subj_param_df, est_hyparms, est_var, loglik)
  names(result) = c("dd_data", "subj_parms", "est_hyparms", "est_var", "log_lik")
  
  return(result)
}
