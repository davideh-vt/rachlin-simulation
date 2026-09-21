library(mvtnorm)

simulate_dataset_hyperboloid <- function(groups,
                                         num_subj,
                                         time_points,
                                         mean_param,
                                         sigma_sq, Sigma){
  
  # Check that the input is valid
  # REMOVED
  
  # Get number of time points
  n_tp <- length(time_points)
  
  # Set up output
  result <- data.frame(matrix(nrow = n_tp*sum(num_subj), ncol = 6))
  names(result) <- c("subj", "true_ln_k", "true_s", "group", "delay", "indiff")
  
  group_index <- c(0, cumsum(num_subj))
  
  for(i in 1:length(groups)){
    n_group <- num_subj[i]  # Number of subjects for this group
    
    # Simulate the true ln_k
    true_param <- rmvnorm(n_group, mean = mean_param[,i], sigma = Sigma)
    
    # Get all the terms that compose indiff_transform
    ln_k_term <- rep(true_param[,1], each = n_tp)
    s_term <- rep(true_param[,2], each = n_tp)
    random_term <- rnorm(n_group*n_tp, mean = 0, sd = sqrt(sigma_sq))
    
    indiff_transform <- ln_k_term + s_term*log(time_points) + random_term
    
    # Get original indiff
    indiff <- 1/(exp(indiff_transform) + 1)
    
    group_df <- data.frame(subj = rep((group_index[i]+1):group_index[i+1], each = n_tp),
                           true_ln_k = rep(true_param[,1], each = n_tp),
                           true_s = rep(true_param[,2], each = n_tp),
                           group = groups[i],
                           delay = rep(time_points, n_group),
                           indiff = indiff)
    
    result[(n_tp*group_index[i]+1):(n_tp*group_index[i+1]),] <- group_df
  }
  return(result)
}
