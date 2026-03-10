read_raw_profile <- function(raw_profile,
                             create_predictions_vector=FALSE){
  
  
  if (create_predictions_vector) {
    raw_profile <- raw_profile[, c("index", "channel_1", "channel_2",
                                   "channel_3", "channel_4", "molw")]
  }
  # Keep relevant columns
  else{
  raw_profile <- raw_profile[, 1:5]}
  
  current_n <- nrow(raw_profile)
  
  # If already >= 5000, truncate to 5000
  if(current_n >= 5000){
    return(raw_profile[1:5000, ])
  }
  if(current_n == 5000){
    return(raw_profile)
  }
  
  # Baseline window (use the last available 500 points if 4000:4500 not valid)
  baseline_start <- max(1, current_n - 500)
  baseline_end   <- current_n
  
  baseline_block <- raw_profile[baseline_start:baseline_end,
                                2:ncol(raw_profile)]
  
  # Estimate Gaussian parameters per channel
  baseline_mean <- colMeans(baseline_block, na.rm = TRUE)
  
  baseline_sd <- apply(baseline_block, 2, sd, na.rm = TRUE)
  
  # Protect against degenerate variance
  baseline_sd[!is.finite(baseline_sd) | baseline_sd <= 0] <- 1e-6
  
  # Number of rows to generate
  n_pad <- 5000 - current_n
  
  # Generate noise per channel
  padding_channels <- sapply(
    seq_along(baseline_mean),
    function(i){
      rnorm(n_pad,
            mean = baseline_mean[i],
            sd   = baseline_sd[i])
    }
  )
  
  padding_channels <- as.matrix(padding_channels)
  
  # Build index column
  padding_index <- seq(current_n + 1, 5000)
  
  padding <- cbind(padding_index, padding_channels)
  colnames(padding) <- colnames(raw_profile)
  # Append
  raw_profile <- rbind(raw_profile, padding)

  # Reset index in first column
  raw_profile[, 1] <- 1:nrow(raw_profile)
  
  return(raw_profile)
}
