create_noisy_idealized_profile <- function(raw_profile, peak_labels, noise_level = 0.3) {
  
  idealized_profile <- data.frame(scan = raw_profile$scan)
  
  # For each channel
  for (channel in 2:ncol(raw_profile)) {
    
    # Initialize with baseline noise
    baseline_noise <- rnorm(nrow(raw_profile), mean = 50, sd = 20)
    baseline_noise[baseline_noise < 0] <- 0
    idealized_signal <- baseline_noise
    
    # Get peaks for this channel
    channel_peaks <- peak_labels[peak_labels$channel == (channel - 1), ]
    
    # For each labeled peak
    for (i in 1:nrow(channel_peaks)) {
      peak_scan <- channel_peaks$scan[i]
      peak_height <- channel_peaks$height[i]
      
      # Add variation to peak characteristics
      actual_height <- peak_height * rnorm(1, mean = 1, sd = 0.1)  # ±10% height variation
      actual_position <- peak_scan + rnorm(1, mean = 0, sd = 0.5)  # slight position jitter
      peak_width <- rnorm(1, mean = 5, sd = 0.5)  # variable peak width
      
      # Create Gaussian peak with noise
      scan_range <- (peak_scan - 15):(peak_scan + 15)
      scan_range <- scan_range[scan_range > 0 & scan_range <= nrow(raw_profile)]
      
      # Gaussian formula
      gaussian_peak <- actual_height * exp(-((scan_range - actual_position)^2) / (2 * peak_width^2))
      
      # Add noise to the peak shape itself
      peak_noise <- gaussian_peak * rnorm(length(gaussian_peak), mean = 0, sd = noise_level)
      gaussian_peak <- gaussian_peak + peak_noise
      gaussian_peak[gaussian_peak < 0] <- 0
      
      # Add to idealized signal
      idealized_signal[scan_range] <- idealized_signal[scan_range] + gaussian_peak
    }
    
    # Add some random small peaks (false positives the model should learn to ignore)
    n_random_peaks <- sample(0:5, 1)
    for(j in 1:n_random_peaks){
      random_scan <- sample(1:nrow(raw_profile), 1)
      random_height <- runif(1, 50, 200)  # small random peaks
      random_width <- rnorm(1, mean = 3, sd = 0.5)
      
      scan_range <- (random_scan - 10):(random_scan + 10)
      scan_range <- scan_range[scan_range > 0 & scan_range <= nrow(raw_profile)]
      
      random_peak <- random_height * exp(-((scan_range - random_scan)^2) / (2 * random_width^2))
      idealized_signal[scan_range] <- idealized_signal[scan_range] + random_peak
    }
    
    idealized_profile[, paste0("channel_", channel-1)] <- idealized_signal
  }
  
  return(idealized_profile)
}