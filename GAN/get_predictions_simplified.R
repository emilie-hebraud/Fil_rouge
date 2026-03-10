get_predictions_simplified <- function(raw_profile, model) {
  
  n_channels     <- 4
  context_window <- 100
  saturation     <- 33000
  
  n_scans    <- nrow(raw_profile)
  startScan  <- 1 + context_window
  endScan    <- n_scans - context_window
  categories <- c("B", "P", "S")
  
  # Initialize predictions
  predictions <- matrix("B", nrow = n_scans, ncol = n_channels)
  
  # Process scan by scan
  for (iScan in startScan:endScan) {
    
    # Build input array: (1, 4, 201, 1)
    input_data <- array(0, dim = c(1, n_channels, 1 + 2*context_window, 1))
    
    scan_range <- (iScan - context_window):(iScan + context_window)
    
    for (idye in 1:n_channels) {
      # Extract window
      window <- raw_profile[scan_range, idye + 1]
      # Normalize: subtract mode, divide by saturation
      mode_val <- getmode(window)
      window   <- (window - mode_val) / saturation
      input_data[1, idye, , 1] <- window
    }
    
    # Convert to tensor and predict
    input_tensor <- tensorflow::tf$constant(input_data, dtype = tensorflow::tf$float32)
    scan_predictions <- model$predict(input_tensor, verbose = 0L)
    
    # Get predicted category for each channel
    for (idye in 1:n_channels) {
      probs     <- scan_predictions[[idye]][1, 1, 1, ]
      predicted <- categories[which.max(probs)]
      predictions[iScan, idye] <- predicted
    }
  }
  
  return(predictions)
}