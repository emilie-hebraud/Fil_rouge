extrapolate_molw <- function(molw_vec) {
  n <- length(molw_vec)
  idx <- seq_along(molw_vec)
  
  first_valid <- min(which(!is.na(molw_vec)))
  last_valid  <- max(which(!is.na(molw_vec)))
  
  valid_idx  <- idx[first_valid:last_valid]
  valid_molw <- molw_vec[first_valid:last_valid]
  fit <- lm(valid_molw ~ valid_idx)
  
  na_idx <- which(is.na(molw_vec))
  molw_vec[na_idx] <- predict(fit, newdata = data.frame(valid_idx = na_idx))
  
  return(molw_vec)
}

pre_process_data <- function(profile_dirs, real=FALSE){
  number_of_directories <- length(profile_dirs)
  print(number_of_directories)
  for (iDir in 1:number_of_directories){
    SampleNames <- list.files(path = profile_dirs[iDir], pattern="\\.csv$")
    NoOfSamples <- length(SampleNames)
    print(NoOfSamples)
    for (iSample in 1:NoOfSamples){
      sampleName <- SampleNames[iSample]
      print(paste("pre processing sample ", sampleName, sep=""))
      
      config <- list(
        features_path = paste(profile_dirs[iDir], sampleName, sep="/"),
        labels_path   = paste(here(), "input_data/peak_positions_detailed.csv", sep="/"),
        channels      = c("channel_1", "channel_3", "channel_4"),
        mw_window     = 10,
        #output_path   = paste(here(), "pre_process_data", sampleName, sep="/")
        output_path   = paste("/home/emilie/Documents/Fil_rouge/Data/Test/Test/pre_processed", sampleName, sep="/")
      )
      
      raw_profile <- read.table(config$features_path, sep=";", header=TRUE)
      raw_profile <- raw_profile[, c("index","channel_1","channel_3","channel_4","molw")]
      
      if(real){
        raw_profile <- raw_profile[raw_profile$molw > 100, ]
      }
      
      current_n <- nrow(raw_profile)
      
      if(current_n >= 5000){
        raw_profile <- raw_profile[1:5000, ]
      } else {
        baseline_start <- max(1, current_n - 500)
        baseline_block <- raw_profile[baseline_start:current_n,
                                      c("channel_1","channel_3","channel_4")]
        baseline_mean <- colMeans(baseline_block, na.rm=TRUE)
        baseline_sd   <- apply(baseline_block, 2, sd, na.rm=TRUE)
        baseline_sd[!is.finite(baseline_sd) | baseline_sd <= 0] <- 1e-6
        
        n_pad <- 5000 - current_n
        padding_channels <- sapply(seq_along(baseline_mean), function(i){
          rnorm(n_pad, mean=baseline_mean[i], sd=baseline_sd[i])
        })
        padding <- data.frame(
          index = seq(current_n + 1, 5000),
          padding_channels,
          molw  = NA
        )
        colnames(padding) <- colnames(raw_profile)
        raw_profile <- rbind(raw_profile, padding)
        raw_profile[, 1] <- seq_len(nrow(raw_profile))
      }
      
      # Extrapolate molw over padded NAs (and any leading NAs)
      raw_profile$molw <- extrapolate_molw(raw_profile$molw)
      
      write.table(raw_profile,
                  file      = config$output_path,
                  sep       = ";",
                  row.names = FALSE,
                  col.names = TRUE)
    }
  }
}