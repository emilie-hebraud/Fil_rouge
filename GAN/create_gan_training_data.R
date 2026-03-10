create_gan_training_data <- function(settings){
  
  profile_dirs <- settings$real_and_smooth_dirs
  
  number_of_directories <- length(profile_dirs)
  
  real_profile_paths <- NULL
  smooth_profiles_paths <- NULL
  
  for (iDir in 1:number_of_directories){
    # Directory containing smooth profiles
    smooth_dir <- paste0(here(), "/pre_process_data/real")
    
    # List input profiles (only real .csv files, no smooth)
    input_files <- list.files(
      path = profile_dirs[iDir],
      pattern = "\\.csv$",
      full.names = FALSE
    )
    # Remove smooth files from input list
    input_files <- input_files[!grepl("_smooth\\.csv$", input_files)]
    
    # List only smooth csv files
    smooth_files <- list.files(
      path = smooth_dir,
      pattern = "_smooth\\.csv$",
      full.names = FALSE
    )
    
    # Extract comparable IDs
    input_ids  <- sub("\\.csv$", "", input_files)
    smooth_ids <- sub("_smooth\\.csv$", "", smooth_files)
    
    # Match input profiles to their smooth counterparts
    match_idx <- match(input_ids, smooth_ids)
    #print(length(match_idx))
    
    # Build paths only when a match exists
    for (iSample in seq_along(input_files)) {
      
      if (!is.na(match_idx[iSample])) {
        
        real_profile_paths <- c(
          real_profile_paths,
          file.path(profile_dirs[iDir], input_files[iSample])
        )
        
        smooth_profiles_paths <- c(
          smooth_profiles_paths,
          file.path(smooth_dir, smooth_files[match_idx[iSample]])
        )
      }
    }
  }
    
  number_sample_with_input_and_output <- length(real_profile_paths)
  #the number of dyes
  number_of_dyes <- settings$number_of_dyes
  #the starting scan point
  startScan <- settings$startScan
  #the number of scan points to be trained / generated
  number_of_scanpoints <- settings$number_of_scanpoints
  #the end scan point
  endScan <- startScan + number_of_scanpoints - 1 #need to minus 1 so that startScan:endScan gives number_of_scanpoints data points
  #now go through and read in files to create input and target arrays
  real_profiles <- array(0, dim = c(number_sample_with_input_and_output, number_of_dyes+1, number_of_scanpoints)) #+1 to include molw
  smooth_profiles <- array(0, dim = c(number_sample_with_input_and_output, number_of_dyes+1, number_of_scanpoints))
  generated_profiles <- array(0, dim = c(number_sample_with_input_and_output, number_of_dyes+1, number_of_scanpoints))
  #saturation only applied to real profiles
  saturation <- settings$saturation
  
  for (iProfile in 1:number_sample_with_input_and_output){
    #loads the real profiles
    real_profile_raw <- read.table(file = real_profile_paths[iProfile], header = TRUE, stringsAsFactors = FALSE, sep=";")
  
    #real_profile_raw <- baseline_raw_profile_simple(real_profile_raw_unbaselined)
    #stores the scans in the input_array
    real_profiles[iProfile,,] <- t(real_profile_raw[startScan:endScan,2:(number_of_dyes + 2)])/saturation #molw added
    #loads and stores the smoothed profile
    smooth_profile_raw <- read.table(file = smooth_profiles_paths[iProfile], header = TRUE, stringsAsFactors = FALSE, sep=";")
    #print(head(smooth_profile_raw))
    smooth_profiles[iProfile,,] <- t(smooth_profile_raw[startScan:endScan,2:(number_of_dyes + 2)])/saturation #molw added
  }
  
  return(list(real_profile_target = real_profiles,
              input_paths = real_profile_paths,
              source_profile_smooth = smooth_profiles))
  
}
