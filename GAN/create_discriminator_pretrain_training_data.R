create_discriminator_pretrain_training_data <- function(settings){
  
  source('baseline_raw_profile_simple.R')
  fake_profiles_dir <- settings$fake_profile_dirs
  profile_dirs <- settings$real_and_smooth_dirs
  number_of_directories <- length(profile_dirs)
  
  real_profile_paths <- NULL
  smooth_profiles_paths <- NULL
  fake_profiles_paths <- NULL
  
  for (iDir in 1:number_of_directories){
    # print(profile_dirs)
    sampleNames <- list.files(path = profile_dirs[iDir], pattern="\\.csv$")
    NoOfSamples <- length(sampleNames)
    
    #the samples with "smooth" in the title
    smooth_sample_indices <- grep(pattern = "smooth", sampleNames)
    #partition the samples
    real_samples <- sampleNames[-smooth_sample_indices]
    smoothed_samples <- sampleNames[smooth_sample_indices]
    fake_samples <- list.files(path = fake_profiles_dir, pattern="\\.csv$")
    #for each sample make sure a fake sample exists and then build up arrays of paths
    for (iSample in 1:length(real_samples)){
      sampleName <- real_samples[iSample]
      sample_name_without_suffix <- substr(x = sampleName, start = 1, stop = nchar(sampleName) - 4)
      corresponding_smooth_sample_index <- grep(pattern = sample_name_without_suffix, x = smoothed_samples)
      corresponding_fake_sample_index <- grep(pattern = sample_name_without_suffix, x = fake_samples)
      if (length(corresponding_fake_sample_index) == 1 & length(corresponding_smooth_sample_index) == 1){
        real_profile_paths <- c(real_profile_paths, paste(profile_dirs[iDir], sampleName, sep="/"))
        smooth_profiles_paths <- c(smooth_profiles_paths, paste(profile_dirs[iDir], smoothed_samples[corresponding_smooth_sample_index], sep="/"))
        fake_profiles_paths <- c(fake_profiles_paths, paste(fake_profiles_dir, fake_samples[corresponding_fake_sample_index], sep="/"))
      }
    }
  }
  
  number_sample_with_input_and_output <- length(real_profile_paths)
  targets <- rep(NA, 2*number_sample_with_input_and_output)
  #the number of dyes
  number_of_dyes <- settings$number_of_dyes
  #the starting scan point
  startScan <- settings$startScan
  #the number of scan points to be trained / generated
  number_of_scanpoints <- settings$number_of_scanpoints
  #the end scan point
  endScan <- startScan + number_of_scanpoints - 1 #need to minus 1 so that startScan:endScan gives number_of_scanpoints data points

  #now go through and read in files to create input and target arrays
  input_profiles <- array(0, dim = c(2*number_sample_with_input_and_output, number_of_dyes, number_of_scanpoints))
  source_profiles <- array(0, dim = c(2*number_sample_with_input_and_output, number_of_dyes, number_of_scanpoints))
  #saturation only applied to real profiles
  saturation <- settings$saturation
  
  for (iProfile in 1:number_sample_with_input_and_output){
    odd_spot <- 2*(iProfile - 1) + 1 #real profiles
    even_spot <- odd_spot + 1 #generated profiles
    #loads the real profiles
    real_profile_raw <- read.table(file = real_profile_paths[iProfile], header = TRUE, stringsAsFactors = FALSE, sep=";")
    
    #conducts simple profile baselining
    real_profile <- baseline_raw_profile_simple(real_profile_raw)
    
    #stores the scans in the input_array
    input_profiles[odd_spot,,] <- t(real_profile[startScan:endScan,2:(number_of_dyes + 1)])/saturation
    targets[odd_spot] <- 1
    #loads the fake profile
    fake_profile <- read.table(file = fake_profiles_paths[iProfile], header = TRUE, stringsAsFactors = FALSE, sep=";")
    #stores the scans in the target_array (note that saturation is not applied here asI save the fakes still normalised)
    input_profiles[even_spot,,] <- t(fake_profile[,2:(number_of_dyes + 1)])/saturation
    targets[even_spot] <- 0
    #stores the smoothed profile
    smooth_profile_raw <- read.table(file = smooth_profiles_paths[iProfile], header = TRUE, stringsAsFactors = FALSE, sep=";")
    #lazy, I store the same thing twice rather than being clever about it
    smooth_profile_formatted <- t(smooth_profile_raw[startScan:endScan,2:(number_of_dyes + 1)])/saturation
    source_profiles[odd_spot,,] <- smooth_profile_formatted
    source_profiles[even_spot,,] <- smooth_profile_formatted
  }
  
  return(list(input_profiles = input_profiles,
              input_paths = real_profile_paths,
              source_profile_smooth = source_profiles,
              targets = targets))
}
