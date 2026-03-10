create_generator_pretrain_training_data <- function(settings){
                  
  profile_dirs <- settings$real_and_smooth_dirs
  number_of_directories <- length(profile_dirs)
  
  input_profile_paths <- NULL
  target_profiles_paths <- NULL
  
  for (iDir in 1:number_of_directories){
    sampleNames <- list.files(path = profile_dirs[iDir], pattern="\\.csv$")
    #print(sampleNames)
    NoOfSamples <- length(sampleNames)
    
    #the samples with "smooth" in the title
    smoothed_sample_indices <- grep(pattern = "smooth", sampleNames)
    #print(smoothed_sample_indices)
    #partition the samples
    non_smoothed_samples <- sampleNames[-smoothed_sample_indices]
    smoothed_samples <- sampleNames[smoothed_sample_indices]
    #for each sample make sure a smoothed sample exists and then build up arrays of paths
    for (iSample in 1:length(non_smoothed_samples)){
      sampleName <- non_smoothed_samples[iSample]
      sample_name_without_suffix <- substr(x = sampleName, start = 1, stop = nchar(sampleName) - 4)
      corresponding_smoothed_sample_index <- grep(pattern = sample_name_without_suffix, x = smoothed_samples)
      if (length(corresponding_smoothed_sample_index) == 1){
        input_profile_paths <- c(input_profile_paths, paste(profile_dirs[iDir], smoothed_samples[corresponding_smoothed_sample_index], sep="/"))
        target_profiles_paths <- c(target_profiles_paths, paste(profile_dirs[iDir], sampleName, sep="/"))
      }
    }
  }
  
  
  number_sample_with_input_and_output <- length(input_profile_paths)
  #the number of dyes
  number_of_dyes <- settings$number_of_dyes
  #the starting scan point
  startScan <- settings$startScan
  #the number of scan points to be trained / generated
  number_of_scanpoints <- settings$number_of_scanpoints
  #the end scan point
  endScan <- startScan + number_of_scanpoints - 1 #need to minus 1 so that startScan:endScan gives number_of_scanpoints data points
  #now go through and read in files to create input and target arrays
  input_profiles <- array(0, dim = c(number_sample_with_input_and_output, number_of_dyes+1, number_of_scanpoints)) # +1 dim (molw)
  target_profiles <- array(0, dim = c(number_sample_with_input_and_output, number_of_dyes+1, number_of_scanpoints))
  
  for (iProfile in 1:number_sample_with_input_and_output){
    #loads the smoothed (input) profile
    smoothed_profile <- read.table(file = input_profile_paths[iProfile], header = TRUE, stringsAsFactors = FALSE, sep=";")
    #print(head(smoothed_profile))
    #print(head(smoothed_profile[startScan:endScan,2:(number_of_dyes+2)]))
    #stores the scans in the input_array
    input_profiles[iProfile,,] <- t(smoothed_profile[startScan:endScan,2:(number_of_dyes+2)]) #stores the molw as well
    #print(head(smoothed_profile[startScan:endScan,2:(number_of_dyes+2)]))
    
    #loads the raw profile
    raw_profile <- read.table(file = target_profiles_paths[iProfile], header = TRUE, stringsAsFactors = FALSE, sep=";")
    #print(dim(raw_profile))
    #print(head(raw_profile))
    
    #conducts simple profile baselining
    #source("baseline_raw_profile_simple.R")
    #baselined_raw_profile <- baseline_raw_profile_simple(raw_profile)
    
    #stores the scans in the target_array
    target_profiles[iProfile,,] <- t(raw_profile[startScan:endScan,2:(number_of_dyes+2)]) #stores molw as well
  }
  #the final step is to apply the saturation threshold
  saturation <- settings$saturation
  return(list(input_profiles = input_profiles/saturation,
              target_profiles = target_profiles/saturation,
              input_paths = input_profile_paths,
              target_paths = target_profiles_paths))
  
}

