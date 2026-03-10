
simulate_profile <- function(profile_simulation_settings){

  #loads libraries
  library('here')
  library('simDNAmixtures')
  #sources
  source('setup_profile_simulation.R')
  source('convert_bp_to_scan.R')
  
  #allele frequencies
  freqs <- read_allele_freqs(paste(here(), "\\input_data\\AlleleFreq\\Australian Caucasian.csv", sep=""))
  
  #create the model settings
  model_settings <- setup_profile_simulation()
  
  #define the sampling parameters
  sampling_parameters <- list(min_template = profile_simulation_settings$min_template, #note that these can be arrays if you want to create specific mixture proportions
                              max_template = profile_simulation_settings$max_template, #note that these can be arrays if you want to create specific mixture proportions
                              degradation_shape = profile_simulation_settings$degradation_shape, 
                              degradation_scale = profile_simulation_settings$degradation_scale) 
  
  #the bounds on the profiles we wish to simulate
  number_of_simulated_profiles <- profile_simulation_settings$number_of_simulated_profiles
  upper_NoC <- profile_simulation_settings$upper_NoC
  lower_NoC <- profile_simulation_settings$lower_NoC
  simulated_profiles <- list()
  NoC_array <- rep(NA, number_of_simulated_profiles)
  
  for (iProfile in 1:number_of_simulated_profiles){
    NoC <- sample(lower_NoC:upper_NoC, 1)
    NoC_array[iProfile] <- NoC
    #sample n mixtures using the log_normal model
    sim_name <- paste("profile_", iProfile, "_NoC_", NoC, sep="")
    simulated_profiles[[sim_name]] <- sample_mixtures(n = 1,
                                                      contributors = sprintf("U%s", 1:NoC),
                                                      freqs = freqs,
                                                      sampling_parameters = sampling_parameters,
                                                      model_settings = model_settings,
                                                      sample_model = sample_log_normal_model
                                                      )
  }
  
  #get one of the mixtures from the simulated_profiles object
  profile_focus <- 1
  
  one_profile <- simulated_profiles[[profile_focus]]$samples[[1]]$mixture
  #need to add in Amelogenin (as the simulator tool doesn't include this locus and it will be obvious when comparing to fake data)
  templates <- rep(NA, NoC_array[profile_focus])
  for (iNoC in 1:NoC_array[profile_focus]){
    templates[iNoC] <- simulated_profiles[[profile_focus]]$parameter_summary[[paste("template", iNoC, sep="")]]
  }
  #choose random male or female for each contrib (1 = male, 0 = female)
  sex_array <- sample(0:1, NoC_array[profile_focus], replace = TRUE)
  x_height <- 0
  y_height <- 0
  for (iNoC in 1:NoC_array[profile_focus]){
    if (sex_array[iNoC] == 0){
      x_height <- x_height + templates[iNoC]*2
      #TODO  perturb using variance param in future
    } else {
      x_height <- x_height + templates[iNoC]
      y_height <- y_height + templates[iNoC]
      #TODO  perturb using variance param in future
    }
  }
  #Add to list
  if (x_height > 0){
    one_profile_addition <- c("AMEL", "X", round(x_height, 0), 106)
    one_profile <- rbind(one_profile, one_profile_addition)
  }
  if (y_height > 0){
    one_profile_addition <- c("AMEL", "Y", round(y_height, 0), 113)
    one_profile <- rbind(one_profile, one_profile_addition)
  }
  
  #adds a scan column to the profile
  one_profile$Scans <- convert_bp_to_scan(one_profile$Size) 
  return(one_profile)
}

extended_profile_input <- function(one_profile, one_simulated_profile){
  #get the annotation
  annotate_mix <- one_simulated_profile$samples[[1]]$annotated_mixture
  #cut down to just what was detected
  annotate_mix <- annotate_mix[annotate_mix$HeightAtOrAboveDetectionThreshold,]
  #compile the parts that can be used as inputs used for ANN reader
  column_names_to_target <- c("HeightAllele", "HeightUncappedBackStutter", "HeightUncappedForwardStutter", "HeightUncappedDoubleBackStutter", "HeightUncappedHalfStutter")
  short_names <- c("A", "S", "F", "S", "H")
  info_for_ANN_reader <- annotate_mix[column_names_to_target]
  #and now the amelogenin bits
  number_amel_alleles <- sum(one_profile$Locus == "AMEL")
  if (number_amel_alleles > 0){
    Amel_matrix <- matrix(one_profile$Height[one_profile$Locus == "AMEL"], ncol = 1)
    Amel_matrix <- cbind(Amel_matrix, matrix(rep(0, 8), ncol=4))
  }
  colnames(Amel_matrix) <- column_names_to_target
  info_for_ANN_reader <- rbind(info_for_ANN_reader, Amel_matrix)
  #now add on to profile
  one_profile_to_return <- cbind(one_profile, info_for_ANN_reader)
  #finally add the ANN call to the final column
  one_profile_to_return <- transform(one_profile_to_return, Assign = short_names[max.col(one_profile_to_return[,6:10])])
}
