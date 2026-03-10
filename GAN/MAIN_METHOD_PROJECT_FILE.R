#clears environment
rm(list=ls())

# Required R packages
# install.packages(c("here", "tidyverse", "ggplot2", "gridExtra", "reticulate"))
# TensorFlow/Keras setup
Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")
library(reticulate)
library(tensorflow)
library(keras)
library(here)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)

#sets booleans to determine what work needs to be done
#convert_EPG_xml_to_csv <- FALSE
pre_process_real_data <- TRUE
pre_process_synthetic_data <- FALSE
create_smoothed_profile_from_real_profile <- TRUE
carry_out_profile_generator_pretraining <- TRUE
carry_out_profile_discriminator_pretraining <- TRUE
create_profile_GAN_training_dataset <- TRUE
carry_out_profile_GAN_training <- TRUE
carry_out_additional_profile_discrimiantor_training_from_model <- FALSE
plot_all_predictions_of_training_profiles <- FALSE
simulate_a_profile <- TRUE

#filenames of Tensorflow models to load if pre-training is not occurring
#profile_generator_model_filename <- "Profile_generator_GAN_trained_model_random_24-09-2023.h5"
#profile_discriminator_model_filename <- "Profile_discriminator_GAN_trained_model_random_24-09-2023.h5"

#-----------------------------------------pre_process_real_data--------------------------------------------------------------------------------
if(pre_process_real_data){
  source("pre_process_data.R")
  #profile_dirs <- c(paste(here(), "input_profiles", sep="/"))
  profile_dirs <- "/home/emilie/Documents/Fil_rouge/Data/Test/Test/Data"
  pre_process_data(profile_dirs, real=TRUE)
}

if(pre_process_synthetic_data){
  source("pre_process_data.R")
  profile_dirs <- c(paste(here(), "synthetic_data/signals", sep="/"))
  pre_process_data(profile_dirs, real=FALSE)
}


#------------------------------------------------------------------------------------------------------------------------------------------
#profile_dirs <- c(paste(here(), "pre_process_data", sep="/"))
profile_dirs <- "/home/emilie/Documents/Fil_rouge/Data/Test/Test/pre_processed"
#fake_profile_dirs <- c(paste(here(), "prediction_generator", sep="/"))
fake_profile_dirs <- "/home/emilie/Documents/Fil_rouge/Data/Test/Test/fake"

#settings that are carried throughout the analysis
settings <- list(number_of_dyes = 3, #the number of dyes in the dataset
                 startScan = 1, #the starting scan point
                 number_of_scanpoints = 5000, #the number of scan points to be trained / generated
                 saturation = 1, #the saturation level of the instrument (for normalisation) - note the real saturation is 33000 but I am using a smaller value so that the prediction scale is larger (and the MSE loss encurs bigger penalties for not detecting peaks)
                 fake_profile_dirs = fake_profile_dirs,
                 real_and_smooth_dirs = profile_dirs,
                 saturation_max = 50000,
                 include_noise_generation = TRUE) 
#----------------------------------------- create_smoothed_profile_from_real_profile  ----------------------------------------------------

if(create_smoothed_profile_from_real_profile){
  #read in source
  source('smooth_profile.R')
  source('compare_raw_v_smoothed_from_scans.R')
  source('create_prediction_vector.R')
  # source('get_predictions_using_ANN.R')
  # the profile reading ANN, initialised as NULL and only loaded if needed
  # profile_reading_ANN <- NULL
  # get number of directories that hold files
  number_of_directories <- length(profile_dirs)
  print(number_of_directories)
  # loop over each directory
  for (iDir in 1:number_of_directories){
    #get names of cvs samples
    SampleNames <- list.files(path = profile_dirs[iDir], pattern="\\.csv$")
    # print(SampleNames)
    # print(head(SampleNames))
    # get number of samples
    NoOfSamples <- length(SampleNames)
    print(NoOfSamples)
    # loop over each sample
    for (iSample in 1:NoOfSamples){
      #gets the sample
      sampleName <- SampleNames[iSample]
      print(paste("pre processing sample ", sampleName, sep=""))
      
      # define config
      config <- list(
        
        # Path to your features CSV
        features_path = paste(profile_dirs[iDir], sampleName, sep="/"),
        
        # Path to your labels CSV
        labels_path    = paste(here(), "input_data/labels.csv", sep="/"), 
        #change path to make it automatic
        
        # Channels to process
        channels       = c("channel_1", "channel_3", "channel_4"),
        
        # Molecular weight window around each peak (in Da)
        # A scan is labelled "M" if its molw is within [mu - mw_window, mu + mw_window]
        mw_window      = 10,
        
        # Output path for the prediction vector
        #output_path    = paste(here(),"label_data", sampleName, sep="/")
        output_path    = paste("/home/emilie/Documents/Fil_rouge/Data/Test/Test/labels", sampleName,sep="/")
      )
      
      raw_profile <- read.table(config$features_path, sep=";", header=TRUE, row.names = NULL,fill = TRUE)
      print(dim(raw_profile))
      
      predictions <- create_predictions_vector(config)
      #plot_predictions_2(raw_profile, predictions)

      # smooths the profile data
      smoothed_profile <- smooth_profile(raw_profile, # the baselined profile
                                         predictions, #the predictions
                                         0) #threshold of 6 means 1 million to 1
     
      print("smoothing done")
      # plots a simple baselined vs smoothed (this shows more what the ANN will learn to emulate)
      sample_name_without_suffix <- substr(x = sampleName, start = 1, stop = nchar(sampleName) - 4)
      print(sample_name_without_suffix)
      print(dim(raw_profile))
    
      compare_raw_v_smoothed_from_scans(raw_profile,
                                        smoothed_profile,
                                        paste(profile_dirs[iDir], sample_name_without_suffix, sep="/")
      )
      # and saves the smoothed csv-format file
      pre_save_name <- paste(profile_dirs[iDir], sample_name_without_suffix, sep="/")
      print(pre_save_name)
      save_name <- paste(pre_save_name, "_smooth.csv", sep="")
      write.table(x = smoothed_profile, file = save_name, col.names = TRUE, row.names = FALSE, sep=";")
    } # next sample
  } # next directory
}


#----------------------------------------- carry_out_profile_generator_pretraining  ----------------------------------------------------
profile_dirs <- c(paste(here(), "pre_process_data", sep="/"))
fake_profile_dirs <- c(paste(here(), "prediction_generator", sep="/"))

#settings that are carried throughout the analysis
settings <- list(number_of_dyes = 3, #the number of dyes in the dataset
                 startScan = 1, #the starting scan point
                 number_of_scanpoints = 5000, #the number of scan points to be trained / generated
                 saturation = 1, #the saturation level of the instrument (for normalisation) - note the real saturation is 33000 but I am using a smaller value so that the prediction scale is larger (and the MSE loss encurs bigger penalties for not detecting peaks)
                 fake_profile_dirs = fake_profile_dirs,
                 real_and_smooth_dirs = profile_dirs,
                 saturation_max = 50000,
                 include_noise_generation = TRUE) 
if (carry_out_profile_generator_pretraining) {
  source("build_profile_generator.R")
  source("create_generator_pretrain_training_data.R")
  source("Profile_generator_pre_training.R")
  
  training_epochs <- 50L
  use_sample_weights <- FALSE
  save_profile_images <- TRUE
  save_profile_csvs <- TRUE
  
  profile_generator <- build_profile_generator(settings)
  print("profile_generator done (architecture)")
  
  profile_generator_training_data <- create_generator_pretrain_training_data(settings)
  print(" generator pre training dataset done")
  
  profile_generator <- profile_generator_pre_training(settings,
                                 profile_generator,
                                 training_epochs,
                                 save_profile_images,
                                 save_profile_csvs,
                                 use_sample_weights,
                                 profile_generator_training_data)
  
  print("profile generator pre training done (training)")
  
  # create output directory if needed
  saved_models_dir <- file.path(here(), "saved_models", sep="/")
  print(saved_models_dir)
  if (!dir.exists(saved_models_dir)) {
    dir.create(saved_models_dir, recursive = TRUE)
  }
  
  # date string (base R)
  date_string <- format(Sys.Date(), "%Y-%m-%d")
  
  # save models
  save_model_hdf5(
    profile_generator,
    file.path(
      saved_models_dir,
      paste0("Profile_generator_pretrained_model_NO_noise_", date_string, ".h5")
    )
  )
}

#----------------------------------------- carry_out_profile_discriminator_pretraining  ----------------------------------------------------

if(carry_out_profile_discriminator_pretraining){
  #build the tensorflow model
  source('build_profile_discriminator.R')
  source('create_discriminator_pretrain_training_data.R')
  source('Profile_discriminator_pre_training.R')
  
  useL1Regularization <- FALSE
  
  profile_discriminator <- build_profile_discriminator(settings, useL1Regularization)
  print("architecture of profile discriminator done")
  
  profile_discriminator_training_data <- create_discriminator_pretrain_training_data(settings)
  print("training data for profile discriminator done")
  
  #assign settings
  number_epochs <- 10L
  
  #carry out the discriminator pre-training
  profile_discriminator <- profile_discriminator_pre_training(settings,
                                                              profile_discriminator,
                                                              number_epochs,
                                                              profile_discriminator_training_data)
  print("training of profile discriminator done")
  
  # directory for saved models
  saved_models_dir <- file.path(here(), "saved_models", sep="/")
  
  # create directory if needed
  if (!dir.exists(saved_models_dir)) {
    dir.create(saved_models_dir, recursive = TRUE)
  }
    
  # save the TensorFlow model
  date_string <- format(Sys.Date(), "%Y-%m-%d")
  
  save_model_hdf5(
    profile_discriminator,
    file.path(
      saved_models_dir,
      paste0(
        "Profile_discriminator_pretrained_model_",
        date_string,
        ".h5"
      )
    )
  )
}
#----------------------------------------- carry_out_profile_GAN_training_no_noise  ----------------------------------------------------
profile_dirs <- c(paste(here(), "pre_process_data/real", sep="/"))
fake_profile_dirs <- c(paste(here(), "prediction_generator", sep="/"))

#settings that are carried throughout the analysis
settings <- list(number_of_dyes = 3, #the number of dyes in the dataset
                 startScan = 1, #the starting scan point
                 number_of_scanpoints = 5000, #the number of scan points to be trained / generated
                 saturation = 1, #the saturation level of the instrument (for normalisation) - note the real saturation is 33000 but I am using a smaller value so that the prediction scale is larger (and the MSE loss encurs bigger penalties for not detecting peaks)
                 fake_profile_dirs = fake_profile_dirs,
                 real_and_smooth_dirs = profile_dirs,
                 saturation_max = 50000,
                 include_noise_generation = TRUE) 
#que sur données réelles !!!!
if (carry_out_profile_GAN_training) {
  
  source("create_gan_training_data.R")
  source("Profile_GAN_trainer.R")
  
  print("Sarting training dataset")
  GAN_training_dataset <- create_gan_training_data(settings)
  print("create GAN dataset done")
  training_epochs <- 30L
  use_sample_weights <- TRUE
  
  # To adapt
  profile_generator_path <- file.path(here(), "saved_models",
                                      "Profile_generator_pretrained_model_NO_noise_2026-02-25.h5")
  # To adapt
  profile_discriminator_path <- file.path(here(), "saved_models",
                                          "Profile_discriminator_pretrained_model_2026-02-25.h5")
  
  
  profile_generator <- load_model_hdf5(profile_generator_path)
  profile_discriminator <- load_model_hdf5(profile_discriminator_path)
  
  # Training GAN
  updated_gen_disc_ANNS <- profile_GAN_trainer(
    settings,
    profile_discriminator,
    profile_generator,
    use_sample_weights,
    training_epochs,
    GAN_training_dataset
  )
  
  profile_generator <- updated_gen_disc_ANNS$updated_profile_generator
  profile_discriminator <- updated_gen_disc_ANNS$updated_profile_discriminator
  
  # Final save
  date_string <- format(Sys.Date(), "%d-%m-%Y")
  
  save_model_hdf5(profile_generator,
                  file.path(here(), "saved_models",
                            paste0("Profile_generator_GAN_trained_model_NO_noise_", date_string, ".h5")))
  
  save_model_hdf5(profile_discriminator,
                  file.path(here(), "saved_models",
                            paste0("Profile_discriminator_GAN_trained_model_NO_noise_", date_string, ".h5")))
}


#---------------------------------plot the results----------------------------------------------------------------------------------

plot_generated_vs_input <- function(smooth_profile_df, 
                                    generated_profile_df, 
                                    number_of_dyes,
                                    savename = NULL) {
  library(ggplot2)
  library(gridExtra)

  # Common y axis limits across both profiles
  plot_max_rfu <- max(
    max(as.numeric(unlist(smooth_profile_df[, 2:(number_of_dyes + 1)]))),
    max(as.numeric(unlist(generated_profile_df[, 2:(number_of_dyes + 1)])))
  )
  plot_min_rfu <- min(
    min(as.numeric(unlist(smooth_profile_df[, 2:(number_of_dyes + 1)]))),
    min(as.numeric(unlist(generated_profile_df[, 2:(number_of_dyes + 1)])))
  )
  
  plot_list <- vector('list', number_of_dyes)
  
  for (dye in 1:number_of_dyes) {
    
    col_name <- paste0("channel_", dye)
    
    p <- ggplot() +
      # raw/smooth input in black
      geom_line(data = smooth_profile_df,
                aes_string(x = "molw", y = col_name),
                col = "black", linewidth = 0.2, alpha = 0.7) +
      # generated output in red
      geom_line(data = generated_profile_df,
                aes_string(x = "molw", y = col_name),
                col = "red", linewidth = 0.2, alpha = 0.7) +
      ylim(plot_min_rfu, plot_max_rfu) +
      ylab("RFU") +
      xlab("") +
      ggtitle(paste0("Channel ", dye)) +
      theme_classic() +
      theme(plot.title = element_text(size = 8),
            axis.text.x = element_text(size = 5))
    
    plot_list[[dye]] <- p
  }
  
  # Add a shared x label at the bottom
  full_plot <- do.call("grid.arrange", c(plot_list, ncol = 1,
                                         bottom = "molecular weight (pdb)",
                                         top = "Black = smooth input   |   Red = generated"))
  
  if (!is.null(savename)) {
    ggsave(savename, full_plot, units = 'px', width = 2000, height = number_of_dyes * 500)
  }
  
  return(full_plot)
}

#---------------------------------- simulate random profiles, make them look real and save them ------------------------------------

if (simulate_a_profile){
  # Load an existing smooth profile as the generator input
  smooth_profile_raw <- read.table(file = "/home/emilie/Documents/Fil_rouge/Code/Rproj/pre_process_data/M1_pl0035_smooth.csv", 
                                 header = TRUE, stringsAsFactors = FALSE, sep=";")
  molw <- smooth_profile_raw$molw
  smooth_profile <- t(smooth_profile_raw[settings$startScan:settings$number_of_scanpoints, 
                                         2:(settings$number_of_dyes + 1)])
  
  profile_generator_path <- file.path(here(), "saved_models",
                                      "Profile_generator_GAN_trained_model_NO_noise_25-02-2026.h5")
  profile_generator <- profile_generator <- load_model_hdf5(profile_generator_path)
  # Numbered input
  if (!settings$include_noise_generation) {
    numbered_input            <- matrix(rep(1:500, settings$number_of_dyes), ncol = 500, byrow = TRUE) / 500
    numbered_input_one_sample <- array(numbered_input, dim = c(1, settings$number_of_dyes, 500, 1))
  }
  
  # Run through generator — each call produces a different realistic-looking profile
  # from the same smooth input (due to noise input variation)
  for (iSim in 1:2) {
    
    if (settings$include_noise_generation) {
      numbered_input_one_sample <- array(runif(settings$number_of_dyes * 500),
                                         dim = c(1, settings$number_of_dyes, 500, 1))
    }
    
    smooth_profile_array  <- array(smooth_profile, 
                                   dim = c(1, settings$number_of_dyes,settings$number_of_scanpoints, 1)) / settings$saturation
    generated_profile_raw <- profile_generator$predict(list(smooth_profile_array, numbered_input_one_sample))
    
    # Format and save
    generated_profile_df <- as.data.frame(
      cbind(seq_len(settings$number_of_scanpoints), 
            t(generated_profile_raw[1, , , 1]) * settings$saturation)
    )
    generated_profile_df <- cbind(generated_profile_df,molw)
    colnames(generated_profile_df) <- c("index", paste0("channel_", seq_len(settings$number_of_dyes)), "molw")
    
    write.table(generated_profile_df, 
              file.path(here(), "simulated_profiles", paste0(iSim, "_generated.csv")),
              row.names = FALSE, sep=";")
    # After saving the CSV, also plot
    smooth_profile_df <- as.data.frame(
      cbind(seq_len(settings$number_of_scanpoints), t(smooth_profile))
    )
    smooth_profile_df <- cbind(smooth_profile_df, molw)
    colnames(smooth_profile_df) <- c("index", paste0("channel_", seq_len(settings$number_of_dyes)), "molw")
    
    plot_generated_vs_input(
      smooth_profile_df    = smooth_profile_df,
      generated_profile_df = generated_profile_df,
      number_of_dyes       = settings$number_of_dyes,
      savename             = file.path(here(), "simulated_profiles", paste0(iSim, "_comparison.jpg"))
    )
  }
}
# ----------------------------------------------------------------------------------------------------------------------------------------