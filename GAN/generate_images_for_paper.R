library('here')
library('tensorflow')
library('keras')
library('tidyverse')

source('baseline_raw_profile.R')
source('baseline_raw_profile_simple.R')
source('peak_detection.R')
source('smooth_profile.R')
source('get_predictions_using_ANN.R')
source('create_gan_training_data.R')

#----------------------------------------------- Figure 1 -------------------------------------------------

plot_Fig1 <- FALSE

if(plot_Fig1){

savename <- paste(here(), "\\paper_Figure1.jpg", sep="")
jpeg(file=savename, width = 4000, height=6000, res=600)
  #set up the 5 panels (b,l,t,r)
  par(mfrow=c(5, 1), mar=c(2, 2, 0, 0), oma=c(2,2,2,2))
  
  startScan <- 4000
  endScan <- 9000
  plot_dye <- 1
  
  #load in DNA profile
  profile_path <- paste(here(), "\\ProvedIT_profiles\\PROVEDIt_2-5-Person Profiles_3500 5sec_GF29cycles\\5 sec\\RD14-0003(020316ADG_5sec)\\A02_RD14-0003-31_32-1;1-M2c-0.03GF-Q2.0_01.5sec.csv", sep="")
  profile_data <- read.csv(file = profile_path, header = TRUE, stringsAsFactors = FALSE)
  
  ################A) the raw profile plot
  baselined_profle <- baseline_raw_profile_simple(profile_data)
  
  yplot_max <- max(baselined_profle[startScan:endScan, plot_dye+1])
  yplot_min <- min(baselined_profle[startScan:endScan, plot_dye+1])
  
  plot(x = baselined_profle$scan[startScan:endScan],
       y = baselined_profle[startScan:endScan, plot_dye+1],
       ylim = c(yplot_min, yplot_max),
       type = 'l',
       xlab = "",
       ylab="")
  text(x = startScan, y = yplot_max*0.9, "A")
  
  ###############B) lowess line
  lowess_smooth <- lowess(y = baselined_profle[,plot_dye+1], x = baselined_profle[,1], f = 0.05)
  
  plot(x = lowess_smooth$x[startScan:endScan],
       y = lowess_smooth$y[startScan:endScan],
       ylim = c(yplot_min, yplot_max),
       type = 'l',
       xlab = "",
       ylab="")
  text(x = startScan, y = yplot_max*0.9, "B")
  
  ##################C) baselined using lowess
  baselines_profle_lowess <- baseline_raw_profile(baselined_profle)
  
  plot(x = baselines_profle_lowess$scan[startScan:endScan],
       y = baselines_profle_lowess[startScan:endScan, plot_dye+1],
       ylim = c(yplot_min, yplot_max),
       type = 'l',
       xlab = "",
       ylab="")
  text(x = startScan, y = yplot_max*0.9, "C")
  
  ###############D) peak centres
  probability_threshold <- 10
  peak_RFU_threshold <- 20
  idye <- plot_dye + 1
  
  #identify peaks at certain thresholds and store position and height
  peak_probs_for_dye <- peak_detection(baselines_profle_lowess[,idye], peak_RFU_threshold)
  
  #screen out non-peaks below probability threshold
  thinned_peak_probs_for_dye <- peak_probs_for_dye
  thinned_peak_probs_for_dye[thinned_peak_probs_for_dye < probability_threshold] <- 0
  #use a min-max function so that each peak is left with one central scan point
  scan_window <- 10
  
  for (iScan in 1:(length(thinned_peak_probs_for_dye) - scan_window)){
    max_prob_in_scan_window <- max(thinned_peak_probs_for_dye[iScan:(iScan + scan_window)])
    indices_of_non_maxes <- which(thinned_peak_probs_for_dye[iScan:(iScan + scan_window)] != max_prob_in_scan_window)
    thinned_peak_probs_for_dye[iScan:(iScan + scan_window)][indices_of_non_maxes] <- 0
  }
  
  #plot(thinned_peak_probs_for_dye[startScan:endScan], type='l')
  
  #screen out values prior ro primer flare
  thinned_peak_probs_for_dye[1:4000] <- 0
  #screen out any peak centres for which profile has below set fluorescence
  thinned_peak_probs_for_dye[baselines_profle_lowess[,idye] < peak_RFU_threshold] <- 0
  
  #screen out any peaks that have a prediction type of P (pull-up) or B (baseline) or X (unassigned)
  profile_reading_ANN <- load_model_hdf5(paste(here(), "\\profile_reading_ANN.h5", sep=""))
  raw_profile_predictions <- cbind(profile_data[,1], get_predictions_using_ANN(profile_data, profile_reading_ANN))
  
  thinned_peak_probs_for_dye[raw_profile_predictions[,idye] == 'P' |
                               raw_profile_predictions[,idye] == 'B' |
                               raw_profile_predictions[,idye] == 'B'] <- 0
  
  peak_centres_at_indices <- which(thinned_peak_probs_for_dye > 0)
  
  plot(x = 0,
       y = 0,
       col = "white",
       xlim = c(startScan, endScan),
       ylim= c(0, 1),
       xlab = "",
       ylab="",
       yaxt='n')
  for (iPeak in 1:length(peak_centres_at_indices)){
    abline(v = peak_centres_at_indices[iPeak], lwd=2)
  }
  text(x = startScan, y = 0.9, "D")
  
  #################E) final smoothed profile
  #smooths the profile data
  smoothed_profile <- smooth_profile(baselines_profle_lowess, # the baselined profile
                                     raw_profile_predictions, #the predictions
                                     10) #threshold of 6 means 1 million to 1
  
  plot(x = smoothed_profile$scan[startScan:endScan],
       y = smoothed_profile[startScan:endScan, plot_dye+1],
       ylim = c(yplot_min, yplot_max),
       type = 'l',
       xlab = "",
       ylab="")
  text(x = startScan, y = yplot_max*0.9, "E")
  
  #adding axis labels
  mtext("Scan point", side=1, line=1, cex=1, col="black", outer=TRUE)
  mtext("fluorescence", side=2, line=1, cex=1, col="black", outer=TRUE)
  #reste defaults
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), oma=c(0,0,0,0))
dev.off()

}

#----------------------------------------------- Figure 5 -------------------------------------------------

plot_Fig5 <- FALSE

if(plot_Fig5){
  
  relaculate_values <- FALSE
  
  if (relaculate_values){

    #path to saved models and loss 
    results_path <- paste(here(), "\\GAN_increments\\1078 profiles success - random\\", sep="")
    
    ######### A) loss plots
    loss_data <- read.csv(file=paste(results_path, "Disc_and_Gen_loss_full_run.csv", sep=""), header = TRUE, stringsAsFactors = FALSE)
    
    fig5_yplot_min <- min(loss_data$disc_loss, loss_data$gen_loss)
    fig5_yplot_max <- max(loss_data$disc_loss, loss_data$gen_loss)
    
    ##### B) discriminator performance boxplots
    
    main_profile_path <- paste(here(), "\\ProvedIT_profiles\\PROVEDIt_2-5-Person Profiles_3500 5sec_GF29cycles\\5 sec\\", sep="")
    original_profile_dirs <- c(paste(main_profile_path, "RD14-0003(020316ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(020516ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(021016ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(021516ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(021716ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(022316ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(022516ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(022616ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(022916ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(030216ADG_5sec)\\", sep=""),
                               paste(main_profile_path, "RD14-0003(030716ADG_5sec)\\", sep="")
    )
    fake_profile_dirs <- paste(here(), "\\prediction_generator\\", sep="")
    
    additonal_profile_dirs_path <- paste(here(), "\\input_profiles\\", sep="")
    additonal_profile_dirs <- c(paste(additonal_profile_dirs_path, "GF_3500_RUN_14975\\", sep=""),
                                paste(additonal_profile_dirs_path, "GF_3500_RUN_14999\\", sep=""),
                                paste(additonal_profile_dirs_path, "GF_3500_RUN_15015\\", sep=""),
                                paste(additonal_profile_dirs_path, "GF_3500_RUN_15015 - NOT YET\\", sep=""),
                                paste(additonal_profile_dirs_path, "GF_3500_RUN_15022\\", sep=""),
                                paste(additonal_profile_dirs_path, "GlobalFiler 3500 mixtures\\", sep=""),
                                paste(additonal_profile_dirs_path, "GlobalFiler 3500 refs\\", sep=""),
                                paste(additonal_profile_dirs_path, "GlobalFiler_3500_2021_lowlevel\\", sep=""),
                                paste(additonal_profile_dirs_path, "GlobalFiler_3500_PROVEDIt\\", sep=""),
                                paste(additonal_profile_dirs_path, "GlobalFiler_3500_PROVEDIt_mixtures\\", sep="")
    )
    
    profile_dirs <- c(original_profile_dirs, additonal_profile_dirs)
    
    settings <- list(number_of_dyes = 6, #the number of dyes in the dataset
                     startScan = 4000, #the starting scan point
                     number_of_scanpoints = 5000, #the number of scan points to be trained / generated
                     saturation = 100, #the saturation level of the instrument (for normalisation) - note the real saturation is 33000 but I am using a smaller value so that the prediction scale is larger (and the MSE loss encurs bigger penalties for not detecting peaks)
                     fake_profile_dirs = fake_profile_dirs,
                     real_and_smooth_dirs = profile_dirs,
                     saturation_max = 33000,
                     fluorescence_column_names = paste("dye", 1:6, sep=""),
                     include_noise_generation = TRUE) 
    
    #creates the GAN training and target data
    GAN_training_dataset <- create_gan_training_data(settings)
    
    number_increments <- 21
    number_profiles <- dim(GAN_training_dataset$real_profile_target)[1]
    
    real_predictions <- matrix(NA, ncol = number_increments, nrow = number_profiles)
    fake_predictions <- matrix(NA, ncol = number_increments, nrow = number_profiles)
    
    for (iInc in 1:number_increments){
      #indicator
      cat(paste("starting round ", iInc, " / ", number_increments, "\n", sep=""))
      #epoch in filename to load
      epoch_increment <- (iInc - 1)*10
      #load discriminator
      dicriminiator_increment_filename <- paste("Profile_discriminator_GAN_trained_model_epoch_", epoch_increment, ".h5", sep="")
      discriminator <-load_model_hdf5(paste(results_path, dicriminiator_increment_filename, sep=""))
    
      #load generator
      generator_increment_filename <- paste("Profile_generator_GAN_trained_model_epoch_", epoch_increment, ".h5", sep="")
      generator <-load_model_hdf5(paste(results_path, generator_increment_filename, sep=""))
      
      #regenerate real and fake profiles
      source('create_discriminator_pretrain_training_data_from_models.R')
      profile_discriminator_training_data <- create_discriminator_pretrain_training_data_from_models(settings,
                                                                                                     discriminator,
                                                                                                     generator,
                                                                                                     settings$include_noise_generation,
                                                                                                     GAN_training_dataset)
      #go through each profile and classify
      for(iProfile in 1:number_profiles){
        even_spot <- iProfile*2 #fake profiles
        odd_spot <- even_spot - 1 #real profiles
        #the conditional smoothed input
        smoothed_profile <- array(profile_discriminator_training_data$src_images_smooth[odd_spot,,,], dim = c(1, 6, 5000, 1))
        #real profile
        real_profile <- array(profile_discriminator_training_data$targets_profiles_real_and_fake[odd_spot,,,], dim = c(1, 6, 5000, 1))
        discriminator_prediction_real <- discriminator %>% predict(list(smoothed_profile, real_profile), verbose = FALSE)
        real_predictions[iProfile, iInc] <- mean(discriminator_prediction_real)
        #fake profile
        fake_profile <- array(profile_discriminator_training_data$targets_profiles_real_and_fake[even_spot,,,], dim = c(1, 6, 5000, 1))
        discriminator_prediction_fake <- discriminator %>% predict(list(smoothed_profile, fake_profile), verbose = FALSE)
        fake_predictions[iProfile, iInc] <- mean(discriminator_prediction_fake)
      }
      
    }
    #save the results (because this takes quite a while to generate)
    write.table(x = real_predictions, file=paste(results_path, "real_prediction_increments.csv", sep=""), sep=",", col.names = TRUE, row.names = FALSE)
    write.table(x = fake_predictions, file=paste(results_path, "generated_prediction_increments.csv", sep=""), sep=",", col.names = TRUE, row.names = FALSE)
  }
  
  if (settgings$include_noise_generation){
    savename <- paste(here(), "\\paper_Figure5 - random.jpg", sep="")
  } else {
    savename <- paste(here(), "\\paper_Figure5 - structured.jpg", sep="")
  }
  jpeg(file=savename, width = 6000, height=4000, res=450)
    par(mfrow=c(2,1), mar=c(3, 4, 0.1, 0.1)) #(b,l,t,r)
    
    plot(x = loss_data$epoch,
         y = loss_data$disc_loss,
         ylim = c(fig5_yplot_min, fig5_yplot_max),
         xlab = "epoch",
         ylab = "loss",
         col="black",
         type='l',
         lwd=2)
    lines(x = loss_data$epoch,
          y = loss_data$gen_loss,
          col="black",
          lwd=2,
          lty=2)
    legend("topright", legend=c("discriminator", "generator"), lty=c(1,2), col="black", lwd=2)
    
    #now create the boxplot
    plot(x = 0, y = 0, col="white",
         xlim=c(1, number_increments),
         ylim=c(0, 1.2),
         ylab="discriminator prediction",
         xlab="",
         xaxt='n',
         yaxt='n')
  
    for (iInc in 1:number_increments){
      boxplot(x = fake_predictions[,iInc], at = iInc - 0.2, col="grey30", yaxt='n', add=TRUE)
      boxplot(x = real_predictions[,iInc], at = iInc + 0.2, col="white", yaxt='n', add=TRUE)
      abline(v = iInc - 0.5, col="grey30")
    }
  
    mtext(side = 1, line = 0.5, text = seq(from = 0, to = number_increments*10, by = 10), at = 1:number_increments)
    mtext(side = 1, line = 1.5, text = "epoch")
    
    mtext(side = 2, line = 1, text = seq(from = 0, to = 1, by = 0.1), at = seq(from = 0, to = 1, by = 0.1), las=2)
      
    legend("topleft", fill=c("grey30", "white"), legend=c("generated", "real"))
    
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  dev.off()
}

#----------------------------------------------- Figure 6 -------------------------------------------------

plot_Fig6 <- FALSE

if(plot_Fig6){
#  profile_path <- paste(here(), "\\ProvedIT_profiles\\PROVEDIt_2-5-Person Profiles_3500 5sec_GF29cycles\\5 sec\\RD14-0003(020516ADG_5sec)\\A06_RD14-0003-44_45_46-1;2;2-M3a-0.075GF-Q0.6_01.5sec.csv", sep="")
#  smooth_profile_path <- paste(here(), "\\ProvedIT_profiles\\PROVEDIt_2-5-Person Profiles_3500 5sec_GF29cycles\\5 sec\\RD14-0003(020516ADG_5sec)\\A06_RD14-0003-44_45_46-1;2;2-M3a-0.075GF-Q0.6_01.5sec_smooth.csv", sep="")
  
  profile_path <- paste(here(), "\\ProvedIT_profiles\\PROVEDIt_2-5-Person Profiles_3500 5sec_GF29cycles\\5 sec\\RD14-0003(020516ADG_5sec)\\D01_RD14-0003-36_37_38-1;2;1-M3a-0.124GF-Q0.7_04.5sec.csv", sep="")
  smooth_profile_path <- paste(here(), "\\ProvedIT_profiles\\PROVEDIt_2-5-Person Profiles_3500 5sec_GF29cycles\\5 sec\\RD14-0003(020516ADG_5sec)\\D01_RD14-0003-36_37_38-1;2;1-M3a-0.124GF-Q0.7_04.5sec_smooth.csv", sep="")
  
  
  
  
  profile <- read.csv(file = profile_path, header = TRUE, stringsAsFactors = FALSE)
  smooth_profile <- read.csv(file = smooth_profile_path, header = TRUE, stringsAsFactors = FALSE)
  
  smooth_profile_for_gen <- t(smooth_profile[settings$startScan:(settings$startScan + settings$number_of_scanpoints - 1),2:(number_of_dyes + 1)])
  smooth_profile_for_gen_array <- array(smooth_profile_for_gen/settings$saturation, dim=c(1, 6, 5000, 1))
  
  generator_increment_filename <- paste(here(), "\\saved_models\\Profile_generator_GAN_trained_model_random_24-09-2023.h5", sep="")
  generator <-load_model_hdf5(generator_increment_filename)
  
  if (!settings$include_noise_generation){
    numbered_input_array <- array(runif(settings$number_of_dyes*500), dim=c(1, settings$number_of_dyes, 500, 1))
  }  else {
    numbered_input <- matrix(rep(1:500, 6), ncol=500, byrow = TRUE)/500
    numbered_input_array <- array(NA, dim=c(1, 6, 500, 1))
  }
  generated_profile <- generator %>% predict(list(smooth_profile_for_gen_array, numbered_input_array))
  
  plot_dye <- 3
  
  savename <- paste(here(), "\\paper_Figure6.jpg", sep="")
  jpeg(file=savename, width = 4000, height=4000, res=450)
    par(mfrow=c(3,1), mar=c(2,0,0,0), oma = c(3,3,0.1,0.1)) #(b,l,t,r)
    
    plot_start <- 4000
    plot_end <- 5000
    
    # the real profile
    plot(profile[plot_start:plot_end, plot_dye + 1],
         yaxt='n',
          ylab="",
         xlab="",
          type='l')
    # the smooth profile
    plot(smooth_profile[plot_start:plot_end, plot_dye + 1],
         yaxt='n',
         ylab="",
         xlab="",
         type='l')
    # the generated profiles
    plot(generated_profile[1, plot_dye, (plot_start - settings$startScan + 1):(plot_end - plot_start), 1],
         yaxt='n',
         ylab="",
         xlab="",
         type='l')
    
    mtext(side = 1, line = 1, text="scan point", outer = TRUE)
    mtext(side = 2, line = 1, text="generated profile (RFU)", at = 0.15, outer = TRUE)
    mtext(side = 2, line = 1, text="smoothed profile (RFU)", at = 0.5, outer = TRUE)
    mtext(side = 2, line = 1, text="real profile (RFU)", at = 0.85, outer = TRUE)
  
  dev.off()
  
}

#----------------------------------------------- Figure 7 -------------------------------------------------

plot_Fig7 <- FALSE

if(plot_Fig7){
  number_of_dyes <- 6
  #load profile
  strong_profile_path <- paste(here(), "\\ProvedIT_profiles\\PROVEDIt_2-5-Person Profiles_3500 5sec_GF29cycles\\5 sec\\RD14-0003(020316ADG_5sec)\\F07_RD14-0003-35_50-1;9-M2a-0.75GF-Q0.7_06.5sec.csv", sep="")
  strong_smooth_profile_path <- paste(here(), "\\ProvedIT_profiles\\PROVEDIt_2-5-Person Profiles_3500 5sec_GF29cycles\\5 sec\\RD14-0003(020316ADG_5sec)\\F07_RD14-0003-35_50-1;9-M2a-0.75GF-Q0.7_06.5sec_smooth.csv", sep="")

  strong_profile <- read.csv(file = strong_profile_path, header = TRUE, stringsAsFactors = FALSE)
  strong_smooth_profile <- read.csv(file = strong_smooth_profile_path, header = TRUE, stringsAsFactors = FALSE)
  
  strong_smooth_profile_for_gen <- t(strong_smooth_profile[settings$startScan:(settings$startScan + settings$number_of_scanpoints - 1),2:(number_of_dyes + 1)])
  strong_smooth_profile_for_gen_array <- array(strong_smooth_profile_for_gen/settings$saturation, dim=c(1, 6, 5000, 1))
  
  generator_increment_filename <- paste(here(), "\\saved_models\\Profile_generator_GAN_trained_model_random_24-09-2023.h5", sep="")
  generator <-load_model_hdf5(generator_increment_filename)
  
  if (settings$include_noise_generation){
    numbered_input_array <- array(runif(number_of_dyes*500), dim=c(1, number_of_dyes, 500, 1))
  } else {
    numbered_input <- matrix(rep(1:500, 6), ncol=500, byrow = TRUE)/500
    numbered_input_array <- array(NA, dim=c(1, 6, 500, 1))
  }
  generated_profile <- generator %>% predict(list(strong_smooth_profile_for_gen_array, numbered_input_array))
  
  savename <- paste(here(), "\\paper_Figure7.jpg", sep="")
  jpeg(file=savename, width = 4000, height=4000, res=450)
    par(mfcol=c(6,3), mar=c(0,0,0,0), oma = c(3,3,3,0.1)) #(b,l,t,r)
    #real (left)
    for (iDye in 1:number_of_dyes){
      plot(strong_profile[4500:5500, iDye+1],
           xaxt='n',
           yaxt='n',
           type='l')
    }
    #smooth (mid)
    for (iDye in 1:number_of_dyes){
      plot(strong_smooth_profile_for_gen_array[1, iDye, 500:1500, 1],
           xaxt='n',
           yaxt='n',
           type='l')
    }
    #generated (right)
    for (iDye in 1:number_of_dyes){
      plot(generated_profile[1, iDye, 500:1500, 1],
           xaxt='n',
           yaxt='n',
           type='l')
    }
    
    mtext("Scan point", side=1, line=1, cex=1, col="black", outer=TRUE)
    mtext("fluorescence", side=2, line=1, cex=1, col="black", outer=TRUE)
    mtext("Original profile", side=3, line=1, at = 0.15, cex=1, col="black", outer=TRUE)
    mtext("smoothed profile", side=3, line=1, at = 0.5, cex=1, col="black", outer=TRUE)
    mtext("generated profile", side=3, line=1, at = 0.85, cex=1, col="black", outer=TRUE)
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), oma=c(0,0,0,0))
  dev.off()
}

#----------------------------------------------- Figure 8 -------------------------------------------------

plot_Fig8 <- FALSE

if(plot_Fig8){ 
  library('lattice')

  profiles <- c("1 - 449753.xml.csv", "2 - 449755.xml.csv", "3 - 449746.xml.csv", "4 - 449756.xml.csv", "5 - 449742.xml.csv")
  
  for (iProfile in 1:length(profiles)){
  
    profile_path <- paste(here(), "\\test_profiles\\", profiles[iProfile], sep="")
    profile_raw_data <- read.csv(file = profile_path, header = TRUE, stringsAsFactors = FALSE)
    
    #profile reading ANN
    profile_reading_ANN <- load_model_hdf5(paste(here(), "\\profile_reading_ANN.h5", sep=""))
    
    #classify
    classification_on_real_profile <- cbind(profile_raw_data[,1], get_predictions_using_ANN(profile_raw_data, profile_reading_ANN))
    
    #smooth, regen and reclassify
    profile_for_gen <- t(profile_raw_data[settings$startScan:(settings$startScan + settings$number_of_scanpoints - 1),2:(number_of_dyes + 1)])
    profile_for_gen_array <- array(profile_for_gen/settings$saturation, dim=c(1, 6, 5000, 1))
    
    generator_increment_filename <- paste(here(), "\\saved_models\\Profile_generator_GAN_trained_model_random_24-09-2023.h5", sep="")
    generator <-load_model_hdf5(generator_increment_filename)
    
    if (settings$include_noise_generation){
      numbered_input_array <- array(runif(settings$number_of_dyes*500), dim=c(1, settings$number_of_dyes, 500, 1))
    } else {
      numbered_input <- matrix(rep(1:500, 6), ncol=500, byrow = TRUE)/500
      numbered_input_array <- array(NA, dim=c(1, 6, 500, 1))
    }
    
    baselined_profle <- baseline_raw_profile_simple(profile_raw_data)
    baselines_profle_lowess <- baseline_raw_profile(baselined_profle)
    smoothed_profile <- smooth_profile(baselines_profle_lowess, # the baselined profile
                                       classification_on_real_profile, #the predictions
                                       10) #threshold of 6 means 1 million to 1
    
    
    smooth_profile_for_gen <- t(smoothed_profile[settings$startScan:(settings$startScan + settings$number_of_scanpoints - 1),2:(number_of_dyes + 1)])
    smooth_profile_for_gen <- array(smooth_profile_for_gen/settings$saturation, dim=c(1, 6, 5000, 1))
    
    generated_profile <- generator %>% predict(list(smooth_profile_for_gen, numbered_input_array))
    
    generated_profile_padded <- rbind(matrix(0, ncol=number_of_dyes, nrow = settings$startScan - 1), t(generated_profile[1,,,1]), matrix(0, ncol=number_of_dyes, nrow = 10000 - (settings$startScan + settings$number_of_scanpoints)))
    scan_number_sequence <- seq(from = 1, to = dim(generated_profile_padded)[1], by = 1)
    generated_profile_for_classifier <- cbind(scan_number_sequence, generated_profile_padded*settings$saturation)
    classification_on_generated_profile <- cbind(scan_number_sequence, get_predictions_using_ANN(generated_profile_for_classifier, profile_reading_ANN))
    
    #create comparison array
    prediction_before_vs_after <- matrix(NA, ncol = 3, nrow = number_of_dyes*settings$number_of_scanpoints)
    colnames(prediction_before_vs_after) <- c("dye", "raw_classify", "gen_classify")
    for (iDye in 1:number_of_dyes){
      row_start <- (iDye - 1)*5000 + 1
      row_end <- iDye*5000
      prediction_before_vs_after[row_start:row_end,1] <- iDye
      prediction_before_vs_after[row_start:row_end,2] <- classification_on_real_profile[settings$startScan:(settings$startScan + settings$number_of_scanpoints - 1),iDye+1]
      prediction_before_vs_after[row_start:row_end,3] <- classification_on_generated_profile[settings$startScan:(settings$startScan + settings$number_of_scanpoints - 1),iDye+1]
    }
    
    
    comaprison_table <- table(prediction_before_vs_after[,2], prediction_before_vs_after[,3])
    number_cats <- dim(comaprison_table)[1]
    cat_names <- row.names(comaprison_table)
    comaprison_table_proportions <- prop.table(comaprison_table, 1)
    
  
    #now colour the table (using lattice) and add text in the cells showing the value (need to set up the coordinate arrays)
    x_array <- rep(1:number_cats, number_cats)
    y_array <- rep(1, number_cats)
    for (i in 2:number_cats){
      y_array <- c(y_array, rep(i, number_cats))
    }
    #then do the plot with a tricky bit of code to get the labels in there
    saveDir <- paste(here(), "\\test_profiles\\", sep="")
    jpeg(file=paste(saveDir, "Figure8_", iProfile,".jpg", sep=""), width = 2000, height=2000, res=300)
      levelplot(comaprison_table_proportions,
#                col.regions = heat.colors(100),
                col.regions = gray(100:0/100),
                panel = function(...) {
                   arg <- list(...)
                   panel.levelplot(...)
                   panel.text(y_array, x_array, t(comaprison_table))
                },
               colorkey=FALSE,
               xlab="",
               ylab="")
    dev.off()
     
  }

  
#  grid.arrange(grobs = lattice_plots, ncol = 2)
  
}

#----------------------------------------------- Figure 9 -------------------------------------------------

plot_Fig9 <- FALSE

if (plot_Fig9){

  generator_increment_filename <- paste(here(), "\\saved_models\\Profile_generator_GAN_trained_model_random_24-09-2023.h5", sep="")
  generator <-load_model_hdf5(generator_increment_filename)
  
  #the number of profiles to simulate
  number_of_profiles_to_simulate <- 100
  
  #load in scripts
  source('simulate_profiles.R')
  source('get_mock_EPG.R')
  source('get_mock_EPG_from_scans.R')
  #profile simulation settings
  profile_simulation_settings <- list(min_template = 20,   # note that these can be arrays if you want to create specific mixture proportions
                                      max_template = 100,       # note that these can be arrays if you want to create specific mixture proportions
                                      degradation_shape = 2.5,  # Gamma distribution from which degradation for each contributor is sampled
                                      degradation_scale = 1e-3, # Gamma distribution from which degradation for each contributor is sampled
                                      number_of_simulated_profiles = 1,
                                      lower_NoC = 2,
                                      upper_NoC = 6) 
  
  #convert it to a scan input format
  plot_EPG <- FALSE
  EPG_plot_savename <- NULL
  #the numbered input
  if (!settings$include_noise_generation){
    numbered_input <- matrix(rep(1:500, 6), ncol=500, byrow = TRUE)/500
    numbered_input_one_sample <- array(numbered_input, dim=c(1, 6, 500, 1))
  }
  
  for(iSim in 1:number_of_profiles_to_simulate){
    #simulate basic profile information in STRmix input format
    simulated_profile <- simulate_profile(profile_simulation_settings)
    #create the mock EPG
    mock_EPG <- get_mock_EPG(simulated_profile, plot_EPG, EPG_plot_savename)
    #generate the profile
    mock_EPG_array <- array(mock_EPG, dim=c(1,settings$number_of_dyes,settings$number_of_scanpoints,1))/settings$saturation
    #the numbered input
    if (settings$include_noise_generation){
      numbered_input_one_sample <- array(runif(settings$number_of_dyes*500), dim=c(1, settings$number_of_dyes, 500, 1))
    }
    generated_profile_from_sim <- generator %>% predict(list(mock_EPG_array, numbered_input_one_sample))
    #format the input for plotting
    simulated_profile_complete_input <- t(generated_profile_from_sim[1,,,1])*settings$saturation
    simulated_profile_complete_input <- cbind(seq(1:dim(simulated_profile_complete_input)[1]), simulated_profile_complete_input)
    colnames(simulated_profile_complete_input) <- c("scan", paste("dye", seq(1:(dim(simulated_profile_complete_input)[2] - 1)), sep=""))
    simulated_profile_complete_input <- as.data.frame(simulated_profile_complete_input)
    #plot the profile
    simulated_profile_complete <- get_mock_EPG_from_scans(simulated_profile_complete_input, 1, settings$number_of_scanpoints)
    #save the profile
    savename <- paste(here(), "\\simulated_profiles\\", iSim, "_simulated_profile.jpg", sep="")
  #  ggsave(savename, simulated_profile_complete, units='px', width=3000, height=5000)
    ggsave(savename, simulated_profile_complete, units='px', width=1500, height=2000)
  }
}

#------------------------------------------ Figure 10--------------------------------------------------

plot_Fig10 <- FALSE

if (plot_Fig10){
  
  generator_increment_filename <- paste(here(), "\\saved_models\\Profile_generator_GAN_trained_model_random_24-09-2023.h5", sep="")
  generator <-load_model_hdf5(generator_increment_filename)
  
  #the number of profiles to simulate
  number_of_profiles_to_simulate <- 20
  
  #load in scripts
  source('simulate_profiles.R')
  source('get_mock_EPG.R')
  source('get_mock_EPG_from_scans.R')
  #profile simulation settings
  profile_simulation_settings <- list(min_template = 20,   # note that these can be arrays if you want to create specific mixture proportions
                                      max_template = 100,       # note that these can be arrays if you want to create specific mixture proportions
                                      degradation_shape = 2.5,  # Gamma distribution from which degradation for each contributor is sampled
                                      degradation_scale = 1e-3, # Gamma distribution from which degradation for each contributor is sampled
                                      number_of_simulated_profiles = 1,
                                      lower_NoC = 2,
                                      upper_NoC = 6) 
  
  #convert it to a scan input format
  plot_EPG <- FALSE
  EPG_plot_savename <- NULL
  #the numbered input
  if (!settings$include_noise_generation){
    numbered_input <- matrix(rep(1:500, 6), ncol=500, byrow = TRUE)/500
    numbered_input_one_sample <- array(numbered_input, dim=c(1, 6, 500, 1))
  }
  
  
  #simulate basic profile information in STRmix input format
  simulated_profile <- simulate_profile(profile_simulation_settings)
  #create the mock EPG
  mock_EPG <- get_mock_EPG(simulated_profile, plot_EPG, EPG_plot_savename)
  #generate the profile
  mock_EPG_array <- array(mock_EPG, dim=c(1,settings$number_of_dyes,settings$number_of_scanpoints,1))/settings$saturation
    
  for(iSim in 1:number_of_profiles_to_simulate){
    #the numbered input
    if (settings$include_noise_generation){
      numbered_input_one_sample <- array(runif(settings$number_of_dyes*500), dim=c(1, settings$number_of_dyes, 500, 1))
    }
    generated_profile_from_sim <- generator %>% predict(list(mock_EPG_array, numbered_input_one_sample))
    #format the input for plotting
    simulated_profile_complete_input <- t(generated_profile_from_sim[1,,,1])*settings$saturation
    simulated_profile_complete_input <- cbind(seq(1:dim(simulated_profile_complete_input)[1]), simulated_profile_complete_input)
    colnames(simulated_profile_complete_input) <- c("scan", paste("dye", seq(1:(dim(simulated_profile_complete_input)[2] - 1)), sep=""))
    simulated_profile_complete_input <- as.data.frame(simulated_profile_complete_input)
    #plot the profile
    simulated_profile_complete <- get_mock_EPG_from_scans(simulated_profile_complete_input, 1, settings$number_of_scanpoints)
    #save the profile
    savename <- paste(here(), "\\simulated_profiles\\same_", iSim, "_simulated_profile.jpg", sep="")
    #  ggsave(savename, simulated_profile_complete, units='px', width=3000, height=5000)
    ggsave(savename, simulated_profile_complete, units='px', width=1500, height=2000)
  }
}

#----------------------------------------------------------------------------------------------------


