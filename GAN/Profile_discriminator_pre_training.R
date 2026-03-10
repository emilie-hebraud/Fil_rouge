profile_discriminator_pre_training <- function(settings,
                                               discriminator,
                                               number_epochs,
                                               training_data){

  #libraries
  library('tensorflow')
  library('keras')
  library('tidyverse')
  
  discriminator_version <- 5
  
  #-------------------------  some global ANN settings ---------------------------
  
  #assigning local variables with the values from the settings list
  number_of_dyes <- settings$number_of_dyes
  startScan <- settings$startScan
  number_scanpoint_in_input_layer <- settings$number_of_scanpoints
  
  #--------------- run the model -----------------------------------------------------
  
  number_of_profiles_to_use <- dim(training_data$input_profiles)[1]
  batch_size <- as.integer(number_of_profiles_to_use/100)
  number_batches <- round(number_of_profiles_to_use/batch_size - 0.5, 0) #-0.5 ensures the round down of number of batches
  
  src_images_smooth <- array(training_data$source_profile_smooth[1:number_of_profiles_to_use,,], dim=c(number_of_profiles_to_use,number_of_dyes,number_scanpoint_in_input_layer,1))
  targets_profiles_real_and_fake <- array(training_data$input_profiles[1:number_of_profiles_to_use,,], dim=c(number_of_profiles_to_use, number_of_dyes, number_scanpoint_in_input_layer,1))
  targets_single_value <- training_data$targets[1:number_of_profiles_to_use]
  
  src_images_smooth <- tf$convert_to_tensor(src_images_smooth, dtype = "float32")
  targets_profiles_real_and_fake <- tf$convert_to_tensor(targets_profiles_real_and_fake, dtype = "float32")
  # convert to tensors
  
  #need to replicate the targets so that they have 20 outputs (the length of the final output layer of the patchGAN)
  number_patchGAN_outputs <- 20
  targets <- array(NA, dim=c(number_of_profiles_to_use, 1, number_patchGAN_outputs, 1))
  for (iOut in 1:number_of_profiles_to_use){
    targets[iOut,1,,1] <- targets_single_value[iOut]
  }
  targets <- tf$convert_to_tensor(targets, dtype="float32")
  
  history <- discriminator$fit(
    list(src_images_smooth, targets_profiles_real_and_fake),
    targets,
    epochs = number_epochs,
    batch_size = batch_size,
    verbose = 1,
    validation_split = 0.1
  )
  
  loss_values <- history$history$loss
  val_loss_values <- history$history$val_loss
  
  jpeg(file=paste(here(), "/training_results/Profile_discriminator_loss_V", 
                  discriminator_version, ".jpg", sep=""), 
       width = 2000, height = 2000, res = 200)
  
  # Replace zeros or negatives (just in case) with a small positive value
  loss_plot <- loss_values
  val_loss_plot <- val_loss_values
  loss_plot[loss_plot <= 0] <- 1e-6
  val_loss_plot[val_loss_plot <= 0] <- 1e-6
  
  # Use actual min/max (do NOT convert to integer)
  yplot_min <- min(loss_plot, val_loss_plot)
  yplot_max <- max(loss_plot, val_loss_plot)
  
  # Plot with log scale
  plot(loss_plot, type='l', col="blue", ylim=c(yplot_min, yplot_max),
       ylab="MeanSquaredError", xlab="Epoch", log='y')
  lines(val_loss_plot, col="dark green")
  legend("topright", col=c("blue", "darkgreen"), lwd=1, legend=c("loss", "val loss"))
  
  dev.off()
  
  all_preds <- discriminator$predict(list(training_data$source_profile_smooth, training_data$input_profiles))
  #average the preds per profile for plotting
  preds_to_plot <- rep(NA, number_of_profiles_to_use)
  for (iOut in 1:number_of_profiles_to_use){
    preds_to_plot[iOut] <- mean(all_preds[iOut,1,,1])
  }
  
  jpeg(file=paste(here(), "/training_results/Profile_discriminator_performance_V", discriminator_version, "_rand_full.jpg", sep=""), width = 2000, height=2000, res=200)
    boxplot(preds_to_plot[training_data$targets == 0], at = 1, xlim=c(0.5, 2.5), ylim=c(0,1), xlab="category", ylab="probability", xaxt='n')
    boxplot(preds_to_plot[training_data$targets == 1], at = 2, add=TRUE)
    mtext(side = 1, line = 1, at = 1:2, text = c("fake", "real"))
  dev.off()

  #returns the model
  return(discriminator)
}