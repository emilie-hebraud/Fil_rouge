create_discriminator_pretrain_training_data_from_models <- function(settings,
                                                                    profile_discriminator,
                                                                    profile_generator,
                                                                    include_noise,
                                                                    GAN_training_dataset){
  
  #libraries
  library('tensorflow')
  library('keras')
  library('tidyverse')
  library('tfdatasets')
  
  #some settings
  number_of_profile_gan <- dim(GAN_training_dataset$real_profile_target)[1]
  number_of_dyes <- settings$number_of_dyes
  number_of_scanpoints <- settings$number_of_scanpoints
  include_noise <- settings$include_noise_generation
  
  #creates the numbered input for the generator
  if (!include_noise){
    numbered_input <- matrix(rep(1:500, 6), ncol=500, byrow = TRUE)/500
    numbered_input_array <- array(numbered_input, dim=c(1, 6, 500, 1))
  }
  # -------------------------------- set up the discriminator only training ----------------------------------
  
  src_images_smooth <- array(NA, dim = c(number_of_profile_gan*2, number_of_dyes, number_of_scanpoints, 1))
  targets_profiles_real_and_fake <- array(NA, dim = c(number_of_profile_gan*2, number_of_dyes, number_of_scanpoints, 1))
  number_patchGAN_outputs <- 20
  discriminator_targets <- array(NA, dim=c(number_of_profile_gan*2, 1, number_patchGAN_outputs, 1))
  ANN_array_input_dim <- c(1,number_of_dyes,number_of_scanpoints,1)
  for (iProfile in 1:number_of_profile_gan){
    odd_spot <- 2*(iProfile - 1) + 1 #real profiles
    even_spot <- odd_spot + 1 #generated profiles
    #smooth profiles
    smooth_profile_formatted <- array(GAN_training_dataset$source_profile_smooth[iProfile,,], dim=ANN_array_input_dim)
    src_images_smooth[odd_spot,,,] <- smooth_profile_formatted
    src_images_smooth[even_spot,,,] <- smooth_profile_formatted
    #real profiles
    targets_profiles_real_and_fake[odd_spot,,,] <- array(GAN_training_dataset$real_profile_target[iProfile,,], dim=ANN_array_input_dim)
    #fake profile
    if (include_noise) {
      numbered_input_array <- array(runif(number_of_dyes*500), dim=c(1, number_of_dyes, 500, 1))
    }
    generated_profile <- profile_generator$predict(list(smooth_profile_formatted, numbered_input_array), verbose = FALSE)

    targets_profiles_real_and_fake[even_spot,,,] <- generated_profile
    #targets
    discriminator_targets[odd_spot,1,,1] <- 1
    discriminator_targets[even_spot,1,,1] <- 0
    #progress
    if (iProfile%%100==0){
      cat(round(100*iProfile/number_of_profile_gan, 0), "% | ")
    }
  }
  
  # -------------------------------------------------------------------------------------------------------------
  
  return(list(src_images_smooth = src_images_smooth, targets_profiles_real_and_fake = targets_profiles_real_and_fake, discriminator_targets = discriminator_targets))
  
}