profile_GAN_trainer <- function(settings,
                                profile_discriminator,
                                profile_generator,
                                use_sample_weights,
                                epochs_of_training,
                                GAN_training_dataset){
  
  #libraries
  library('here')
  library('tensorflow')
  library('keras')
  library('tidyverse')
  library('tfdatasets')

  #pass in the Tensorflow objects
  discriminator <- profile_discriminator
  generator <- profile_generator
  
  #creates the discriminator only training dataset (as this is used later to train just the disciminator at points through the GAN training)
  source('create_discriminator_pretrain_training_data_from_models.R')
  include_noise <- FALSE
  discriminator_only_training_dataset <- create_discriminator_pretrain_training_data_from_models(settings,
                                                                                                 discriminator,
                                                                                                 generator,
                                                                                                 include_noise,
                                                                                                 GAN_training_dataset)
  # discrimintor_only_training_dataset_dims <- dim(discrimintor_only_training_dataset$targets_profiles_real_and_fake)
  #some settings
  number_of_profile_gan <- dim(GAN_training_dataset$real_profile_target)[1]
  
  number_of_dyes <- settings$number_of_dyes
  number_of_scanpoints <- settings$number_of_scanpoints
  noise_multiplier <- settings$saturation_max/settings$saturation
  include_noise <- settings$include_noise_generation
  
  #profiles per batch
  profiles_per_batch <- 8
  #creates the numbered input for the generator
  if(!include_noise){
    numbered_input <- matrix(rep(1:500, number_of_dyes), ncol=500, byrow = TRUE)/500
    numbered_input_one_sample <- array(numbered_input, dim=c(1, number_of_dyes, 500, 1))
    numbered_input_array <- array(NA, dim=c(profiles_per_batch, number_of_dyes, 500, 1))
    
    for(iBatch in 1:profiles_per_batch){
      numbered_input_array[iBatch, , , 1] <- numbered_input
    }
  }
  # -------------------------- functions ------------------------------------------------
  
  #a function to plot the output
  plot_incremental_GAN_profiles <- function(epoch, starting_epoch = 50) {
    starting_epoch <- starting_epoch
    plot_profile <- grep(pattern="M1_pl10036", x = GAN_training_dataset$input_paths) #this is just the profile I am using
    plot_dye <- 1
    
    real_profile_plot <- GAN_training_dataset$real_profile_target[plot_profile,1:number_of_dyes+1,]
    real_profile_array <- array(real_profile_plot, dim=c(1,number_of_dyes,number_of_scanpoints,1))
    smooth_profile_plot <- GAN_training_dataset$source_profile_smooth[plot_profile,1:number_of_dyes+1,]
    molw <- GAN_training_dataset$source_profile_smooth[plot_profile, 4, ]
    
    if(include_noise){
      numbered_input_one_sample <- array(runif(number_of_dyes*500), dim=c(1, number_of_dyes, 500, 1))
    }
    
    #convert to array
    smooth_profile_array <- array(smooth_profile_plot, dim=c(1,number_of_dyes,number_of_scanpoints,1))
    generated_profile_for_plot <- generator$predict(list(smooth_profile_array, numbered_input_one_sample))
    
    generated_profile_plot <- generated_profile_for_plot[1, , , 1]
    disc_value_on_generated_profile <- discriminator$predict(list(smooth_profile_array, generated_profile_for_plot))
    #calculates the discriminator value for the fake profile
    savename <- paste(here(), "/GAN_increments/training_profile_", plot_profile, "_epoch_", (starting_epoch+epoch), "_NO_noise.jpg", sep="")
    # print(savename)
    
    #calculates the discriminator value for the real profile
    disc_value_on_real_profile <- discriminator$predict(list(smooth_profile_array, real_profile_array))
      
    print(dim(real_profile_plot))
    slice_real <- real_profile_plot[plot_dye, ]  # dim(slice_real) = 4 x 5000
    slice_smooth <- smooth_profile_plot[plot_dye , ]
    slice_generated <- generated_profile_plot[plot_dye, ]
     
     # Calculer les min/max pour le plot
     plot_ymax <- max(slice_real, slice_smooth, slice_generated)
     plot_ymin <- min(slice_real, slice_smooth, slice_generated)
     
    #builds the plot
    jpeg(file=savename, width = 4000, height=4000, res=400)
      par(mfrow=c(3,1))
      #the smooth profile
      plot(x = molw,
           y = smooth_profile_plot[plot_dye,],
           ylim=c(-1000, 10000),
           xlab="molw (pdb)",
           ylab="rfu",
           main=paste("input profile (smoothed) : epoch ", (starting_epoch+epoch), sep=""),
           type='l',
           col="black")
      #the generated profile
      plot(x = molw,
           y = generated_profile_plot[plot_dye,],
           ylim=c(-1000, 10000),
           xlab="molw (pdb)",
           ylab="rfu",
           main=paste("fake profile (from ANN) : disciminator ", mean(disc_value_on_generated_profile), sep=""),
           type='l',
           col="black")
      #the real profile
      plot(x = molw,
           y = real_profile_plot[plot_dye,],
           ylim=c(-1000, 10000),
           xlab="molw (pdb)",
           ylab="rfu",
           main=paste("real profile (from lab) : disciminator ", mean(disc_value_on_real_profile), sep=""),
           type='l',
           col="black")
      par(mfrow=c(1,1))
    dev.off()
  }
  
  #-------------------------------------------------------------------------------------
  
  #code adapted from: https://blogs.rstudio.com/ai/posts/2018-09-20-eager-pix2pix/
  
  # discriminator_optimizer <- tf$compat$v1$train$AdamOptimizer(5e-5, beta1 = 0.5, beta2 = 0.999)
  # generator_optimizer <- tf$compat$v1$train$AdamOptimizer(5e-5, beta1 = 0.5, beta2 = 0.999)
  
  generator_optimizer <- tf$keras$optimizers$Adam(1e-4, beta_1 = 0.5)
  discriminator_optimizer <- tf$keras$optimizers$Adam(1e-4, beta_1 = 0.5)
  
  generator_loss <- function(disc_judgment, generated_output, target, input_image) {
    gen_lambda <- 1 # 100 was the value chosen by the authors of the pix2pix paper (but I am not scaling the same way)
    #discriminator component
    labels <- tf$ones_like(disc_judgment)
    #add noise to labels
    labels <- labels + 0.05 * tf$random$uniform(tf$shape(labels))
    #calculate loss
    gan_loss <- tf$compat$v1$losses$sigmoid_cross_entropy(labels, disc_judgment)
    #generator component
    if (use_sample_weights){
      sample_weights <- (input_image*10 + 0.9) #weight the peaks as more important than the baseline
    } else {
      sample_weights <- 1
    }
    l1_loss <- tf$reduce_mean(tf$abs(target - generated_output)*sample_weights)
    #total loss
    total_loss <- gan_loss + (gen_lambda * l1_loss)
    return(total_loss)
  }
  
  discriminator_loss <- function(real_output, generated_output) {
    #labels
    real_labels <- tf$ones_like(real_output)
    #add noise to labels
    real_labels <- real_labels + 0.05 * tf$random$uniform(tf$shape(real_labels))
    #calculate loss
    real_loss <- tf$compat$v1$losses$sigmoid_cross_entropy(multi_class_labels = real_labels, logits = real_output)
    
    #labels
    gen_labels <- tf$zeros_like(generated_output)
    #add noise to labels
    gen_labels <- gen_labels + 0.05 * tf$random$uniform(tf$shape(gen_labels))
    #calculate loss
    generated_loss <- tf$compat$v1$losses$sigmoid_cross_entropy(multi_class_labels = gen_labels, logits = generated_output)
    real_loss + generated_loss
  }
  
  
  train <- function(GAN_training_dataset, num_epochs) {
    #tallies to hold the training info
    disc_loss_tally <- NULL
    gen_loss_tally <- NULL
    #number of profiles
    number_of_profile_gan <- dim(GAN_training_dataset$real_profile_target)[1]
    if (profiles_per_batch > number_of_profile_gan){
        profiles_per_batch <- number_of_profile_gan
      }
    number_batches <- number_of_profile_gan/profiles_per_batch 
    #cycle through epochs
    for (epoch in 1:num_epochs) {
      total_loss_gen <- 0
      total_loss_disc <- 0
      
      if (number_batches %% 1 > 0){
        number_batches <- round(number_batches + 0.5, 0)  #the addition of 0.5 to make sure all samples included
      }

      #cycle through batches
      for (iBatch in 1:number_batches){
        #starting indicator
        cat("Epoch ", epoch, "/", num_epochs, 
            ": Batch ", iBatch, "/", number_batches)
        #set up profiles in batch
        batch_start <- (iBatch - 1)*profiles_per_batch + 1
        batch_end <- iBatch*profiles_per_batch
        if (batch_end <= number_of_profile_gan){
          #normal batch
          profiles_to_include <- c(seq(from = batch_start, to = batch_end, by = 1))
        } else {
          #batches is too small (the last one) is just filled with the required number of profiles from the beginning again
          numeric_batch_size <- number_of_profile_gan - batch_start + 1
          profiles_to_include <- c(seq(from = 1, to = profiles_per_batch - numeric_batch_size, by = 1), seq(from = batch_start, to = number_of_profile_gan, by = 1))
        }
        #get inputs and targets for batch
        input_image <- array(GAN_training_dataset$source_profile_smooth[profiles_to_include,,], dim=c(profiles_per_batch,number_of_dyes,number_of_scanpoints,1))
        # print("input image class")
        # print(class(input_image))
        
        #get target image
        target_image <- array(GAN_training_dataset$real_profile_target[profiles_to_include,,], dim=c(profiles_per_batch,number_of_dyes,number_of_scanpoints,1))
        # print("target image class")
        # print(class(target_image))
        
        if(include_noise){
          numbered_input_array <- array(runif(number_of_dyes*500*profiles_per_batch), dim=c(profiles_per_batch, number_of_dyes, 500, 1))
        }
        # convert array to tensor
        input_image_tensor <- tf$convert_to_tensor(input_image, dtype = "float32")
        target_image_tensor <- tf$convert_to_tensor(target_image, dtype = "float32")
        
        numbered_tensor <- tf$convert_to_tensor(
          numbered_input_array[1:profiles_per_batch,,,],
          dtype = "float32"
        )
        # print("input image/ target image tensor class")
        # print(class(input_image_tensor))
        # print(class(numbered_tensor))
        
        with(tf$GradientTape() %as% gen_tape, {
          with(tf$GradientTape() %as% disc_tape, {
            #generates a fake profile
            gen_output <- generator(list(input_image_tensor,numbered_tensor))
            
            # print("gen output class")
            # print(class(gen_output))
            
            #gets the discriminator output for the real image
            disc_real_output <- discriminator((list(input_image_tensor, target_image_tensor)))

            # print("disc real output class")
            # print(class(disc_real_output))
            
            #gets the discriminator output for the generated image
            disc_generated_output <- discriminator((list(input_image_tensor, gen_output)))
            
            # print("disc gen output class")
            # print(class(disc_generated_output))
            
            #calculate the losses
            gen_loss <- generator_loss(disc_generated_output, gen_output, target_image_tensor, input_image_tensor)
            disc_loss <- discriminator_loss(disc_real_output, disc_generated_output)
          })
        })
        #calculate the gradients
        generator_gradients <- gen_tape$gradient(gen_loss, generator$trainable_weights)
        discriminator_gradients <- disc_tape$gradient(disc_loss, discriminator$trainable_weights)
        
        #apply the gradients in order to train the ANNs
        # generator_optimizer$apply_gradients(transpose(list(generator_gradients, generator$variables)))
        # discriminator_optimizer$apply_gradients(transpose(list(discriminator_gradients, discriminator$variables)))
        
        # cat("Gradients générateur:", length(generator_gradients), "\n")
        # cat("Variables générateur:", length(generator$trainable_variables), "\n")
        # cat("Gradients discriminateur:", length(discriminator_gradients), "\n")
        # cat("Variables discriminateur:", length(discriminator$trainable_variables), "\n")
        
        generator_optimizer$apply_gradients(
          Map(list, generator_gradients, generator$trainable_weights)
        )
        discriminator_optimizer$apply_gradients(
          Map(list, discriminator_gradients, discriminator$trainable_weights)
        )
        # print("passed optimization")
        #update tallies
        total_loss_gen <- total_loss_gen + as.numeric(gen_loss)
        total_loss_disc <- total_loss_disc + as.numeric(disc_loss)
        cat("| gen loss = (", as.numeric(gen_loss), 
            ") : disc loss = (", as.numeric(disc_loss), 
            ")\n")
     
        #every 10 iterations store the predictions so more training can be done
        # if(epoch%%10 == 0){
        #   #places the generated DNA profiles into the even positions of the discriminator dataset in their correct position
        #   profiles_to_replace <- seq(from = batch_start, to = min(batch_end, number_of_profile_gan), by = 1)*2
        #   gen_output_start <- 1
        #   gen_output_end <- min(batch_end, number_of_profile_gan) - batch_start + 1
        #   discrimintor_only_training_dataset$src_images_smooth[profiles_to_replace,,,] <- array(gen_output[gen_output_start:gen_output_end,,,], dim=dim(gen_output[gen_output_start:gen_output_end,,,]))
        # }
        
      } #next batch
      
      #add to tally (scaled by the number of batches)
      disc_loss_tally <- c(disc_loss_tally, total_loss_disc/number_batches)
      gen_loss_tally <- c(gen_loss_tally, total_loss_gen/number_batches)
      
      #every 10 iterations
      if(epoch%%10 == 0){
        #save the current models
        date_string <- format(today(), "%d-%m-%Y")
        save_model_hdf5(generator, paste(here(), "/saved_models/Profile_generator_GAN_trained_model_epoch_", epoch, "_", date_string, ".h5", sep=""))
        save_model_hdf5(discriminator, paste(here(), "/saved_models/Profile_discriminator_GAN_trained_model_epoch_", epoch, "_", date_string, ".h5", sep=""))
        #do some additional discriminator training
        # discriminator %>% fit(
        #   list(discrimintor_only_training_dataset$src_images_smooth, discrimintor_only_training_dataset$targets_profiles_real_and_fake),
        #   discrimintor_only_training_dataset$discriminator_targets,
        #   epochs = 100,
        #   batch_size = number_of_profile_gan,
        #   view_metrics = TRUE,
        #   validation_split = 0.1
        # )
        #save the loss values
        write.table(x = disc_loss_tally, file = paste(here(), "/GAN_increments/Disc_loss_epoch_", epoch, ".csv", sep=""), sep=";")
        write.table(x = gen_loss_tally, file = paste(here(), "/GAN_increments/Gen_loss_epoch_", epoch, ".csv", sep=""), sep=";")
      }

      #print progress to screen
      cat("Epoch ", epoch, "\n")
      cat("Generator loss: ", total_loss_gen/number_batches, "\n")
      cat("Discriminator loss: ", total_loss_disc/number_batches, "\n\n")
      #save the output after each epoch graphically
      plot_incremental_GAN_profiles(epoch)
    } #next epoch
    return(list(disc_loss = disc_loss_tally, gen_loss = gen_loss_tally))
  }
  
  
  # ------------------------------------ carry out some training ----------------------------------------
  #the starting model performance
  plot_incremental_GAN_profiles(0)
  #trains the model for some epochs
  training <- train(GAN_training_dataset, epochs_of_training)
  
  # ----------------------------------- plot the results ----------------------------------------------------

  #plots and saves the losses
  disc_loss_tally <- training$disc_loss
  gen_loss_tally <- training$gen_loss
  
  date_string <- format(today(), "%d-%m-%Y")
  savename <- paste(here(), "/training_results/GAN_losses_", date_string, ".jpg", sep="")
  jpeg(file=savename, width = 4000, height=4000, res=400)
    par(mfrow=c(2,1))
    yplot_max <- max(disc_loss_tally)
    yplot_min <- min(disc_loss_tally)
    
    plot(x = 1:length(disc_loss_tally),
         y = disc_loss_tally,
         ylim = c(0,5000),
         ylab = "loss",
         xlab = "epoch",
         type = 'l',
         main = "discriminator total loss",
         col="black",
         lwd = 2)
    
    yplot_max <- max(gen_loss_tally)
    yplot_min <- min(gen_loss_tally)
    
    plot(x = 1:length(gen_loss_tally),
         y = gen_loss_tally,
         ylim = c(0, 5000),
         ylab = "loss",
         xlab = "epoch",
         type = 'l',
         main = "generator total loss",
         col="black",
         lwd = 2)
    par(mfrow=c(1,1))
  dev.off()
  
  
  # -----------------------------------------------------------------------------------------------------------------------------------------------

  #return the updated models
  return(list(updated_profile_generator = generator, updated_profile_discriminator = discriminator, profile_dloss = disc_loss_tally, profile_gloss = gen_loss_tally))
  
}