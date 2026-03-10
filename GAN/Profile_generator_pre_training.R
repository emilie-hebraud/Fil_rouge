profile_generator_pre_training <- function(settings,
                                           complete_model,
                                           number_epochs,
                                           save_profile_images,
                                           save_profile_csvs,
                                           train_with_initial_sample_weights,
                                           training_data){
  
  #libraries
  library('tensorflow')
  library('keras')
  library('tidyverse')

  #-------------------------  some global ANN settings ---------------------------
  #assigning local variables with the values from the settings list
  number_of_dyes <- settings$number_of_dyes
  startScan <- settings$startScan
  number_scanpoint_in_input_layer <- settings$number_of_scanpoints
  saturation <- settings$saturation
  profile_dirs <- settings$real_and_smooth_dirs
  #the end scan point
  endScan <- startScan + number_scanpoint_in_input_layer - 1 #need to minus 1 so that startScan:endScan gives number_scanpoint_in_input_layer data points
  
  
  number_of_profiles <- dim(training_data$input_profiles)[1]
  
  
  
  number_of_profiles_to_use <- dim(training_data$input_profiles)[1]
  batch_size <- as.integer(number_of_profiles_to_use/100, 0)
  number_batches <- round(number_of_profiles_to_use/batch_size - 0.5, 0)
  
  ######################################################################################
  
  # formats the training inputs
  input_profiles_smooth <- array(training_data$input_profiles[1:number_of_profiles_to_use, 1:number_of_dyes, ], 
                                 dim=c(number_of_profiles_to_use, number_of_dyes, number_scanpoint_in_input_layer, 1))
  
  targets_profiles_real <- array(training_data$target_profiles[1:number_of_profiles_to_use, 1:number_of_dyes, ], 
                                 dim=c(number_of_profiles_to_use, number_of_dyes, number_scanpoint_in_input_layer, 1))
  
  molw_profile <- training_data$input_profiles[1:number_of_profiles_to_use, 4, ]
  channel_names <- c("channel_1","channel_3","channel_4")
  # convert arrays into tensors
  input_profiles_smooth <- tf$convert_to_tensor(input_profiles_smooth, dtype = "float32")
  targets_profiles_real <- tf$convert_to_tensor(targets_profiles_real, dtype = "float32")
  
  # the numbered input array to assist with sequential training
  numbered_input <- matrix(rep(1:500, number_of_dyes), ncol=500, byrow = TRUE)/500
  numbered_input_array <- array(NA, dim=c(number_of_profiles_to_use, number_of_dyes, 500, 1))
  
  for(iBatch in 1:number_of_profiles_to_use){
    numbered_input_array[iBatch, , , 1] <- numbered_input
  }
  
  numbered_input_array <- tf$convert_to_tensor(
    numbered_input_array[1:number_of_profiles_to_use,,,],
    dtype = "float32")
  
  range(input_profiles_smooth)
  range(targets_profiles_real)
  
  loss_values <- NULL
  val_loss_values <- NULL
  
  if(train_with_initial_sample_weights){
    #fits the model but with the fixed random values generated upfront
    
    history <- complete_model$fit(
      list(input_profiles_smooth, numbered_input_array),
      targets_profiles_real,
      epochs = number_epochs,
      batch_size = batch_size,
      verbose = 1,
      validation_split = 0.1,
      sample_weight = input_profiles_smooth*100 + 1 #this just forces the classifier to concentrate on getting the peaks right
    )
    
    
    loss_values <- c(loss_values, history$history$loss)
    val_loss_values <- c(val_loss_values, history$history$val_loss)
  }
  #now trains without the emphasis on peaks
  history <- complete_model$fit(
    list(input_profiles_smooth, numbered_input_array),
    targets_profiles_real,
    epochs = number_epochs,
    batch_size = batch_size,
    verbose = 1,
    validation_split = 0.1
  )
  
  
  loss_values <- c(loss_values, history$history$loss)
  val_loss_values <- c(val_loss_values, history$history$val_loss)
  
  #plots the loss
  jpeg(file=paste(here(), "/training_results/ANN_generator_pre-training_loss.jpg", sep=""), width = 2000, height=2000, res=200)
    yplot_max <- as.integer(max(loss_values, val_loss_values))
    yplot_min <- as.integer(min(loss_values, val_loss_values))
    plot(loss_values, type='l', col="blue", ylim=c(yplot_min,yplot_max), ylab="MeanSquaredError", xlab="epoch", log='y')
    lines(val_loss_values, col="dark green")
    abline(v = number_epochs, col = "grey", lty = 2)
    legend("topright", col=c("blue", "dark green"), lwd=1, legend=c("loss", "val loss"))
  dev.off()
  #########################################################################################
  
  #-------------------------------------------------------------------------------------
  #shows the performance of the ANN in making the smoothed data from the training set look realistic
  if(save_profile_images | save_profile_csvs){
    
    #get predictions from training dataset
#    all_preds <- complete_model %>% predict(input_profiles_smooth)
    all_preds <- complete_model$predict(list(input_profiles_smooth, numbered_input_array))
    
    for (iProfile in 1:dim(all_preds)[1]){
      #the profile you wish to view
      pred_plot_profile <- iProfile
      
      #set saturation (leave as 1 to have raw ANN output)
      saturation <- settings$saturation
      #gets the data for that profile
      ANN_predicted_profile_data_plot <- all_preds[pred_plot_profile, , , 1]*saturation
      input_profiles_smooth_plot <- input_profiles_smooth[pred_plot_profile, , , 1]*saturation
      targets_profiles_real_plot <- targets_profiles_real[pred_plot_profile, , , 1]*saturation
      
      # get the molw for this specific profile
      molw <- molw_profile[pred_plot_profile,]
      
      
      
      library(tensorflow)
      # Conversion explicite en float32 (sécurité)
      input_tensor  <- tf$cast(input_profiles_smooth, tf$float32)
      target_tensor <- tf$cast(targets_profiles_real, tf$float32)
      
      #get the number of dyes
      number_of_dyes <- dim(ANN_predicted_profile_data_plot)[1]
      
      
      #-------converts ANN output and training data input to dataframes--------
      #transpose data
      ANN_predicted_profile_data_plot.df <- t(ANN_predicted_profile_data_plot)
      #add scans
      ANN_predicted_profile_data_plot.df <- cbind(1:5000, ANN_predicted_profile_data_plot.df)
      # add molw as a column
      ANN_predicted_profile_data_plot.df <- cbind(ANN_predicted_profile_data_plot.df, molw)
      #add dye names
      colnames(ANN_predicted_profile_data_plot.df) <- c("index", channel_names, "molw")
      #convert to dataframe
      ANN_predicted_profile_data_plot.df <- as.data.frame(ANN_predicted_profile_data_plot.df)
      
      #transpose data
      #re conversion tensor -> array
      input_profiles_smooth_plot <- as.array(input_profiles_smooth_plot[1:number_of_dyes, ])
      input_profiles_smooth_plot.df <- t(input_profiles_smooth_plot)
      #add scans
      input_profiles_smooth_plot.df <- cbind(1:5000, input_profiles_smooth_plot.df)
      #add molw
      input_profiles_smooth_plot.df <- cbind(input_profiles_smooth_plot.df, molw)
      #add dye names
      colnames(input_profiles_smooth_plot.df) <- c("index", channel_names, "molw")
      #convert to dataframe
      input_profiles_smooth_plot.df <- as.data.frame(input_profiles_smooth_plot.df)
      # print(head(input_profiles_smooth_plot.df))
      
      profile_name <- basename(training_data$input_paths[iProfile])
      profile_name_without_suffix <- substr(x = profile_name, start = 1, stop = nchar(profile_name) - 4)
      #calls the comparison plotting method
      source('compare_raw_v_smoothed_from_scans.R')
      save_name <- paste(here(), "/pre_process_data/", profile_name_without_suffix, iProfile, sep="")
      #mettre légendes sur graph généré, reel; smooth avec des couleurs différentes
      # print(save_name)
      
      compare_raw_v_smoothed_from_scans(ANN_predicted_profile_data_plot.df,
                                        input_profiles_smooth_plot.df,
                                        save_name,
                                        1,
                                        5000,
                                        1,
                                        5000,
                                        save_profile_images,
                                        generator=TRUE)
      
      if (save_profile_csvs){
        #saves predicted profile
        write.table(x = ANN_predicted_profile_data_plot.df,
                    file = paste(here(), "/prediction_generator/", profile_name_without_suffix, "_fake.csv", sep=""),
                    sep=";",
                    col.names = TRUE,
                    row.names = FALSE)
      }
    }
   
  }
  #---------------------------------------------------------------------------------------------------------------------------
   
  #returns the model
  return(complete_model)
}
