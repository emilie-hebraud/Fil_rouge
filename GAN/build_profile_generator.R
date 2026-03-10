# based on https://blogs.rstudio.com/ai/posts/2018-09-20-eager-pix2pix/

build_profile_generator <- function(settings){
  #libraries
  library('tensorflow')
  library('keras')
  library('tidyverse')
  
  # ------------------ setting up the profile ANN --------------------------
  number_of_dyes <- settings$number_of_dyes
  number_scanpoint_in_input_layer <- settings$number_of_scanpoints
  print(number_scanpoint_in_input_layer)
  number_channels <- 1
  layer1_number_of_filters <- 256
  layer2_number_of_filters <- 128
  
  layer1_input_shape <- c(number_of_dyes, number_scanpoint_in_input_layer, number_channels)
  
  #build the ANN
  profile_inputs <- layer_input(shape=layer1_input_shape, name="profile_input_layer")
  
  profile_downsize <- profile_inputs %>% 
    layer_conv_2d(filters = 16, strides=c(1, 10), kernel_size = c(1, 10), activation = "relu", padding = 'same', name = "downsize1")
  
  profile_horizontal_model <- profile_downsize %>% 
    layer_conv_2d(filters = 32, strides=c(1, 1), kernel_size = c(1, 5), activation = "relu", padding = 'same', name = "profile_horizontal_model_L1") %>%
    layer_batch_normalization()%>%
    layer_conv_2d(filters = 32, strides=c(1, 1), kernel_size = c(1, 5), activation = "relu", padding = 'same', name = "profile_horizontal_model_L2") %>%
    layer_batch_normalization()
    
  profile_horizontal_model_downsize2 <- profile_horizontal_model %>% 
    layer_conv_2d(filters = 64, strides=c(1, 10), kernel_size = c(1, 10), activation = "relu", padding = 'same', name = "downsize2_L1") %>%
    layer_batch_normalization() %>%
    layer_conv_2d(filters = 64, strides=c(1, 1), kernel_size = c(1, 5), activation = "relu", padding = 'same', name = "downsize2_L2") %>%
    layer_batch_normalization()
  
  profile_horizontal_model_downsize_centre <- profile_horizontal_model_downsize2 %>% 
    layer_conv_2d(filters = 128, strides=c(1, 50), kernel_size = c(1, 50), activation = "relu", padding = 'same', name = "downsize3") %>%
    layer_batch_normalization()
  
  profile_horizontal_upsize2 <- profile_horizontal_model_downsize_centre %>%
    layer_upsampling_2d(size = c(1, 50)) %>%
    {layer_concatenate(inputs = list(profile_horizontal_model_downsize2, .), axis = 3)} %>%
    layer_conv_2d(filters = 64, kernel_size = c(1, 50), padding = "same", activation = "relu") %>%
    layer_batch_normalization()
  
  profile_horizontal_upsize1 <- profile_horizontal_upsize2 %>%
    layer_upsampling_2d(size = c(1, 10)) %>%
    {layer_concatenate(inputs = list(profile_horizontal_model, .), axis = 3)} %>%
    layer_conv_2d(filters = 32, kernel_size = c(1, 10), padding = "same", activation = "relu") %>%
    layer_batch_normalization()
  
  profile_vertical_model <- profile_downsize %>% 
    layer_conv_2d(filters = layer1_number_of_filters, strides=c(1, 1), kernel_size = c(number_of_dyes, 1), activation = "relu", padding = 'valid', name = "profile_vertical_model_L1") %>%
    layer_batch_normalization()%>%
    layer_conv_2d(filters = layer2_number_of_filters, strides=c(1, 1), kernel_size = c(1, 1), activation = "relu", padding = 'valid', name = "profile_vertical_model_L2") %>%
    layer_batch_normalization()%>%
    layer_conv_2d(filters = 32, strides=c(1, 1), kernel_size = c(1, 1), activation = "relu", padding = 'same', name = "profile_vertical_model_L3") %>%
    layer_batch_normalization()
  
  
  #hold numbers 0 to 1
  number_input_layer <- layer_input(shape=c(number_of_dyes,500, 1), name="numbered_input_layer")
  
  # to change if number_of_dyes value changes
  duplicate_vertical_layer <- layer_concatenate(inputs = list(profile_vertical_model,
                                                              profile_vertical_model,
                                                              profile_vertical_model),
                                                axis = 1, name = "vertical_stack_creation")
  
#  conatenate_500_layers <- layer_concatenate(inputs = list(profile_horizontal_upsize1, duplicate_vertical_layer), axis = 3,  name = "final_concat")
  conatenate_500_layers <- layer_concatenate(inputs = list(profile_horizontal_upsize1, duplicate_vertical_layer, number_input_layer), axis = 3,  name = "final_concat")
  
  
  profile_500_input_CNN <- conatenate_500_layers %>%
    layer_conv_2d(filters = layer1_number_of_filters, strides=c(1, 1), kernel_size = c(1, 1), activation = "relu", padding = 'same', name = "profile_input_model_L1") %>%
    layer_batch_normalization()%>%
    layer_conv_2d(filters = layer2_number_of_filters, strides=c(1, 1), kernel_size = c(1, 1), activation = "relu", padding = 'same', name = "profile_input_model_L2") %>%
    layer_batch_normalization()%>%
    layer_conv_2d(filters = 32, strides=c(1, 1), kernel_size = c(1, 1), activation = "relu", padding = 'same', name = "profile_input_model_L3") %>%
    layer_batch_normalization()
  
  final_upsample <- profile_500_input_CNN %>%
    layer_upsampling_2d(size = c(1,10)) 
  
  output_layer <- layer_concatenate(inputs = list(profile_inputs, final_upsample), axis = 3) %>%
    layer_conv_2d(filters = layer1_number_of_filters, strides=c(1, 1), kernel_size = c(1, 5), activation = "relu", padding = 'same', name = "profile_output_L1") %>%
    layer_batch_normalization()%>%
    layer_spatial_dropout_2d(rate = 0.25) %>%
    layer_conv_2d(filters = layer2_number_of_filters, strides=c(1, 1), kernel_size = c(1, 5), activation = "relu", padding = 'same', name = "profile_output_L2") %>%
    layer_batch_normalization()%>%
    layer_spatial_dropout_2d(rate = 0.25) %>%
    layer_conv_2d(filters = 1, strides=c(1, 1), kernel_size = c(1, 5), 
                  padding = 'same', name = "profile_output_L3",
                  activation = "linear")  # Pas d'activation pour régression
  
  #puts it all together in the final model
#  complete_model <- keras_model(inputs = profile_inputs, outputs = output_layer)
  complete_model <- keras_model(inputs = list(profile_inputs, number_input_layer), outputs = output_layer)
  
  #summarises the model in a text
  summary(complete_model)
  
  # ---------------------------------------------------------------------
  
  #compile the model
  complete_model$compile(
    optimizer = optimizer_adam(learning_rate = 2e-4, beta_1 = 0.5),
    loss = "MeanSquaredError"
  )
  
  return(complete_model)
}