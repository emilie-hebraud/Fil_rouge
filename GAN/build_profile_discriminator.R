build_profile_discriminator <- function(settings, useL1Regularization = TRUE){
  
  #libraries
  library('tensorflow')
  library('keras')
  library('tidyverse')
  
  number_of_dyes <- settings$number_of_dyes
  startScan <- settings$startScan
  number_scanpoint_in_input_layer <- settings$number_of_scanpoints
  print(paste("Building discriminator with", number_of_dyes, "dyes and", number_scanpoint_in_input_layer, "scanpoints"))
  
  # source image input [number_of_dyes x number_scanpoint_in_input_layer x 1]
  in_src_image <- layer_input(shape=c(number_of_dyes, number_scanpoint_in_input_layer, 1), name="input_src_layer")
  # target image input [number_of_dyes x number_scanpoint_in_input_layer x 1]
  in_target_image <- layer_input(shape=c(number_of_dyes, number_scanpoint_in_input_layer, 1), name="input_target_layer")
  # concatenate images channel-wise [number_of_dyes x number_scanpoint_in_input_layer x 2]
  merged_inputs <- layer_concatenate(list(in_src_image, in_target_image), axis = 3, name = "merged_input")
  
  if (useL1Regularization){
    # C64 - First convolution layer
    # Input: [number_of_dyes x number_scanpoint_in_input_layer x 2]
    # Output: [number_of_dyes x (number_scanpoint_in_input_layer/4) x 64]
    discrim_out <- merged_inputs %>%
      layer_conv_2d(
        filters = 64, 
        kernel_size = c(1, 10), 
        strides = c(1, 4), 
        padding = 'same', 
        kernel_regularizer =  keras$regularizers$L1(0.005), 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      ) %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # C128 - Second convolution layer
      # Output: [number_of_dyes x (number_scanpoint_in_input_layer/16) x 128]
      layer_conv_2d(
        filters = 128, 
        kernel_size = c(1, 10), 
        strides = c(1, 4), 
        padding = 'same', 
        kernel_regularizer =  keras$regularizers$L1(0.005), 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      ) %>%
      layer_batch_normalization() %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # C256 - Third convolution layer
      # Output: [number_of_dyes x (number_scanpoint_in_input_layer/64) x 256]
      layer_conv_2d(
        filters = 256, 
        kernel_size = c(1, 10), 
        strides = c(1, 4), 
        padding = 'same', 
        kernel_regularizer =  keras$regularizers$L1(0.005), 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      ) %>%
      layer_batch_normalization() %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # C512 - Fourth convolution layer
      # Output: [number_of_dyes x (number_scanpoint_in_input_layer/256) x 512]
      layer_conv_2d(
        filters = 512, 
        kernel_size = c(1, 10), 
        strides = c(1, 4), 
        padding = 'same', 
        kernel_regularizer =  keras$regularizers$L1(0.005), 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      ) %>%
      layer_batch_normalization() %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # Second last output layer - reduce spatial dimension across dyes
      # CRITICAL FIX: Use kernel_size matching number_of_dyes (3 instead of hardcoded 6)
      # Output: [1 x (number_scanpoint_in_input_layer/256) x 512]
      layer_conv_2d(
        filters = 512, 
        kernel_size = c(number_of_dyes, 1), 
        padding = 'valid'
      ) %>%
      layer_batch_normalization() %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # Patch output - final classification layer
      # Output: [1 x (number_scanpoint_in_input_layer/256) x 1]
      layer_conv_2d(
        filters = 1, 
        kernel_size = c(1, 1), 
        padding = 'same', 
        kernel_regularizer =  keras$regularizers$L1(0.005), 
        activation = 'sigmoid', 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      )
    
  } else {
    # C64 - First convolution layer (no L1 regularization version)
    discrim_out <- merged_inputs %>%
      layer_conv_2d(
        filters = 64, 
        kernel_size = c(1, 10), 
        strides = c(1, 4), 
        padding = 'same', 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      ) %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # C128 - Second convolution layer
      layer_conv_2d(
        filters = 128, 
        kernel_size = c(1, 10), 
        strides = c(1, 4), 
        padding = 'same', 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      ) %>%
      layer_batch_normalization() %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # C256 - Third convolution layer
      layer_conv_2d(
        filters = 256, 
        kernel_size = c(1, 10), 
        strides = c(1, 4), 
        padding = 'same', 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      ) %>%
      layer_batch_normalization() %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # C512 - Fourth convolution layer
      layer_conv_2d(
        filters = 512, 
        kernel_size = c(1, 10), 
        strides = c(1, 4), 
        padding = 'same', 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      ) %>%
      layer_batch_normalization() %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # Second last output layer - CRITICAL FIX: kernel_size matches number_of_dyes
      layer_conv_2d(
        filters = 512, 
        kernel_size = c(number_of_dyes, 1), 
        padding = 'valid'
      ) %>%
      layer_batch_normalization() %>%
      layer_activation_leaky_relu(alpha = 0.2) %>%
      
      # Patch output - final classification
      layer_conv_2d(
        filters = 1, 
        kernel_size = c(1, 1), 
        padding = 'same', 
        activation = 'sigmoid', 
        kernel_initializer = initializer_random_normal(mean = 0.0, stddev = 0.02)
      )
  }
  
  # define model
  discriminator <- keras_model(inputs = list(in_src_image, in_target_image), outputs = discrim_out)
  
  # summarises the model in a text
  summary(discriminator)
  
  # compile the model
  discriminator$compile(
    optimizer = optimizer_adam(learning_rate = 2e-4, beta_1 = 0.5),
    loss = "MeanSquaredError"
  )
  
  return(discriminator)
  
}