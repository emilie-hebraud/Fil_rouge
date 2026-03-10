build_simplified_profile_reading_ANN <- function(
    n_channels    = 4,      # your 4 dye channels
    context_window = 100,   # same as original
    n_categories  = 3       # B, P, S
) {
  
  library(keras)
  library(tensorflow)
  
  # Input: (batch, n_channels, 1 + 2*context_window, 1)
  # = (batch, 4, 201, 1)
  input_layer <- layer_input(
    shape = c(n_channels, 1 + 2 * context_window, 1),
    name  = "input_layer"
  )
  
  # ── Shared layers ──────────────────────────────────────────
  
  # layer1: large conv across the context window
  # Original: 512 x 1 x 201 x 6 → we adapt for 4 channels
  shared_1 <- input_layer %>%
    layer_conv_2d(
      filters     = 512,
      kernel_size = c(1, 201),   # across full context window
      padding     = "valid",
      activation  = "relu",
      name        = "layer1"
    ) %>%
    layer_dropout(rate = 0.3)
  
  # layer2: 1x1 conv to compress features
  # Original: 256 x 512 x 1 x 1
  shared_2 <- shared_1 %>%
    layer_conv_2d(
      filters     = 256,
      kernel_size = c(1, 1),
      padding     = "valid",
      activation  = "relu",
      name        = "layer2"
    ) %>%
    layer_dropout(rate = 0.3)
  
  # ── Per-channel branches ────────────────────────────────────
  # Original had per-dye branches (FAM, VIC, NED, TAZ, LIZ, SID)
  # You have 4 channels → 4 branches
  
  channel_names <- c("channel_1", "channel_2", "channel_3", "channel_4")
  outputs <- vector("list", n_channels)
  
  for (i in 1:n_channels) {
    
    # layer3: per-channel intermediate layer
    # Original: 128 x 256 x 1 x 1
    branch <- shared_2 %>%
      layer_conv_2d(
        filters     = 128,
        kernel_size = c(1, 1),
        padding     = "valid",
        activation  = "relu",
        name        = paste0("layer3_", channel_names[i])
      ) %>%
      layer_dropout(rate = 0.3)
    
    # output: per-channel classification
    # Original FAM/VIC/NED: 5 categories → you use 3 (B, P, S)
    outputs[[i]] <- branch %>%
      layer_conv_2d(
        filters     = n_categories,
        kernel_size = c(1, 1),
        padding     = "valid",
        activation  = "softmax",
        name        = paste0("output_", channel_names[i])
      )
  }
  
  # Build model
  model <- keras_model(
    inputs  = input_layer,
    outputs = outputs,
    name    = "simplified_profile_reading_ANN"
  )
  
  # Compile 
  model$compile(
    optimizer = optimizer_adam(),
    loss = rep("categorical_crossentropy", n_channels),
    metrics = rep("accuracy", n_channels)
  )
  
  return(model)
}