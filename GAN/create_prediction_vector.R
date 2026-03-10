# =============================================================================
# create_predictions_vector.R
#
# Creates a predictions vector from features and labels files.
#
# For each scan position (index) and each channel:
#   - "1" if within +/- 10 Da (molecular weight) of a main peak
#   - "0" otherwise (background)
#
# Input files:
#   features : index;channel_1;channel_2;channel_3;channel_4;channel_5;molw
#   labels   : plant_id,plant_id_str,multiplex,channel,marker_id,peak_kind,
#              peak_index,mu_pb,sigma_pb,peak_width_pb,amplitude,amp_class,
#              plant_peak_width_max_pb
#
# Link between files:
#   channel  = channel index (1-4) extracted from column name
#   plant_id = extracted from feature filename: *_pl{plant_id}*
#   multiplex= extracted from feature filename: *M{multiplex}_pl*
# =============================================================================

# HELPERS — extract plant_id and multiplex from filename

extract_plant_id <- function(filepath) {
  # Expects pattern: ..._pl{plant_id}...  e.g. features_M1_pl0001.csv
  fname <- tools::file_path_sans_ext(basename(filepath))
  as.integer(sub(".*_pl([0-9]+).*", "\\1", fname))
}

extract_multiplex <- function(filepath) {
  # Expects pattern: ...M{multiplex}_pl...  e.g. features_M1_pl0001.csv
  fname <- tools::file_path_sans_ext(basename(filepath))
  as.integer(sub(".*M([0-9]+)_pl.*", "\\1", fname))
}

# MAIN FUNCTION
create_predictions_vector <- function(cfg) {
  
  # Load data
  cat("Loading features from:", cfg$features_path, "\n")
  features <- read.table(cfg$features_path, sep = ";", header = TRUE,
                       stringsAsFactors = FALSE, row.names=NULL, fill=TRUE)
  cat("Loading labels from:  ", cfg$labels_path, "\n")
  # Read as raw lines to inspect around the problem
  raw <- readLines(cfg$labels_path)
  
  labels <- read.table(cfg$labels_path,
                       sep = ";",
                       header = TRUE,
                       dec = ".",
                       quote = "",
                       comment.char = "",
                       stringsAsFactors = FALSE,
                       row.names=NULL)
  # Extract plant_id and multiplex from filename 
  plant_idx <- extract_plant_id(cfg$features_path)
  multiplex  <- extract_multiplex(cfg$features_path)
  cat("Detected plant_id:", plant_idx, "| multiplex:", multiplex, "\n")
  
  # Initialise predictions with "0" (background) 
  n_scans     <- nrow(features)
  # print(n_scans)
  predictions <- data.frame(index = features$index,
                            molw  = features$molw,
                            stringsAsFactors = FALSE)
  
  for (ch in cfg$channels) {
    predictions[[ch]] <- rep("0", n_scans)
  }
  
  # Fill predictions channel by channel
  for (ch in cfg$channels) {
    
    ch_idx <- as.integer(sub("channel_", "", ch))   # "channel_2" → 2
    
    # Filter labels for this plant / multiplex / channel
    df_peaks <- labels[
      labels$channel   == ch_idx   &
        labels$plant_id  == plant_idx &
        labels$multiplex == multiplex,
    ]
    #print(df_peaks)
    cat(sprintf("  Channel %d: %d peaks found in labels\n",
                ch_idx, nrow(df_peaks)))
    
    if (nrow(df_peaks) == 0) next
    
    # For each scan position, determine label
    molw_vec <- features$molw   # molecular weight at each scan
    for (i in seq_len(nrow(df_peaks))) {
      
      mu        <- df_peaks$mu_pb[i]
      peak_kind <- tolower(trimws(df_peaks$peak_kind[i]))
      # print(class(df_peaks$mu_pb))
      # print(df_peaks$mu_pb)
      # print(class(features$molw))
      
      # ADD THIS
      label <- switch(peak_kind,
                      "main"      = "1",
                      "spurious"  = "0",
                      "0"  # fallback for unknown peak_kind
      )
      
      
      # Scan indices within the molecular weight window
      cat(sprintf("    mu=%.2f | molw range=[%.2f, %.2f]\n", 
                  mu, min(molw_vec, na.rm=TRUE), max(molw_vec, na.rm=TRUE)))
      in_window <- which(
        molw_vec >= (mu - cfg$mw_window) &
          molw_vec <= (mu + cfg$mw_window)
      )
      # print(in_window)
      for (idx in in_window) {
        if (predictions[[ch]][idx] != "1") {   # don't overwrite a main peak
          predictions[[ch]][idx] <- label
        }
      }
    }   # closes for (i in ...)
  }     # closes for (ch in cfg$channels)
  
  # Save AFTER both loops
  write.table(predictions,
              file      = cfg$output_path,
              sep       = ";",
              row.names = FALSE,
              col.names = TRUE,
              quote     = FALSE)
  cat("\nPredictions saved to:", cfg$output_path, "\n")
  return(predictions)
}         # closes create_predictions_vector


# QUICK VISUALISATION
plot_predictions_2 <- function(features, predictions,
                             channels = c("channel_1",
                                          "channel_3", "channel_4")) {
  library(ggplot2)
  library(gridExtra)
  
  # Helper: merge contiguous same-label scans into single rectangles
  get_label_regions <- function(df, target_label) {
    idx <- which(df$label == target_label)
    if (length(idx) == 0) return(data.frame(xmin = numeric(0), 
                                            xmax = numeric(0)))
    breaks  <- c(0, which(diff(idx) > 1), length(idx))
    data.frame(
      xmin = df$molw[idx[breaks[-length(breaks)] + 1]],
      xmax = df$molw[idx[breaks[-1]]]
    )
  }
  
  dye_cols  <- c("blue", "darkgreen", "goldenrod3", "red")
  plot_list <- vector("list", length(channels))
  
  for (i in seq_along(channels)) {
    ch <- channels[i]
    df <- data.frame(
      molw = features$molw,
      index  = features$index,
      signal = features[[ch]],
      label  = predictions[[ch]],
      stringsAsFactors = FALSE
    )
    
    # Get contiguous regions per label
    main_regions <- get_label_regions(df, "1")
    bg_regions <- get_label_regions(df, "0")
    
    # Dummy data for legend (uses aes fill mapping)
    legend_df <- data.frame(
      label = c("Main", "Background"),
      fill  = c("steelblue", "orange")
    )
    
    p <- ggplot(df, aes(x = molw, y = signal))
    
    # 2. Contiguous background bands
    if (nrow(main_regions) > 0) {
      p <- p + geom_rect(data = main_regions,
                         aes(xmin = xmin, xmax = xmax,
                             ymin = -Inf, ymax = Inf),
                         fill = "steelblue", alpha = 0.15,
                         inherit.aes = FALSE)
    }
    
    if (nrow(bg_regions) > 0) {
      p <- p + geom_rect(data = bg_regions,
                         aes(xmin = xmin, xmax = xmax,
                             ymin = -Inf, ymax = Inf),
                         fill = "orange", alpha = 0.15,
                         inherit.aes = FALSE)
    }
    
    # 3. Filled ribbon under curve — only within labelled regions
    #    Force continuity by setting non-matching labels to NA
    df_M        <- df
    df_M$signal <- ifelse(df$label == "1", df$signal, NA)
    
    df_S        <- df
    df_S$signal <- ifelse(df$label == "0", df$signal, NA)
    
    p <- p +
      geom_ribbon(data = df_M,
                  aes(ymin = 0, ymax = signal),
                  fill = "steelblue", alpha = 0.5,
                  na.rm = TRUE) +
      geom_ribbon(data = df_S,
                  aes(ymin = 0, ymax = signal),
                  fill = "orange", alpha = 0.5,
                  na.rm = TRUE)
    
    # 4. Signal line on top
    p <- p + geom_line(colour = dye_cols[i], linewidth = 0.4)
    
    # 5. Manual legend via annotate
    p <- p +
      annotate("rect", xmin = -Inf, xmax = -Inf,   # invisible — just for legend
               ymin = -Inf, ymax = -Inf,
               fill = NA, colour = NA) +
      labs(title   = paste(ch),
           x       = "Molecular weight (pb)",
           y       = "Signal",
           caption = "Blue = Main peak   |   Orange = Background") +
      theme_classic() +
      theme(plot.title   = element_text(size = 9, face = "bold"),
            plot.caption = element_text(size = 7, colour = "grey40"))
    
    plot_list[[i]] <- p
  }
  
  do.call("grid.arrange", c(plot_list, ncol = 1))
}

# Single file
#predictions <- create_predictions_vector(config)
#
# Visualise
#features <- read.table(config$features_path, sep=";", header=TRUE)
# Check molw range after filtering
#print(range(features$molw, na.rm=TRUE))

# Check if any rows are near the expected peak
#plot_predictions_2(features, predictions)