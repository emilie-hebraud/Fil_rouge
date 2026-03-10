get_mock_EPG_from_scans <- function(profile_data, startScan = 1, endScan = 4961){
  
  library('ggplot2')
  library('gridExtra')
  
  # ==================== ADAPT THIS SECTION TO YOUR KIT ====================
  # Define your locus names based on your DNA profiling kit
  # Example for a 4-dye system - MODIFY ACCORDING TO YOUR ACTUAL KIT
  
  # Option 1: If you know your kit's loci, replace these with your actual locus names
  locus_names <- c(
    "Locus1", "Locus2", 
    "Locus3", "Locus4", 
    "Locus5","Locus6", 
    "Locus7", "Locus8"        
  )
  
  # Define which dye/channel each locus belongs to
  # Use channel names that match your data columns
  all_dyes <- c(
    "blue", "blue",           # Loci in channel 1
    "green", "green",        # Loci in channel 2
    "yellow", "yellow",     # Loci in channel 3
    "red", "red"            # Loci in channel 4
  )
  
  # ==================== END OF ADAPTATION SECTION ====================

  # Get unique dyes (channels)
  unique_dyes <- unique(all_dyes)
  number_of_dyes <- length(unique_dyes)
  
  # Calculate plot limits
  plot_max_rfu <- max(as.numeric(unlist(profile_data[startScan:endScan, 2:(dim(profile_data)[2])])), na.rm = TRUE)
  plot_min_rfu <- min(as.numeric(unlist(profile_data[startScan:endScan, 2:(dim(profile_data)[2])])), na.rm = TRUE)
  
  # Create list to hold plots
  mockEPGPlot <- vector('list', number_of_dyes)
  
  # Generate a plot for each dye/channel
  for (dye in 1:number_of_dyes) {
    
    current_dye <- unique_dyes[dye]
    loci_in_dyeset <- locus_names[all_dyes == current_dye]
    
    # Get the correct column name from your data
    # Your data has columns: scan, channel_1, channel_2, channel_3, channel_4
    column_name <- paste("dye_", dye, sep = "_")
    
    # Check if column exists, if not use generic dye naming
    if (!column_name %in% colnames(profile_data)) {
      column_name <- paste("dye_", dye, sep = "")
    }
    
    mockEPGPlot[[dye]] <- ggplot(data = profile_data[startScan:endScan, ]) +
      geom_line(aes_string(x = "scan", y = column_name), 
                col = current_dye, 
                linewidth = 0.1) +
      theme(legend.position = "none", 
            axis.text.x = element_text(size = 5)) +
      ylim(plot_min_rfu, plot_max_rfu) +
      xlab("") +
      ylab("RFU") +
      theme_classic()
    
  }
  
  # Combine all plots
  full_EPG_plot <- do.call("grid.arrange", c(mockEPGPlot, ncol = 1))
  
  # Return the plot
  return(full_EPG_plot)
}
