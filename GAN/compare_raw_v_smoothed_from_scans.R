compare_raw_v_smoothed_from_scans <- function(profile_data,
                                              smoothed_profile,
                                              sample_name_with_path,
                                              startScan_profile = 1,
                                              endScan_profile = 5000,
                                              startScan_smoothed = 1,
                                              endScan_smoothed = 5000,
                                              save_jpg = TRUE,
                                              startScan_real = 1,
                                              endScan_real = 5000,
                                              generator = FALSE){
  
  library(ggplot2)
  library(patchwork)
  
  final_column <- dim(profile_data)[2]-1
  plot_max_rfu <- max(max(as.numeric(unlist(profile_data[startScan_profile:endScan_profile, 2:final_column]))),
                      max(as.numeric(unlist(smoothed_profile[startScan_profile:endScan_profile, 2:final_column]))))
  plot_min_rfu <- min(min(as.numeric(unlist(profile_data[startScan_profile:endScan_profile, 2:final_column]))),
                      min(as.numeric(unlist(smoothed_profile[startScan_profile:endScan_profile, 2:final_column]))))
  
  mockEPGPlot <- list()
  
  for (dye in c(1,3,4))#channel_2 is missing
    {
    print(dye)
    plot_df <- data.frame(index = numeric(0), value = numeric(0), label = character(0))
    
    if(!is.null(profile_data)){
      profile_label <- ifelse(exists("generator") && generator, "generated", "raw")
      plot_df <- rbind(plot_df, data.frame(
        index = profile_data[startScan_profile:endScan_profile, "molw"],
        value = profile_data[startScan_profile:endScan_profile, paste0("channel_", dye)],
        label = paste0(profile_label, " - channel_", dye)
      ))
    }
    
    if(!is.null(smoothed_profile)){
      plot_df <- rbind(plot_df, data.frame(
        index = smoothed_profile[startScan_smoothed:endScan_smoothed, "molw"],
        value = smoothed_profile[startScan_smoothed:endScan_smoothed, paste0("channel_", dye)],
        label = paste0("smoothed - channel_", dye)
      ))
    }
    
    plot_df$label <- factor(plot_df$label, levels = unique(plot_df$label))
    
    p <- ggplot(plot_df, aes(x = index, y = value, color = label)) +
      geom_line(linewidth = 0.1) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 5)) +
      ylim(plot_min_rfu, plot_max_rfu) +
      xlab("") +
      ylab("RFU") +
      labs(color = "Signal type & channel")
    
    mockEPGPlot[[length(mockEPGPlot) + 1]] <- p
  }
  
  # Combine with patchwork
  full_EPG_plot <- wrap_plots(mockEPGPlot, ncol = 1)
  
  if(save_jpg){
    save_name <- paste0(sample_name_with_path, ".jpg")
    ggsave(save_name, full_EPG_plot, width = 15, height = 8, dpi = 600, units = "in")
  }
}