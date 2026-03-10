smooth_profile <- function(raw_profile,
                           raw_profile_predictions,
                           probability_threshold){

  source('peak_detection.R')
  
  #the number of dyes
  number_dyes <- dim(raw_profile)[2] - 1
  #rfu threshold below which peaks will be screened out
  peak_RFU_threshold <- 0
  #makes a copy of the profile to return
  return_profile <- raw_profile
  channel_name <- c("channel_1", "channel_3", "channel_4")
  #goes through each dye
  for (idye in 2:number_dyes){
    
    print(paste("starting dye ", (idye - 1), sep=""))
    cat(sprintf("Dye %d: length=%d, NAs=%d, range=[%.1f, %.1f]\n",
                idye, 
                length(raw_profile[,idye]),
                sum(is.na(raw_profile[,idye])),
                min(raw_profile[,idye], na.rm=TRUE),
                max(raw_profile[,idye], na.rm=TRUE)))
    
    #identify peaks at certain thresholds and store position and height
    peak_probs_for_dye <- peak_detection(raw_profile[,idye], peak_RFU_threshold)
    
    #screen out non-peaks below probability threshold
    thinned_peak_probs_for_dye <- peak_probs_for_dye
    thinned_peak_probs_for_dye[thinned_peak_probs_for_dye < probability_threshold] <- 0
    #use a min-max function so that each peak is left with one central scan point
    scan_window <- 10

    for (iScan in 1:(length(thinned_peak_probs_for_dye) - scan_window)){
      max_prob_in_scan_window <- max(thinned_peak_probs_for_dye[iScan:(iScan + scan_window)])
      indices_of_non_maxes <- which(thinned_peak_probs_for_dye[iScan:(iScan + scan_window)] != max_prob_in_scan_window)
      thinned_peak_probs_for_dye[iScan:(iScan + scan_window)][indices_of_non_maxes] <- 0
    }
    
    
    #screen out values prior ro primer flare
    # thinned_peak_probs_for_dye[1:4000] <- 0
    #screen out any peak centres for which profile has below set fluorescence
    # thinned_peak_probs_for_dye[raw_profile[,idye] < peak_RFU_threshold] <- 0
    ch_name <- channel_name[idye - 1]
    thinned_peak_probs_for_dye[raw_profile_predictions[[ch_name]] == '0'] <- 0
    
    #now get heights for each peak position
    peak_centres_at_indices <- which(thinned_peak_probs_for_dye > 0)
    heights_of_peak_centres <- raw_profile[peak_centres_at_indices,idye]
    
    print(paste("found ", length(peak_centres_at_indices), "peaks", sep=""))
    #draws the peaks
    plot_y <- rep(0, length( raw_profile[,idye]))
    if (length(peak_centres_at_indices) > 0){
      #now create smoothed profile from these centres
      peak_width <- 4
      plot_x <- 1:length( raw_profile[,idye])
      for (peak in 1:length(peak_centres_at_indices)){
        current_size <- peak_centres_at_indices[peak]
        current_height <- heights_of_peak_centres[peak]
        #multiplier so that rfu = height of graph peaks
        peak_height_add_array <- dnorm(plot_x, mean = current_size, sd = peak_width)
        multiplier <- current_height/max(peak_height_add_array, na.rm = TRUE)
        plot_y <- plot_y + multiplier*peak_height_add_array
      }
    }
    #round to integer
    rounded_plot_y <- round(plot_y, 0)
    #add to the return profile
    return_profile[[idye]] <- rounded_plot_y
    
    print(paste("finished dye ", (idye - 1), sep=""))
  }
  #returns profile
  return(return_profile)
}

