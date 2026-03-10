library('here')
library('ggplot2')
library('gridExtra')

profile_dirs <- c(paste(additonal_profile_dirs_path, "GF_3500_RUN_14975\\", sep=""),
                  paste(additonal_profile_dirs_path, "GF_3500_RUN_14999\\", sep=""),
                  paste(additonal_profile_dirs_path, "GF_3500_RUN_15015\\", sep=""),
                  paste(additonal_profile_dirs_path, "GF_3500_RUN_15022\\", sep=""),
                  paste(additonal_profile_dirs_path, "GlobalFiler 3500 mixtures\\", sep=""),
                  paste(additonal_profile_dirs_path, "GlobalFiler 3500 refs\\", sep=""),
                  paste(additonal_profile_dirs_path, "GlobalFiler_3500_2021_lowlevel\\", sep=""),
                  paste(additonal_profile_dirs_path, "GlobalFiler_3500_PROVEDIt\\", sep=""),
                  paste(additonal_profile_dirs_path, "GlobalFiler_3500_PROVEDIt_mixtures\\", sep="")
)



#read in source
source('baseline_raw_profile.R')
source('baseline_raw_profile_simple.R')
source('smooth_profile.R')
source('compare_raw_v_smoothed_from_scans.R')

startScan_profile = 4000
endScan_profile = 9000
startScan_smoothed = 4000
endScan_smoothed = 9000
startScan_real = 4000
endScan_real = 9000

#get number of directories that hold files
number_of_directories <- length(profile_dirs)
#loop over each directory
for (iDir in 1:number_of_directories){
  #get names of cvs samples
  SampleNames <- list.files(path = profile_dirs[iDir], pattern="\\.csv")
  #get number of samples
  NoOfSamples <- length(SampleNames)
  #loop over each sample
  for (iSample in 1:NoOfSamples){
    #gets the sample
    sampleName <- SampleNames[iSample]

    #only work on non-smoothed profiles
    got_smooth <- grep(pattern="smooth", x = sampleName)
    got_jpg <- grep(pattern="jpg", x = sampleName)
    if(length(got_smooth) == 0 & length(got_smooth) == 0) {
      
      #gets the raw profile path
      raw_profile_path <- paste(profile_dirs[iDir], sampleName, sep="")
      
      #reads in the raw profile data
      raw_profile <- read.csv(file = raw_profile_path, sep=",", header = TRUE, stringsAsFactors = FALSE)
      
      #creates the simple baselined raw profile
      profile_data <- baseline_raw_profile_simple(raw_profile)
      
      #loads the smoothe profile
      sample_name_without_suffix <- substr(x = sampleName, start = 1, stop = nchar(sampleName) - 4)
      smooth_profile_path <- paste(profile_dirs[iDir], sample_name_without_suffix, "_smooth.csv", sep="")
      smooth_profile <- read.csv(file = smooth_profile_path, sep=",", header = TRUE, stringsAsFactors = FALSE)
  
      #plots a simple baselined vs smoothed (this shows more what the ANN will learn to emulate)
      final_column <- (dim(profile_data)[2])
      plot_max_rfu = 300
      plot_min_rfu = min(min(as.numeric(unlist(profile_data[startScan_profile:endScan_profile, 2:final_column]))),
                         min(as.numeric(unlist(smooth_profile[startScan_profile:endScan_profile, 2:final_column]))))
      
      number_of_dyes <- dim(smooth_profile)[2] - 1
      
      mockEPGPlot <- vector('list', length(number_of_dyes))
      for (dye in 1:number_of_dyes){
        p <- ggplot()
        p <- p + geom_line(data = profile_data[startScan_profile:endScan_profile,], aes_string(x = "scan", y = paste("dye", dye, sep="")), col = "orange", linewidth = 0.1)
        p <- p + geom_line(data = smooth_profile[startScan_smoothed:endScan_smoothed,], aes_string(x = "scan", y = paste("dye", dye, sep="")), col = "black", linewidth = 0.1)
        p <- p + theme(legend.position = "none", axis.text.x = element_text(size = 5))
        p <- p + ylim(plot_min_rfu, plot_max_rfu)
        p <- p + xlab("")
        p <- p + ylab("RFU")
        p <- p + theme_classic()
        mockEPGPlot[[dye]] <- p
      }
      full_EPG_plot <- do.call("grid.arrange", c(mockEPGPlot, ncol = 1))
      save_name <- paste(profile_dirs[iDir], sample_name_without_suffix, "_zoom.jpg", sep="")
      ggsave(save_name, full_EPG_plot, width = 15, height = 8, dpi = 600, units = "in")
    }
    
    
  } # next sample
} # next directory
