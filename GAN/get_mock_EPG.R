get_mock_EPG <- function(all_EPG_data, plot_EPG = FALSE, EPG_plot_savename = "EPG.jpg", context_reach = 0){
  
  source('convert_bp_to_scan.R')
  
  #locus names
  GF_locus_names <- c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "Yindel",     "AMEL",       "D8S1179",     "D21S11",    "D18S51",     "DYS391",     "D2S441", "D19S433", "TH01", "FGA", "D22S1045", "D5S818", "D13S317", "D7S820", "SE33", "size",    "D10S1248", "D1S1656", "D12S391", "D2S1338")
  all_GF_Dyes <- c("blue",     "blue",  "blue",    "blue",   "blue", "dark green", "dark green", "dark green", "dark green", "dark green", "dark green", "black", "black",  "black", "black", "red",      "red",    "red",     "red",    "red",  "orange", "purple",   "purple",  "purple",  "purple")
  locus_bp_size <- c(130,       180,    250,       305,      360,     80,          100,          145,           220,          290,          380,          85,       150,    200,      260,     110,        160,      220,       280,      360,    100,      110,        185,       240,       320)
  #matrix that links loci to dyes
  GFcolourLinkMatrix <- matrix(ncol=3, nrow=length(GF_locus_names))
  GFcolourLinkMatrix[,1] <- GF_locus_names
  GFcolourLinkMatrix[,2] <- all_GF_Dyes
  locus_scan_size <- convert_bp_to_scan(locus_bp_size)
  GFcolourLinkMatrix[,3] <- locus_scan_size
  colnames(GFcolourLinkMatrix) <- c("loci", "dye", "size")
  GFcolourLinkMatrix.df <- data.frame(GFcolourLinkMatrix)
  unique_dyes <- unique(all_GF_Dyes)
  number_of_dyes <- length(unique_dyes)
  #sets up x values
  startScan <- 1 - context_reach
  endScan <- 5000 + context_reach - 1
  plot_x <- seq(from = startScan, to = endScan, by = 1)
    
  #now adds in peaks for the ILS
  ILS_GF_peaks_bp <- c(60,80,100,114,120,140,160,180,200,214,220,240,250,260,280,300,314,320,340,360,380,400,414,420,440,460,480,500)
  ILS_GF_peaks_scans <- convert_bp_to_scan(ILS_GF_peaks_bp)
  number_ILS_peaks <- length(ILS_GF_peaks_bp)
  ILS_matirx_colnames <- c("Locus", "Allele", "Height", "Size", "Scans")
  ILS_add_matrix <- matrix(NA, ncol=length(ILS_matirx_colnames), nrow = number_ILS_peaks)
  colnames(ILS_add_matrix) <- ILS_matirx_colnames
  for (iILS_peak in 1:number_ILS_peaks){
    ILS_add_matrix[iILS_peak, 1] <- "size"
    ILS_add_matrix[iILS_peak, 2] <- ILS_GF_peaks_bp[iILS_peak]
    ILS_add_matrix[iILS_peak, 3] <- round(rnorm(n = 1, mean = 500, sd = 10), 0) #just some height around 500 with a little random noise (this may need to be more realistic to fool an AI discriminator)
    ILS_add_matrix[iILS_peak, 4] <- ILS_GF_peaks_bp[iILS_peak]
    ILS_add_matrix[iILS_peak, 5] <- ILS_GF_peaks_scans[iILS_peak]
  }
  
  all_EPG_data <- rbind(all_EPG_data, ILS_add_matrix)
  
  #the return array
  return_array <- matrix(ncol = number_of_dyes, nrow = (endScan - startScan + 1))
  colnames(return_array) <- unique_dyes
  
  #variable for peak width
  peak_width = 2
  #object to hold all plots
  mockEPGPlot <- vector('list', length(number_of_dyes))
  for (dye in 1:number_of_dyes){
    current_dye <- unique_dyes[dye]
    loci_in_dyeset <- GF_locus_names[all_GF_Dyes==current_dye]
    number_loci_in_dyeset <- length(loci_in_dyeset)
    
    #starts dataframe for plotting
    dye.df <- NULL
    #list of alleles
    allele_array <- NULL
    #list of peak scans
    scans_array <- NULL
    #sets up proflie y vaues
    plot_y <- rep(0, length(plot_x))
    for (locus in 1:number_loci_in_dyeset){
      current_locus <- loci_in_dyeset[locus]
      peaks_at_profile_and_locus <- all_EPG_data[all_EPG_data$Locus==current_locus,]$Allele
      #if there are no peaks at locus then leave all as zeros
      if(length(peaks_at_profile_and_locus)>0){
        #get data
        peak_heights_at_profile_and_locus <- as.numeric(all_EPG_data[all_EPG_data$Locus==current_locus,]$Height)
        peak_sizes_at_profile_and_locus <- as.numeric(all_EPG_data[all_EPG_data$Locus==current_locus,]$Scans)
        #add to tallies
        allele_array <- c(allele_array, peaks_at_profile_and_locus)
        scans_array <- c(scans_array, peak_sizes_at_profile_and_locus)
        
        for (peak in 1:length(peaks_at_profile_and_locus)){
          current_size <- peak_sizes_at_profile_and_locus[peak]
          if (current_size <= endScan){
            current_height <- peak_heights_at_profile_and_locus[peak]
            #multiplier so that rfu = height of graph peaks
            peak_height_add_array <- dnorm(plot_x, mean=as.numeric(current_size), sd=peak_width)
            multiplier <- current_height/max(peak_height_add_array)
            plot_y <- plot_y + multiplier*peak_height_add_array
          }
        }
      }
    }
    dye.df <- rbind(dye.df, data.frame(x = plot_x, y = plot_y))
    #stores the values
    return_array[,dye] <- as.numeric(dye.df$y)
    
    #if you want to plot the EPG
    if(plot_EPG){
      #plotting
      plot_max_rfu <- max(as.numeric(all_EPG_data$Height))
      ggplotbreaks <- scans_array
      ggplotlabels <- allele_array
      locus_label.df <- data.frame(ll_Locus = loci_in_dyeset, ll_Size = as.numeric(GFcolourLinkMatrix.df$size[GFcolourLinkMatrix.df$dye==current_dye]), ll_y=rep(plot_max_rfu, number_loci_in_dyeset))
      
      mockEPGPlot[[dye]] <- ggplot(data=dye.df) +
        geom_line(aes(x=x, y=y), col=current_dye, size=0.1) +
        theme(legend.position = "none", axis.text.x = element_text(size = 5)) +
        scale_x_continuous(breaks = as.numeric(ggplotbreaks), labels = ggplotlabels)+
        ylim(0, plot_max_rfu)+
        xlab("") +
        ylab("RFU")+
        geom_label(data=locus_label.df, aes(x=ll_Size, y=ll_y, label=ll_Locus))
    }
  }
  if (plot_EPG){
    save_name <- paste(here(), "\\EPGs\\", EPG_plot_savename, sep="")
    full_EPG_plot <- arrangeGrob(mockEPGPlot[[1]], mockEPGPlot[[2]], mockEPGPlot[[3]], mockEPGPlot[[4]],mockEPGPlot[[5]],mockEPGPlot[[6]], ncol=1)
    #print to screen
    grid::grid.draw(full_EPG_plot)
    #saves
    ggsave(save_name, full_EPG_plot, width = 15, height = 8, dpi = 600, units = "in")
  }
  return(t(return_array))
}