get_mock_EPG_from_scans <- function(profile_data, startScan = 1, endScan = 4961){

  library('ggplot2')
  library('gridExtra')
  source('convert_bp_to_scan.R')
  
  # GF_locus_names <- c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "Yindel", "AMEL", "D8S1179", "D21S11", "D18S51", "DYS391", "D2S441", "D19S433", "TH01", "FGA", "D22S1045", "D5S818", "D13S317", "D7S820", "SE33", "size standard", "D10S1248", "D1S1656", "D12S391", "D2S1338")
  #matrix that links loci to dyes
  GFcolourLinkMatrix <- matrix(ncol=3, nrow=length(GF_locus_names))
  GFcolourLinkMatrix[,1] <- GF_locus_names
  all_GF_Dyes <- c("blue", "darkgreen","purple", "red")
  GFcolourLinkMatrix[,2] <- all_GF_Dyes
  locus_bp_size <- c(130,       180,    250,       305,      360,     80,          100,          145,          220,          290,          380,          85,       150,       200,      260,      110,        160,      220,       280,      360,    200,      110,        185,       240,       320)
  locus_scan_size <- convert_bp_to_scan(locus_bp_size)
  GFcolourLinkMatrix[,3] <- locus_scan_size
  colnames(GFcolourLinkMatrix) <- c("loci", "dye_", "size")
  GFcolourLinkMatrix.df <- data.frame(GFcolourLinkMatrix)
  unique_dyes <- unique(all_GF_Dyes)
  number_of_dyes <- length(unique_dyes)
  

  
  plot_max_rfu = max(as.numeric(unlist(profile_data[startScan:endScan, 2:(dim(profile_data)[2])])))
  plot_min_rfu = min(as.numeric(unlist(profile_data[startScan:endScan, 2:(dim(profile_data)[2])])))
  
  mockEPGPlot <- vector('list', length(number_of_dyes))
  
  for (dye in 1:number_of_dyes){
  
    current_dye <- unique_dyes[dye]
    loci_in_dyeset <- GF_locus_names[all_GF_Dyes==current_dye]
    
    
    mockEPGPlot[[dye]] <- ggplot(data = profile_data[startScan:endScan,]) +
                          geom_line(aes_string(x = "index", y = paste("channel_", dye, sep="")), col = current_dye, linewidth = 0.1) +
                          theme(legend.position = "none", axis.text.x = element_text(size = 5)) +
                          ylim(plot_min_rfu, plot_max_rfu) +
                          xlab("") +
                          ylab("RFU")
      
  }
  #put all the plots together
  full_EPG_plot <- do.call("grid.arrange", c(mockEPGPlot, ncol = 1))
  #return the plot
  return(full_EPG_plot)
  
}