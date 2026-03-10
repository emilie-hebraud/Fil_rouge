convert_bp_to_scan <- function(bp_array){
  
  #the size is generated in base pairs we we want it to be in scan points so that it matches the data in the EPG format
  bp_to_scan_mean <- 11.20068182
  
  #just  a fixed linear conversion
  Scans <- round(as.numeric(bp_array) * bp_to_scan_mean, 0) + 3500
  
  return(Scans)
}