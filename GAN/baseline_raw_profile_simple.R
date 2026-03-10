baseline_raw_profile_simple <- function(raw_profile){
  source('getmode.R')
  
  number_dyes <- 4
  mode_startscan <- 1
  mode_endcan <- 5000
  
  return_profile <- raw_profile[, 1:(number_dyes+1)]
  
  #very basic baselining function
  for (idye in 2:(number_dyes+1)){
    # cat("\n=== Dye", idye-1, "===\n")
    values_subset <- raw_profile[[idye]][mode_startscan:mode_endcan]
    # cat("Range input:", range(values_subset, na.rm=TRUE), "\n")
    # cat("Any NA in input:", any(is.na(values_subset)), "\n")
    
    mode <- getmode(values_subset)
    # cat("Mode calculated:", mode, "\n")
    # cat("Mode is NA:", is.na(mode), "\n")
    
    result <- raw_profile[[idye]] - mode
    # cat("Range after subtraction:", range(result, na.rm=TRUE), "\n")
    # cat("Any NA after subtraction:", any(is.na(result)), "\n")
    
    return_profile[[idye]] <- result
  }
  
  return(return_profile)
}