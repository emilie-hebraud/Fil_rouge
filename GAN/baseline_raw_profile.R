baseline_raw_profile <- function(raw_profile){

  column_names <- colnames(raw_profile)
  #the columns with the fluorescence info are those that have names "dye#", whereas those with "dye#GT" hold classifications
  columns_with_dye <- grep(pattern="channel", x = column_names)
  columns_with_GT <- grep(pattern="GT", x = column_names)
  
  columns_holding_fluorescence_info <- NULL
  for (iCol in 1:length(columns_with_dye)){
    col_index_in_dye_array <- columns_with_dye[iCol]
     if(sum(col_index_in_dye_array == columns_with_GT) == 0){
       columns_holding_fluorescence_info <- c(columns_holding_fluorescence_info, columns_with_dye[iCol])
     }
  }
  
  number_dyes <- length(columns_holding_fluorescence_info)
  return_profile <- raw_profile[,c(1, columns_holding_fluorescence_info)]

  #a bit of a more sophisticated function that uses lowess
  for (idye in 1:number_dyes){
    dye_column_index <- columns_holding_fluorescence_info[idye]
    #gets the baseline
    lowess_smooth <- lowess(y = raw_profile[,dye_column_index], x = raw_profile[,1], f = 0.05)
    #removes it from the raw profile
    return_profile[[dye_column_index]] <- raw_profile[[dye_column_index]] - lowess_smooth$y
  }
  
  
  return(return_profile)
  
}
