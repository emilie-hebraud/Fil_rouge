convert_xml_to_csv <- function(folder_holding_xmls){

  #loads required library
  library("XML")
  
  #set up working dir
  workingDir <- folder_holding_xmls
  
  #for loop across all xml files in worwing directory
  SampleNames <- list.files(path=workingDir, pattern="\\.xml")
  NoOfSamples <- length(SampleNames)
  
  for (s in 1:NoOfSamples){
    
    #input filename
    inputFilename <- SampleNames[s]
    #inputFilename minus the ".xml" suffix
    saveFilename <- substr(inputFilename,1,nchar(inputFilename)-4)
    
    #save data as datframe
    inputdata <- xmlParse(paste(workingDir, inputFilename, sep=""))
    
    #convert to list
    xml_data <- xmlToList(inputdata)
    
    #data from dyes
    scans <- length(as.numeric(xml_data[[20]][[1]][[2]]))
    dyes <- length(xml_data[[20]])
    #big matrix that will hold all the informaiton : [dye][scan+1]
    CompiledMatrix <- matrix(NA, nrow = scans, ncol = dyes+1)
    
    for (dye in 1:dyes){
      CompiledMatrix[,dye+1] <- as.numeric(xml_data[[20]][[dye]][[2]])
    }
    
    #names columns in Matirx
    colnames(CompiledMatrix) <- c("scan", sprintf("dye%s", seq(1:dyes)))
    #adds scan number ot first column
    CompiledMatrix[,1] <- seq(1:scans)
    
    #saves as table
    write.csv(CompiledMatrix, file=paste(workingDir, saveFilename, ".csv", sep=""), row.names=FALSE, col.names=TRUE)
  }
}