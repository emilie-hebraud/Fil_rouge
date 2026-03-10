peak_detection <- function(baselined_profile, peak_RFU_threshold){
  
  source('sumLogProb.R')
  source('getmode.R')
  
  normdist <- function(x, SD) { #currently on a ln scale
    return(-0.5 * (x / SD)^2 - log(SD * sqrt(2.0 * pi)))
  }
  
  #the window over which the peak detection is carried out
  scanRange <- 20
  
  #parameters to trial
  sigmaMax <- 7
  sigmaMin <- 2
  stepFactor <- 0.5
  sigmaPrior <- 1.0 / (stepFactor * (sigmaMax - sigmaMin))
  meanPrior <- 1.0 / (scanRange - 1)
  #standard deviation of last 1000 scans
#  NoiseSD <- sd(baselined_profile[(length(baselined_profile)-1000):length(baselined_profile)]) #this is not working smoothly so I am going to use magic number of 3.8
  NoiseSD <- 3.8
  
  #array to hold the peak label probs
  peakProb <- rep(0, length(baselined_profile))
  for (iScan in (scanRange / 2):(length(baselined_profile) - scanRange / 2)) {
      NumeratorProb <- 0
      DenominatorProbNoPeak <- 0
      DenominatorProbPeak <- 0
      #set up region of interest
      ROI <- baselined_profile[(iScan - scanRange / 2):(iScan + scanRange / 2)]

      #baseline ROI based on mode
      ROI <- ROI - getmode(ROI)
      
      #check to see if the central peak is greater than the threshold, otherwise skip it
      if (ROI[round(scanRange/2, 0)] > peak_RFU_threshold){
      
        #no peak
        no_peak_tally = 0;
        for (iWindow in 1:scanRange) {
          no_peak_tally <- no_peak_tally + normdist(ROI[iWindow], NoiseSD)
        }
        DenominatorProbNoPeak <- no_peak_tally
        
        #peak at different POI
        for (mean in 1:scanRange) {
          A <- ROI[mean]
          sigma_trial <- sigmaMin
          while (sigma_trial <= sigmaMax) {
            peak_tally <- 0
            if (mean < scanRange / 2 - 1 | mean > scanRange / 2 + 1) { #not the central peak
              for (window in 1:scanRange) {
                expected <- A * exp(-0.5 * ((window - mean) / sigma_trial)^2)
                peak_tally <- peak_tally + normdist(expected - ROI[window], NoiseSD)
              }
              if (mean == 1 & sigma_trial == sigmaMin) { #first one
                DenominatorProbPeak <- peak_tally
              } else {
                DenominatorProbPeak <- sumLogProb(DenominatorProbPeak, log(meanPrior*((scanRange - 3)/scanRange) * sigmaPrior) + peak_tally) #in log scale
              }
            }
            sigma_trial <- sigma_trial + 1 / stepFactor
          }
        }
        
        #peak at POI
        for (mean in (round(scanRange / 2, 0) - 1):(round(scanRange / 2, 0) + 1)) {
          A <- ROI[mean]
          sigma_trial <- sigmaMin
          while (sigma_trial <= sigmaMax) {
            tally <- 0
            A <- ROI[mean]
            for(window in 1:scanRange) {
              expected <- A * exp(-0.5 * ((window - mean) / sigma_trial)^2)
              tally <- tally + normdist(expected - ROI[window], NoiseSD)
            }
            if (sigma_trial == sigmaMin) {
              NumeratorProb <- tally
            } else {
              NumeratorProb <- sumLogProb(NumeratorProb, log(1/3 * sigmaPrior) + tally) #in log scale
            }
            sigma_trial <- sigma_trial + 1 / stepFactor
          }
        }
        #calculate the posterior peak probability
        peakProb[iScan] <- NumeratorProb - sumLogProb(log(0.5) + DenominatorProbNoPeak, log(0.5) + DenominatorProbPeak)
        
      } # end not > threshold peak
  }
  
  return(peakProb)
}


