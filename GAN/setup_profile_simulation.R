setup_profile_simulation <- function(){
  
  #load the size regression
  size_regression <- read_size_regression(paste(here(), "\\input_data\\size_regression\\GlobalFiler_SizeRegression.csv", sep=""))
  #locus names
  GF_locus_names <- c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "Yindel", "D8S1179", "D21S11", "D18S51", "DYS391", "D2S441", "D19S433", "TH01", "FGA", "D22S1045", "D5S818", "D13S317", "D7S820", "SE33", "D10S1248", "D1S1656", "D12S391", "D2S1338")
  #number of loci
  num_loci <- length(GF_locus_names)
  #repeat length for eahc locus
  GF_repeat_lengths <- c(D3S1358 = 4, vWA = 4, D16S539 = 4, CSF1PO = 4, TPOX = 4, D8S1179 = 4,
                         D21S11 = 4, D18S51 = 4, D2S441 = 4, D19S433 = 4, TH01 = 4, FGA = 4,
                         D22S1045 = 3, D5S818 = 4, D13S317 = 4, D7S820 = 4, SE33 = 4,
                         D10S1248 = 4, D1S1656 = 4, D12S391 = 4, D2S1338 = 4)
  
  #defining the stutter models
  
  #------back stutter------------
  # prepare stutter regression
  filename_bs_regression <- paste(here(), "\\input_data\\Stutter\\GlobalFiler_Back_Stutter.txt", sep="")
  bs_regression <- read_stutter_regression(filename_bs_regression)
  # prepare exceptions, i.e. where does the regression not apply?
  filename_bs_exceptions <- paste(here(), "\\input_data\\Stutter\\GlobalFiler_Stutter_Exceptions_3500.csv", sep="")
  bs_exceptions <- read_stutter_exceptions(filename_bs_exceptions)
  # prepare a stutter type
  backstutter <- stutter_type(name = "BackStutter", delta = -1,
                              stutter_regression = bs_regression,
                              stutter_exceptions = bs_exceptions)
  
  #------ forward stutter------------
  # prepare stutter regression
  filename_fs_regression <- paste(here(), "\\input_data\\Stutter\\GlobalFiler_Forward_Stutter.txt", sep="")
  fs_regression <- read_stutter_regression(filename_fs_regression)
  # prepare a stutter type
  forwardstutter <- stutter_type(name = "ForwardStutter", delta = +1,
                                 stutter_regression = fs_regression)
  
  #------ double back stutter------------
  # prepare stutter regression
  filename_dbs_regression <- paste(here(), "\\input_data\\Stutter\\GlobalFiler_Double_Back_Stutter_3500.txt", sep="")
  dbs_regression <- read_stutter_regression(filename_dbs_regression)
  # prepare a stutter type
  double.backstutter <- stutter_type(name = "DoubleBackStutter", delta = -2,
                                     stutter_regression = dbs_regression)
  
  #------ half stutter------------
  # prepare stutter regression
  filename_hs_regression <- paste(here(), "\\input_data\\Stutter\\GlobalFiler_Half_Stutter_3500.txt", sep="")
  hs_regression <- read_stutter_regression(filename_hs_regression)
  # prepare a stutter type
  halfstutter <- stutter_type(name = "HalfStutter",
                              delta = c(-1,2),
                              stutter_regression = hs_regression,
                              repeat_length_by_marker = GF_repeat_lengths,
                              applies_to_all_loci = FALSE,
                              applies_to_loci = c("SE33", "D1S1656"))
  
  #---------------------- set up model settings -------------------------------------
  
  dection_threshold <- setNames(rep(50, length(GF_locus_names)), GF_locus_names)
  
  # we create one allele specific stutter model with all stutter types
  
  allele_specific_stutter_model <- allele_specific_stutter_model(list(BackStutter = backstutter,
                                                                      ForwardStutter = forwardstutter,
                                                                      DoubleBackStutter = double.backstutter,
                                                                      HalfStutter = halfstutter), size_regression)
  
  # create the stutter variability model
  stutter_variability_model <- list(BackStutter = list(k2_prior = c(3.2662914731086,5.0150354973372),
                                                       inversely_proportional_to_parent = TRUE,
                                                       max_stutter_ratio = 0.3),
                                    ForwardStutter = list(k2_prior = c(4.5119034898389,1.2734535288209),
                                                          inversely_proportional_to_parent = TRUE,
                                                          max_stutter_ratio = 0.15),
                                    DoubleBackStutter = list(k2_prior = c(4.378255634018,1.4728182274552),
                                                             inversely_proportional_to_parent = TRUE,
                                                             max_stutter_ratio = 0.05),
                                    HalfStutter = list(k2_prior = c(3.3546980566024,1.7580535950464),
                                                       inversely_proportional_to_parent = TRUE,
                                                       max_stutter_ratio = 0.1)
  )
  #put them together in the model settings
  model_settings <- list(locus_names = GF_locus_names,
                         degradation_parameter_cap = 0.01,
                         c2_prior = c(8.4770608137371, 1.3193914013074),
                         LSAE_variance_prior = 0.0162351082765,
                         detection_threshold = dection_threshold,
                         size_regression = size_regression,
                         stutter_model = allele_specific_stutter_model,
                         stutter_variability = stutter_variability_model)
  
  return(model_settings)
  
}