
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")

library(rhdf5)
path <- "/home/emilie/Documents/Fil_rouge/Code/Rproj/mmc1/EPG_classification_Keras_MHCNN_model.h5"
h5ls(path)
mydata <- h5read(path, )