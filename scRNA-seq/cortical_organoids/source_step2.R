#Normalize samples ##############################################################################################

#Load libraries
library(Seurat)
library(tidyverse)

#Create directory to store results of normalized samples
dir.create("2_results_norm", showWarnings=T)

#Create function to normalize samples
normalizeData <- function(input, output){

  #Load filtered data
  df = LoadSeuratRds(input)

  #Normalize data using standard normalization
  df = NormalizeData(df, scale.factor=10000)
  df = FindVariableFeatures(df, selection.method="vst", nfeatures=2000)
  df = ScaleData(df, features=rownames(df))

  #Get cell cycle scores per cell (it should be done after normalization)
  df = CellCycleScoring(df, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=F)
  df$CC.Difference = df$S.Score - df$G2M.Score

  #Save normalized Seurat object
  SaveSeuratRds(df, output)

}

#Run pipeline using standard normalization -----------------------------------------------------------------

#List input files
input_files = list.files(path="1_results_QC_filt", pattern=".rds", full.names=T)

#For each input file
for(input_file in input_files){

  #Set output file
  output_file = gsub("1_results_QC_filt/filt_", "2_results_norm/norm_", input_file)

  #Run normalization
  normalizeData(input_file, output_file)

}

