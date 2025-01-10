#Run QC filtering for each sample ##############################################################################################

#Load libraries
library(Seurat)
library(tidyverse)
library(scDblFinder)

#Create directory
dir.create("1_results_QC_filt", showWarnings=T)

#Create function to read unfiltered counts, create Seurat object and run QC filtering for each sample
run_qc <- function(folder, sample, mycolor){

  #Set filenames
  file_mtx = paste0(folder, "/counts_unfiltered/cells_x_genes.mtx")
  file_features = paste0(folder, "/counts_unfiltered/cells_x_genes.genes.names.txt")
  file_cells = paste0(folder, "/counts_unfiltered/cells_x_genes.barcodes.txt")

  #Create Seurat object, keeping genes expressed in at least 3 cells and cells with at least 200 genes expressed
  exp_matrix = ReadMtx(mtx=file_mtx, features=file_features, cells=file_cells, feature.column=1, mtx.transpose=T)
  df = CreateSeuratObject(counts=exp_matrix, project=sample, min.cells=3, min.features=200)
 
  #Get basic statistics before QC filtering
  beforeQC.number.cells=nrow(df@meta.data)
  beforeQC.median.genes.per.cell=median(df@meta.data$nFeature_RNA)
  beforeQC.median.UMIs.per.cell=median(df@meta.data$nCount_RNA)
  beforeQC.detected.genes=nrow(df[['RNA']]@layers$counts)

  #Calculate % of expressed mitochondrial genes per cell (indicator of cell quality) and add to metadata
  df[["percent.mt"]] = PercentageFeatureSet(df, pattern="^MT-")

  #Plot QC metrics as violin plots
  pdf(paste0("1_results_QC_filt/plot_metrics_beforeQC_", sample, ".pdf"), width=5, height=5)
  print(VlnPlot(df, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, cols=mycolor, alpha=0.1))
  garbage = dev.off()

  #Filter out cells with deviating from the range of 200 to 3,000 features and those harbouring over 5% of reads from mitochondrial genes
  metadata = df@meta.data
  metadata = metadata %>% mutate(out_mito = ifelse(percent.mt > 5, "TRUE", "FALSE")) %>%
                          mutate(out_feature = ifelse(nFeature_RNA < 200 | nFeature_RNA > 3000, "TRUE", "FALSE"))
  df$out_mito = metadata$out_mito
  df$out_feature = metadata$out_feature
  df$quality = ifelse(df$out_mito == TRUE | df$out_feature == TRUE, "Low-Quality", "High-Quality")
  df = subset(df, subset = quality == "High-Quality")
  
  #Get basic statistics after excluding low-quality cells
  afterQC.number.cells=nrow(df@meta.data)
  afterQC.median.genes.per.cell=median(df@meta.data$nFeature_RNA)
  afterQC.median.UMIs.per.cell=median(df@meta.data$nCount_RNA)
  afterQC.detected.genes=nrow(df[['RNA']]@layers$counts)
  
  #Identify doublets using scDblFinder and remove them from data
  set.seed(123)
  dbl = scDblFinder(as.SingleCellExperiment(df), returnType = 'table') %>% as.data.frame() %>% filter(type == 'real')
  dbl %>% dplyr::count(class)
  df = df[, dbl %>% dplyr::filter(class == "singlet") %>% rownames()]

  #Visualize QC metrics after subsetting
  pdf(paste0("1_results_QC_filt/plot_metrics_afterQC_", sample, ".pdf"), width=5, height=5)
  print(VlnPlot(df, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, cols=mycolor, alpha=0.3))
  garbage = dev.off()

  #Get basic statistics after doublets removal
  afterDoublets.number.cells=nrow(df@meta.data)
  afterDoublets.median.genes.per.cell=median(df@meta.data$nFeature_RNA)
  afterDoublets.median.UMIs.per.cell=median(df@meta.data$nCount_RNA)
  afterDoublets.detected.genes=nrow(df[['RNA']]@layers$counts)
  
  #Combine basic statistics
  number.cells = data.frame(metric="Number of cells", beforeQC=beforeQC.number.cells, afterQC=afterQC.number.cells, afterDoublets=afterDoublets.number.cells)
  median.genes.per.cell = data.frame(metric="Median genes per cell", beforeQC=beforeQC.median.genes.per.cell, afterQC=afterQC.median.genes.per.cell, afterDoublets=afterDoublets.median.genes.per.cell)
  median.UMIs.per.cell = data.frame(metric="Median UMIs per cell", beforeQC=beforeQC.median.UMIs.per.cell, afterQC=afterQC.median.UMIs.per.cell, afterDoublets=afterDoublets.median.UMIs.per.cell)
  detected.genes = data.frame(metric="Detected genes", beforeQC=beforeQC.detected.genes, afterQC=afterQC.detected.genes, afterDoublets=afterDoublets.detected.genes)
  summary = rbind(number.cells, median.genes.per.cell, median.UMIs.per.cell, detected.genes)

  #Write statistics to file
  write.table(summary, paste0("1_results_QC_filt/summary_", sample, ".txt"), row.names=F, col.names=T, quote=F, sep="\t")
  
  #Save filtered data to disk
  SaveSeuratRds(df, paste0("1_results_QC_filt/filt_", sample, ".rds"))

}

#Run QC filtering for each sample and save intermediate Seurat objects to disk
run_qc("230222-03-GM0216s1-organoid-C15-ctlr-10k-M", "C15-CTL", "lightgray")
run_qc("230222-04-GM0216s2-organoid-C15-10uM-10k-M", "C15-10uM", "lightgray")
run_qc("230222-05-GM0216s3-organoid-C15-30uM-10k-M", "C15-30uM", "lightgray")
run_qc("230222-06-GM0216s4-organoid-C27-ctlr-10k-M", "C27-CTL", "lightgray")
run_qc("230222-07-GM0216s5-organoid-C27-10uM-10k-M", "C27-10uM", "lightgray")
run_qc("230222-08-GM0216s6-organoid-C27-30uM-10k-M", "C27-30uM", "lightgray")
run_qc("230303-01-GM0302s1-organoid-ASC-CTRL-M", "ASC-CTL", "lightgray")
run_qc("230303-02-GM0302s2-organoid-ASC-10uM-M", "ASC-10uM", "lightgray")
run_qc("230303-03-GM0302s3-organoid-ASC-30uM-M", "ASC-30uM", "lightgray")
run_qc("230303-04-GM0302s4-organoid-C136-CTRL-M", "C136-CTL", "lightgray")
run_qc("230303-05-GM0302s5-organoid-C136-10uM-M", "C136-10uM", "lightgray")
run_qc("230303-06-GM0302s6-organoid-C136-30uM-M", "C136-30uM", "lightgray")
