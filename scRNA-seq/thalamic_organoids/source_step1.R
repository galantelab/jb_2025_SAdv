#Run QC filtering for each sample ############################################################################################################

#Load libraries
library(Seurat)
library(tidyverse)
library(scDblFinder)

#Create directory
dir.create("1_results_QC_filt", showWarnings=T)

#Create function to read counts, create Seurat object and run QC filtering for each sample
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
    #Distribution of number of expressed genes (nFeature_RNA), number of UMIs (nCount_RNA), and % mitochondrial genes expressed (percent.mt) per cell
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
  afterDoublets.number.cells=nrow(df@meta.data) #number of cells
  afterDoublets.median.genes.per.cell=median(df@meta.data$nFeature_RNA) #number of expressed genes per cell
  afterDoublets.median.UMIs.per.cell=median(df@meta.data$nCount_RNA) #number of UMIs per cell
  afterDoublets.detected.genes=nrow(df[['RNA']]@layers$counts) #expressed genes
  
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
run_qc("221021-04-GM1011S4-C136-CTL1-10k-M", "C136-CTL_1", "lightgray")
run_qc("221103-04-GM1011S4-C136-CTL1-10k-M", "C136-CTL_2", "lightgray")
run_qc("221021-05-GM1011S5-C136-10uM-10k-M", "C136-10uM_1", "lightgray")
run_qc("221103-05-GM1011S5-C136-10uM-10k-M", "C136-10uM_2", "lightgray")
run_qc("221021-01-GM1011S1-ASC1-13-CTL1-10k-M", "ASC-CTL_1", "lightgray")
run_qc("221103-01-GM1011S1-ASC1-13-CTL1-10k-M", "ASC-CTL_2", "lightgray")
run_qc("221021-02-GM1011S2-ASC1-13-10uM-10k-M", "ASC-10uM_1", "lightgray")
run_qc("221103-02-GM1011S2-ASC1-13-10uM-10k-M", "ASC-10uM_2", "lightgray")

