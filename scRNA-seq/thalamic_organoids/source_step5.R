#Make plots (part 1) ###############################################################################################

#Load libraries
library(Seurat)
library(tidyverse)

#Create directory to store plots
dir.create("4_plots", showWarnings=T)

#Load mapped Seurat data with reductions
df = LoadSeuratRds("3_integrated_samples/integrated_with_reductions.rds")

#Subset Seurat object
df = subset(df, idents = c("ASC-10uM_1", "ASC-10uM_2", "ASC-CTL_1", "ASC-CTL_2", "C136-10uM_1", "C136-10uM_2", "C136-CTL_1", "C136-CTL_2"))

#Create column with final sample names in metadata
metadata = df@meta.data
metadata = metadata %>% rownames_to_column(var="cell.id")
samples = read.delim("map_samples.txt", header=F)
colnames(samples) = c("orig.ident", "sample")
metadata = merge(metadata, samples, by="orig.ident") %>% column_to_rownames(var="cell.id")
df@meta.data = metadata

#Set order of samples
df@meta.data$sample = factor(df@meta.data$sample, levels=c("NOVA1 hu/hu CTL", "NOVA1 hu/hu 10uM", "NOVA1 ar/ar CTL", "NOVA1 ar/ar 10uM"))

#Plot UMAP grouped by samples
pdf("4_plots/plot_UMAP_samples.pdf", width=6, height=4)
print(DimPlot(df, reduction="umap", group.by="sample", label=F) + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2"))
garbage = dev.off()

#Plot UMAP splitted by samples
pdf("4_plots/plot_UMAP_splitBySample.pdf", width=4.2, height=4)
print(DimPlot(df, reduction="umap", group.by="sample", split.by="sample", ncol=2, label=F) + NoLegend() + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2"))
garbage = dev.off()

#Plot UMAP grouped by predicted cell types
pdf("4_plots/plot_UMAP_predictedCellTypes.pdf", width=6, height=4)
print(DimPlot(df, reduction="umap", group.by="predicted.id", label=F) + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2"))
garbage = dev.off()

#Plot UMAP grouped by predicted cell types and splitted by samples
pdf("4_plots/plot_UMAP_predictedCellTypes_splitBySample.pdf", width=6, height=4)
print(DimPlot(df, reduction="umap", group.by="predicted.id", split.by="sample", ncol=2, label=F) + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2"))
garbage = dev.off()

#Set default assay
DefaultAssay(df) <- "RNA"

#Plot UMAP with LMNA expression (all samples together)
pdf("4_plots/plot_UMAP_FOXP2.pdf", width=6, height=4)
print(FeaturePlot(df, features="FOXP2") + xlab("UMAP 1") + ylab("UMAP 2"))
garbage = dev.off()

#Plot UMAP with LMNA expression (splitted by samples)
pdf("4_plots/plot_UMAP_FOXP2_splitBySample.pdf", width=16, height=4)
print(FeaturePlot(df, features="FOXP2", split.by="sample"))
garbage = dev.off()
