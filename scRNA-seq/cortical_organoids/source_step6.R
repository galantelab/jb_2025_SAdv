#Make plots ###############################################################################################

#Load libraries
library(Seurat)
library(tidyverse)

#Create directory to store plots
dir.create("4_plots_mixed_1.5_2M", showWarnings=T)

#Load mapped Seurat data with reductions
df = LoadSeuratRds("3_integrated_samples/integrated_mixed_1.5_2M_with_reductions.rds")

#Create column with final sample groups in metadata
metadata = df@meta.data
metadata = metadata %>% rownames_to_column(var="cell.id")
groups = read.delim("sample_groups.txt", header=T)
metadata = merge(metadata, groups, by.x="orig.ident", by.y="sample") %>% column_to_rownames(var="cell.id")
df@meta.data = metadata

#Set order of groups
df@meta.data$group = factor(df@meta.data$group, levels=c("NOVA1 hu/hu CTL", "NOVA1 hu/hu 10uM", "NOVA1 hu/hu 30uM", "NOVA1 ar/ar CTL", "NOVA1 ar/ar 10uM", "NOVA1 ar/ar 30uM"))

#Plot UMAP grouped by groups
pdf("4_plots_mixed_1.5_2M/plot_UMAP_groups.pdf", width=6, height=4)
print(DimPlot(df, reduction="umap", group.by="group", label=F) + ggtitle(""))
garbage = dev.off()

#Plot UMAP splitted by group
pdf("4_plots_mixed_1.5_2M/plot_UMAP_splitByGroup.pdf", width=8, height=4)
print(DimPlot(df, reduction="umap", group.by="group", split.by="group", ncol=3, label=F) + ggtitle(""))
garbage = dev.off()

#Plot UMAP grouped by predicted cell types
pdf("4_plots_mixed_1.5_2M/plot_UMAP_predictedCellTypes.pdf", width=6, height=4)
print(DimPlot(df, reduction="umap", group.by="predicted.id", label=F) + ggtitle(""))
garbage = dev.off()

#Plot UMAP grouped by predicted cell types and splitted by groups
pdf("4_plots_mixed_1.5_2M/plot_UMAP_predictedCellTypes_splitByGroup.pdf", width=8, height=4)
print(DimPlot(df, reduction="umap", group.by="predicted.id", split.by="group", ncol=3, label=F) + ggtitle(""))
garbage = dev.off()
