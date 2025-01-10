#Make plots (part 2) ###############################################################################################

#Libraries
library(Seurat)
library(tidyverse)
library(RColorBrewer)

#Load mapped Seurat data with reductions
df = LoadSeuratRds("3_integrated_samples/integrated_with_reductions.rds")

#Create column with samples in metadata
df$sample = gsub("_1|_2", "", df$orig.ident)

#Change sample names
df$sample <- gsub("ASC", "NOVA1-hu/hu", df$sample)
df$sample <- gsub("C136", "NOVA1-ar/ar", df$sample)
df_subset <- subset(df, subset = sample == "NOVA1-ar/ar-CTL" | sample == "NOVA1-ar/ar-10uM" | sample == "NOVA1-hu/hu-CTL" | sample == "NOVA1-hu/hu-10uM" )

## Plot proportion / cell number composition per sample
# Prep
ggData = data.frame(prop.table(table(df_subset$predicted.id, df_subset$sample), margin = 2))
colnames(ggData) = c("predicted.id", "sample", "value")
nClust <- length(unique(ggData$predicted.id))
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)

# Control vs Control
ggData = ggData[ggData$sample == "NOVA1-ar/ar-CTL" | ggData$sample == "NOVA1-hu/hu-CTL",]

p1 <- ggplot(ggData, aes(sample, value, fill = predicted.id)) +
  geom_col() + xlab("Sample") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colCls) + coord_flip() +
  labs(fill = "Cell Type")

png(filename="4.1_rmercuri_plots/plot_celltype_control_vs_control.png", type="cairo", units ="in", width=6, height=4, pointsize=12, res=96)
p1
garbage = dev.off()

# C136
ggData = data.frame(prop.table(table(df_subset$predicted.id, df_subset$sample), margin = 2))
colnames(ggData) = c("predicted.id", "sample", "value")
ggData = ggData[ggData$sample == "NOVA1-ar/ar-CTL" | ggData$sample == "NOVA1-ar/ar-10uM",]
ggData$sample <- factor(ggData$sample, levels = c("NOVA1-ar/ar-10uM", "NOVA1-ar/ar-CTL"))

p1 <- ggplot(ggData, aes(sample, value, fill = predicted.id)) +
  geom_col() + xlab("Condition") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colCls) + coord_flip() +
  scale_x_discrete(labels=c("10uM","CTL")) +
  labs(fill = "Cell Type") + 
  ggtitle("NOVA1-ar/ar") +
  theme(plot.title = element_text(hjust = 0.5))

png(filename="4.1_rmercuri_plots/plot_celltype_C136.png", type="cairo", units ="in", width=6, height=4, pointsize=12, res=96)
p1
garbage = dev.off()

# ASC
ggData = data.frame(prop.table(table(df_subset$predicted.id, df_subset$sample), margin = 2))
colnames(ggData) = c("predicted.id", "sample", "value")
ggData = ggData[ggData$sample == "NOVA1-hu/hu-CTL" | ggData$sample == "NOVA1-hu/hu-10uM",]
ggData$sample <- factor(ggData$sample, levels = c("NOVA1-hu/hu-10uM", "NOVA1-hu/hu-CTL"))

p1 <- ggplot(ggData, aes(sample, value, fill = predicted.id)) +
  geom_col() + xlab("Condition") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colCls) + coord_flip() +
  scale_x_discrete(labels=c("10uM","CTL")) +
  labs(fill = "Cell Type") + 
  ggtitle("NOVA1-hu/hu") +
  theme(plot.title = element_text(hjust = 0.5))

png(filename="4.1_rmercuri_plots/plot_celltype_ASC.png", type="cairo", units ="in", width=6, height=4, pointsize=12, res=96)
p1
garbage = dev.off()

## uMAPs
#Set default assay
DefaultAssay(df_subset) <- "RNA"

df_subset$sample = factor(df_subset$sample, levels = c("NOVA1-hu/hu-CTL", "NOVA1-hu/hu-10uM","NOVA1-ar/ar-CTL", "NOVA1-ar/ar-10uM"))

#FOXP2 with max.cutoff=q50
png(filename="4.1_rmercuri_plots/plot_UMAP_FOXP2_splitBySample.png", type="cairo", units ="in", width=6, height=6, pointsize=12, res=96)
print(FeaturePlot(df_subset, features="FOXP2", split.by="sample", ncol=1, max.cutoff = "q50"))+
    patchwork::plot_layout(ncol = 2, nrow = 2)
garbage = dev.off()

#LXH1 with max.cutoff=q50
png(filename="4.1_rmercuri_plots/plot_UMAP_LHX1_splitBySample.png", type="cairo", units ="in", width=6, height=6, pointsize=12, res=96)
print(FeaturePlot(df_subset, features="LHX1", split.by="sample", ncol=1, max.cutoff = "q50"))+
    patchwork::plot_layout(ncol = 2, nrow = 2)
garbage = dev.off()

#TCF7L2 with max.cutoff=q50
png(filename="4.1_rmercuri_plots/plot_UMAP_TCF7L2_splitBySample.png", type="cairo", units ="in", width=6, height=6, pointsize=12, res=96)
print(FeaturePlot(df_subset, features="TCF7L2", split.by="sample", ncol=1, max.cutoff = "q50"))+
    patchwork::plot_layout(ncol = 2, nrow = 2)
garbage = dev.off()

##Bubble plot
colGEX = c("grey85", brewer.pal(7, "Reds"))
genes.to.plot <- c("FOXP2", "GBX2", "LHX1", "LHX5", "OTX2", "TCF7L2")

#By cells
p1 <- DotPlot(df_subset, group.by = "predicted.id", features = genes.to.plot, scale = FALSE) + 
  coord_flip() + scale_color_gradientn(colours = colGEX) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  xlab("") + ylab("")
  
png(filename="4.1_rmercuri_plots/plot_markers_splitBycells.png", type="cairo", units ="in", width=8, height=4, pointsize=12, res=96)
p1
garbage = dev.off()

#By cells from NOVA1-hu/hu
df_subset$cluster_join <- paste(df_subset$predicted.id, df_subset$sample, sep = "-")
p1 = DotPlot(df_subset, features = genes.to.plot, group.by = "cluster_join", scale = FALSE) +
  coord_flip() + scale_color_gradientn(colours = colGEX) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(limits = c("EN1-NOVA1-hu/hu-CTL","EN1-NOVA1-hu/hu-10uM","EN2-NOVA1-hu/hu-CTL","EN2-NOVA1-hu/hu-10uM","IN1-NOVA1-hu/hu-CTL","IN1-NOVA1-hu/hu-10uM","IN2-NOVA1-hu/hu-CTL","IN2-NOVA1-hu/hu-10uM","IN3-NOVA1-hu/hu-CTL","IN3-NOVA1-hu/hu-10uM","IPC2-NOVA1-hu/hu-CTL","IPC2-NOVA1-hu/hu-10uM","MES-NOVA1-hu/hu-CTL","MES-NOVA1-hu/hu-10uM","RG1-NOVA1-hu/hu-CTL","RG1-NOVA1-hu/hu-10uM","RG2-NOVA1-hu/hu-CTL","RG2-NOVA1-hu/hu-10uM"), labels=c("EN1-CTL","EN1-10uM","EN2-CTL","EN2-10uM","IN1-CTL","IN1-10uM","IN2-CTL","IN2-10uM","IN3-CTL","IN3-10uM","IPC2-CTL","IPC2-10uM","MES-CTL","MES-10uM","RG1-CTL","RG1-10uM","RG2-CTL","RG2-10uM")) +
    ggtitle("NOVA1-hu/hu") +
  xlab("") + ylab("")

png(filename="4.1_rmercuri_plots/plot_markers_splitByASC.png", type="cairo", units ="in", width=10, height=4, pointsize=12, res=96)
p1
garbage = dev.off()


#By cells from NOVA1-ar/ar
df_subset$cluster_join <- paste(df_subset$predicted.id, df_subset$sample, sep = "-")
p1 = DotPlot(df_subset, features = genes.to.plot, group.by = "cluster_join", scale = FALSE) +
  coord_flip() + scale_color_gradientn(colours = colGEX) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(limits = c("EN1-NOVA1-ar/ar-CTL","EN1-NOVA1-ar/ar-10uM","EN2-NOVA1-ar/ar-CTL","EN2-NOVA1-ar/ar-10uM","IN1-NOVA1-ar/ar-CTL","IN1-NOVA1-ar/ar-10uM","IN2-NOVA1-ar/ar-CTL","IN2-NOVA1-ar/ar-10uM","IN3-NOVA1-ar/ar-CTL","IN3-NOVA1-ar/ar-10uM","IPC2-NOVA1-ar/ar-CTL","IPC2-NOVA1-ar/ar-10uM","MES-NOVA1-ar/ar-CTL","MES-NOVA1-ar/ar-10uM","RG1-NOVA1-ar/ar-CTL","RG1-NOVA1-ar/ar-10uM","RG2-NOVA1-ar/ar-CTL","RG2-NOVA1-ar/ar-10uM"), labels=c("EN1-CTL","EN1-10uM","EN2-CTL","EN2-10uM","IN1-CTL","IN1-10uM","IN2-CTL","IN2-10uM","IN3-CTL","IN3-10uM","IPC2-CTL","IPC2-10uM","MES-CTL","MES-10uM","RG1-CTL","RG1-10uM","RG2-CTL","RG2-10uM")) +
  ggtitle("NOVA1-ar/ar") +
  xlab("") + ylab("")

png(filename="4.1_rmercuri_plots/pplot_markers_splitByC136.png", type="cairo", units ="in", width=10, height=4, pointsize=12, res=96)
p1
garbage = dev.off()

#By cells from control samples
p1 = DotPlot(df_subset, features = genes.to.plot, group.by = "cluster_join", scale = FALSE) +
  coord_flip() + scale_color_gradientn(colours = colGEX) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(limits = c("EN1-NOVA1-ar/ar-CTL","EN1-NOVA1-hu/hu-CTL","EN2-NOVA1-ar/ar-CTL","EN2-NOVA1-hu/hu-CTL","IN1-NOVA1-ar/ar-CTL","IN1-NOVA1-hu/hu-CTL","IN2-NOVA1-ar/ar-CTL","IN2-NOVA1-hu/hu-CTL","IN3-NOVA1-ar/ar-CTL","IN3-NOVA1-hu/hu-CTL","IPC2-NOVA1-ar/ar-CTL","IPC2-NOVA1-hu/hu-CTL","MES-NOVA1-ar/ar-CTL","MES-NOVA1-hu/hu-CTL","RG1-NOVA1-ar/ar-CTL","RG1-NOVA1-hu/hu-CTL","RG2-NOVA1-ar/ar-CTL","RG2-NOVA1-hu/hu-CTL"), labels=c("EN1-NOVA1-ar/ar","EN1-NOVA1-hu/hu","EN2-NOVA1-ar/ar","EN2-NOVA1-hu/hu","IN1-NOVA1-ar/ar","IN1-NOVA1-hu/hu","IN2-NOVA1-ar/ar","IN2-NOVA1-hu/hu","IN3-NOVA1-ar/ar","IN3-NOVA1-hu/hu","IPC2-NOVA1-ar/ar","IPC2-NOVA1-hu/hu","MES-NOVA1-ar/ar","MES-NOVA1-hu/hu","RG1-NOVA1-ar/ar","RG1-NOVA1-hu/hu","RG2-NOVA1-ar/ar","RG2-NOVA1-hu/hu")) +
  ggtitle("CTL Condition") +
  xlab("") + ylab("") 
  
png(filename="4.1_rmercuri_plots/plot_markers_splitByCTL.png", type="cairo", units ="in", width=10, height=4, pointsize=12, res=96)
p1
garbage = dev.off()

#By samples
p1 <- DotPlot(df_subset, group.by = "sample", features = genes.to.plot, scale = FALSE) + 
  coord_flip() + scale_color_gradientn(colours = colGEX) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  xlab("") + ylab("")
  
png(filename="4.1_rmercuri_plots/plot_markers_splitBysamples.png", type="cairo", units ="in", width=6, height=5, pointsize=12, res=96)
p1
garbage = dev.off()
