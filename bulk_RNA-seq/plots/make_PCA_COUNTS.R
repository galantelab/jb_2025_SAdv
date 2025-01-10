#!/usr/bin/env Rscript

# FFERREIRA 06/11/2024
# BEPE Collabs - Nina Neanderthals

## Makes PCA plots from COUNTS data generated in 'st1_TXimport.R'
### Apply filters: (i) protein-coding only; (ii) |L2FC| > 1 | 2

################
# 0. SETS UP ENV
################

# Loads LIBs
library(ggplot2)
library(DESeq2)

# Sets WD
wd <- getwd()
setwd(wd)

# Gets CLI-ARGs
args <- commandArgs(trailingOnly = T)

# Checks CLI-ARGs
if (length(args) != 8)
{
  stop(paste("8 CLI-ARGs must be supplied: ",
             "(I) INDIR;",
             "(II) OUTDIR;",
             "(III) CONDITION to be tested;",
             "(IV) CONTROL group;",
             "(V) Control SAMPLE SIZE;",
             "(VI) TXImport type;",
             "(VII) L2FC filter;",
             "(VIII) Gene type filter.",
             sep = "\n"), call. = F)
}

# Processes CLI-ARGs
inDir <- args[1]  # TXimport/
outDir <- args[2]  # Plots/
group1 <- args[3]  #e.g., Condition; not using '-' to avoid R message - see below
group2 <- args[4]  #e.g., Control
nCtrl <- as.numeric(args[5])  #e.g., 2 for all
txi <- args[6]  #e.g., COUNTS | TPM
l2fc_filter <- args[7]  #e.g., |L2FC| < [1,2] | PASS
gene_filter <- args[8]  #e.g., protein_coding | many others

# Manual ARGs
## GROUPS: ASC-10, ASC-30, ASC-CTRL, C136-10, C136-30, C136-CTRL,
## GROUPS: NOVA1-C15-10, NOVA1-C15-30, NOVA1-C15-CTRL, NOVA1-C27-10, NOVA1-C27-30, NOVA1-C27-CTRL
group1 <- "NOVA1-ArAr-CTRL"  # LEFT (condition; not using '-' to avoid R message - see below)
group2 <- "NOVA1-HuHu-CTRL"  # RIGHT (control)
inDir <- "TXimport"  # tximport/
outDir <- "Plots"  # figures/
nCtrl <- 4  # It seems for all comparisons there were just duplicates
txi <- "COUNTS"  #e.g., COUNTS | TPM
l2fc_filter <- "PASS"  #e.g., |L2FC| < 2* | PASS; * Here 1, cause Nina's data is more specific
gene_filter <- "protein_coding"  #e.g., protein_coding | many others

##########################
# 1. PREPARES DATA FOR PCA
##########################

# Reads TXimport raw data
df <- as.data.frame(read.delim(file = paste(inDir, "/all_", txi, "_batch2-", group1, "_vs_", group2, ".txt", sep = ""),
                               header = T, check.names = F, dec = ".", sep = "\t", na.strings = c(NA,"NA","")))

# Keeps only genes with counts >= 5 in any sample
## PS: more than half of the genes are excluded here
### https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-DESeq2/#:~:text=Filter%20raw%20counts,through%20and%20will%20be%20faster.
### In the above link it says the recommendation is to remove genes with less than 10, not 5 reads/counts...
df <- df[apply(df, 1, function(x){any(x >= 5)}),]

# Rounds gene expression (floats --> integers)
## PS: This is important because the "DESeqDataSetFromMatrix (DESeq2)" R function only accepts integers!
df <- round(x = df, digits = 0)

# Reads DESeq2 filtered data
## L2FC info
map <- as.data.frame(read.delim(file = paste(wd, "/Inputs/DEGs_", group1, "_vs_", group2, ".ALL.final.tsv", sep = ""),
                                header = T, check.names = F, dec = ".", sep = "\t", na.strings = c(NA,"NA","")))

# Reads conversion table file MAPping IDs data
## ENSG --> Symbol
if (gene_filter == "all")
{
  id2symbol <- as.data.frame(read.delim(file = paste(wd, "/map_gene_IDs.tsv", sep = ""), header = F,
                                        check.names = F, dec = ".", sep = "\t", na.strings = c(NA,"NA","")))
} else {  # PC-only
  id2symbol <- as.data.frame(read.delim(file = paste(wd, "/map_gene_IDs.PC.tsv", sep = ""), header = F,
                                        check.names = F, dec = ".", sep = "\t", na.strings = c(NA,"NA","")))
}

# Renames id2symbol COLs
colnames(id2symbol) <- c("ID", "Name")

# Adds annotation (symbol) info while keeping only the PC-genes, if previously selected
data <- merge(x = map, y = id2symbol, by.x = "Gene", by.y = "ID", all = F, sort = T)

# Creates "L2FC.Mark" COL for filtering later
## PS: changed criteria here to use an acceptable amount of genes!
data[["L2FC.Mark"]] <- ifelse(test = abs(data[["L2FC"]]) < 1, yes = "|L2FC| < 1", no = "PASS")

# Filters data for |L2FC| >= [1,2] only
## PS: many other genes are removed here, again!
metadata <- data[data[["L2FC.Mark"]] == l2fc_filter,]

# Removes ROWs with any NA
## PS: some other genes are removed here as well!
metadata <- na.omit(metadata)

# Renames COLs
colnames(metadata) <- c("Gene.ID", "AVG.L2.COUNTS", "L2FC", "FDR", "Gene.Symbol", "L2FC.Mark")

# Adds symbolic "Type" COL to metadata
metadata[["Type"]] <- gene_filter

# Gets the IDs of the selected genes before calling plotPCA
idx.meta <- metadata[["Gene.ID"]]

# Creates vector of conditions (must be as "factor" to avoid WARN)
## PS: here, for the CTRL group (REF), the samples are always the last ones
condition <- factor(x = c(rep(group1, (length(colnames(df))) - nCtrl), rep(group2, nCtrl)), levels = c(group1, group2))

# Converts to DF and sets ROW names matching COL names from original DF
## PS: It should be in the same order of input data!
col_data <- as.data.frame(condition)
rownames(col_data) <- colnames(df)

################
# 2. RUNS DESEQ2
################

# Creates "DESeqDataSet" object to be used in "DESeq2" R functions below
## PS: ALWAYS let DESeq2 run using ALL GENES --> https://support.bioconductor.org/p/106025/#106032
dataset <- DESeqDataSetFromMatrix(countData = df, colData = col_data, design = ~condition)

# Applies "DESeq2" normalization
dds <- estimateSizeFactors(dataset)
dds <- estimateDispersions(dds)

# Log2-transforms data
rld <- vst(dds, blind = F)
colnames(rld) <- colnames(df)

##################
# 3. PLOTS FIGURES
##################

# Plots PCA
## |L2FC| > 2 (PASS) & Protein-Coding

### PDF
pdf(file = paste(outDir, paste("PCA_", group1, "_vs_", group2, ".L2FC-PC.pdf", sep = ""), sep = "/"),
    width = 8, height = 10, onefile = T, title = "PCAs", bg = "transparent",
    fg = "black", pointsize = 12, pagecentre = T, useKerning = T, compress = T)

## Directs PCA to OBJ
myPCA <- plotPCA(object = rld[rownames(rld) %in% idx.meta,]) + ggtitle("") + theme_classic() +
  scale_colour_brewer(palette = "Accent", direction = -1) +
  theme(legend.text = element_text(size = 14, face = "bold"), legend.title = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "plain", color = "black")) +
  scale_color_manual(values = c("#EAC435", "#345995"), breaks = c(group1, group2),
                     labels = c(paste(group1, "(Condition)", sep = " "), paste(group2, "(Control)", sep = " ")))

# Plots PCA
print(myPCA)

# Devs off (closes OUTFILE)
dev.off()

### PNG
png(filename = paste(outDir, paste("PCA_", group1, "_vs_", group2, ".L2FC-PC.png", sep = ""), sep = "/"),
    width = 20, height = 25, res = 300, units = "cm",
    family = "ArialMT", title = "Neanderthals NOVA1", bg = "white", pointsize = 12)

## Directs PCA to OBJ
myPCA <- plotPCA(object = rld[rownames(rld) %in% idx.meta,]) + ggtitle("") + theme_classic() +
  scale_colour_brewer(palette = "Accent", direction = -1) +
  theme(legend.text = element_text(size = 14, face = "bold"), legend.title = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "plain", color = "black")) +
  scale_color_manual(values = c("#EAC435", "#345995"), breaks = c(group1, group2),
                     labels = c(paste(group1, "(Condition)", sep = " "), paste(group2, "(Control)", sep = " ")))

# Plots PCA
print(myPCA)

# Devs off (closes OUTFILE)
dev.off()
