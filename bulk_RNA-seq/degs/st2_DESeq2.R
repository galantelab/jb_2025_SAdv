#!/usr/bin/env Rscript

# FFERREIRA 03/24/2024
# Nina Project - Neanderthal

## Runs differential gene expression (DEG) analysis
### PS: it should have more than 1 sample by group for comparison (at least 2x2), otherwise it throws an ERROR
#### NICE TUTORIALS / EXPLANATIONS
# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html (DESeq2 - normalization)
# https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html (DESeq2 - DE analyses)

################
# 0. SETS UP ENV
################

# Loads LIBs
cat(paste("Loading libraries...\n", sep = ""))
library(DESeq2)
library(dplyr)
#library(icesTAF) # Only necessary to perform CLI tasks from R
cat(paste("\tDone!\n", sep = ""))

# Sets WD
cat(paste("Setting WD...\n", sep = ""))
wd <- getwd()
setwd(wd)
cat(paste("\tDone!\n", sep = ""))

# Gets CLI-ARGs
cat(paste("Parsing ARGs...\n", sep = ""))
args <- commandArgs(trailingOnly = T)

# Checks CLI-ARGs
if (length(args) != 6)
{
  stop(paste("6 CLI-ARGs must be supplied: ",
             "(I) CONDITION to be tested;",
             "(II) CONTROL group;",
             "(III) INDIR (after TXimport);",
             "(IV) OUTDIR (of DESeq2);",
             "(V) Control SAMPLE SIZE;",
             "(VI) My BATCH.",
             sep = "\n"), call. = F)
}

# Processes CLI-ARGs
group1 <- args[1]  #e.g., "ASC-10" (condition; not using '-' to avoid R message - see below)
group2 <- args[2]  #e.g., "ASC-CTRL" (control)
inDir <- args[3]  #e.g., "tximport"
outDir <- args[4]  #e.g., "deseq2"
nCtrl <- as.numeric(args[5])  #e.g., "5" (ls XXX/ | wc -l)
myBatch <- args[6]  #e.g., "batch1-2"

# Manual ARGs
# inDir <- "tximport"
# outDir <- "deseq2"
# group1 <- "PVT"  # condition (not using '-' to avoid R message - see below)
# group2 <- "Shaker"  # control
# nCtrl <- 5  # sample size of control group --> 1 (test) | 5 (all)
# myBatch <- "batch1-2"  #e.g., "batch1-2"

# TESTS
cat(paste("\tBatch: ", myBatch, "\n", sep = ""))
cat(paste("\tG1: ", group1, "\n", sep = ""))
cat(paste("\tG2: ", group2, "\n", sep = ""))
cat(paste("\tINDIR: ", inDir, "\n", sep = ""))
cat(paste("\tOUTDIR: ", outDir, "\n", sep = ""))
cat(paste("\tSample Size Ctrl: ", nCtrl, "\n", sep = ""))

cat(paste("\tDone!\n", sep = ""))

#############################
# 1. PREPARES DATA FOR DESEQ2
#############################

# Comparison
cat(paste("Checking if comparison is right...\n", sep = ""))
myComp <- paste(myBatch, "-", group1, "_vs_", group2, sep = "")

# Prints comparison --> IMPORTANT CHECKPOINT!
myComp
cat(paste("\tDone!\n", sep = ""))

# Reads input COUNTS data (from TXimport)
## "all_COUNTS_batch1-2-ASC-10_vs_ASC-CTRL.txt"
cat(paste("Reading and editing INfile from 'tximport'. ATTENTION here to 'match IDs from SHELL'...\n", sep = ""))
data <- as.data.frame(read.delim(file = paste(inDir, "/all_COUNTS_", myComp, ".txt", sep = ""),
                                 header = T, sep = "\t", dec = ".", check.names = F, na.strings = c(NA,"NA","")))

# Renames COLs: Sample IDs --> Groups
# Prints control IDs --> IMPORTANT CHECKPOINT!
## "PS:" ALWAYS create DIRs which are "COND_vs_CTRL"! Plus, make SYMLINK NAMES such that CTRL GROUP appears ALWAYS LAST after LS CMD!
control_IDs <- names(data[,c((ncol(data) - (nCtrl - 1)):ncol(data))])
control_IDs
colnames(data) <- ifelse(test = names(data) %in% control_IDs, yes = group2, no = group1)

# Relabels the COLs based on groups (adding a number TAG to turn them into unique IDs --> later use in "rownames")
# PS: Order: Condition (G1) --> Control (G2)!
df <- cbind(data[,grep(pattern = group1, x = colnames(data), ignore.case = F, fixed = T, value = F, invert = F)],
            data[,grep(pattern = group2, x = colnames(data), ignore.case = F, fixed = T, value = F, invert = F)])

# Rounds gene COUNTS (floats --> integers)
## PS: This is important because the "DESeqDataSetFromMatrix" "DESeq2" R function only accepts integers!
df <- round(x = df, digits = 0)
cat(paste("\tDone!\n", sep = ""))

# Creates vector of conditions (must be as "factor" to avoid WARN)
cat(paste("Establishing 'condition' formula to be 'tested in fact'...\n", sep = ""))
condition <- factor(x = c(rep(x = group1, table(names(data))[[group1]]),
                          rep(x = group2, table(names(data))[[group2]])),
                    levels = c(group1, group2))

# Prints CONDITION for FORMULA --> IMPORTANT CHECKPOINT!
condition

# Converts to DF and sets ROW names matching COL names from original DF
## PS: It should be in the same order of input data!
col_data <- as.data.frame(condition)
rownames(col_data) <- colnames(df)
cat(paste("\tDone!\n", sep = ""))

################
# 2. RUNS DESEQ2
################

cat(paste("Running 'DESeq2' itself...\n", sep = ""))

#################################### MESSAGE ####################################
# converting counts to integer mode
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
################################################################################

# Creates "DESeqDataSet" OBJ to be used in "DESeq" R FUN below
dataset <- DESeqDataSetFromMatrix(countData = df, colData = col_data, design = ~condition)

# Filters out low counts (only to speed up DESeq2 processing) by removing rows that have only "0 | 1 read"
dataset <- dataset[rowSums(counts(dataset)) > 1,]

# Applies "DESeq" R FUN to find DEGs
dds <- DESeq(object = dataset)
cat(paste("\tDone!\n", sep = ""))

#############################
# 3. EDITS AND WRITES RESULTS
#############################

cat(paste("ATTENTION: Checking results. 1st group appearing is 'condition'. 2nd one is the 'control'...\n", sep = ""))

# Extracts results from DESeq2 - contrast PAR specifies the order of comparison
## pAdjustMethod = "BH" ("fdr" seems to be an alias. controls FDR) --> less conservative corrections (?p.adjust)
## lfcThreshold = 0 (we apply quality filters later in "st3_Filters.R")
## PS: contrast = a CHAR vector with exactly 3 elements: the name of a factor in the design formula ("condition" here),
## the name of the NUMERATOR (condition) level for the FC (G1 = "tumor_general" here),
## and the name of the DENOMINATOR (control) level for the FC (G2 = "normal" here) (simplest case)
### https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
### The level given last is the base level (control) for the comparison
res <- results(object = dds, contrast = c("condition", group1, group2))

# PS: Sees the ordering of conditions comparison --> IMPORTANT CHECKPOINT!
## For example, if you see: "log2FoldChange  results log2 fold change (MAP): group1 vs group2"
## Then it means that positive log2FoldChange is up-regulated in group1
## and negative log2FoldChange is down-regulated in group1
order <- mcols(res,use.names = T)
order$description  # checks order of groups, statistical test used, correction method, etc...

# Re-orders results by Adjusted P-Value (BH: Benjamini & Hochberg, 1995 - "BH" or its alias "fdr")
res_ordered <- as.data.frame(res[order(res$padj),])
cat(paste("\tDone!\n", sep = ""))

# Writes final DEG DESeq2 results to outfile
cat(paste("Writing OUT files...\n", sep = ""))
write.table(x = res_ordered, file = paste(outDir, "/result_DESeq2_", myComp, ".tsv", sep = ""),
            quote = F, sep = "\t", row.names = T, col.names = T, na = c("NA", "", NA))
cat(paste("\tDone!\n", sep = ""))
