#!/usr/bin/env Rscript

# FFERREIRA 03/24/2024
# Nina Project - Neanderthal

# Filters Differentially Expressed Genes (DEGs) by:
## 1) |L2FC| > 1 (mark if > 2)
## 2) FDR < 0.05

################
# 0. SETS UP ENV
################

# Sets WD
cat(paste("Setting WD...\n", sep = ""))
wd <- getwd()
setwd(wd)
cat(paste("\tDone!\n", sep = ""))

# General OBJs
cat(paste("Parsing ARGs...\n", sep = ""))
sep <- "\t"
nas <- c("NA", "", NA)
l2fc <- 1
fdr <- 0.05

# Gets CLI-ARGs
args <- commandArgs(trailingOnly = T)

# Checks CLI-ARGs
if (length(args) != 4)
{
  stop(paste("4 CLI-ARGs must be supplied: ",
             "(I) MAP file (ENSG --> Symbol & Type);",
             "(II) IODIR;",
             "(III) COMPARISON;",
             "(IV) BATCH.",
             sep = "\n"), call. = F)
}

# Processes CLI-ARGs
inFile <- args[1]  #e.g., "map_pcgene_IDs.txt" --> "map_gene_IDs.txt" (OLD / it has duplicated gene IDs)
IODir <- args[2]  #e.g., "deseq2"
myComp <- args[3]  #e.g., ASC-10_vs_ASC-CTRL | Others
myBatch <- args[4]  #e.g., "batch1-2" | batch1 | batch2

# Manual ARGs
# inFile <- "map_pcgene_IDs.txt"  # MAP file (ENSG --> Infos - symbol and type). "../map_gene_IDs.txt" has duplicated gene IDs, so it does not work
# IODir <- "deseq2"  # IODIR
# myComp <- "PVT_vs_Shaker"  #e.g., PVT_vs_Shaker | Others
# myBatch <- batch1-2  #e.g., "batch1-2" | batch1 | batch2

cat(paste("\tDone!\n", sep = ""))

#############################
# 1. EDITS AND WRITES RESULTS
#############################

# Reads MAP file between Ensembl Gene IDs and Symbols / Types
## PS: First COL must be set as "row.names"
cat(paste("Reading and editing INfile from 'DESeq2'. ATTENTION: Check if file is from comparison '", myComp, "'...\n",
          sep = ""))
map <- as.data.frame(read.delim(file = inFile, row.names = 1, header = F, sep = sep, na.strings = nas))

# Lists raw DESeq2 result files
files <- list.files(path = IODir, pattern = paste("result_DESeq2_", myBatch, "-", myComp, ".tsv", sep = ""))
files
# file <- files
cat(paste("\tDone!\n", sep = ""))

# Iterates over list of files
for (file in files)
{
  # Reads raw DESeq2 result file, one by one
  ## PS: No need to set "row.names" here!
  deseq2 <- as.data.frame(read.delim(file = paste(IODir, file, sep = "/"), header = T, sep = sep, na.strings = nas))
  
  # Maps gene symbols & types (MAP file) to raw results data (DESeq2 file)
  ## Sorts the OUT DF by the "by" COL
  ## PS: It is not interesting here to use the "all = T" ARG...
  data <- merge(x = map, y = deseq2, by = 0, sort = T)
  
  # Selects relevant COLs
  ## PS: fixed values for now!
  data <- data[,c(1:3,5,9)]
  
  # Renames COLs
  ## Ensembl Gene IDs and their respective Symbols, Types, Log2 Fold Change, and adjusted p-val (BH | FDR)
  colnames(data) <- c("ENSG", "Symbol", "Type", "L2FC", "FDR")
  
  # Removes NA values for FDR since they are not relevant for us in downstream analyses
  cat(paste("Applying filters: (i) |L2FC| > ", l2fc, " (mark if > 2) and FDR < ", fdr, "...\n", sep = ""))
  data <- data[!is.na(data[["FDR"]]),]
  
  # Filters only DEGs
  ## PS: |L2FC| > 1 & FDR < .05
  data <- data[((abs(data[["L2FC"]]) > l2fc) & data[["FDR"]] < fdr),]
  
  # Marks genes with |L2FC| > 2, it may be useful in the future!
  data[["L2FC.Mark"]] <- ifelse(test = (abs(data[["L2FC"]]) > 2), yes = "PASS", no = "|L2FC| <= 2")
  cat(paste("\tDone!\n", sep = ""))
  
  # Separates UP/DOWN-regulated genes
  up <- data[((data[["L2FC"]] > l2fc) & data[["FDR"]] < fdr),]  # UP
  down <- data[((data[["L2FC"]] < (-l2fc)) & data[["FDR"]] < fdr),]  # DOWN
  
  # Orders DFs by FDR (BH adjusted p-value)
  up <- up[order(up[["FDR"]]),]
  down <- down[order(down[["FDR"]]),]
  
  # Writes final (edited and filtered) DESeq2 results (UP/DOWN-regulated genes, separated) to files
  cat(paste("Writing OUT files...\n", sep = ""))
  write.table(x = data,
              file = paste(IODir, gsub(pattern = "result_DESeq2", replacement = "genes_ALL",
                                       x = file, ignore.case = F, fixed = T), sep = "/"),
              row.names = F, col.names = T, quote = F, sep = sep, na = nas)  # ALL
  write.table(x = up,
              file = paste(IODir, gsub(pattern = "result_DESeq2", replacement = "genes_UP",
                                       x = file, ignore.case = F, fixed = T), sep = "/"),
              row.names = F, col.names = T, quote = F, sep = sep, na = nas)  # UP
  write.table(x = down,
              file = paste(IODir, gsub(pattern = "result_DESeq2", replacement = "genes_DOWN",
                                       x = file, ignore.case = F, fixed = T), sep = "/"),
              row.names = F, col.names = T, quote = F, sep = sep, na = nas)  # DOWN
  cat(paste("\tDone!\n", sep = ""))
}
