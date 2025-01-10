#!/usr/bin/env Rscript

# FFERREIRA 03/24/2024
# Nina Project - Neanderthal

## Aggregates gene-level COUNTS / get TPM values

################
# 0. SETS UP ENV
################

# Loads libraries
cat(paste("Loading libraries...\n", sep = ""))
library(tximport)
library(readr)
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
if (length(args) != 5)
{
  stop(paste("5 CLI-ARGs must be supplied: ",
             "(I) Transcript to Gene conversion table;",
             "(II) INDIR;",
             "(III) OUTDIR;",
             "(IV) COMPARISON;",
             "(V) BATCH.",
             sep = "\n"), call. = F)
}

# Processes CLI-ARGs
## Use ID2Gene.txt (Kallisto ID) | ENST-2-ENSG.txt (use 'ignoreAfterBar = T' then)
trans2gene <- args[1]  # ../ID2Gene.txt (Transcript --> Gene conversion table) | ENST-2-ENSG.txt (OLD)
inDir <- args[2]  # ../kallisto/batch/Condition_vs_Control (INDIR with Kallisto files)
outDir <- args[3]  # tximport (OUTDIR where outputs will be placed)
myComp <- args[4]  # Lead_vs_Nothing | Others (my comparisons)
myBatch <- args[5]  # batch1-2 | batch1 | batch2 (batch used)

# Manual ARGs
## Use ID2Gene.txt (Kallisto ID) | ENST-2-ENSG.txt (use 'ignoreAfterBar = T' then)
# trans2gene <- "../ID2Gene.txt"
# myComp <- "PVT_vs_Shaker"
# inDir <- paste(wd, "/../kallisto/", myBatch, "/", myComp, sep = "")
# outDir <- "tximport"
# myBatch <- "batch1-2"

cat(paste("\tDone!\n", sep = ""))

###############################
# 1. PREPARES DATA FOR TXIMPORT
###############################

# Gets files to be processed
## If there is no 'H5' files available, then use the regular 'TSV' ones
cat(paste("Reading and editing files...\n", sep = ""))
files <- list.files(path = inDir, pattern = "abundance.tsv", recursive = T)

# Gets sample IDs
samples <- unlist(strsplit(x = files, split = "/"))[c(seq(from = 2, to = (length(files)*3), by = 3))]
names(files) <- samples

# Prints named files --> IMPORTANT CHECKPOINT!
files

# Reads ID2GENE file
tx2gene <- as.data.frame(read.delim(file = trans2gene, header = F, sep = "\t", na.strings = c(NA,"NA",""), dec = "."))

# Changes WD only to get TXIMPORT to work
setwd(inDir)
cat(paste("\tDone!\n", sep = ""))

##################
# 2. RUNS TXIMPORT
##################

# TXImports files
## TXimport can read TSV files, but it is slower...
## LINK: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto
# txi.kallisto <- tximport(files[1], type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
cat(paste("Running 'tximport' itself...\n", sep = ""))
txi.kallisto <- tximport(files, tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM",  # TXimport ARGs
                         type = "kallisto")  # Kallisto specific ARGs
cat(paste("\tDone!\n", sep = ""))

# Changes WD back to original
setwd(wd)

# Extracts Gene COUNTS & TPM
counts <- txi.kallisto$counts
tpm <- txi.kallisto$abundance

# Rounds Gene COUNTS (floats to integers)
## Not doing right now... I do not see the point in that
#counts <- round(x = counts, digits = 0)

###################
# 3. WRITES RESULTS
###################

# Writes results to files
## Change manually the OUT filenames
cat(paste("Writing OUT files...\n", sep = ""))
write.table(x = counts,
            file = paste(outDir, "/all_COUNTS_", myBatch, "-", myComp, ".txt", sep = ""),
            eol = "\n", na = c(NA,"NA",""), dec = ".",
            row.names = T, col.names = T, quote = F, sep = "\t")
write.table(x = tpm,
            file = paste(outDir, "/all_TPM_", myBatch, "-", myComp, ".txt", sep = ""),
            eol = "\n", na = c(NA,"NA",""), dec = ".",
            row.names = T, col.names = T, quote = F, sep = "\t")
cat(paste("\tDone!\n", sep = ""))
