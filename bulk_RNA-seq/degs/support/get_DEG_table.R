# FFERREIRA 12/26/2024
# Sanford Consortium - UCSD
# Prepares table to export as SUP table

# Loads LIBs
library(tidyverse)

# Sets WD
wd <- getwd()
setwd(wd)

# OBJs to read INFILE
g1 <- "NOVA1-ArAr-CTRL"  # LEFT
g2 <- "NOVA1-HuHu-CTRL"  # RIGHT
ext <- ".tsv"
sep <- "\t"
dec <- "."
nas <- c(NA,"NA","")
l2fc <- 1
fdr <- .05

# Reads MAP file
map <- as.data.frame(read.delim(file = paste(wd, "/map_gene_IDs", ext, sep = ""), header = F,
                                check.names = F, dec = dec, sep = sep, na.strings = nas))

# Renames MAP COLs
colnames(map) <- c("ID", "Name")

# Reads expression raw results after DESeq2
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/09_sleuth.html (output explained)
degs <- as.data.frame(read.delim(file = paste(wd, "/result_DESeq2_batch2-", g1, "_vs_", g2, ext, sep = ""),
                                 header = T, row.names = 1, check.names = F, dec = dec, sep = sep, na.strings = nas))

# Adds annotation info
data <- merge(x = map, y = degs, by.x = 1, by.y = 0, all = F, sort = T)

# Removes ROWs with any NA
data <- na.omit(data)

# Renames COLs
colnames(data) <- c("Gene.ID", "Gene.Symbol", "Base.Mean", "L2FC", "L2FCSE", "Statistic", "P-Value", "FDR")

# Changes FDR values for better visualization in Volcano Plot
data[["1-FDR"]] <- 1 - data$FDR

# Categorizes DEGs in UP and DOWN
data <- mutate(data, DEG = ifelse(test = (data[["L2FC"]] > l2fc & data[["FDR"]] < fdr), yes = "UP",
                                  no = ifelse(test = (data[["L2FC"]] < (-l2fc) & data[["FDR"]] < fdr),
                                              yes = "DOWN", no = "NS")))

# Order DEG COL
data[["DEG"]] <- factor(x = data[["DEG"]], levels = c("UP", "DOWN", "NS"))

# Sorts data frame
data <- data[with(data, order(data[["DEG"]], data[["L2FC"]], data[["FDR"]], data[["Base.Mean"]], data[["Gene.Symbol"]],
                              decreasing = F, na.last = T, method = "radix")),]

# Writes final SUP table
write.table(x = data, file = paste(wd, "/final_DEG_table_", g1, "_vs_", g2, ext, sep = ""),
            sep = sep, na = nas, dec = dec, quote = F, row.names = F, col.names = T)
