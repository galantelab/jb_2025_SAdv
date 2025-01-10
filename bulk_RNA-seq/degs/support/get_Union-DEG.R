# FFERREIRA 12/26/2024
# Sanford Consortium - UCSD
# Prepares table to export as SUP table

# Loads LIBs
library(tidyverse)

# Sets WD
wd <- getwd()
setwd(wd)

# OBJs to read INFILE
g1 <- "10"  # Treatment
g2 <- "CTRL"  # Control
c1 <- "C136"  # cell lineage #1
c2 <- "NOVA1-C15"  # cell lineage #2
batch <- "batch2"  # good batch used
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
degs1 <- as.data.frame(read.delim(file = paste(wd, "/result_DESeq2_", batch, "-",
                                               c1, "-", g1, "_vs_", c1, "-", g2, ext, sep = ""),
                                  header = T, row.names = 1, check.names = F, dec = dec, sep = sep, na.strings = nas))
degs2 <- as.data.frame(read.delim(file = paste(wd, "/result_DESeq2_", batch, "-",
                                               c2, "-", g1, "_vs_", c2, "-", g2, ext, sep = ""),
                                  header = T, row.names = 1, check.names = F, dec = dec, sep = sep, na.strings = nas))

# Adds annotation info
data1 <- merge(x = map, y = degs1, by.x = 1, by.y = 0, all = F, sort = T)
data2 <- merge(x = map, y = degs2, by.x = 1, by.y = 0, all = F, sort = T)

# Removes ROWs with any NA
data1 <- na.omit(data1)
data2 <- na.omit(data2)

# Renames COLs
colnames(data1) <- c("Gene.ID", "Gene.Symbol", "Base.Mean", "L2FC", "L2FCSE", "Statistic", "P-Value", "FDR")
colnames(data2) <- c("Gene.ID", "Gene.Symbol", "Base.Mean", "L2FC", "L2FCSE", "Statistic", "P-Value", "FDR")

# Selects only DEGs
data1 <- data1[(abs(data1[["L2FC"]]) > l2fc & data1[["FDR"]] < fdr),]
data2 <- data2[(abs(data2[["L2FC"]]) > l2fc & data2[["FDR"]] < fdr),]

# Checks redundancy
sort(intersect(x = data1$Gene.Symbol, y = data2$Gene.Symbol))

# Combines the data frames
combined_data <- bind_rows(data1, data2)

# Selects rows with the highest absolute L2FC for each gene
data <- combined_data %>% group_by(Gene.Symbol) %>% slice_max(order_by = abs(L2FC), n = 1) %>% ungroup()

# Changes FDR values for better visualization in Volcano Plot
data[["1-FDR"]] <- 1 - data[["FDR"]]

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
## ALL
write.table(x = data, file = paste(wd, "/final_DEG_table_NOVA1-ArAr-", g1, "_vs_", g2, ext, sep = ""),
            sep = sep, na = nas, dec = dec, quote = F, row.names = F, col.names = T)

## UP
write.table(x = data[data[["DEG"]] == "UP",],
            file = paste(wd, "/final_DEG-UP_table_NOVA1-ArAr-", g1, "_vs_", g2, ext, sep = ""),
            sep = sep, na = nas, dec = dec, quote = F, row.names = F, col.names = T)

## DOWN
write.table(x = data[data[["DEG"]] == "DOWN",],
            file = paste(wd, "/final_DEG-DOWN_table_NOVA1-ArAr-", g1, "_vs_", g2, ext, sep = ""),
            sep = sep, na = nas, dec = dec, quote = F, row.names = F, col.names = T)
