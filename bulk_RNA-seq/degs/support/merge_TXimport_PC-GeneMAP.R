# FFERREIRA 03/24/2024
# Nina Project - Neanderthal

# Merges the MAP file for PC-Genes only created in the previous R code with TXimport data

# Sets WD
wd <- getwd()
setwd(wd)

# General OBJs
sep <- "\t"
nas <- c("NA", "", NA)
myComps <- c("ASC-10_vs_ASC-30", "ASC-10_vs_ASC-CTRL", "ASC-30_vs_ASC-CTRL",
             "C136-10_vs_C136-30", "C136-10_vs_C136-CTRL", "C136-30_vs_C136-CTRL",
             "NOVA1-C15-10_vs_NOVA1-C15-30", "NOVA1-C15-10_vs_NOVA1-C15-CTRL", "NOVA1-C15-30_vs_NOVA1-C15-CTRL",
             "NOVA1-C27-10_vs_NOVA1-C27-30", "NOVA1-C27-10_vs_NOVA1-C27-CTRL", "NOVA1-C27-30_vs_NOVA1-C27-CTRL")
myBatches <- c("batch1-2", "batch1", "batch2")

# Reads MAP file between Ensembl Gene IDs and Symbols / Types
map <- as.data.frame(read.delim(file = "../map_pcgene_IDs.txt", header = F, sep = sep, na.strings = nas))

for (myBatch in myBatches)
{
  for (myComp in myComps)
  {
    # Reads TXimport data
    ## PS: First COL must be set as "row.names"
    counts <- as.data.frame(read.delim(file = paste("all_COUNTS_", myBatch, "-", myComp, ".txt", sep = ""),
                                       row.names = 1, header = T, sep = sep, na.strings = nas))
    tpm <- as.data.frame(read.delim(file = paste("all_TPM_", myBatch, "-", myComp, ".txt", sep = ""),
                                    row.names = 1, header = T, sep = sep, na.strings = nas))
    
    # Reads conversion table
    conversion <- as.data.frame(read.delim(file = paste("mapID2GROUP.", myBatch, ".tsv", sep = ""),
                                           header = T, sep = sep, na.strings = nas))
    
    # Merges MAP file with TXimport files
    merge.counts <- merge(x = map, y = counts, by.x = "V1", by.y = 0, all = T)
    merge.tpm <- merge(x = map, y = tpm, by.x = "V1", by.y = 0, all = T)
    
    # Renames COLs
    names(merge.counts) <- c("ID", "Symbol", "Type", names(counts)[1:ncol(counts)])
    names(merge.tpm) <- c("ID", "Symbol", "Type", names(tpm)[1:ncol(tpm)])
    
    # Fills in the GAPs in 'Symbol' and 'Type' COLs for non-pc genes
    merge.counts$Symbol <- ifelse(test = is.na(merge.counts$Symbol), yes = merge.counts$ID, no = merge.counts$Symbol)
    merge.counts$Type <- ifelse(test = is.na(merge.counts$Type), yes = "non_coding", no = merge.counts$Type)
    merge.tpm$Symbol <- ifelse(test = is.na(merge.tpm$Symbol), yes = merge.tpm$ID, no = merge.tpm$Symbol)
    merge.tpm$Type <- ifelse(test = is.na(merge.tpm$Type), yes = "non_coding", no = merge.tpm$Type)
    
    # Use conversion table to get more significant group names for comparison
    myIDs <- names(merge.counts)[4:ncol(merge.counts)]
    conversion <- conversion[conversion$ID %in% myIDs,]
    myNewIDs <- c(rep(names(table(conversion$Group)[1]), table(conversion$Group)[[1]]),
                  rep(names(table(conversion$Group)[2]), table(conversion$Group)[[2]]))
    
    # Renames COLs again
    names(merge.counts) <- c("ID", "Symbol", "Type", myNewIDs)
    names(merge.tpm) <- c("ID", "Symbol", "Type", myNewIDs)
    
    # Sorts tables
    merge.counts <- merge.counts[with(merge.counts, order(merge.counts$Type, merge.counts$Symbol,
                                                          decreasing = T, na.last = T, method = "radix")),]
    merge.tpm <- merge.tpm[with(merge.tpm, order(merge.tpm$Type, merge.tpm$Symbol,
                                                 decreasing = T, na.last = T, method = "radix")),]
    
    # Removes DUPs
    ## I had already done this before
    # map <- map[!duplicated(map),]
    
    # Writes final IN tables for Degust (good for 1 x 1 comparisons)
    write.table(x = merge.counts, file = paste("counts_", myBatch, "-", myComp, ".annotated.sorted.tsv", sep = ""),
                row.names = F, col.names = T, quote = F, sep = sep, na = nas)
    write.table(x = merge.tpm, file = paste("tpm_", myBatch, "-", myComp, ".annotated.sorted.tsv", sep = ""),
                row.names = F, col.names = T, quote = F, sep = sep, na = nas)
  }
}
