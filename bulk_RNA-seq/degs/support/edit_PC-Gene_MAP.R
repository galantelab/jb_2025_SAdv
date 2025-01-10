# FFERREIRA 11/27/2023
# Livia Project - NASA

# Edits MAP files so that we have only unique gene appearences

# Sets WS
wd <- getwd()
setwd(wd)

# General OBJs
sep <- "\t"
nas <- c("NA", "", NA)

# Reads MAP file between Ensembl Gene IDs and Symbols / Types
## PS: First COL must be set as "row.names"
map <- as.data.frame(read.delim(file = paste(wd, "/../map_gene_IDs.txt", sep = ""),
                                header = F, sep = sep, na.strings = nas))

# Sorts table
map <- map[with(map, order(map$V1, map$V3,
                           decreasing = F, na.last = T, method = "radix")),]

# Removes DUPs
## I had already done this before
# map <- map[!duplicated(map),]

# Selects only protein-coding (PC) genes
map.pc <- map[map$V3 == "protein_coding",]

# Writes edited MAP
write.table(x = map.pc, file = "map_pcgene_IDs.txt", row.names = F, col.names = F, quote = F, sep = sep, na = nas)
