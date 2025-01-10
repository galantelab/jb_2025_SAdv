#Run GO-BP enrichment analyses ###################################################################

#Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

#List files
files = list.files(path=".", pattern="exclusive_")

#Create empty dataframe to store results
final = data.frame(GO.id=character(), GO.term=character(), gene.ratio=character(), bg.ratio=character(), rich.factor=numeric(), fold.enrichment=numeric(), z.score=numeric(), p.value=numeric(), p.adjust=numeric(), q.value=numeric(), list.genes=character(), n.genes=numeric(), direction=character(), comparison=character(), cell.population=character(), direction=character(), comparison=character())

#For each file
for(file in files){

  #Extract direction (up/down) and comparison
  direction = unlist(strsplit(file, "_"))[2]
  comparison = gsub(".txt", "", paste(unlist(strsplit(file, "_"))[3:7], collapse="_"))

  #Read file
  df = read.delim(file, header=T)

  #Run enrichment for each cell population
  for(cell_pop in c("aRG", "Cortical hem", "Newborn CFuPN", "PN")){

    cat(comparison, direction, cell_pop, "\n")

    #Get gene symbols and convert them to Entrez ids
    gene_symbols = df[df$cell.population==cell_pop,]$gene.symbol
    gene_entrez = bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
    gene_entrez = gene_entrez$ENTREZID

    #Run enrichment analysis for GO-BP, with default p-value and adjusted p-value (FDR) thresholds
    go_bp_enrich = enrichGO(gene=gene_entrez, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.2, minGSSize=5, maxGSSize=2000, pAdjustMethod="BH", readable=T)
    res = as.data.frame(go_bp_enrich)
    rownames(res) = NULL

    #If there is any enrichment
    if(nrow(res) >0 ){

      #Add cell type, direction (up/down) and comparison
      res$cell.population = cell_pop
      res$direction = direction
      res$comparison = comparison

      #Rename columns
      colnames(res) = c("GO.id", "GO.term", "gene.ratio", "bg.ratio", "rich.factor", "fold.enrichment", "z.score", "p.value", "p.adjust", "q.value", "list.genes", "n.genes", "cell.population", "direction", "comparison")

      #Format genes list
      res$list.genes = gsub("/", " ", res$list.genes)

      #Add result to final dataframe
      final = rbind(final, res)

    }

  }

}

#Reorder dataframe columns
final = final %>% select("comparison", "direction", "cell.population", "GO.id", "GO.term", "gene.ratio", "bg.ratio", "rich.factor", "fold.enrichment", "z.score", "p.value", "p.adjust", "q.value", "n.genes", "list.genes")

#Reorder dataframe rows
final = final[order(final$comparison, final$direction, final$cell.population, final$p.adjust),]

#Write result to file
write.table(final, "enrichment_GO-BP.txt", row.names=F, col.names=T, quote=F, sep="\t")
