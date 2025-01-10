#Run differential gene expression analyses between cell populations ####################################################

#Load libraries
library(Seurat)
library(tidyverse)

#Create directory to store plots
dir.create("5_diff_exp_mixed_1.5_2M", showWarnings=T)

#Load mapped Seurat data with reductions
df = LoadSeuratRds("3_integrated_samples/integrated_mixed_1.5_2M_with_reductions.rds")

#Create column with final grouped samples in metadata
metadata = df@meta.data
metadata = metadata %>% rownames_to_column(var="cell.id")
groups = read.delim("sample_groups.txt", header=T)
metadata = merge(metadata, groups, by.x="orig.ident", by.y="sample") %>% column_to_rownames(var="cell.id") %>% rename("sample"="group")
df@meta.data = metadata

#Set order of grouped samples
df@meta.data$sample = factor(df@meta.data$sample, levels=c("NOVA1 hu/hu CTL", "NOVA1 hu/hu 10uM", "NOVA1 hu/hu 30uM", "NOVA1 ar/ar CTL", "NOVA1 ar/ar 10uM", "NOVA1 ar/ar 30uM"))

#Read map between gene symbols and types
map_genes = read.table("map_genes.txt", header=F)
colnames(map_genes) = c("gene.symbol", "gene.type")
map_genes = map_genes[!duplicated(map_genes$gene.symbol),]

#Create function to run differential gene expression analysis (group1 x group2)
run_diff <- function(data, group1, group2, comparison){

  #Set identity class to make comparisons
  Idents(data) = "sample"

  #Create empty dataframe to store results
  final = data.frame(gene.symbol=character(), gene.type=character(), avg.log2FC=numeric(), pct.1=numeric(), pct.2=numeric(), pct.1_pct.2=numeric(), pval.adj=numeric(), cell.population=character())

  #For each cell population
  for(cell_population in c("aRG", "Newborn CFuPN", "Cortical hem", "PN")){

    #Subset cell population
    curr_df = subset(data, predicted.id == cell_population)

    #Run FindMarkers
    diff = suppressWarnings(FindMarkers(object=curr_df, assay="integrated", ident.1=group1, ident.2=group2, slot="data", logfc.threshold=0.1, min.pct=0.01))

    #Format data frame
    diff = diff %>% select(-"p_val") %>% rename("avg.log2FC"="avg_log2FC", "pval.adj"="p_val_adj")

    #Add difference between pct 1 and pct2, and reorder columns
    diff$pct.1_pct.2 = diff$pct.1 - diff$pct.2
    diff = diff %>% select("avg.log2FC", "pct.1", "pct.2", "pct.1_pct.2", "pval.adj")

    #Add column with gene type
    diff = merge(map_genes, diff, by.x="gene.symbol", by.y=0)

    #Add column with cell population
    diff$cell.population = cell_population

    #Add result to final dataframe
    final = rbind(final, diff)    

  }

  #Write result to file
  write.table(final, paste0("5_diff_exp_mixed_1.5_2M/unfilt_", comparison, ".txt"), row.names=F, col.names=T, quote=F, sep="\t")

}

#Call function to run differential gene expression analyses between groups ------------------------

  #NOVA1 ar/ar CTL x NOVA1 hu/hu CTL
  run_diff(df, "NOVA1 ar/ar CTL", "NOVA1 hu/hu CTL", "ar_CTL_vs_hu_CTL")

  #NOVA1 hu/hu 10uM x NOVA1 hu/hu CTL
  run_diff(df, "NOVA1 hu/hu 10uM", "NOVA1 hu/hu CTL", "hu_10uM_vs_hu_CTL")

  #NOVA1 hu/hu 30uM x NOVA1 hu/hu CTL
  run_diff(df, "NOVA1 hu/hu 30uM", "NOVA1 hu/hu CTL", "hu_30uM_vs_hu_CTL")

  #NOVA1 ar/ar 10uM x NOVA1 ar/ar CTL
  run_diff(df, "NOVA1 ar/ar 10uM", "NOVA1 ar/ar CTL", "ar_10uM_vs_ar_CTL")

  #NOVA1 ar/ar 30uM x NOVA1 ar/ar CTL
  run_diff(df, "NOVA1 ar/ar 30uM", "NOVA1 ar/ar CTL", "ar_30uM_vs_ar_CTL")
