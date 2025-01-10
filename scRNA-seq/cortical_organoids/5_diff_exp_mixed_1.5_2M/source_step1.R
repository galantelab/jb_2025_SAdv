#Filter differentially expressed protein-coding genes (|log2FC| > 0.1 & adj.pvalue < 0.05) ######################################

#List files
files = list.files(path=".", pattern="unfilt")

#For each file
for(file in files){

  #Read file
  df = read.delim(file, header=T)

  #Filter diff exp protein-coding genes with Bonferroni adjusted p-value < 0.05
  df = df[df$gene.type=="protein_coding" & df$pval.adj<0.05,]

  #Separate up/down-regulated genes
  up = df[df$avg.log2FC > 0,]
  down = df[df$avg.log2FC < 0,]

  #Reoder dataframes by cell population and pval.adj
  up = up[order(up$cell.population, up$pval.adj),]
  down = down[order(down$cell.population, down$pval.adj),]

  #Save result to file
  write.table(up, gsub("unfilt", "up", file), row.names=F, col.names=T, quote=F, sep="\t")
  write.table(down, gsub("unfilt", "down", file), row.names=F, col.names=T, quote=F, sep="\t")

}
