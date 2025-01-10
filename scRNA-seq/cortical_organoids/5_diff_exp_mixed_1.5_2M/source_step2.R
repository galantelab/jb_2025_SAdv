#Plot Venn diagrams with intersections between cell populations ################################################################

#Load library
library(VennDiagram)

#List files
up_files = list.files(path=".", pattern="^up_")
down_files = list.files(path=".", pattern="^down_")
files = c(up_files, down_files)

#For each file
for(file in files){

  #Read file
  df = read.delim(file, header=T)

  #Split genes by cell population
  list_genes = split(df$gene.symbol, df$cell.population)

  #Make Venn diagram
  pdf(paste0("venn_", gsub(".txt", ".pdf", file)), width=5.5, height=5.5)
  venn.plot = venn.diagram(x=list_genes, filename=NULL, output=T, 
    category.names=gsub(" ", "\n", names(list_genes)),
    fill=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6"), 
    alpha=0.3, cat.cex=1.5, cex=1.2, cat.dist=c(0.2, 0.23, 0.13, 0.08),
    lty="blank", show.plot=T)
  grid.draw(venn.plot)
  dev.off()

  #Filter only exclusive genes per cell population and write to file
  exclusive_genes = df[!duplicated(df$gene.symbol) & !duplicated(df$gene.symbol, fromLast=T), ]
  write.table(exclusive_genes, paste0("exclusive_", file), row.names=F, col.names=T, quote=F, sep="\t")

}
