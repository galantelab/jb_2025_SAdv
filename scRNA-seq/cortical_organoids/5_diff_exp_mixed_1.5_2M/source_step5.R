#Make bubbleplots without REVIGO filtering ########################################################################

#Load libraries
library(ggplot2)
library(tidyverse)

#Read final results
df = read.delim("brain_enrichment_GO-BP.txt", header=T)

#Select relevant columns
df = df %>% select(comparison, direction, cell.population, GO.term, fold.enrichment, p.adjust, n.genes)

#Add group info: cell population + direction
df$group = paste(df$cell.population, df$direction, sep=" ")

#Reorder dataframe
df$direction = factor(df$direction, levels=c("up", "down"))
df$cell.population = factor(df$cell.population, levels=c("aRG", "Cortical hem", "Newborn CFuPN", "PN"))
df = df[order(df$direction, df$cell.population),]
df$group = factor(df$group, levels=df[!duplicated(df$group),]$group)

#For each comparison
for(comparison in df[!duplicated(df$comparison),]$comparison){

  #Filter current comparison
  tmp = df[df$comparison==comparison,]

  #Keep terms occurring in more than one group
  tmp$count.GO.term = ave(tmp$GO.term, tmp$GO.term, FUN = length)
  tmp$count.GO.term = as.numeric(tmp$count.GO.term)
  to_keep1 = tmp[tmp$count.GO.term > 1,]

  #For the remaining terms, filter top 10 most enriched terms per group (cell population + direction)
  to_check = tmp[tmp$count.GO.term == 1,]
  to_keep2 = data.frame(comparison=character(), direction=character(), cell.population=character(), GO.term=character(), fold.enrichment=numeric(), p.adjust=numeric(), group=character(), count.GO.term=numeric())
  for(group in to_check[!duplicated(to_check$group),]$group){

    #Filter current group
    tmp2 = to_check[to_check$group==group,]

    #Reorder dataframe by fold enrichment and add top 10 terms to final dataframe
    #to_keep2 = rbind(to_keep2, head(tmp2[order(tmp2$fold.enrichment, decreasing=T),],10))

    #Filter GO terms with n.genes >= 5 and sort by fold.enrichment
    filt1 = tmp2[tmp2$n.genes >= 5, ]
    filt1 = filt1[order(-filt1$fold.enrichment), ]

    #Filter GO terms with n.genes < 5 and sort by fold.enrichment
    filt2 = tmp2[tmp2$n.genes < 5, ]
    filt2 = filt2[order(-filt2$fold.enrichment), ]

    #Combine the two datasets, prioritizing terms with n.genes >= 5
    top_terms = rbind(filt1, filt2)

    #Add top 5 GO terms to final dataframe
    to_keep2 = rbind(to_keep2, head(top_terms, 5))

  }
  tmp = rbind(to_keep1, to_keep2)

  #Sort the dataframe by group_count (descending) and then by GO.term (optional, alphabetically)
  tmp = tmp[order(-tmp$count.GO.term, tmp$count.GO.term==1, tmp$group, -tmp$fold.enrichment), ]
  tmp$GO.term = factor(tmp$GO.term, levels=rev(tmp[!duplicated(tmp$GO.term),]$GO.term))

  #Create the bubble plot
  ggplot(tmp, aes(x=group, y=GO.term, size=fold.enrichment, color=p.adjust)) + geom_point() + theme_classic() + theme(axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10)) + xlab("") + ylab("") + scale_color_gradient(low="steelblue4", high="steelblue1", limits=c(0, 0.05)) + scale_size_continuous(range=c(1, 6), limits = c(2, 110)) + scale_x_discrete(drop = FALSE)
  if(comparison=="ar_10uM_vs_ar_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_blue.pdf"), width=8, height=4) }
  if(comparison=="ar_30uM_vs_ar_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_blue.pdf"), width=8, height=6) }
  if(comparison=="hu_10uM_vs_hu_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_blue.pdf"), width=8, height=8) }
  if(comparison=="hu_30uM_vs_hu_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_blue.pdf"), width=8, height=8.5) }
  if(comparison=="ar_CTL_vs_hu_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_blue.pdf"), width=8, height=6) }

  ggplot(tmp, aes(x=group, y=GO.term, size=fold.enrichment, color=p.adjust)) + geom_point() + theme_classic() + theme(axis.text.x=element_text(size=10, angle=45, hjust=1), axis.text.y=element_text(size=10)) + xlab("") + ylab("") + scale_color_gradient(low="firebrick4", high="firebrick1", limits=c(0, 0.05)) + scale_size_continuous(range=c(1, 6), limits = c(2, 110)) + scale_x_discrete(drop = FALSE)
  if(comparison=="ar_10uM_vs_ar_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_red.pdf"), width=8, height=4) }
  if(comparison=="ar_30uM_vs_ar_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_red.pdf"), width=8, height=6) }
  if(comparison=="hu_10uM_vs_hu_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_red.pdf"), width=8, height=8) }
  if(comparison=="hu_30uM_vs_hu_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_red.pdf"), width=8, height=8.5) }
  if(comparison=="ar_CTL_vs_hu_CTL"){ ggsave(paste0("bubbleplots_", comparison, "_red.pdf"), width=8, height=6) }

}

