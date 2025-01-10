# FFERREIRA 05/26/2024
# Sanford Consortium - UCSD
# DEG results from bulk "full" normal brain organoids - Volcano Plot

################
# 0. SETS UP ENV
################

# Loads LIBs
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(tidyverse)

# Figure CONFIG
dpi <- 1000
#formats <- c("jpeg","pdf","png","svg","tiff")
formats <- c("pdf", "png")
options(scipen = 999) # DISABLE SCIENTIFIC NOTATION

# Sets WD
wd <- getwd()
setwd(wd)

# Sets theme to be used in plots
## theme_bw: The classic dark-on-light ggplot2 theme. May work better for presentations displayed with a projector
### WARN: Adds the extra X/Y lines that I do not like!
## theme_classic: A classic-looking theme, with x and y axis lines and no gridlines
### WARN: Removes the gridlines that are sometimes useful!
mytheme <- #theme_bw(base_family = "ArialMT", base_size = 12) +
  theme_classic(base_family = "ArialMT", base_size = 12) +
  theme(panel.grid.major.y = element_line(size = .05, colour = "#2F3437"),
        panel.grid.minor.y = element_line(size = .025, colour = "#2F3437"),
        text = element_text(family = "ArialMT", face = "plain", size = 10, colour = "#2F3437"),
        axis.title = element_text(family = "ArialMT", size = 14,
                                  vjust = 1, hjust = .5, face = "bold", colour = "#090A0B"),
        axis.text = element_text(family = "ArialMT", size = 12,
                                 face = "bold", hjust = 0.5, vjust = .5, colour = "#2F3437"),
        plot.title = element_text(family = "ArialMT", size = 18,
                                  vjust = 2, hjust = 0.5, face = "bold", lineheight = 1, colour = "#090A0B"),
        plot.subtitle = element_text(family = "ArialMT", size = 16,
                                     vjust = 2, hjust = 0.5, face = "italic", lineheight = 1, colour = "#090A0B"),
        legend.text = element_text(family = "ArialMT",
                                   colour = "#2F3437", size = 10, face = "plain"),
        legend.title = element_text(family = "ArialMT",
                                    colour = "#2F3437", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.spacing.x = unit(-.1, "cm"),
        axis.ticks = element_line(size = .5, colour = "#2F3437"),
        axis.line = element_line(size = .5, colour = "#2F3437"))

##################
# 1. PREPARES DATA
##################

# OBJs to read INFILE
## GROUPS: ASC-10, ASC-30, ASC-CTRL, C136-10, C136-30, C136-CTRL,
## GROUPS: NOVA1-C15-10, NOVA1-C15-30, NOVA1-C15-CTRL, NOVA1-C27-10, NOVA1-C27-30, NOVA1-C27-CTRL
g1 <- "NOVA1-C27-10"  # LEFT
g2 <- "NOVA1-C27-CTRL"  # RIGHT
ext <- ".tsv"
sep <- "\t"
dec <- "."
nas <- c(NA,"NA","")
l2fc <- 1
fdr <- .05

# Reads MAP file
map <- as.data.frame(read.delim(file = paste(wd, "/map_gene_IDs", ext, sep = ""), header = F,
                                check.names = F, dec = dec, sep = sep, na.strings = nas))

# Removes and renames MAP COLs
# map <- map[,c(1,2)]
colnames(map) <- c("ID", "Name")

# Reads expression raw results after DESeq2
## PS: after filtering for surfaceome genes only

# https://hbctraining.github.io/DGE_workshop_salmon/lessons/09_sleuth.html (output explained)
data <- as.data.frame(read.delim(file = paste(wd, "/Inputs/DEGs_", g1, "_vs_", g2, ".ALL.final", ext, sep = ""),
                                 header = T, check.names = F, dec = dec, sep = sep, na.strings = nas))

# Keeps only relevant COLs
# data <- data[,c(1,2,6)]

# Adds annotation info
data <- merge(x = data, y = map, by = 1, all = F, sort = T)

# Removes ROWs with any NA
data <- na.omit(data)

# Renames COLs
colnames(data) <- c("Gene.ID", "AVG.L2.COUNTS", "L2FC", "FDR", "Gene.Symbol")

# Iterates over selection options to get at most 50 IDs to plot as labels
# IDs <- rep("test",51)
# selection <- 1
# while (length(IDs) > 50)
# {
#   IDs <- data[((abs(data[["L2FC"]]) > (l2fc * selection)) & data[["FDR"]] < fdr),][["Gene.ID"]]
#   print(paste("Selection: ", selection, sep = ""))
#   print(paste("'IDs' length: ", length(IDs), sep = ""))
#   selection <- selection + 1
# }

# Creates threshold variable
data <- mutate(data, threshold = ifelse(test = (data[["L2FC"]] > l2fc & data[["FDR"]] < fdr), yes = "UP",
                                        no = ifelse(test = (data[["L2FC"]] < (-l2fc) & data[["FDR"]] < fdr),
                                                    yes = "DOWN", no = "NS")))

# Names specific genes that would not be named due to threshold
data$Label <- F
# data$Label <- ""
# data$Label <- ifelse(test = data[["Gene.ID"]] %in% IDs, yes = T, no = F)

# Additional layer of filtering for labeling
## PS: maybe this does not make sense anymore
# data$Label <- ifelse(data[["Label"]] & abs(data[["L2FC"]]) > (l2fc * selection), yes = T, no = F)

# Changes FDR values for better visualization in Volcano Plot
data$FDRPlot <- 1 - data$FDR

############################
# 2 - GENERATES VOLCANO PLOT
############################

# Makes Volcano Plot
## -log10(FDR)
VolcanoPlot <- ggplot(data = data, aes(x = L2FC, y = FDRPlot)) +
  geom_point(aes(col = threshold), size = 1.5, alpha = .3) +
  scale_color_manual(values = c("UP" = "#1465AC", "DOWN" = "#B31B21", "NS" = "darkgray"),
                     name = "DEGs (NOVA1-C27 10uM vs. Control)",
                     breaks = c("DOWN", "UP"), #breaks = c("DOWN", "UP", "NS"),
                     # labels = c(paste("Down: ", table(data$threshold)["DOWN"][[1]], sep = ""),
                     #            paste("Up: ", table(data$threshold)["UP"][[1]], sep = ""),
                     #            paste("NS: ", table(data$threshold)["NS"][[1]], sep = "")),
                     labels = c(paste("Down: ", table(data$threshold)["DOWN"][[1]], sep = ""),
                                paste("Up: ", table(data$threshold)["UP"][[1]], sep = ""))) +
  xlab("Log2 Fold Change (L2FC)") + ylab("Adjusted P-Value (1-FDR)") +
  geom_hline(yintercept = (1 - fdr), linetype = "dashed", col = "#1C211E") +
  geom_vline(xintercept = c((-l2fc), l2fc), linetype = "dashed", col = "#1C211E") +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = .1)) + mytheme +
  scale_x_continuous(breaks = seq(from = round(min(data[["L2FC"]]),0), to = round(max(data[["L2FC"]]),0), by = 1)) +
  geom_label_repel(aes(label = ifelse(data$Label, as.character(data[["Gene.Symbol"]]),"")), size = 2,
                   box.padding = unit(2.0, "lines"), label.padding = unit(0.2, "lines"), fill = "#FFFFFF",
                   color = "#090A0B", segment.size  = .1, segment.color = "#1C211E", max.overlaps = 100)

# Prints FIG
print(x = VolcanoPlot, newpage = T)

# ELSEVIER DPI = 1000
for (format in formats)
{
  ggsave(filename = paste("VolcanoPlot_", g1, "_vs_", g2, ".DEG.", toupper(format), sep = ""),
         plot = VolcanoPlot, device = format, path = paste(wd, "Plots", "PAPER", sep = "/"),
         scale = 1, units = "in", limitsize = T, dpi = dpi, width = 7, height = 5)
}

# Saves UP & DOWN DEGs for further Enrichement Analyses
## Numbers here are different from the pipeline because here we have all genes
## There we have only the protein-coding ones

### UP
write.table(x = data[data$threshold == "UP", c(1,5,3,4,2)],
            file = paste(wd, "/Plots/PAPER/input_", g1, "_vs_", g2, ".UP", ext, sep = ""),
            sep = sep, na = nas, dec = dec, quote = F, row.names = F, col.names = T)

### DOWN
write.table(x = data[data$threshold == "DOWN", c(1,5,3,4,2)],
            file = paste(wd, "/Plots/PAPER/input_", g1, "_vs_", g2, ".DOWN", ext, sep = ""),
            sep = sep, na = nas, dec = dec, quote = F, row.names = F, col.names = T)
