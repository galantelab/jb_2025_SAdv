# FFERREIRA 12/11/2024
# BEPE - Nina Colaboration / Neanderthals project

# MAKE BUBBLE PLOTS FOR GENE ONTOLOGY BIOLOGICAL PROCESS (GO:BP)
## ShinyGO + REVIGO + chatGPT

# COMPARISON: NOVA1-ArAr-CTRL vs NOVA1-HuHu-CTRL

# Loads LIBs
library("RColorBrewer")
library("ggplot2")
library("ggtext")
library("tidyverse")
library("ggnewscale")
library("forcats")
library("patchwork")
library("Cairo")
library("extrafont")

# Imports fonts into R (I may need to do this once)
font_import()
loadfonts(device = "pdf")  

# DISABLE SCIENTIFIC NOTATION
options(scipen = 999)

# FIGURE CONFIG
dpi <- 1000
formats <- c("png", "pdf")
#formats <- c("jpeg","pdf","png","svg","tiff")

# Sets the WD
wd <- getwd()
setwd(wd)

# Sets theme (plots)
mytheme <- theme_bw(base_family = "ArialMT", base_size = 12) +
  theme_classic(base_family = "ArialMT", base_size = 12) +
  theme(panel.grid = element_blank(),
        text = element_text(family = "ArialMT", face = "plain", size = 12, colour = "black"),
        axis.title = element_text(family = "ArialMT", size = 24,
                                  vjust = 1, hjust = .5, face = "bold", colour = "black"),
        axis.text = element_text(family = "ArialMT", size = 14,
                                 face = "plain", hjust = .5, vjust = .5, colour = "black"),
        plot.title = element_text(family = "ArialMT", size = 30,
                                  vjust = 2, hjust = 0.5, face = "bold", lineheight = 1, colour = "black"),
        plot.subtitle = element_text(family = "ArialMT", size = 28,
                                     vjust = 2, hjust = 0.5, face = "italic", lineheight = 1, colour = "black"),
        legend.text = element_text(family = "ArialMT", colour = "black", size = 10, face = "plain"),
        legend.title = element_text(family = "ArialMT", colour = "black", size = 12, face = "plain"),
        axis.ticks = element_line(size = .25, colour = "black"),
        axis.line = element_line(size = .25, colour = "black"))

# Reads inputs
up <- as.data.frame(read.delim(file = paste(wd, "/final_NOVA1-ArAr-30uM_vs_CTRL.UP.tsv", sep = ""),
                               header = T, sep = "\t", dec = ".", na.strings = c(NA,"NA","")))
down <- as.data.frame(read.delim(file = paste(wd, "/final_NOVA1-ArAr-30uM_vs_CTRL.DOWN.tsv", sep = ""),
                                 header = T, sep = "\t", dec = ".", na.strings = c(NA,"NA","")))

# Edits inputs
up[["Category"]] <- "UP"
down[["Category"]] <- "DOWN"

# Merges inputs
tmp_data <- rbind(up, down)

# Reorders and renames COLs
tmp_data <- tmp_data[,c(3,2,1,5,4)]
colnames(tmp_data) <- c("Pathway", "Fold", "FDR", "Category", "Genes")

# Removes GO IDs
tmp_data[["Pathway"]] <- gsub(pattern = "^GO:\\d+\\s+", replacement = "", x = tmp_data[["Pathway"]])

# Removes extra spaces at the end of pathway names
tmp_data[["Pathway"]] <- trimws(x = tmp_data[["Pathway"]], which = "right")

# Reads chatGPT info
upGPT <- as.data.frame(read.delim(file = paste(wd, "/chatGPT_NOVA1-ArAr-30uM_vs_CTRL.UP.tsv", sep = ""),
                                  header = T, sep = "\t", dec = ".", na.strings = c(NA,"NA","")))
downGPT <- as.data.frame(read.delim(file = paste(wd, "/chatGPT_NOVA1-ArAr-30uM_vs_CTRL.DOWN.tsv", sep = ""),
                                    header = T, sep = "\t", dec = ".", na.strings = c(NA,"NA","")))

# Edits inputs
upGPT[["Category"]] <- "UP"
downGPT[["Category"]] <- "DOWN"

# Merges inputs
tmp_data2 <- rbind(upGPT, downGPT)

# Reorders and renames COLs
tmp_data2 <- tmp_data2[,c(2,3,4)]
colnames(tmp_data2) <- c("Pathway", "Group", "Category")

# Merges TMP tables
final_data <- merge(x = tmp_data2, y = tmp_data, by = c("Pathway", "Category"), all = F, sort = T)

# Establishes Y-axis colors based on GO Terms
final_data[["myYColors"]] <- NA
final_data[["myYColors"]] <- ifelse(test = final_data[["Group"]] == "Neural Tissue and Organ Development",
                                    yes = "#A49E28", no = final_data[["myYColors"]])
final_data[["myYColors"]] <- ifelse(test = final_data[["Group"]] == "Neurodevelopment and Neuronal Differentiation",
                                    yes = "#C82F04", no = final_data[["myYColors"]])
final_data[["myYColors"]] <- ifelse(test = final_data[["Group"]] == "Synaptic and Neural Communication",
                                    yes = "#854488", no = final_data[["myYColors"]])
final_data[["myYColors"]] <- ifelse(test = final_data[["Group"]] == "Neurological Responses and Disorders",
                                    yes = "#1688B6", no = final_data[["myYColors"]])
final_data[["myYColors"]] <- ifelse(test = final_data[["Group"]] == "Neurodevelopment and Glial Activation",
                                    yes = "#0C7C59", no = final_data[["myYColors"]])
final_data[["myYColors"]] <- ifelse(test = final_data[["Group"]] == "Cell Adhesion, Motility, and Morphogenesis",
                                    yes = "#FF7733", no = final_data[["myYColors"]])

# Sorts table
final_data <- final_data[with(final_data, order(final_data[["Group"]], tolower(final_data[["Pathway"]]),
                                                decreasing = F, na.last = T, method = "radix")),]

# Adjusts levels
final_data[["Pathway"]] <- factor(final_data[["Pathway"]], levels = unique(final_data[["Pathway"]]))

# Writes Bubble Plot table
write.table(x = final_data, file = paste(wd, "bubblePlot_table.tsv", sep = "/"),
            quote = F, sep = "\t", na = c(NA,"NA",""), dec = ".", col.names = T, row.names = F)

# Ensures the factor rule stands for the subsets as well
down_data <- filter(final_data, Category == "DOWN") %>%
  mutate(Pathway = factor(Pathway, levels = levels(final_data$Pathway)))
up_data <- filter(final_data, Category == "UP") %>%
  mutate(Pathway = factor(Pathway, levels = levels(final_data$Pathway)))

# Makes plot
p <- ggplot(final_data, aes(fill = Group)) +
  geom_point(data = down_data,
             aes(x = Category, y = Pathway, size = Fold, color = FDR),
             alpha = .8, shape = 19, position = position_dodge(width = 0)) +
  scale_colour_gradientn(name = "FDR (DOWN)",
                         breaks = c(0.001, 0.01, 0.02),
                         labels = c("≤ 0.001", "≤ 0.01", "≤ 0.02"),
                         colors = c("#CC0007", "#FF999C")) +
  new_scale_colour() +
  geom_point(data = up_data,
             aes(x = Category, y = Pathway, size = Fold, color = FDR),
             alpha = .8, shape = 19, position = position_dodge(width = 0)) +
  scale_colour_gradientn(name = "FDR (UP)",
                         breaks = c(0.01, 0.025, 0.04),
                         labels = c("≤ 0.01", "≤ 0.025", "≤ 0.04"),
                         colors = c("#476785", "#BDCDDB")) +
  scale_y_discrete(drop = F) + scale_x_discrete(labels = c("Down", "Up")) +
  scale_size(range = c(1, 10), limits = c(1, 21), breaks = seq(2, 20, by = 4),
             labels = seq(2, 20, by = 4), name = "Fold Enrichment") +
  scale_fill_manual(name = "GO:BP Groups", breaks = final_data$Group,
                    values = final_data$myYColors, labels = final_data$Group) + mytheme +
  theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(colour = final_data[!duplicated(final_data[,c(1,6)]),][["myYColors"]]),
        legend.spacing = unit(1, "cm"), legend.position = "right", legend.direction = "vertical",
        legend.box = "vertical", aspect.ratio = 5, legend.margin = margin(0, 0, 0, 1, unit = "cm"),
        legend.box.margin = margin(0, 0, 0, 1, unit = "cm"),
        text = element_text(family = "ArialMT")) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 6), order = 1),
         size = guide_legend(label.position = "right", title.position = "top", title.hjust = 0, order = 2),
         color = guide_colorbar(title.position = "top", title.hjust = 0, ncol = 2, order = 3))

# Prints Boxplot
print(p)

# Saves FIG
for (format in formats)
{
  if (format == "pdf")
  {
    ggsave(filename = paste("BubblePlot_GO-BP_NOVA1-ArAr-30uM_vs_CTRL.v1.", format, sep = ""),
           plot = p, device = CairoPDF, path = wd, scale = 1, units = "in", limitsize = T,
           dpi = dpi, width = 12, height = 8)
  } else {
    ggsave(filename = paste("BubblePlot_GO-BP_NOVA1-ArAr-30uM_vs_CTRL.v1.", format, sep = ""),
           plot = p, device = format, path = wd, scale = 1, units = "in", limitsize = T,
           dpi = dpi, width = 12, height = 8)
  }
}
