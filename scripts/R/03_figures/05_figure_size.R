##Date: 21/07/2015
##Author: Zhipeng
## This script is used to make plot for size distribution for all small RNA data.
## Reads the consolidated count tables from the strategy output directory
## resolved by the shared bootstrap (init.R) and plots, for each annotation
## class, the read-size distribution faceted as samples (rows) x read flavor
## "unique reads"/"all reads" (columns), in the same style as Figure_04.
##
## Input tables (step 02):
##   tables/Table_02a_annotation_count_unique_reads.csv   (unique reads)
##   tables/Table_02b_annotation_count_all_reads.csv      (all reads)
##
## Outputs (step 05):
##   figures/Figure_05a.<class>_size_barplot.pdf/.png              (counts)
##   figures/Figure_05b.<class>_size_barplot.percentage.pdf/.png   (percentage)

library(ggplot2)
library(scales)
library(data.table)
rm(list = ls())

source("../lib/init.R")
dir.tab = file.path(out.dir, "tables")
dir.fig = file.path(out.dir, "figures")
dir.create(dir.fig, recursive = TRUE, showWarnings = FALSE)

###consolidated count tables from step 02 (unique + all reads)
t2a.file = file.path(dir.tab, "Table_02a_annotation_count_unique_reads.csv")
t2b.file = file.path(dir.tab, "Table_02b_annotation_count_all_reads.csv")
if(!file.exists(t2a.file) || !file.exists(t2b.file))
      stop("no Table_02a/Table_02b annotation count tables found in ", dir.tab, " (run step 02 first)")
unique.df = as.data.table(read.csv(t2a.file))
all.df = as.data.table(read.csv(t2b.file))
unique.df[, flavor := "unique reads"]
all.df[, flavor := "all reads"]
all.read.df = rbind(unique.df, all.df)
all.read.df[, flavor := factor(flavor, levels = c("unique reads", "all reads"))]

class.list = c("read", "matmiRNA", "piRNA", "snoRNA", "tRNA")
size.levels = sort(unique(as.integer(all.read.df$item[all.read.df$category %in%
                        paste0(class.list, ".size")])))

## x-axis tick labels every 5 nt on a continuous axis to avoid overstacking
x.breaks = seq(ceiling(min(size.levels) / 5) * 5, max(size.levels), by = 5)
if (x.breaks[1] > min(size.levels)) x.breaks = c(min(size.levels), x.breaks)

sample.num = length(unique(all.read.df$sample))
fig.d = fig.dims(sample.num, 2, per.h = 4.1)
fig.width  = fig.d["width"]
fig.height = fig.d["height"]
for(j in 1:length(class.list)){
      single.class.df = all.read.df[all.read.df$category == paste0(class.list[j], ".size"), ]
      if(nrow(single.class.df) == 0) next
      single.class.df[, per := Freq/sum(Freq), by = .(sample, flavor)]
      single.class.df$item = as.numeric(as.character(single.class.df$item))
      p.sample.class.barplot = ggplot(data = single.class.df, aes(x = item, y = Freq)) +
		   geom_col(width = 0.8) +
		   labs(title = paste0(class.list[j], " read-size distribution (counts)"),
		        x = "Read size (nt)", y = "Count") +
		   scale_y_continuous(labels = comma) +
		   scale_x_continuous(breaks = x.breaks) +
		   theme_bw() + small.font() +
		   facet_grid(sample ~ flavor, scales = "free_y")
      ggsave(p.sample.class.barplot, file = file.path(dir.fig, paste0("Figure_05a.", class.list[j], "_size_barplot.pdf")), width = fig.width, height = fig.height)
      ggsave(p.sample.class.barplot, file = file.path(dir.fig, paste0("Figure_05a.", class.list[j], "_size_barplot.png")), width = fig.width, height = fig.height, dpi = 300)

      p.sample.class.percent.plot = ggplot(data = single.class.df, aes(x = item, y = per)) +
		   geom_col(width = 0.8) +
		   labs(title = paste0(class.list[j], " read-size distribution (percentage)"),
		        x = "Read size (nt)", y = "Percentage") +
		   scale_y_continuous(labels = percent) +
		   scale_x_continuous(breaks = x.breaks) +
		   theme_bw() + small.font() +
		   facet_grid(sample ~ flavor)
      ggsave(p.sample.class.percent.plot, file = file.path(dir.fig, paste0("Figure_05b.", class.list[j], "_size_barplot.percentage.pdf")), width = fig.width, height = fig.height)
      ggsave(p.sample.class.percent.plot, file = file.path(dir.fig, paste0("Figure_05b.", class.list[j], "_size_barplot.percentage.png")), width = fig.width, height = fig.height, dpi = 300)
}
