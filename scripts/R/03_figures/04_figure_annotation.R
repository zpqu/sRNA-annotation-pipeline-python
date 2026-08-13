##Date: 21/07/2015
##Author: Zhipeng
## This script is used to make plots for small RNA annotation results (step 3).
## It reads the two consolidated count tables produced by step 02 from the
## strategy output directory resolved by the shared bootstrap (init.R):
##   Table_02a_annotation_count_unique_reads.csv   (unique reads)
##   Table_02b_annotation_count_all_reads.csv      (all reads)
## For each class (read, matmiRNA, snoRNA, piRNA, tRNA) a count barplot (log10
## y-axis) and a percentage barplot are drawn as separate figures, using a 2D
## facet grid (samples on the y-axis rows, read flavor unique/all on the columns).
## Colour is used only to differentiate sense, antisense (AS.*) and
## other/intergenic categories.
## Output figures follow the step naming rule: Figure_04a, Figure_04b, ...

library(ggplot2)
library(scales)
rm(list = ls())

source("../lib/init.R")
dir.tab = file.path(out.dir, "tables")
dir.fig = file.path(out.dir, "figures")
dir.create(dir.fig, recursive = TRUE, showWarnings = FALSE)

###read consolidated count tables from step 2
unique.tab = read.csv(file.path(dir.tab, "Table_02a_annotation_count_unique_reads.csv"))
all.tab = read.csv(file.path(dir.tab, "Table_02b_annotation_count_all_reads.csv"))
unique.tab$flavor = "unique reads"
all.tab$flavor = "all reads"
all.df = rbind(unique.tab, all.tab)

class.map = c(read = "read.annotation",
              matmiRNA = "matmiRNA.annotation",
              snoRNA = "snoRNA.annotation",
              piRNA = "piRNA.annotation",
              tRNA = "tRNA.annotation")
class.list = names(class.map)

read.levels = c("matmiRNA", "piRNA", "snoRNA", "tRNA", "RM", "refGene.NM.exon",
                "refGene.NM.intron", "lincRNA.exon", "AS.matmiRNA", "AS.piRNA",
                "AS.snoRNA", "AS.tRNA", "AS.RM", "AS.refGene.NM.exon",
                "AS.refGene.NM.intron", "AS.lincRNA.exon", "other")
region.levels = c("CDS", "5UTR", "3UTR", "intron", "up1k", "down1k", "RM",
                  "AS.CDS", "AS.5UTR", "AS.3UTR", "AS.intron", "AS.up1k",
                  "AS.down1k", "AS.RM", "intergenic")

###sense/antisense/other orientation
orient.of = function(lv){
      factor(ifelse(grepl("^AS\\.", lv), "antisense",
             ifelse(lv %in% c("intergenic", "other"), "other", "sense")),
             levels = c("sense", "antisense", "other"))
}
orient.cols = c(sense = "#1f78b4", antisense = "#e31a1c", other = "grey50")

rot.theme = theme(axis.text.x = element_text(angle = 45, hjust = 1))

for(j in 1:length(class.list)){
      class = class.list[j]
      cat.name = class.map[[class]]
      levels.use = if(class == "read") read.levels else region.levels
      sub.dt = all.df[all.df$category == cat.name, ]
      if(nrow(sub.dt) == 0) next
      sub.dt$item = factor(sub.dt$item, levels = levels.use)
      sub.dt$group = orient.of(sub.dt$item)
      sample.num = length(unique(sub.dt$sample))
      fd = fig.dims(sample.num, 2, per.h = 4.1)
      fig.width = fd["width"]
      fig.height = fd["height"]

      ##percentage panel: complete all item x sample x flavor combos (zero-filled)
      all.combo = expand.grid(sample = unique(sub.dt$sample),
                              item = levels(sub.dt$item),
                              flavor = c("unique reads", "all reads"),
                              stringsAsFactors = FALSE)
      sub.full = merge(sub.dt[, c("sample", "item", "flavor", "Freq")], all.combo,
                       by = c("sample", "item", "flavor"), all = TRUE)
      sub.full$Freq[is.na(sub.full$Freq)] = 0
      sub.full$item = factor(sub.full$item, levels = levels(sub.dt$item))
      sub.full$group = orient.of(sub.full$item)
      sub.full$per = ave(sub.full$Freq, sub.full$sample, sub.full$flavor,
                         FUN = function(x) x / sum(x))

      p.count = ggplot(sub.dt, aes(x = item, y = Freq, fill = group)) +
            geom_bar(stat = "identity", width = 0.7) +
            scale_y_log10(labels = comma) +
            scale_fill_manual(values = orient.cols, name = "Orientation") +
            facet_grid(sample ~ flavor) +
            labs(title = paste0(class, " annotation - count"),
                 y = "Count (log10 scale)", x = NULL) +
            theme_bw() + rot.theme + small.font()

      p.per = ggplot(sub.full, aes(x = item, y = per, fill = group)) +
            geom_bar(stat = "identity", width = 0.7) +
            scale_y_continuous(labels = percent) +
            scale_fill_manual(values = orient.cols, name = "Orientation") +
            facet_grid(sample ~ flavor) +
            labs(title = paste0(class, " annotation - percentage"),
                 y = "Percentage", x = NULL) +
            theme_bw() + rot.theme + small.font()

      save.plot = function(gg, base) {
            pdf(paste0(base, ".pdf"), width = fig.width, height = fig.height)
            print(gg)
            dev.off()
            png(paste0(base, ".png"), width = fig.width * 300, height = fig.height * 300, res = 300)
            print(gg)
            dev.off()
            print(paste("saved:", base, Sys.time()))
      }
      save.plot(p.count, file.path(dir.fig, paste0("Figure_04", letters[j], ".", class, "_annotation_count_barplot")))
      save.plot(p.per, file.path(dir.fig, paste0("Figure_04", letters[j], ".", class, "_annotation_percentage_barplot")))
}
