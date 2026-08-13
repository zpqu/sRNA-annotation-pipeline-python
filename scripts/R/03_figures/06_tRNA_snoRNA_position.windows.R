##Date: 25/07/2015
##Author: Zhipeng
## This script tests whether piRNA-annotated reads are degradation products of
## tRNA/snoRNA genes (step 06). Reads come from the step-2 annotated objects
## (rdata/*.bam.annotated.gr.RData in the strategy output directory); only reads
## annotated as piRNA (type "piRNA", sense to a piRNA locus) are used.
## With the piRNA > snoRNA annotation priority (step 02), piRNA-annotated reads
## may overlap snoRNA genes on the same strand, so their position profiles over
## tRNA and snoRNA gene bodies reveal whether they look like degradation
## fragments of those genes.
## For each tRNA/snoRNA gene a 21 bp window is slid in 1 bp steps along the gene
## body, and the mean read count per window position is plotted per sample.
## The sense panel counts piRNA reads on the same strand as the gene, the
## antisense panel piRNA reads on the opposite strand (AS.piRNA reads are not
## used for the position profiles).
## Both read flavors are counted:
##   all reads    : weighted by 'count' (redundant reads)
##   unique reads : each non-redundant read counted once
## A summary figure reports the % of piRNA reads that overlap snoRNA
## (sense/antisense), tRNA (sense/antisense) genes or neither; only reads
## annotated as piRNA are included (AS.piRNA reads are excluded).
##
## Outputs (step 06):
##   figures/Figure_06.piRNA_vs_tRNA_pos_barplot(.pdf/.png) + _AS
##   figures/Figure_06.piRNA_vs_snoRNA_pos_barplot(.pdf/.png) + _AS
##   figures/Figure_06.piRNA_overlap_summary(.pdf/.png)
##   tables/Table_06_piRNA_overlap_summary.csv
##   rdata/piRNA_on_tRNA.20bp.dis.all.df.RData, piRNA_on_snoRNA.20bp.dis.all.df.RData

suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(data.table))
rm(list = ls())

source("../lib/init.R")
dir.tab = file.path(out.dir, "tables")
dir.rdata = file.path(out.dir, "rdata")
dir.fig = file.path(out.dir, "figures")
dir.create(dir.fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir.tab, recursive = TRUE, showWarnings = FALSE)

## ---- sliding 21 bp windows (1 bp step) along each gene body -----------------
sliding.windows = function(gr){
      win = gr[integer(0)]
      for(i in seq_along(gr)){
            g = gr[i]
            len = width(g)
            if(len <= 20) next
            if(as.character(strand(g)) == "+"){
                  starts = start(g):(end(g) - 20)
            }else{
                  starts = (end(g) - 20):start(g)
            }
            pos = seq_along(starts)
            win = c(win, GRanges(seqnames = seqnames(g),
                                 ranges = IRanges(start = starts, end = starts + 20),
                                 strand = strand(g),
                                 position = pos,
                                 gene_id = values(g)$gene_id))
      }
      win
}

## ---- sum the weighted count of reads overlapping each window -----------------
weighted.counts = function(wins, reads, cnt, ignore.strand = TRUE){
      out = integer(length(wins))
      if(length(reads) > 0){
            ol = findOverlaps(wins, reads, type = "any", ignore.strand = ignore.strand)
            cs = rowsum(cnt[subjectHits(ol)], queryHits(ol))
            out[as.integer(rownames(cs))] = cs[, 1]
      }
      out
}

## ---- format tRNAs -----------------------------------------------------------
load(file.path(db.dir, "tRNA.gr.RData"))
tRNA.20bp.gr = sliding.windows(tRNA.gr)
save(tRNA.20bp.gr, file = file.path(dir.rdata, "tRNA.20bp.gr.RData"))
print(paste("tRNA windows:", length(tRNA.20bp.gr)))

## ---- format snoRNAs (genes longer than 20 bp) --------------------------------
load(file.path(db.dir, "snoRNA.gr.RData"))
snoRNA.gr = snoRNA.gr[width(snoRNA.gr) > 20]
snoRNA.20bp.gr = sliding.windows(snoRNA.gr)
save(snoRNA.20bp.gr, file = file.path(dir.rdata, "snoRNA.20bp.gr.RData"))
print(paste("snoRNA windows:", length(snoRNA.20bp.gr)))

## ---- opposite-strand copies of the gene windows (for antisense counts) --------
tRNA.anti.20bp.gr = tRNA.20bp.gr
strand(tRNA.anti.20bp.gr) = ifelse(as.character(strand(tRNA.anti.20bp.gr)) == "+", "-", "+")
snoRNA.anti.20bp.gr = snoRNA.20bp.gr
strand(snoRNA.anti.20bp.gr) = ifelse(as.character(strand(snoRNA.anti.20bp.gr)) == "+", "-", "+")

## ---- piRNA reads per window, per sample, for both read flavors ----------------
piRNA_on_tRNA.dis.all.df = NULL
piRNA_on_snoRNA.dis.all.df = NULL
piRNA.overlap.tab = NULL
files = list.files(path = dir.rdata, pattern = ".bam.annotated.gr.RData$")
if(length(files) == 0) stop("no *.bam.annotated.gr.RData found in ", dir.rdata, " (run step 02 first)")
for(i in seq(along = files)){
      file.name = file.path(dir.rdata, files[i])
      sample.name = files[i]
      sample.name = gsub("\\.bam\\.annotated\\.gr\\.RData", "", sample.name)
      sample.name = gsub("\\.(bwa|bowtie2)$", "", sample.name)
      print(paste("Now is processing ...", sample.name, Sys.time()))
      load(file.name)
      cnt = mcols(reads.bam.annotated.gr)$count
      pirna_sense.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "piRNA"]
      pirna.cnt = cnt[reads.bam.annotated.gr$type == "piRNA"]

      ## piRNA reads (sense to a piRNA locus) on tRNA windows; sense = same
      ## strand as the gene, antisense = opposite strand
      if(length(tRNA.20bp.gr) > 0){
            tRNA.base = as.data.frame(tRNA.20bp.gr)
            tRNA.all = tRNA.base
            tRNA.all$sense = weighted.counts(tRNA.20bp.gr, pirna_sense.gr, pirna.cnt,
                                                       ignore.strand = FALSE)
            tRNA.all$antisense = weighted.counts(tRNA.anti.20bp.gr, pirna_sense.gr, pirna.cnt,
                                                           ignore.strand = FALSE)
            tRNA.all$sample = sample.name
            tRNA.all$flavor = "all reads"
            tRNA.uniq = tRNA.base
            tRNA.uniq$sense = weighted.counts(tRNA.20bp.gr, pirna_sense.gr,
                                                        rep(1, length(pirna_sense.gr)),
                                                        ignore.strand = FALSE)
            tRNA.uniq$antisense = weighted.counts(tRNA.anti.20bp.gr, pirna_sense.gr,
                                                            rep(1, length(pirna_sense.gr)),
                                                            ignore.strand = FALSE)
            tRNA.uniq$sample = sample.name
            tRNA.uniq$flavor = "unique reads"
            piRNA_on_tRNA.dis.all.df = rbind(piRNA_on_tRNA.dis.all.df, tRNA.all, tRNA.uniq)
      }

      ## piRNA reads (sense to a piRNA locus) on snoRNA windows; sense = same
      ## strand as the gene, antisense = opposite strand
      if(length(snoRNA.20bp.gr) > 0){
            snoRNA.base = as.data.frame(snoRNA.20bp.gr)
            snoRNA.all = snoRNA.base
            snoRNA.all$sense = weighted.counts(snoRNA.20bp.gr, pirna_sense.gr, pirna.cnt,
                                                        ignore.strand = FALSE)
            snoRNA.all$antisense = weighted.counts(snoRNA.anti.20bp.gr, pirna_sense.gr, pirna.cnt,
                                                            ignore.strand = FALSE)
            snoRNA.all$sample = sample.name
            snoRNA.all$flavor = "all reads"
            snoRNA.uniq = snoRNA.base
            snoRNA.uniq$sense = weighted.counts(snoRNA.20bp.gr, pirna_sense.gr,
                                                        rep(1, length(pirna_sense.gr)),
                                                        ignore.strand = FALSE)
            snoRNA.uniq$antisense = weighted.counts(snoRNA.anti.20bp.gr, pirna_sense.gr,
                                                            rep(1, length(pirna_sense.gr)),
                                                            ignore.strand = FALSE)
            snoRNA.uniq$sample = sample.name
            snoRNA.uniq$flavor = "unique reads"
            piRNA_on_snoRNA.dis.all.df = rbind(piRNA_on_snoRNA.dis.all.df, snoRNA.all, snoRNA.uniq)
      }

      ## ---- piRNA overlap summary (snoRNA/tRNA sense or antisense, or neither) ----
      cat = rep("none", length(pirna_sense.gr))
      if(length(pirna_sense.gr) > 0){
            sno.sense = unique(queryHits(findOverlaps(pirna_sense.gr, snoRNA.gr, ignore.strand = FALSE)))
            cat[sno.sense] = "snoRNA_sense"
            sno.anti = setdiff(unique(queryHits(findOverlaps(pirna_sense.gr, snoRNA.gr, ignore.strand = TRUE))), sno.sense)
            cat[sno.anti] = "snoRNA_antisense"
            remain = which(cat == "none")
            if(length(remain)){
                  trna.sense = unique(queryHits(findOverlaps(pirna_sense.gr[remain], tRNA.gr, ignore.strand = FALSE)))
                  cat[remain[trna.sense]] = "tRNA_sense"
            }
            remain = which(cat == "none")
            if(length(remain)){
                  trna.anti = unique(queryHits(findOverlaps(pirna_sense.gr[remain], tRNA.gr, ignore.strand = TRUE)))
                  cat[remain[trna.anti]] = "tRNA_antisense"
            }
      }
      ov = data.table(sample = sample.name, category = factor(cat, levels = c("snoRNA_sense", "snoRNA_antisense",
                        "tRNA_sense", "tRNA_antisense", "none")),
                      count = as.numeric(pirna_sense.gr$count))
      ov.sum = ov[, .(n_unique = as.numeric(.N), n_reads = base::sum(count)), by = .(sample, category)]
      ov.sum[, pct_unique := round(100 * n_unique / base::sum(n_unique), 2)]
      ov.sum[, pct_reads := round(100 * n_reads / base::sum(n_reads), 2)]
      piRNA.overlap.tab = rbind(piRNA.overlap.tab, ov.sum)
}

save(piRNA_on_tRNA.dis.all.df, file = file.path(dir.rdata, "piRNA_on_tRNA.20bp.dis.all.df.RData"))
save(piRNA_on_snoRNA.dis.all.df, file = file.path(dir.rdata, "piRNA_on_snoRNA.20bp.dis.all.df.RData"))
write.csv(piRNA.overlap.tab, file = file.path(dir.tab, "Table_06_piRNA_overlap_summary.csv"), row.names = FALSE)

## ---- plot position profiles (piRNA reads on gene windows) ---------------------
position.plot = function(dis.all.df, fig.base, fig.title){
      if(is.null(dis.all.df) || nrow(dis.all.df) == 0) return(NULL)
      dis.all.df$flavor = factor(dis.all.df$flavor, levels = c("unique reads", "all reads"))
      sample.num = length(unique(dis.all.df$sample))
      fd = fig.dims(sample.num, 2, per.h = 4.1)
      fig.width  = fd["width"]
      fig.height = fd["height"]
      test.sense.df = aggregate(sense ~ position + sample + flavor, data = dis.all.df, mean)
      test.antisense.df = aggregate(antisense ~ position + sample + flavor, data = dis.all.df, mean)

      p.sense = ggplot(data = test.sense.df, aes(x = position, y = sense)) +
            geom_bar(stat = "identity") +
            labs(title = paste0(fig.title, " - piRNA reads, sense position distribution (20 bp windows)"),
                 x = "", y = "Mean count") +
            scale_y_continuous(labels = comma) +
            theme_bw() + small.font() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            facet_grid(sample ~ flavor, scales = "free_y")
      ggsave(p.sense, file = paste0(fig.base, ".pdf"), width = fig.width, height = fig.height)
      ggsave(p.sense, file = paste0(fig.base, ".png"), width = fig.width, height = fig.height, dpi = 300)

      p.antisense = ggplot(data = test.antisense.df, aes(x = position, y = antisense)) +
            geom_bar(stat = "identity") +
            labs(title = paste0(fig.title, " - piRNA reads, antisense position distribution (20 bp windows)"),
                 x = "", y = "Mean count") +
            scale_y_continuous(labels = comma) +
            theme_bw() + small.font() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            facet_grid(sample ~ flavor, scales = "free_y")
      ggsave(p.antisense, file = paste0(fig.base, "_AS.pdf"), width = fig.width, height = fig.height)
      ggsave(p.antisense, file = paste0(fig.base, "_AS.png"), width = fig.width, height = fig.height, dpi = 300)
      print(paste("saved:", fig.base, Sys.time()))
}

position.plot(piRNA_on_tRNA.dis.all.df, file.path(dir.fig, "Figure_06.piRNA_vs_tRNA_pos_barplot"), "tRNA")
position.plot(piRNA_on_snoRNA.dis.all.df, file.path(dir.fig, "Figure_06.piRNA_vs_snoRNA_pos_barplot"), "snoRNA")

## ---- plot piRNA overlap summary (stacked % bar, samples x flavor) -------------
ov.long = melt(piRNA.overlap.tab[, .(sample, category, n_unique, n_reads)],
               id.vars = c("sample", "category"),
               measure.vars = c("n_unique", "n_reads"),
               variable.name = "flavor", value.name = "n",
               variable.factor = FALSE)
ov.long[flavor == "n_unique", flavor := "unique reads"]
ov.long[flavor == "n_reads", flavor := "all reads"]
ov.long[, flavor := factor(flavor, levels = c("unique reads", "all reads"))]
ov.long[, pct := 100 * n / base::sum(n), by = .(sample, flavor)]
ov.long[, category := factor(category, levels = c("snoRNA_sense", "snoRNA_antisense",
                                                   "tRNA_sense", "tRNA_antisense", "none"))]
cat.cols = c(snoRNA_sense = "#1f78b4", snoRNA_antisense = "#a6cee3",
             tRNA_sense = "#e31a1c", tRNA_antisense = "#fb9a99",
             none = "grey70")
sample.num = length(unique(ov.long$sample))
fd = fig.dims(sample.num, 2, per.h = 2.7)
p.summary = ggplot(ov.long, aes(x = sample, y = pct, fill = category)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = cat.cols, name = "piRNA overlap") +
      labs(title = "Percentage of piRNA reads overlapping snoRNA/tRNA genes",
           x = NULL, y = "% of piRNA reads") +
      theme_bw() + small.font() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_grid(. ~ flavor)
ggsave(p.summary, file = file.path(dir.fig, "Figure_06.piRNA_overlap_summary.pdf"),
       width = fd["width"], height = fd["height"])
ggsave(p.summary, file = file.path(dir.fig, "Figure_06.piRNA_overlap_summary.png"),
       width = fd["width"], height = fd["height"], dpi = 300)
print("piRNA overlap summary figure saved")
