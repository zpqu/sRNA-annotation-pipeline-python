##Date: 2026
##Author: Zhipeng
## Step 03: abundance structure of the current strategy's annotation (out.dir).
## For each genome-feature category we compute the per-locus read abundance and
## describe how reads are distributed across loci:
##   matmiRNA : per mature miRNA locus (miRBase Name, up to 2110 loci)
##   snoRNA   : per snoRNA gene (refGene gene_name, 135 loci)
##   tRNA     : per tRNA gene (<genome>_tRNAs gene_id, 1137 loci for mm39)
##   piRNA    : per read position (chr:start:end:strand) -- the piRNA feature set
##              is ~83 M predicted features (mm39), so per-read loci are used instead
## Contents:
##   1) Table_03a: per-category summary statistics (n loci, reads, median/mean log10,
##                 Gini, top-1/5/10 % read shares, top locus, cross-sample tests)
##   2) Table_03b: long-format per-locus abundance (category, locus, reads) per sample
##   3) Figure_03a: per-locus log10 read-distribution (violin + boxplot)
##   4) Figure_03b: rank-abundance (Whittaker) curves per category
##   5) Figure_03c: Lorenz curves (matmiRNA) -- skewness visualisation
## Outputs are written to <out.dir>/tables and <out.dir>/figures (init.R).
## All three figures use facet_wrap(~sample) so each sample is shown in its own panel.

library(ggplot2)
library(scales)
library(data.table)
library(GenomicRanges)

out = function(...) cat(sprintf(...), "\n")

source("../lib/init.R")
dir.out = out.dir
dir.create(file.path(dir.out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dir.out, "figures"), recursive = TRUE, showWarnings = FALSE)

## ---- per-locus read abundance for one category --------------------------------------
locus.tab = data.table()
for(s in samples){
  load(paste0(dir.out, "/rdata/", s, ".bam.annotated.gr.RData"))
  for(cat in c("matmiRNA", "snoRNA", "tRNA")){
    g = reads.bam.annotated.gr[reads.bam.annotated.gr$type == cat]
    if(length(g) == 0) next
    load(file.path(db.dir, paste0(cat, ".gr.RData")))
    feat = get(paste0(cat, ".gr"))
    name.col = switch(cat, matmiRNA = "Name", snoRNA = "gene_name", tRNA = "gene_id")
    ol = findOverlaps(g, feat, type = "any")
    if(length(ol) == 0) next
    dt = data.table(sample = s, category = cat,
                    locus = mcols(feat)[[name.col]][subjectHits(ol)],
                    count = as.numeric(g$count[queryHits(ol)]))
    dt = dt[, .(n_reads = base::sum(count)), by = .(sample, category, locus)]
    locus.tab = rbind(locus.tab, dt)
  }
  ## piRNA: per read position
  g = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "piRNA"]
  if(length(g) > 0){
    dt = data.table(sample = s, category = "piRNA",
                    locus = paste(seqnames(g), start(g), end(g), strand(g)),
                    n_reads = as.numeric(g$count))
    locus.tab = rbind(locus.tab, dt)
  }
}

## ---- Gini coefficient (on per-locus read counts, sorted ascending) -------------------
gini = function(x){
  x = as.numeric(x)
  if(length(x) < 2) return(NA_real_)
  x = sort(x)
  n = length(x)
  if(base::sum(x) == 0) return(NA_real_)
  (2*base::sum(seq_len(n)*x))/(n*base::sum(x)) - (n + 1)/n
}

## ---- Table_03a: summary statistics + cross-sample tests ---------------------------------
sum.tab = data.table()
for(cat in c("matmiRNA", "snoRNA", "tRNA", "piRNA")){
  sub = locus.tab[category == cat]
  ## cross-sample statistics on loci present in both samples (paired Spearman /
  ## Wilcoxon are only defined for exactly two samples; otherwise left NA)
  rho = NA_real_
  w = NULL
  if (length(samples) == 2){
    cs = dcast(sub, locus ~ sample, value.var = "n_reads", fill = 0)
    if (all(samples %in% names(cs))){
      setnames(cs, samples, c("v1", "v2"))
      cs[, v1 := log10(v1 + 1)][, v2 := log10(v2 + 1)]
      rho = suppressWarnings(cor(cs$v1, cs$v2, method = "spearman"))
      w = suppressWarnings(wilcox.test(cs$v1, cs$v2, paired = TRUE))
    }
  }
  for(s in samples){
    x = sub[sample == s]
    if(nrow(x) == 0) next
    tot = base::sum(x$n_reads)
    setorder(x, -n_reads)
    sum.tab = rbind(sum.tab, data.table(
      category = cat, sample = s,
      n_loci = nrow(x),
      total_reads = tot,
      mean_log10 = round(mean(log10(x$n_reads + 1)), 3),
      median_log10 = round(median(log10(x$n_reads + 1)), 3),
      gini = round(gini(x$n_reads), 3),
      top1_pct = round(100*x$n_reads[1]/tot, 2),
      top5_pct = round(100*base::sum(x$n_reads[1:min(5, nrow(x))])/tot, 2),
      top10_pct = round(100*base::sum(x$n_reads[1:min(10, nrow(x))])/tot, 2),
      top_locus = x$locus[1],
      spearman_cross_sample = if (length(samples) == 2) round(rho, 3) else NA_real_,
      wilcoxon_p_2sample = if (length(samples) == 2) format.pval(w$p.value, digits = 3) else NA_character_))
  }
}
write.csv(sum.tab, paste0(dir.out, "/tables/Table_03a_category_abundance_summary.csv"), row.names = FALSE)
out("== Table_03a: per-category abundance summary ==")
print(sum.tab, row.names = FALSE)

## ---- Table_03b: long-format per-locus abundance ------------------------------------------
write.csv(locus.tab[order(category, sample, -n_reads)],
          paste0(dir.out, "/tables/Table_03b_per_locus_abundance.csv"), row.names = FALSE)

## ---- Kruskal-Wallis across categories (per sample; descriptive, units differ) -----------
for(s in samples){
  sub = locus.tab[sample == s][, log10(n_reads + 1), by = category]
  kw = kruskal.test(V1 ~ category, data = sub)
  out("%s Kruskal-Wallis across categories (matmiRNA/snoRNA/tRNA loci + piRNA read-positions): p = %s",
      s, format.pval(kw$p.value, digits = 3))
}

## ---- Figure 1: per-locus log10 read distribution ----------------------------------------
plot.dt = locus.tab[, log10_reads := log10(n_reads + 1)]
plot.dt[, category := factor(category, levels = c("matmiRNA", "piRNA", "tRNA", "snoRNA"))]
cat.cols = c(matmiRNA = "#1f78b4", piRNA = "#33a02c", tRNA = "#e31a1c", snoRNA = "#ff7f00")
f.ncol = if (length(samples) <= 8) 2 else 4
p1 = ggplot(plot.dt, aes(x = category, y = log10_reads, fill = category)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.2) +
  facet_wrap(~sample, ncol = f.ncol) +
  labs(title = "per-locus read abundance by annotation category (current strategy)",
       x = NULL, y = expression(log[10](reads + 1))) +
  scale_fill_manual(values = cat.cols) +
  theme_bw() + small.font()
d1 = fig.dims(length(samples), f.ncol, per.h = 4.1)
ggsave(p1, file = paste0(dir.out, "/figures/Figure_03a_per_locus_distribution.pdf"),
       width = d1["width"], height = d1["height"])
ggsave(p1, file = paste0(dir.out, "/figures/Figure_03a_per_locus_distribution.png"),
       width = d1["width"], height = d1["height"], dpi = 300)

## ---- Figure 2: rank-abundance (Whittaker) curves ------------------------------------------
ra = locus.tab[, .(n_reads = base::sum(n_reads)), by = .(sample, category, locus)]
ra[, rank := frank(-n_reads, ties.method = "first"), by = .(sample, category)]
ra[, n_loci := .N, by = .(sample, category)]
ra[, rank := rank + rnorm(.N, 0, 0.08)]  ## jitter ties
ra[, category := factor(category, levels = c("matmiRNA", "piRNA", "tRNA", "snoRNA"))]
p2 = ggplot(ra, aes(x = rank, y = n_reads, color = category)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  facet_wrap(~sample, scales = "free_x", ncol = f.ncol) +
  labs(title = "rank-abundance curves per annotation category",
       x = "rank of locus (log10)", y = "reads (log10)") +
  scale_color_manual(values = cat.cols) +
  theme_bw() + small.font()
d2 = fig.dims(length(samples), f.ncol, per.h = 4.1)
ggsave(p2, file = paste0(dir.out, "/figures/Figure_03b_rank_abundance.pdf"),
       width = d2["width"], height = d2["height"])
ggsave(p2, file = paste0(dir.out, "/figures/Figure_03b_rank_abundance.png"),
       width = d2["width"], height = d2["height"], dpi = 300)

## ---- Figure 3: Lorenz curves for matmiRNA ------------------------------------------------
lorenz = function(x){
  x = sort(as.numeric(x))
  n = length(x)
  c(0, cumsum(x)/base::sum(x), 1)[1:(n + 1)]
}
lz = data.table()
for(s in samples){
  x = locus.tab[sample == s & category == "matmiRNA"]$n_reads
  x = sort(x, decreasing = TRUE)
  n = length(x)
  lz = rbind(lz, data.table(sample = s, i = c(0, seq_len(n))/n,
                            cum_share = c(0, cumsum(x)/base::sum(x))))
}
p3 = ggplot(lz, aes(i, cum_share)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(color = "steelblue", linewidth = 0.6) +
  coord_equal() +
  facet_wrap(~sample, ncol = f.ncol) +
  labs(title = "Lorenz curves - mature miRNA loci (the more bowed, the more skewed)",
       x = "cumulative fraction of loci", y = "cumulative fraction of reads") +
  theme_bw() + small.font()
d3 = fig.dims(length(samples), f.ncol, per.h = 4.1)
ggsave(p3, file = paste0(dir.out, "/figures/Figure_03c_lorenz_matmiRNA.pdf"),
       width = d3["width"], height = d3["height"])
ggsave(p3, file = paste0(dir.out, "/figures/Figure_03c_lorenz_matmiRNA.png"),
       width = d3["width"], height = d3["height"], dpi = 300)

## ---- console summary ----------------------------------------------------------------------
out("== per-category dominance (all reads, both samples) ==")
for(cat in c("matmiRNA", "piRNA", "tRNA", "snoRNA")){
  for(s in samples){
    r = sum.tab[category == cat & sample == s]
    if(nrow(r) == 0) next
    out("  %-9s %-15s n_loci=%-7.0f total=%-11.0f Gini=%-5.2f top1=%-5.1f%% top5=%-5.1f%% top10=%-5.1f%% top=%-22s",
        cat, s, r$n_loci, r$total_reads, r$gini, r$top1_pct, r$top5_pct, r$top10_pct, r$top_locus)
  }
}
out("Abundance analysis written to %s/tables/ and %s/figures/", dir.out, dir.out)
