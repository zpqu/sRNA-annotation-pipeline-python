##Date: 16/07/2015
##Author: Zhipeng
## This script is used to format small RNA data as rdata
## Update (31/07/2026): also build a non-redundant reads GRanges where reads
## with exactly the same coordinates (chr, start, end, strand) appear only
## once, store their frequency as a "count" metadata column, plot the read
## size distribution against the count distribution/density, and export the
## underlying summaries as tables.
##
## Outputs (step 01, written under out.base = output/ or output/comparison/ in
## comparison mode; the shared pre-processing stage runs once per pipeline run):
##   <out.base>/rdata/<sample>.bam.gr.RData            all reads
##   <out.base>/rdata/<sample>.bam.unique.gr.RData     non-redundant reads with 'count'
##   <out.base>/figures/Figure_01.read_size_vs_count.pdf/.png   (faceted by sample)
##   <out.base>/tables/Table_01a_sample_summary.csv
##   <out.base>/tables/Table_01b_read_size_distribution.csv
##   <out.base>/tables/Table_01c_read_count_distribution.csv
##   <out.base>/tables/Table_01d_read_size_vs_count.csv

library(rtracklayer)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
rm(list = ls())

source("../lib/init.R")
dir.rdata = file.path(out.base, "rdata")
dir.fig = file.path(out.base, "figures")
dir.tab = file.path(out.base, "tables")
if(!dir.exists(dir.rdata)) dir.create(dir.rdata, recursive = TRUE)
if(!dir.exists(dir.fig)) dir.create(dir.fig, recursive = TRUE)
if(!dir.exists(dir.tab)) dir.create(dir.tab, recursive = TRUE)

files = list.files(path = "../../../bams/", pattern = ".bam$")
summ.list = list()
size.list = list()
count.list = list()
bin.list = list()
dens.list = list()

for(i in seq(along = files)){
  file.name = paste("../../../bams/", files[i], sep = "");
  sam = files[i];
  print(paste("Now is processing ...", sam, "at", Sys.time()))
  sam = gsub("\\.bam", "", sam)
  sam = gsub("\\_R1\\.trimmed\\_NEB\\.bwa\\.hg19", "", sam)
  sample.label = gsub("\\.(bwa|bowtie2)$", "", sam)
  reads.bam = readGAlignments(file.name, use.names = T)
  reads.strand = gsub("chr", "", seqnames(reads.bam))
  reads.strand = gsub("Chr", "", reads.strand)
  reads.strand = paste("chr", reads.strand, sep = "")
  reads.bam.gr = GRanges(seqnames = reads.strand,
			    ranges = ranges(reads.bam),
			    strand = strand(reads.bam))
  rm(reads.bam)
  gc()

  ## ----- non-redundant reads --------------------------------------------------
  ## reads with identical chr/start/end/strand are shown once; their frequency
  ## is recorded in the 'count' metadata column of the GRanges object.
  read.key = paste(seqnames(reads.bam.gr), start(reads.bam.gr),
                   end(reads.bam.gr), strand(reads.bam.gr), sep = ":")
  read.uniq = unique(read.key)
  read.count = tabulate(match(read.key, read.uniq), nbins = length(read.uniq))
  is.first = !duplicated(read.key)
  reads.bam.unique.gr = reads.bam.gr[is.first]
  mcols(reads.bam.unique.gr)$count = read.count
  rm(read.key, read.uniq, read.count, is.first)
  gc()

  print(paste("total reads:", length(reads.bam.gr),
              " unique reads:", length(reads.bam.unique.gr),
              " count range:", paste(range(mcols(reads.bam.unique.gr)$count), collapse = "-"),
              "at", Sys.time()))

  ## ----- save -----------------------------------------------------------------
  gr.name = file.path(dir.rdata, paste0(sam, ".bam.gr.RData"))
  save(reads.bam.gr, file = gr.name)
  gr.unique.name = file.path(dir.rdata, paste0(sam, ".bam.unique.gr.RData"))
  save(reads.bam.unique.gr, file = gr.unique.name)

  ## ----- plot/table data -------------------------------------------------------
  plot.dt = data.table(width = width(reads.bam.unique.gr),
                       count = mcols(reads.bam.unique.gr)$count)
  size.dt = plot.dt[, .(n_unique = .N, n_reads = sum(count)), by = width]
  cdist.dt = plot.dt[, .(n_unique = .N), by = .(lc = round(log2(count)))]
  bin.dt = plot.dt[, .(n_unique = .N), by = .(width, lc = round(log2(count)))]

  ## ----- sample summary ---------------------------------------------------------
  tb = sort(table(plot.dt$width), decreasing = TRUE)
  wb = plot.dt[, .(n_reads = sum(count)), by = width][order(-n_reads)]
  summ.list[[i]] = data.table(
    sample = sample.label,
    total_reads = length(reads.bam.gr),
    unique_reads = length(reads.bam.unique.gr),
    singleton_reads = sum(plot.dt$count == 1),
    pct_singletons = round(100 * sum(plot.dt$count == 1) / nrow(plot.dt), 2),
    median_count = median(plot.dt$count),
    max_count = max(plot.dt$count),
    dominant_size_by_unique_nt = as.integer(names(tb)[1]),
    dominant_size_by_unique_n = as.integer(tb[1]),
    dominant_size_by_reads_nt = wb$width[1],
    dominant_size_by_reads_n = wb$n_reads[1],
    top5_sizes_nt = paste(names(tb)[1:5], collapse = ";")
  )
  size.list[[i]] = copy(size.dt)[, sample := sample.label]
  count.list[[i]] = copy(cdist.dt)[, sample := sample.label]
  bin.list[[i]] = copy(bin.dt)[, sample := sample.label]

  ## ---- density data (subsampled, one entry per unique read) ----------------------
  set.seed(123)
  n.plot = nrow(plot.dt)
  dens.list[[i]] = data.table(sample = sample.label,
                              lc = log2(plot.dt$count[base::sample(n.plot, min(n.plot, 1e6))]))

  rm(reads.bam.gr, reads.bam.unique.gr, plot.dt, size.dt, cdist.dt, bin.dt)
  gc()
}

## ----- combined tables ------------------------------------------------------------
summ.tab = rbindlist(summ.list)
size.tab = rbindlist(size.list)
count.tab = rbindlist(count.list)
bin.tab = rbindlist(bin.list)
dens.tab = rbindlist(dens.list)

write.csv(summ.tab, file.path(dir.tab, "Table_01a_sample_summary.csv"), row.names = FALSE)
write.csv(size.tab, file.path(dir.tab, "Table_01b_read_size_distribution.csv"), row.names = FALSE)
write.csv(count.tab, file.path(dir.tab, "Table_01c_read_count_distribution.csv"), row.names = FALSE)
write.csv(bin.tab, file.path(dir.tab, "Table_01d_read_size_vs_count.csv"), row.names = FALSE)

## ----- size vs count distribution/density figures (samples as rows, 3 plot types as columns)
x.breaks = seq(floor(min(size.tab$width) / 5) * 5, max(size.tab$width), by = 5)
rot.theme = theme(axis.text.x = element_text(angle = 45, hjust = 1))
n.samp = length(unique(size.tab$sample))

p.size = ggplot(size.tab, aes(x = width, y = n_unique)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(breaks = x.breaks) +
  facet_wrap(~sample, ncol = 1, scales = "free_x") +
  labs(title = "Read size distribution (unique reads)",
       x = "Read size (nt)", y = "Number of unique reads") +
  theme_bw() + rot.theme + small.font()

p.count = ggplot(dens.tab, aes(x = lc)) +
  geom_density(fill = "grey70", alpha = 0.6) +
  facet_wrap(~sample, ncol = 1, scales = "free_y") +
  labs(title = "Count distribution (unique reads)",
       x = "log2(count)", y = "Density") +
  theme_bw() + rot.theme + small.font()

p.2d = ggplot(bin.tab, aes(x = width, y = lc, fill = n_unique)) +
  geom_tile() +
  scale_fill_gradient(low = "grey90", high = "darkred", trans = "log10",
                      name = "unique reads") +
  scale_x_continuous(breaks = x.breaks) +
  facet_wrap(~sample, ncol = 1, scales = "free") +
  labs(title = "Read size vs count density",
       x = "Read size (nt)", y = "log2(count)") +
  theme_bw() + rot.theme + small.font()

per.h.mm = 69
total.height.mm = min(n.samp * per.h.mm, 220)
fig.width  = 180 / 25.4
fig.height = total.height.mm / 25.4
fig.base = file.path(dir.fig, "Figure_01.read_size_vs_count")
pdf(paste0(fig.base, ".pdf"), width = fig.width, height = fig.height)
print(p.size | p.count | p.2d)
dev.off()
png(paste0(fig.base, ".png"), width = round(fig.width * 300), height = round(fig.height * 300), res = 300)
print(p.size | p.count | p.2d)
dev.off()

print(paste0("Tables saved to ", dir.tab, ": Table_01a_sample_summary.csv, Table_01b_read_size_distribution.csv, Table_01c_read_count_distribution.csv, Table_01d_read_size_vs_count.csv"))
