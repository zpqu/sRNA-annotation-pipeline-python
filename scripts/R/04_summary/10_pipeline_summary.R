##Date: 2026
##Author: Zhipeng
## Step 10: end-of-pipeline summary. After all annotation/analysis steps have
## finished, this script
##   1) lists every output file (with a description and size), grouped by
##      output directory / producing step,
##   2) summarises the main results (sample stats, read sizes, annotation
##      composition, top mature miRNAs, per-category abundance),
##   3) runs a set of sanity checks (totals consistent, tables consistent,
##      logs error-free) and reports PASS / FAIL / WARN per check.
## Layout-aware:
##   single-strategy run   -> summarizes out.base = output/
##   comparison run        -> summarizes out.base = output/comparison/ including
##                            the per-strategy subfolders (fully_contained/union/any).
##
## Run from scripts/R/04_summary/ (same convention as the other pipeline
## scripts) or via run_smallRNA_annotation.sh.
##
## Outputs (step 10):
##   <out.base>/pipeline_summary.md        this report
##   <out.base>/tables/Table_10_sanity_checks.csv

suppressMessages({library(data.table); library(GenomicRanges)})

source("../lib/init.R")

## -----------------------------------------------------------------------------
## small helpers
## -----------------------------------------------------------------------------
out = function(...) cat(sprintf(...), "\n")
rep.lines = character()
rep = function(...) rep.lines <<- c(rep.lines, sprintf(...))

read.tab = function(f){
  if (file.exists(f)) fread(f) else NULL
}
fmt = function(x) format(x, big.mark = ",", scientific = FALSE)
pct2 = function(x) round(100 * x, 2)

## active strategies: in comparison mode all strategy dirs that exist; otherwise
## the current (single) strategy resolved by config.
strat.names = if (is.comparison) names(strategy.dir) else cur.strategy
strat.path = function(st){
  if (is.comparison) file.path(out.base, strategy.dir[[st]]) else out.dir
}
strat.dir.exists = function(st) dir.exists(strat.path(st))
active.strats = strat.names[sapply(strat.names, strat.dir.exists)]
st.label = function(st) if (is.comparison) st else ""

rep("## Small-RNA annotation pipeline - end-of-pipeline summary\n")
rep("- Reference genome : `%s`", genome)
rep("- Feature DB       : `%s` (built by scripts/R/00_build_DB/)", db.dir)
rep("- Samples          : %s", paste(samples, collapse = ", "))
rep("- Strategy mode    : `%s`%s", strategy.mode,
    ifelse(is.comparison, "", paste0(" (", cur.strategy, ")")))
rep("- Output base      : `%s`", out.base)
rep("- Strategies summarized: %s", paste(sprintf("`%s`", active.strats), collapse = ", "))
rep("- Step 03 (abundance) in active strategies: %s",
    paste(vapply(active.strats, function(st)
      sprintf("`%s` %s", st, ifelse(file.exists(file.path(strat.path(st), "tables/Table_03a_category_abundance_summary.csv")), "run", "not run")),
      character(1)), collapse = "; "))
rep("- Generated        : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
rep("")

## -----------------------------------------------------------------------------
## 1. Output-file inventory
## -----------------------------------------------------------------------------
rep("## 1. Output files\n")
rep("Each output file produced by the pipeline is listed below with a short description and its size on disk.\n")

SAM <- paste(sprintf("(%s)", samples), collapse = "|")
desc.rules <- list(
  ## step 01 tables
  list(re = "^Table_01a_sample_summary\\.csv$", d = "Per-sample summary: total reads, unique reads, singletons, dominant read sizes"),
  list(re = "^Table_01b_read_size_distribution\\.csv$", d = "Read size distribution (n_unique and n_reads per read size)"),
  list(re = "^Table_01c_read_count_distribution\\.csv$", d = "Read count (expression) distribution, log2 bins"),
  list(re = "^Table_01d_read_size_vs_count\\.csv$", d = "2D read-size x read-count grid"),
  ## step 02 consolidated count tables
  list(re = "^Table_02a_annotation_count_unique_reads\\.csv$", d = "Annotation composition, unique reads (each unique read counted once)"),
  list(re = "^Table_02b_annotation_count_all_reads\\.csv$", d = "Annotation composition, all reads (weighted by read frequency)"),
  ## step 02 per-unique-read annotation tables
  list(re = sprintf("^Table_02_(%s)_unique_reads_annotation\\.csv$", SAM), d = "Per-unique-read annotation + abundance for this sample (sorted by count descending)"),
  ## step s01 (overlap-rule comparison, comparison mode only)
  list(re = "^Table_s01a_overlap_rule_composition\\.csv$", d = "Overlap-rule (any/within/union) annotation-composition comparison"),
  list(re = "^Table_s01b_overlap_rule_category_totals\\.csv$", d = "Overlap-rule per-category read totals"),
  list(re = "^Table_s01c_read_movement_contained_vs_any\\.csv$", d = "Read movement between fully-contained and any strategies"),
  list(re = "^Table_s01d_read_movement_contained_vs_union\\.csv$", d = "Read movement between fully-contained and union strategies"),
  list(re = "^Table_s01e_mature_miRNA_expression_strategies\\.csv$", d = "Top mature miRNA expression under each overlap-rule strategy"),
  ## step 03 (abundance)
  list(re = "^Table_03a_category_abundance_summary\\.csv$", d = "Per-category abundance statistics (n loci, total reads, Gini, top-1/5/10 %, top locus, cross-sample tests)"),
  list(re = "^Table_03b_per_locus_abundance\\.csv$", d = "Per-locus read abundance, long format"),
  ## step 10 (this script)
  list(re = "^Table_10_sanity_checks\\.csv$", d = "Pipeline sanity-check results (step 10)"),
  list(re = "^Table_06_piRNA_overlap_summary\\.csv$", d = "piRNA reads overlapping snoRNA/tRNA genes (sense/antisense/none), unique + all reads (step 06)"),
  ## rdata
  list(re = sprintf("^(%s)\\.bam\\.gr\\.RData$", SAM), d = "All aligned reads (GRanges)"),
  list(re = sprintf("^(%s)\\.bam\\.unique\\.gr\\.RData$", SAM), d = "Non-redundant unique reads with 'count' (expression) column"),
  list(re = sprintf("^(%s)\\.bam\\.annotated\\.gr\\.RData$", SAM), d = "Annotated reads (category 'type', gene context 'region', 'count', 'feature.id', 'n_features')"),
  list(re = "^(snoRNA|tRNA)\\.20bp\\.gr\\.RData$", d = "20 bp sliding windows along each gene body (step 06)"),
  list(re = "^piRNA_on_(snoRNA|tRNA)\\.20bp\\.dis\\.all\\.df\\.RData$", d = "piRNA read counts per window position, unique + all reads (step 06)"),
  ## figures
  list(re = "^Figure_01\\.read_size_vs_count\\.(pdf|png)$", d = "Read size / count distribution figure, faceted by sample (step 01)"),
  list(re = "^Figure_04[a-e]\\..*\\.(pdf|png)$", d = "Annotation count / percentage barplots (step 04)"),
  list(re = "^Figure_05a\\..*_size_barplot\\.(pdf|png)$", d = "Per-class read-size barplots, unique + all reads (step 05)"),
  list(re = "^Figure_05b\\..*_size_barplot\\.percentage\\.(pdf|png)$", d = "Per-class read-size percentage barplots, unique + all reads (step 05)"),
  list(re = "^Figure_06\\.piRNA_vs_(tRNA|snoRNA)_pos_barplot(_AS)?\\.(pdf|png)$", d = "Position-distribution of piRNA reads over tRNA/snoRNA gene windows; _AS = opposite-strand (antisense) panel, unique + all reads (step 06)"),
  list(re = "^Figure_06\\.piRNA_overlap_summary\\.(pdf|png)$", d = "Percentage of piRNA reads overlapping snoRNA/tRNA genes (sense/antisense/none) (step 06)"),
  list(re = "^Figure_s01a_overlap_rules_composition\\.(pdf|png)$", d = "Overlap-rule annotation composition (step s01)"),
  list(re = "^Figure_s01b_overlap_rules_category_size\\.(pdf|png)$", d = "Overlap-rule per-category read-size distributions (step s01)"),
  list(re = "^Figure_03a_per_locus_distribution\\.(pdf|png)$", d = "Per-locus log10 abundance distribution per category (step 03)"),
  list(re = "^Figure_03b_rank_abundance\\.(pdf|png)$", d = "Rank-abundance (Whittaker) curves per category (step 03)"),
  list(re = "^Figure_03c_lorenz_matmiRNA\\.(pdf|png)$", d = "Lorenz curves for mature-miRNA loci (step 03)")
)
describe.file = function(f){
  for (r in desc.rules) if (grepl(r$re, f)) return(r$d)
  "See producing step for details"
}

dir.map <- list(
  "tables"  = "Step 01: read-preprocessing tables (shared)",
  "rdata"   = "Step 01: read GRanges objects (shared)",
  "figures" = "Steps 01, 03-06: read-size, annotation, abundance and position-distribution figures"
)
if (is.comparison)
  for (st in active.strats){
    sub = file.path(strategy.dir[[st]])
    dir.map[[file.path(sub, "tables")]]  = sprintf("Steps 02, 03, 10: %s-strategy count tables, abundance and this summary", st)
    dir.map[[file.path(sub, "rdata")]]   = sprintf("Steps 02, 06: %s-strategy annotated-read GRanges and window objects", st)
    dir.map[[file.path(sub, "figures")]] = sprintf("Steps 03-06: %s-strategy abundance, annotation, size and position figures", st)
  }
n.files = 0L
for (d in names(dir.map)){
  full = file.path(out.base, d)
  if (!dir.exists(full)) { rep("### `%s`\n_(not run in this analysis)_\n", d); next }
  files = list.files(full)
  if (length(files) == 0) { rep("### `%s`\n_(empty)_\n", d); next }
  rep("### `%s`\n\n_%s_\n", d, dir.map[[d]])
  for (f in files){
    sz = format(file.info(file.path(full, f))$size, big.mark = ",")
    rep("- `%s` (%s B) - %s", f, sz, describe.file(f))
    n.files = n.files + 1L
  }
  rep("")
}
if (is.comparison && length(active.strats) < length(strat.names))
  rep("_(strategies not found: %s)_\n",
      paste(setdiff(strat.names, active.strats), collapse = ", "))

rep("### Ancillary inputs\n")
rep("- `scripts/R/lib/init.R` - shared bootstrap: genome/strategy from the environment, samples auto-detected from `bams/`, output paths derived")
rep("- `DB/rdata_%s/` - feature GRanges DB (%s objects), built by step 00", genome, length(list.files(db.dir, pattern = "\\.RData$")))
rep("- `scripts/logs/step*.log` - console log of each step")
rep("")

## -----------------------------------------------------------------------------
## 2. Results summary
## -----------------------------------------------------------------------------
rep("## 2. Results summary\n")

## ---- 2.1 sample-level -------------------------------------------------------
t1a = read.tab(file.path(out.base, "tables/Table_01a_sample_summary.csv"))
if (!is.null(t1a)){
  rep("### 2.1 Sample-level read statistics (Table_01a)\n")
  rep("| sample | total reads | unique reads | singletons | %% singletons | median count | max count | dominant size by reads |")
  rep("|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(t1a)))
    rep("| %s | %s | %s | %s | %.1f | %s | %s | %d nt |",
        t1a$sample[i], fmt(t1a$total_reads[i]), fmt(t1a$unique_reads[i]),
        fmt(t1a$singleton_reads[i]), t1a$pct_singletons[i],
        fmt(t1a$median_count[i]), fmt(t1a$max_count[i]), t1a$dominant_size_by_reads_nt[i])
  rep("")
}

## ---- 2.2 read-size composition ----------------------------------------------
t1b = read.tab(file.path(out.base, "tables/Table_01b_read_size_distribution.csv"))
if (!is.null(t1b)){
  rep("### 2.2 Read-size composition (Table_01b)\n")
  rep("| sample | reads 20-22 nt | %% of all reads | dominant size |")
  rep("|---|---|---|---|")
  for (s in samples){
    x = t1b[sample == s]
    if (nrow(x) == 0) next
    f = pct2(sum(x$n_reads[x$width %in% 20:22]) / sum(x$n_reads))
    rep("| %s | %s | %.2f | %d nt |", s, fmt(sum(x$n_reads[x$width %in% 20:22])), f,
        x$width[which.max(x$n_reads)])
  }
  rep("")
}

## ---- 2.3 annotation composition (all reads) ---------------------------------
comp.md = function(t2, tag){
  x = t2[category == "read.annotation"]
  if (nrow(x) == 0) return(invisible())
  tot = x[, .(tot = sum(Freq)), by = sample]
  x = merge(x, tot, by = "sample")
  x[, pct := round(100 * Freq / tot, 2)]
  setorder(x, sample, -Freq)
  rep("### %s\n", tag)
  rep("| sample | category | reads | %% of reads |")
  rep("|---|---|---|---|")
  for (i in seq_len(nrow(x)))
    rep("| %s | %s | %s | %.2f |", x$sample[i], x$item[i], fmt(x$Freq[i]), x$pct[i])
  rep("")
}

t2b.list = t2a.list = list()
for (st in active.strats){
  t2b.list[[st]] = read.tab(file.path(strat.path(st), "tables/Table_02b_annotation_count_all_reads.csv"))
  t2a.list[[st]] = read.tab(file.path(strat.path(st), "tables/Table_02a_annotation_count_unique_reads.csv"))
  if (!is.null(t2b.list[[st]]))
    comp.md(t2b.list[[st]], sprintf("2.3a Annotation composition, all reads (strategy `%s`)", st))
  if (!is.null(t2a.list[[st]]))
    comp.md(t2a.list[[st]], sprintf("2.3b Annotation composition, unique reads (strategy `%s`)", st))
}

## ---- 2.4 mature-miRNA read sizes (20-22 nt) ---------------------------------
rep("### 2.4 Mature-miRNA read sizes (20-22 nt; step 02 size categories)\n")
rep("| strategy | sample | matmiRNA reads | matmiRNA reads 20-22 nt | %% |")
rep("|---|---|---|---|---|")
for (st in active.strats){
  t2b = t2b.list[[st]]
  if (is.null(t2b)) next
  for (s in samples){
    mt = t2b[category == "matmiRNA.size" & sample == s]
    if (nrow(mt) == 0) next
    mt[, w := as.integer(item)]
    f = pct2(sum(mt$Freq[mt$w %in% 20:22]) / sum(mt$Freq))
    rep("| %s | %s | %s | %s | %.2f |", st, s, fmt(sum(mt$Freq)),
        fmt(sum(mt$Freq[mt$w %in% 20:22])), f)
  }
}
rep("")

## ---- 2.5 top mature miRNAs ---------------------------------------------------
for (st in active.strats){
  csvs = file.path(strat.path(st), "tables", paste0("Table_02_", samples, "_unique_reads_annotation.csv"))
  if (!all(file.exists(csvs))) next
  t2m = rbindlist(lapply(csvs, function(csv){
    dt = fread(csv)
    dt[category == "matmiRNA" & !is.na(feature_id),
       .(n_unique = .N, n_reads = sum(count)), by = .(sample, name = feature_id)]
  }))
  if (nrow(t2m) == 0) next
  rep("### 2.5 Top expressed mature miRNAs - strategy `%s` (derived from `Table_02_<sample>_unique_reads_annotation.csv`)\n", st)
  for (s in samples){
    x = t2m[sample == s][order(-n_reads)][1:min(10, .N)]
    rep("**%s**\n", s)
    rep("| rank | mature miRNA | n_unique | n_reads |")
    rep("|---|---|---|---|")
    for (i in seq_len(nrow(x)))
      rep("| %d | %s | %s | %s |", i, x$name[i], fmt(x$n_unique[i]), fmt(x$n_reads[i]))
    rep("")
  }
}

## ---- 2.6 per-category abundance (step 03, if run) ------------------------------
for (st in active.strats){
  t5a = read.tab(file.path(strat.path(st), "tables/Table_03a_category_abundance_summary.csv"))
  if (is.null(t5a)) next
  rep("### 2.6 Per-category abundance - strategy `%s` (Table_03a, step 03)\n", st)
  rep("| category | sample | n_loci | total reads | Gini | top1 %% | top5 %% | top10 %% | top locus |")
  rep("|---|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(t5a)))
    rep("| %s | %s | %s | %s | %.3f | %.1f | %.1f | %.1f | %s |",
        t5a$category[i], t5a$sample[i], fmt(t5a$n_loci[i]), fmt(t5a$total_reads[i]),
        t5a$gini[i], t5a$top1_pct[i], t5a$top5_pct[i], t5a$top10_pct[i], t5a$top_locus[i])
  rep("")
}

## ---- 2.7 per-unique-read annotation table ------------------------------------
for (st in active.strats){
  rep("### 2.7 Per-unique-read annotation + abundance table - strategy `%s` (step 02)\n", st)
  rep("| sample | unique reads | total reads | categories | top feature (reads) |")
  rep("|---|---|---|---|---|")
  for (s in samples){
    csv = file.path(strat.path(st), "tables", paste0("Table_02_", s, "_unique_reads_annotation.csv"))
    if (!file.exists(csv)) next
    dt = fread(csv)
    top = dt[which.max(count)]
    rep("| %s | %s | %s | %d | %s | %s |",
        s, fmt(nrow(dt)), fmt(sum(dt$count)), uniqueN(dt$category),
        ifelse(is.na(top$feature_id), "other", top$feature_id), fmt(top$count))
  }
  rep("")
}

## -----------------------------------------------------------------------------
## 3. Sanity checks
## -----------------------------------------------------------------------------
rep("## 3. Sanity checks\n")
rep("Each check returns PASS, FAIL or WARN (WARN = input not available, e.g. optional step not run).\n")

checks = data.table(id = character(), description = character(), status = character(), detail = character())
check = function(id, desc, ok, detail = ""){
  st = if (is.na(ok)) "WARN" else if (ok) "PASS" else "FAIL"
  checks <<- rbind(checks, data.table(id = id, description = desc, status = st, detail = as.character(detail)))
}

total.reads = if (!is.null(t1a)) setNames(as.numeric(t1a$total_reads), t1a$sample) else NULL
uniq.reads  = if (!is.null(t1a)) setNames(as.numeric(t1a$unique_reads), t1a$sample) else NULL

for (st in active.strats){
  t2b = t2b.list[[st]]
  t2a = t2a.list[[st]]
  sdir = strat.path(st)
  tag = ifelse(st.label(st) == "", "", paste0(st.label(st), "."))
  for (s in samples){
    ## totals
    if (!is.null(t2b)){
      v = t2b[category == "read.annotation" & sample == s, sum(Freq)]
      check(paste0("s1_", tag, s), sprintf("sum of all category counts equals total reads (%s, strategy %s)", s, st),
            !is.null(total.reads) && isTRUE(all.equal(v, as.numeric(total.reads[s]))),
            sprintf("observed %s, expected %s", fmt(v), fmt(total.reads[s])))
    }
    if (!is.null(t2a)){
      v = t2a[category == "read.annotation" & sample == s, sum(Freq)]
      check(paste0("s2_", tag, s), sprintf("sum of unique-read category counts equals unique reads (%s, strategy %s)", s, st),
            !is.null(uniq.reads) && isTRUE(all.equal(v, as.numeric(uniq.reads[s]))),
            sprintf("observed %s, expected %s", fmt(v), fmt(uniq.reads[s])))
    }
    ## unique-reads annotation CSV
    csv = file.path(sdir, "tables", paste0("Table_02_", s, "_unique_reads_annotation.csv"))
    if (file.exists(csv)){
      dt = fread(csv)
      check(paste0("s3_", tag, s), sprintf("unique-reads table row count equals unique reads (%s, strategy %s)", s, st),
            isTRUE(all.equal(nrow(dt), as.numeric(uniq.reads[s]))),
            sprintf("observed %s, expected %s", fmt(nrow(dt)), fmt(uniq.reads[s])))
      check(paste0("s4_", tag, s), sprintf("unique-reads table sum(count) equals total reads (%s, strategy %s)", s, st),
            isTRUE(all.equal(sum(dt$count), as.numeric(total.reads[s]))),
            sprintf("observed %s, expected %s", fmt(sum(dt$count)), fmt(total.reads[s])))
      if (!is.null(t2b)){
        v = t2b[category == "read.annotation" & item == "matmiRNA" & sample == s, Freq]
        check(paste0("s5_", tag, s), sprintf("Table_02b matmiRNA equals sum of matmiRNA reads in unique-reads table (%s, strategy %s)", s, st),
              isTRUE(all.equal(v, sum(dt$count[dt$category == "matmiRNA"]))),
              sprintf("observed %s, expected %s", fmt(sum(dt$count[dt$category == "matmiRNA"])), fmt(v)))
      }
      check(paste0("s6_", tag, s), sprintf("unique-reads table sorted by count descending, all counts >= 1 (%s, strategy %s)", s, st),
            dt$count[1] == max(dt$count) && all(dt$count >= 1), "")
      check(paste0("s7_", tag, s), sprintf("sum(cpm) ~ 1e6 in unique-reads table (%s, strategy %s)", s, st),
            isTRUE(abs(sum(dt$cpm) - 1e6) / 1e6 < 1e-3),
            sprintf("observed %s", fmt(round(sum(dt$cpm))))
      )
    }
    ## matmiRNA 20-22 nt fraction (union/any rules deliberately add long reads
    ## spanning a miRNA locus, so the expected fraction is lower than for
    ## fully-contained)
    if (!is.null(t2b)){
      mt = t2b[category == "matmiRNA.size" & sample == s]
      if (nrow(mt) > 0){
        mt[, w := as.integer(item)]
        f = sum(mt$Freq[mt$w %in% 20:22]) / sum(mt$Freq)
        thresh = if (st == "fully-contained") 0.90 else 0.85
        check(paste0("s9_", tag, s), sprintf("matmiRNA reads are predominantly 20-22 nt (>= %.2f) (%s, strategy %s)", thresh, s, st),
              f >= thresh, sprintf("observed %.4f", f))
      }
    }
    ## n_unique <= n_reads per item
    if (!is.null(t2a) && !is.null(t2b)){
      a = t2a[category == "read.annotation" & sample == s]
      b = t2b[category == "read.annotation" & sample == s]
      m = merge(a, b, by = "item", suffixes = c(".u", ".r"))
      ok = all(m$Freq.u <= m$Freq.r)
      check(paste0("s10_", tag, s), sprintf("unique-read count <= all-read count for every category (%s, strategy %s)", s, st), ok, "")
    }
  }
  ## Table_03a vs Table_02b consistency. Table_03a sums per-locus abundance, so a
  ## matmiRNA read overlapping n loci contributes n x count (n = n_features in
  ## the unique-reads table). The expected total is therefore the sum over
  ## matmiRNA reads of count*n_features (equals Table_02b matmiRNA when every
  ## matmiRNA read overlaps exactly one locus, as for fully-contained/union).
  t5a = read.tab(file.path(sdir, "tables/Table_03a_category_abundance_summary.csv"))
  if (!is.null(t5a) && !is.null(t2b)){
    for (s in samples){
      v2 = t2b[category == "read.annotation" & item == "matmiRNA" & sample == s, Freq]
      v5 = t5a[category == "matmiRNA" & sample == s, total_reads]
      csv = file.path(sdir, "tables", paste0("Table_02_", s, "_unique_reads_annotation.csv"))
      exp.reads = v2
      if (file.exists(csv)){
        dt = fread(csv)
        mat = dt[category == "matmiRNA"]
        exp.reads = sum(mat$count * pmax(mat$n_features, 1))
      }
      if (length(v5) == 1)
        check(paste0("s13_", tag, s), sprintf("Table_03a matmiRNA total equals locus-assignment sum from unique-reads table (%s, strategy %s)", s, st),
              isTRUE(all.equal(as.numeric(v5), as.numeric(exp.reads))),
              sprintf("observed %s, expected %s (Table_02b matmiRNA: %s)", fmt(v5), fmt(exp.reads), fmt(v2)))
    }
  }
}

## read-size consistency (shared step 1)
if (!is.null(t1b)){
  for (s in samples){
    x = t1b[sample == s]
    ok = isTRUE(all.equal(sum(x$n_reads), as.numeric(total.reads[s]))) &&
         isTRUE(all.equal(sum(x$n_unique), as.numeric(uniq.reads[s])))
    check(paste0("s8_", s), sprintf("Table_01b read sizes sum to total/unique reads (%s)", s), ok, "")
  }
}

## feature DB
db.files = list.files(db.dir, pattern = "\\.RData$")
check("s14_db", "feature DB has the expected 24 objects", length(db.files) >= 24,
      sprintf("observed %d objects", length(db.files)))
if (file.exists(file.path(db.dir, "matmiRNA.gr.RData"))){
  e = new.env(); load(file.path(db.dir, "matmiRNA.gr.RData"), envir = e)
  n.mir = length(e$matmiRNA.gr)
  if (genome == "mm39")
    check("s15_db", "mature miRNA loci in DB == 2110 (mm39)", n.mir == 2110, sprintf("observed %d", n.mir))
  else
    check("s15_db", sprintf("mature miRNA loci in DB (reported for %s)", genome), NA, sprintf("observed %d", n.mir))
}

## logs: all step*.log present and error-free
log.dir = file.path(project.root, "scripts", "logs")
logs = list.files(log.dir, pattern = "^step.*\\.log$")
if (length(logs) == 0){
  check("s16_logs", "all step logs present", FALSE, paste("no step*.log found in", log.dir))
}else{
  check("s16_logs", "all step logs present", TRUE, paste(length(logs), "logs found"))
  bad = sapply(logs, function(l){
    txt = readLines(file.path(log.dir, l), warn = FALSE)
    any(grepl("^Error|Execution halted", txt))
  })
  check("s17_logs", "no fatal errors in step logs", !any(bad),
        ifelse(any(bad), paste("errors in:", paste(logs[bad], collapse = ", ")), ""))
}

## ---- sanity-check summary table ---------------------------------------------
rep("### Sanity-check results\n")
rep("| id | description | status | detail |")
rep("|---|---|---|---|")
for (i in seq_len(nrow(checks)))
  rep("| %s | %s | **%s** | %s |", checks$id[i], checks$description[i], checks$status[i],
      ifelse(checks$detail[i] == "", "-", checks$detail[i]))
rep("")
rep("**Overall**: %d PASS, %d FAIL, %d WARN",
    sum(checks$status == "PASS"), sum(checks$status == "FAIL"), sum(checks$status == "WARN"))

## -----------------------------------------------------------------------------
## write outputs
## -----------------------------------------------------------------------------
report.file = file.path(out.base, "pipeline_summary.md")
writeLines(rep.lines, report.file)
fwrite(checks, file.path(out.base, "tables/Table_10_sanity_checks.csv"))

cat(sprintf("Pipeline summary written to %s\n", report.file))
cat(sprintf("Sanity checks: %d PASS, %d FAIL, %d WARN; %d files described\n",
            sum(checks$status == "PASS"), sum(checks$status == "FAIL"),
            sum(checks$status == "WARN"), n.files))
