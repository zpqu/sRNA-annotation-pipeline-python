##Date: 2026
##Author: Zhipeng
## Shared bootstrap sourced by every pipeline R script. The reference genome and
## annotation strategy come from the
## environment, which the runner (scripts/run_smallRNA_annotation.sh) sets from
## its CLI arguments; sample labels are auto-detected from the bams/ folder
## using the same derivation as step 01. All path objects derive from the
## project root, which is resolved from the working directory the script is run
## from (scripts/R/<stage>/, as set by the runner).
##
##   genome     : assembly ID used to locate the feature DB and raw files
##                (e.g. "mm39", "hg38"), from SMALLRNA_GENOME.
##   strategy   : sense-annotation rule, from SMALLRNA_STRATEGY
##                (or SMALLRNA_SUBSTRATEGY inside a comparison run):
##                  "fully-contained"  read fully contained in the feature
##                  "union"            read within feature OR feature within read
##                  "any"              any overlap (>= 1 bp) between read and feature
##                  "comparison"       runs all three rules and then the overlap-rule
##                                     comparison analysis (step s01)
##   samples    : auto-detected from bams/*.bam, sorted for a deterministic order.
##   chr.style  : chromosome naming convention enforced on the feature DB build,
##                from SMALLRNA_CHR_STYLE (default "chr").

genome  <- Sys.getenv("SMALLRNA_GENOME", unset = "mm39")  # env var overrides the default
chr.style <- Sys.getenv("SMALLRNA_CHR_STYLE", unset = "chr")

## project root, resolved regardless of the working directory the script is run from
project.root <- if (basename(getwd()) %in% c("00_build_DB", "lib", "01_preprocess",
                                              "02_annotation", "03_figures", "04_summary")) {
  "../../.."
} else if (file.exists("bams")) {
  "."
} else {
  stop("cannot locate project root from working directory ", getwd())
}

## publication-style figure sizing: fixed 180 mm width, max 220 mm height.
## Returns c(width, height) in inches.  Height is computed from content and
## capped at the maximum; width is always 7.09 in (180 mm).
##   n      : number of facet panels (typically the number of samples)
##   ncol   : number of facet columns
##   per.h  : target panel height (inches)
##   stack  : vertically stacked ggplots sharing the same facet layout
fig.dims = function(n, ncol, per.h = 3.2, stack = 1){
  width.mm  = 180
  height.mm = 220
  width  = width.mm / 25.4          # 7.09 in
  ncol   = max(1L, min(ncol, n))
  rows   = ceiling(n / ncol)
  height = min(stack * rows * per.h, height.mm / 25.4)   # capped at 8.66 in
  c(width = width, height = height)
}

## publication theme: smaller fonts suitable for 180 mm figures
## Defined as a function so ggplot2 only needs to be loaded when the caller uses it.
small.font = function(){
  theme(
    plot.title       = element_text(size = 8),
    axis.title       = element_text(size = 8),
    axis.text        = element_text(size = 7),
    strip.text       = element_text(size = 8),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.25, "cm")
  )
}

## ---- annotation strategy -----------------------------------------------------
canon.strategy = function(x){
  x = tolower(trimws(x))
  x = gsub("[ _]+", "-", x)
  if (grepl("union", x)) return("union")
  if (grepl("any", x))   return("any")
  if (grepl("compar", x)) return("comparison")
  "fully-contained"
}
strategy.mode <- canon.strategy(Sys.getenv("SMALLRNA_STRATEGY", unset = "fully-contained"))
sub.strategy  <- Sys.getenv("SMALLRNA_SUBSTRATEGY", unset = NA_character_)
is.comparison <- strategy.mode == "comparison"
cur.strategy  <- if (!is.na(sub.strategy)) canon.strategy(sub.strategy) else strategy.mode

## directories (relative to project.root)
db.dir    <- file.path(project.root, "DB", paste0("rdata_", genome))         # built feature GRanges
orig.dir  <- file.path(project.root, "DB", paste0("original_data_", genome)) # raw feature files for the DB build

## output layout:
##   single-strategy run (fully-contained | union | any) -> output/
##   comparison run                                       -> output/comparison/ with one
##      subfolder per strategy (fully_contained/, union/, any/) plus shared step-01
##      results and the comparison-analysis tables/figures at the comparison root.
strategy.dir <- c("fully-contained" = "fully_contained", "union" = "union", "any" = "any")
out.base    <- file.path(project.root, if (is.comparison) "output/comparison" else "output")
out.dir     <- if (is.comparison && !is.na(sub.strategy)) {
                 file.path(out.base, unname(strategy.dir[cur.strategy]))
               } else {
                 out.base
               }

## sample labels: auto-detected from bams/ in the same way step 01 derives them
## from the BAM file names (strip ".bam" then a trailing ".bwa"/".bowtie2").
## Sorted for a deterministic display order across figure scripts.
bam.files <- list.files(file.path(project.root, "bams"), pattern = "\\.bam$")
samples <- sort(unique(sub("\\.(bwa|bowtie2)$", "", sub("\\.bam$", "", bam.files))))
if (length(samples) == 0)
  warning("no BAM files found in ", file.path(project.root, "bams"),
          "; 'samples' is empty (step 01 and later steps will fail)")

stopifnot(dir.exists(file.path(project.root, "DB")))
