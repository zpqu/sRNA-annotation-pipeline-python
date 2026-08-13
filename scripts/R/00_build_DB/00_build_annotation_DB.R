##Author: Zhipeng
## This script rebuilds the small-RNA annotation genomic-feature DB for the
## reference genome selected via the environment (assembly ID `genome`, e.g.
## "mm39" or "hg38"). Raw annotation files are read from
## DB/original_data_<genome>/ (refGene.gtf, miRNA.gff3, tRNAs.bed, RM.bed,
## piRBase.bed, piRNAdb.gtf) and one .RData per feature object is written to
## DB/rdata_<genome>/.
##
## The set of GRanges / GRangesList objects mirrors the old mm10 DB, but the
## genome prefix is dropped: object/file names are generic (matmiRNA.gr,
## refGene.NM.exon.grl, ...), so the annotation scripts are genome-agnostic.
##
## Run from scripts/R/00_build_DB/ (same convention as the other pipeline
## scripts), or invoke directly: Rscript scripts/R/00_build_DB/00_build_annotation_DB.R

suppressMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(data.table)
})

source("../lib/init.R")
in.dir  <- orig.dir
out.dir <- db.dir
stopifnot(dir.exists(in.dir))
dir.create(out.dir, showWarnings = FALSE)

msg <- function(...) cat(sprintf("[0_build_%s] ", genome), paste0(...), "\n", sep = "")
save.one <- function(obj) {
  obj.name <- deparse(substitute(obj))
  save(list = obj.name, file = file.path(out.dir, paste0(obj.name, ".RData")))
  msg("saved ", obj.name, " (n=", if (is(obj, "GRangesList")) length(unlist(obj)) else length(obj), ")")
}

## keep the same metadata column order as the old mm10 refGene-derived objects
refGene.cols <- c("source", "type", "score", "phase", "gene_id", "transcript_id",
                  "exon_number", "exon_id", "gene_name")
reorder.mcols <- function(gr) {
  mcols(gr) <- mcols(gr)[, refGene.cols, drop = FALSE]
  gr
}

##----------------------------------------------------------------------------
## chromosome-name normalization (applied to every raw source so that read
## coordinates (normalized to "chr" in step 1) and feature coordinates always
## agree, regardless of the source convention: chr1/1/Chr1 -> chr1, MT -> chrM)
##----------------------------------------------------------------------------
normalize.seqnames <- function(x) {
  x <- gsub("^Chr", "chr", x)
  if (chr.style == "chr") {
    x <- ifelse(grepl("^chr", x), x, paste0("chr", x))
    x <- ifelse(x == "chrMT", "chrM", x)
    x
  } else if (chr.style == "") {
    gsub("^chr", "", x)
  } else {
    x
  }
}
## keep only primary chromosomes / scaffolds (drop alts, randoms, haplotypes)
primary.chr <- function(x) {
  !grepl("_alt|_random|_fix|_hap|chrUn|^Un|^HSCHR", x)
}
## subset to primary chromosomes, then rename seqlevels via a named old->new map
normalize.gr <- function(gr) {
  lev.old <- seqlevels(gr)
  keep <- primary.chr(lev.old)
  gr <- gr[as.character(seqnames(gr)) %in% lev.old[keep]]
  lev.new <- normalize.seqnames(lev.old[keep])
  names(lev.new) <- lev.old[keep]
  seqlevels(gr) <- lev.new
  gr
}

##----------------------------------------------------------------------------
## 1. refGene GTF (UCSC) -> protein-coding NM_ + noncoding NR_ features
##----------------------------------------------------------------------------
msg("importing refGene.gtf ...")
gtf <- import(file.path(in.dir, "refGene.gtf"), format = "gtf")
gtf <- normalize.gr(gtf)
tx_id     <- gtf$transcript_id
gene_name <- gtf$gene_name
is_NM     <- grepl("^NM_", tx_id)
is_NR     <- grepl("^NR_", tx_id)
is_snord  <- !is.na(gene_name) & grepl("^Snor", gene_name, ignore.case = TRUE)  # Snord*/Snora* (mouse), SNORD*/SNORA* (human)

## --- NM_ (protein coding): exon / CDS / 5UTR / 3UTR ---
refGene.NM.exon.gr <- reorder.mcols(gtf[is_NM & gtf$type == "exon"])
refGene.NM.CDS.gr  <- reorder.mcols(gtf[is_NM & gtf$type == "CDS"])
refGene.NM.5UTR.gr <- reorder.mcols(gtf[is_NM & gtf$type == "5UTR"])
refGene.NM.3UTR.gr <- reorder.mcols(gtf[is_NM & gtf$type == "3UTR"])

## split by transcript (only transcripts carrying the feature)
refGene.NM.exon.grl <- split(refGene.NM.exon.gr, refGene.NM.exon.gr$transcript_id)
refGene.NM.CDS.grl  <- split(refGene.NM.CDS.gr,  refGene.NM.CDS.gr$transcript_id)
refGene.NM.5UTR.grl <- split(refGene.NM.5UTR.gr, refGene.NM.5UTR.gr$transcript_id)
refGene.NM.3UTR.grl <- split(refGene.NM.3UTR.gr, refGene.NM.3UTR.gr$transcript_id)

## --- NM_ transcript spans / introns / +/-1kb flanks ---
exons_by_tx <- split(refGene.NM.exon.gr, refGene.NM.exon.gr$transcript_id)
tx_ranges   <- range(exons_by_tx)                    # transcript span per NM_ transcript
intron_list <- setdiff(tx_ranges, exons_by_tx)       # introns per NM_ transcript

tx_gr   <- unlist(tx_ranges, use.names = TRUE)       # one range per transcript
tx_ids  <- names(tx_gr)
tx_lvls <- names(tx_ranges)
names(tx_gr) <- NULL
mcols(tx_gr) <- NULL

refGene.NM.intron.grl <- intron_list            # GRangesList, one element per NM_ transcript

refGene.NM.mRNA.gr  <- tx_gr
refGene.NM.mRNA.grl <- split(tx_gr, factor(tx_ids, levels = tx_lvls))

refGene.NM.up1k.gr  <- trim(promoters(tx_gr, upstream = 1000, downstream = 0))
refGene.NM.up1k.grl <- split(refGene.NM.up1k.gr, factor(tx_ids, levels = tx_lvls))

refGene.NM.down1k.gr <- tx_gr
s_plus  <- as.logical(strand(tx_gr) == "+")
s_minus <- as.logical(strand(tx_gr) == "-")
start(refGene.NM.down1k.gr)[s_plus]  <- end(tx_gr)[s_plus] + 1L
end(refGene.NM.down1k.gr)[s_plus]    <- end(tx_gr)[s_plus] + 1000L
end(refGene.NM.down1k.gr)[s_minus]   <- start(tx_gr)[s_minus] - 1L
start(refGene.NM.down1k.gr)[s_minus] <- start(tx_gr)[s_minus] - 1000L
refGene.NM.down1k.gr <- trim(refGene.NM.down1k.gr)
refGene.NM.down1k.grl <- split(refGene.NM.down1k.gr, factor(tx_ids, levels = tx_lvls))

## --- NR_ (noncoding): snoRNAs (Snor*) and other ncRNA ("lincRNA" in old DB) ---
snoRNA.gr       <- reorder.mcols(gtf[is_NR & is_snord & gtf$type == "exon"])
refGene.NR.exon.gr <- reorder.mcols(gtf[is_NR & gtf$type == "exon"])
lincRNA.exon.gr  <- reorder.mcols(gtf[is_NR & !is_snord & gtf$type == "exon"])
lincRNA.exon.grl <- split(lincRNA.exon.gr, lincRNA.exon.gr$transcript_id)

##----------------------------------------------------------------------------
## 2. miRNA GFF3 (miRBase v22) -> primary transcripts + mature miRNAs
##----------------------------------------------------------------------------
msg("importing miRNA.gff3 ...")
mir <- import(file.path(in.dir, "miRNA.gff3"), format = "gff3")
mir <- normalize.gr(mir)
primiRNA.gr <- mir[mir$type == "miRNA_primary_transcript"]
matmiRNA.gr <- mir[mir$type == "miRNA"]

##----------------------------------------------------------------------------
## 3. tRNAs (tRNAscan-SE BED12)
##----------------------------------------------------------------------------
msg("importing tRNAs.bed ...")
tRNA.bed <- fread(file.path(in.dir, "tRNAs.bed"), sep = "\t", header = FALSE,
                  select = 1:6, col.names = c("chr", "start", "end", "name", "score", "strand"))
keep <- primary.chr(tRNA.bed$chr)
tRNA.bed <- tRNA.bed[keep]
tRNA.bed[, chr := normalize.seqnames(chr)]
tRNA.gr <- GRanges(seqnames = tRNA.bed$chr,
                   ranges = IRanges(tRNA.bed$start + 1L, tRNA.bed$end),
                   strand = tRNA.bed$strand,
                   source = paste0(genome, "_tRNAs"), type = "exon",
                   score = tRNA.bed$score, phase = NA_integer_,
                   gene_id = tRNA.bed$name, transcript_id = tRNA.bed$name)

##----------------------------------------------------------------------------
## 4. RepeatMasker
##----------------------------------------------------------------------------
msg("importing RM.bed ...")
RM.bed <- fread(file.path(in.dir, "RM.bed"), sep = "\t", header = FALSE,
                select = 1:6, col.names = c("chr", "start", "end", "name", "score", "strand"))
keep <- primary.chr(RM.bed$chr)
RM.bed <- RM.bed[keep]
RM.bed[, chr := normalize.seqnames(chr)]
RM.gr <- GRanges(seqnames = RM.bed$chr,
                 ranges = IRanges(RM.bed$start + 1L, RM.bed$end),
                 strand = RM.bed$strand,
                 source = paste0(genome, "_rmsk"), type = "exon",
                 score = RM.bed$score, phase = NA_integer_,
                 gene_id = RM.bed$name, transcript_id = RM.bed$name)

##----------------------------------------------------------------------------
## 5. piRNAs: piRBase (reference) + non-redundant piRNAdb
##    piRBase is the reference set. piRNAdb entries whose coordinates
##    (chr, start, end, strand) exactly match a piRBase entry are dropped;
##    the remaining piRNAdb loci are appended to the piRBase set.
##----------------------------------------------------------------------------
msg("importing piRBase.bed (this is the big one) ...")
piRNA.bed <- fread(file.path(in.dir, "piRBase.bed"), sep = "\t", header = FALSE,
                   select = 1:6,
                   col.names = c("chr", "start", "end", "name", "score", "strand"),
                   colClasses = list(character = c(1, 4, 6), integer = c(2, 3, 5)))
keep <- primary.chr(piRNA.bed$chr)
piRNA.bed <- piRNA.bed[keep]
piRNA.bed[, chr := normalize.seqnames(chr)]
piRNA.gr <- GRanges(seqnames = piRNA.bed$chr,
                    ranges = IRanges(piRNA.bed$start + 1L, piRNA.bed$end),
                    strand = piRNA.bed$strand,
                    name = piRNA.bed$name, score = piRNA.bed$score)

msg("importing piRNAdb.gtf ...")
piRNAdb.gr <- import(file.path(in.dir, "piRNAdb.gtf"), format = "gtf")
piRNAdb.gr <- normalize.gr(piRNAdb.gr)
## piRNAdb GTF is 1-based inclusive, piRBase GRanges is 1-based -> comparable directly
exact.hits <- findOverlaps(piRNAdb.gr, piRNA.gr, type = "equal")
piRNAdb.new.gr <- piRNAdb.gr[setdiff(seq_along(piRNAdb.gr), queryHits(exact.hits))]
## drop loci duplicated within piRNAdb itself (same coordinates, different codes)
dup.key <- duplicated(paste(seqnames(piRNAdb.new.gr), start(piRNAdb.new.gr),
                            end(piRNAdb.new.gr), strand(piRNAdb.new.gr)))
piRNAdb.new.gr <- piRNAdb.new.gr[!dup.key]
mcols(piRNAdb.new.gr) <- DataFrame(name = piRNAdb.new.gr$piRNA_code, score = NA_integer_)
msg("piRNAdb non-redundant loci added: ", length(piRNAdb.new.gr))
piRNA.gr <- c(piRNA.gr, piRNAdb.new.gr)

##----------------------------------------------------------------------------
## save everything (one .RData per object, mirroring DB/rdata layout)
##----------------------------------------------------------------------------
msg("saving objects to ", out.dir)
save.one(primiRNA.gr)
save.one(matmiRNA.gr)
save.one(snoRNA.gr)
save.one(piRNA.gr)
save.one(tRNA.gr)
save.one(RM.gr)
save.one(refGene.NM.exon.gr)
save.one(refGene.NM.exon.grl)
save.one(refGene.NM.CDS.gr)
save.one(refGene.NM.CDS.grl)
save.one(refGene.NM.5UTR.gr)
save.one(refGene.NM.5UTR.grl)
save.one(refGene.NM.3UTR.gr)
save.one(refGene.NM.3UTR.grl)
save.one(refGene.NM.intron.grl)
save.one(refGene.NM.mRNA.gr)
save.one(refGene.NM.mRNA.grl)
save.one(refGene.NM.up1k.gr)
save.one(refGene.NM.up1k.grl)
save.one(refGene.NM.down1k.gr)
save.one(refGene.NM.down1k.grl)
save.one(refGene.NR.exon.gr)
save.one(lincRNA.exon.gr)
save.one(lincRNA.exon.grl)

msg("done. ", length(list.files(out.dir, pattern = ".RData$")), " objects written to ", out.dir)
