##Date: 21/07/2015
##Author: Zhipeng
## Step 2: small-RNA annotation, mature-miRNA strategy, parameterised by the
## sense-annotation rule. The miRNA category is annotated with MATURE miRNA
## loci (matmiRNA.gr, miRBase v22, ~22 nt). The sense rule is controlled by
## config/bootstrap (init.R: SMALLRNA_STRATEGY / SMALLRNA_SUBSTRATEGY):
##   fully-contained : a read is annotated to a feature only if the read is
##                      fully contained in the feature (findOverlaps type="within")
##   union           : read within feature OR feature within read (union
##                      containment; captures long reads spanning a small feature)
##   any             : any overlap (>= 1 bp) between read and feature
## Antisense annotation is always type="any" with ignore.strand = TRUE.
## Priority: matmiRNA>piRNA>snoRNA>tRNA>RM>refGene_NM_exon>refGene_NM_intron>lincRNA
## >antisense_matmiRNA>antisense_piRNA>antisense_snoRNA>antisense_tRNA>antisense_RM
## >antisense_NM_exon>antisense_NM_intron>antisense_lincRNA
## NOTE: piRNA is deliberately ranked above snoRNA so that reads overlapping both a
## piRNA locus and a snoRNA gene are annotated as piRNA - this lets downstream
## analyses (step 06) test whether piRNA-annotated reads are degradation products
## of snoRNA/tRNA genes.
## Reads are annotated as non-redundant unique reads (chr/start/end/strand).
## Outputs are written to the directory resolved by the shared bootstrap
## (scripts/R/lib/init.R; output/ for a single-strategy run,
## output/comparison/<strategy>/ inside a comparison run):
##   tables/Table_02a_annotation_count_unique_reads.csv
##   tables/Table_02b_annotation_count_all_reads.csv
##   tables/Table_02_<sample>_unique_reads_annotation.csv  (per-read annotation + abundance)
##   rdata/<sample>.bam.annotated.gr.RData

library(GenomicFeatures)
library(data.table)
rm(list = ls())

source("../lib/init.R")
sense.rule <- if (cur.strategy == "union") "union" else if (cur.strategy == "any") "any" else "within"
dir.tab  <- file.path(out.dir, "tables")
dir.rdata <- file.path(out.dir, "rdata")
dir.create(dir.tab, recursive = TRUE, showWarnings = FALSE)
dir.create(dir.rdata, recursive = TRUE, showWarnings = FALSE)
print(sprintf("annotation strategy (sense rule): %s (config cur.strategy=%s)", sense.rule, cur.strategy))

###load genomic features (mature miRNA)
load(file.path(db.dir, "matmiRNA.gr.RData"))
load(file.path(db.dir, "snoRNA.gr.RData"))
load(file.path(db.dir, "piRNA.gr.RData"))
load(file.path(db.dir, "tRNA.gr.RData"))
load(file.path(db.dir, "RM.gr.RData"))
load(file.path(db.dir, "refGene.NM.exon.grl.RData"))
load(file.path(db.dir, "refGene.NM.intron.grl.RData"))
load(file.path(db.dir, "lincRNA.exon.grl.RData"))

#load additional gene features
load(file.path(db.dir, "refGene.NM.CDS.grl.RData"))
load(file.path(db.dir, "refGene.NM.5UTR.grl.RData"))
load(file.path(db.dir, "refGene.NM.3UTR.grl.RData"))
load(file.path(db.dir, "refGene.NM.up1k.grl.RData"))
load(file.path(db.dir, "refGene.NM.down1k.grl.RData"))

####make genomic feature list or gene feature list
genomicFeature.id = c("matmiRNA.gr", "piRNA.gr", "snoRNA.gr",
		  "tRNA.gr", "RM.gr", "refGene.NM.exon.grl",
		  "refGene.NM.intron.grl",
		  "lincRNA.exon.grl")
geneFeature.list = c("refGene.NM.CDS.grl", "refGene.NM.5UTR.grl",
	       "refGene.NM.3UTR.grl", "refGene.NM.intron.grl",
	       "refGene.NM.up1k.grl", "refGene.NM.down1k.grl", "RM.gr")

### helper: id vector for a genomic feature, aligned to the subject index used
### by findOverlaps on that feature (element index for GRangesList features).
feature.id.vector = function(feature.name, tx2gn){
  f = get(feature.name)
  if(is(f, "GRangesList")){
    u = unlist(f)
    if(!is.null(u$transcript_id) && !is.null(u$gene_name)){
      gn = as.character(u$gene_name); gi = as.character(u$gene_id)
      gn[is.na(gn) | gn == ""] = gi[is.na(gi) | gi == ""]
      dt = unique(data.table(tx = as.character(u$transcript_id), id = gn))
      return(dt$id[match(names(f), dt$tx)])
    }
    return(tx2gn$gene_name[match(names(f), tx2gn$transcript_id)])
  }
  id.col = switch(feature.name,
    matmiRNA.gr = "Name", snoRNA.gr = "gene_id", piRNA.gr = "name",
    tRNA.gr = "gene_id", RM.gr = "gene_id")
  as.character(mcols(f)[[id.col]])
}
## transcript -> gene-name map (for intron grl, whose unlisted ranges carry no mcols)
tx2gn = unique(data.table(transcript_id = as.character(unlist(refGene.NM.exon.grl)$transcript_id),
                          gene_name = as.character(unlist(refGene.NM.exon.grl)$gene_name)))
fids.all = lapply(genomicFeature.id, feature.id.vector, tx2gn = tx2gn)
names(fids.all) = genomicFeature.id

### helper: for a list of query/subject hit pairs, return the specific feature id
### (and the number of overlapping features) for each hit read, aligned to
### unique(qh). piRNA uses the highest-scoring locus among the overlapping
### piRNAs; other features use an arbitrary (first) hit. For the union rule the
### "feature within read" direction (feature query, read subject) is added to the
### candidate pairs before the arbitrary/score-based selection.
feature.hit.info = function(qh, sh, feature.name, fids, rule, ignore.strand = FALSE){
  if(length(qh) == 0) return(list(feature.id = character(0), n_features = integer(0)))
  feat = get(feature.name)
  if(rule == "union"){
    ol2 = findOverlaps(feat, test.bam.gr, type = "within")
    qh = c(qh, subjectHits(ol2)); sh = c(sh, queryHits(ol2))
  }
  ## distinct (read, feature) pairs: for the union rule the same pair can be
  ## recorded twice (e.g. a read exactly coinciding with a feature), which would
  ## otherwise inflate n_features.
  keep.pair = !duplicated(paste(qh, sh, sep = ":"))
  qh = qh[keep.pair]; sh = sh[keep.pair]
  rows = unique(qh)
  nf = tabulate(qh, nbins = length(test.bam.gr))[rows]
  if(feature.name == "piRNA.gr"){
    ord = order(qh, -as.numeric(mcols(feat)$score[sh]))
    keep = !duplicated(qh[ord]); map = sh[ord][keep]; names(map) = qh[ord][keep]
    fid = fids[map[as.character(rows)]]
  }else{
    if(rule == "union"){
      keep = !duplicated(qh); map = sh[keep]; names(map) = qh[keep]
      fid = fids[map[as.character(rows)]]
    }else{
      ol = findOverlaps(test.bam.gr[rows], feat, type = rule, ignore.strand = ignore.strand, select = "arbitrary")
      fid = fids[ol]
    }
  }
  list(feature.id = fid, n_features = nf)
}

### sense-hit query/subject pairs for the configured rule
sense.hits = function(reads.gr, feature.gr){
  if(sense.rule == "within"){
    ol = findOverlaps(reads.gr, feature.gr, type = "within")
    return(list(qh = queryHits(ol), sh = subjectHits(ol)))
  }
  if(sense.rule == "union"){
    ol1 = findOverlaps(reads.gr, feature.gr, type = "within")
    ol2 = findOverlaps(feature.gr, reads.gr, type = "within")
    return(list(qh = c(queryHits(ol1), subjectHits(ol2)),
                sh = c(subjectHits(ol1), queryHits(ol2))))
  }
  ol = findOverlaps(reads.gr, feature.gr, type = "any")
  return(list(qh = queryHits(ol), sh = subjectHits(ol)))
}

###
source("mygeneFeature.R")

### count per group: n_unique (each unique read once) and n_reads (weighted by count)
count.by.2 = function(group, counts){
  dt = data.table(Var1 = as.character(group), count = as.numeric(counts))
  dt = dt[, .(n_unique = .N, n_reads = sum(count)), by = Var1]
  return(dt)
}

### append group counts to the two consolidated tables (unique-based and all-reads-based)
unique.all.tab = data.table()
reads.all.tab = data.table()
add.counts = function(category, group, counts){
  dt = count.by.2(group, counts)
  unique.all.tab <<- rbind(unique.all.tab,
                           data.table(sample = sample, category = category,
                                      item = dt$Var1, Freq = dt$n_unique))
  reads.all.tab <<- rbind(reads.all.tab,
                          data.table(sample = sample, category = category,
                                     item = dt$Var1, Freq = dt$n_reads))
}

dir.step1 = file.path(out.base, "rdata")
files = list.files(path = dir.step1, pattern = ".bam.unique.gr.RData$")
if(length(files) == 0) stop("no *.bam.unique.gr.RData found in ", dir.step1, " (run step 1 first)")
for(i in seq(along = files)){
      file.name = file.path(dir.step1, files[i])
      sample = files[i]
      sample = gsub("\\.bam\\.unique\\.gr\\.RData", "", sample)
      sample = gsub("\\.(bwa|bowtie2)$", "", sample)

      print(paste("Now is processing ...", sample, Sys.time()))
      load(file.name)
      test.all.id = c(1:length(reads.bam.unique.gr))

      test.bam.gr = reads.bam.unique.gr
      test.bam.gr$type = NULL
      test.bam.gr$region = NULL
      ## position key for removing already-annotated reads. Names are NOT unique:
      ## multi-mapping reads share one name across loci (step 1 collapses by
      ## chr/start/end/strand), so removal must use this key, not read names.
      test.bam.gr$rm.key = paste(seqnames(test.bam.gr), start(test.bam.gr),
                                 end(test.bam.gr), strand(test.bam.gr), sep = ":")
      test.bam.new.gr = test.bam.gr
      test.bam.new.gr = test.bam.new.gr[-(1:length(test.bam.new.gr))]
      test.other.gr = NULL
      ###miRNA[0]>piRNA[1]>snoRNA[2]>tRNA[3]>RM_refGeneNM[4]>RM[5]>refGene_NM_exon[6]>refGene_NM_intron[7]>lincRNA[8]
      ### >antisense_matmiRNA>antisense_piRNA>antisense_snoRNA>antisense_tRNA>antisense_RM_refGeneNM>antisense_RM>antisense_NM_exon>antisense_NM_intron>anitsense_lincRNA

      ##sense
      for(j in 1:length(genomicFeature.id)){
      	    feature.name = genomicFeature.id[j]
	    feature.id = gsub("\\.grl", "", feature.name)
	    feature.id = gsub("\\.gr", "", feature.id)

	    test.bam.hit.gr = NULL
	    hs = sense.hits(test.bam.gr, get(feature.name))
	    test.bam.hit.gr = test.bam.gr[unique(hs$qh)]
	    test.bam.hit.gr$region = rep("NA", length(test.bam.hit.gr))
	    hit.info = feature.hit.info(hs$qh, hs$sh, feature.name, fids.all[[feature.name]], sense.rule)
	    test.bam.hit.gr$feature.id = hit.info$feature.id
	    test.bam.hit.gr$n_features = hit.info$n_features
	    if(feature.id == "matmiRNA"){
	    		  test.bam.hit.gr = mygeneFeature(bam = "test.bam.hit.gr", type = sense.rule)
			  add.counts("matmiRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("matmiRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "snoRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr", type = sense.rule)
			  add.counts("snoRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("snoRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "piRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr", type = sense.rule)
			  add.counts("piRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("piRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "tRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr", type = sense.rule)
			  add.counts("tRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("tRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    test.bam.hit.gr$type = rep(feature.id, length(test.bam.hit.gr))
	    test.bam.new.gr = c(test.bam.new.gr, test.bam.hit.gr)
	    if(length(unique(hs$qh)) == 0){
	     test.bam.gr = test.bam.gr
	    }else{
	     test.bam.gr = test.bam.gr[!test.bam.gr$rm.key %in% test.bam.gr$rm.key[unique(hs$qh)]]
	    }
      }

      ##antisense (always any-overlap, ignore strand)
      for(k in 1:length(genomicFeature.id)){
            feature.name = genomicFeature.id[k]
            feature.id = gsub("\\.grl", "", feature.name)
            feature.id = gsub("\\.gr", "", feature.id)

            test.bam.hit.gr = NULL
            test.bam.ol = findOverlaps(test.bam.gr, get(feature.name), type = "any", ignore.strand = T)
            test.bam.hit.gr = test.bam.gr[unique(queryHits(test.bam.ol))]
	    test.bam.hit.gr$region = rep("NA", length(test.bam.hit.gr))
	    hit.info = feature.hit.info(queryHits(test.bam.ol), subjectHits(test.bam.ol), feature.name,
	                                fids.all[[feature.name]], "any", ignore.strand = T)
	    test.bam.hit.gr$feature.id = hit.info$feature.id
	    test.bam.hit.gr$n_features = hit.info$n_features
            test.bam.hit.gr$type = rep(paste("AS.", feature.id, sep = ""), length(test.bam.hit.gr))
            test.bam.new.gr = c(test.bam.new.gr, test.bam.hit.gr)
            if(length(unique(queryHits(test.bam.ol))) == 0){
	        test.bam.gr = test.bam.gr
            }else{
	    	test.bam.gr = test.bam.gr[!test.bam.gr$rm.key %in% test.bam.gr$rm.key[unique(queryHits(test.bam.ol))]]
            }

      }

      ##other
      test.bam.other.gr = test.bam.gr
      test.bam.other.gr$region = rep("NA", length(test.bam.other.gr))
      test.bam.other.gr$type = rep("other", length(test.bam.other.gr))
      test.bam.other.gr$feature.id = rep(NA_character_, length(test.bam.other.gr))
      test.bam.other.gr$n_features = rep(NA_integer_, length(test.bam.other.gr))
      test.bam.new.gr = c(test.bam.new.gr, test.bam.other.gr)
      add.counts("read.annotation", test.bam.new.gr$type, test.bam.new.gr$count)
      add.counts("read.size", width(test.bam.new.gr), test.bam.new.gr$count)
      reads.bam.annotated.gr = test.bam.new.gr
      reads.bam.annotated.gr$rm.key = NULL
      save(reads.bam.annotated.gr, file = file.path(dir.rdata, paste(sample, ".bam.annotated.gr.RData", sep = "")))

      ### per-read annotation + abundance table (unique reads, sorted by expression)
      sample.dt = data.table(
        sample = sample,
        read_id = seq_along(reads.bam.annotated.gr),
        chr = as.character(seqnames(reads.bam.annotated.gr)),
        start = start(reads.bam.annotated.gr),
        end = end(reads.bam.annotated.gr),
        strand = as.character(strand(reads.bam.annotated.gr)),
        size = width(reads.bam.annotated.gr),
        count = as.numeric(reads.bam.annotated.gr$count),
        cpm = round(as.numeric(reads.bam.annotated.gr$count) /
                      sum(as.numeric(reads.bam.annotated.gr$count)) * 1e6, 3),
        category = as.character(reads.bam.annotated.gr$type),
        feature_id = as.character(reads.bam.annotated.gr$feature.id),
        gene_context = as.character(reads.bam.annotated.gr$region),
        n_features = as.integer(reads.bam.annotated.gr$n_features))
      setorder(sample.dt, -count, chr, start, end)
      fwrite(sample.dt, file = file.path(dir.tab, paste0("Table_02_", sample, "_unique_reads_annotation.csv")))
}

## ----- consolidated count tables (unique reads vs all reads) -------------------
write.csv(unique.all.tab, file.path(dir.tab, "Table_02a_annotation_count_unique_reads.csv"), row.names = FALSE)
write.csv(reads.all.tab, file.path(dir.tab, "Table_02b_annotation_count_all_reads.csv"), row.names = FALSE)
print(paste("Tables saved to", dir.tab, ": Table_02a_annotation_count_unique_reads.csv, Table_02b_annotation_count_all_reads.csv, Table_02_<sample>_unique_reads_annotation.csv"))
