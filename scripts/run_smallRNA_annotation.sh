#!/usr/bin/env bash
## Small-RNA annotation pipeline runner.
## Usage: run_smallRNA_annotation.sh [GENOME] [STRATEGY]
##
##   GENOME   reference genome assembly, default mm39 (also SMALLRNA_GENOME).
##            The feature DB for the genome is auto-detected
##            (DB/rdata_<GENOME>/ with the 24 feature objects) and, if missing
##            or incomplete, rebuilt automatically via step 00
##            (scripts/R/00_build_DB/00_build_annotation_DB.R). Set
##            SMALLRNA_FORCE_REBUILD_DB=1 to force a rebuild of an existing DB
##            (e.g. after changing the DB build script).
##   STRATEGY sense-annotation rule:
##              fully-contained  (default) read fully contained in the feature
##              union                     read within feature OR feature within read
##              any                       any overlap (>= 1 bp) between read and feature
##              comparison                runs fully-contained, union and any (each
##                                        through steps 02/03/04/05/06) and then
##                                        performs the overlap-rule comparison
##                                        (step s01) and the end-of-pipeline
##                                        summary (step 10)
##
## Output layout (see scripts/R/lib/init.R):
##   single strategy -> output/          (tables/figures/rdata)
##   comparison      -> output/comparison/  (shared step-01 at the root, one
##                                           subfolder per strategy + comparison
##                                           analysis and summary at the root)

usage() {
    cat <<'EOF'
Small-RNA annotation pipeline runner.

Usage:
  run_smallRNA_annotation.sh [GENOME] [STRATEGY]
  run_smallRNA_annotation.sh -h | --help

Arguments:
  GENOME   reference genome assembly, default mm39 (also SMALLRNA_GENOME).
           The feature DB for the genome is auto-detected (DB/rdata_<GENOME>/
           with the 24 feature objects) and, if missing or incomplete, rebuilt
           automatically via step 00 (scripts/R/00_build_DB/00_build_annotation_DB.R).
  STRATEGY sense-annotation rule (also SMALLRNA_STRATEGY):
             fully-contained  (default) read fully contained in the feature
             union                     read within feature OR feature within read
             any                       any overlap (>= 1 bp) between read and feature
             comparison                runs fully-contained, union and any (each
                                       through steps 02/03/04/05/06), then the
                                       overlap-rule comparison (step s01) and the
                                       end-of-pipeline summary (step 10)

Output layout (see scripts/R/lib/init.R):
  single strategy -> output/            (tables/figures/rdata)
  comparison      -> output/comparison/ (shared step-01 at the root, one subfolder
                                         per strategy (fully_contained/, union/,
                                         any/) plus the comparison analysis and
                                         summary at the root)

Environment:
  SMALLRNA_GENOME              same as [GENOME]
  SMALLRNA_STRATEGY            same as [STRATEGY]
  SMALLRNA_FORCE_REBUILD_DB=1  force a feature-DB rebuild even if it is complete
  SMALLRNA_SUBSTRATEGY         set by the runner inside comparison mode (do not set)

Examples:
  run_smallRNA_annotation.sh                 # mm39, fully-contained
  run_smallRNA_annotation.sh hg38 union      # human, union rule
  run_smallRNA_annotation.sh hg38 comparison # all three rules + comparison
EOF
}

set -e
for arg in "$@"; do
    case "$arg" in
        -h|--help|help) usage; exit 0 ;;
    esac
done
GENOME="${1:-${SMALLRNA_GENOME:-mm39}}"
STRATEGY="${2:-${SMALLRNA_STRATEGY:-fully-contained}}"
case "$STRATEGY" in
  fully-contained|fully_contained|within) STRATEGY="fully-contained" ;;
  union)        STRATEGY="union" ;;
  any)          STRATEGY="any" ;;
  comparison|compar) STRATEGY="comparison" ;;
  *) echo "ERROR: unknown STRATEGY '$STRATEGY' (use fully-contained|union|any|comparison)"; exit 1 ;;
esac
export SMALLRNA_GENOME=$GENOME
export SMALLRNA_STRATEGY=$STRATEGY
echo -ne "Reference genome: $GENOME, strategy: $STRATEGY. The whole pipeline starts at: "
date

## clean previous results (output base for this run + logs)
BASE=../output
if [ "$STRATEGY" == "comparison" ]; then BASE=../output/comparison; fi
rm -rf "$BASE"
mkdir -p "$BASE/tables" "$BASE/rdata" "$BASE/figures"
rm -rf logs
mkdir -p logs

## run one pipeline script from its stage directory, logging to scripts/logs/
run_R() {
    local stage_dir="$1" script="$2" log="$3"
    echo -ne "$log starts at: "
    date
    ( cd "$stage_dir" && Rscript "$script" > "../../logs/$log" 2>&1 )
    echo -ne "$log finished at: "
    date
}

## feature DB: auto-detect (DB/rdata_<GENOME>/ with >= 24 feature objects) and,
## if missing or incomplete, rebuild it via step 00 (logs to logs/step00.log).
DB_DIR="../DB/rdata_${GENOME}"
DB_ORIG="../DB/original_data_${GENOME}"
DB_N_OBJECTS=24
db_built() {
    [ -d "$DB_DIR" ] && [ "$(ls -1 "$DB_DIR"/*.RData 2>/dev/null | wc -l)" -ge "$DB_N_OBJECTS" ]
}
if [ "${SMALLRNA_FORCE_REBUILD_DB:-0}" = "1" ] || ! db_built; then
    echo "Feature DB $DB_DIR is missing or incomplete; checking raw files ..."
    missing=""
    for f in refGene.gtf miRNA.gff3 tRNAs.bed RM.bed piRBase.bed piRNAdb.gtf; do
        if [ ! -f "$DB_ORIG/$f" ]; then missing="$missing $f"; fi
    done
    if [ -n "$missing" ]; then
        echo "ERROR: cannot build feature DB for $GENOME: raw files missing in $DB_ORIG:"
        for f in $missing; do echo "  - $f"; done
        exit 1
    fi
    echo "Building feature DB for $GENOME ..."
    run_R R/00_build_DB 00_build_annotation_DB.R step00.log
    if ! db_built; then
        echo "ERROR: feature DB build failed for $GENOME (see logs/step00.log)"
        exit 1
    fi
    echo "Feature DB for $GENOME is ready."
fi

if [ "$STRATEGY" == "comparison" ]; then
    ## step 01: shared pre-processing (runs once)
    run_R R/01_preprocess 01_format_bam.R step01.log

    ## steps 02/03/04/05/06 for each strategy
    for SUB in fully-contained union any; do
        export SMALLRNA_SUBSTRATEGY=$SUB
        echo "===== strategy: $SUB ====="
        run_R R/02_annotation 02_annotation_smallRNA.R "step02_${SUB}.log"
        run_R R/03_figures 03_abundance_by_category.R "step03_${SUB}.log"
        run_R R/03_figures 04_figure_annotation.R "step04_${SUB}.log"
        run_R R/03_figures 05_figure_size.R "step05_${SUB}.log"
        run_R R/03_figures 06_tRNA_snoRNA_position.windows.R "step06_${SUB}.log"
    done
    unset SMALLRNA_SUBSTRATEGY

    ## step s01: overlap-rule comparison analysis
    run_R R/03_figures s01_compare_overlap_rules.R steps01.log
else
    ## step 01
    run_R R/01_preprocess 01_format_bam.R step01.log

    ## step 02
    run_R R/02_annotation 02_annotation_smallRNA.R step02.log

    ## steps 03/04/05/06
    run_R R/03_figures 03_abundance_by_category.R step03.log
    run_R R/03_figures 04_figure_annotation.R step04.log
    run_R R/03_figures 05_figure_size.R step05.log
    run_R R/03_figures 06_tRNA_snoRNA_position.windows.R step06.log
fi

## step 10: end-of-pipeline summary
run_R R/04_summary 10_pipeline_summary.R step10.log

echo -ne "The whole pipeline ends at: "
date
