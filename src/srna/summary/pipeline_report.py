"""Step 10: end-of-pipeline summary report.

Port of ``10_pipeline_summary.R``.

Writes ``<out.base>/pipeline_summary.md`` (output-file inventory, results
summary, sanity checks) and ``<out.base>/tables/Table_10_sanity_checks.csv``.
Single-strategy runs summarise ``output/``; comparison runs summarise
``output/comparison/`` including the per-strategy sub-folders.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from srna.config import Strategy
from srna.features.feature_store import FeatureStore
from srna.paths import PathResolver

_EXPECTED_MATMIRNA = 2110


def _fmt(x) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "-"
    return f"{int(x):,}" if isinstance(x, (int, np.integer)) else f"{x:,.0f}"


def _pct(x: float) -> float:
    return round(100 * x, 2)


class _Report:
    def __init__(self) -> None:
        self.lines: list[str] = []
        self.checks: list[dict] = []

    def rep(self, fmt: str, *args) -> None:
        self.lines.append(fmt % args if args else fmt)

    def check(self, check_id: str, desc: str, ok, detail: str = "") -> None:
        status = "WARN" if ok is None else ("PASS" if ok else "FAIL")
        self.checks.append(
            {"id": check_id, "description": desc, "status": status, "detail": detail}
        )


def _read_tab(path: Path) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path)
    return None


def _active_strategies(
    resolver: PathResolver,
) -> tuple[list[Strategy], callable]:
    """Return the strategies to summarise and the per-strategy output dir."""
    config = resolver.config
    if config.strategy.is_comparison:
        strat = [Strategy.FULLY_CONTAINED, Strategy.UNION, Strategy.ANY]
        return [s for s in strat if resolver.strategy_dir(s).exists()], resolver.strategy_dir
    return [config.strategy], lambda s: resolver.out_dir


def run_pipeline_summary(
    resolver: PathResolver,
    samples: list[str],
    store: FeatureStore,
) -> None:
    """Write the pipeline summary report + sanity-check table."""
    out_base = resolver.output_base
    out_base.mkdir(parents=True, exist_ok=True)
    r = _Report()
    strategies, strat_path = _active_strategies(resolver)
    strat_names = [s.value for s in strategies]

    r.rep("## Small-RNA annotation pipeline - end-of-pipeline summary\n")
    r.rep("- Reference genome : `%s`", resolver.config.genome)
    r.rep("- Feature DB       : `%s`", resolver.db_cache_dir)
    r.rep("- Samples          : %s", ", ".join(samples))
    r.rep("- Strategy mode    : `%s`", resolver.config.strategy.value)
    r.rep("- Output base      : `%s`", out_base)
    r.rep("- Strategies summarized: %s", ", ".join(f"`{s}`" for s in strat_names))
    r.rep("- Generated        : %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    r.rep("")

    # ---- 1. output-file inventory -------------------------------------------
    r.rep("## 1. Output files\n")
    dir_map = {
        "tables": "Step 01: read-preprocessing tables (shared)",
        "rdata": "Step 01: read GRanges objects (shared)",
        "figures": (
            "Steps 01, 03-06: read-size, annotation, abundance and " "position-distribution figures"
        ),
    }
    if resolver.config.strategy.is_comparison:
        for st in strategies:
            sub = st.dir_name
            dir_map[f"{sub}/tables"] = f"Steps 02, 03, 10: {st.value}-strategy count tables"
            dir_map[f"{sub}/rdata"] = f"Steps 02, 06: {st.value}-strategy annotated-read objects"
            dir_map[f"{sub}/figures"] = f"Steps 03-06: {st.value}-strategy figures"
    n_files = 0
    for d, desc in dir_map.items():
        full = out_base / d
        if not full.exists():
            r.rep("### `%s`\n_(not run in this analysis)_\n", d)
            continue
        files = sorted(full.iterdir())
        if not files:
            r.rep("### `%s`\n_(empty)_\n", d)
            continue
        r.rep("### `%s`\n\n_%s_\n", d, desc)
        for f in files:
            if f.is_dir():
                continue
            sz = f"{f.stat().st_size:,}"
            r.rep("- `%s` (%s B) - %s", f.name, sz, describe(f.name))
            n_files += 1
        r.rep("")
    r.rep("")

    # ---- 2. results summary ---------------------------------------------------
    r.rep("## 2. Results summary\n")
    _summary_section(r, resolver, out_base, samples, strategies, strat_path, store)

    # ---- 3. sanity checks ------------------------------------------------------
    _sanity_checks(r, resolver, out_base, samples, strategies, strat_path, store)
    r.rep("### Sanity-check results\n")
    r.rep("| id | description | status | detail |")
    r.rep("|---|---|---|---|")
    for c in r.checks:
        r.rep(
            "| %s | %s | **%s** | %s |",
            c["id"],
            c["description"],
            c["status"],
            c["detail"] if c["detail"] else "-",
        )
    r.rep("")
    n_pass = sum(1 for c in r.checks if c["status"] == "PASS")
    n_fail = sum(1 for c in r.checks if c["status"] == "FAIL")
    n_warn = sum(1 for c in r.checks if c["status"] == "WARN")
    r.rep("**Overall**: %d PASS, %d FAIL, %d WARN", n_pass, n_fail, n_warn)

    report_file = out_base / "pipeline_summary.md"
    report_file.write_text("\n".join(r.lines), encoding="utf-8")
    pd.DataFrame(r.checks).to_csv(out_base / "tables" / "Table_10_sanity_checks.csv", index=False)
    logger.info(
        "step 10: summary written to {} ({} PASS, {} FAIL, {} WARN, {} files)",
        report_file,
        n_pass,
        n_fail,
        n_warn,
        n_files,
    )


def _describe_rules() -> list[tuple[str, str]]:
    return [
        (
            r"^Table_01a_sample_summary\.csv$",
            "Per-sample summary: total reads, unique reads, singletons, dominant read sizes",
        ),
        (
            r"^Table_01b_read_size_distribution\.csv$",
            "Read size distribution (n_unique and n_reads per read size)",
        ),
        (
            r"^Table_01c_read_count_distribution\.csv$",
            "Read count (expression) distribution, log2 bins",
        ),
        (r"^Table_01d_read_size_vs_count\.csv$", "2D read-size x read-count grid"),
        (r"^Table_02a_annotation_count_unique_reads\.csv$", "Annotation composition, unique reads"),
        (
            r"^Table_02b_annotation_count_all_reads\.csv$",
            "Annotation composition, all reads (weighted by read frequency)",
        ),
        (
            r"^Table_02_.*_unique_reads_annotation\.csv$",
            "Per-unique-read annotation + abundance for this sample",
        ),
        (
            r"^Table_s01a_overlap_rule_composition\.csv$",
            "Overlap-rule annotation-composition comparison",
        ),
        (
            r"^Table_s01b_overlap_rule_category_totals\.csv$",
            "Overlap-rule per-category read totals",
        ),
        (
            r"^Table_s01c_read_movement_contained_vs_any\.csv$",
            "Read movement between fully-contained and any strategies",
        ),
        (
            r"^Table_s01d_read_movement_contained_vs_union\.csv$",
            "Read movement between fully-contained and union strategies",
        ),
        (
            r"^Table_s01e_mature_miRNA_expression_strategies\.csv$",
            "Top mature miRNA expression under each overlap-rule strategy",
        ),
        (r"^Table_03a_category_abundance_summary\.csv$", "Per-category abundance statistics"),
        (r"^Table_03b_per_locus_abundance\.csv$", "Per-locus read abundance, long format"),
        (
            r"^Table_06_piRNA_overlap_summary\.csv$",
            "piRNA reads overlapping snoRNA/tRNA genes (sense/antisense/none)",
        ),
        (r"^Table_10_sanity_checks\.csv$", "Pipeline sanity-check results (step 10)"),
        (r"^Figure_01\.", "Read size / count distribution figure, faceted by sample (step 01)"),
        (
            r"^Figure_03[abc]_",
            "Per-locus abundance distribution / rank-abundance / Lorenz curves (step 03)",
        ),
        (r"^Figure_04[a-e]\.", "Annotation count / percentage barplots (step 04)"),
        (r"^Figure_05[ab]\.", "Per-class read-size barplots, unique + all reads (step 05)"),
        (r"^Figure_06\.", "piRNA-vs-tRNA/snoRNA position profiles and overlap summary (step 06)"),
        (r"^Figure_s01[ab]_", "Overlap-rule comparison figures (step s01)"),
    ]


def describe(fname: str) -> str:
    import re

    for pattern, desc in _describe_rules():
        if re.match(pattern, fname):
            return desc
    return "See producing step for details"


def _summary_section(
    r: _Report,
    resolver: PathResolver,
    out_base: Path,
    samples: list[str],
    strategies: list[Strategy],
    strat_path: callable,
    store: FeatureStore,
) -> None:
    t1a = _read_tab(out_base / "tables" / "Table_01a_sample_summary.csv")
    t1b = _read_tab(out_base / "tables" / "Table_01b_read_size_distribution.csv")
    if t1a is not None:
        r.rep("### 2.1 Sample-level read statistics (Table_01a)\n")
        r.rep(
            "| sample | total reads | unique reads | singletons | %% singletons | "
            "median count | max count | dominant size by reads |"
        )
        r.rep("|---|---|---|---|---|---|---|---|")
        for row in t1a.itertuples(index=False):
            r.rep(
                "| %s | %s | %s | %s | %.1f | %s | %s | %d nt |",
                row.sample,
                _fmt(row.total_reads),
                _fmt(row.unique_reads),
                _fmt(row.singleton_reads),
                row.pct_singletons,
                _fmt(row.median_count),
                _fmt(row.max_count),
                row.dominant_size_by_reads_nt,
            )
        r.rep("")
    if t1b is not None:
        r.rep("### 2.2 Read-size composition (Table_01b)\n")
        r.rep("| sample | reads 20-22 nt | %% of all reads | dominant size |")
        r.rep("|---|---|---|---|")
        for s in samples:
            x = t1b[t1b["sample"] == s]
            if x.empty:
                continue
            reads_2022 = int(x.loc[x["width"].isin([20, 21, 22]), "n_reads"].sum())
            f = _pct(reads_2022 / x["n_reads"].sum())
            dom = int(x.loc[x["n_reads"].idxmax(), "width"])
            r.rep("| %s | %s | %.2f | %d nt |", s, _fmt(reads_2022), f, dom)
        r.rep("")

    for st in strategies:
        sdir = strat_path(st)
        t2a = _read_tab(sdir / "tables" / "Table_02a_annotation_count_unique_reads.csv")
        t2b = _read_tab(sdir / "tables" / "Table_02b_annotation_count_all_reads.csv")
        if t2b is not None:
            _comp_md(r, t2b, f"2.3a Annotation composition, all reads (strategy `{st.value}`)")
        if t2a is not None:
            _comp_md(r, t2a, f"2.3b Annotation composition, unique reads (strategy `{st.value}`)")

    for st in strategies:
        sdir = strat_path(st)
        t2b = _read_tab(sdir / "tables" / "Table_02b_annotation_count_all_reads.csv")
        if t2b is None:
            continue
        r.rep("### 2.4 Mature-miRNA read sizes (20-22 nt) - strategy `%s`\n", st.value)
        r.rep("| sample | matmiRNA reads | matmiRNA reads 20-22 nt | %% |")
        r.rep("|---|---|---|---|")
        for s in samples:
            mt = t2b[(t2b["category"] == "matmiRNA.size") & (t2b["sample"] == s)]
            if mt.empty:
                continue
            mt = mt.copy()
            mt["w"] = mt["item"].astype(int)
            tot = int(mt["Freq"].sum())
            w2022 = int(mt.loc[mt["w"].isin([20, 21, 22]), "Freq"].sum())
            r.rep("| %s | %s | %s | %.2f |", s, _fmt(tot), _fmt(w2022), _pct(w2022 / tot))
        r.rep("")

    for st in strategies:
        sdir = strat_path(st)
        t2m_rows: list[pd.DataFrame] = []
        for s in samples:
            per = _read_tab(sdir / "tables" / f"Table_02_{s}_unique_reads_annotation.csv")
            if per is None:
                continue
            mat = per[per["category"] == "matmiRNA"]
            if mat.empty:
                continue
            mat = mat.copy()
            mat["name"] = mat["feature_id"].astype(str)
            agg = (
                mat.groupby(["name"], sort=False)
                .agg(n_unique=("count", "size"), n_reads=("count", "sum"))
                .reset_index()
            )
            agg["sample"] = s
            t2m_rows.append(agg)
        if not t2m_rows:
            continue
        t2m = pd.concat(t2m_rows, ignore_index=True)
        r.rep("### 2.5 Top expressed mature miRNAs - strategy `%s`\n", st.value)
        for s in samples:
            x = t2m[t2m["sample"] == s].sort_values("n_reads", ascending=False).head(10)
            r.rep("**%s**\n", s)
            r.rep("| rank | mature miRNA | n_unique | n_reads |")
            r.rep("|---|---|---|---|")
            for i, row in enumerate(x.itertuples(index=False), start=1):
                r.rep("| %d | %s | %s | %s |", i, row.name, _fmt(row.n_unique), _fmt(row.n_reads))
            r.rep("")

    for st in strategies:
        sdir = strat_path(st)
        t5a = _read_tab(sdir / "tables" / "Table_03a_category_abundance_summary.csv")
        if t5a is None:
            continue
        r.rep("### 2.6 Per-category abundance - strategy `%s` (Table_03a)\n", st.value)
        r.rep(
            "| category | sample | n_loci | total reads | Gini | top1 %% | top5 %% | top10 %% | "
            "top locus |"
        )
        r.rep("|---|---|---|---|---|---|---|---|---|")
        for row in t5a.itertuples(index=False):
            r.rep(
                "| %s | %s | %s | %s | %.3f | %.1f | %.1f | %.1f | %s |",
                row.category,
                row.sample,
                _fmt(row.n_loci),
                _fmt(row.total_reads),
                row.gini,
                row.top1_pct,
                row.top5_pct,
                row.top10_pct,
                row.top_locus,
            )
        r.rep("")

    for st in strategies:
        sdir = strat_path(st)
        r.rep("### 2.7 Per-unique-read annotation + abundance table - strategy `%s`\n", st.value)
        r.rep("| sample | unique reads | total reads | categories | top feature (reads) |")
        r.rep("|---|---|---|---|---|")
        for s in samples:
            per = _read_tab(sdir / "tables" / f"Table_02_{s}_unique_reads_annotation.csv")
            if per is None:
                continue
            top = per.loc[per["count"].idxmax()]
            r.rep(
                "| %s | %s | %s | %d | %s | %s |",
                s,
                _fmt(len(per)),
                _fmt(int(per["count"].sum())),
                per["category"].nunique(),
                top["feature_id"],
                _fmt(top["count"]),
            )
        r.rep("")


def _comp_md(r: _Report, t2: pd.DataFrame, tag: str) -> None:
    x = t2[t2["category"] == "read.annotation"].copy()
    if x.empty:
        return
    x["pct"] = round(100 * x["Freq"] / x.groupby("sample")["Freq"].transform("sum"), 2)
    x = x.sort_values(["sample", "Freq"], ascending=[True, False])
    r.rep("### %s\n", tag)
    r.rep("| sample | category | reads | %% of reads |")
    r.rep("|---|---|---|---|")
    for row in x.itertuples(index=False):
        r.rep("| %s | %s | %s | %.2f |", row.sample, row.item, _fmt(row.Freq), row.pct)
    r.rep("")


def _sanity_checks(
    r: _Report,
    resolver: PathResolver,
    out_base: Path,
    samples: list[str],
    strategies: list[Strategy],
    strat_path: callable,
    store: FeatureStore,
) -> None:
    t1a = _read_tab(out_base / "tables" / "Table_01a_sample_summary.csv")
    t1b = _read_tab(out_base / "tables" / "Table_01b_read_size_distribution.csv")
    total_reads = (
        dict(zip(t1a["sample"], t1a["total_reads"], strict=False)) if t1a is not None else {}
    )
    uniq_reads = (
        dict(zip(t1a["sample"], t1a["unique_reads"], strict=False)) if t1a is not None else {}
    )

    t2b_list: dict[str, pd.DataFrame] = {}
    t2a_list: dict[str, pd.DataFrame] = {}
    for st in strategies:
        sdir = strat_path(st)
        t2b_list[st.value] = _read_tab(sdir / "tables" / "Table_02b_annotation_count_all_reads.csv")
        t2a_list[st.value] = _read_tab(
            sdir / "tables" / "Table_02a_annotation_count_unique_reads.csv"
        )

    for st in strategies:
        tag = "" if not resolver.config.strategy.is_comparison else st.value + "."
        t2a = t2a_list[st.value]
        t2b = t2b_list[st.value]
        sdir = strat_path(st)
        for s in samples:
            if t2b is not None:
                v = t2b[(t2b["category"] == "read.annotation") & (t2b["sample"] == s)]["Freq"].sum()
                ok = s in total_reads and v == total_reads[s]
                r.check(
                    f"s1_{tag}{s}",
                    f"sum of all category counts equals total reads ({s}, strategy {st.value})",
                    ok,
                    f"observed {_fmt(v)}, expected {_fmt(total_reads.get(s))}",
                )
            if t2a is not None:
                v = t2a[(t2a["category"] == "read.annotation") & (t2a["sample"] == s)]["Freq"].sum()
                ok = s in uniq_reads and v == uniq_reads[s]
                r.check(
                    f"s2_{tag}{s}",
                    f"sum of unique-read category counts equals unique reads "
                    f"({s}, strategy {st.value})",
                    ok,
                    f"observed {_fmt(v)}, expected {_fmt(uniq_reads.get(s))}",
                )
            csv = sdir / "tables" / f"Table_02_{s}_unique_reads_annotation.csv"
            per = _read_tab(csv)
            if per is not None:
                ok = s in uniq_reads and len(per) == uniq_reads[s]
                r.check(
                    f"s3_{tag}{s}",
                    f"unique-reads table row count equals unique reads ({s}, strategy {st.value})",
                    ok,
                    f"observed {_fmt(len(per))}, expected {_fmt(uniq_reads.get(s))}",
                )
                ok = s in total_reads and int(per["count"].sum()) == total_reads[s]
                r.check(
                    f"s4_{tag}{s}",
                    f"unique-reads table sum(count) equals total reads "
                    f"({s}, strategy {st.value})",
                    ok,
                    f"observed {_fmt(int(per['count'].sum()))}, "
                    f"expected {_fmt(total_reads.get(s))}",
                )
                if t2b is not None:
                    v = t2b[
                        (t2b["category"] == "read.annotation")
                        & (t2b["item"] == "matmiRNA")
                        & (t2b["sample"] == s)
                    ]["Freq"]
                    v = int(v.sum()) if len(v) else 0
                    mat_sum = int(per.loc[per["category"] == "matmiRNA", "count"].sum())
                    r.check(
                        f"s5_{tag}{s}",
                        f"Table_02b matmiRNA equals sum of matmiRNA reads "
                        f"in unique-reads table ({s}, strategy {st.value})",
                        v == mat_sum,
                        f"observed {_fmt(mat_sum)}, expected {_fmt(v)}",
                    )
                counts = per["count"].to_numpy()
                ok = counts[0] == counts.max() and bool((counts >= 1).all())
                r.check(
                    f"s6_{tag}{s}",
                    f"unique-reads table sorted by count descending, "
                    f"all counts >= 1 ({s}, strategy {st.value})",
                    ok,
                    "",
                )
                cpm_sum = float(per["cpm"].sum())
                r.check(
                    f"s7_{tag}{s}",
                    f"sum(cpm) ~ 1e6 in unique-reads table ({s}, strategy {st.value})",
                    abs(cpm_sum - 1e6) / 1e6 < 1e-3,
                    f"observed {_fmt(cpm_sum)}",
                )
            if t2b is not None:
                mt = t2b[(t2b["category"] == "matmiRNA.size") & (t2b["sample"] == s)]
                if not mt.empty:
                    mt = mt.copy()
                    mt["w"] = mt["item"].astype(int)
                    f = mt.loc[mt["w"].isin([20, 21, 22]), "Freq"].sum() / mt["Freq"].sum()
                    thresh = 0.90 if st.value == "fully-contained" else 0.85
                    r.check(
                        f"s9_{tag}{s}",
                        f"matmiRNA reads predominantly 20-22 nt (>= {thresh}) "
                        f"({s}, strategy {st.value})",
                        f >= thresh,
                        f"observed {f:.4f}",
                    )
            if t2a is not None and t2b is not None:
                a = t2a[(t2a["category"] == "read.annotation") & (t2a["sample"] == s)]
                b = t2b[(t2b["category"] == "read.annotation") & (t2b["sample"] == s)]
                m = a.merge(b, on="item", suffixes=(".u", ".r"))
                ok = bool((m["Freq.u"] <= m["Freq.r"]).all())
                r.check(
                    f"s10_{tag}{s}",
                    f"unique-read count <= all-read count for every category "
                    f"({s}, strategy {st.value})",
                    ok,
                    "",
                )
        t5a = _read_tab(sdir / "tables" / "Table_03a_category_abundance_summary.csv")
        if t5a is not None and t2b is not None:
            for s in samples:
                v2_row = t2b[
                    (t2b["category"] == "read.annotation")
                    & (t2b["item"] == "matmiRNA")
                    & (t2b["sample"] == s)
                ]["Freq"]
                v2 = int(v2_row.sum()) if len(v2_row) else 0
                v5_row = t5a[(t5a["category"] == "matmiRNA") & (t5a["sample"] == s)]["total_reads"]
                per = _read_tab(sdir / "tables" / f"Table_02_{s}_unique_reads_annotation.csv")
                exp = v2
                if per is not None:
                    mat = per[per["category"] == "matmiRNA"]
                    exp = (
                        int((mat["count"] * mat["n_features"].fillna(1).clip(lower=1)).sum())
                        if not mat.empty
                        else 0
                    )
                if len(v5_row) == 1:
                    v5 = int(v5_row.iloc[0])
                    r.check(
                        f"s13_{tag}{s}",
                        f"Table_03a matmiRNA total equals locus-assignment sum "
                        f"from unique-reads table ({s}, strategy {st.value})",
                        v5 == exp,
                        f"observed {_fmt(v5)}, expected {_fmt(exp)} "
                        f"(Table_02b matmiRNA: {_fmt(v2)})",
                    )

    if t1b is not None:
        for s in samples:
            x = t1b[t1b["sample"] == s]
            ok = bool(x["n_reads"].sum() == total_reads.get(s)) and bool(
                x["n_unique"].sum() == uniq_reads.get(s)
            )
            r.check(f"s8_{s}", f"Table_01b read sizes sum to total/unique reads ({s})", ok, "")

    n_mat = len(store.table("matmiRNA"))
    r.check(
        "s15_db",
        f"mature miRNA loci in DB == {_EXPECTED_MATMIRNA} (mm39)",
        n_mat == _EXPECTED_MATMIRNA,
        f"observed {n_mat}",
    )

    log_dir = resolver.log_dir
    logs = sorted(log_dir.glob("*.log")) if log_dir.exists() else []
    if not logs:
        r.check("s16_logs", "all step logs present", False, f"no step logs found in {log_dir}")
    else:
        r.check("s16_logs", "all step logs present", True, f"{len(logs)} logs found")
        bad = []
        for log in logs:
            text = log.read_text(encoding="utf-8", errors="ignore")
            if "ERROR" in text or "Traceback" in text:
                bad.append(log.name)
        r.check(
            "s17_logs",
            "no fatal errors in step logs",
            not bad,
            f"errors in: {', '.join(bad)}" if bad else "",
        )
