"""Compare ARIBA per-cohort recovery to Bacotype's complete-vs-SR penetrance.

Reads every ``*.report.tsv`` under ``<run-dir>/reports/``, tallies per-locus
presence, gene-level concordance, and recovery quality (kb assembled + mean
contig depth), then — if Bacotype's ``complete_vs_sr_genomes/penetrance/
<cohort>.csv`` is available — joins to produce the comparison tables that
answer:

1. Did ARIBA find more than Bacotype's existing short-read MAG detection?
2. How does ARIBA's rate compare to the complete-genome ground truth?
3. When ARIBA detects a locus, is the full gene set there and is the
   assembly deep enough that we should trust the call?

Headers use **ARIBA** (our reads-based call), **sr-MAG** (Bacotype's
assembly-based call on the same short reads), and **Complete** (Bacotype
ground truth from complete genomes).

The Bacotype penetrance CSV is keyed by clonal group. For sublineages
(SL<N>), we fall back to CG<N>.csv when present — SL/CG pairs are often
near-equivalent cohorts.

The "mean depth" column uses ARIBA's own ``ctg_cov`` value (mean per-base
depth on the local assembly) — that's exactly the
``reads × read_length / ctg_len`` quantity, computed from the actual BAM,
so we don't re-estimate it.

Usage
-----
    pixi run -e dev python -m mag_rescue.tl.assess_recovery \
        --ariba-run-dir <RDS>/.../mag_rescue/kleb_virulence/CG39 \
        --cohort CG39 \
        --bacotype-dir <RDS>/.../complete_vs_sr_genomes
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger("assess_recovery")


# --------------------------------------------------------------------------
# Locus registry
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class LocusSpec:
    """One locus's identity for tally + Bacotype-feature lookup."""

    name: str
    n_total_genes: int
    bacotype_feature: str
    gene_match_re: re.Pattern[str]


DB_LOCI: dict[str, list[LocusSpec]] = {
    "kleb_virulence": [
        LocusSpec("ybt", 11, "Yersiniabactin_bsc", re.compile(r"^ybt:")),
        LocusSpec("iuc", 5, "Aerobactin_bsc", re.compile(r"^iuc:")),
        LocusSpec("iro", 4, "Salmochelin_bsc", re.compile(r"^iro:")),
        LocusSpec("clb", 15, "Colibactin_bsc", re.compile(r"^clb:")),
        # Bacotype tracks the rmp ADC operon and the rmpA2 paralog separately.
        # Our metadata clusters both under "rmp:..." so we split by gene name.
        LocusSpec("rmpADC", 3, "RmpADC_bsc", re.compile(r"^rmp:rmp[ACD]$")),
        LocusSpec("rmpa2", 1, "rmpA2_bsc", re.compile(r"^rmp:rmpA2$")),
    ],
}


# --------------------------------------------------------------------------
# Reading ARIBA reports
# --------------------------------------------------------------------------


# ARIBA report.tsv column positions (1-indexed in docs; here 0-indexed).
_COL_READS = 5
_COL_REF_LEN = 7
_COL_REF_BASE_ASSEMBLED = 8
_COL_PC_IDENT = 9
_COL_CTG_LEN = 11
_COL_CTG_COV = 12
_COL_FREE_TEXT = 30


@dataclass(frozen=True)
class ReportRow:
    """One row from an ARIBA report.tsv (only the columns we use)."""

    free_text: str
    ref_len: int
    ref_base_assembled: int
    pc_ident: float
    ctg_len: int
    ctg_cov: float
    reads: int


def _safe_int(s: str) -> int | None:
    try:
        return int(s)
    except ValueError:
        return None


def _safe_float(s: str) -> float | None:
    try:
        return float(s)
    except ValueError:
        return None


def _read_report(report_path: Path) -> list[ReportRow]:
    """Return typed rows from one ARIBA report. Header + bad rows are skipped."""
    rows: list[ReportRow] = []
    with report_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= _COL_FREE_TEXT:
                continue
            reads = _safe_int(cols[_COL_READS])
            ref_len = _safe_int(cols[_COL_REF_LEN])
            ref_base_assembled = _safe_int(cols[_COL_REF_BASE_ASSEMBLED])
            pc_ident = _safe_float(cols[_COL_PC_IDENT])
            ctg_len = _safe_int(cols[_COL_CTG_LEN])
            ctg_cov = _safe_float(cols[_COL_CTG_COV])
            if None in (reads, ref_len, ref_base_assembled, pc_ident, ctg_len, ctg_cov):
                continue
            rows.append(
                ReportRow(
                    free_text=cols[_COL_FREE_TEXT],
                    ref_len=ref_len,
                    ref_base_assembled=ref_base_assembled,
                    pc_ident=pc_ident,
                    ctg_len=ctg_len,
                    ctg_cov=ctg_cov,
                    reads=reads,
                )
            )
    return rows


def _per_sample_locus_rows(reports_dir: Path, loci: list[LocusSpec]) -> dict[str, dict[str, list[ReportRow]]]:
    """Return ``{accession: {locus_name: [ReportRow, ...]}}``."""
    out: dict[str, dict[str, list[ReportRow]]] = {}
    for report in sorted(reports_dir.glob("*.report.tsv")):
        acc = report.name.removesuffix(".report.tsv")
        rows = _read_report(report)
        by_locus: dict[str, list[ReportRow]] = {spec.name: [] for spec in loci}
        for row in rows:
            for spec in loci:
                if spec.gene_match_re.match(row.free_text):
                    by_locus[spec.name].append(row)
                    break
        out[acc] = by_locus
    return out


# --------------------------------------------------------------------------
# Reading Bacotype's penetrance CSV
# --------------------------------------------------------------------------


def _find_penetrance_csv(bacotype_dir: Path, cohort: str) -> Path | None:
    """Find Bacotype's penetrance/<cohort>.csv with an SL→CG fallback.

    Bacotype's complete_vs_sr_genomes/penetrance/ is keyed by clonal group.
    For sublineages (SL<N>) we fall back to penetrance/CG<N>.csv when present
    — they're not identical cohorts but for many SLs (SL23/CG23, SL39/CG39,
    etc.) the dominant CG covers the SL.
    """
    direct = bacotype_dir / "penetrance" / f"{cohort}.csv"
    if direct.is_file():
        return direct
    if cohort.startswith("SL"):
        alt = bacotype_dir / "penetrance" / f"CG{cohort[2:]}.csv"
        if alt.is_file():
            logger.warning("no penetrance/%s.csv; falling back to %s (SL→CG substitution)", cohort, alt.name)
            return alt
    return None


def _read_bacotype_penetrance(csv_path: Path) -> dict[str, dict]:
    """Parse Bacotype's per-cohort penetrance CSV keyed by feature name."""
    rows: dict[str, dict] = {}
    if not csv_path.is_file():
        return rows
    with csv_path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows[row["feature"]] = row
    return rows


# --------------------------------------------------------------------------
# Tabulate
# --------------------------------------------------------------------------


@dataclass
class LocusRow:
    """Per-locus tally for one cohort, including the Bacotype comparison."""

    name: str
    n_total_genes: int
    # presence
    ariba_n: int  # samples with >=1 gene of locus
    ariba_pct: float
    # concordance
    ariba_full_n: int
    ariba_full_pct: float  # of detected, fraction with all genes
    ariba_half_or_more_n: int
    ariba_half_or_more_pct: float  # of detected, fraction with >= ceil(n_total/2) genes
    ariba_completeness_hist: list[tuple[int, int]]  # [(n_genes, n_samples), ...]
    # recovery quality (mean over positive samples)
    mean_kb_recovered: float
    mean_pct_ref_recovered: float
    mean_ctg_cov: float
    # Bacotype comparison
    srmag_pct: float | None = None
    complete_pct: float | None = None
    n_complete: int | None = None
    n_srmag: int | None = None


def _aggregate_locus(spec: LocusSpec, per_sample: dict[str, dict[str, list[ReportRow]]]) -> dict:
    """Compute per-locus aggregates across all samples (positives only for quality)."""
    n_total = len(per_sample)
    half_threshold = math.ceil(spec.n_total_genes / 2)

    gene_counts: list[int] = []
    kb_recovered: list[float] = []
    pct_recovered: list[float] = []
    ctg_covs: list[float] = []

    for acc in per_sample:
        rows = per_sample[acc][spec.name]
        if not rows:
            continue
        unique_genes = {r.free_text for r in rows}
        gene_counts.append(len(unique_genes))
        sum_assembled = sum(r.ref_base_assembled for r in rows)
        sum_ref_len = sum(r.ref_len for r in rows)
        kb_recovered.append(sum_assembled / 1000.0)
        if sum_ref_len > 0:
            pct_recovered.append(sum_assembled / sum_ref_len)
        non_zero_cov = [r.ctg_cov for r in rows if r.ctg_cov > 0]
        if non_zero_cov:
            ctg_covs.append(sum(non_zero_cov) / len(non_zero_cov))

    n_present = len(gene_counts)
    n_full = sum(1 for c in gene_counts if c == spec.n_total_genes)
    n_half_or_more = sum(1 for c in gene_counts if c >= half_threshold)
    hist = sorted(Counter(gene_counts).items())

    return {
        "n_total": n_total,
        "n_present": n_present,
        "n_full": n_full,
        "n_half_or_more": n_half_or_more,
        "completeness_hist": hist,
        "mean_kb_recovered": (sum(kb_recovered) / len(kb_recovered)) if kb_recovered else 0.0,
        "mean_pct_ref_recovered": (sum(pct_recovered) / len(pct_recovered)) if pct_recovered else 0.0,
        "mean_ctg_cov": (sum(ctg_covs) / len(ctg_covs)) if ctg_covs else 0.0,
    }


def assess(
    *,
    ariba_run_dir: Path,
    db: str,
    bacotype_penetrance_csv: Path | None,
) -> tuple[int, list[LocusRow]]:
    """Tally per-locus stats for one cohort. Returns (n_samples_total, rows)."""
    reports_dir = ariba_run_dir / "reports"
    if not reports_dir.is_dir():
        logger.error("reports dir missing: %s", reports_dir)
        sys.exit(1)
    loci = DB_LOCI[db]
    per_sample = _per_sample_locus_rows(reports_dir, loci)
    n_total = len(per_sample)
    if n_total == 0:
        logger.error("no *.report.tsv in %s", reports_dir)
        sys.exit(1)

    bacotype = _read_bacotype_penetrance(bacotype_penetrance_csv) if bacotype_penetrance_csv else {}

    rows: list[LocusRow] = []
    for spec in loci:
        agg = _aggregate_locus(spec, per_sample)
        b = bacotype.get(spec.bacotype_feature, {})
        rows.append(
            LocusRow(
                name=spec.name,
                n_total_genes=spec.n_total_genes,
                ariba_n=agg["n_present"],
                ariba_pct=agg["n_present"] / agg["n_total"],
                ariba_full_n=agg["n_full"],
                ariba_full_pct=(agg["n_full"] / agg["n_present"]) if agg["n_present"] > 0 else 0.0,
                ariba_half_or_more_n=agg["n_half_or_more"],
                ariba_half_or_more_pct=((agg["n_half_or_more"] / agg["n_present"]) if agg["n_present"] > 0 else 0.0),
                ariba_completeness_hist=agg["completeness_hist"],
                mean_kb_recovered=agg["mean_kb_recovered"],
                mean_pct_ref_recovered=agg["mean_pct_ref_recovered"],
                mean_ctg_cov=agg["mean_ctg_cov"],
                srmag_pct=float(b["sr_penetrance"]) if "sr_penetrance" in b else None,
                complete_pct=float(b["complete_penetrance"]) if "complete_penetrance" in b else None,
                n_complete=int(b["n_complete"]) if "n_complete" in b else None,
                n_srmag=int(b["n_sr"]) if "n_sr" in b else None,
            )
        )
    return n_total, rows


# --------------------------------------------------------------------------
# Render
# --------------------------------------------------------------------------


def _fmt_pct(x: float | None) -> str:
    return "n/a" if x is None else f"{x:.2%}"


def render_markdown(*, cohort: str, n_total: int, rows: list[LocusRow], penetrance_source: str | None = None) -> str:
    """Build the human-readable markdown summary."""
    has_bact = any(r.srmag_pct is not None for r in rows)
    out: list[str] = []
    out.append(f"# {cohort} — ARIBA recovery assessment\n")
    out.append(f"**ARIBA samples processed:** {n_total}\n")
    if has_bact:
        any_n_sr = next((r.n_srmag for r in rows if r.n_srmag is not None), None)
        any_n_c = next((r.n_complete for r in rows if r.n_complete is not None), None)
        if any_n_sr or any_n_c:
            src = f" (from {penetrance_source})" if penetrance_source else ""
            out.append(f"**Bacotype reference{src}:** n_complete={any_n_c}, n_sr-MAG={any_n_sr}\n")

    out.append("\n## Per-locus presence (% of samples with ≥1 gene)\n")
    if has_bact:
        out.append("| Locus | Complete | sr-MAG | ARIBA | Δ |")
        out.append("|---|---|---|---|---|")
        for r in rows:
            delta = ""
            if r.srmag_pct is not None:
                d = r.ariba_pct - r.srmag_pct
                if r.srmag_pct > 0:
                    ratio = r.ariba_pct / r.srmag_pct
                    delta = f"{d:+.2%} ({ratio:.1f}×)"
                else:
                    delta = f"{d:+.2%} (∞)"
            out.append(
                f"| {r.name} "
                f"| {_fmt_pct(r.complete_pct)} ({r.n_complete or '?'}) "
                f"| {_fmt_pct(r.srmag_pct)} "
                f"| **{_fmt_pct(r.ariba_pct)}** ({r.ariba_n}/{n_total}) "
                f"| {delta} |"
            )
    else:
        out.append("| Locus | ARIBA | n |")
        out.append("|---|---|---|")
        for r in rows:
            out.append(f"| {r.name} | **{_fmt_pct(r.ariba_pct)}** | {r.ariba_n}/{n_total} |")

    out.append("\n## Locus concordance & recovery quality (positives only)\n")
    out.append("| Locus | n positive | % full concordance | % ≥50% concordance | mean kb recovered | mean depth (×) |")
    out.append("|---|---|---|---|---|---|")
    for r in rows:
        if r.ariba_n == 0:
            out.append(f"| {r.name} | 0 | — | — | — | — |")
            continue
        out.append(
            f"| {r.name} "
            f"| {r.ariba_n} "
            f"| {r.ariba_full_pct:.1%} "
            f"| {r.ariba_half_or_more_pct:.1%} "
            f"| {r.mean_kb_recovered:.1f} "
            f"| {r.mean_ctg_cov:.0f} |"
        )
    return "\n".join(out) + "\n"


def render_tsv(*, cohort: str, n_total: int, rows: list[LocusRow]) -> str:
    """Machine-readable TSV with every metric, plus the completeness histogram."""
    header = [
        "cohort",
        "n_ariba",
        "locus",
        "n_total_genes",
        "ariba_n",
        "ariba_pct",
        "ariba_full_n",
        "ariba_full_frac",
        "ariba_half_or_more_n",
        "ariba_half_or_more_frac",
        "mean_kb_recovered",
        "mean_pct_ref_recovered",
        "mean_ctg_cov",
        "completeness_hist",
        "srmag_pct",
        "complete_pct",
        "n_srmag",
        "n_complete",
    ]
    lines = ["\t".join(header)]
    for r in rows:
        hist = ";".join(f"{c}:{n}" for c, n in r.ariba_completeness_hist)
        lines.append(
            "\t".join(
                str(x)
                for x in [
                    cohort,
                    n_total,
                    r.name,
                    r.n_total_genes,
                    r.ariba_n,
                    f"{r.ariba_pct:.6f}",
                    r.ariba_full_n,
                    f"{r.ariba_full_pct:.6f}",
                    r.ariba_half_or_more_n,
                    f"{r.ariba_half_or_more_pct:.6f}",
                    f"{r.mean_kb_recovered:.4f}",
                    f"{r.mean_pct_ref_recovered:.6f}",
                    f"{r.mean_ctg_cov:.2f}",
                    hist,
                    "" if r.srmag_pct is None else f"{r.srmag_pct:.6f}",
                    "" if r.complete_pct is None else f"{r.complete_pct:.6f}",
                    "" if r.n_srmag is None else r.n_srmag,
                    "" if r.n_complete is None else r.n_complete,
                ]
            )
        )
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main() -> None:
    """CLI entrypoint."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument(
        "--ariba-run-dir",
        type=Path,
        required=True,
        help="The cohort run dir containing reports/. e.g. <RDS>/.../mag_rescue/kleb_virulence/CG39",
    )
    ap.add_argument("--cohort", required=True, help="Cohort name (e.g. CG39, SL23). Used to find Bacotype CSVs.")
    ap.add_argument("--db", default="kleb_virulence", help="DB name in DB_LOCI registry.")
    ap.add_argument(
        "--bacotype-dir",
        type=Path,
        default=None,
        help="Path to Bacotype's complete_vs_sr_genomes/ dir. Comparison tables omitted if not given.",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Where to write <cohort>_recovery.{md,tsv}. Defaults to <ariba-run-dir>/assessment/.",
    )
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")

    bact_csv = _find_penetrance_csv(args.bacotype_dir, args.cohort) if args.bacotype_dir else None
    penetrance_source = bact_csv.name if bact_csv else None

    n_total, rows = assess(ariba_run_dir=args.ariba_run_dir, db=args.db, bacotype_penetrance_csv=bact_csv)
    md = render_markdown(cohort=args.cohort, n_total=n_total, rows=rows, penetrance_source=penetrance_source)
    tsv = render_tsv(cohort=args.cohort, n_total=n_total, rows=rows)

    out_dir = args.output_dir or (args.ariba_run_dir / "assessment")
    out_dir.mkdir(parents=True, exist_ok=True)
    md_path = out_dir / f"{args.cohort}_recovery.md"
    tsv_path = out_dir / f"{args.cohort}_recovery.tsv"
    md_path.write_text(md)
    tsv_path.write_text(tsv)

    logger.info("wrote %s", md_path)
    logger.info("wrote %s", tsv_path)
    sys.stdout.write(md)


if __name__ == "__main__":
    main()
