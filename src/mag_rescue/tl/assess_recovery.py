"""Compare ARIBA per-cohort recovery to Bacotype's complete-vs-SR penetrance.

Reads every ``*.report.tsv`` under ``<run-dir>/reports/``, tallies per-locus
presence and per-locus gene completeness, then — if Bacotype's
``complete_vs_sr_genomes/penetrance/<cohort>.csv`` is available — produces
the comparison tables that answer:

1. Did ARIBA find more than Bacotype's existing SR detection? (rescue check)
2. How does ARIBA's SR rate compare to the complete-genome ground truth?
3. When ARIBA detects a locus, is the full gene set there? (call quality)

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
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger("assess_recovery")


# --------------------------------------------------------------------------
# Locus registry: maps our metadata.tsv free_text labels to Bacotype features.
#
# Each entry: (locus_name, n_total_genes, bacotype_feature_name, gene_match_regex)
# gene_match_regex matches against the report's free_text column (last col).
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


def _free_text_labels(report_path: Path) -> set[str]:
    """Return the unique ``free_text`` values from one report (last col)."""
    labels: set[str] = set()
    with report_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols:
                labels.add(cols[-1])
    return labels


def _per_sample_locus_genes(reports_dir: Path, loci: list[LocusSpec]) -> dict[str, dict[str, set[str]]]:
    """Return ``{accession: {locus_name: {gene_label, ...}}}``.

    A "gene_label" is a ``<locus>:<gene>`` free-text entry. Empty inner sets
    mean the locus wasn't detected for that sample.
    """
    out: dict[str, dict[str, set[str]]] = {}
    for report in sorted(reports_dir.glob("*.report.tsv")):
        acc = report.name.removesuffix(".report.tsv")
        labels = _free_text_labels(report)
        per_locus = {l.name: {lab for lab in labels if l.gene_match_re.match(lab)} for l in loci}
        out[acc] = per_locus
    return out


# --------------------------------------------------------------------------
# Reading Bacotype's penetrance CSV
# --------------------------------------------------------------------------


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
    aribasr_n: int  # samples with >=1 gene of locus
    aribasr_pct: float
    aribasr_full_n: int  # samples with all genes (n_total_genes)
    aribasr_full_pct: float  # of detected, fraction with all genes
    aribasr_completeness_hist: list[tuple[int, int]]  # [(n_genes, n_samples), ...]
    bacotype_sr_pct: float | None = None  # penetrance as fraction (0..1)
    bacotype_complete_pct: float | None = None
    bacotype_n_complete: int | None = None
    bacotype_n_sr: int | None = None


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
    per_sample = _per_sample_locus_genes(reports_dir, loci)
    n_total = len(per_sample)
    if n_total == 0:
        logger.error("no *.report.tsv in %s", reports_dir)
        sys.exit(1)

    bacotype = _read_bacotype_penetrance(bacotype_penetrance_csv) if bacotype_penetrance_csv else {}

    rows: list[LocusRow] = []
    for spec in loci:
        gene_counts = [len(per_sample[acc][spec.name]) for acc in per_sample if per_sample[acc][spec.name]]
        n_present = len(gene_counts)
        n_full = sum(1 for c in gene_counts if c == spec.n_total_genes)
        hist = sorted(Counter(gene_counts).items())

        b = bacotype.get(spec.bacotype_feature, {})
        rows.append(
            LocusRow(
                name=spec.name,
                aribasr_n=n_present,
                aribasr_pct=n_present / n_total,
                aribasr_full_n=n_full,
                aribasr_full_pct=(n_full / n_present) if n_present > 0 else 0.0,
                aribasr_completeness_hist=hist,
                bacotype_sr_pct=float(b["sr_penetrance"]) if "sr_penetrance" in b else None,
                bacotype_complete_pct=float(b["complete_penetrance"]) if "complete_penetrance" in b else None,
                bacotype_n_complete=int(b["n_complete"]) if "n_complete" in b else None,
                bacotype_n_sr=int(b["n_sr"]) if "n_sr" in b else None,
            )
        )
    return n_total, rows


# --------------------------------------------------------------------------
# Render
# --------------------------------------------------------------------------


def _fmt_pct(x: float | None) -> str:
    return "n/a" if x is None else f"{x:.2%}"


def render_markdown(*, cohort: str, n_total: int, rows: list[LocusRow]) -> str:
    """Build the human-readable markdown summary."""
    has_bact = any(r.bacotype_sr_pct is not None for r in rows)
    out: list[str] = []
    out.append(f"# {cohort} — ARIBA recovery assessment\n")
    out.append(f"**ARIBA samples processed:** {n_total}\n")
    if has_bact:
        any_n_sr = next((r.bacotype_n_sr for r in rows if r.bacotype_n_sr is not None), None)
        any_n_c = next((r.bacotype_n_complete for r in rows if r.bacotype_n_complete is not None), None)
        if any_n_sr or any_n_c:
            out.append(f"**Bacotype reference:** n_complete={any_n_c}, n_sr={any_n_sr}\n")

    out.append("\n## Per-locus presence (% of samples with ≥1 gene)\n")
    if has_bact:
        out.append("| Locus | Complete | Bacotype SR | ARIBA SR | Δ (ARIBA − Bacotype SR) |")
        out.append("|---|---|---|---|---|")
        for r in rows:
            delta = ""
            if r.bacotype_sr_pct is not None:
                d = r.aribasr_pct - r.bacotype_sr_pct
                ratio = (r.aribasr_pct / r.bacotype_sr_pct) if r.bacotype_sr_pct > 0 else float("inf")
                delta = f"{d:+.2%} ({ratio:.1f}×)" if r.bacotype_sr_pct > 0 else f"{d:+.2%} (∞)"
            out.append(
                f"| {r.name} "
                f"| {_fmt_pct(r.bacotype_complete_pct)} ({r.bacotype_n_complete or '?'}) "
                f"| {_fmt_pct(r.bacotype_sr_pct)} "
                f"| **{_fmt_pct(r.aribasr_pct)}** ({r.aribasr_n}/{n_total}) "
                f"| {delta} |"
            )
    else:
        out.append("| Locus | ARIBA SR | n |")
        out.append("|---|---|---|")
        for r in rows:
            out.append(f"| {r.name} | **{_fmt_pct(r.aribasr_pct)}** | {r.aribasr_n}/{n_total} |")

    out.append("\n## Call consistency (when locus is detected, how many of its genes are present?)\n")
    out.append("| Locus | Total genes | Full locus | Partial breakdown (n_genes → n_samples) |")
    out.append("|---|---|---|---|")
    for r in rows:
        if r.aribasr_n == 0:
            out.append(f"| {r.name} | {DB_LOCI_NTOTAL[r.name]} | — | not detected |")
            continue
        partial = ", ".join(f"{c} → {n}" for c, n in r.aribasr_completeness_hist if c != DB_LOCI_NTOTAL[r.name])
        partial = partial or "(none partial)"
        out.append(
            f"| {r.name} | {DB_LOCI_NTOTAL[r.name]} "
            f"| {r.aribasr_full_n}/{r.aribasr_n} ({r.aribasr_full_pct:.1%}) "
            f"| {partial} |"
        )
    return "\n".join(out) + "\n"


# Built lazily after the dataclass — avoids repeating the gene counts.
DB_LOCI_NTOTAL: dict[str, int] = {}


def render_tsv(*, cohort: str, n_total: int, rows: list[LocusRow]) -> str:
    """Machine-readable TSV with all the same data plus completeness histogram as JSON-ish."""
    header = [
        "cohort",
        "n_aribasr",
        "locus",
        "n_total_genes",
        "ariba_sr_n",
        "ariba_sr_pct",
        "ariba_sr_full_n",
        "ariba_sr_full_frac_of_detected",
        "ariba_sr_completeness_hist",
        "bacotype_sr_pct",
        "bacotype_complete_pct",
        "bacotype_n_sr",
        "bacotype_n_complete",
    ]
    lines = ["\t".join(header)]
    for r in rows:
        hist = ";".join(f"{c}:{n}" for c, n in r.aribasr_completeness_hist)
        lines.append(
            "\t".join(
                str(x)
                for x in [
                    cohort,
                    n_total,
                    r.name,
                    DB_LOCI_NTOTAL[r.name],
                    r.aribasr_n,
                    f"{r.aribasr_pct:.6f}",
                    r.aribasr_full_n,
                    f"{r.aribasr_full_pct:.6f}",
                    hist,
                    "" if r.bacotype_sr_pct is None else f"{r.bacotype_sr_pct:.6f}",
                    "" if r.bacotype_complete_pct is None else f"{r.bacotype_complete_pct:.6f}",
                    "" if r.bacotype_n_sr is None else r.bacotype_n_sr,
                    "" if r.bacotype_n_complete is None else r.bacotype_n_complete,
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

    # Populate DB_LOCI_NTOTAL once per run.
    DB_LOCI_NTOTAL.clear()
    DB_LOCI_NTOTAL.update({l.name: l.n_total_genes for l in DB_LOCI[args.db]})

    bact_csv = (args.bacotype_dir / "penetrance" / f"{args.cohort}.csv") if args.bacotype_dir else None
    if bact_csv and not bact_csv.is_file():
        logger.warning("bacotype penetrance file missing: %s — emitting ARIBA-only tables", bact_csv)
        bact_csv = None

    n_total, rows = assess(ariba_run_dir=args.ariba_run_dir, db=args.db, bacotype_penetrance_csv=bact_csv)
    md = render_markdown(cohort=args.cohort, n_total=n_total, rows=rows)
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
