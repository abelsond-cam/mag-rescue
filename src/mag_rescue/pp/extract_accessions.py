"""One-shot accession-list builder for the ARIBA Slurm array.

Reads Bacotype's curated metadata TSV, filters to the kleb short-read subset
(``kpsc_final_list==True AND is_refseq==False``), parses ``fastq_ftp`` URLs,
and writes a slim 5-column TSV the array job consumes line-by-line. Every
row that does NOT make it into the inclusion list is written to a sidecar
skipped TSV with a ``reason`` column, so deferred edge cases stay visible.

Usage
-----
    pixi run -e dev python -m mag_rescue.pp.extract_accessions \
        --metadata <bacotype_metadata.tsv> \
        --outdir <run-accessions-dir> \
        --version v1
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger("extract_accessions")

# Column names we read from Bacotype's metadata.tsv.
COL_KPSC = "kpsc_final_list"
COL_REFSEQ = "is_refseq"
COL_RUN_ACC = "run_accession"
COL_FASTQ_FTP = "fastq_ftp"
COL_FASTQ_MD5 = "fastq_md5"
COL_PLATFORM = "metadata.runs.instrument.platform"
COL_SUBLINEAGE = "Sublineage"
COL_CLONAL_GROUP = "Clonal group"

OUTPUT_HEADER = ["run_accession", "r1_url", "r2_url", "r1_md5", "r2_md5"]
SKIPPED_HEADER = [
    "run_accession",
    "reason",
    "platform",
    "fastq_ftp",
    "fastq_md5",
]


def _https_url(ena_path: str) -> str:
    """Convert an ENA bare-host path to an HTTPS URL.

    Bacotype's ``fastq_ftp`` column omits the scheme: e.g.
    ``ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/SRR000001/SRR000001_1.fastq.gz``.
    ENA serves the same paths over HTTPS, which is faster and firewall-friendlier
    than passive FTP from CSD3.
    """
    p = ena_path.strip()
    if p.startswith(("http://", "https://", "ftp://")):
        return p
    return f"https://{p}"


def _split_semicolon(field: str) -> list[str]:
    """Split a metadata field on ``;`` (ENA convention), dropping empties."""
    if not field:
        return []
    return [tok for tok in (t.strip() for t in field.split(";")) if tok]


def _pick_paired(urls: list[str], md5s: list[str]) -> tuple[str, str, str, str] | None:
    """Pick the (R1, R2) pair from ENA's url/md5 lists.

    Returns ``(r1_url, r2_url, r1_md5, r2_md5)`` on success, or ``None`` if no
    valid pair is found. Handles two cases:

    * 2 items → take in order (assumes ENA's standard order: R1, R2).
    * 3 items → keep the entries whose URLs end in ``_1.fastq.gz`` /
      ``_2.fastq.gz`` (drop the orphan single-end). Lists must align by index.

    For other counts (0, 1, 4+), returns ``None`` — caller logs to skipped.tsv.
    """
    if len(urls) != len(md5s):
        return None
    if len(urls) == 2:
        return _https_url(urls[0]), _https_url(urls[1]), md5s[0], md5s[1]
    if len(urls) == 3:
        idx_r1 = next((i for i, u in enumerate(urls) if u.endswith("_1.fastq.gz")), None)
        idx_r2 = next((i for i, u in enumerate(urls) if u.endswith("_2.fastq.gz")), None)
        if idx_r1 is None or idx_r2 is None:
            return None
        return _https_url(urls[idx_r1]), _https_url(urls[idx_r2]), md5s[idx_r1], md5s[idx_r2]
    return None


def _classify_row(row: dict[str, str]) -> tuple[list[str], str | None]:
    """Decide whether a row is includable; if not, return the reason.

    Returns ``(output_row, None)`` for include, or ``([...], reason)`` for skip.
    """
    acc = row.get(COL_RUN_ACC, "").strip()
    platform = row.get(COL_PLATFORM, "").strip()
    fastq_ftp = row.get(COL_FASTQ_FTP, "").strip()
    fastq_md5 = row.get(COL_FASTQ_MD5, "").strip()

    if not acc:
        return [], "no_run_accession"
    if not fastq_ftp:
        return [], "no_fastq_ftp"
    if "ILLUMINA" not in platform:
        return [], "non_illumina"
    if platform.strip() != "ILLUMINA":
        # Multi-run rows like "ILLUMINA || OXFORD_NANOPORE" — defer per plan.
        return [], "multi_platform_or_multi_run"

    urls = _split_semicolon(fastq_ftp)
    md5s = _split_semicolon(fastq_md5)
    pair = _pick_paired(urls, md5s)
    if pair is None:
        return [], f"unsupported_url_count_{len(urls)}"

    r1_url, r2_url, r1_md5, r2_md5 = pair
    return [acc, r1_url, r2_url, r1_md5, r2_md5], None


def _filter_to_kleb_short_reads(
    metadata_path: Path,
    *,
    sublineage: str | None = None,
    clonal_group: str | None = None,
) -> tuple[list[list[str]], list[list[str]]]:
    """Stream the metadata TSV and partition into included + skipped rows.

    Optional ``sublineage`` and ``clonal_group`` further restrict the inclusion
    set. Rows that pass kpsc/refseq but fail the SL/CG filter are simply not
    written anywhere — they're not the same kind of "skip" as a malformed row.
    """
    required = (COL_KPSC, COL_REFSEQ, COL_RUN_ACC, COL_FASTQ_FTP, COL_FASTQ_MD5, COL_PLATFORM)
    if sublineage is not None:
        required = (*required, COL_SUBLINEAGE)
    if clonal_group is not None:
        required = (*required, COL_CLONAL_GROUP)

    included: list[list[str]] = []
    skipped: list[list[str]] = []
    with metadata_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for col in required:
            if col not in reader.fieldnames:
                raise SystemExit(f"metadata TSV missing required column: {col}")
        for row in reader:
            if row.get(COL_KPSC) != "True" or row.get(COL_REFSEQ) != "False":
                continue
            if sublineage is not None and row.get(COL_SUBLINEAGE, "").strip() != sublineage:
                continue
            if clonal_group is not None and row.get(COL_CLONAL_GROUP, "").strip() != clonal_group:
                continue
            out_row, reason = _classify_row(row)
            if reason is None:
                included.append(out_row)
            else:
                skipped.append(
                    [
                        row.get(COL_RUN_ACC, "").strip(),
                        reason,
                        row.get(COL_PLATFORM, "").strip(),
                        row.get(COL_FASTQ_FTP, "").strip(),
                        row.get(COL_FASTQ_MD5, "").strip(),
                    ]
                )
    return included, skipped


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _write_summary(path: Path, included: list[list[str]], skipped: list[list[str]], metadata_path: Path) -> None:
    counts = Counter(row[1] for row in skipped)
    lines = [
        f"source: {metadata_path}",
        f"included: {len(included)}",
        f"skipped:  {len(skipped)}",
        "skip reasons:",
    ]
    for reason, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        lines.append(f"  {reason:32s} {n}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    """CLI entrypoint."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--metadata", type=Path, required=True, help="Bacotype's curated metadata TSV.")
    ap.add_argument("--outdir", type=Path, required=True, help="Where to write the three output files.")
    ap.add_argument("--version", default="v1", help="Filename suffix tag (default: v1).")
    ap.add_argument("--sublineage", default=None, help="Optional: restrict to one Sublineage (e.g. SL23).")
    ap.add_argument("--clonal-group", default=None, help="Optional: restrict to one Clonal group (e.g. CG39).")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")

    if not args.metadata.is_file():
        logger.error("metadata not found: %s", args.metadata)
        sys.exit(1)

    logger.info(
        "Reading %s%s%s",
        args.metadata,
        f"  [Sublineage={args.sublineage}]" if args.sublineage else "",
        f"  [Clonal group={args.clonal_group}]" if args.clonal_group else "",
    )
    included, skipped = _filter_to_kleb_short_reads(
        args.metadata,
        sublineage=args.sublineage,
        clonal_group=args.clonal_group,
    )

    base = f"kleb_short_reads_{args.version}"
    out_tsv = args.outdir / f"{base}.tsv"
    skipped_tsv = args.outdir / f"{base}.skipped.tsv"
    summary_txt = args.outdir / f"{base}.summary.txt"

    _write_tsv(out_tsv, OUTPUT_HEADER, included)
    _write_tsv(skipped_tsv, SKIPPED_HEADER, skipped)
    _write_summary(summary_txt, included, skipped, args.metadata)

    logger.info("Wrote %s (%d rows)", out_tsv, len(included))
    logger.info("Wrote %s (%d rows)", skipped_tsv, len(skipped))
    logger.info("Wrote %s", summary_txt)


if __name__ == "__main__":
    main()
