"""Per-sample ARIBA worker (one Slurm array task = one accession).

Idempotent: exits 0 if the accession's report.tsv already exists. Downloads
fastqs from ENA over HTTPS, verifies md5, runs ``ariba run`` against the
selected DB's prepareref output, copies a slim set of summary files to RDS,
deletes scratch.

Usage
-----
    pixi run -e dev python -m mag_rescue.pp.run_ariba \
        --db kleb_virulence \
        --accession SRR... --r1-url ... --r2-url ... --r1-md5 ... --r2-md5 ... \
        --workdir $SLURM_TMPDIR/<acc> \
        --run-dir <RDS>/.../mag_rescue/kleb_virulence/all \
        --threads $SLURM_CPUS_PER_TASK [--detailed]
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import shutil
import subprocess
import sys
from pathlib import Path

logger = logging.getLogger("run_ariba")

# Output filenames ARIBA produces in its outdir we care to keep.
ARIBA_REPORT = "report.tsv"
ARIBA_DETAILED_FILES = ("assembled_seqs.fa.gz", "assemblies.fa.gz")

# Repo root (..../mag_rescue/mag-rescue/), used to resolve the prepareref dir.
REPO_ROOT = Path(__file__).resolve().parents[3]


def _md5sum(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for buf in iter(lambda: fh.read(chunk), b""):
            h.update(buf)
    return h.hexdigest()


def _verify_md5(path: Path, expected: str) -> bool:
    got = _md5sum(path)
    if got != expected:
        logger.error("md5 mismatch for %s: expected %s, got %s", path, expected, got)
        return False
    return True


def _curl_download(url: str, dest: Path, retries: int = 5, retry_delay: int = 30) -> int:
    """Download ``url`` to ``dest`` via curl. Returns curl's exit code."""
    cmd = [
        "curl",
        "-fsSL",
        "--retry",
        str(retries),
        "--retry-delay",
        str(retry_delay),
        "-o",
        str(dest),
        url,
    ]
    logger.info("curl -> %s", dest.name)
    res = subprocess.run(cmd, check=False)
    return res.returncode


def _atomic_copy(src: Path, dest: Path) -> None:
    """Copy ``src`` to ``dest`` via a ``.tmp`` and atomic rename."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    shutil.copy2(src, tmp)
    tmp.replace(dest)


def _prepareref_dir(db: str) -> Path:
    """Resolve the ARIBA prepareref output dir for a given DB name."""
    return REPO_ROOT / "refs" / db / "prepareref_out"


def _run(
    *,
    db: str,
    accession: str,
    r1_url: str,
    r2_url: str,
    r1_md5: str,
    r2_md5: str,
    workdir: Path,
    run_dir: Path,
    threads: int,
    detailed: bool,
    ariba_timeout_min: int,
) -> int:
    """Execute the per-sample pipeline. Returns the process exit code."""
    reports_dir = run_dir / "reports"
    sample_logs_dir = run_dir / "sample_logs"
    report_dest = reports_dir / f"{accession}.report.tsv"

    if report_dest.exists():
        logger.info("already done — skipping (%s)", report_dest.relative_to(run_dir))
        return 0

    prep = _prepareref_dir(db)
    if not prep.is_dir():
        logger.error("prepareref dir missing: %s — build the ref DB first", prep)
        return 2

    workdir.mkdir(parents=True, exist_ok=True)
    sample_logs_dir.mkdir(parents=True, exist_ok=True)

    r1 = workdir / f"{accession}_1.fastq.gz"
    r2 = workdir / f"{accession}_2.fastq.gz"

    for url, dest, expected in ((r1_url, r1, r1_md5), (r2_url, r2, r2_md5)):
        rc = _curl_download(url, dest)
        if rc != 0:
            logger.error("curl failed (rc=%d) for %s", rc, url)
            return 3
        if not _verify_md5(dest, expected):
            return 4

    ariba_out = workdir / "ariba_out"
    if ariba_out.exists():
        shutil.rmtree(ariba_out)
    ariba_cmd = [
        "timeout",
        f"{ariba_timeout_min}m",
        "ariba",
        "run",
        "--threads",
        str(threads),
        str(prep),
        str(r1),
        str(r2),
        str(ariba_out),
    ]
    logger.info("ariba run (threads=%d, timeout=%dm)", threads, ariba_timeout_min)
    res = subprocess.run(ariba_cmd, check=False)
    if res.returncode != 0:
        logger.error("ariba run failed (rc=%d)", res.returncode)
        return 5

    report_src = ariba_out / ARIBA_REPORT
    if not report_src.is_file():
        logger.error("ariba completed but %s is missing", report_src)
        return 6
    _atomic_copy(report_src, report_dest)
    logger.info("wrote %s", report_dest.relative_to(run_dir))

    if detailed:
        for fname in ARIBA_DETAILED_FILES:
            src = ariba_out / fname
            if not src.is_file():
                logger.warning("--detailed requested but %s missing — skipping", fname)
                continue
            dest = run_dir / fname.split(".", 1)[0] / f"{accession}.{fname.split('.', 1)[1]}"
            # e.g. assembled_seqs.fa.gz -> <run_dir>/assembled_seqs/<acc>.fa.gz
            _atomic_copy(src, dest)
            logger.info("wrote %s", dest.relative_to(run_dir))

    shutil.rmtree(workdir, ignore_errors=True)
    return 0


def main() -> None:
    """CLI entrypoint."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--db", required=True, help="DB name under refs/<db>/prepareref_out/.")
    ap.add_argument("--accession", required=True)
    ap.add_argument("--r1-url", required=True)
    ap.add_argument("--r2-url", required=True)
    ap.add_argument("--r1-md5", required=True)
    ap.add_argument("--r2-md5", required=True)
    ap.add_argument("--workdir", type=Path, required=True, help="Per-task scratch dir (e.g. $SLURM_TMPDIR/<acc>).")
    ap.add_argument("--run-dir", type=Path, required=True, help="<RDS>/.../mag_rescue/<db>/<run-name>/")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--detailed", action="store_true", help="Also keep assembled_seqs.fa.gz + assemblies.fa.gz.")
    ap.add_argument("--ariba-timeout-min", type=int, default=90, help="Per-sample ariba run wall-time cap.")
    args = ap.parse_args()

    log_path = args.run_dir / "sample_logs" / f"{args.accession}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )

    rc = _run(
        db=args.db,
        accession=args.accession,
        r1_url=args.r1_url,
        r2_url=args.r2_url,
        r1_md5=args.r1_md5,
        r2_md5=args.r2_md5,
        workdir=args.workdir,
        run_dir=args.run_dir,
        threads=args.threads,
        detailed=args.detailed,
        ariba_timeout_min=args.ariba_timeout_min,
    )
    sys.exit(rc)


if __name__ == "__main__":
    main()
