"""Slurm orchestrator for the per-sample ARIBA array.

Builds the sbatch command for a `kleb_short_reads_<ver>.tsv` list file and
submits it (or prints it via ``--dry-run``). Runs on the login node — does
not itself process samples.

Subcommands
-----------
submit    — sbatch a fresh array (optionally capped via ``--n-samples``).
status    — tally per-state counts + cross-check report.tsv presence.
retry     — sbatch only the indices that failed in a prior job.

Usage
-----
    pixi run -e dev python -m mag_rescue.pp.parallel_ariba submit \
        --db kleb_virulence --run-name all \
        --mag-rescue-root <RDS>/processed/mag_rescue \
        --repo-dir ~/workspace/mag-rescue \
        [--n-samples 10] [--subset-metadata <vip-list>] \
        [--concurrency 100] [--dry-run]

    pixi run -e dev python -m mag_rescue.pp.parallel_ariba status --job-id 29154321
    pixi run -e dev python -m mag_rescue.pp.parallel_ariba retry  --job-id 29154321 \
        --db kleb_virulence --run-name all \
        --mag-rescue-root <RDS>/processed/mag_rescue \
        --repo-dir ~/workspace/mag-rescue
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger("parallel_ariba")

DEFAULT_CONCURRENCY = 100
DEFAULT_LIST_FILENAME = "kleb_short_reads_v1.tsv"
SBATCH_SCRIPT_RELPATH = "slurm_scripts/ariba_array.sh"

FAILED_STATES = {"FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED"}


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------


def _run_dir(mag_rescue_root: Path, db: str, run_name: str) -> Path:
    return mag_rescue_root / db / run_name


def _list_path(mag_rescue_root: Path, db: str, run_name: str, list_filename: str) -> Path:
    return _run_dir(mag_rescue_root, db, run_name) / "accessions" / list_filename


def _slurm_logs_dir(mag_rescue_root: Path) -> Path:
    return mag_rescue_root / "slurm_logs"


# ---------------------------------------------------------------------------
# sbatch command construction
# ---------------------------------------------------------------------------


def _count_data_rows(list_path: Path) -> int:
    """Count rows in the list TSV, excluding the header."""
    with list_path.open() as fh:
        next(fh, None)  # header
        return sum(1 for _ in fh)


def _build_sbatch_cmd(
    *,
    repo_dir: Path,
    list_path: Path,
    db: str,
    run_dir: Path,
    slurm_logs_dir: Path,
    array_spec: str,
    ariba_sif: Path,
    subset_metadata: Path | None,
) -> list[str]:
    """Compose the sbatch command. Returns the argv list (subprocess-ready)."""
    sbatch_script = repo_dir / SBATCH_SCRIPT_RELPATH
    output_pat = slurm_logs_dir / "mag-ariba_%A_%a.out"
    error_pat = slurm_logs_dir / "mag-ariba_%A_%a.err"
    export = (
        f"ALL,LIST_FILE={list_path},DB={db},RUN_DIR={run_dir},"
        f"REPO_DIR={repo_dir},ARIBA_SIF={ariba_sif},"
        f"SUBSET_METADATA={subset_metadata or ''}"
    )
    return [
        "sbatch",
        f"--array={array_spec}",
        f"--output={output_pat}",
        f"--error={error_pat}",
        f"--export={export}",
        str(sbatch_script),
    ]


def _array_spec_range(n_total: int, n_samples: int, concurrency: int) -> str:
    """Build a `1-N%M` array spec, capped to ``n_samples`` if positive."""
    n = n_total if n_samples <= 0 else min(n_samples, n_total)
    return f"1-{n}%{concurrency}"


def _array_spec_indices(indices: list[int], concurrency: int) -> str:
    """Build a Slurm array spec from a list of indices, collapsing consecutive runs.

    Consecutive indices become `N-M` ranges (e.g. `2272-2277`) instead of being
    enumerated individually. CSD3 transient curl failures during high-concurrency
    ENA fetches tend to cluster in long runs, so the unfolded comma-list can
    blow past Slurm's argv length limit (~4 KB) — at which point sbatch returns
    "Pathname of a file, directory or other parameter too long". Range collapse
    typically halves the spec length and keeps the same task set.
    """
    if not indices:
        return f"%{concurrency}"
    sorted_idx = sorted(set(indices))
    runs: list[str] = []
    start = prev = sorted_idx[0]
    for i in sorted_idx[1:]:
        if i == prev + 1:
            prev = i
            continue
        runs.append(f"{start}-{prev}" if prev > start else str(start))
        start = prev = i
    runs.append(f"{start}-{prev}" if prev > start else str(start))
    return f"{','.join(runs)}%{concurrency}"


def _parse_jobid(stdout: str) -> str | None:
    """Pull `Submitted batch job <N>` out of sbatch's stdout."""
    m = re.search(r"Submitted batch job (\d+)", stdout)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# sacct queries
# ---------------------------------------------------------------------------


def _sacct_states(job_id: str) -> list[tuple[int, str]]:
    """Return [(array_index, state), ...] for the given job, via sacct.

    Filters to array tasks (rows with ``JobID == <job_id>_<index>``); skips
    .batch / .extern child rows and the parent stem.
    """
    res = subprocess.run(
        ["sacct", "-j", job_id, "--format=JobID,State", "-X", "--parsable2", "--noheader"],
        check=False,
        capture_output=True,
        text=True,
    )
    if res.returncode != 0:
        logger.error("sacct failed: %s", res.stderr.strip())
        return []
    out: list[tuple[int, str]] = []
    pat = re.compile(rf"^{re.escape(job_id)}_(\d+)$")
    for line in res.stdout.splitlines():
        if not line.strip():
            continue
        jobid, state = line.split("|", 1)
        m = pat.match(jobid)
        if not m:
            continue
        out.append((int(m.group(1)), state.split()[0]))  # state may have suffix
    return out


def _failed_indices(job_id: str) -> list[int]:
    return sorted(idx for idx, st in _sacct_states(job_id) if st in FAILED_STATES)


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------


def cmd_submit(args: argparse.Namespace) -> int:
    """Sbatch a fresh array."""
    list_path = _list_path(args.mag_rescue_root, args.db, args.run_name, args.list_filename)
    if not list_path.is_file():
        logger.error("list file missing: %s", list_path)
        return 1
    if not args.ariba_sif.is_file():
        logger.error("ariba SIF missing: %s", args.ariba_sif)
        return 1
    n_total = _count_data_rows(list_path)
    array_spec = _array_spec_range(n_total, args.n_samples, args.concurrency)
    cmd = _build_sbatch_cmd(
        repo_dir=args.repo_dir,
        list_path=list_path,
        db=args.db,
        run_dir=_run_dir(args.mag_rescue_root, args.db, args.run_name),
        slurm_logs_dir=_slurm_logs_dir(args.mag_rescue_root),
        array_spec=array_spec,
        ariba_sif=args.ariba_sif,
        subset_metadata=args.subset_metadata,
    )
    logger.info("array spec: %s  (n_total=%d, n_samples=%d)", array_spec, n_total, args.n_samples)
    logger.info("sbatch cmd:\n  %s", " \\\n  ".join(cmd))

    if args.dry_run:
        logger.info("--dry-run set; not submitting")
        return 0
    _slurm_logs_dir(args.mag_rescue_root).mkdir(parents=True, exist_ok=True)
    res = subprocess.run(cmd, check=False, capture_output=True, text=True)
    sys.stdout.write(res.stdout)
    sys.stderr.write(res.stderr)
    if res.returncode != 0:
        return res.returncode
    job_id = _parse_jobid(res.stdout)
    if job_id:
        logger.info("submitted job %s — track with `parallel_ariba status --job-id %s`", job_id, job_id)
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    """Tally per-state counts for a Slurm array, cross-check report.tsv presence."""
    states = _sacct_states(args.job_id)
    if not states:
        logger.error("no array tasks found for job %s", args.job_id)
        return 1
    counts = Counter(st for _, st in states)
    n = sum(counts.values())
    print(f"job {args.job_id}: {n} array tasks")
    for state, c in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {state:12s} {c}")

    if args.run_dir is not None:
        reports = args.run_dir / "reports"
        if reports.is_dir():
            n_reports = sum(1 for _ in reports.glob("*.report.tsv"))
            print(f"  reports/   {n_reports} files in {reports}")
    return 0


def cmd_retry(args: argparse.Namespace) -> int:
    """Re-submit only the indices that failed in a prior job."""
    indices = _failed_indices(args.job_id)
    if not indices:
        logger.info("no failed indices for job %s — nothing to retry", args.job_id)
        return 0
    list_path = _list_path(args.mag_rescue_root, args.db, args.run_name, args.list_filename)
    if not list_path.is_file():
        logger.error("list file missing: %s", list_path)
        return 1
    array_spec = _array_spec_indices(indices, args.concurrency)
    cmd = _build_sbatch_cmd(
        repo_dir=args.repo_dir,
        list_path=list_path,
        db=args.db,
        run_dir=_run_dir(args.mag_rescue_root, args.db, args.run_name),
        slurm_logs_dir=_slurm_logs_dir(args.mag_rescue_root),
        array_spec=array_spec,
        ariba_sif=args.ariba_sif,
        subset_metadata=args.subset_metadata,
    )
    logger.info("retry %d failed indices from job %s", len(indices), args.job_id)
    logger.info("array spec: %s", array_spec)
    logger.info("sbatch cmd:\n  %s", " \\\n  ".join(cmd))
    if args.dry_run:
        return 0
    res = subprocess.run(cmd, check=False, capture_output=True, text=True)
    sys.stdout.write(res.stdout)
    sys.stderr.write(res.stderr)
    return res.returncode


# ---------------------------------------------------------------------------
# argparse
# ---------------------------------------------------------------------------


def _add_shared_run_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--db", default="kleb_virulence", help="DB name (default: kleb_virulence).")
    p.add_argument("--run-name", default="all", help="Run name under <db>/ (default: all).")
    p.add_argument("--mag-rescue-root", type=Path, required=True, help="<RDS>/processed/mag_rescue/")
    p.add_argument("--repo-dir", type=Path, required=True, help="Path to the cloned mag-rescue repo on HPC.")
    p.add_argument("--ariba-sif", type=Path, required=True, help="Path to the ariba apptainer container.")
    p.add_argument("--list-filename", default=DEFAULT_LIST_FILENAME)
    p.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    p.add_argument("--subset-metadata", type=Path, default=None)
    p.add_argument("--dry-run", action="store_true")


def main() -> None:
    """CLI entrypoint."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_submit = sub.add_parser("submit", help="Submit a new array.")
    _add_shared_run_args(p_submit)
    p_submit.add_argument("--n-samples", type=int, default=-1, help="Cap to first N (default -1 = all).")
    p_submit.set_defaults(func=cmd_submit)

    p_status = sub.add_parser("status", help="Tally per-state counts for a job.")
    p_status.add_argument("--job-id", required=True)
    p_status.add_argument("--run-dir", type=Path, default=None, help="Optional: count reports/ files too.")
    p_status.set_defaults(func=cmd_status)

    p_retry = sub.add_parser("retry", help="Re-submit failed indices from a prior job.")
    _add_shared_run_args(p_retry)
    p_retry.add_argument("--job-id", required=True)
    p_retry.set_defaults(func=cmd_retry)

    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")

    if args.cmd in {"submit", "retry"} and shutil.which("sbatch") is None:
        logger.warning("`sbatch` not on PATH — submit/retry will fail. (Use --dry-run to compose only.)")

    sys.exit(args.func(args))


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# Helpers exposed for testing
# ---------------------------------------------------------------------------


def _write_list_for_test(path: Path, n_data_rows: int) -> Path:
    """Test helper: create a list TSV with header + n data rows."""
    rows = [["run_accession", "r1_url", "r2_url", "r1_md5", "r2_md5"]]
    for i in range(n_data_rows):
        rows.append([f"SRR{i:06d}", f"https://x/{i}_1.fq.gz", f"https://x/{i}_2.fq.gz", "a", "b"])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    return path
